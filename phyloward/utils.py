import argparse
import subprocess
import sys
import warnings
from functools import wraps
from collections import OrderedDict

import re


def coroutine_simple(func):  # coroutine decorator
    @wraps(func)
    def start(*args, **kwargs):
        cr = func(*args, **kwargs)
        next(cr)
        return cr
    return start


@coroutine_simple
def progress(total=None, head_message=''):  # progress indicator in CLI
    if not total or total < 0:
        total = float('inf')
        string = '[{prog}] '
    else:
        total = int(total)
        digits = len(str(total))
        string = '[{{prog:{digits}d}}/{tot:{digits}d}] '.format(digits=digits, tot=total)

    i = 0
    while i < total:
        try:
            print(head_message, string.format(prog=i),
                  sep='', end='                    \r', file=sys.stderr, flush=True)
            i += 1
            _ = yield
        except StopIteration:
            break
    print(head_message, string.format(prog=i), '                  ', sep='', file=sys.stderr)
    yield  # prevents raising StopIteration as soon as finished


# --- FASTA file IO ---


class InvalidFastaFormat(ValueError):
    pass


class SimpleFastaReader(object):
    def __init__(self, handle):
        if hasattr(handle, 'readable') and handle.readable():
            self.handle = handle
        else:
            raise TypeError

    def __iter__(self):
        self.handle.seek(0)
        label = None
        frags = []
        for line in self.handle:
            if line.isspace() or line.startswith(';'):
                continue
            if line.startswith('>'):
                if label and isinstance(label, str):
                    yield (label, ''.join(frags))
                label = line[1:].split(maxsplit=1)[0]  # description is ignored
                frags = []
            elif label:
                frags.append(line.strip())
        if label and isinstance(label, str):
            yield (label, ''.join(frags))
        else:
            raise InvalidFastaFormat
        return


class SimpleFastaWriter(object):
    def __init__(self, d, width=None):
        self.data = OrderedDict(d)
        self.width = width

    @staticmethod
    def _split_lines(s, length=None):
        if length is None:
            return s
        return '\n'.join(s[i:i + length] for i in range(0, len(s), length))

    def __iter__(self):
        for label, s in self.data.items():
            # if re.search(r'\s', label):
            #     label = '\'' + label + '\''
            yield label, (s if s else '-')

    def __str__(self):
        frags = ['>{0}\n{1}\n'.format(label, self._split_lines(s, self.width))
                 for label, s in self]
        return ''.join(frags)


# --- Abstract class for external programs ---

class _BaseExternalPipe(object):
    PROG_NAME = NotImplemented
    exec_path = NotImplemented
    supports_stdin = NotImplemented  # TODO cls property?

    def __init__(self):
        self._stdin = None
        self._stdout = None
        self._stderr = None
        self._args_log = None

    @property
    def stdin(self):
        return self._stdin

    @stdin.setter
    def stdin(self, s):
        if not self.supports_stdin and s:
            warnings.warn('This program does not support stdin.')
        self._stdout = None
        self._stdin = str(s)

    @property
    def stdout(self):
        return self._stdout

    @property
    def stderr(self):
        return self._stderr

    @property
    def args_log(self):
        return self._args_log

    @staticmethod
    def _split_args(args):
        if args is None:
            return []
        try:
            return args.split()
        except AttributeError:
            return list(args)

    # Get version info from out or err
    @classmethod
    def get_version(cls):
        # out, err = cls._process_executor(args='#version_arg')
        raise NotImplementedError('Subclass should implement this')

    # Default option parameters
    def _options(self, *args, **kwargs):
        options = [NotImplemented, ]
        raise NotImplementedError('Subclass should implement this')
        return options

    def run(self, override_args=None, cwd=None):
        if override_args:
            options = self._split_args(override_args)
        else:
            options = self._options()
        self._args_log = options
        self._stdout, self._stderr = \
            self._execute_process(self.stdin, options=options, cwd=cwd)
        return self._stdout

    @classmethod
    def _execute_process(cls, stdin=None, options=None, cwd=None, timeout=None):
        if options is None:
            options = []
        try:  # args is string
            options = options.split()
        except AttributeError:  # args is list
            pass

        args = [cls.exec_path] + options
        proc = subprocess.Popen(args,
                                stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                universal_newlines=True, cwd=cwd)
        try:
            out, err = proc.communicate(input=stdin, timeout=timeout)  # type: str
        except subprocess.TimeoutExpired:
            proc.kill()
            raise TimeoutError('Timeout Error from {}!\n'.format(cls.PROG_NAME))
        if proc.returncode != 0:
            raise RuntimeError('Runtime Error during executing {}...\n'.format(cls.PROG_NAME) + err)
        return out, err


# --- Argument parsing classes overridden ---

class AliasedAction(argparse._StoreAction):
    def __init__(self, *args, **kwargs):
        self.choice_alias = kwargs.get('alias_map', {})
        try:
            del kwargs['alias_map']
        except KeyError:
            pass
        super().__init__(*args, **kwargs)


class PatchedHelpFormatter(argparse.HelpFormatter):
    def add_usage(self, usage, actions, groups, prefix=None):
        return super().add_usage(usage, actions, groups, prefix)

    def _format_usage(self, usage, actions, groups, prefix):
        prefix = 'Usage: '
        # prefix = ''
        return super()._format_usage(usage, actions, groups, prefix)

    def _metavar_formatter(self, action, default_metavar):
        return super()._metavar_formatter(action, default_metavar)

    def _format_action_invocation(self, action):
        if not action.option_strings:
            default = self._get_default_metavar_for_positional(action)
            metavar, = self._metavar_formatter(action, default)(1)
            return metavar
        else:
            parts = []
            # Patched: -s/--long  # Original: -s, --long
            if action.nargs == 0:
                parts.extend(action.option_strings)
                result = '/'.join(parts)
            # Patched: -s/--long ARGS  # Original: -s ARGS, --long ARGS
            else:
                default = self._get_default_metavar_for_optional(action)
                args_string = self._format_args(action, default)
                for option_string in action.option_strings:
                    parts.append('%s' % (option_string))
                result = '/'.join(parts) + ' ' + args_string
            return result


class PatchedArgParser(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        if not kwargs.get('formatter_class'):
            kwargs['formatter_class'] = PatchedHelpFormatter
        super().__init__(*args, **kwargs)
        self.subparser_dest = None

    def parse_args(self, args=None, namespace=None):
        args, argv = self.parse_known_args(args, namespace)
        if argv:
            msg = argparse._('unrecognized arguments: %s')
            if args.__getattribute__(self.subparser_dest) in self.subparser_names():
                subparser = self.get_subcommand_parser(args.__getattribute__(self.subparser_dest))
                subparser.error(msg % ' '.join(argv))
            self.error(msg % ' '.join(argv))
        return args

    def add_subparsers(self, **kwargs):
        if self._subparsers is not None:
            self.error(argparse._('cannot have multiple subparser arguments'))
        if not kwargs.get('parser_class'):
            kwargs['parser_class'] = self.__class__
        if kwargs.get('dest'):
            self.subparser_dest = kwargs['dest']
        return super().add_subparsers(**kwargs)

    def error(self, message):
        self.print_help()
        sys.stderr.write('\n{prog}: error: {message}\n'.format(prog=self.prog, message=message))
        sys.exit(2)

    def _get_value(self, action, arg_string):
        result = super()._get_value(action, arg_string)
        try:
            aliased_result = action.choice_alias.get(result, result)
        except (AttributeError, NameError):
            aliased_result = result
        return aliased_result

    def _check_value(self, action, value):
        try:
            aliased_value = action.choice_alias.get(value, value)
        except (AttributeError, NameError):
            aliased_value = value
        if action.choices is not None and aliased_value not in action.choices:
            args = {'value': value,
                    'choices': ', '.join(map(repr, action.choices))}
            msg = argparse._('invalid choice: %(value)r (choose from %(choices)s)')
            raise argparse.ArgumentError(action, msg % args)

    def add_argument(self, *args, **kwargs):
        return super().add_argument(*args, **kwargs)

    def get_subcommand_parser(self, name):
        ag = self._subparsers  # type: argparse._ArgumentGroup
        action = ag._group_actions[0]
        parser = action._name_parser_map[name]  # type: argparse.ArgumentParser
        return parser

    def subparser_names(self):
        return self._subparsers._group_actions[0]._name_parser_map.keys()


# class Spinner(object):
#     busy = False
#     delay = 0.125
#
#     @staticmethod
#     def spinning_cursor():
#         while True:
#             for cursor in '|/-\\':
#                 yield cursor
#
#     def __init__(self, delay=None):
#         self.spinner_generator = self.spinning_cursor()
#         self.end_message = ' '
#         if delay and float(delay):
#             self.delay = delay
#         self.start()
#
#         def _handle_keyboard_interrupt(signal, frame):
#             self.stop()
#             raise KeyboardInterrupt
#
#         import signal
#         signal.signal(signal.SIGINT, _handle_keyboard_interrupt)
#
#     def spinner_task(self):
#         sys.stderr.write(' ')
#         while self.busy:
#             sys.stderr.write('\b')
#             sys.stderr.flush()
#             sys.stderr.write(next(self.spinner_generator))
#             sys.stderr.flush()
#             time.sleep(self.delay)
#
#     def start(self):
#         self.busy = True
#         import threading
#         threading.Thread(target=self.spinner_task).start()
#
#     def stop(self):
#         self.busy = False
#         # time.sleep(self.delay)
#
#     def end(self, message=' '):
#         self.stop()
#         sys.stderr.write('\b')
#         sys.stderr.write(message)
#         sys.stderr.write('\n')
#         sys.stderr.flush()
