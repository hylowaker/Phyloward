#!/usr/bin/env python3
"""
Command line python script for Phyloward
"""
import concurrent.futures.thread
import json
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

from phyloward import __version__ as ver, extractor, aligner
from phyloward.utils import PatchedArgParser, AliasedAction, coroutine_simple, progress

PACKAGE_NAME = 'phyloward'


def _print_welcome():
    welcome = (
        '------------------------------------------------\n'
        ' Phyloward version {ver}\n'
        ' Copyright (c) 2018 JaeHeung Han\n'
        '------------------------------------------------\n'
    ).format(ver=ver)
    print(welcome, file=sys.stderr)


def _argument_parser():
    # --------------- Argparse: subcommand ---------------
    usage = PACKAGE_NAME + ' {extract,align} ...'
    parser = PatchedArgParser(prog=PACKAGE_NAME, usage=usage)  # TODO better help
    parser.add_argument('-v', '--version', action='store_true')
    subparsers = parser.add_subparsers(title='Subcommands', dest='command', prog=PACKAGE_NAME)

    usage_extract = PACKAGE_NAME + ' extract [options..] <INPUT> <OUT_DIR>'
    parser_extract = subparsers.add_parser('extract', usage=usage_extract,
                                           help='[Step 1] Search core genes from genome')  # type: PatchedArgParser

    usage_align = PACKAGE_NAME + ' align [-t] [options..] <DIRECTORY> <RESULT_DIR>'
    parser_align = subparsers.add_parser('align', usage=usage_align,
                                         help='[Step 2] Align each core gene and concatenate to infer phylogeny')  # type: PatchedArgParser

    # --------------- Argparse: extract ---------------
    parser_extract.add_argument('input', metavar='INPUT',
                                help='A single genome FASTA file '
                                     'OR a directory containing multiple FASTA files')
    parser_extract.add_argument(metavar='OUT_DIR', dest='out', nargs='?',
                                help='Directory where extracted core gene JSON files '
                                     'should be stored')
    parser_extract.add_argument('-a', '--archaea', action='store_true',
                                help='ARCHAEA mode')
    parser_extract.add_argument('-c', '--custom', metavar='HMM',
                                help='Use custom HMM profile')
    parser_extract.add_argument('-l', '--label', metavar='S',
                                help='Name of the genome (Ignored if INPUT is directory)')
    parser_extract.add_argument('-p', '--process', metavar='N', type=int,
                                help='Numbers of simultaneous workers (default: Use all available cores)')
    # parser_extract.add_argument('--table', metavar='N', default='11',
    #                             help='Translation table for Prodigal (default: 11)')
    # TODO metadata args

    # --------------- Argparse: align ---------------
    parser_align.add_argument(metavar='DIRECTORY', dest='input_dir', nargs='?',
                              help='Directory which contains extracted core gene JSON files')
    parser_align.add_argument(metavar='RESULT_DIR', dest='result_dir', nargs='?',
                              help='Directory where result files to be stored. '
                                   'If not given, result JSON string is printed as stdout.')
    parser_align.add_argument('-t', '--tree', action='store_true',
                              help='Include phylogenetic trees in the result. '
                                   'The tree is represented as Newick format, process by FastTree.')
    parser_align.add_argument('-g', '--full-gene', action='store_true',
                              help='Include whole gene sequence instead of trimmed domain region. '
                                   'Note that the result may include duplicated sequences if you enable this option.')
    _choice_alias = {'nt': 'nucleotide', 'nucl': 'nucleotide',
                     'amino': 'protein', 'aa': 'protein', 'prot': 'protein'}
    parser_align.add_argument('-s', '--seq-type', metavar='S', action=AliasedAction,
                              choices=['nucleotide', 'protein', 'codon', 'codon12'],
                              alias_map=_choice_alias,
                              default='aa', help='Type of aligned sequence to use on tree building. '
                                                 'Choose from \'nucleotide\', \'protein\', \'codon\', \'codon12\'.'
                                                 ' (default: protein)')
    parser_align.add_argument('-f', '--filter', metavar='D', type=float, default=0.,
                              help='Remove gap-rich columns from concatenated sequence. '
                                   'The threshold is a percentage of gaps in a column. '
                                   'The higher the value, more columns will be filtered out. '
                                   '(range: 0~100, default: 0)')
    parser_align.add_argument('-p', '--process', metavar='N', type=int,
                              help='Numbers of simultaneous workers (default: Use all available cores)')
    parser_align.add_argument('--gsi-threshold', metavar='D', type=float, default=95.,
                              help='# experimental')

    return parser


parser = None
args = None
argv = None  # for test/debugging


# =========================================================================

def main():
    global parser
    global args
    parser = _argument_parser()
    args = parser.parse_args(argv)
    if args.version:
        print('Phyloward pipeline version:', ver, file=sys.stderr)
        return

    if 'extract' == args.command:
        _pipe_extract()
    elif 'align' == args.command:
        _pipe_align()
    else:
        _print_welcome()
        parser.print_help()

# =========================================================================


def _extract_single(file, label=None):
    try:
        # Main Part
        cg_extracted = extractor.extract_core_genes(
            file, is_archaea=args.archaea, profile=args.custom, label=label,
            # codon_table=args.table
        )
    except RuntimeError:
        raise

    return cg_extracted


@coroutine_simple
def _writes_extracted(directory=None, indent=4):
    try:
        os.makedirs(directory, exist_ok=True)
    except TypeError:
        print('# OUT_DIR argument is not provided. Writing to stdout instead...',
              file=sys.stderr)
    except FileExistsError as e:
        print('# ERROR:', e, file=sys.stderr)
        sys.exit(1)

    count = 0
    cg_extracted_list = []
    while True:
        try:
            cg_extracted = yield  # type: extractor.ExtractedCoreGenes
            if not cg_extracted:
                break
        except StopIteration:
            break

        if directory:
            _base = os.path.basename(cg_extracted.info['extracted_from'])  # abspath :--> basename
            _name = os.path.splitext(_base)[0] + '_extracted.json'
            target = os.path.join(directory, _name)
            with open(target, 'w') as fh:
                print(cg_extracted.as_json(indent), file=fh)
        else:
            cg_extracted_list.append(cg_extracted)
        count += 1

    if directory:
        full_dir = os.path.abspath(directory)
        if count > 1:
            print('Extracted JSON files are saved in "{}" directory.'.format(full_dir),
                  file=sys.stderr)
        elif count == 1:
            print('Extracted JSON file is saved in "{}" directory.'.format(full_dir),
                  file=sys.stderr)
    else:
        if len(cg_extracted_list) > 1:
            print(json.dumps(cg_extracted_list, indent=indent))
        elif len(cg_extracted_list) == 1:
            print(json.dumps(cg_extracted_list[0], indent=indent))
        else:
            print('{}')


def _pipe_extract():
    args_input = args.input
    outdir = args.out
    workers = args.process if args.process else (os.cpu_count() if os.cpu_count() else 1)
    pref_dom = 'ARCHAEAL' if args.archaea else 'BACTERIAL'
    print('# Custom HMM profile: ' + args.custom, file=sys.stderr)

    writer = _writes_extracted(outdir, indent=4)

    is_multiple = os.path.isdir(args_input)
    if is_multiple:
        # process multiple fasta files in a directory
        input_dir = os.path.abspath(args_input)
        _inputs = (x for x in os.scandir(input_dir) if x.is_file())
        _inputs = sorted(_inputs, key=lambda x: x.name)
        _inputs = [x.path for x in _inputs]
        input_list = [x for x in _inputs
                      if x.endswith('.fasta') or x.endswith('.fna') or x.endswith('fa')]
        if not input_list:
            print('# ERROR: No FASTA file were found in "{}" directory'.format(args_input),
                  file=sys.stderr)
            sys.exit(1)
        reporter = progress(len(input_list),
                            'Extracting {} core genes from {} files... '
                            .format(pref_dom, len(input_list)))
    elif os.path.exists(args_input):
        # extract from single fasta file
        input_list = [args_input]
        print('Extracting {} core genes from file "{}"...'
              .format(pref_dom, os.path.basename(args_input)),
              file=sys.stderr)
        reporter = None
    else:
        print('# ERROR: Cannot access "{}"'.format(args_input), file=sys.stderr)
        sys.exit(2)

    with ThreadPoolExecutor(workers) as pool:
        futures = [pool.submit(_extract_single, file) for file in input_list]
        try:
            for future in as_completed(futures):
                cg_extracted = future.result()
                writer.send(cg_extracted)  # If args.out, writes a file. Else, does nothing
                try:
                    reporter.send('')
                except (AttributeError, StopIteration):
                    pass
        except (Exception, KeyboardInterrupt):
            pool._threads.clear()
            concurrent.futures.thread._threads_queues.clear()
            raise

    try:
        writer.throw(StopIteration)
    except StopIteration:
        pass


# =========================================================================


def _pipe_align():
    # set parameters ------------------------------------------------------
    input_dir = args.input_dir
    result_dir = args.result_dir
    domain_only = not args.full_gene
    seq_type = args.seq_type
    workers = args.process if args.process else (os.cpu_count() if os.cpu_count() else 1)
    cutoff = args.filter / 100

    print('# Draw tree: {}'.format('YES' if args.tree else 'NO'), file=sys.stderr)
    print('# Use whole gene: {}'.format('NO' if domain_only else 'YES'), file=sys.stderr)
    print('# Sequence type: {}'.format(seq_type), file=sys.stderr)

    if not 0 <= cutoff <= 1.:
        print('# ERROR: filter parameter must be in range [0, 100]', file=sys.stderr)
        sys.exit(2)
    gsi_threshold = args.gsi_threshold / 100

    if sys.stdin.isatty() and not input_dir:  # no input arguments --> print help
        parser_align = parser.get_subcommand_parser('align')
        parser_align.error('the following arguments are required: DIRECTORY')
        return

    # get inputs ----------------------------------------------------------
    if not sys.stdin.isatty():  # alternative input: read from stdin pipe
        print('# WARNING: Reading from stdin. This feature is experimental.', file=sys.stderr)
        cg_extracted_list = extractor.extracted_list_from_json(raw_string=sys.stdin.read())
        cg_aligned = aligner.import_from_extracted_list(cg_extracted_list)

    elif not os.path.isdir(input_dir):  # alternative input: read from joined JSON file
        infile = input_dir
        print('# WARNING! "{}" is not a directory'.format(infile), file=sys.stderr)
        try:
            cg_aligned = aligner.import_from_files([infile, ])
            print('# Assuming "{}" contains multiple extracted results...'.format(infile),
                  file=sys.stderr)
        except (TypeError, ValueError):
            print('# ERROR: Cannot read input files', file=sys.stderr)
            sys.exit(2)

    else:  # Normally goes this route. Read from files
        extracted_dir = os.path.abspath(input_dir)
        _inputs = sorted((x for x in os.scandir(extracted_dir)
                          if x.is_file() and x.name.endswith('.json')),
                         key=lambda x: x.name)
        input_list = [x.path for x in _inputs]
        if not input_list:
            print('# ERROR: Directory "{}" is empty'.format(input_dir), file=sys.stderr)
            sys.exit(2)
        # cg_extracted_list = aligner._extracted_list_from_json_files(input_list)
        cg_aligned = aligner.import_from_files(input_list)

    if not result_dir:
        print('# RESULT_DIR argument is not provided. Writing to stdout instead...',
              file=sys.stderr)

    # Main align process! -------------------------------------------------
    domain_msg = 'ARCHAEAL' if cg_aligned.is_archaeal else 'BACTERIAL'
    print('# Assuming all {} inputs are {}...'.format(cg_aligned.number_of_genomes(), domain_msg),
          file=sys.stderr)
    isnt = True if seq_type in ('nucleotide', 'nt', 'nucl') else False
    cg_aligned.align(domain_only=domain_only, nucleotide=isnt, cutoff=cutoff, workers=workers, report_progress=progress)  # time-consuming
    info_alignment_program = cg_aligned.msa_program_info  # comment
    cg_aligned.add_comment('arguments', ' '.join(sys.argv))
    cg_aligned.add_comment('superkingdom', 'Archaea' if cg_aligned.is_archaeal else 'Bacteria')
    cg_aligned.add_comment('alignment_program', str(info_alignment_program))

    # fill trees ----------------------------------------------------------
    if args.tree:  # Without '--tree' option, result only includes alignments.
        # workers = 3  # TODO raxml debug test
        cg_aligned.fill_trees(seq_type, workers=workers, report_progress=progress)
        cg_aligned.add_comment('phylogeny_program', str(cg_aligned.tree_program_info))
        cg_aligned.calculate_gene_support_index(threshold=gsi_threshold)

    # write final results -------------------------------------------------
    cg_aligned.add_comment('created_on', time.strftime('%b %d %Y %H:%M:%S'))
    cg_aligned.add_comment('{}_pipeline_version'.format(PACKAGE_NAME), ver)
    if not result_dir:  # print result as stdout
        print(cg_aligned.result_as_json(indent=4))
    else:
        _dump_aligned_results(cg_aligned, result_dir)
        print('Result files are saved in "{}" directory'.format(os.path.abspath(result_dir)),
              file=sys.stderr)
    # TODO done message


def _dump_aligned_results(result, outdir, label=None):
    result = result  # type: aligner.CoreGenesNewickTree
    os.makedirs(outdir, exist_ok=True)
    outdir = os.path.abspath(outdir)
    # prefix = label if label else os.path.basename(outdir)
    currdir = os.getcwd()
    os.chdir(outdir)

    # Full result JSON
    with open('FullResult.json', 'w') as fh:
        print(result.result_as_json(indent=4), file=fh)

    # Alignment FASTA (concatenated core genes)
    with open('Alignment_CoreGenes.fasta', 'w') as fh:
        print(result.get_fasta(), file=fh, end='')

    # Alignment FASTA (each core gene)
    for gene in result.genes():
        assert gene not in result._excluded
        try:
            s = result.get_fasta(gene)
        except KeyError:
            continue
        with open('Alignment_' + gene + '.fasta', 'w') as fh:  # beware characters \/:*?"<>|
            print(s, file=fh, end='')

    # logfile
    with open('Log.log', 'w') as fh:
        print('Genomes included in the analysis.', file=fh)
        print('uid', 'label', sep='\t', file=fh)
        for uid, label in result.uid_label.items():
            print(uid, label, sep='\t', file=fh)
        print(file=fh)
        print('The length of alignments', file=fh)
        for gene, seqdict in result._data.items():
            try:
                len_ = seqdict.length_of_alignment()
            except aligner.NotAlignedError:
                len_ = '#excluded'
            print(gene, len_, sep='\t', file=fh)

    if not args.tree:  # stop here if not drawing trees
        os.chdir(currdir)
        return

    # Newick tree (concatenated core genes)
    with open('Tree.nwk', 'w') as fh:
        print(result.get_tree(), file=fh)

    # Newick tree GSI (concatenated core genes)
    with open('Tree(GSI-marked).nwk', 'w') as fh:
        print(result.get_tree_gsi(), file=fh)

    # Newick tree (each core gene)
    for gene in result.genes():
        assert gene not in result._excluded
        try:
            s = result.get_tree(gene)
        except KeyError:
            continue
        with open('Tree_' + gene + '.nwk', 'w') as fh:
            print(s, file=fh)

    os.chdir(currdir)

# =========================================================================


if __name__ == '__main__':
    main()
