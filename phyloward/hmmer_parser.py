# Modified by JaeHeung Han <hylowaker@gmail.com>
# The original code is written by Wibowo Arindrarto.
# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its license.
# ------------------------------------------------------------------
# Biopython License Agreement
#
# Permission to use, copy, modify, and distribute this software and its
# documentation with or without modifications and for any purpose and
# without fee is hereby granted, provided that any copyright notices
# appear in all copies and that both those copyright notices and
# this permission notice appear in supporting documentation, and that
# the names of the contributors or copyright holders not be used
# in advertising or publicity pertaining to distribution of the software
# without specific prior permission.
#
# THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM
# ALL WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
# CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL,
# INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING
# FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT,
# NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION
# WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.


import re
from collections import namedtuple

# precompile regex patterns
_RE_PROGRAM = re.compile(r'^# (\w*hmm\w+) :: .*$')  # Program name capture
_RE_VERSION = re.compile(r'# \w+ ([\w+.]+) .*; http.*$')  # Version string capture
_RE_OPT = re.compile(r'^# (.+):\s+(.+)$')  # Option string capture
_QUERY_ID_LEN = re.compile(r'^Query:\s*(.*)\s+\[\w=(\d+)\]')  # Query id and length capture
_DOMHIT_PREFX = re.compile(r'^\s+\d')


QueryResult = namedtuple('QueryResult', ['id', 'desc', 'hits'])


class HmmerParser(object):
    def __init__(self, stream):
        self.handle = stream  # TODO handle non-blocking live stream
        self._line = None
        self._read_forward()
        self.metadata = self._parse_preamble()
        # TODO line number?

    def __iter__(self):
        self.handle.seek(0)
        self._line = None
        for qresult in self._parse_query_generator():
            yield qresult

    def _read_forward(self):
        while True:
            # Read next non-whitespace line
            line = self.handle.readline()
            if (line and line.strip()) or (not line):
                self._line = line
                return

    def _read_next_until(self, bool_func):
        while True:
            self._read_forward()
            if bool_func(self._line) or not self._line:
                return

    def _parse_preamble(self):
        metadata = {}
        # Bool flag for storing state ~ whether we are parsing the option lines or not
        has_opts = False
        while True:
            # No pound sign means we've left the preamble
            if not self._line.startswith('#'):
                break
            # Dashes could either mean we are entering or leaving the options
            # Section ~ so it's a switch for the has_opts flag
            elif '- - -' in self._line:
                if not has_opts:
                    # If flag is false, that means we're entering opts so switch the flag accordingly
                    has_opts = True
                else:
                    # If flag is true, that means we've reached the end of opts so we can break out of the function
                    break
            elif not has_opts:
                regx = re.search(_RE_PROGRAM, self._line)
                if regx:
                    metadata['program'] = regx.group(1)
                regx = re.search(_RE_VERSION, self._line)
                if regx:
                    metadata['version'] = regx.group(1)
            elif has_opts:
                regx = re.search(_RE_OPT, self._line)
                # if target in regx.group(1), then we store the key as target
                if 'target' in regx.group(1):
                    metadata['target'] = regx.group(2).strip()
                else:
                    metadata[regx.group(1)] = regx.group(2)

            self._read_forward()

        return metadata

    def _parse_query_generator(self):
        """Parse a HMMER3 query block (PRIVATE)."""
        self._read_next_until(lambda line: line.startswith('Query:'))

        while self._line:
            regx = re.search(_QUERY_ID_LEN, self._line)
            while not regx:
                self._read_forward()
                regx = re.search(_QUERY_ID_LEN, self._line)

            # Get query id and length
            query = regx.group(1).strip()
            # Store qresult attributes  # TODO make use of additional info?
            # qresult_attrs = {
            #     # 'name': query,
            #     'seq_len': int(regx.group(2)),
            #     'program': self._meta.get('program'),
            #     'version': self._meta.get('version'),
            #     'target': self._meta.get('target'),
            # }

            # Get description and accession, if they exist
            description = '<unknown description>'  # placeholder
            while not self._line.startswith('Scores for '):
                self._read_forward()

                if self._line.startswith('Description:'):
                    description = self._line.strip().split(' ', 1)[1].strip()
                    # qresult_attrs['description'] = description

            # --- Parse the query hits ---
            hit_list = self._parse_hit(query, description)
            while self._line and not self._line.startswith('//'):
                # self._read_next_until(lambda s: s.startswith('Internal pipeline'))  # Pass all alignments info
                self._read_forward()

            pass  # TODO ?
            qresult = QueryResult(id=query, desc=description, hits=hit_list)
            yield qresult

            self._read_forward()
            while self._line.startswith('#'):
                self._read_forward()
            if self._line.startswith('[ok]'):
                break

    def _parse_hit(self, query, qdesc=''):
        """Parse a HMMER3 hit block, beginning with the hit table (PRIVATE)."""
        # get to the end of the hit table delimiter and read one more line
        self._read_next_until(lambda line: line.startswith('    ------- ------ -----'))
        self._read_forward()

        # Assume every hit is in inclusion threshold until
        # the inclusion threshold line is encountered
        is_included = True

        # Parse the hit table
        hit_attr_list = []
        while True:
            if not self._line:
                return []

            if self._line.startswith('  ------ inclusion'):
                is_included = False
                self._read_forward()
            # If there are no hits, then there are no hsps
            # so we forward-read until 'Internal pipeline..'
            elif self._line.startswith('   [No hits detected that satisfy '
                                       'reporting'):
                while True:
                    self._read_forward()
                    if self._line.startswith('Internal pipeline'):
                        assert len(hit_attr_list) == 0
                        return []
            elif self._line.startswith('Domain annotation for each '):
                self._read_forward()
                i = 0
                while True:
                    if self._line.startswith('>>'):
                        hit_id = self._line.split()[1]
                        assert hit_id == hit_attr_list[i]['id']
                        self._read_forward()
                        if self._line.startswith('   [No individual domains'):
                            hit_attr_list[i]['domain_hit'] = None
                        else: 
                            self._read_forward()
                            maxscore = -1.
                            while True:
                                self._read_forward()
                                if not re.match(_DOMHIT_PREFX, self._line):
                                    break
                                sl = self._line.split()
                                assert sl[1] == '!'
                                score = float(sl[2])
                                if score > maxscore:
                                    # sl[9] 1-based, sl[10] inclusive
                                    hit_attr_list[i]['domain_hit'] = (int(sl[9]) - 1, int(sl[10]))
                                    maxscore = score
                        i += 1
                    else:
                        self._read_forward()
                    if self._line.startswith('Internal pipeline'):
                        break

                return hit_attr_list

            # Entering hit results row
            # Parse the columns into a list
            row = [x for x in self._line.strip().split(maxsplit=9)]
            # If there's no description, set it to an empty string
            if len(row) < 10:
                row.append('')
                assert len(row) == 10
            hit_attrs = {
                'id': row[8],
                'query_id': query,
                'evalue': float(row[0]),
                'bitscore': float(row[1]),
                # 'bias': float(row[2]),
                # 'domain_exp_num': float(row[6]),
                # 'domain_obs_num': int(row[7]),
                # 'description': row[9],
                'protein': None,
                'nucleotide': None,
                'domain_hit': None,
                'evalue_domain': float(row[3]),
                'bitscore_domain': float(row[4]),
                'is_included': is_included,
            }
            hit_attr_list.append(hit_attrs)

            self._read_forward()
