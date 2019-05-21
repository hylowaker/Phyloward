import json
import os
import os.path
import re
import time
import uuid
import warnings
from collections import namedtuple, OrderedDict

from phyloward import global_config, __version__ as progver, __name__ as progname
from phyloward.hmmer_parser import HmmerParser
from phyloward.utils import SimpleFastaReader, _BaseExternalPipe

_file_abspath = os.path.dirname(os.path.abspath(__file__))
PATH_PROFILE_HMM_ARCH = os.path.join(_file_abspath, 'profile_arch.hmm')
PATH_PROFILE_HMM_BACT = os.path.join(_file_abspath, 'profile_bact.hmm')

_ssu_bac = 'RF00177'  # test
_ssu_arc = 'RF01959'
_lsu_bac = 'RF02541'
_lsu_arc = 'RF02540'


# ========= Wrapper classes for external programs =========

class ProdigalPipe(_BaseExternalPipe):  # Tested on Prodigal v2.6.2
    """Predict prokaryotic genes with Prodigal"""
    PROG_NAME = 'Prodigal'
    exec_path = global_config.prodigal_path
    _supports_stdin = True
    _TempOutPath = namedtuple('_TempOutPath', ['protein', 'nucleotide'])

    def __init__(self, input=None, *, codon_table='11'):
        super().__init__()
        self._input = input   # file name
        self._codon_table = str(codon_table) if codon_table is not None else '11'
        self.temp_path = self._locate_temp()

    def _locate_temp(self):
        import tempfile
        tempdir = tempfile.gettempdir()
        temp_pro = os.path.join(tempdir, str(uuid.uuid1()) + '.faa.tmp')
        temp_nuc = os.path.join(tempdir, str(uuid.uuid1()) + '.fna.tmp')
        return self._TempOutPath(temp_pro, temp_nuc)

    def _delete_temp(self):
        try:
            os.remove(self.temp_path.nucleotide)
            os.remove(self.temp_path.protein)
        except OSError:
            pass

    def _args_as_list(self):
        raise NotImplementedError

    @property
    def input(self):
        return self._input

    @property
    def codon_table(self):
        return self._codon_table

    @codon_table.setter
    def codon_table(self, value):
        self._codon_table = str(value)  # Foolproof

    @classmethod
    def get_version(cls):
        _, err = cls._execute_process(options='-v')
        regx = re.search(r'Prodigal\s+V(\d+\.\d+\.?\d*).+', err.strip())
        version = regx.group(1) if regx else None
        return version

    def _options(self):
        options = ['-a', self.temp_path.protein,
                   '-d', self.temp_path.nucleotide,
                   '-g', self._codon_table,
                   '-q']
        if self._input is not None:  # else: read from stdin
            options += ['-i', self._input]
        return options


class HmmsearchPipe(_BaseExternalPipe):
    """Select core genes from predicted genes"""
    PROG_NAME = 'hmmsearch'
    exec_path = global_config.hmmsearch_path
    _supports_stdin = False

    def __init__(self, hmmfile, seqdb):
        super().__init__()
        self.profile = hmmfile
        self.seqdb = seqdb
        self.threads = 4
        self._cut_arg = '--cut_tc'

    def get_version(self):
        out, _ = self._execute_process(options='-h')
        regx = re.search(r'^# HMMER (\d+\.\d+\w*)', out.splitlines()[1])
        version = regx.group(1) if regx else None
        return version

    def set_cut_threshold(self, cutoff):
        if cutoff.casefold in ('trusted', 'tc'):
            self._cut_arg = '--cut_tc'
        elif cutoff.casefold in ('gathering', 'ga'):
            self._cut_arg = '--cut_ga'
        elif cutoff.casefold in ('noise', 'nc'):
            self._cut_arg = '--cut_nc'
        elif cutoff is None or cutoff.casefold == 'none':
            self._cut_arg = ''
        else:
            raise ValueError('cutoff parameter must be one of {\'tc\', \'ga\', \'nc\', None}')

    def _options(self):
        options = [self._cut_arg,
                   '--noali',
                   '--cpu', str(self.threads),
                   self.profile,
                   self.seqdb]
        return options


# ========= Data class for extracted genes =========


class _Sequence:
    __slots__ = ('_protein', '_nucleotide', '_domain_start', '_domain_end')

    def __init__(self, prot, nucl):
        self._protein = prot if prot else ''  # type: str
        self._nucleotide = nucl if nucl else ''  # type: str
        self._domain_start = None
        self._domain_end = None
    
    # @classmethod  #
    # def from_extracted_core_gene(cls, extracted, gene):
    #     if not isinstance(extracted, ExtractedCoreGenes):
    #         raise TypeError
    #     aa_seq = extracted.best_hit_sequence(gene)  # None-able
    #     nt_seq = extracted.best_hit_sequence(gene, nucleotide=True)
    #     seq = cls(aa_seq, nt_seq)
    #     return seq
    
    def set_domain_region(self, start, end):
        if not isinstance(start, int) or not isinstance(end, int):
            raise TypeError
        self._domain_start = start  # inclusive
        self._domain_end = end  # exclusive
    
    @property
    def protein(self):
        return self._protein
    
    @property
    def protein_domain(self):
        a = self._domain_start
        b = self._domain_end
        if a is not None and b is not None:
            return self._protein[a:b]
        return ''

    @property
    def nucleotide(self):
        return self._nucleotide
    
    @property
    def nucleotide_domain(self):
        a = self._domain_start
        b = self._domain_end
        if a is not None and b is not None:
            return self._nucleotide[a*3:b*3]
        return ''

    def __bool__(self):
        return bool(self._nucleotide)

    def __str__(self):
        if self.protein:
            return self._protein
        else:
            return self._nucleotide


class ExtractedCoreGenes:
    def __init__(self, info, data):
        self.info = info    # Meta-information (dictionary)
        self.data = data    # Actual core gene data (dictionary)

    def __bool__(self):
        pipeline_name = self.info and self.info.get('program')
        if not pipeline_name or not pipeline_name.startswith('phyloward'):
            return False
        return bool(self.info and self.data)

    @property
    def is_archaeal(self):
        domain = self.info['domain']
        if domain.casefold() == 'archaea':
            return True
        elif domain.casefold() == 'bacteria':
            return False
        raise ValueError('The domain(superkingdom) is neither Archaea nor Bacteria!')

    @classmethod
    def from_dict(cls, d):
        return cls(d.get('info'), d.get('data'))

    def get_label(self):
        return self.info['label']

    def get_uid(self):
        return self.info.get('uid')

    def get_hit_list(self, gene):
        return self.data[gene]

    def gene_set(self, *, include_empty=False):
        if include_empty:
            return set(self.data.keys())
        else:
            return set(gene for gene in self.data.keys() if self.get_hit_list(gene))

    def best_hit_sequence(self, gene):
        """Return hit Sequence object which has highest e-value"""
        hits = self.get_hit_list(gene)  # assumes already sorted by evalue!
        for hit in hits:
            if hit.get('is_included'):
                seq = _Sequence(hit['protein'], hit['nucleotide'])
                try: 
                    hit_start, hit_end = hit['domain_hit'].split(':')
                    seq.set_domain_region(int(hit_start), int(hit_end))
                except AttributeError:
                    pass
                return seq
        return _Sequence('', '')

    def as_dict(self):
        return OrderedDict(info=self.info, data=self.data)

    def as_json(self, indent=None):
        import json
        return json.dumps(self.as_dict(), indent=indent)

    @classmethod
    def from_json(cls, s):
        import json
        d = json.loads(s, object_pairs_hook=OrderedDict)
        try:  # check extracted gene json format validity
            if not d['info']['program'].startswith('phyloward'):
                return None
        except (KeyError, AttributeError):
            return None

        return cls(info=d['info'], data=d['data'])



# ========= Module functions =========


def extract_core_genes(file, *, is_archaea=False, profile=None, **kwargs):
    """Search core genes from a genome file

    :param file: A FASTA file with genome sequences.
    :param archaea: If True, extract archaeal core genes instead of bacterial ones
    :param kwargs: Metadata (ex. label)
    :return: ExtractedCoreGenes object
    """
    try:   # file is IO-like
        s = file.read()
        prodigal = ProdigalPipe()
        prodigal.stdin = s
    except AttributeError:   # file is str (file name)
        prodigal = ProdigalPipe(file)
    if kwargs.get('codon_table'):
        prodigal.codon_table = kwargs['codon_table']

    # ------ Prodigal ------
    try:
        prodigal.run()
    except RuntimeError:
        raise
    tmpaa, tmpnt = prodigal.temp_path

    # ------ hmmsearch ------
    if profile is None:
        profile = PATH_PROFILE_HMM_ARCH if is_archaea else PATH_PROFILE_HMM_BACT
    profile = os.path.abspath(profile)

    hmmsearch = HmmsearchPipe(profile, tmpaa)
    try:
        rawver = hmmsearch.get_version()
        intver = int(rawver.split('.')[0])
        if intver != 3:
            warnings.warn('ERROR! This version({}) of HMMER package is not supported.'
                          ' Please use version 3.0 or later.'.format(rawver))
        hmmsearch.run()
    except (FileNotFoundError, PermissionError):
        raise
    except RuntimeError:
        raise

    # --- Parse hmmsearch result ---
    from io import StringIO
    with StringIO(hmmsearch.stdout) as f:
        hmm_result = _parse_hmmsearch_output(f)  # TODO as non-blocking stream

    # ------ Take protein & nucleotide sequence from Prodigal ------
    with open(tmpaa) as fh_aa, \
            open(tmpnt) as fh_nt:
        aa_map = {label: seq[:-1] if seq.endswith('*') else seq
                  for label, seq in SimpleFastaReader(fh_aa)}  # trim '*' char in fasta
        nt_map = {label: seq for label, seq in SimpleFastaReader(fh_nt)}

    try:  # Removing temporary files by prodigal
        os.remove(tmpaa)
        os.remove(tmpnt)
    except OSError:
        pass

    # --- Merge alternative model hits ---  # Caution! manipulating mutable objects!
    marked_to_delete = []
    for gene, hits in hmm_result.items():
        if gene.endswith('*'):  # hit is from alternative model
            marked_to_delete.append(gene)
            original = gene[:-1]
            # for hit in hits:
            #     hit['query_id'] = original  # reset query_id with original gene name
            hmm_result[original].extend(hits)
            hmm_result[original].sort(key=lambda x: x['evalue'])

    for each in marked_to_delete:
        del hmm_result[each]

    # ------ Add protein & nucleotide sequence to hmmresult ------  
    n_cores, n_hits = 0, 0
    for _gene, hits in hmm_result.items():  # Caution! manipulating mutable objects!
        if hits:
            n_cores += 1
        check_duplicate = set()
        for hit in hits:
            hit_id = hit['id']
            s_aa, s_nt = aa_map[hit_id], nt_map[hit_id]
            hit['protein'] = s_aa
            hit['nucleotide'] = s_nt
            if (s_aa, s_nt) in check_duplicate:  # if identical seq exists, use only one.
                hit['is_included'] = False
            check_duplicate.add((s_aa, s_nt))
            if hit.get('is_included'):
                n_hits += 1

    # ------ Metadata ------
    meta = OrderedDict()
    try:
        meta.update(kwargs if kwargs else {})
    except TypeError as e:
        warnings.warn(e)

    if not meta.get('label'):
        # file name as default label
        meta['label'] = os.path.splitext(os.path.basename(file))[0]
    meta['label'] = re.sub(r'\s+', '_', meta['label'])  # replace whitespaces
    meta['uid'] = uuid.uuid1().int >> 64
    meta['domain'] = 'Bacteria' if not is_archaea else 'Archaea'
    meta['hmm_profile'] = profile
    meta['program'] = progname + ' ' + progver
    meta['created_on'] = time.strftime('%b %d %Y %H:%M:%S')
    meta['extracted_from'] = file
    meta['prodigal_ver'] = prodigal.get_version()
    meta['hmmer_ver'] = hmmsearch.get_version()
    meta['n_core_genes'] = n_cores
    meta['n_hits'] = n_hits

    extracted = ExtractedCoreGenes(info=meta, data=hmm_result)
    return extracted


def _parse_hmmsearch_output(stream):
    hmmparser = HmmerParser(stream)
    result = OrderedDict((qresult.id, qresult.hits) for qresult in hmmparser)
    # id (core gene) => hits
    return result


def _is_valid(d):
    try:  # check extracted gene json format validity
        if not d['info']['program'].startswith('phyloward'):
            return False
        if d['info']['domain'].casefold() not in ('bacteria', 'archaea'):
            return False
    except (KeyError, AttributeError):
        return False
    return True


def extracted_list_from_json(*files, raw_string=None):
    """Return a list of ExtractedCoreGenes objects from JSON

    :param files: JSON files to read from
    :param raw_string: Raw JSON string. If defined, 'files' parameter will be ignored.
    :return: List of ExtractedCoreGenes objects
    """
    if raw_string:
        d = json.loads(raw_string, object_pairs_hook=OrderedDict)
        if isinstance(d, dict):
            return [d, ]
        return [ExtractedCoreGenes.from_dict(each) for each in d if _is_valid(each)]

    def _iter_read(files):
        for f in files:
            with open(f) as fh:
                yield fh.read()

    extracted_list = []
    for s in _iter_read(files):
        d = json.loads(s, object_pairs_hook=OrderedDict)
        if isinstance(d, dict) and _is_valid(d):
            extracted_list.append(ExtractedCoreGenes.from_dict(d))
        elif isinstance(d, list):
            extracted_list.extend([ExtractedCoreGenes.from_dict(x) for x in d if _is_valid(x)])
    extracted_list = [ex for ex in extracted_list if ex]
    return extracted_list
