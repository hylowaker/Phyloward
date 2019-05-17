import concurrent.futures.thread
import io
import json
import os
import re
import tempfile
import uuid
import warnings
from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor, as_completed

from phyloward import extractor
from phyloward import global_config
from phyloward.branchanalysis import BranchAnalysis
from phyloward.utils import SimpleFastaReader, SimpleFastaWriter, _BaseExternalPipe


# TODO zZ and label replacer!


# ========= Wrapper classes for external programs =========

# TODO different aligner

class MafftPipe(_BaseExternalPipe):  # Tested on MAFFT v7.313
    """Build multiple alignments of each core gene"""
    PROG_NAME = 'MAFFT'
    exec_path = global_config.mafft_path
    supports_stdin = False

    def __init__(self, input=''):
        super().__init__()
        self._input = input  # fasta file  (doesn't support stdin)
        self.threads = 1

    @property
    def input(self):
        return self._input

    @classmethod
    def get_version(cls):
        _, err = cls._execute_process(options='--version')  # type: str
        regx = re.search(r'^v(\d+\.\d+)', err)
        version = regx.group(1) if regx else None
        return version

    def _options(self):
        options = ['--auto', '--thread', str(self.threads),
                   self._input]
        # options = ['--globalpair', '--maxiterate', '1000',
        #            '--thread', str(self.threads),
        #            self._input]  # test
        return options


class _TreePipeBase(_BaseExternalPipe):
    def __init__(self, input=None):
        super().__init__()
        self._input = input
        # self._seq_type = 'protein'
        self.is_nucl = False  # TODO property
        self.model = None

    # @property
    # def seq_type(self):
    #     return self._seq_type
    #
    # @seq_type.setter
    # def seq_type(self, value):
    #     if value in ('protein', 'prot', 'pro', 'aa'):
    #         self._seq_type = 'protein'
    #     elif value in ('nucleotide', 'nucl', 'nuc', 'nt'):
    #         self._seq_type = 'nucleotide'
    #     else:
    #         raise ValueError

    @classmethod
    def get_version(cls):
        raise NotImplementedError('Subclass should implement this')

    def _options(self, *args, **kwargs):
        raise NotImplementedError('Subclass should implement this')

    @property
    def input(self):
        return self._input

    def produce_tree(self):
        raise NotImplementedError('Subclass should implement this')


class FastTreePipe(_TreePipeBase):  # Tested on FastTree 2.1.10
    PROG_NAME = 'FastTree'
    exec_path = global_config.fasttree_path
    supports_stdin = True

    def __init__(self, input=None):
        super().__init__(input)

    @classmethod
    def get_version(cls):
        _, err = cls._execute_process(options='-expert')
        regx = re.search(r'FastTree\s?(?:version)?\s?(\d+\.+\d+\.?\d*)', err.splitlines()[0])
        version = regx.group(1) if regx else None
        return version

    def _options(self):
        options = []
        # if self._seq_type == 'nucleotide':
        if self.is_nucl:
            options.append('-nt')
            options.append('-gtr')
        if self._input:
            options.append(self._input)  # test
        return options

    def produce_tree(self):
        if self._stdout is not None:
            s = self._stdout.strip()  # type: str
            if s.startswith('(') and s.endswith(';'):
                return s
            else:
                raise ValueError('Something has gone wrong with newick tree')
        else:
            self.run()
            return self.produce_tree()


class RAxMLPipe(_TreePipeBase):  # Tested on RAxML 8.2.12
    PROG_NAME = 'RAxML'
    exec_path = global_config.raxml_path
    supports_stdin = False

    def __init__(self, input):
        super().__init__(input)
        self.run_id = os.path.basename(self._input)

    @classmethod
    def get_version(cls):
        out, _ = cls._execute_process(options='-h')
        regx = re.search(r'RAxML\s*(?:version)?\s*(\d+\.\d+\.?\d*)', out)
        version = regx.group(1) if regx else None
        return version

    def _options(self):
        options = []
        options += ['-s', self._input]
        options += ['-n', self.run_id]
        options += ['-p', '123']
        # assert self._seq_type in ('protein', 'nucleotide')
        # if self._seq_type == 'protein':
        if self.is_nucl:
            options += ['-m', 'GTRCAT']
        else:
            options += ['-m', 'PROTCATJTT']
        return options

    def produce_tree(self):
        # I hate RAxML's interface!
        tempdir = tempfile.gettempdir()
        self.run(cwd=tempdir)
        tempf = os.path.join(tempdir, 'RAxML_bestTree.{}'.format(self.run_id))
        with open(tempf) as fh:
            result = fh.read().strip()
        assert result.startswith('(') and result.endswith(';')
        try:
            _j = os.path.join
            os.remove(_j(tempdir, 'RAxML_bestTree.{}'.format(self.run_id)))
            os.remove(_j(tempdir, 'RAxML_info.{}'.format(self.run_id)))
            os.remove(_j(tempdir, 'RAxML_log.{}'.format(self.run_id)))
            os.remove(_j(tempdir, 'RAxML_parsimonyTree.{}'.format(self.run_id)))
            os.remove(_j(tempdir, 'RAxML_result.{}'.format(self.run_id)))
        except OSError:
            pass
        try:
            os.remove(self._input + '.reduced')
        except OSError:
            pass

        return result


# ========= Data class for aligned genes =========


CAT = 'Concatenated'
msa_engine = MafftPipe
tree_engine = FastTreePipe
# tree_engine = RAxMLPipe


class NotEnoughGenesExtracted(ValueError):
    pass


class NotAlignedError(ValueError):
    pass


class _SequenceAligned(extractor._Sequence):  # TODO Need refactor: Do not inherit
    __slots__ = ('is_aligned')

    def __init__(self, prot, nucl, *, aligned=False):
        super().__init__(prot, nucl)
        self.is_aligned = aligned

    def get_aligned_codon(self, *, omit_third=False):
        if not self.protein:
            raise ValueError('No protein seqeunce exists')
        if not self.is_aligned:
            raise NotAlignedError

        if omit_third:
            return ''.join(ch for i, ch in enumerate(self.nucleotide) if i % 3 != 2)
        else:
            return self.nucleotide


class _SequenceBunch:
    def __init__(self):
        self._data = OrderedDict()
        self._is_aligned = False
        self._is_aligned_by_domain = False  # TODO
        self._is_nucl_alignment = False
        self.is_filtered = False
        self._fasta_seq = None
        self._fasta_seq_width = None
        self.newick_tree = None

    def __getitem__(self, item):
        return self._data[item]

    def __setitem__(self, key, value):
        self._data[key] = value

    # def __getattr__(self, item):
    #     return self._data.__getattribute__(item)

    def __len__(self):
        return len(self._data)

    def __iter__(self):
        self._data.__iter__()

    @property
    def is_aligned(self):
        return self._is_aligned

    @is_aligned.setter
    def is_aligned(self, b):
        self._is_aligned = b
        for _label, seq in self.items():
            assert isinstance(seq, _SequenceAligned)
            seq.is_aligned = b

    def keys(self):
        return self._data.keys()

    def values(self):
        return self._data.values()

    def items(self):
        return self._data.items()

    def align(self, domain_only=True, nucleotide=False):  # TODO Refactor domain alignment. it's spaghetti code for now :(
        tempdir = tempfile.gettempdir()
        tempf_align = os.path.join(tempdir, str(uuid.uuid1()) + '.fasta.tmp')  # TODO reduce disk io?

        if domain_only:
            d = dict((label, seq.nucleotide_domain) if nucleotide else (label, seq.protein_domain)
                     for label, seq in self.items())
        else:
            d = dict((label, seq.nucleotide) if nucleotide else (label, seq.protein)
                     for label, seq in self.items())

        if sum(bool(s) for s in d.values()) < 3:
            raise NotEnoughGenesExtracted('Extracted genes are less than 3')
        elif not any(d.values()):
            raise NotEnoughGenesExtracted

        # --- Aligning process ---
        with open(tempf_align, 'w') as fh:  # writes temporary file
            print(SimpleFastaWriter(d), file=fh, end='')
        try:
            msa_program = msa_engine(input=tempf_align)
            msa_program.run()
        except (FileNotFoundError, PermissionError) as e:
            raise e
        except RuntimeError as e:
            raise e

        with io.StringIO(msa_program.stdout) as fp:
            aligned_map = {label: s for label, s in SimpleFastaReader(fp)}

        # --- Assign aligned Sequence values ---
        for label, s in aligned_map.items():
            if not nucleotide:  # also fill codon alignments
                nt = list(self[label].nucleotide_domain if domain_only else self[label].nucleotide)
                for i, ch in enumerate(s):
                    if ch == '-':
                        nt[i * 3:i * 3] = list('---')
                last_three = ''.join(nt[-3:]).casefold()
                if last_three not in ('taa', 'tag', 'tga'):
                    nt += list('---')
                assert len(nt) == 3 * len(s) + 3
                self[label] = _SequenceAligned(s, ''.join(nt), aligned=True)
            else:
                self[label] = _SequenceAligned('', s, aligned=True)  # TODO fix info loss!
        self._is_aligned = True
        self._is_aligned_by_domain = True if domain_only else False
        self._is_nucl_alignment = True if nucleotide else False

        try:  # removes temporary file
            os.remove(tempf_align)
        except OSError:
            pass

    def _get_alignment_map(self, seq_type):
        if not self._is_aligned:
            raise NotAlignedError

        result = OrderedDict()
        for label, seq in self.items():
            if not seq.is_aligned:
                raise NotAlignedError

            if seq_type in ('amino', 'protein', 'prot', 'pro', 'aa'):
                result[label] = seq.protein
            elif seq_type in ('nucleotide', 'nucl', 'nuc', 'nt'):
                result[label] = seq.nucleotide
            elif seq_type == 'codon':
                result[label] = seq.get_aligned_codon()
            elif seq_type == 'codon12':
                result[label] = seq.get_aligned_codon(omit_third=True)
            else:
                raise ValueError('Invalid sequence_type: {}'.format(seq_type))
        return result

    def as_fasta(self, seq_type=None, width=None):
        if seq_type is None:
            seq_type = 'protein' if not self._is_nucl_alignment else 'nucleotide'
        if self._is_nucl_alignment \
                and seq_type in ('protein', 'prot', 'pro', 'aa', 'codon', 'codon12'):
            raise ValueError('Cannot get protein alignment since the data is nucleotide alignment.')

        if self._fasta_seq is None or width != self._fasta_seq_width:
            self._fasta_seq = str(SimpleFastaWriter(self._get_alignment_map(seq_type), width))
        return self._fasta_seq

    def length_of_alignment(self, seq_type=None):
        if not self._is_aligned:
            raise NotAlignedError

        it = iter(SimpleFastaReader(io.StringIO(self.as_fasta(seq_type))))
        _, seq0 = next(it)
        _, seq1 = next(it)
        assert len(seq0) == len(seq1)
        return len(seq0)

    def build_tree(self, seq_type=None, args=None):
        if not self._is_aligned:
            raise NotAlignedError

        if tree_engine.supports_stdin:  # FastTree
            tree_builder = tree_engine()
            tree_builder.stdin = self.as_fasta(seq_type)
        else:  # stupid RAxML
            tempdir = tempfile.gettempdir()
            tempf = os.path.join(tempdir, str(uuid.uuid1().int) + '.temp.fasta')
            with open(tempf, 'w') as fh:
                print(self.as_fasta(seq_type), file=fh, end='')
            tree_builder = tree_engine(input=tempf)

        if self._is_nucl_alignment:
            tree_builder.seq_type = 'nucleotide'
        else:
            tree_builder.seq_type = 'protein'

        nwk_tree = tree_builder.produce_tree()
        assert nwk_tree
        self.newick_tree = nwk_tree  # type: str
        return nwk_tree

    # def get_newick_tree(self, seq_type=None):
    #     if not self._nwk_tree:
    #         self._build_tree(seq_type)
    #     return self._nwk_tree


class AlignedCoreGenes():
    def __init__(self):
        self._data = OrderedDict()
        self._is_archaeal = False
        self._core_genes = None
        self.is_aligned = False
        self._is_domain_only = False
        self._is_nucleotide_alignment = False
        self._excluded = []
        try:
            self.msa_program_info = (msa_engine.PROG_NAME, msa_engine.get_version())
        except (FileNotFoundError, RuntimeError):
            self.msa_program_info = None

        self._genomes = set()
        self._comments = {}
        self.uid_label = {}

    def __getitem__(self, key):
        return self._data[key]
    
    def __setitem__(self, key, value):
        self._data[key] = value
    
    def __iter__(self):
        raise NotImplementedError  # TODO implement this

    @property
    def is_archaeal(self):
        return self._is_archaeal

    def genes(self, scope='included'):
        if scope == 'all':
            cores = self._core_genes[:]
        elif scope == 'included':
            cores = [gene for gene in self._core_genes if gene not in self._excluded]
        elif scope == 'excluded':
            cores = self._excluded[:]
        else:
            raise ValueError('Parameter \'scope\' must be one of '
                             '(\'all\', \'included\', \'excluded\')')
        return cores

    @classmethod
    def from_extracted_list(cls, extracted_list):
        cg_aligned = cls()

        # --- Infer phylogenetic domain (superkingdom) ---
        is_all_archaea = all(each.is_archaeal for each in extracted_list)
        is_all_bacteria = all(not each.is_archaeal for each in extracted_list)
        if is_all_archaea:
            cg_aligned._is_archaeal = True
        elif is_all_bacteria:
            cg_aligned._is_archaeal = False
        else:
            raise ValueError('All genomes must belong to same domain(superkingdom).')

        # --- are profiles identical? ---
        if len(set(each.info['hmm_profile'] for each in extracted_list)) != 1:
            raise ValueError('HMM profile not identical!')
        cg_aligned._comments['hmm_profile'] = extracted_list[0].info['hmm_profile']

        # -- Mode: Full sequence or Best 1 domain
        # if mode == 'full':
        #     best_domain = False
        # elif 'domain' in mode:
        #     best_domain = True
        # else:
        #     raise ValueError('Mode must be either \'full\' or \'domain\'')

        # --- Set un-aligned Sequence values ---
        core_genes = extracted_list[0].gene_set(include_empty=True)  # type: set
        for gene in sorted(core_genes):
            seqdict = _SequenceBunch()
            for each in extracted_list:
                if not each:
                    continue
                assert isinstance(each, extractor.ExtractedCoreGenes)
                if each.gene_set(include_empty=True) != core_genes:
                    raise ValueError('Extracted core gene information may be corrupted.')
                seq = each.best_hit_sequence(gene)
                # label_genome = re.sub(r'\s+', '_', each.get_label())  # replace whitespaces
                label_genome = each.get_label()
                seqdict[label_genome] = seq
                cg_aligned._genomes.add(label_genome)
                cg_aligned.uid_label[each.get_uid()] = each.get_label()
            cg_aligned[gene] = seqdict

        cg_aligned._core_genes = sorted(core_genes)
        assert cg_aligned._data.keys() == core_genes
        return cg_aligned

    # def _determine_core_genes(self):
    #     pass

    def _concatenate_aligned(self):
        assert self.is_aligned

        seqdict_to_concat = OrderedDict()  # genome name => list of aliged sequences to be joined
        seqdict_joined = _SequenceBunch()
        for gene, seqdict in self._data.items():
            if gene in self._excluded:
                continue
            assert self._genomes == seqdict.keys()
            for label_genome, seq in seqdict.items():
                # Append (aligned) sequence of each gene
                assert seq.is_aligned
                seqdict_to_concat.setdefault(label_genome, []).append(seq)

            for label, seqlist in seqdict_to_concat.items():
                joined_aaseq = ''.join(seq.protein for seq in seqlist)
                joined_ntseq = ''.join(seq.nucleotide for seq in seqlist)
                joined_seq = _SequenceAligned(joined_aaseq, joined_ntseq, aligned=True)
                seqdict_joined[label] = joined_seq
        seqdict_joined.is_aligned = True
        self[CAT] = seqdict_joined
        self._data.move_to_end(CAT, last=False)

    def _filter_concatenated(self, cutoff=0.):
        assert self.is_aligned
        assert self._data.get(CAT)

        if not cutoff:
            return
        elif not 0. <= cutoff <= 1.:
            raise ValueError('cutoff value out of range [0.0, 1.0]')

        def _c(seq):
            assert seq.is_filtered is False
            return seq

        seqdict = self[CAT]  # type: _SequenceBunch
        seqmatrix_aa = [list(_c(seq).protein) for seq in seqdict.values()]
        seqmatrix_nt = [list(_c(seq).nucleotide) for seq in seqdict.values()]

        col_length = len(seqdict)
        threshold = col_length * (1. - cutoff)

        i = 0
        excluded_aa = []
        while True:
            try:
                column_aa = [s[i] for s in seqmatrix_aa]
            except IndexError:
                break
            if cutoff == 1. and '-' in column_aa:
                excluded_aa.append(i)
            elif column_aa.count('-') > threshold:
                excluded_aa.append(i)
            i += 1

        i = 0
        excluded_nt = []
        while True:
            try:
                column_nt = [s[i] for s in seqmatrix_nt]
            except IndexError:
                break
            if cutoff == 1. and '-' in column_nt:
                excluded_nt.append(i)
            elif column_nt.count('-') > threshold:
                excluded_nt.append(i)
            i += 1

        # assign new filtered string to each Sequence object
        for seq in seqdict.values():
            seq.protein = ''.join([ch for i, ch in enumerate(seq.protein) if i not in excluded_aa])
            seq.nucleotide = ''.join([ch for i, ch in enumerate(seq.nucleotide) if i not in excluded_nt])
            seq.is_filtered = True

    def _align_each(self, gene):
        seqdict = self[gene]  # type: _SequenceBunch
        try:
            seqdict.align(domain_only=self._is_domain_only, nucleotide=self._is_nucleotide_alignment)
        except NotEnoughGenesExtracted:
            self._excluded.append(gene)

    def align(self, *, domain_only=True, nucleotide=False, cutoff=0., workers=1, report_progress=None):
        """Dummy Docstring"""
        self._is_domain_only = domain_only
        self._is_nucleotide_alignment = nucleotide

        progress_reporter = None
        if report_progress:
            try:
                progress_reporter = report_progress(len(self._data), head_message='Align : ')
            except TypeError:
                progress_reporter = report_progress(len(self._data))
            if not hasattr(progress_reporter, 'send'):
                raise AttributeError

        with ThreadPoolExecutor(workers) as pool:
            futures = [pool.submit(self._align_each, gene) for gene in self._core_genes]
            try:
                for future in as_completed(futures):
                    future.result()
                    try:
                        progress_reporter.send('')
                    except StopIteration:
                        break
            except (Exception, KeyboardInterrupt):
                pool._threads.clear()
                concurrent.futures.thread._threads_queues.clear()
                raise

        self.is_aligned = True
        self._concatenate_aligned()
        self._filter_concatenated(cutoff=cutoff)

    def number_of_genomes(self):
        return len(self._genomes)

    def add_comment(self, field, value):
        self._comments[field] = value

    def get_fasta(self, gene=None, seq_type=None):
        if gene in self._excluded:
            warnings.warn('The gene {} is excluded from analysis.'.format(gene))
            return None

        seqdict = self[gene if gene else CAT]  # type: _SequenceBunch
        return seqdict.as_fasta(seq_type)

    def as_dict_with_fasta(self, seq_type=None):
        gene_to_fasta_map = OrderedDict()
        for gene in self._data.keys():  # if aligned, includes CAT
            gene_to_fasta_map[gene] = None
            if gene in self._excluded:
                continue
            gene_to_fasta_map[gene] = self.get_fasta(gene, seq_type)
        return gene_to_fasta_map

    # noinspection PyTypeChecker
    def result_as_json(self, seq_type=None, *, indent=None):
        data = OrderedDict()
        data['__comment__'] = self._comments
        for gene, seqdict in self._data.items():
            data[gene] = OrderedDict()
            if gene in self._excluded:
                data[gene]['__excluded__'] = 1
                continue
            data[gene]['fasta'] = seqdict.as_fasta(seq_type)
            # data[gene]['tree'] = seqdict.newick_tree
        return json.dumps(data, indent=indent)


class CoreGenesNewickTree(AlignedCoreGenes):
    def __init__(self):
        super().__init__()
        self._tree_gsi_marked = None
        try:
            self.tree_program_info = (tree_engine.PROG_NAME, tree_engine.get_version())
        except (FileNotFoundError, RuntimeError):
            self.tree_program_info = None

    def _draw_each_tree(self, gene=None, seq_type=None, args=None):
        if gene in self._excluded:
            return
        seqdict = self[gene if gene else CAT]  # type: _SequenceBunch
        seqdict.build_tree(seq_type)

    def fill_trees(self, seq_type=None, *, args=None, workers=1, report_progress=None):
        """TODO docstring"""
        progress_reporter = None
        if report_progress:
            try:
                progress_reporter = report_progress(len(self._data), head_message='Tree  : ')
            except TypeError:
                progress_reporter = report_progress(len(self._data))
            if not hasattr(progress_reporter, 'send'):
                raise AttributeError

        with ThreadPoolExecutor(workers) as pool:
            futures = [pool.submit(self._draw_each_tree, gene, seq_type, args) for gene in self._data]
            try:
                for future in as_completed(futures):
                    future.result()
                    try:
                        progress_reporter.send('')
                    except (AttributeError, StopIteration):
                        pass
            except (Exception, KeyboardInterrupt):
                pool._threads.clear()
                concurrent.futures.thread._threads_queues.clear()
                raise

    def calculate_gene_support_index(self, threshold=.95):
        n_genome = self.number_of_genomes()
        all_gene_trees = [self.get_tree(gene) for gene in self._data
                          if gene != CAT and gene not in self._excluded]
        assert all(isinstance(tree, str) for tree in all_gene_trees)
        ba = BranchAnalysis(all_gene_trees)
        marked = ba.mark_tree(self.get_tree(),
                              is_rooted=False, is_core=True,
                              missing_tolerance=-1, mismatch_tolerance=int(n_genome*(1.-threshold)))
        self._tree_gsi_marked = marked
        return marked

    def get_tree(self, gene=None):
        if gene in self._excluded:
            warnings.warn('The gene {} is excluded from analysis.'.format(gene))
            return None

        seqdict = self[gene if gene else CAT]  # type: _SequenceBunch
        return seqdict.newick_tree  # TODO check exists

    def get_tree_gsi(self):
        return self._tree_gsi_marked  # TODO check exists

    # noinspection PyTypeChecker
    def result_as_json(self, seq_type=None, *, indent=None):
        data = OrderedDict()
        data['__comment__'] = self._comments
        for gene, seqdict in self._data.items():
            data[gene] = OrderedDict()
            if gene in self._excluded:
                data[gene]['__excluded__'] = 1
                continue
            data[gene]['fasta'] = seqdict.as_fasta(seq_type)
            data[gene]['tree'] = seqdict.newick_tree
        data[CAT]['tree_gsi_marked'] = self._tree_gsi_marked
        return json.dumps(data, indent=indent)


# ========= Module functions =========


def import_from_files(file_list):
    extracted_list = extractor.extracted_list_from_json(*file_list)
    return CoreGenesNewickTree.from_extracted_list(extracted_list)


def import_from_extracted_list(extracted_list):
    return CoreGenesNewickTree.from_extracted_list(extracted_list)
