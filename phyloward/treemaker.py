# import concurrent.futures.thread
# import io
# import os
# import re
# import tempfile
# import uuid
# from collections import OrderedDict
# from concurrent.futures import ThreadPoolExecutor, as_completed
#
# from phyloward import global_config
# from phyloward.aligner import CAT, AlignedCoreGenes, _SequenceBunch
# from phyloward.branchanalysis import BranchAnalysis
# from phyloward.utils import _BaseExternalPipe, SimpleFastaWriter, SimpleFastaReader
#
#
# # =========
#
# class _TreePipeBase(_BaseExternalPipe):
#     def __init__(self, input=None):
#         super().__init__()
#         self._input = input
#         self._seq_type = 'protein'
#         self.model = None
#
#     @property
#     def seq_type(self):
#         return self._seq_type
#
#     @seq_type.setter
#     def seq_type(self, value):
#         if value in ('protein', 'prot', 'pro', 'aa'):
#             self._seq_type = 'protein'
#         elif value in ('nucleotide', 'nucl', 'nuc', 'nt'):
#             self._seq_type = 'nucleotide'
#         else:
#             raise ValueError
#
#     def get_version(self):
#         raise NotImplementedError('Subclass should implement this')
#
#     def _options(self, *args, **kwargs):
#         raise NotImplementedError('Subclass should implement this')
#
#     @property
#     def input(self):
#         return self._input
#
#     def produce_tree(self):
#         raise NotImplementedError('Subclass should implement this')
#
#
# class FastTreePipe(_TreePipeBase):  # Tested on FastTree 2.1.10
#     PROG_NAME = 'FastTree'
#     exec_path = global_config.fasttree_path
#     _supports_stdin = True
#
#     def __init__(self, input=None):
#         super().__init__(input)
#
#     def get_version(self):
#         _, err = self._execute_process(options='-expert')
#         regx = re.search(r'FastTree\s?(?:version)?\s?(\d+\.+\d+\.?\d*)', err.splitlines()[0])
#         version = regx.group(1) if regx else None
#         return version
#
#     def _options(self):
#         options = []
#         if self._seq_type == 'nucleotide':
#             options.append('-nt')
#             options.append('-gtr')
#         # TODO other options
#         if self._input:
#             options.append(self._input)  # test
#         return options
#
#     def produce_tree(self):
#         if self._stdout is not None:
#             s = self._stdout.strip()  # type: str
#             if s.startswith('(') and s.endswith(';'):
#                 return s
#             else:
#                 return None
#         else:
#             self.run()
#             return self.produce_tree()
#
#
# class RAxMLPipe(_TreePipeBase):  # Tested on RAxML 8.2.12
#     PROG_NAME = 'RAxML'
#     exec_path = global_config.raxml_path
#     _supports_stdin = False
#
#     def __init__(self, input):
#         super().__init__(input)
#         self.run_id = os.path.basename(self._input)
#
#     def get_version(self):
#         out, _ = self._execute_process(options='-h')
#         regx = re.search(r'RAxML\s*(?:version)?\s*(\d+\.\d+\.?\d*)', out)
#         version = regx.group(1) if regx else None
#         return version
#
#     def _options(self):
#         options = []
#         options += ['-s', self._input]
#         options += ['-n', self.run_id]
#         options += ['-p', '123']
#         assert self._seq_type in ('protein', 'nucleotide')
#         if self._seq_type == 'protein':
#             options += ['-m', 'PROTCATJTT']
#         elif self._seq_type == 'nucleotide':
#             options += ['-m', 'GTRCAT']
#         # TODO custom options
#         return options
#
#     def produce_tree(self):
#         # I hate RAxML's interface!
#         tempdir = tempfile.gettempdir()
#         self.run(cwd=tempdir)
#         with open('RAxML_bestTree.{}'.format(self.run_id)) as fh:
#             result = fh.read().strip()
#         assert result.startswith('(') and result.endswith(';')
#         try:
#             _j = os.path.join
#             os.remove(_j(tempdir, 'RAxML_bestTree.{}'.format(self.run_id)))
#             os.remove(_j(tempdir, 'RAxML_info.{}'.format(self.run_id)))
#             os.remove(_j(tempdir, 'RAxML_log.{}'.format(self.run_id)))
#             os.remove(_j(tempdir, 'RAxML_parsimonyTree.{}'.format(self.run_id)))
#             os.remove(_j(tempdir, 'RAxML_result.{}'.format(self.run_id)))
#         except OSError:
#             pass
#         try:
#             os.remove(self._input + '.reduced')
#         except OSError:
#             pass
#
#         return result
#
#
# # =========
#
#
# tree_engine = FastTreePipe
#
#
# class _SequenceAlignmentTree(_SequenceBunch):
#     def __init__(self):
#         super().__init__()
#         self._nwk_tree = None
#
#     def build_tree(self, args=None):
#         if tree_engine.supports_stdin:  # FastTree
#             tree_builder = tree_engine()
#             tree_builder.stdin = self.as_fasta()
#         else:  # stupid RAxML
#             tempdir = tempfile.gettempdir()
#             tempf = os.path.join(tempdir, str(uuid.uuid1().int) + '.temp.fasta')
#             with open(tempf, 'w') as fh:
#                 print(SimpleFastaWriter(self.as_fasta()), file=fh, end='')
#             tree_builder = tree_engine(input=tempf)
#
#         # if self.is_
#         # TODO !!!@!
#
#     def get_newick_tree(self):
#         if not self._nwk_tree:
#             self.build_tree()
#         return self._nwk_tree
#
#
# class CoreGenesNewickTree2(AlignedCoreGenes):
#     def __init__(self, cg_aligned):
#         super().__init__()
#         if isinstance(cg_aligned, AlignedCoreGenes):
#             self.__dict__.update(cg_aligned.__dict__)
#
#         try:
#             self.tree_program_info = (tree_engine.PROG_NAME, tree_engine().get_version())
#         except (FileNotFoundError, RuntimeError):
#             self.tree_program_info = None
#
#     def _draw_tree(self, gene=None, args=None):
#         if gene is None:
#             gene = CAT
#
#         if tree_engine.supports_stdin:  # FastTree
#             tree_builder = tree_engine()
#             tree_builder.stdin = self.get_fasta(gene)
#         else:  # stupid RAxML
#             tempdir = tempfile.gettempdir()
#             tempf = os.path.join(tempdir, str(uuid.uuid1().int) + '.temp.fasta')
#             with open(tempf, 'w') as fh:
#                 print(SimpleFastaWriter(self[gene]['fasta']), file=fh, end='')
#             tree_builder = tree_engine(input=tempf)
#
#         if self._is_nucleotide_alignment:
#             tree_builder.seq_type = 'protein'
#         else:
#             tree_builder.seq_type = 'nucleotide'
#         # TODO custom model
#         nwk_tree = tree_builder.produce_tree()  # type: str
#         assert nwk_tree
#         self[gene]['tree'] = nwk_tree
#         return nwk_tree
#
#     def fill_trees(self, *, args=None, workers=1, report_progress=None):
#         """TODO docstring"""
#         progress_reporter = None
#         if report_progress:
#             try:
#                 progress_reporter = report_progress(len(self), head_message='Tree  : ')
#             except TypeError:
#                 progress_reporter = report_progress(len(self))
#             if not hasattr(progress_reporter, 'send'):
#                 raise AttributeError
#
#         with ThreadPoolExecutor(workers) as pool:
#             futures = [pool.submit(self._draw_tree, gene, args) for gene in self]
#             try:
#                 for future in as_completed(futures):
#                     try:
#                         future.result()
#                     except Exception:
#                         raise
#                     try:
#                         progress_reporter.send('')
#                     except (AttributeError, StopIteration):
#                         pass
#             except KeyboardInterrupt:
#                 pool._threads.clear()
#                 concurrent.futures.thread._threads_queues.clear()
#                 raise
#
#     def add_trees(self):
#         pass
#
#     def get_tree(self, gene=None):
#         tree = self[gene if gene else CAT]['tree']  # type: str
#         if not tree:
#             return self._draw_tree(gene)
#         else:
#             return tree
#
#
#
#
#
#
# class CoreGenesNewickTree(OrderedDict):
#     def __init__(self, data, seq_type='protein'):
#         super().__init__()
#         if isinstance(data, AlignedCoreGenes):
#             if not data.is_aligned:
#                 raise ValueError
#             d = data.as_dict_with_fasta(seq_type)
#             self.uid_label = data.uid_label
#         elif isinstance(data, dict):  # testing purpose
#             d = data
#             self.uid_label = None
#         else:
#             raise ValueError
#
#         for gene, fasta in d.items():  # includes 'Concatenated'
#             self[gene] = {}
#             self[gene]['fasta'] = fasta
#             self[gene]['tree'] = None
#
#         if seq_type in ('protein', 'prot', 'pro', 'aa'):
#             self.sequence_type = 'protein'
#             self.is_protein_seq = True
#         elif seq_type in ('nucleotide', 'nucl', 'nuc', 'nt'):
#             self.sequence_type = 'nucleotide'
#             self.is_protein_seq = False
#         elif seq_type in ('codon', 'codon12'):
#             self.sequence_type = seq_type
#             self.is_protein_seq = False
#         else:
#             raise ValueError('Invalid sequence_type: {}'.format(seq_type))
#
#         self.__total = len(self)
#         self._n_genome = 0
#         self.tree_engine = tree_engine
#         try:
#             self.tree_program_info = (self.tree_engine.PROG_NAME, self.tree_engine().get_version())
#         except (FileNotFoundError, RuntimeError):
#             self.tree_program_info = None
#
#     @classmethod
#     def from_json(cls, data, seq_type='protein'):
#         import json
#         try:
#             d = json.loads(data)
#         except TypeError:
#             d = json.load(data)
#
#         return cls(d, seq_type)
#
#     def _draw_tree(self, gene=None, args=None):
#         if gene is None:
#             gene = CAT
#
#         if tree_engine.supports_stdin:  # FastTree
#             tree_builder = self.tree_engine()
#             tree_builder.stdin = self[gene]['fasta']
#         else:  # stupid RAxML
#             tempdir = tempfile.gettempdir()
#             tempf = os.path.join(tempdir, str(uuid.uuid1().int) + '.temp.fasta')
#             with open(tempf, 'w') as fh:
#                 print(SimpleFastaWriter(self[gene]['fasta']), file=fh, end='')
#             tree_builder = self.tree_engine(input=tempf)
#
#         if self.is_protein_seq:
#             tree_builder.seq_type = 'protein'
#         else:
#             tree_builder.seq_type = 'nucleotide'
#         # TODO custom model
#         nwk_tree = tree_builder.produce_tree()  # type: str
#         assert nwk_tree
#         self[gene]['tree'] = nwk_tree
#         return nwk_tree
#
#     def fill_trees(self, *, args=None, workers=1, report_progress=None):
#         """TODO docstring"""
#         progress_reporter = None
#         if report_progress:
#             try:
#                 progress_reporter = report_progress(self.__total, head_message='Tree  : ')
#             except TypeError:
#                 progress_reporter = report_progress(self.__total)
#             if not hasattr(progress_reporter, 'send'):
#                 raise AttributeError
#
#         with ThreadPoolExecutor(workers) as pool:
#             futures = [pool.submit(self._draw_tree, gene, args) for gene in self]
#             try:
#                 for future in as_completed(futures):
#                     try:
#                         future.result()
#                     except Exception:
#                         raise
#                     try:
#                         progress_reporter.send('')
#                     except (AttributeError, StopIteration):
#                         pass
#             except KeyboardInterrupt:
#                 pool._threads.clear()
#                 concurrent.futures.thread._threads_queues.clear()
#                 raise
#
#     def add_comment(self, field, value):
#         self.setdefault('__comment__', {})[field] = value
#         self.move_to_end('__comment__', last=False)
#
#     def get_tree(self, gene=None):
#         tree = self[gene if gene else CAT]['tree']  # type: str
#         if not tree:
#             return self._draw_tree(gene)
#         else:
#             return tree
#
#     def __iter__(self):
#         for k in super().__iter__():
#             if k.startswith('__'):
#                 continue
#             yield k
#
#     def get_fasta(self, gene=None):
#         return self[gene if gene else CAT]['fasta']
#
#     def number_of_genomes(self):
#         if not self._n_genome:
#             self._n_genome = sum(1 for s in self.get_fasta().splitlines() if s.startswith('>'))
#         return self._n_genome
#
#     def calculate_gene_support_index(self, threshold=.95):
#         n_genome = self.number_of_genomes()
#         all_gene_trees = [self.get_tree(gene) for gene in self if gene != CAT]
#         ba = BranchAnalysis(all_gene_trees)
#         marked = ba.mark_tree(self.get_tree(),
#                               is_rooted=False, is_core=True,
#                               missing_tolerance=-1, mismatch_tolerance=int(n_genome*(1.-threshold)))
#         self[CAT]['tree_gsi_marked'] = marked
#         return marked
#
#     def get_tree_gsi(self):
#         return self[CAT].get('tree_gsi_marked')
#
#     def as_json(self, indent=None):
#         import json
#         return json.dumps(self, indent=indent)
#
#
# def import_data(data, sequence_type='protein'):
#     try:  # data is AlignedCoreGenes (or dict)
#         return CoreGenesNewickTree(data, sequence_type)
#     except ValueError:  # data is JSON string (for test purpose)
#         return CoreGenesNewickTree.from_json(data, sequence_type)
