# import sys
# import json
# from collections import OrderedDict
#
#
# def from_legacy(s):
#     """Convert legacy UBCG bcg format to phyloward bacterial JSON format"""
#     d = json.loads(s)
#
#     uid = int(d[0]['uid'])
#     label = d[1]['label']
#     # _acc = d[2]['accession']
#     # _tax = d[3]['taxon_name']
#     # _ncb = d[4]['ncbi_name']
#     # _str = d[5]['strain_name']
#     # _typ = d[6]['strain_type']
#     # _prp = d[7]['strain_property']
#     # _txn = d[8]['taxonomy']
#     ubcg_genenum_ver = d[9]['UBCG_target_gene_number|version']
#     n_ubcg = d[10]['n_ubcg']
#     n_genes = d[11]['n_genes']
#     # _nparalog = d[12]['n_paralog_ubcg']
#     # _struc = d[13]['data_structure']
#
#     new_info = OrderedDict()
#     new_info['label'] = label
#     new_info['uid'] = uid
#     new_info['domain'] = 'Bacteria'
#     new_info['program'] = 'phyloward' + ' ##Converted from legacy bcg file (UBCG [{}])'.format(ubcg_genenum_ver)
#     new_info['created_on'] = None
#     new_info['n_core_genes'] = n_ubcg
#     new_info['n_hits'] = n_genes
#
#     data = d[14]['data']
#     new_data = OrderedDict()
#     for gene, array in data.items():
#         new_data[gene] = []
#         for hit in array[1:]:
#             nt_seq = hit[1]
#             aa_seq = hit[2]
#             evalue = float(hit[3])
#             new_data[gene].append({
#                 'id': None,
#                 'evalue': evalue,
#                 'bitscore': None,
#                 'is_included': True,
#                 'protein': aa_seq,
#                 'nucleotide': nt_seq
#             })
#
#     newd = OrderedDict(info=new_info, data=new_data)
#     result = json.dumps(newd, indent=4)
#     return result
#
#
# def to_legacy(s):
#     """Convert phyloward bacterial JSON format to legacy UBCG bcg format"""
#     d = json.loads(s)
#     if d['info']['domain'].casefold() != 'bacteria':
#         raise ValueError('Not bacteria. Cannot convert to bcg format.')
#
#     t = target = []
#     t.append({'uid': str(d['info']['uid'])})  # 0
#     t.append({'label': d['info']['label']})
#     t.append({'accession': None})
#     t.append({'taxon_name': None})
#     t.append({'ncbi_name': None})
#     t.append({'strain_name': None})  # 5
#     t.append({'strain_type': None})
#     t.append({'strain_property': None})
#     t.append({'taxonomy': None})
#     from phyloward import __version__ as ver, __name__ as name, CORE_GENES_BACT as cgl
#     t.append({'UBCG_target_gene_number|version': '{}|##converted from {} {}'.format(len(cgl), name, ver)})
#     t.append({'n_ubcg': d['info']['n_core_genes']})  # 10
#     t.append({'n_genes': d['info']['n_hits']})
#     t.append({'n_paralog_ubcg': d['info']['n_hits'] - d['info']['n_core_genes']})
#     t.append({"data_structure": {"gene_name": ["n_genes", ["feature_index", "dna", "protein", "evalue"]]}})  # 13
#
#     data = {}
#     for gene, hits in d['data'].items():
#         data[gene] = [sum(1 for hit in hits if hit['is_included']), ]
#         for hit in hits:
#             if not hit['is_included']:
#                 continue
#             nt_seq = hit['nucleotide']
#             aa_seq = hit['protein']
#             evalue = hit['evalue']
#             data[gene].append([None, nt_seq, aa_seq, str(evalue).upper()])
#
#     t.append({'data': data})
#     result = json.dumps(t)
#     return result
#
#
# if __name__ == '__main__':
#     print('Convert extracted JSON to legacy BCG format (give -l option to convert opposite)', file=sys.stderr)
#     if sys.argv[1] in ('-l', '--legacy'):
#         f = sys.argv[2]
#         with open(f) as fh:
#             print(from_legacy(fh.read()))
#     else:
#         f = sys.argv[1]
#         with open(f) as fh:
#             print(to_legacy(fh.read()))
