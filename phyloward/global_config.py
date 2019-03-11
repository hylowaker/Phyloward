import json
import os

if os.name == 'posix':
    env_paths = os.environ.get('PATH').split(':')
elif os.name == 'nt':
    env_path = os.environ.get('PATH').split(';')
conf = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'config.json')

# default path settings
try:
    with open(conf) as fh:
        conf = json.load(fh)
except OSError:
    raise
    conf = {
        "prodigal_path": "prodigal",
        "hmmsearch_path": "hmmsearch",
        "mafft_path": "mafft",
        "fasttree_path": "FastTree",
        "raxml_path": "raxmlHPC",
    }

for k in ['prodigal_path', 'hmmsearch_path', 'mafft_path', 'fasttree_path', 'raxml_path']:
    if os.path.isfile(conf[k]):
        continue
    for env_path in env_paths:
        fpath = os.path.join(os.path.expanduser(env_path), conf[k])
        if os.path.isfile(fpath):
            conf[k] = fpath
            break

prodigal_path = conf['prodigal_path']
hmmsearch_path = conf['hmmsearch_path']
mafft_path = conf['mafft_path']
fasttree_path = conf['fasttree_path']
raxml_path = conf['raxml_path']
