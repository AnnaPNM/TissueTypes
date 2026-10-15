import pandas as pd
import numpy as np
import requests

c = pd.read_csv('data/counts_hidden.csv', index_col=0)
m = pd.read_csv('data/meta_tissues_hidden.csv').rename(columns={'Unnamed: 0':'sample'})
m = m.set_index('sample').loc[c.columns]
y = m['Tissue_Number'].values

lib = c.sum(axis=0)
cpm = c.div(lib, axis=1) * 1e6
logcpm = np.log2(cpm + 1)

def strip_ver(g): return g.split('.')[0]

all_markers = {}
for t in sorted(set(y)):
    mask = y == t
    mean_t = logcpm.loc[:, mask].mean(axis=1)
    mean_o = logcpm.loc[:, ~mask].mean(axis=1)
    diff = (mean_t - mean_o).sort_values(ascending=False)
    all_markers[t] = diff.head(30).index.tolist()

# collect all unique gene ids
all_ids = set()
for v in all_markers.values():
    all_ids.update(strip_ver(g) for g in v)
all_ids = sorted(all_ids)

# batch query mygene
symbols = {}
for i in range(0, len(all_ids), 100):
    batch = all_ids[i:i+100]
    r = requests.post('https://mygene.info/v3/query', json={'q': batch, 'scopes':'ensembl.gene','fields':'symbol','species':'human'}, timeout=60)
    for d in r.json():
        symbols[d.get('query')] = d.get('symbol')

for t, genes in all_markers.items():
    print(f'=== Tissue {t} ===')
    for g in genes:
        print(f'  {g} -> {symbols.get(strip_ver(g))}')
    print()
