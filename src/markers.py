import pandas as pd
import numpy as np

c = pd.read_csv('data/counts_hidden.csv', index_col=0)
m = pd.read_csv('data/meta_tissues_hidden.csv').rename(columns={'Unnamed: 0':'sample'})
m = m.set_index('sample').loc[c.columns]
y = m['Tissue_Number'].values

lib = c.sum(axis=0)
cpm = c.div(lib, axis=1) * 1e6
logcpm = np.log2(cpm + 1)

for t in sorted(set(y)):
    mask = y == t
    mean_t = logcpm.loc[:, mask].mean(axis=1)
    mean_o = logcpm.loc[:, ~mask].mean(axis=1)
    diff = (mean_t - mean_o).sort_values(ascending=False)
    print(f'=== Tissue {t} (n={mask.sum()}) top 25 markers ===')
    print(diff.head(25).index.tolist())
    print()
