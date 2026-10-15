import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline

c = pd.read_csv('data/counts_hidden.csv', index_col=0)
m = pd.read_csv('data/meta_tissues_hidden.csv').rename(columns={'Unnamed: 0':'sample'})
m = m.set_index('sample').loc[c.columns]
y = m['Tissue_Number'].values

# CPM + log2 normalization
lib = c.sum(axis=0)
cpm = c.div(lib, axis=1) * 1e6
logcpm = np.log2(cpm + 1)

# feature selection: top variance genes
var = logcpm.var(axis=1).sort_values(ascending=False)
top = var.head(2000).index
X = logcpm.loc[top].T.values

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
rf = RandomForestClassifier(n_estimators=300, random_state=42, n_jobs=-1)
lr = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000, random_state=42))

for name, model in [('RF', rf), ('LR', lr)]:
    scores = cross_val_score(model, X, y, cv=cv, scoring='accuracy')
    print(f'{name}: acc={scores.mean():.4f} +/- {scores.std():.4f}')

# fit RF on all data, get feature importances
rf.fit(X, y)
imp = pd.Series(rf.feature_importances_, index=top)
imp = imp.sort_values(ascending=False)
print('\nTop 40 marker genes:')
print(imp.head(40).index.tolist())
