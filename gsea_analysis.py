"""
Gen Seti Zenginlestirme Analizi (GSEA) - Genc Fare Karacigeri, ISS-Terminal
FLT vs GC karsilastirmasi icin, GO Biyolojik Surec kategorilerine karsi.

Kullanilan gen seti veritabani:
ELTEbioinformatics/GMT_files_for_mulea (GO_BP_Mus_musculus_EnsemblID.gmt)
https://github.com/ELTEbioinformatics/GMT_files_for_mulea
"""
import gseapy as gp
import pandas as pd
import numpy as np
import re

df = pd.read_csv('GLDS-379_rna_seq_Normalized_Counts_rRNArm_GLbulkRNAseq.csv')
df = df.rename(columns={df.columns[0]: 'gene_id'}).set_index('gene_id')

pattern = re.compile(r'RR8_LVR_([A-Z]+)_(ISS-T|LAR)_([A-Z]+)_([A-Za-z0-9]+)')
meta = []
for c in df.columns:
    m = pattern.match(c)
    group, timepoint, age, sid = m.groups()
    meta.append({'col': c, 'group': group, 'timepoint': timepoint, 'age': age})
meta_df = pd.DataFrame(meta).set_index('col')

mean_expr_all = df.mean(axis=1)
df_filt = df[mean_expr_all >= 10.0].copy()
df_log = np.log2(df_filt + 1)

def get_cols(group, timepoint, age):
    sel = meta_df[(meta_df.group == group) & (meta_df.timepoint == timepoint) & (meta_df.age == age)]
    return list(sel.index)

flt_cols = get_cols('FLT', 'ISS-T', 'YNG')
gc_cols = get_cols('GC', 'ISS-T', 'YNG')
A = df_log[flt_cols].values
B = df_log[gc_cols].values
var_A = A.var(axis=1, ddof=1)
var_B = B.var(axis=1, ddof=1)

eps = 1e-9
log_var_ratio = np.log2((var_A + eps) / (var_B + eps))

rank_df = pd.DataFrame({'gene': df_filt.index.astype(str), 'score': log_var_ratio})
rank_df = rank_df.dropna().sort_values('score', ascending=False).reset_index(drop=True)

pre_res = gp.prerank(
    rnk=rank_df,
    gene_sets='GO_BP_Mus_musculus_EnsemblID_clean.gmt',
    min_size=15,
    max_size=500,
    permutation_num=1000,
    outdir=None,
    seed=42,
    threads=4,
    no_plot=True,
)

res = pre_res.res2d.sort_values('FDR q-val')
res.to_csv('gsea_results.csv', index=False)
print(res[['Term', 'ES', 'NES', 'NOM p-val', 'FDR q-val']].head(20).to_string())
