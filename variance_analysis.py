"""
RR-8 / RRRM-1 Liver Transcriptional Variance Analysis
GLDS-379 (OSD-379) - NASA GeneLab

Question: Does spaceflight (FLT) increase per-gene expression variance
relative to Ground Control (GC) in mouse liver, and does this differ
by age (YNG/OLD) or persist after return to Earth (LAR)?
"""
import pandas as pd
import numpy as np
from scipy import stats
import re
import json

np.random.seed(42)

# ---------- Load & parse ----------
df = pd.read_csv('/mnt/user-data/uploads/GLDS-379_rna_seq_Normalized_Counts_rRNArm_GLbulkRNAseq.csv')
df = df.rename(columns={df.columns[0]: 'gene_id'})
df = df.set_index('gene_id')

pattern = re.compile(r'RR8_LVR_([A-Z]+)_(ISS-T|LAR)_([A-Z]+)_([A-Za-z0-9]+)')
meta = []
for c in df.columns:
    m = pattern.match(c)
    group, timepoint, age, sid = m.groups()
    meta.append({'col': c, 'group': group, 'timepoint': timepoint, 'age': age})
meta_df = pd.DataFrame(meta).set_index('col')

# ---------- Expression filter ----------
# Remove genes with near-zero mean expression across ALL samples (avoids the
# division-by-near-zero-variance artifact seen previously with contaminating
# skeletal-muscle transcripts in a liver dataset).
mean_expr_all = df.mean(axis=1)
MIN_MEAN_EXPR = 10.0
df_filt = df[mean_expr_all >= MIN_MEAN_EXPR].copy()

print(f"Genes before filter: {df.shape[0]}")
print(f"Genes after mean-expression filter (>= {MIN_MEAN_EXPR}): {df_filt.shape[0]}")
print()

def get_cols(group, timepoint, age):
    sel = meta_df[(meta_df.group == group) & (meta_df.timepoint == timepoint) & (meta_df.age == age)]
    return list(sel.index)

def levene_variance_test(dataA_cols, dataB_cols, label_A, label_B, n_perm=1000, seed=1):
    """
    Per-gene Levene's test (center='median', robust to non-normality) comparing
    variance of group A vs group B across all genes in df_filt.
    Returns per-gene results + a permutation test on the summary statistic
    (fraction of genes with higher variance in A than B).
    """
    A = df_filt[dataA_cols].values
    B = df_filt[dataB_cols].values
    n_genes = A.shape[0]

    levene_p = np.empty(n_genes)
    var_A = A.var(axis=1, ddof=1)
    var_B = B.var(axis=1, ddof=1)

    for i in range(n_genes):
        try:
            stat, p = stats.levene(A[i, :], B[i, :], center='median')
        except Exception:
            p = np.nan
        levene_p[i] = p

    frac_higher_in_A = np.mean(var_A > var_B)
    sig_mask = levene_p < 0.05
    n_sig = int(np.nansum(sig_mask))
    frac_sig_higher_in_A = np.mean(var_A[sig_mask] > var_B[sig_mask]) if n_sig > 0 else np.nan

    # Wilcoxon signed-rank test on log-variance ratio (paired, gene-by-gene),
    # testing whether the genome-wide shift in variance is directional.
    eps = 1e-9
    log_ratio = np.log2((var_A + eps) / (var_B + eps))
    log_ratio_finite = log_ratio[np.isfinite(log_ratio)]
    try:
        wstat, wp = stats.wilcoxon(log_ratio_finite)
    except Exception:
        wstat, wp = np.nan, np.nan

    # Permutation test: pool A+B samples, reshuffle into groups of the same
    # sizes as A and B, recompute frac_higher_in_A_perm, build null distribution.
    pooled_cols = list(dataA_cols) + list(dataB_cols)
    nA = len(dataA_cols)
    pooled = df_filt[pooled_cols].values
    rng = np.random.default_rng(seed)
    perm_stats = np.empty(n_perm)
    for p_i in range(n_perm):
        idx = rng.permutation(pooled.shape[1])
        a_idx = idx[:nA]
        b_idx = idx[nA:]
        pv_a = pooled[:, a_idx].var(axis=1, ddof=1)
        pv_b = pooled[:, b_idx].var(axis=1, ddof=1)
        perm_stats[p_i] = np.mean(pv_a > pv_b)

    observed = frac_higher_in_A
    perm_p = (np.sum(perm_stats >= observed) + 1) / (n_perm + 1)  # one-sided

    result = {
        'comparison': f'{label_A} (n={len(dataA_cols)}) vs {label_B} (n={len(dataB_cols)})',
        'n_genes_tested': n_genes,
        'frac_genes_higher_var_in_A': round(float(frac_higher_in_A), 4),
        'n_genes_levene_p<0.05': n_sig,
        'frac_sig_genes_higher_var_in_A': round(float(frac_sig_higher_in_A), 4) if n_sig > 0 else None,
        'permutation_null_mean': round(float(perm_stats.mean()), 4),
        'permutation_null_std': round(float(perm_stats.std()), 4),
        'permutation_p_value': round(float(perm_p), 5),
        'wilcoxon_signed_rank_p': float(wp),
    }
    return result, levene_p, var_A, var_B

# ---------- Main comparisons: FLT vs GC, ISS-Terminal, by age ----------
results = {}

for age in ['YNG', 'OLD']:
    flt_cols = get_cols('FLT', 'ISS-T', age)
    gc_cols = get_cols('GC', 'ISS-T', age)
    res, levene_p, var_flt, var_gc = levene_variance_test(flt_cols, gc_cols, 'FLT', 'GC', n_perm=1000, seed=1)
    results[f'FLT_vs_GC_ISS-T_{age}'] = res
    # save top genes for young group
    if age == 'YNG':
        top_idx = np.argsort(levene_p)
        top_genes_yng = df_filt.index[top_idx[:20]]
        top_p_yng = levene_p[top_idx[:20]]
        top_ratio_yng = (var_flt[top_idx[:20]] / (var_gc[top_idx[:20]] + 1e-9))

# ---------- Follow-up: FLT vs GC at LAR (post-return) ----------
for age in ['YNG', 'OLD']:
    flt_cols = get_cols('FLT', 'LAR', age)
    gc_cols = get_cols('GC', 'LAR', age)
    res, _, _, _ = levene_variance_test(flt_cols, gc_cols, 'FLT', 'GC', n_perm=1000, seed=2)
    results[f'FLT_vs_GC_LAR_{age}'] = res

# ---------- Negative control: GC vs VIV at ISS-T (both non-flight) ----------
for age in ['YNG', 'OLD']:
    gc_cols = get_cols('GC', 'ISS-T', age)
    viv_cols = get_cols('VIV', 'ISS-T', age)
    res, _, _, _ = levene_variance_test(gc_cols, viv_cols, 'GC', 'VIV', n_perm=1000, seed=3)
    results[f'NEGCTRL_GC_vs_VIV_ISS-T_{age}'] = res

print(json.dumps(results, indent=2))

# Save top gene table (young, ISS-T, FLT vs GC) for the report
top_df = pd.DataFrame({
    'gene_id': top_genes_yng,
    'levene_p': top_p_yng,
    'var_ratio_FLT_over_GC': top_ratio_yng
})
top_df.to_csv('/home/claude/top_variance_genes_YNG_ISS-T.csv', index=False)

with open('/home/claude/variance_results_summary.json', 'w') as f:
    json.dump(results, f, indent=2)

print()
print("Top 20 variance-shifted genes (YNG, ISS-T, FLT vs GC) saved.")
