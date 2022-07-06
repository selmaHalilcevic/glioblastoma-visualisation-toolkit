import os
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
plt.rcParams['figure.dpi'] = 500
# plt.rcParams['font.size'] = 12



def get_pvalues_significance_annotations(r_df, n):
    significance_level = 0.05
    r_df = r_df.loc[:, ~r_df.columns.duplicated()]
    pvalues_df = r_df.copy()
    for col in pvalues_df.columns:
        l = list(pvalues_df[col])
        pvalues = []
        for r in l:
            if r == 1.0 or r == -1.0:
                pvalues.append(0.0)
            else:
                t_statistic = (r * np.sqrt(n - 2)) / np.sqrt(1 - r * r)
                pval = stats.t.sf(np.abs(t_statistic), n - 2) * 2  # degrees of freedom = n - 2  # two-sided p value
                pvalues.append(pval)
        pvalues_df[col] = pvalues

    temp = r_df.to_numpy()

    annot = [[f"{temp[row][col]:.3f}"
              + ('\n★' if val <= significance_level else '')
              for col, val in enumerate(column)] for row, column in enumerate(pvalues_df.to_numpy())]

    return r_df, annot


def get_heatmaps(workspace_path, data, img_name, img_size=(10, 6), vmin=-1, vmax=1,
                 title=None):
    plt.figure(figsize=img_size)
    n = data.shape[0]  # sample size
    r_df = data.corr(method='pearson')
    _, annot = get_pvalues_significance_annotations(r_df, n)
    cmap = sns.diverging_palette(230, 0, 90, 60, as_cmap=True)
    ax = sns.heatmap(r_df, cmap=cmap, annot=annot, vmin=vmin, vmax=vmax, fmt='',
                     cbar_kws={
                         'label': f'Pearson\'s correlations coefficient for n = {n} samples, ★ marks p-value $\leq$ 0.05'})
    ax.set_xticklabels(
        ax.get_xticklabels(),
        rotation=45,
        horizontalalignment='right'
    )
    ax.set_yticklabels(
        ax.get_yticklabels(),
        rotation=0
    )
    if title != None: plt.title(f'{title} - Correlations')
    plt.savefig(os.path.join(workspace_path, 'results', 'correlations', f'{img_name}.png'), bbox_inches='tight')
    plt.close()


def get_dummy_variables(df, column_name, new_columns, mutation):
    df_dummies = pd.get_dummies(df[column_name])
    df_new = pd.concat([df, df_dummies], axis=1)
    for col in new_columns:
        df_new = df_new.rename(columns={col: mutation + " - " + col})
    return df_new


def split_mutations_to_dummy_variables(df):
    prev_cols = df.columns
    df_new = df
    mutations = ['Co-deletion 1p/19q', 'ATRX', 'p53', 'MGMT methylation', 'IDH R132H', 'IDH1', 'BRAF',
                 'EGFR', 'TERT']
    for m in mutations:
        df_new = get_dummy_variables(df_new, column_name=m,
                                     new_columns=list(set(list(df_new[m]))),
                                     mutation=m)
    new_cols = df_new.columns
    n = [c for c in new_cols if c not in prev_cols]

    return df_new, n


def split_location_to_dummy_vars(df, col_name_loc, col_name_contact):
    prev_cols = df.columns
    values_side = []
    for side in df[col_name_loc]:
        if side == "r" or side == "R":
            values_side.append('right')
        elif side == "l" or side == "L":
            values_side.append('left')
        elif side == "r,l" or side == 'l,r' or side == 'R,L' or side == 'L,R':
            values_side.append('right and left')
        elif side == 'bill':
            values_side.append('bill')
    df[col_name_loc] = values_side

    if df[col_name_loc].isnull().values.any() or df[col_name_contact].isnull().values.any():
        df = df.dropna()
    df_side = get_dummy_variables(df, column_name=col_name_loc,
                                  new_columns=list(set(list(df[col_name_loc]))),
                                  mutation='Tumor Side')

    new_cols_side = [c for c in df_side.columns if c not in prev_cols and 'bill' not in c]
    df = df[prev_cols]
    values_side = []
    for contact in df[col_name_contact]:
        if contact == "cortex":
            values_side.append('no')
        elif contact == "ventricle":
            values_side.append('yes')
        elif contact == "cortex+ventricle" or contact == 'ventricle+cortex':
            values_side.append('yes')
        elif contact == 'none':
            values_side.append('none')
    df[col_name_contact] = values_side
    df_contact = get_dummy_variables(df, column_name=col_name_contact,
                                     new_columns=list(set(list(df[col_name_contact]))),
                                     mutation='Ventricle Contact')

    new_cols_contact = [c for c in df_contact.columns if c not in prev_cols and 'none' not in c]
    df[new_cols_side] = df_side[new_cols_side]
    df[new_cols_contact] = df_contact[new_cols_contact]
    return df, new_cols_side + new_cols_contact


def plot_mutations_with_combinations(wp, df, genes, filename, img_size, add_cols=None, col_map='Blues',
                                     add_title=False):
    mutations = ['Co-deletion 1p/19q', 'ATRX', 'p53', 'MGMT methylation', 'EGFR', 'BRAF', 'TERT']
    for m in mutations:
        data = []
        for col in add_cols:
            if 'other mutation' in col:
                continue
            if m in col:
                data.append(col)
        if len(data) == 0:
            continue
        n_with_sample_no = []
        for old_col_name in data:
            total_samples = (df[old_col_name] == 1).sum()
            new_col_name = old_col_name + f" (n = {total_samples})"
            df = df.rename(columns={old_col_name: new_col_name})
            n_with_sample_no.append(new_col_name)
        data = ['Overall survival'] + n_with_sample_no + genes
        mut = re.sub("[ ?.!/;:]", '_', m)
        f = re.sub("[ ?.!;:_]", ' ', filename)
        title = None if not add_title else f'{f.capitalize()} and {mut} mutation'
        get_heatmaps(wp, df[data],
                     img_name=os.path.join('mutations and gene markers combinations', f'{filename}_{mut}'),
                     img_size=img_size, vmin=-1, vmax=1, title=title)

    mutations = ['IDH R132H', 'IDH1']
    data = []
    for m in mutations:
        for col in add_cols:
            if 'other mutation' in col:
                continue
            if m in col:
                if np.all((df[col] == 0)):
                    continue
                else:
                    data.append(col)
    n_with_sample_no = []
    for old_col_name in data:
        total_samples = (df[old_col_name] == 1).sum()
        new_col_name = old_col_name + f" (n = {total_samples})"
        df = df.rename(columns={old_col_name: new_col_name})
        n_with_sample_no.append(new_col_name)
    data = ['Overall survival'] + n_with_sample_no + genes
    f = re.sub("[ ?.!;:_]", ' ', filename)
    title = None if not add_title else f'{f.capitalize()} and IDH1 mutation'
    get_heatmaps(wp, df[data],
                 img_name=os.path.join('mutations and gene markers combinations', f'{filename}_IDH1_mutations'),
                 img_size=img_size, vmin=-1, vmax=1, title=title)


def plot_mutations_with_location(df, mutation_cols=None, location_cols=None, wp=None, col_map='Blues', img_size=None,
                                 add_title=False):
    loc_col_with_sample_no = []
    for old_col_name in location_cols:
        total_samples = (df[old_col_name] == 1).sum()
        new_col_name = old_col_name + f" (n = {total_samples})"
        df[new_col_name] = df[old_col_name]
        loc_col_with_sample_no.append(new_col_name)

    mutations = ['Co-deletion 1p/19q', 'ATRX', 'p53', 'MGMT methylation', 'EGFR', 'BRAF', 'TERT']
    for m in mutations:
        data = []
        for col in mutation_cols:
            if 'other mutation' in col:
                continue
            if m in col:
                data.append(col)
        if len(data) == 0:
            continue
        n_with_sample_no = []
        for old_col_name in data:
            total_samples = (df[old_col_name] == 1).sum()
            new_col_name = old_col_name + f" (n = {total_samples})"
            df = df.rename(columns={old_col_name: new_col_name})
            n_with_sample_no.append(new_col_name)
        data = ['Overall survival'] + n_with_sample_no + loc_col_with_sample_no
        mut = re.sub("[ ?.!/;:]", '_', m)
        title = None if not add_title else f'Tumor location and {mut} mutation'
        get_heatmaps(wp, df[data],
                     img_name=os.path.join('mutations and tumor location', f'tumor_location_{mut}'),
                     img_size=img_size, vmin=-1, vmax=1, title=title)

    # Plot IDH mutations
    mutations = ['IDH R132H', 'IDH1']
    data = []
    for m in mutations:
        for col in mutation_cols:
            if 'other mutation' in col:
                continue
            if m in col:
                if np.all((df[col] == 0)):
                    continue
                else:
                    data.append(col)
    n_with_sample_no = []
    for old_col_name in data:
        total_samples = (df[old_col_name] == 1).sum()
        new_col_name = old_col_name + f" (n = {total_samples})"
        df = df.rename(columns={old_col_name: new_col_name})
        n_with_sample_no.append(new_col_name)
    data = ['Overall survival'] + n_with_sample_no + loc_col_with_sample_no
    title = None if not add_title else f'Tumor location and IDH1 mutation'
    get_heatmaps(wp, df[data],
                 img_name=os.path.join('mutations and tumor location', f'tumor_location_IDH1_mutations'),
                 img_size=img_size, vmin=-1, vmax=1, title=title)


def plot_combinations(workspace_path, data, h_diff, v_diff, img_name,
                      c_map='Blues', annot_=True, img_size=(10, 6), multiple=False, color=['black'], add_title=False):
    plt.figure(figsize=img_size)
    n = data.shape[0]  # sample size
    r_df = data.corr(method='pearson')
    r_df, annot = get_pvalues_significance_annotations(r_df, n)
    cmap = sns.diverging_palette(230, 0, 90, 60, as_cmap=True)
    ax = sns.heatmap(r_df, cmap=cmap, annot=annot, fmt='', vmin=-1, vmax=1,
                     cbar_kws={
                         'label': f'Pearson\'s correlations coefficient for n = {n} samples, ★ marks p-value $\leq$ 0.05'})
    if not multiple:
        ax.hlines([h_diff[0]], *ax.get_xlim(), colors=[color[0]], linestyles='dashdot')
        ax.vlines([v_diff[0]], *ax.get_ylim(), colors=[color[0]], linestyles='dashdot')
        ax.set_xticklabels(
            ax.get_xticklabels(),
            rotation=45,
            horizontalalignment='right'
        )
        ax.set_yticklabels(
            ax.get_yticklabels(),
            rotation=0
        )
    else:
        if len(h_diff) == len(v_diff):
            color = color * len(h_diff) if len(color) == 1 else color
            for i in range(0, len(h_diff)):
                ax.hlines([h_diff[i]], *ax.get_xlim(), colors=[color[i]], linestyles='dashdot')
                ax.vlines([v_diff[i]], *ax.get_ylim(), colors=[color[i]], linestyles='dashdot')
                ax.set_xticklabels(
                    ax.get_xticklabels(),
                    rotation=45,
                    horizontalalignment='right'
                )
                ax.set_yticklabels(
                    ax.get_yticklabels(),
                    rotation=0
                )
    f = re.sub("[ ?.!;:_]", ' ', img_name)
    title = None if not add_title else f'{f.capitalize()}'
    plt.savefig(os.path.join(workspace_path, 'results', 'correlations', 'gene markers combinations', f'{img_name}.png'),
                bbox_inches='tight')
    plt.close()
