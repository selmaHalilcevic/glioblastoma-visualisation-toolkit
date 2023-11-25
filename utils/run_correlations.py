import os
import re

import pandas as pd

from utils.correlations import plot_mutations_with_combinations, split_mutations_to_dummy_variables, \
    get_heatmaps, plot_combinations, split_location_to_dummy_vars, plot_mutations_with_location
from utils.variables import combinations, groups, subtypes


def plot_mutations(df, wp):
    df_new, cols_with_mutations_dummy_vars = split_mutations_to_dummy_variables(df)
    df_loc, cols_location = split_location_to_dummy_vars(df_new, 'Side of tumor L, R', 'Ventricle contact - combined')

    loc_col_with_sample_no = []
    for old_col_name in cols_location:
        total_samples = (df_loc[old_col_name] == 1).sum()
        new_col_name = old_col_name + f" (n = {total_samples})"
        df_loc[new_col_name] = df_loc[old_col_name]
        loc_col_with_sample_no.append(new_col_name)

    plot_mutations_with_location(df_loc, mutation_cols=cols_with_mutations_dummy_vars,
                                 location_cols=cols_location, wp=wp, img_size=(15, 12), add_title=False)
    for key in combinations:
        plot_mutations_with_combinations(wp, df=df_new, genes=combinations[key][0],
                                         filename=key, img_size=combinations[key][1],
                                         add_cols=cols_with_mutations_dummy_vars, col_map='Blues', add_title=False)
        data = ['Overall survival'] + loc_col_with_sample_no + combinations[key][0]
        f = re.sub("[_]", ' ', key).capitalize()
        get_heatmaps(wp, df_loc[data],
                     img_name=os.path.join('gene combinations and tumor location', f'tumor_location_{key}'),
                     img_size=combinations[key][1],
                     vmin=-1, vmax=1)# , title=f'Tumor location and {f}')


def plot_corr_combinations(df):
    # Plot correlations between each gene in a combination
    for key in combinations:
        f = re.sub("[_]", ' ', key).capitalize()
        get_heatmaps(workspace_path, data=df[['Overall survival'] + combinations[key][0]],
                     img_name=os.path.join('gene combinations', key),
                     img_size=combinations[key][1], vmin=-1)# , title=f)

    # Correlations between two combinations
    for group in groups:
        data = df[combinations[group[0]][0] + combinations[group[1]][0]]
        h_diff = [len(combinations[group[0]][0])]
        v_diff = h_diff
        name = f"{group[0]}_LEFT_and_{group[1]}_RIGHT"
        plot_combinations(workspace_path, data, h_diff, v_diff,
                          img_name=name, img_size=group[2], add_title=False)

    # Plot correlations between gbm subtypes and overall survival
    data = ['Overall survival']
    [diff, prev, count, names] = [[], 1, 0, ""]
    for key in subtypes:
        names += key
        count += 1
        data = data + subtypes[key]
        if count == len(subtypes): break
        diff.append(len(subtypes[key]) + prev)
        prev = len(subtypes[key]) + 1
        names += "_"
    plot_combinations(workspace_path, data=df[data], h_diff=diff, v_diff=diff, img_name=f'tcga_subtypes_{names}',
                      img_size=(18, 14), multiple=True, color=['black'], add_title=False)


def main(workspace='../'):
    global workspace_path
    workspace_path = workspace

    pato = pd.read_csv(os.path.join(workspace_path, 'data', 'pathohistological_dataset.csv'), sep=';')
    gene_df = pd.read_csv(os.path.join(workspace_path, 'data', 'gene_expression_dataset.csv'), sep=';')
    df = pd.merge(pato, gene_df, on='SampleName')
    plot_mutations(df, workspace_path)
    plot_corr_combinations(df)


if __name__ == "__main__":
    main()
