import os
import warnings

import pandas as pd

warnings.filterwarnings("ignore")

from utils.kaplan_meier import plot_kaplan_meier, plot_location, plot_contact, \
    plot_subtypes, plot_combinations, plot_gender
from utils.variables import genes, combinations

workspace_path = "./"


def plot_genes(gene_name, T, E, vals, wp, show_censors=True, show_title=True, at_risk=True):
    quartiles = {0.25: 'Q1', 0.50: 'Q2', 0.75: 'Q3'}
    img_name = os.path.join('All genes', gene_name)
    for q in quartiles:
        plot_kaplan_meier(vals, T, E, percentile=q, var_name=gene_name, img_name=f"{img_name}_{quartiles[q]}",
                          workspace_path=wp, add_threshold=True, show_censors=show_censors,
                          show_plot_title=show_title, at_risk=at_risk, threshold_three_decimals=True)


def plot_age(T, E, vals, wp, show_censors=True, show_title=True, at_risk=False):
    quartiles = {0.25: 'Q1', 0.50: 'Q2', 0.75: 'Q3'}
    img_name = os.path.join('Age, gender, KPS, location, subtypes', 'age_at_the_time_of_1st_diagnosis')
    for q in quartiles:
        plot_kaplan_meier(vals, T, E, percentile=q, var_name='Age', img_name=f"{img_name}_{quartiles[q]}",
                          workspace_path=wp, add_threshold=True, show_censors=show_censors,
                          show_plot_title=show_title, at_risk=at_risk)


def plot_kps(T, E, post_kps, pre_kps, wp, show_censors=True, show_title=True, at_risk=False):
    quartiles = {0.25: 'Q1', 0.50: 'Q2', 0.75: 'Q3'}
    for q in quartiles:
        plot_kaplan_meier(pre_kps, T, E, percentile=q, var_name='Pre-KPS',
                          img_name=os.path.join('Age, gender, KPS, location, subtypes',
                                                f"pre-operative-kps-{quartiles[q]}"),
                          workspace_path=wp, add_threshold=True, show_censors=show_censors,
                          show_plot_title=show_title, at_risk=at_risk)
        plot_kaplan_meier(post_kps, T, E, percentile=q, var_name='Post-KPS',
                          img_name=os.path.join('Age, gender, KPS, location, subtypes',
                                                f"post-operative-kps-{quartiles[q]}"),
                          workspace_path=wp, add_threshold=True, show_censors=show_censors,
                          show_plot_title=show_title, at_risk=at_risk)


def main(workspace='../'):
    global workspace_path
    workspace_path = workspace

    df_gene = pd.read_csv(os.path.join(workspace_path, 'data', 'gene_expression_dataset.csv'), sep=';')
    df_clinical = pd.read_csv(os.path.join(workspace_path, 'data', 'clinical_dataset.csv'), sep=';')
    df_pato = pd.read_csv(os.path.join(workspace_path, 'data', 'pathohistological_dataset.csv'), sep=';')
    df = df_gene[['SampleName', 'Subtype']].merge(df_clinical[['SampleName', 'Overall survival', 'Alive']],
                                                  on='SampleName')
    df = df.dropna()
    T = df['Overall survival']
    plot_subtypes(T, E=df['Alive'], groups=df['Subtype'], show_title=True, at_risk=True, workspace_path=workspace_path)

    df = df_pato[
        ['SampleName', 'Side of tumor L, R', 'Ventricle contact (Cortex, cortex+ventricle, ventricle, none)']].merge(
        df_clinical[['SampleName', 'Overall survival', 'Alive']], on='SampleName')
    # df = df.dropna()
    T = df['Overall survival']
    plot_location(T, E=df['Alive'], vals=df['Side of tumor L, R'], workspace_path=workspace_path, show_title=True)
    plot_contact(T, E=df['Alive'], groups=df['Ventricle contact (Cortex, cortex+ventricle, ventricle, none)'],
                 workspace_path=workspace_path, show_title=True)

    plot_age(T=df_clinical['Overall survival'], E=df_clinical['Alive'],
             vals=df_clinical['Age at time of 1st diagnosis'], show_title=True, at_risk=True, wp=workspace_path)
    plot_kps(T=df_clinical['Overall survival'], E=df_clinical['Alive'], post_kps=df_clinical['Post operative KPS'],
             pre_kps=df_clinical['Pre-operative KPS'],
             wp=workspace_path, show_title=True, at_risk=True)
    plot_gender(T=df_clinical['Overall survival'], E=df_clinical['Alive'], vals=df_clinical['Gender (m/f)'],
                workspace_path=workspace_path, show_censors=True, at_risk=True, show_title=True)

    # TODO: COMBINATIONS

    df = df_clinical[['SampleName', 'Overall survival', 'Alive']].merge(
        df_gene, on='SampleName'
    )
    # df = df.dropna()

    # PLOTS ALL GENES
    for gene in genes:
        plot_genes(gene, T=df['Overall survival'], E=df['Alive'], vals=df[gene], wp=workspace_path, show_censors=True,
                   show_title=True, at_risk=True)

    # Plots all genes thresholded by quartiles in one pdf file
    plot_combinations(T=df['Overall survival'], genes=genes, df=df, rows=2, cols=3,
                      workspace_path=workspace_path,
                      filename='all_genes', show_censors=True, at_risk=False)

    # Plots every combination to a single pdf file
    for key in combinations:
        genes_ = combinations[key][0]
        plot_combinations(T=df['Overall survival'], genes=genes_, df=df, rows=2, cols=3,
                          workspace_path=workspace_path,
                          filename=key, show_censors=True, at_risk=True)


if __name__ == "__main__":
    main()
