import os
import statistics
import warnings

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.plotting import add_at_risk_counts
from lifelines.statistics import logrank_test, pairwise_logrank_test
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.offsetbox import AnchoredText
plt.rcParams['figure.dpi'] = 500

warnings.filterwarnings("ignore")


def get_groups(E, T, percentile):
    [threshold] = E.quantile([percentile])
    i1 = (E <= threshold)
    i2 = (E > threshold)

    median_i1 = statistics.median(list(T[i1]))
    median_i2 = statistics.median(list(T[i2]))

    return i1, i2, "{:.3f}".format(threshold), "{:.3f}".format(median_i1), "{:.3f}".format(median_i2)


def get_survival_median(kmf_model):
    median = kmf_model.median_survival_time_
    temp = kmf_model.survival_function_
    median_idx = list(temp.index).index(median)
    y_max = list(temp[temp.columns[0]])[median_idx]

    return median, y_max


def check_p_value(p_value):
    if p_value < 0.001:
        return "<0.001"
    else:
        return "{:.3f}".format(p_value)


def plot_kaplan_meier(vals, T, E, percentile, var_name, img_name, workspace_path, add_threshold=True, show_censors=True,
                      show_plot_title=False, at_risk=False, threshold_three_decimals=False):
    matplotlib.rcParams.update({"font.size": 9})
    i1, i2, threshold, median_i1, median_i2 = get_groups(vals, T, percentile)
    fig, ax = plt.subplots()
    kmf1 = KaplanMeierFitter()
    if add_threshold and threshold_three_decimals:
        label_1 = f"{var_name} $\leq$ {threshold}"
    elif add_threshold and not threshold_three_decimals:
        label_1 = f"{var_name} $\leq$ {float(threshold)}"
    else:
        label_1 = var_name
    kmf1.fit(durations=T[i1], event_observed=E[i1],
             label=f'{label_1} (n = {len(T[i1])}, median OS = {float(median_i1)})')
    median1, y_max1 = get_survival_median(kmf1)
    ax = kmf1.plot_survival_function(ax=ax, show_censors=show_censors, censor_styles={'ms': 6, 'marker': 's'},
                                     loc=slice(0., 36.))
    if add_threshold and threshold_three_decimals:
        label_2 = f"{var_name} > {threshold}"
    elif add_threshold and not threshold_three_decimals:
        label_2 = f"{var_name} > {float(threshold)}"
    else:
        label_2 = var_name
    kmf2 = KaplanMeierFitter()
    kmf2.fit(durations=T[i2], event_observed=E[i2],
             label=f'{label_2} (n = {len(T[i2])}, median OS = {float(median_i2)})')
    median2, y_max2 = get_survival_median(kmf2)
    ax = kmf2.plot_survival_function(ax=ax, show_censors=show_censors, censor_styles={'ms': 6, 'marker': 's'},
                                     loc=slice(0., 36.))
    plt.xticks(np.arange(0, 36, step=6))
    plt.ylabel("Survival probability")
    plt.xlabel("Time in months")

    if at_risk: add_at_risk_counts(kmf1, kmf2, labels=[label_1, label_2], ax=ax)

    plt.vlines(x=median1, ymin=0, ymax=y_max1,
               colors='black',
               linestyle='dotted')
    plt.hlines(y=y_max1, xmin=0, xmax=median1,
               colors='black',
               linestyle='dotted')
    plt.vlines(x=median2, ymin=0, ymax=y_max2,
               colors='black',
               linestyle='dotted')
    plt.hlines(y=y_max2, xmin=0, xmax=median2,
               colors='black',
               linestyle='dotted')
    results = logrank_test(durations_A=T[i1], durations_B=T[i2], event_observed_A=E[i1], event_observed_B=E[i2])
    p_value = check_p_value(results.p_value)
    at = AnchoredText(f"p-value: {p_value}",
                      loc='center right', prop=dict(size=12), frameon=True)
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax.add_artist(at)
    q = ""
    if percentile == 0.25:
        q = 1
    elif percentile == 0.50:
        q = 2
    elif percentile == 0.75:
        q = 3
    if show_plot_title: plt.title(f"{var_name}, Q{q}")
    plt.savefig(os.path.join(workspace_path, 'results', 'kaplan_meier_curves', f"{img_name}.png"), bbox_inches='tight')
    plt.close(fig)
    matplotlib.rcParams.update({"font.size": 12})


def plot_subtypes(T, E, groups, show_censors=True, at_risk=True, show_title=False, workspace_path=None):
    matplotlib.rcParams.update({"font.size": 9})
    subtypes = ['CL', 'MES', 'MIX', 'PN']
    i_CL = (groups == 'CL')
    i_MES = (groups == 'MES')
    i_MIX = (groups == 'MIX')
    i_PN = (groups == 'PN')

    fig, ax = plt.subplots()  # figsize=(7,9))

    kmf1 = KaplanMeierFitter()
    kmf1.fit(durations=T[i_CL], event_observed=E[i_CL],
             label=f'CL (n = {len(T[i_CL])}, median OS = {statistics.median(list(T[i_CL]))})')
    ax = kmf1.plot_survival_function(ax=ax, show_censors=show_censors, censor_styles={'ms': 6, 'marker': 's'},
                                     loc=slice(0., 36.), ci_alpha=0.1, ci_show=True)

    kmf2 = KaplanMeierFitter()
    kmf2.fit(durations=T[i_MES], event_observed=E[i_MES],
             label=f'MES (n = {len(T[i_MES])}, median OS = {statistics.median(list(T[i_MES]))})')
    ax = kmf2.plot_survival_function(ax=ax, show_censors=True, censor_styles={'ms': 6, 'marker': 's'},
                                     loc=slice(0., 36.), ci_alpha=0.1, ci_show=True)

    kmf3 = KaplanMeierFitter()
    kmf3.fit(durations=T[i_MIX], event_observed=E[i_MIX],
             label=f'MIX (n = {len(T[i_MIX])}, median OS = {statistics.median(list(T[i_MIX]))})')
    ax = kmf3.plot_survival_function(ax=ax, show_censors=show_censors, censor_styles={'ms': 6, 'marker': 's'},
                                     loc=slice(0., 36.), ci_alpha=0.1, ci_show=True)

    kmf4 = KaplanMeierFitter()
    kmf4.fit(durations=T[i_PN], event_observed=E[i_PN],
             label=f'PN (n = {len(T[i_PN])}, median OS = {statistics.median(list(T[i_PN]))})')
    ax = kmf4.plot_survival_function(ax=ax, show_censors=True, censor_styles={'ms': 6, 'marker': 's'},
                                     loc=slice(0., 36.), ci_alpha=0.1, ci_show=True)
    plt.xticks(np.arange(0, 36, step=6))
    plt.ylabel("Survival probability")
    plt.xlabel("Time in months")
    matplotlib.rcParams.update({"font.size": 7})
    if at_risk: add_at_risk_counts(kmf1, kmf2, kmf3, kmf4, labels=subtypes, ax=ax,
                                   rows_to_show=['At risk', 'Censored', 'Events'])
    matplotlib.rcParams.update({"font.size": 9})

    results2 = pairwise_logrank_test(event_durations=T, groups=groups, event_observed=E)
    results2_df = results2.summary.p.reset_index(drop=False)
    at = AnchoredText(f"P-VALUES:\n{split_multiple_p_values(results2_df)}",
                      loc='center right', prop=dict(size=11), frameon=True)

    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax.add_artist(at)
    if show_title: plt.title('Kaplan-Meier Survival Prediction Based on Four TCGA Subtypes')
    plt.savefig(os.path.join(workspace_path, 'results', 'kaplan_meier_curves',
                             'Age, gender, KPS, location, subtypes', "four_tcga_subtypes.png"), bbox_inches='tight')
    plt.close(fig)
    matplotlib.rcParams.update({"font.size": 9})

    fig2, ax2 = plt.subplots(figsize=(5, 6))
    ax2 = kmf1.plot_survival_function(ax=ax2, show_censors=show_censors, censor_styles={'ms': 6, 'marker': 's'},
                                      loc=slice(0., 36.), ci_alpha=0.1, ci_show=True)
    ax2 = kmf3.plot_survival_function(ax=ax2, show_censors=show_censors, censor_styles={'ms': 6, 'marker': 's'},
                                      loc=slice(0., 36.), ci_alpha=0.1, ci_show=True)

    plt.xticks(np.arange(0, 36, step=6))
    plt.ylabel("Survival probability")
    plt.xlabel("Time in months")
    matplotlib.rcParams.update({"font.size": 8})
    if at_risk: add_at_risk_counts(kmf1, kmf3, labels=['CL', 'MIX'], ax=ax2,
                                   rows_to_show=['At risk', 'Censored', 'Events'])
    matplotlib.rcParams.update({"font.size": 9})
    plt.tight_layout()
    results = logrank_test(durations_A=T[i_CL], durations_B=T[i_MIX],
                           event_observed_A=E[i_CL], event_observed_B=E[i_MIX])
    p_value = check_p_value(results.p_value)
    at = AnchoredText(f"p-value: {p_value}",
                      loc='center right', prop=dict(size=12), frameon=True)
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax2.add_artist(at)
    if show_title: plt.title('CL and MIX TCGA Subtypes')

    plt.savefig(os.path.join(workspace_path, 'results', 'kaplan_meier_curves',
                             'Age, gender, KPS, location, subtypes', "CL_and_MIX_tcga_subtypes.png"),
                bbox_inches='tight')
    plt.close(fig)
    matplotlib.rcParams.update({"font.size": 12})


def helper(m, n, axarr):
    if m > 1:
        arr_ij = [(x, y) for x, y in np.ndindex(axarr.shape)]
        return [axarr[index] for index in arr_ij]  # subplots
    else:
        arr_ij = list(zip([0] * n, list(range(n))))
        return [axarr[index[1]] for index in arr_ij]  # subplots


def plot_combinations(T, genes, df, rows, cols, workspace_path, filename, show_censors=False, at_risk=False):
    matplotlib.rcParams.update({"font.size": 5})
    if at_risk:
        m, n = 1, 3
    else:
        m, n = rows, cols
    with PdfPages(os.path.join(workspace_path, 'results', 'kaplan_meier_curves', 'combinations',
                               f'{filename}.pdf')) as pdf:
        f, axarr = plt.subplots(m, n, sharex="col", sharey="row")
        plt.subplots_adjust(wspace=0.05)
        subplots = helper(m, n, axarr)

        splot_index = 0
        for s, splot in enumerate(subplots):
            last_row = m * n - s < n + 1
            first_in_row = s % n == 0
            if last_row:
                splot.set_xlabel("Time in months")
            if first_in_row:
                splot.set_ylabel("Survival probability")

        for gene in genes:
            E = df['Alive']

            # TODO: 25%
            i1, i2, threshold, median_i1, median_i2 = get_groups(df[gene], T, 0.25)
            kmf1 = KaplanMeierFitter()
            label_1 = f'$\leq$ {threshold}'
            kmf1.fit(durations=T[i1], event_observed=E[i1],
                     label=f'{label_1} (n = {len(T[i1])}, median = {float(median_i1)})')
            median1, y_max1 = get_survival_median(kmf1)
            subplots[splot_index] = kmf1.plot_survival_function(ax=subplots[splot_index], show_censors=show_censors,
                                                                censor_styles={'ms': 5, 'marker': '+'},
                                                                loc=slice(0., 36.))

            label_2 = f'> {threshold}'
            kmf2 = KaplanMeierFitter()
            kmf2.fit(durations=T[i2], event_observed=E[i2],
                     label=f'{label_2} (n = {len(T[i2])}, median = {float(median_i2)})')
            median2, y_max2 = get_survival_median(kmf2)
            subplots[splot_index] = kmf2.plot_survival_function(ax=subplots[splot_index], show_censors=show_censors,
                                                                censor_styles={'ms': 5, 'marker': '+'},
                                                                loc=slice(0., 36.))

            plt.xticks(np.arange(0, 36, step=6))
            plt.ylabel("Survival probability")
            plt.xlabel("Time in months")

            if at_risk: add_at_risk_counts(kmf1, kmf2, labels=[label_1, label_2], ax=subplots[splot_index])

            subplots[splot_index].vlines(x=median1, ymin=0, ymax=y_max1,
                                         colors='black',
                                         linestyle='dotted')
            subplots[splot_index].hlines(y=y_max1, xmin=0, xmax=median1,
                                         colors='black',
                                         linestyle='dotted')
            subplots[splot_index].vlines(x=median2, ymin=0, ymax=y_max2,
                                         colors='black',
                                         linestyle='dotted')
            subplots[splot_index].hlines(y=y_max2, xmin=0, xmax=median2,
                                         colors='black',
                                         linestyle='dotted')
            results = logrank_test(durations_A=T[i1], durations_B=T[i2], event_observed_A=E[i1],
                                   event_observed_B=E[i2])
            p_value = check_p_value(results.p_value)
            plt.sca(subplots[splot_index])
            subplots[splot_index].plot()
            subplots[splot_index].set_title(f"{gene}, Q1")
            at = AnchoredText(f"p-value: {p_value}",
                              loc='center right', prop=dict(size=6), frameon=True)
            at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
            subplots[splot_index].add_artist(at)
            plt.xticks(np.arange(0, 36, step=6))
            plt.ylabel("Survival probability")
            plt.xlabel("Time in months")
            splot_index += 1
            if splot_index == m * n:
                pdf.savefig()
                plt.close(f)
                f, axarr = plt.subplots(m, n, sharex="col", sharey="row")
                plt.subplots_adjust(wspace=0.2)
                subplots = helper(m, n, axarr)
                splot_index = 0
                for s, splot in enumerate(subplots):
                    last_row = (m * n) - s < n + 1
                    first_in_row = s % n == 0
                    if last_row:
                        splot.set_xlabel("Time in months")
                    if first_in_row:
                        splot.set_ylabel("Survival probability")

            # TODO: 50%
            i1, i2, threshold, median_i1, median_i2 = get_groups(df[gene], T, 0.50)
            kmf1 = KaplanMeierFitter()
            label_1 = f'$\leq$ {threshold}'
            kmf1.fit(durations=T[i1], event_observed=E[i1],
                     label=f'{label_1} (n = {len(T[i1])}, median = {float(median_i1)})')
            median1, y_max1 = get_survival_median(kmf1)
            subplots[splot_index] = kmf1.plot_survival_function(ax=subplots[splot_index], show_censors=show_censors,
                                                                censor_styles={'ms': 5, 'marker': '+'},
                                                                loc=slice(0., 36.))

            label_2 = f'> {threshold}'
            kmf2 = KaplanMeierFitter()
            kmf2.fit(durations=T[i2], event_observed=E[i2],
                     label=f'{label_2} (n = {len(T[i2])}, median = {float(median_i2)})')
            median2, y_max2 = get_survival_median(kmf2)
            subplots[splot_index] = kmf2.plot_survival_function(ax=subplots[splot_index], show_censors=show_censors,
                                                                censor_styles={'ms': 5, 'marker': '+'},
                                                                loc=slice(0., 36.))

            plt.xticks(np.arange(0, 36, step=6))
            plt.ylabel("Survival probability")
            plt.xlabel("Time in months")

            if at_risk: add_at_risk_counts(kmf1, kmf2, labels=[label_1, label_2], ax=subplots[splot_index])

            subplots[splot_index].vlines(x=median1, ymin=0, ymax=y_max1,
                                         colors='black',
                                         linestyle='dotted')
            subplots[splot_index].hlines(y=y_max1, xmin=0, xmax=median1,
                                         colors='black',
                                         linestyle='dotted')
            subplots[splot_index].vlines(x=median2, ymin=0, ymax=y_max2,
                                         colors='black',
                                         linestyle='dotted')
            subplots[splot_index].hlines(y=y_max2, xmin=0, xmax=median2,
                                         colors='black',
                                         linestyle='dotted')
            results = logrank_test(durations_A=T[i1], durations_B=T[i2], event_observed_A=E[i1],
                                   event_observed_B=E[i2])
            p_value = check_p_value(results.p_value)
            plt.sca(subplots[splot_index])
            subplots[splot_index].plot()
            subplots[splot_index].set_title(f"{gene}, Q2")
            at = AnchoredText(f"p-value: {p_value}",
                              loc='center right', prop=dict(size=6), frameon=True)
            at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
            subplots[splot_index].add_artist(at)
            plt.xticks(np.arange(0, 36, step=6))
            plt.ylabel("Survival probability")
            plt.xlabel("Time in months")
            # plt.show()
            splot_index += 1
            if splot_index == m * n:
                pdf.savefig()
                plt.close(f)
                f, axarr = plt.subplots(m, n, sharex="col", sharey="row")
                arr_ij = [(x, y) for x, y in np.ndindex(axarr.shape)]
                plt.subplots_adjust(wspace=0.2)
                subplots = [axarr[index] for index in arr_ij]
                splot_index = 0
                for s, splot in enumerate(subplots):
                    splot.set_ylim(0, 1.12)
                    splot.set_xlim(0, 100)
                    last_row = (m * n) - s < n + 1
                    first_in_row = s % n == 0
                    if last_row:
                        splot.set_xlabel("Time in months")
                    if first_in_row:
                        splot.set_ylabel("Survival probability")

            # TODO: 75%
            i1, i2, threshold, median_i1, median_i2 = get_groups(df[gene], T, 0.75)
            kmf1 = KaplanMeierFitter()
            label_1 = f'$\leq$ {threshold}'
            kmf1.fit(durations=T[i1], event_observed=E[i1],
                     label=f'{label_1} (n = {len(T[i1])}, median = {float(median_i1)})')
            median1, y_max1 = get_survival_median(kmf1)
            subplots[splot_index] = kmf1.plot_survival_function(ax=subplots[splot_index], show_censors=show_censors,
                                                                censor_styles={'ms': 5, 'marker': '+'},
                                                                loc=slice(0., 36.))

            label_2 = f'> {threshold}'
            kmf2 = KaplanMeierFitter()
            kmf2.fit(durations=T[i2], event_observed=E[i2],
                     label=f'{label_2} (n = {len(T[i2])}, median = {float(median_i2)})')
            median2, y_max2 = get_survival_median(kmf2)
            subplots[splot_index] = kmf2.plot_survival_function(ax=subplots[splot_index], show_censors=show_censors,
                                                                censor_styles={'ms': 5, 'marker': '+'},
                                                                loc=slice(0., 36.))

            plt.xticks(np.arange(0, 36, step=6))
            plt.ylabel("Survival probability")
            plt.xlabel("Time in months")

            if at_risk: add_at_risk_counts(kmf1, kmf2, labels=[label_1, label_2], ax=subplots[splot_index])

            subplots[splot_index].vlines(x=median1, ymin=0, ymax=y_max1,
                                         colors='black',
                                         linestyle='dotted')
            subplots[splot_index].hlines(y=y_max1, xmin=0, xmax=median1,
                                         colors='black',
                                         linestyle='dotted')
            subplots[splot_index].vlines(x=median2, ymin=0, ymax=y_max2,
                                         colors='black',
                                         linestyle='dotted')
            subplots[splot_index].hlines(y=y_max2, xmin=0, xmax=median2,
                                         colors='black',
                                         linestyle='dotted')
            results = logrank_test(durations_A=T[i1], durations_B=T[i2], event_observed_A=E[i1],
                                   event_observed_B=E[i2])
            p_value = check_p_value(results.p_value)
            plt.sca(subplots[splot_index])
            subplots[splot_index].plot()
            subplots[splot_index].set_title(f"{gene}, Q3")
            at = AnchoredText(f"p-value: {p_value}",
                              loc='center right', prop=dict(size=6), frameon=True)
            at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
            subplots[splot_index].add_artist(at)
            plt.ylabel("Survival probability")
            plt.xlabel("Time in months")
            splot_index += 1
            if splot_index == m * n:
                pdf.savefig()
                plt.close(f)
                f, axarr = plt.subplots(m, n, sharex="col", sharey="row")
                plt.subplots_adjust(wspace=0.2)
                subplots = helper(m, n, axarr)
                splot_index = 0
                for s, splot in enumerate(subplots):
                    last_row = (m * n) - s < n + 1
                    first_in_row = s % n == 0
                    if last_row:
                        splot.set_xlabel("Time in months")
                    if first_in_row:
                        splot.set_ylabel("Survival probability")

        if (len(genes) % 2 != 0) and m == 2:
            pdf.savefig()
            plt.close(f)
        elif (len(genes) % 3 != 0) and m == 3:
            pdf.savefig()
            plt.close(f)

    matplotlib.rcParams.update({"font.size": 12})


def split_multiple_p_values(df_pvalue):
    s = ""
    i = 0
    l0 = list(df_pvalue['level_0'])
    l1 = list(df_pvalue['level_1'])
    for p in df_pvalue['p']:
        s += f'{l0[i]} & {l1[i]} = {check_p_value(p)}'
        if (i + 1) != len(l0):
            s += '\n'
        i += 1
    return s


def plot_location(T, E, vals, workspace_path, at_risk=True, show_title=False):
    matplotlib.rcParams.update({"font.size": 9})
    values = []
    for side in vals:
        if side == "r" or side == "R":
            values.append('Right')
        elif side == "l" or side == "L":
            values.append('Left')
        elif side == "r,l" or side == 'l,r' or side == 'R,L' or side == 'L,R':
            values.append('Left and right')
        elif side == 'bill':
            values.append('bill')
    temp = pd.DataFrame(values, columns=['Events'])
    groups = temp['Events']
    i_r = (groups == 'Right')
    i_l = (groups == 'Left')
    i_lr = (groups == 'Left and right')

    groups_ = list(groups[i_r]) + list(groups[i_l]) + list(groups[i_lr])
    T_ = list(T[i_r]) + list(T[i_l]) + list(T[i_lr])
    E_ = list(E[i_r]) + list(E[i_l]) + list(E[i_lr])

    fig, ax = plt.subplots()  # figsize=(5,5))  #figsize=(7,9))
    kmf1 = KaplanMeierFitter()
    kmf1.fit(durations=T[i_r], event_observed=E[i_r],
             label=f'Right (n = {len(T[i_r])}, median OS = {statistics.median(list(T[i_r]))})')
    ax = kmf1.plot_survival_function(ax=ax, show_censors=True, censor_styles={'ms': 6, 'marker': 's'},
                                     loc=slice(0., 36.), ci_alpha=0.1)

    kmf2 = KaplanMeierFitter()
    kmf2.fit(durations=T[i_l], event_observed=E[i_l],
             label=f'Left (n = {len(T[i_l])}, median OS = {statistics.median(list(T[i_l]))})')
    ax = kmf2.plot_survival_function(ax=ax, show_censors=True, censor_styles={'ms': 6, 'marker': 's'},
                                     loc=slice(0., 36.), ci_alpha=0.1)

    kmf3 = KaplanMeierFitter()
    kmf3.fit(durations=T[i_lr], event_observed=E[i_lr],
             label=f'Left and right (n = {len(T[i_lr])}, median OS = {statistics.median(list(T[i_lr]))})')
    ax = kmf3.plot_survival_function(ax=ax, show_censors=True, censor_styles={'ms': 6, 'marker': 's'},
                                     loc=slice(0., 36.), ci_alpha=0.1)

    plt.xticks(np.arange(0, 36, step=6))
    plt.ylabel("Survival probability")
    plt.xlabel("Time in months")

    if at_risk: add_at_risk_counts(kmf1, kmf2, kmf3, labels=['RIGHT', 'LEFT', 'LEFT AND RIGHT'], ax=ax)

    results2 = pairwise_logrank_test(event_durations=T_, groups=groups_, event_observed=E_)
    results2_df = results2.summary.p.reset_index(drop=False)
    for col in ['level_0', 'level_1']:
        t = []
        for item in results2_df[col]:
            if item == 'Left':
                t.append('L')
            elif item == 'Right':
                t.append('R')
            elif item == 'Left and right':
                t.append('L, R')
            elif item == 'bill':
                t.append('bill')
        results2_df[col] = t
    at = AnchoredText(f"P-VALUES:\n{split_multiple_p_values(results2_df)}",
                      loc='center right', prop=dict(size=12), frameon=True)

    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax.add_artist(at)
    if show_title: plt.title('Tumor Location')
    plt.savefig(os.path.join(workspace_path, 'results', 'kaplan_meier_curves',
                             'Age, gender, KPS, location, subtypes', "location.png"), bbox_inches='tight')
    plt.close(fig)
    matplotlib.rcParams.update({"font.size": 12})


def plot_contact(T, E, groups, workspace_path, at_risk=True, show_title=None):
    matplotlib.rcParams.update({"font.size": 9})
    i_cortex = (groups == 'cortex')
    i_cv = (groups == 'cortex+ventricle')
    i_ventricle = (groups == 'ventricle')
    i_together = i_ventricle + i_cv

    groups_ = list(groups[i_cortex]) + list(groups[i_cv]) + list(groups[i_ventricle])
    T_ = list(T[i_cortex]) + list(T[i_cv]) + list(T[i_ventricle])
    E_ = list(E[i_cortex]) + list(E[i_cv]) + list(E[i_ventricle])

    fig, ax = plt.subplots()  # figsize=(7, 9))
    kmf1 = KaplanMeierFitter()
    kmf1.fit(durations=T[i_cortex], event_observed=E[i_cortex],
             label=f'Cortex (n = {len(T[i_cortex])}, median OS = {statistics.median(list(T[i_cortex]))})')
    ax = kmf1.plot_survival_function(ax=ax, show_censors=True, censor_styles={'ms': 6, 'marker': 's'},
                                     loc=slice(0., 36.), ci_alpha=0.1)

    kmf2 = KaplanMeierFitter()
    kmf2.fit(durations=T[i_ventricle], event_observed=E[i_ventricle],
             label=f'Ventricle (n = {len(T[i_ventricle])}, median OS = {statistics.median(list(T[i_ventricle]))})')
    ax = kmf2.plot_survival_function(ax=ax, show_censors=True, censor_styles={'ms': 6, 'marker': 's'},
                                     loc=slice(0., 36.), ci_alpha=0.1)

    kmf3 = KaplanMeierFitter()
    kmf3.fit(durations=T[i_cv], event_observed=E[i_cv],
             label=f'Cortex + Ventricle (n = {len(T[i_cv])}, median OS = {statistics.median(list(T[i_cv]))})')
    ax = kmf3.plot_survival_function(ax=ax, show_censors=True, censor_styles={'ms': 6, 'marker': 's'},
                                     loc=slice(0., 36.), ci_alpha=0.1)

    plt.xticks(np.arange(0, 36, step=6))
    plt.ylabel("Survival probability")
    plt.xlabel("Time in months")

    if at_risk: add_at_risk_counts(kmf1, kmf2, kmf3, labels=['CORTEX', 'VENTRICLE', 'CORTEX + VENTRICLE'], ax=ax)

    results2 = pairwise_logrank_test(event_durations=T_, groups=groups_, event_observed=E_)
    results2_df = results2.summary.p.reset_index(drop=False)
    for col in ['level_0', 'level_1']:
        t = []
        for item in results2_df[col]:
            if item == 'cortex':
                t.append('C')
            elif item == 'ventricle':
                t.append('V')
            elif item == 'cortex+ventricle':
                t.append('C + V')
            elif item == 'none':
                t.append('none')
        results2_df[col] = t
    at = AnchoredText(f"P-VALUES:\n{split_multiple_p_values(results2_df)}",
                      loc='center right', prop=dict(size=11), frameon=True)

    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax.add_artist(at)
    if show_title: plt.title('Ventricle contact')
    plt.savefig(os.path.join(workspace_path, 'results', 'kaplan_meier_curves',
                             'Age, gender, KPS, location, subtypes', "ventricle_contact.png"), bbox_inches='tight')
    plt.close(fig)

    fig, ax = plt.subplots()  # figsize=(5,5))
    kmf1 = KaplanMeierFitter()
    kmf1.fit(durations=T[i_together], event_observed=E[i_together],
             label=f'yes (n = {len(T[i_together])}, median OS = {statistics.median(list(T[i_together]))})')
    ax = kmf1.plot_survival_function(ax=ax, show_censors=True, censor_styles={'ms': 6, 'marker': 's'},
                                     loc=slice(0., 36.), ci_alpha=0.1)

    kmf2 = KaplanMeierFitter()
    kmf2.fit(durations=T[i_cortex], event_observed=E[i_cortex],
             label=f'no (n = {len(T[i_cortex])}, median OS = {statistics.median(list(T[i_cortex]))})')
    ax = kmf2.plot_survival_function(ax=ax, show_censors=True, censor_styles={'ms': 6, 'marker': 's'},
                                     loc=slice(0., 36.), ci_alpha=0.1)
    plt.xticks(np.arange(0, 36, step=6))
    plt.ylabel("Survival probability")
    plt.xlabel("Time in months")

    if at_risk: add_at_risk_counts(kmf1, kmf2, labels=['VENTRICLE CONTACT - YES', 'VENTRICLE CONTACT - NO'], ax=ax)
    results = logrank_test(durations_A=T[i_together], durations_B=T[i_cortex],
                           event_observed_A=E[i_together], event_observed_B=E[i_cortex])
    p_value = check_p_value(results.p_value)
    at = AnchoredText(f"p-value: {p_value}",
                      loc='center right', prop=dict(size=12), frameon=True)

    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax.add_artist(at)
    if show_title: plt.title('Ventricle contact')
    plt.savefig(os.path.join(workspace_path, 'results', 'kaplan_meier_curves',
                             'Age, gender, KPS, location, subtypes',
                             "ventricle contact - ventricle and cortex_ventricle combined.png"),
                bbox_inches='tight')
    plt.close(fig)

    matplotlib.rcParams.update({"font.size": 12})


def plot_gender(T, E, vals, workspace_path, show_censors=True, at_risk=False, show_title=False):
    matplotlib.rcParams.update({"font.size": 9})
    values = []
    for gender in vals:
        if gender == "F" or gender == "f":
            values.append('f')
        elif gender == "M" or gender == "m":
            values.append('m')

    temp = pd.DataFrame(values, columns=['Events'])
    groups = temp['Events']

    i1 = (groups == 'f')  # Female
    i2 = (groups == 'm')  # Male

    fig, ax = plt.subplots()
    kmf1 = KaplanMeierFitter()
    kmf1.fit(durations=T[i1], event_observed=E[i1],
             label=f'Female (n = {len(T[i1])}, median OS = {statistics.median(list(T[i1]))})')
    median1, y_max1 = get_survival_median(kmf1)
    ax = kmf1.plot_survival_function(ax=ax, show_censors=show_censors, censor_styles={'ms': 6, 'marker': 's'},
                                     loc=slice(0., 36.))

    kmf2 = KaplanMeierFitter()
    kmf2.fit(durations=T[i2], event_observed=E[i2],
             label=f'Male (n = {len(T[i2])}, median OS = {statistics.median(list(T[i2]))})')
    median2, y_max2 = get_survival_median(kmf2)
    ax = kmf2.plot_survival_function(ax=ax, show_censors=True, censor_styles={'ms': 6, 'marker': 's'},
                                     loc=slice(0., 36.))

    plt.xticks(np.arange(0, 36, step=6))
    plt.ylabel("Survival probability")
    plt.xlabel("Time in months")

    if at_risk: add_at_risk_counts(kmf1, kmf2, labels=['Female', 'Male'], ax=ax)
    plt.vlines(x=median1, ymin=0, ymax=y_max1,
               colors='black',
               linestyle='dotted')
    plt.hlines(y=y_max1, xmin=0, xmax=median1,
               colors='black',
               linestyle='dotted')
    plt.vlines(x=median2, ymin=0, ymax=y_max2,
               colors='black',
               linestyle='dotted')
    plt.hlines(y=y_max2, xmin=0, xmax=median2,
               colors='black',
               linestyle='dotted')
    results = logrank_test(durations_A=T[i1], durations_B=T[i2], event_observed_A=E[i1], event_observed_B=E[i2])
    p_value = check_p_value(results.p_value)
    at = AnchoredText(f"p-value: {p_value}",
                      loc='center right', prop=dict(size=12), frameon=True)
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax.add_artist(at)
    if show_title: plt.title('Gender')
    plt.savefig(os.path.join(workspace_path, 'results', 'kaplan_meier_curves', 'Age, gender, KPS, location, subtypes',
                             "gender.png"), bbox_inches='tight')
    plt.close(fig)
    matplotlib.rcParams.update({"font.size": 12})
