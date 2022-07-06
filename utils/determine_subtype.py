import os
import pickle
import warnings

warnings.filterwarnings('always')
import graphviz
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi'] = 300
import numpy as np
import pandas as pd
import plotly.express as px
import seaborn as sns
from dtreeviz.trees import dtreeviz
from sklearn import tree
from sklearn.decomposition import PCA
from sklearn.metrics import make_scorer, roc_auc_score
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score, cross_validate
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.tree import DecisionTreeClassifier, _tree
from sklearn import metrics

from utils.variables import biomarkers_for_subtype_definition as features, subtype_target_names

warnings.filterwarnings("ignore")

workspace_path = "../"


def plot_first_three_components(components, total_explained_variance, labels):
    l = {'0': 'PC 1', '1': 'PC 2', '2': 'PC 3'}
    fig = px.scatter_3d(
        components, x=0, y=1, z=2,
        size=abs(components[:, 3]) if components.shape[1] == 4 else None,
        color=labels,
        title=f'Total Explained Variance: {total_explained_variance:.3f}%',
        labels=l.update({'3': 'PC 4'}) if components.shape[1] == 4 else l
    )
    fig.show()


def get_linear_combinations(pca_comp, index):
    cols = [f'PC {i}' for i in range(1, len(pca_comp) + 1)]
    df = pd.DataFrame(pca_comp.T, columns=cols, index=index)
    f = open(os.path.join(workspace_path, 'results', 'determine subtypes', 'pca_results.txt'), 'w')
    f.write("Linear combinations:\n\n")
    for component in cols:
        s = f"{component} = "
        count = 0
        for item in sorted(list(zip(df[component], index)), reverse=True):
            if count == 0:
                s += f"{float('{:.4f}'.format(item[0]))} * {item[1]}"
            else:
                if item[0] >= 0:
                    s += f" + {float('{:.4f}'.format(item[0]))} * {item[1]}"
                else:
                    s += f" - {abs(float('{:.4f}'.format(item[0])))} * {item[1]}"
            count += 1
        f.write(f"{s}\n")
    f.write(
        '\n......................................................................................................\n\n')
    f.write('Sorted loadings by absolute value for the first three principal components\n\n')
    for i in range(1, 4):
        vals = df[f'PC {i}'].abs().sort_values(ascending=False).round(3).to_string()
        f.write(f'\nPC {i} loadings:\n\n')
        f.write(vals)
        f.write("\n")
    f.write(
        '\n......................................................................................................\n\n')
    f.close()

    fig, ax = plt.subplots(figsize=(10, 6))
    cmap = sns.diverging_palette(230, 0, 90, 60, as_cmap=True)

    temp = np.transpose(pca_comp)
    annot = [["{:.3f}".format(val)
              for col, val in enumerate(column)] for row, column in enumerate(temp)]

    ax = sns.heatmap(df, cmap=cmap, annot=annot, fmt='')
    ax.set_xticklabels(
        ax.get_xticklabels(),
        rotation=45,
        horizontalalignment='right'
    )
    ax.set_yticklabels(
        ax.get_yticklabels(),
        rotation=0
    )

    plt.title("Linear Combination Weights - Loadings")
    plt.savefig(
        os.path.join(workspace_path, 'results', 'determine subtypes',
                     f'PC_loadings_for_{len(pca_comp)}_components.png'),
        bbox_inches='tight')
    plt.close(fig)


def get_correlations_between_variables_and_pca(pca_components, pca_expl_var, index):
    loadings = pca_components.T * np.sqrt(pca_expl_var)
    cols = [f'PC {i}' for i in range(1, len(pca_components) + 1)]
    corr_matrix = pd.DataFrame(loadings, columns=cols, index=index)
    plt.figure(figsize=(10, 6))
    cmap = sns.diverging_palette(230, 0, 90, 60, as_cmap=True)
    temp = corr_matrix.to_numpy()
    annot = [["{:.3f}".format(val)
              for col, val in enumerate(column)] for row, column in enumerate(temp)]
    ax = sns.heatmap(corr_matrix, cmap=cmap, annot=annot, fmt='')
    ax.set_xticklabels(
        ax.get_xticklabels(),
        rotation=45,
        horizontalalignment='right'
    )
    ax.set_yticklabels(
        ax.get_yticklabels(),
        rotation=0
    )
    plt.title("Correlations Between Genes and Principal Components")
    plt.savefig(os.path.join(workspace_path, 'results', 'determine subtypes',
                             f'correlations_between_genes_and_principal_components_{len(cols)}'),
                bbox_inches='tight')
    plt.close()


def get_rules(tree, feature_names, class_names):
    tree_ = tree.tree_
    feature_name = [
        feature_names[i] if i != _tree.TREE_UNDEFINED else "undefined!"
        for i in tree_.feature
    ]

    paths = []
    path = []

    leaf_content = dict()

    def recurse(node, path, paths):

        if tree_.feature[node] != _tree.TREE_UNDEFINED:
            name = feature_name[node]
            threshold = tree_.threshold[node]
            p1, p2 = list(path), list(path)
            p1 += [f"({name} <= {np.round(threshold, 3)})"]
            recurse(tree_.children_left[node], p1, paths)
            p2 += [f"({name} > {np.round(threshold, 3)})"]
            recurse(tree_.children_right[node], p2, paths)
        else:
            path += [(tree_.value[node], tree_.n_node_samples[node])]
            paths += [path]

    recurse(0, path, paths)

    # sort by samples count
    samples_count = [p[-1][1] for p in paths]
    ii = list(np.argsort(samples_count))
    paths = [paths[i] for i in reversed(ii)]

    rules = []
    for path in paths:
        rule = "if "

        for p in path[:-1]:
            if rule != "if ":
                rule += " and "
            rule += str(p)
        rule += " then "
        if class_names is None:
            rule += "response: " + str(np.round(path[-1][0][0][0], 3))
        else:
            classes = path[-1][0][0]
            l = np.argmax(classes)
            leaf_content.update({class_names[l]: [path[-1][1], list(
                zip(class_names, classes))]})  # {leaf class : zip(classes, samples)}
            rule += f"class: {class_names[l]} (proba: {np.round(100.0 * classes[l] / np.sum(classes), 2)}%)"
        rule += f" | based on {path[-1][1]:,} samples"
        rules += [rule]

    return rules, leaf_content


def laplace_smoothing(k, n, c, alpha=1.0):
    """
    Perform Laplace smoothing for each leaf in decision tree

    :parameter
    k : The number of training cases of the class label classified by the leaf
    n : The total number of training cases classified by the leaf
    c : The number of classes.
    alpha : float, default=1.0
        Laplace smoothing parameter.
    :return:
    p_l : Laplace smoothed probability
    """

    p_l = (k + alpha) / (n + c)

    return float("{:.3f}".format(p_l))


def decision_tree(X, y, feature_names, subtypes, depth=None):
    com = sorted(list(set(list(zip(LabelEncoder().fit_transform(y), list(subtypes))))))
    target_names = [i[1] for i in com]

    clf = DecisionTreeClassifier(criterion="entropy", splitter="best", max_depth=depth,
                                 random_state=0)
    kfold = StratifiedKFold(n_splits=10, shuffle=True, random_state=0)
    print('Measure\t| mean \t | std')
    scores = cross_val_score(clf, X, y, scoring='accuracy', cv=kfold, n_jobs=-1)
    print('Accuracy: %.3f (%.3f)' % (np.mean(scores), np.std(scores)))

    scores = cross_val_score(clf, X, y, scoring='precision_macro', cv=kfold, n_jobs=-1)
    print('Precision (macro): %.3f (%.3f)' % (np.mean(scores), np.std(scores)))

    scores = cross_val_score(clf, X, y, scoring='recall_macro', cv=kfold, n_jobs=-1)
    print('Recall (macro): %.3f (%.3f)' % (np.mean(scores), np.std(scores)))

    scores = cross_val_score(clf, X, y, scoring='f1_macro', cv=kfold, n_jobs=-1)
    print('F1-score (macro): %.3f (%.3f)' % (np.mean(scores), np.std(scores)))

    scores = cross_val_score(clf, X, y, scoring='precision_weighted', cv=kfold, n_jobs=-1)
    print('Precision (weighted): %.3f (%.3f)' % (np.mean(scores), np.std(scores)))

    scores = cross_val_score(clf, X, y, scoring='recall_weighted', cv=kfold, n_jobs=-1)
    print('Recall (weighted): %.3f (%.3f)' % (np.mean(scores), np.std(scores)))

    scores = cross_val_score(clf, X, y, scoring='f1_weighted', cv=kfold, n_jobs=-1)
    print('F1-score (weighted): %.3f (%.3f)' % (np.mean(scores), np.std(scores)))

    clf.fit(X, y)
    print('auc macro', roc_auc_score(y, clf.predict_proba(X), multi_class='ovr', average='macro'))
    print('auc weighted', roc_auc_score(y, clf.predict_proba(X), multi_class='ovr', average='weighted'))

    modelname = 'predict_subtypes_dt_model.sav'
    pickle.dump(clf, open(os.path.join(workspace_path, modelname), 'wb'))

    text_representation = tree.export_text(clf, feature_names=feature_names)
    print(text_representation)

    rules, leaf_content = get_rules(clf, feature_names=feature_names, class_names=target_names)
    for r in rules:
        print(r)

    # Perform Laplace Smoothing    # laplace_smoothing(k, n, c, alpha=1.0)
    leaves = dict()
    for leaf_class in leaf_content:
        info = leaf_content[leaf_class]
        c = 4
        n = info[0]
        content = dict()
        content.update({'n': n})
        for t in info[1]:
            k = t[1]
            prob = laplace_smoothing(k, n, c)
            content.update({f'{t[0]}-k': k})
            content.update({f'{t[0]}-prob': prob})

        leaves.update({leaf_class: content})

    leaves_prob_df = pd.DataFrame.from_dict(leaves, orient='index')
    leaves_prob_df = leaves_prob_df.reset_index(drop=False).rename(columns={"index": "Leaf Class"})

    leaves_prob_df.to_csv(os.path.join(workspace_path, 'results',
                                       'decision_tree_leaf_class_probabilities_after_laplace_smoothing.csv'),
                          sep=';', index=False)

    dot_data = tree.export_graphviz(clf, out_file=None,
                                    feature_names=feature_names,
                                    class_names=target_names,
                                    filled=True,
                                    max_depth=depth,
                                    rounded=True)

    # Draw a graph
    graph = graphviz.Source(dot_data, format="png",
                            filename=os.path.join(workspace_path, 'results', 'determine subtypes',
                                                  'decision_tree_rules_graph'))
    graph.view()

    # Draw a decision tree with histograms
    viz = dtreeviz(clf, X, y,
                   target_name="TCGA SUBTYPE",
                   feature_names=feature_names,
                   class_names=target_names,
                   fancy=True  # False = To remove histograms
                   )
    viz.save(os.path.join(workspace_path, 'results', 'determine subtypes', 'decision_tree_rules_histograms.svg'))
    viz.view()


def perform_subtype_definition(data, true_labels, no_components, no_components_interest,
                               plot_exp_var_ratio, plot_pca_3D, depth=2):

    data_scaled = StandardScaler().fit_transform(data)
    pca = PCA(n_components=no_components)
    pca.fit(data_scaled)
    print("Explained variance ratio for 3 components: ", pca.explained_variance_ratio_[:3].sum() * 100)
    print("Explained variance ratio for 4 components: ", pca.explained_variance_ratio_[:4].sum() * 100)


    if plot_exp_var_ratio:
        fig, ax = plt.subplots()
        plt.bar(range(1, pca.n_components_ + 1), pca.explained_variance_ratio_, alpha=0.5, align='center',
                label='individual explained variance ratio')
        plt.step(range(1, pca.n_components_ + 1), np.cumsum(pca.explained_variance_ratio_), where='mid',
                 label='cumulative explained variance ratio')
        plt.ylabel('Explained variance ratio')
        plt.xlabel('Principal components')
        plt.legend(loc='best')
        plt.savefig(os.path.join(workspace_path, 'results', 'determine subtypes', 'explained_variance_ratio.png'),
                    bbox_inches='tight')
        # plt.show()
        plt.close(fig)



    pca = PCA(n_components=no_components_interest)
    pca.fit(data_scaled)
    data_scaled = pca.transform(data_scaled)
    get_linear_combinations(pca.components_, data.columns)
    get_correlations_between_variables_and_pca(pca.components_, pca.explained_variance_, index=data.columns)
    feature_names = [f"PC {i + 1}" for i in range(no_components_interest)]
    decision_tree(X=data_scaled, y=LabelEncoder().fit_transform(true_labels), feature_names=feature_names,
                  subtypes=true_labels, depth=depth)

    if plot_pca_3D:
        total_explained_variance = pca.explained_variance_ratio_[:no_components_interest].sum() * 100
        plot_first_three_components(data_scaled, total_explained_variance, true_labels)


def main(workspace='../'):
    global workspace_path
    workspace_path = workspace

    data = pd.read_csv(os.path.join(workspace_path, 'data', 'gene_expression_dataset.csv'), sep=';')
    data = data.dropna()

    target = data.Subtype
    perform_subtype_definition(data[features], true_labels=target, no_components=len(features),
                               no_components_interest=3, plot_exp_var_ratio=True,
                               plot_pca_3D=True, depth=2)


if __name__ == "__main__":
    main()
