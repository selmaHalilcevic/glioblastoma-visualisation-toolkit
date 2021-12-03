import os
import pickle
import warnings

import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from utils.variables import biomarkers_for_subtype_definition as features, subtype_target_names

warnings.filterwarnings("ignore")


def perform_pca(data, n_components=3):
    data_scaled = StandardScaler().fit_transform(data)
    pca = PCA(n_components)
    pca.fit(data_scaled)
    data_scaled = pca.transform(data_scaled)
    return data_scaled


def predict_subtypes(data, model):
    data_scaled = perform_pca(data, n_components=3)
    predicted = model.predict(data_scaled)
    labels = []
    for n in predicted:
        labels.append(subtype_target_names.get(n))
    return labels


def check_extension(path, extension):
    if path[-3:] == extension:
        return path
    else:
        return f"{path}.{extension}"


def main(workspace='../', data_path='',
         model_path='./predict_subtypes_dt_model.sav'):
    global workspace_path
    workspace_path = workspace

    df = pd.read_csv(check_extension(data_path, 'csv'), sep=';')
    model = pickle.load(open(check_extension(model_path, 'sav'), 'rb'))
    predicted_labels = predict_subtypes(df[features], model)
    df['Predicted Subtype'] = predicted_labels
    df.to_csv(os.path.join(workspace_path, 'results', 'determine subtypes', 'predicted subtypes.csv'), sep=';',
              index=False)


if __name__ == "__main__":
    main()
