# Toolkit for the Visualisation of the Gliblastoma Data

## 1. Requirements

Toolkit is developed in Python programming language using ``` python 3.8 ``` interpreter. All needed packages can be found in ``` requirements.txt ``` file, and we list the main needed:
- ``` python 3.8 ```
- ``` pandas == 1.2.3 ```
- ``` numpy == 1.19.4 ```
- ``` lifelines == 0.25.11 ```
- ``` matplotlib == 3.4.1 ```
- ``` seaborn == 0.11.1 ```


## 2. Create workspace
```console
python create_workspace.py --workspace_path workspace-dir-path
```

## 3. Data preparation


## 4. Correlations
```console
python get_correlations.py --workspace_path workspace-dir-path
```

## 5. Kaplan-Meier survival curves
```console
python get_survival_prediction.py --workspace_path workspace-dir-path
```

## 6. TCGA glioblastoma subtype definition and prediciton
```console
python predict_subtypes.py --workspace_path workspace-dir-path --dataset-path dataset --model-path model
```
