# Toolkit for the Visualisation of the Gliblastoma Data

This toolkit is created for the purpose of master thesis at Faculty of Information and Computer Science, University of Ljubljana by Selma Halilčević.

## 1. Requirements

Toolkit is developed in Python programming language using ``` python 3.8 ``` interpreter. All needed packages can be found in ``` requirements.txt ``` file, and we list the main needed:
- ``` python 3.8 ```
- ``` pandas == 1.2.3 ```
- ``` numpy == 1.19.4 ```
- ``` lifelines == 0.25.11 ```
- ``` matplotlib == 3.4.1 ```
- ``` seaborn == 0.11.1 ```


## 2. Create workspace

After cloning of the repository go to the repository folder and execute the following command with workspace directory path stated in ``` workspace-dir-path ``` argument. 
```console
python create_workspace.py --workspace_path workspace-dir-path
```
Workspace must be created in the cloned repository folder. Hence workspace path must lead to the cloned repository folder. This command creates visualisation workspace for glioblastoma data. It creates data and results folders needed for the script execution.

## 3. Data preparation

After the workspace is created, user must prepare datasets to be used by the toolkit. Datasets are stored in created ``` data ``` folder and three datasets are needed:
- ``` clinical_dataset.csv ```
- ``` pathohistological_dataset.csv ```
- ``` gene_expressions_dataset.csv ```

Every dataset must contain feature ``` SampleName ```, which is a unique identifier for glioblastoma samples and it is a feature through which we merge datasets in toolkit. Every dataset uses semicolon, ``` ;```, as a delimeter.

### Clinical dataset
Clinical dataset contains neuroclinical information of glioblastoma patients. Variables that should be in the clinical dataset are listed below in a table.
| Variable                      | Type              | Additional information                                                                                                    |
|-------------------------------|-------------------|---------------------------------------------------------------------------------------------------------------------------|
| SampleName                    | categorical       | Unique identifier.                                                                                                        |
| Gender (m/f)                  | categorical       |                                                                                                                           |
| Date of birth                 | categorical       |                                                                                                                           |
| Date of Diagnosis (first MRI) | categorical       |                                                                                                                           |
| death_date                    | categorical       |                                                                                                                           |
| Dead                          | numerical; binary | Dead (1-dead, 0-censored) at the time <br>the study was finished.                                                         |
| Overall survival              | numerical; int    | Time in months between the date of the first <br>diagnosis using magnetic resonance imagining (MRI) <br>and a death date. |
| Age at time of 1st diagnosis  | numerical; int    |                                                                                                                           |
| Pre-operative KPS             | numerical; int    | Karnofsky performance status (KPS)                                                                                        |
| Post operative KPS            | numerical; int    | Karnofsky performance status (KPS)                                                                                        |
| Tumor side                    | categorical       |                                                                                                                           |
| Ventricle contact             | categorical       |                                                                                                                           |

### Pathohistological dataset
Pathohistological dataset contains data about mutation tests performed on glioblastoma tissue samples. Variables that should be in the clinical dataset are listed below in a table. Each variable that represents one of the mutations or molecular alterations, can be identified, not identified or not tested.
| Variable                  |
|---------------------------|
| Molecular characteristics |
| IDH1                      |
| IDH R132H                 |
| EGFR                      |
| BRAF                      |
| TERT                      |
| ATRX                      |
| p53                       |
| Co-deletion 1p/19q        |
| MGMT methylation          |

### Gene expressions dataset
Gene expressions dataset consists of a variable SampleName to reference samples, 40 gene expression variables and corresponding molecular or TCGA subtype, variable Subtype.

|                             | 40 gene variables and values of Subtype variable:                                                                                                                                                                                                                                         |
|-----------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Gene expression biomarkers: | TUFM, TRPV1, TRIM28, THBS1, TGFB1, STMN4, SOX10, S100A4, P2RX7, KCNF1,<br>GPR55, ERBB3, DAB2, CST7, COL1A2, COL1A, CNR1, CD9, CCR5, CCL5, <br>ACSBG1, VIM, STAT3, SOX2, SNAI1, PROM1, OLIG2, OCT4, NOTCH, NF-KB,<br>ID1, GFAP, CHI3L1, CEBPA, CDH1, CD44, CD15, CCR3, BETATUBULIN, ALYREF |
| Subtype:                    | CL, MES, MIX, PN                                                                                                                                                                                                                                                                          |

## 4. Correlations
To obtain correlations between molecular characterics, tumour location, gene expression biomarkers, and glioblastoma subtypes run the following command:

```console
python get_correlations.py --workspace_path workspace-dir-path
```

This creates heatmaps in ``` results/correlations ``` folder. Correlation heatmaps are divided to different subfolders based on correlation results among molecular characteristics, biomarkers or glioblastoma subtypes.

## 5. Kaplan-Meier survival curves
We use Python library ``` lifelines ``` to perform Kaplan-Meier survival analysis. To obtain Kaplan-Meier survival curves run the following command.
```console
python get_survival_prediction.py --workspace_path workspace-dir-path
```
After the execution of above command Kaplan-Meier survival curves for age at the time of first diagnosis, gender, Karnofsky performance status, subtypes and biomarkers can be found in ``` results/kaplan_meier_curves ```.

## 6. TCGA glioblastoma subtype definition and prediciton

This toolkit comes with already defined decision tree model for glioblastoma subtype prediction. Model predicts four glioblastoma subtypes, classical, mesenchymal, proneural and MIX. MIX subtype is defined when all 15 gene expression biomarkers needed for definition of other three subtypes are expressed. To predict subtype on new dataset use following command, where the needed arguments are path to the working workspace, path to the new dataset and path to the model, which is optional.

```console
python predict_subtypes.py --workspace_path workspace-dir-path --dataset-path dataset --model-path model
```
After the execution of above command, ``` .csv ``` file will be created with results in the ``` results ``` folder.

## Repository structure
    .
    ├── data                                # Contains datasets needed for the toolikit. Created after the creation of working workspace.
    ├── results                             # Contains all results produced by the toolikit. Created after the creation of working workspace.
    ├── utils                               # Contains python files for performing visualisation of correlations, Kaplan-Meier survival and subtype prediction
    │   ├── correlations.py                 # Contains functions needed for execution of correlations
    │   ├── determine_subtype.py            # Determines subtypes. Can be run from IDE.
    │   ├── kaplan_meier.py                 # Contains functions needed for execution of Kaplan-Meier survival analysis
    │   ├── run_correlations.py             # Runs  all correlations. Can be run from IDE.
    │   ├── run_kaplan_meier.py             # Runs all Kaplan-Meier survival curves. Can be run from IDE.
    │   ├── run_subtype_prediction.py       # Runs glioblastoma subtype definition
    │   └── variables.py                    # Contains all variables and gene expression biomarkers combinations.
    ├── LICENCE                   
    ├── README.md                    
    ├── create_workspace.py                 # Creates working workspace.
    ├── get_correlations.py                 # Runs correlations from terminal.
    ├── get_survival_prediction.py          # Runs Kaplan-Meier survival analysis from terminal.
    ├── predict_subtypes.py                 # Runs glioblastoma subtype prediction from terminal.
    ├── predict_subtypes_dt_model.sav       # Decision tree model for glioblastoma subtype prediction.
    └── requirements.txt                    # Requirements file with all needed packages to be installed.
    
## Reference

If you use this toolkit, please cite as follows:

``` Reference for citation to be added later. ```
