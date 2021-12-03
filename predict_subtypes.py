import argparse
import os

from utils.run_subtype_prediction import main as main_


def perform_subtypes_prediction(workspace_path, dataset_path, model_path):
    folder = os.path.join(workspace_path, 'results', 'determine subtypes')
    if not os.path.exists(folder):
        os.mkdir(folder)

    print("Subtype predicition started.")
    main_(workspace_path, dataset_path, model_path)
    print("Subtype predicition is finished.")


def main():
    parser = argparse.ArgumentParser(description='TCGA Glioblastoma Subtypes Prediction')

    parser.add_argument("--workspace_path", help="Path to the visualisation workspace", required=True, action='store')
    parser.add_argument("--dataset_path", help="Path to the new datasest to be used for prediction", required=True, action='store')
    parser.add_argument("--model_path", help="Path to the decision tree model", required=False, action='store')

    args = parser.parse_args()
    perform_subtypes_prediction(args.workspace_path, args.dataset_path, args.model_path)


if __name__ == "__main__":
    main()
