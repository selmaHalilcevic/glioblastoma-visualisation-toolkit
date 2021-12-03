import argparse
import os

from utils.run_kaplan_meier import main as km_main


def perform_survival_prediction(workspace_path):
    kaplan_meier_curves = os.path.join(workspace_path, 'results', 'kaplan_meier_curves')
    if not os.path.exists(kaplan_meier_curves):
        os.mkdir(kaplan_meier_curves)

    folders = ['Age, gender, KPS, location, subtypes', 'All genes', 'combinations']
    for folder in folders:
        f = os.path.join(kaplan_meier_curves, folder)
        if not os.path.exists(f):
            os.mkdir(f)

    print("Performing Kaplan-Meier survival prediction. It may take a while...")
    km_main(workspace_path)
    print("Kaplan-Meier survival prediction obtained.")


def main():
    parser = argparse.ArgumentParser(description='Survival Prediction')

    parser.add_argument("--workspace_path", help="Path to the visualisation workspace", required=True, action='store')

    args = parser.parse_args()

    perform_survival_prediction(args.workspace_path)


if __name__ == "__main__":
    main()
