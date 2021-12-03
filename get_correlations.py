import argparse
import os

from utils.run_correlations import main as corr_main


def perform_correlations(workspace_path):
    correlations = os.path.join(workspace_path, 'results', 'correlations')
    if not os.path.exists(correlations):
        os.mkdir(correlations)

    folders = ['gene combinations', 'gene combinations and tumor location',
               'gene markers combinations', 'mutations and gene markers combinations',
               'mutations and tumor location']
    for folder in folders:
        f = os.path.join(correlations, folder)
        if not os.path.exists(f):
            os.mkdir(f)

    print("Obtaining correlations. It may take a while...")
    corr_main(workspace_path)
    print("Correlations obtained successfully.")


def main():
    parser = argparse.ArgumentParser(description='Correlations')

    parser.add_argument("--workspace_path", help="Path to the visualisation workspace", required=True, action='store')

    args = parser.parse_args()

    perform_correlations(args.workspace_path)


if __name__ == "__main__":
    main()
