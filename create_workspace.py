import argparse
import os


def create_workspace(workspace_path):
    if not os.path.exists(workspace_path):
        print('Directory given as a workspace does not exist. Please create it manually.')
        exit(-1)

    if not os.path.exists(os.path.join(workspace_path, 'results')):
        os.mkdir(os.path.join(workspace_path, 'results'))

    if not os.path.exists(os.path.join(workspace_path, 'data')):
        os.mkdir(os.path.join(workspace_path, 'data'))

    print('Workspace has been created successfully.')


def main():
    # create workspace mora napraviti results i sve ostale podmape i data mapu

    parser = argparse.ArgumentParser(description='Glioblastoma gene expression biomarkers visualisation workspace.')

    parser.add_argument("--workspace_path", help="Path to the visualisation workspace", required=True, action='store')

    args = parser.parse_args()

    create_workspace(args.workspace_path)


if __name__ == "__main__":
    main()