import os
import csv
import sys
import pandas
import pathlib
import logging

logger = logging.getLogger()
logger.addHandler(logging.StreamHandler(sys.stdout))
logger.addHandler(logging.FileHandler("../logs/ids2csv.log"))
logger.setLevel(logging.INFO)


def get_ids(file_path, column_name, delimiter):
    # Reading in file
    df = pandas.read_csv(file_path, delimiter=delimiter)

    # Extract IDS
    ids = df[column_name]

    return ids


def find_files(project_dir, ids):
    # Initialize a dictionary to store the results
    results = {}

    # Loop through each id
    for id in ids:
        dir = pathlib.Path(project_dir)
        files = [str(file) for file in dir.glob(f"**/*{id}*.f*q.gz")]
        files = sorted(files)

        if len(files) == 0:
            logger.warning(f"ID {id}: No reads found.")
            continue

        if len(files) == 1:
            logger.warning("ID {id}: Your read files are unpaired.")
            continue

        if len(files) > 2:
            logger.warning("ID {id}: You have too many files with the same ID.")
            continue

        r1_path = files[0]
        r2_path = files[1]
        # this gives relative file path for full path:  str(files[0].resolve())

        if r1_path and r2_path:
            # Process ID
            results[id] = (r1_path, r2_path)

    return results


def write_to_csv(results, output_filename):
    with open(output_filename, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["ID", "forward", "reverse"])
        for subfolder, (r1_path, r2_path) in results.items():
            writer.writerow([subfolder, r1_path, r2_path])


# Example usage if inside /workflow/scripts:
# python ids2csv.py ../../.test/ ../../config/bold_example_lab_input.tsv "Process ID" "\t"
if __name__ == "__main__":
    if len(sys.argv) != 5:
        print(
            "Usage: python script.py <project_dir> <sample_info_file> <column_name> <file_delimiter>"
        )
        sys.exit(1)

    project_dir = sys.argv[1]
    output_filename = os.getcwd() + "samples_out.csv"
    file_path = sys.argv[2]
    column_name = sys.argv[3]
    delimiter = sys.argv[4]
    ids = get_ids(file_path, column_name, delimiter)
    # print("IDs", ids)
    files_info = find_files(project_dir, ids)
    write_to_csv(files_info, output_filename)
    logger.info(f"CSV file '{output_filename}' created successfully.")
