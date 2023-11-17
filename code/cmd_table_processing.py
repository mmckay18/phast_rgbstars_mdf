import argparse
import os
import re

import pandas as pd
from tqdm import tqdm


def main(filename):
    data_frames = []
    # Initialize variables
    table_start_pattern = (
        r"# Zini     MH  (.+)"  # Regular expression pattern to identify table start
    )

    # Open the file for reading
    table_num = 0
    with open(filename, "r") as file:
        for line in tqdm(file, desc="Processing file", unit="line"):
            line = file.readline()
            match = re.match(table_start_pattern, line)
            if match:
                table_num += 1
                stripped_line = line.strip()
                stripped_line = stripped_line.replace("\t", " ")
                words = stripped_line.split()
                space_separated_line = " ".join(words[1:])
                comma_separated_line = space_separated_line.replace(" ", ",")
                list_of_colname = comma_separated_line.split(",")

                # Make new dataframe
                # Make new dataframe if this is a new table
                if table_num > len(data_frames):
                    print(f"{table_num} New dataframe created")
                    df = pd.DataFrame(columns=list_of_colname)
                    data_frames.append(df)
                else:
                    df = data_frames[table_num]

            elif line[0] == "#":
                pass

            else:
                stripped_line = line.strip()
                stripped_line = stripped_line.replace("\t", " ")
                words = stripped_line.split()
                space_separated_line = " ".join(words)
                comma_separated_line = space_separated_line.replace(" ", ",")
                list_of_data = comma_separated_line.split(",")
                list_of_data = [float(i) for i in list_of_data]

                # Add data to dataframe
                df.loc[len(df)] = list_of_data


 

    os.makedirs(
        f"/Users/mmckay/phd_projects/analysis_routine/DATA/CMD37_csvs/{filename.split('/')[-1].split('.')[0]}",
        exist_ok=True,
    )
    for i, df in enumerate(data_frames):
        df.to_csv(f"/Users/mmckay/phd_projects/analysis_routine/DATA/CMD37_csvs/{filename.split('/')[-1].split('.')[0]}/{i}_{filename.split('/')[-1].split('.')[0]}.csv", index=False)

    # Merge isochrones together
    merged_df = pd.concat(data_frames)
    merged_df.to_csv(
        f"/Users/mmckay/phd_projects/analysis_routine/DATA/CMD37_csvs/{filename.split('/')[-1].split('.')[0]}/merged_{filename.split('/')[-1].split('.')[0]}.csv",
        index=False,
    )
    print("Merged Dataframe Complete")

    
    



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", help="path to the input file")
    args = parser.parse_args()
    main(args.filename)
