import os
import pandas as pd

# Assuming data_folder is already defined as in the original code
data_folder = "/content/drive/MyDrive/20.440 Project/GSE75688/Subtype_DEGs/Annotated"

dfs = {}
for filename in os.listdir(data_folder):
    if filename.endswith(".csv"):  # or any other file extension you need to read
        filepath = os.path.join(data_folder, filename)
        try:
            df_name = filename[:-4]  # remove the .csv extension
            dfs[df_name] = pd.read_csv(filepath)
            print(f"Successfully read {filename} into DataFrame '{df_name}'")
        except pd.errors.ParserError:
            print(f"Error parsing {filename}. Skipping.")
        except Exception as e:
            print(f"An error occurred while reading {filename}: {e}")