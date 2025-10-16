import pandas as pd 
import pyranges as pr 
import os


class BedSample:
    def __init__(self, name : str, condition=None, modification=None, dataframe=None):
        self.name = name
        self.condition = condition
        self.modification = modification
        self.dataframe = dataframe


def load_bed_files(bedfiles: list[str], conditions=None, modifications=None, sep="\t") -> dict[str, BedSample]:
    """
    Loads multiple .bed files containing RNA modifications called via modkit and returns a dictionary of names mapped to a BedSample.

    bedfiles:list[string] : list with the path to each .bed file
    conditions:list[string] : conditions for sample proscessing
    modifications:list[string] : modifications called in the .bed file
    sep:string, optional: Column seperator (default is tab)

    returns a dictionary: mapping of filename -> BedSample 
    """

    if not bedfiles:
        raise ValueError("No .bed files were provided. Please select at least one file.")

    samples = {}

    for path in bedfiles:
        file_name = os.path.splitext(os.path.basename(path))[0]
        file_name_lower = file_name.lower()
        bedsample = BedSample(name=file_name)

        # Detect condition
        if conditions:
            for con in conditions:
                con_lower = str(con).lower()
                if con_lower in file_name_lower:
                    bedsample.condition = con
                    break
        
        # Detect modification
        if modifications:
            for mod in modifications: 
                mod_lower = str(mod).lower()
                if mod_lower in file_name_lower:
                    bedsample.modification = mod
                    break

        # Try loading the dataframe
        try: 
            df = pd.read_csv(path, sep=sep, header=None)
            bedsample.dataframe = df 
        except FileNotFoundError:
            print(f"WARNING: File not found: {path}")
            continue
        except pd.errors.EmptyDataError:
            print(f"WARNING: File is empty: {path}")
            continue
        except Exception as e:
            print(f"WARNING: Error reading path: {path}")
            continue

        # Maybe use something different than file_name
        samples[file_name] = bedsample

    return samples 

def filter_bedfiles(samples: dict[str, BedSample], modification_frequency=10, coverage=20) -> dict[str, BedSample]:
    """
    Applies filter to modification frequency and coverage of multiple .bed files

    samples: dict[str, BedSample]: dictionary of loaded .bed files 
    modification_frequency: Discarding positions with a modification frequency < 10 (default)
    coverage: Discarding positions with a coverage < 20 (default)

    returns a dictionary: mapping of filename -> filtered BedSample
    """

    if not samples: 
        raise ValueError("No files were provided for filtering.")
    
    # Defining necessary columns 
    rename_map = {0: "Chromosome", 1: "Start", 2:"End", 3: "Modified Base Code and Motif", 4: "Score",
                     5: "Strand", 6: "X", 7: "X", 8: "X", 9: "X", 10: "Percent Modified",
                     11: "X", 12: "X", 13: "X", 14: "X", 15: "X", 16: "X", 17:"X"}
    
    for key in samples: 
        df = samples[key].dataframe
        df.rename(columns=rename_map, inplace=True)
        df.drop(columns=[col for col in df if col == "X"], inplace=True)

        filtered_df = df[(df["Score"] >= coverage) & (df["Percent Modified"] >= modification_frequency)]
        samples[key].dataframe = filtered_df

    return samples
    




