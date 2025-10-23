import pandas as pd 
import pyranges as pr
from dataclasses import dataclass 
from typing import Optional
import os


@dataclass
class Filterconfig: 
    """"Configuration for filtering RNA modification data"""
    min_coverage: int = 20
    min_modification_frequency: float = 10.0

    def __post_init__(self): 
        """Validate filter parameters"""
        if not 0 <= self.min_modification_frequency <= 100:
            raise ValueError("Modification frequency must be between 0 and 100")
        if self.min_coverage < 0:
            raise ValueError("Coverage must be non-negative")
        
@dataclass
class FilterStatistics:
    """Statistics about filtered data"""
    original_count: int 
    filtered_count: int 
    removed_count: int 
    removal_percentage: float

@dataclass
class BedSample:
    """Represents a single .bed file with RNA modification data"""
    name: str
    condition: Optional[str] = None
    modification: Optional[str] = None
    dataframe: Optional[pd.DataFrame] = None
    filter_stats: Optional[FilterStatistics] = None 

    def validate(self) -> bool: 
        if self.dataframe is None: 
            return False 
        if self.dataframe.empty:
            raise ValueError(f"Sample {self.name} has empty dataframe")
        return True
        
    def to_pyranges(self) -> pr.PyRanges:
        """Convert dataframe to PyRanges object for genomic operations"""
        if self.dataframe is None:
            raise ValueError(f"Sample {self.name} has no dataframe")
        return pr.PyRanges(self.dataframe)

class ModificationPipeline:
    """Main pipeline for loading and processing RNA modification data"""

    # Functions applied to pipeline 
    # ...
    MODKIT_COLUMNS = {
        0: "Chromosome",
        1: "Start",
        2: "End",
        3: "Modified_Base_Code",
        4: "Score",
        5: "Strand",
        6: "Start_dup",  
        7: "End_dup",    
        8: "Color",
        9: "Coverage",
        10: "Percent_Modified",
        11: "Nmod",
        12: "Ncanonical",
        13: "Nother",
        14: "Ndelete",
        15: "Nfail",
        16: "Ndiff",
        17: "Nnocall"
    }

    def __init__(self):
        self.samples: dict[str, BedSample] = {}
        self.reference: Optional[pd.DataFrame] = None
        self.config: Filterconfig = Filterconfig()

    def load_bed_files(self, bedfiles: list[str], conditions: Optional[list[str]] = None, modifications: Optional[list[str]] = None, sep: str = "\t") -> dict[str, BedSample]:
        """
        Loads multiple .bed files containing RNA modifications from modkit.

        bedfiles:list[string] : List with the path to each .bed file
        conditions:list[string] : Conditions to detect from filenames 
        modifications:list[string] : Modification types called in the .bed files
        sep:string, optional: Column seperator (default is tab)

        Returns a dictionary: mapping of filename -> BedSample objects 
        """

        if not bedfiles:
            raise ValueError("No .bed files were provided. Please select at leat one file.")
        
        samples = {}
        failed_files = []

        for path in bedfiles:
            full_path = os.path.abspath(path)
            if not os.path.exists(full_path):
                failed_files.append(path)
                print(f"File not found: {path}")
                continue 

            file_name = os.path.splitext(os.path.basename(path))[0]
            file_name_lower = file_name.lower()
            bedsample = BedSample(name=file_name)

            # Detect conditions 
            if conditions:
                for con in conditions:
                    if con.lower() in file_name_lower:
                        bedsample.condition = con
                        break
            
            # Detect modification
            if modifications:
                for mod in modifications:
                    if mod.lower() in file_name_lower:
                        bedsample.modification = mod
                        break

            try: 
                df = pd.read_csv(path, sep=sep, header=None)

                if df.empty:
                    failed_files.append(path)
                    print(f"File is empty: {path}")
                    continue

                df = self.standardize_columns(df)
                bedsample.dataframe = df

                samples[file_name] = bedsample
                print(f"Succsessfully loaded: {file_name}")

            except pd.errors.EmptyDataError:
                failed_files.append(str(path))
                print(f"File is empty: {path}")
            except Exception as e:
                failed_files.append(str(path))
                print(f"Error reading {path}: {e}")  
        
        if not samples:
            raise ValueError(f"No valid .bed files were loaded. Failed files {failed_files}")
        if failed_files:
            print(f"Successfully loded {len(samples)}/{len(bedfiles)} files. Failed: {failed_files}")
        
        self.samples = samples
        return samples 
    
    def standardize_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Standardize column names for .bed files

        df: pd.DataFrame: Raw dataframe from .bed file

        Returns dataframe with standardized column names
        """

        n_cols = len(df.columns)

        rename_map = {}
        for i in range(n_cols):
            if i in self.MODKIT_COLUMNS:
                rename_map[i] = self.MODKIT_COLUMNS[i]
        df = df.rename(columns=rename_map)

        essential_cols = ["Chromosome", "Start", "End", "Score", "Percent_Modified"]
        available_essential = [col for col in essential_cols if col in df.columns]

        return df[available_essential]
    
    def filter_samples(self, config: Optional[Filterconfig] = None) -> dict[str, BedSample]:
        """
        Applies coverage and modification frequency filters to all samples

        config: Filterconfig: Filter configuration. default if not provided

        returns: dict[str, BedSample]: dictionary of filtered samples with statistics
        """

        if not self.samples:
            raise ValueError("No samples loaded. Call load_bed_files() first")
        
        if config:
            self.config = config

        for name, sample in self.samples.items():
            if not sample.validate():
                print(f"Skipping invalid sample: {name}")
                continue

            df = sample.dataframe
            original_count = len(df)

            # Apply filters
            filtered_df = df[
                (df["Score"] >= self.config.min_coverage) &
                (df["Percent_Modified"] >= self.config.min_modification_frequency)
            ]

            filtered_count = len(filtered_df)
            removed_count = original_count - filtered_count 
            removal_percentage = (removed_count / original_count * 100) if original_count > 0 else 0

            # Store statistics
            sample.filter_stats = FilterStatistics(
                original_count=original_count,
                filtered_count=filtered_count,
                removed_count=removed_count,
                removal_percentage=removal_percentage
            )

            sample.dataframe = filtered_df
            print(f"{name}: {original_count} -> {filtered_count} positions ({removal_percentage:.1f}% removed)")
        
        return self.samples

    def load_reference(self, genome_path: str) -> pd.DataFrame:
        """
        Loads reference genome annotation file

        genome_path: str: Path to reference genome annotation file

        Returns: pd.DataFrame
        """

        if not genome_path: 
            raise ValueError("No reference genome annotation file was provided.")
        if not os.path.exists(genome_path):
            raise FileNotFoundError(f"Reference file not found: {genome_path}")
        
        try:
            df = pd.read_csv(genome_path, sep=",")

            if df.empty:
                raise ValueError(f"Reference file is empty: {genome_path}")

            column_map = {
                "Chromosome/scaffold name": "Chromosome",
                "Gene start (bp)": "Start",
                "Gene end (bp)": "End"
            }
        
            df = df.rename(columns=column_map)

            self.reference = df 
            print(f"Successfully loaded reference: {genome_path}")

            return df

        except Exception as e:
            raise ValueError(f"Error loading reference file: {e}")