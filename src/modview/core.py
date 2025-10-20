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
        if self.DataFrame.empty:
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
