# pyensembl_operations.py
import os
import numpy as np
from pyensembl import EnsemblRelease
from scripts.globals import PYENSEMBL_CACHE_DIR

def initialize_genome(grch: int = 37) -> EnsemblRelease:
    """Initialize and return the genome assembly."""
    if grch not in [37, 38]:
        raise ValueError("grch must be either 37 or 38")
    
    os.environ['PYENSEMBL_CACHE_DIR'] = PYENSEMBL_CACHE_DIR
    ens_release = 75 if grch == 37 else 111
    assembly = EnsemblRelease(ens_release)
    assembly.download()
    assembly.index()
    return assembly

def cached_pyensembl_call(locus, assembly, canonical_only, function_name):
    """Call a function from the assembly object with the given locus."""
    chrom, pos = locus.split(':')
    func = getattr(assembly, function_name)
    result = func(chrom, int(pos))
    return result[0] if len(result) > 0 else np.nan
