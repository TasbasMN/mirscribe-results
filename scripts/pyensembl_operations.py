from scripts.globals import PYENSEMBL_CACHE_DIR
from functools import lru_cache
import numpy as np

def import_pyensembl(grch):
    if grch not in [37, 38]:
        raise ValueError("grch must be either 37 or 38")

    ens_release = 75 if grch == 37 else 111
    
    import os

    from pyensembl import EnsemblRelease
    os.environ['PYENSEMBL_CACHE_DIR'] = PYENSEMBL_CACHE_DIR
    assembly = EnsemblRelease(ens_release)
    assembly.download()
    assembly.index()
    
    return assembly


def generate_transcript_id_and_gene_name_columns(df, assembly):
    
    df['transcript_id'] = df.apply(lambda x: assembly.transcript_ids_at_locus(x['chr'], x['pos']), axis=1)
    # df["transcript_name"] = df.apply(lambda x: assembly.transcript_names_at_locus(x['chr'], x['pos']), axis=1)
    df["gene_name"] = df.apply(lambda x: assembly.gene_names_at_locus(x['chr'], x['pos']), axis=1)
    
    

@lru_cache(maxsize=None)
def cached_pyensembl_call(locus, assembly, canonical_only, function_name):
    """
    Call a function from the assembly object with the given locus and parameters.

    Args:
        locus (str): A string in the format "chr:pos" representing the chromosome and position.
        assembly (object): An object containing functions to be called.
        canonical_only (bool): If True, return only the first result. If False, return all results.
        function_name (str): The name of the function to call from the assembly object.

    Returns:
        str: The result of the function call, formatted as a string.
            If the result is an empty list, returns np.nan.
            If canonical_only is True, returns the first element of the result list.
            If canonical_only is False, returns a comma-separated string of all elements in the result list.

    Example:
        >>> cached_pyensembl_call('chr1:12345', assembly_obj, True, 'gene_ids_at_locus')
        'ENSGxxxxxxxxx'
        >>> cached_pyensembl_call('chr2:67890', assembly_obj, False, 'transcript_ids_at_locus')
        'ENSTxxxxxxxxx,ENSTxxxxxxxxx,ENSTxxxxxxxxx'
    """
    chrom, pos = locus.split(':')
    
    # Get the function from the assembly object using its name
    func = getattr(assembly, function_name)
    
    result = func(chrom, int(pos))
    
    if len(result) == 0:
        return np.nan
    elif canonical_only:
        return result[0]
    else:
        return ','.join(result)




def get_gene_names(df, assembly, canonical_only=True):
    
    unique_loci = df['locus'].unique()

    gene_name_dict = {
        locus: cached_pyensembl_call(locus, assembly, canonical_only, 'gene_names_at_locus')
        for locus in unique_loci
    }
    
    df['gene_name'] = df['locus'].map(gene_name_dict)

    return df

def get_transcript_ids(df, assembly, canonical_only=False):

    # Get unique locus values
    unique_loci = df['locus'].unique()

    transcript_id_dict = {
        locus: cached_pyensembl_call(locus, assembly, canonical_only, "transcript_ids_at_locus")
        for locus in unique_loci
    }
    
    if canonical_only:
        df['transcript_id'] = df['locus'].map(transcript_id_dict)
    
    else:
        # Create a new DataFrame with expanded transcript IDs
        expanded_df = df['locus'].map(transcript_id_dict).str.split(',', expand=True).stack().reset_index(level=1, drop=True).to_frame('transcript_id')
        
        # Merge the expanded DataFrame with the original DataFrame
        df = df.join(expanded_df)
   

    return df


def get_gene_biotypes(df, assembly):

    # Get unique locus values
    unique_genes = df['gene_id'].unique()

    biotype_dict = {
        gene: 'NA' if gene == 'NA' else assembly.gene_by_id(gene).biotype
        for gene in unique_genes
    }
    # Create a new column 'biotype' by mapping the 'locus' column to the biotype dictionary
    df['biotype'] = df['gene_id'].map(biotype_dict)
    
    return df


