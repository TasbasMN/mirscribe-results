# globals.py
from socket import gethostname

def get_pyensembl_cache_location():
    hostname = gethostname()
    if hostname == "nazo":
        return "/home/nazif/thesis/data"
    elif hostname == "Minerva":
        return "/home/yamak/Code/nazif/data"
    else:
        return "/truba/home/mtasbas/data"

PYENSEMBL_CACHE_DIR = get_pyensembl_cache_location()
GRCH37_DIR = "data/fasta/grch37"
MUTSIG_PROBABILITIES = "data/mutsig_probabilities/probabilities.csv"
