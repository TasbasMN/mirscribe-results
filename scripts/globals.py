from socket import gethostname

def get_pyensembl_cache_location():
    hostname = gethostname()
    if hostname == "nazo":
        return "/home/nazif/thesis/data"
    elif hostname == "Minerva":
        return "/home/yamak/Code/nazif/data"
    elif hostname == 'nazos-MacBook-Air.local':
        return "/Users/nazo/thesis/data"
    else:
        return "/truba/home/mtasbas/data"


    
PYENSEMBL_CACHE_DIR = get_pyensembl_cache_location()

GRCH37_DIR = "data/fasta/grch37"
MIRNA_COORDS_DIR = "data/mirna_coordinates"
TA_SPS_CSV = "data/ta_sps/ta_sps.csv"
MIRNA_CSV = "data/mirna/mirna.csv"
XGB_MODEL = "misc/models/model_with_no_close_proximity.json"

NUCLEOTIDE_OFFSET = 30


AWK_SCRIPT_PATH = "scripts/rnaduplex_to_csv.awk"
MUTSIG_PROBABILITIES = "data/mutsig_probabilities/probabilities.csv"
PROBABILITIES_LATEST = "data/mutsig_probabilities/latest_probabilities/probabilities_december.csv"
CLINPATH_FILE = "data/clinpath/clinpath.csv"

