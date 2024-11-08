import warnings
from mgeasysim.config import config
import pandas as pd
from rapidfuzz import process, fuzz
import glob

GTDB_BASE = config.get('database', 'path')
if GTDB_BASE is not None:
    global GTDB_MD
    global OUTPUT

    GTDB_MD = glob.glob(f'{GTDB_BASE}*metadata*.tsv')[0]
    OUTPUT = config.get('locations', 'output')



# Function to find top match from slist for each item in qlist
def find_top_matches(qlist, slist):
    matches = []
    for query in qlist:
        # Extract the best match from slist for the current query
        match, score, _ = process.extractOne(query, slist)
        matches.append(match)
    return pd.Series(matches, index=qlist, name='top_match')

def get_matching_gtdb(taxfile, search_col='species'):
    """_summary_
    taxfile: a file containing your list of desired taxa for the simulation
    search_col: the column to search for matches: species is a special key 
    that will search on the last field in the "Gtdb_taxonomy" column
    """
    gtdb_md = pd.read_csv(GTDB_MD, sep='\t', index_col='accession')
    gtdb_md_rep = gtdb_md[gtdb_md.gtdb_representative.eq('t')]
    acc2genbank = gtdb_md.set_index('accession')['ncbi_genbank_assembly_accession']

    if search_col == 'species':
        gtdb_md['species'] = gtdb_md.gtdb_taxonomy.apply(lambda x: x.split(';')[-1])

    with open(taxfile, 'r') as handle:
        qstrs = handle.readlines()

    if gtdb_md.shape[0] != len(gtdb_md[search_col].unique()):
        warnings.Warn('Warning: search column is not unique', UserWarning)

    scol2acc = pd.Series(index=gtdb_md[search_col], data=gtdb_md.index)
                         
    matches = pd.merge(find_top_matches(qstrs, gtdb_md[search_col].values).rename_axis(index='query'),
                       scol2acc)
    matches[''] = matches.map(scol2acc)

    return matches, match_accs, acc2genbank

def get_strains(acc_list):
    """For a list of GTDB accs, get a random extra genome from that species cluster
    """
    gtdb_md = pd.read_csv(GTDB_MD, sep='\t', index_col='accession')
    alts = []

# def download_genomes(genbanks):


# tax_matches, tax_accs, a2g = get_matching_gtdb(taxfile)










