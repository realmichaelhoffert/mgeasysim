"""
Module to construct a community given a list of taxon queries
"""

import warnings
import subprocess
import os
import glob

import pandas as pd

import numpy as np
from tqdm import tqdm
from Bio import SeqIO

from mgeasysim import config as cf
from mgeasysim.utils import *

verboseprint = lambda v, _str: print(_str) if v else None 

def get_matching_gtdb(taxfile, logger, search_col='species', verbose=True):
    """_summary_
    taxfile: a file containing your list of desired taxa for the simulation
    search_col: the column to search for matches: species is a special key 
    that will search on the last field in the "Gtdb_taxonomy" column
    """
    # load the GTDB md
    logger.info('Loading GTDB md...')
    verboseprint(verbose, 'Loading GTDB md...')
    gtdb_md = pd.read_csv(cf.GTDB_MD, sep='\t', index_col='accession')
    rep2genomes = gtdb_md.groupby('gtdb_genome_representative').apply(lambda x: list(x.index.unique()))
    
    # only consider representative genomes
    acc2genbank = gtdb_md['ncbi_genbank_assembly_accession']

    # get and save taxon counts
    taxon_count = gtdb_md.groupby('gtdb_taxonomy').count()['ambiguous_bases']
    gtdb_md['gtdb_taxonomy'].map(taxon_count).to_pickle(os.path.join(cf.OUTPUT, 'cluster_sizes.pkl'))

    gtdb_md = gtdb_md[gtdb_md.gtdb_representative.eq('t')]
    

    # special case for species-level search
    if search_col == 'species':
        gtdb_md['species'] = gtdb_md.gtdb_taxonomy.apply(lambda x: x.split(';')[-1])
    verboseprint(verbose, 'Reading queries...')
    
    # open file of queries
    with open(taxfile, 'r') as handle:
        qstrs = [s.strip() for s in handle.readlines()]

    # warn that the column being matched is not unique (query may have multiple matches)
    if gtdb_md.shape[0] != len(gtdb_md[search_col].unique()):
        warnings.warn('Warning: search column is not unique', UserWarning)

    # get the accession of each item in the search column
    # again, search_col should be unique
    scol2acc = pd.Series(index=gtdb_md[search_col], data=gtdb_md.index)

    # make df with query, match, match_acc columns
    logger.info('Finding matches...')
    verboseprint(verbose, 'Finding matches...')
    matches = find_top_matches(qstrs, gtdb_md[search_col].values)

    matches = pd.merge(matches.rename_axis(index='qstr').reset_index(),
                       scol2acc.rename('top_match_accession'), 
                       how='left', left_on='top_match', right_index=True)

    # get an alternate genome and the genbank id for each genome
    matches['top_match_genbank'] = matches['top_match_accession'].map(acc2genbank)
    matches['top_match_alt'] = matches['top_match_accession'].apply(lambda x: np.random.choice(rep2genomes.loc[x]))
    matches['alt_genbank'] = matches['top_match_alt'].map(acc2genbank)
    
    return matches

def download_genomes(genbanks, logger, verbose=True):

    verboseprint(verbose, 'Writing genome list')
    logger.info('Writing genome list')
    with open(os.path.join(cf.OUTPUT, 'genbanklist.txt'), 'w') as handle:
        handle.write('\n'.join(genbanks))
    
    verboseprint(verbose, 'Downloading')
    command1 = [f'datasets download genome accession --inputfile', 
            os.path.join(cf.OUTPUT, 'genbanklist.txt'),
                '--filename', os.path.join(cf.OUTPUT, 'genomes_dataset.zip')]
    
    command2 = ['unzip', '-q', '-o',
                os.path.join(cf.OUTPUT, 'genomes_dataset.zip'),
                '-d',
                os.path.join(cf.OUTPUT, '')]
    
    command1 = ' '.join(command1)
    logger.info(command1)
    run_command(command1, logger, verbose, 'Downloading NCBI dataset failed')

    command2 = ' '.join(command2)
    logger.info(command2)
    run_command(command2, logger, verbose, 'Unzipping NCBI dataset failed')



def rename_files(genbanks, logger, verbose):
    """
    Function to rename files downloaded from NCBI, which have long filenames
    """
    genome2file = get_genome2file()

    # create copies of files for multithreading read simulation
    for genome in genbanks:
        s = genome2file[genome]
        new_filename = '/'.join(s.split('/')[:-1]) + '/' + genome + '_genomic.fna'
        command = ' '.join(['cp', s, new_filename])
        if not os.path.exists(new_filename):
            run_command(command, logger, verbose, f'{command} Copying {s} failed.')

def add_mashdist(matches, logger, verbose=True):
    """
    Function to add a column with ANI between each genome
    and its alternate strain genome use MASH
    matches: the matching dataframe
    returns:
    matches with additional column alt_mashdist
    """
    matches['alt_mashdist'] = 'none'
    verboseprint(verbose, 'Adding mashdists to matches file')
    logger.info('Adding mashdists to matches file')
    for index, row in tqdm(matches.dropna().iterrows()):
        file1 = os.path.join(cf.OUTPUT, 
                            'ncbi_dataset/data', 
                            row['top_match_genbank'],
                            row['top_match_genbank'] + '_genomic.fna')
        file2 = os.path.join(cf.OUTPUT, 
                            'ncbi_dataset/data', 
                            row['alt_genbank'],
                            row['alt_genbank'] + '_genomic.fna')
        # command = f'mash {file1} {file2}'
        command = ' '.join(['mash', 'dist', file1, file2])
        logger.info(command)
        result = subprocess.run(command, capture_output=True, text=True, shell=True)
        output = result.stdout
        d = float(output.split()[2])
        matches.loc[index, 'alt_mashdist'] = 1 - d

    return matches

def distribution(num_bacteria, a = 0.05):
    # Sample from log-normal distribution
    samples = np.random.power(a, num_bacteria)
    
    # Normalize to sum to one
    relative_abundances = samples / samples.sum()
    
    return relative_abundances

def generate_abundances(num_species, min_abundance=1e-8, max_abundance=0.5, exponent=0.5):
    # Generate random values from a Pareto distribution
    pareto_samples = np.random.pareto(exponent, size=num_species)
    
    # Normalize the samples to sum to 1
    scaled_samples = pareto_samples / np.sum(pareto_samples)
    
    # Clip values to enforce the minimum and maximum abundance
    scaled_samples = np.clip(scaled_samples, min_abundance, max_abundance)
    
    # Re-normalize to make the abundances sum to 1 after clipping
    final_abundances = scaled_samples / np.sum(scaled_samples)
    
    return final_abundances

def generate_simulations(logger, matches, n_sims, n_species, power_a, n_strains):
    simulations = []
    logger.info('Making simulations')
    for i in range(n_sims):
        # simulate data
        abundances = generate_abundances(n_species, exponent=power_a)
        simulation = pd.DataFrame(index=matches['top_match_accession'].sample(n_species, replace=False), 
                                                                      data=abundances, 
                                                                      columns=['abun'])
        # column to indicate which strains have duplicates
        simulation['strain_present'] = 0

        taxon_counts = pd.read_pickle(os.path.join(cf.OUTPUT, 'cluster_sizes.pkl'))
        simulation['n_genomes'] = taxon_counts.reindex(simulation.index)

        # only get strains from GTDB genomes with more than 1 genome in the cluster 
        strains_possible = simulation[simulation['n_genomes'] > 1]
        
        # if there's only 1 thing with a genome in the cluster, only 1 strain is possible
        alts = strains_possible.sample(np.min([n_strains, len(strains_possible)]), replace=False)
        logger.info(f'Number of alternate straiins: {len(alts)}')

        # make data
        alt_data = pd.DataFrame(index=[matches[matches.top_match_accession.eq(i)]['top_match_alt'].values[0] for i in alts.index], 
                                columns=['abun', 'strain_present'])
        
        # make each alternate species either 2x or 1/2 as abundant
        
        for j, (index, row) in enumerate(alt_data.iterrows()):
            
            # print(alternate)
            alternate = index
            original = alts.index[j]
            val = simulation.loc[original, 'abun'] / np.random.choice([0.5, 2])
            alt_data.loc[alternate, 'abun'] = val
            alt_data.loc[alternate, 'strain_present'] = 2
            simulation.loc[original, 'strain_present'] = 1
            logger.info('Original: {original}, alt: {alt}, abun: {val}')

        # add to whole data
        simulation = pd.concat([simulation, alt_data])
        simulation['simid'] = i
        simulation['abun'] = simulation['abun'] / simulation['abun'].sum()
        simulations.append(simulation.reset_index())
        
    return pd.concat(simulations)










