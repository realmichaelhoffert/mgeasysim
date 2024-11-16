# simulate module
import subprocess
import os
import logging
import glob

import numpy as np
import pandas as pd

from mgeasysim import config as cf
from mgeasysim.utils import *

from concurrent.futures import ThreadPoolExecutor

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
    

def simulate(simulation_data, N_READS, n_threads, acc2genbank, genome_lengths, genome2file, verbose=True):

    logger = setup_logging_for_function('simulate')
    logger.info(f'Establishing base dir at {cf.OUTPUT} - simulations')

    BASE_PATH = os.path.join(cf.OUTPUT, 'simulations')
    os.makedirs(BASE_PATH, exist_ok=True)
    
    # for each simulation
    total_reads = 0
    for sim in simulation_data['simid'].unique():
        logger.info('-'*50 + f'\nSIMULATION NUMBER: {sim}')
        
        # make a directory for the simulation
        logger.info('Creating directories')
        if not os.path.exists(f'{BASE_PATH}/sim{sim}'):
            os.mkdir(f'{BASE_PATH}/sim{sim}')

        # make subdirs
        
        if not os.path.exists(f'{BASE_PATH}/sim{sim}/fqs'):
            os.mkdir(f'{BASE_PATH}/sim{sim}/fqs')
            os.mkdir(f'{BASE_PATH}/sim{sim}/dbs')
            os.mkdir(f'{BASE_PATH}/sim{sim}/genomes')
        else:
            logger.info('Directories already exist, continuing')


        os.makedirs(f'{BASE_PATH}/sim{sim}/genomes/sylph_db', exist_ok=True)
 
        
        read_simulations = []
        already_exists = 0
        # for each genome included in the current simulation
        for i, (index, row) in enumerate(simulation_data[simulation_data['simid'].eq(sim)].iterrows()):

            # get n reads and fold coverage 
            n_reads = int(np.ceil(row['abun'] * N_READS))
            fold_coverage = f'{(n_reads * 150 * 2) / genome_lengths.loc[row["index"]]:0.20f}'
            logger.info(str(row['index']))
            logger.info('\t\n'.join([str(i), 'NREADS: ', str(n_reads), ' ABUN:', str(row['abun']), 'FOLD COV', str(fold_coverage)]))
            
            total_reads += n_reads
            file = genome2file.loc[acc2genbank.loc[row['index']]]

            # copy genome to genomes folder for default db
            command = f'cp {file} {BASE_PATH}/sim{sim}/genomes/sylph_db/{row["index"]}_genomic.fna'
            run_command(command, 
                        logger, 
                        verbose=verbose, 
                        error_message=f'Copy of genome {row["index"]} to new db folder failed')
            # result = subprocess.run(command.split(' '), shell=False, capture_output=True)
            
            # simulate reads
            outfile = f'{BASE_PATH}/sim{sim}/fqs/{row['index']}'
            if not os.path.exists(f'{outfile}_R1.fq.gz'):
                command1 = f'art_illumina -ss HS25 -i {file} -l 150 -f {fold_coverage} -d {row['index']} -m 300 -s 5 -o {outfile}_R -p'
                command2 = f'gzip {outfile}_R*'
                read_simulations.append(f'{command1} && {command2}')
                logger.info(f'Added command: {command1} && {command2}')
            else:
                logger.info(f'File already exists: {outfile}_R1.fq.gz')
                already_exists += 1
                
        logger.info(' '.join(['Number of samples to simulate: ', str(len(read_simulations))]))
        logger.info(f'Already exists: {already_exists}')

        # run in parallel 
        if len(read_simulations) > 0:
            em1 = 'Art illumina / gzip command failed for file: '
            with ThreadPoolExecutor(max_workers=n_threads) as executor:
                # execute with custom error message
                executor.map(lambda args: run_command(*args), [(cmd, 
                                                                logger, 
                                                                verbose, 
                                                                em1 + cmd.split(' -i ')[-1].split(' -l 150')[0]) for cmd in read_simulations])
        
        logger.info('Dealing with files...')
        
        command = f'cat {BASE_PATH}/sim{sim}/fqs/*_R1.fq.gz > {BASE_PATH}/sim{sim}/sim{sim}_R1.fq.gz'
        logger.info(command)
        run_command(command, logger, verbose, f"Generating {sim} R1 failed")

        command = f'cat {BASE_PATH}/sim{sim}/fqs/*_R2.fq.gz > {BASE_PATH}/sim{sim}/sim{sim}_R2.fq.gz'
        logger.info(command)
        run_command(command, logger, verbose, f"Generating {sim} R2 failed")
        # process = subprocess.Popen(command, 
        #                            stdout=subprocess.PIPE, 
        #                            stderr=subprocess.PIPE, 
        #                            shell=True)
        # process.wait()
        # tdout, stderr = process.communicate()

        # Check for errors and print stderr if the command failed
        # if process.returncode != 0:
            
        #     print("An error occurred:")
        #     print(stderr)
        #     raise ValueError('blah')
        # else:
        #     print("Command succeeded:")

        # logger.info(command)
        # p = subprocess.Popen(command.split(' '), shell=False)
        # p.wait()        

def download_contaminants(logger, verbose=True):
    
    current_dir = os.path.dirname(__file__)
    contam_path = os.path.join(current_dir, '../files', 'contaminants.txt')
    logger.info('Making contaminant directory: ' + os.path.join(cf.OUTPUT, 'contaminants'))
    os.makedirs(os.path.join(cf.OUTPUT, 'contaminants'), exist_ok=True)
    
    command1 = [f'datasets download genome accession --inputfile', 
                contam_path,
                '--filename', os.path.join(cf.OUTPUT, 'contaminants', 'genomes_dataset.zip')]

    command2 = ['unzip', '-q', '-o',
                os.path.join(cf.OUTPUT, 'contaminants', 'genomes_dataset.zip'),
                '-d',
                os.path.join(cf.OUTPUT, 'contaminants')]

    command1 = ' '.join(command1)
    logger.info('Getting contaminants:')
    logger.info(command1)
    run_command(command1, logger, verbose, f"Retrieving contaminants from NCBI failed")

    command2 = ' '.join(command2)
    logger.info('unzipping contaminants:')
    logger.info(command2)
    run_command(command2, logger, verbose, f"Unzipping contaminants NCBI dataset failed")
        
def generate_alt_databases(simulation_data, logger, genome2file, acc2genbank, verbose=True):
    
    logger.info('Generating alternate databases...')
    alt_dbs = ['both_strains', 
                'primary_strain', 
                'secondary_strain', 
                'neither_strain', 
                'contam', 
                'incomp']
    BASE_PATH = os.path.join(cf.OUTPUT, 'simulations')

    for sim in simulation_data['simid'].unique():
        logger.info(f'For {sim}:')

        # generate 
        test_sim_data = simulation_data[simulation_data['simid'].eq(sim)]
        genome_dbs = pd.DataFrame(index=test_sim_data['index'].unique(), columns=alt_dbs)

        # genome_dbs is a dataframe that tells us what to do with each genome for each incomplete db
        for col in genome_dbs.columns:
            if col == 'both_strains':
                genome_dbs[col] = 1
            elif col == 'primary_strain':
                genome_dbs[col] = test_sim_data.set_index('index').reindex(genome_dbs.index)['strain_present'].apply(lambda x: 1 if x < 2 else 0)
            elif col == 'secondary_strain':
                genome_dbs[col] = test_sim_data.set_index('index').reindex(genome_dbs.index)['strain_present'].apply(lambda x: 1 if x != 1 else 0)
            elif col == 'neither_strain':
                genome_dbs[col] = test_sim_data.set_index('index').reindex(genome_dbs.index)['strain_present'].apply(lambda x: 1 if x == 0 else 0)
            elif col == 'incomp':
                genome_dbs[col] = pd.Series(index=genome_dbs.index, data=np.random.beta(a=10, b=2, size=len(genome_dbs.index)))
            elif col == 'contam':
                genome_dbs[col] = create_sparse_series(index=genome_dbs.index)
            else:
                ''
        
        genome_dbs.to_csv(f'{BASE_PATH}/sim{sim}/dbs/genome_dbs.csv', sep='\t')
        
        # for each database
        for col in genome_dbs:
            count = 0
            if not os.path.exists(f'{BASE_PATH}/sim{sim}/genomes/{col}'):
                os.mkdir(f'{BASE_PATH}/sim{sim}/genomes/{col}')
            
            # for each genome
            for genome, val in genome_dbs[col].items():
                
                # locate the file
                genome_file = genome2file.loc[acc2genbank.loc[genome]]
                outfile = f'{BASE_PATH}/sim{sim}/genomes/{col}/{genome}_{col}.fna'
                
                if col != 'contam':
                    if val > 0:
                        count += 1
                        
                    if val < 1:
                        DecompleteGenome(genome_file, 0.75, outfile, to_file=True)
                    else:
                        run_command(f'cp {genome_file} {outfile}', logger, verbose, f"Copying {genome_file} failed")
                else:
                    count += 1
                    # simulate contamination
                    if val == 1:

                        if not os.path.exists(os.path.join(cf.OUTPUT, 'contaminants')):
                            download_contaminants(logger)
                        contam_path = os.path.join(cf.OUTPUT, 'contaminants')    
                        contaminant_files = glob.glob(f'{contam_path}/ncbi_dataset/data/*/*.fna')
                        ContaminateGenome(genome_file, np.random.choice(contaminant_files), outfile)
                    else:
                        subprocess.run(f'cp {genome_file} {outfile}', shell=True)
                        
            logger.info(f' Simulation {sim}: {count} genomes written for db {col}')


def run_sylph(simdata, n_threads, genome2file, acc2genbank, alt_dbs = False, containment=100, verbose=True):

    BASE_PATH = os.path.join(cf.OUTPUT, 'simulations')
    logger = setup_logging_for_function('sylph')

    logger.info('Running sylph...')
    ### RUN SYLPH
    # ----------------------------------------------------------------------
    if alt_dbs:
        dbs = ['sylph_db']
        # set list of dbs
        # code to generate alt dbs

        alt_dbs = ['both_strains', 
                'primary_strain', 
                'secondary_strain', 
                'neither_strain', 
                'contam', 
                'incomp']
        generate_alt_databases(simdata, logger, genome2file, acc2genbank, verbose)
        dbs = dbs + alt_dbs
    else:
        dbs = ['sylph_db']

    for sim in simdata['simid'].unique():

        for _db in dbs:
            logger.info(f'Sketching {_db}...')
            command = f'sylph sketch {BASE_PATH}/sim{sim}/genomes/{_db}/*.fna -o {BASE_PATH}/sim{sim}/dbs/{_db} -t {n_threads} -c {containment}'
            run_command(command, logger, verbose, f'Generating sylph db for {sim} and {_db} failed')
            
            logger.info(f'Profiling with {_db} to output {BASE_PATH}/sim{sim}/sim{sim}_{_db}_profile.tsv')
            command = f'sylph profile {BASE_PATH}/sim{sim}/dbs/{_db}.syldb -1 {BASE_PATH}/sim{sim}/sim{sim}_R1.fq.gz -2 {BASE_PATH}/sim{sim}/sim{sim}_R2.fq.gz -c {containment} -o {BASE_PATH}/sim{sim}/sim{sim}_{_db}_profile.tsv -t {n_threads}'
            run_command(command, logger, verbose, f'Generating sylph profile for {sim} and {_db} failed')
            logger.info('\n')

