# simulate module
import subprocess
import os
import logging
import glob

import numpy as np
import pandas as pd

from mgeasysim import config as cf
from concurrent.futures import ThreadPoolExecutor

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def setup_logging_for_function(func_name):
    """Configure logging for a specific function based on its name."""
    log_filename = os.path.join(cf.OUTPUT, f"{func_name}.log")
    # Create a logger for the specific function
    logger = logging.getLogger(func_name)
    logger.setLevel(logging.INFO)
    
    # Remove any existing handlers to prevent duplication
    if logger.hasHandlers():
        logger.handlers.clear()
    
    # Set up file handler
    file_handler = logging.FileHandler(log_filename)
    file_handler.setLevel(logging.INFO)
    
    # Define log format
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    
    # Add handler to logger
    logger.addHandler(file_handler)
    return logger

# Function to run a command
def run_command(command, logger):
    try:
        result = subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
        logger.info(f"Success: {command}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error: {command}")
        logger.error(e.stderr)

def simulate(simulation_data, N_READS, n_threads, acc2genbank, genome_lengths, genome2file):

    logger = setup_logging_for_function('simulate')

    BASE_PATH = os.path.join(cf.OUTPUT, 'simulations')
    os.makedirs(BASE_PATH, exist_ok=True)
    
    # for each simulation
    total_reads = 0
    for sim in simulation_data['simid'].unique():
        logger.info(f'SIMULATION NUMBER: {sim}')
        
        # make a directory for the simulation
        if not os.path.exists(f'{BASE_PATH}/sim{sim}'):
            os.mkdir(f'{BASE_PATH}/sim{sim}')

        # make subdirs
        if not os.path.exists(f'{BASE_PATH}/sim{sim}/fqs'):
            os.mkdir(f'{BASE_PATH}/sim{sim}/fqs')
            os.mkdir(f'{BASE_PATH}/sim{sim}/dbs')
            os.mkdir(f'{BASE_PATH}/sim{sim}/genomes')


        os.makedirs(f'{BASE_PATH}/sim{sim}/genomes/sylph_db', exist_ok=True)
 
        #
        read_simulations = []
        # for each genome included in the current simulation
        for i, (index, row) in enumerate(simulation_data[simulation_data['simid'].eq(sim)].iterrows()):

            # get n reads and fold coverage 
            n_reads = int(np.ceil(row['abun'] * N_READS))
            fold_coverage = f'{(n_reads * 150 * 2) / genome_lengths.loc[row["index"]]:0.20f}'
            logger.info(str(row['index']))
            logger.info(' '.join([str(i), 'NREADS: ', str(n_reads), ' ABUN:', str(row['abun']), 'FOLD COV', str(fold_coverage)]))
            
            total_reads += n_reads
            file = genome2file.loc[acc2genbank.loc[row['index']]]
            # copy genome to genomes folder for default db
            command = f'cp {file} {BASE_PATH}/sim{sim}/genomes/sylph_db/{row["index"]}_genomic.fna'
            print(command)
            result = subprocess.run(command.split(' '), shell=False, capture_output=True)
            
            # simulate reads
            outfile = f'{BASE_PATH}/sim{sim}/fqs/{row['index']}'
            if not os.path.exists(f'{outfile}_R1.fq.gz'):
                command1 = f'art_illumina -ss HS25 -i {file} -l 150 -f {fold_coverage} -d {row['index']} -m 300 -s 5 -o {outfile}_R -p'
                command2 = f'gzip {outfile}_R*'
                # print(command)
                # result = subprocess.run(command.split(' '), shell=False, capture_output=True)
                read_simulations.append(f'{command1} && {command2}')
                logger.info(f'Added command: {command1} && {command2}')
            else:
                logger.info(f'File already exists: {outfile}_R1.fq.gz')
                
        logger.info(' '.join(['Number of samples to simulate: ', str(len(read_simulations))]))
        # run in parallel 
        if len(read_simulations) > 0:
            with ThreadPoolExecutor(max_workers=n_threads) as executor:
                executor.map(lambda args: run_command(*args), [(cmd, logger) for cmd in read_simulations])
        
        logger.info('Dealing with files...')
        command = f'cat {BASE_PATH}/sim{sim}/fqs/*_R1.fq.gz > {BASE_PATH}/sim{sim}/sim{sim}_R1.fq.gz'
        # Run the command with stderr redirected to subprocess.PIPE
        print(command)
        process = subprocess.Popen(command, 
                                   stdout=subprocess.PIPE, 
                                   stderr=subprocess.PIPE, 
                                   shell=True)
        process.wait()
        stdout, stderr = process.communicate()

        # Check for errors and print stderr if the command failed
        if process.returncode != 0:
            
            print("An error occurred:")
            print(stderr)
            raise ValueError('blah')
        else:
            print("Command succeeded:")

        logger.info(command)
        # p = subprocess.Popen(command.split(' '), shell=False)
        # p.wait()

        command = f'cat {BASE_PATH}/sim{sim}/fqs/*_R2.fq.gz > {BASE_PATH}/sim{sim}/sim{sim}_R2.fq.gz'
        logger.info(command)
        p = subprocess.Popen(command, shell=True)
        p.wait()

# decompletion code
def DecompleteGenome(genome, completeness, outfile, to_file=False):
    # for each sequence
    # get length of sequence to remove
    # pull it out and make two new records
    # save records to file
    with open(genome, 'r') as handle:
        records = [r for r in SeqIO.parse(handle, 'fasta')]
        
    new_records = []
    for r in records:
        # get length of record
        r_len = len(r.seq)
        # get a chunk to remove
        chunk = int(np.floor(r_len * (1 - completeness)))
        # determine start location
        start = np.random.randint(1, r_len - chunk - 1)
        # get chunks
        chunk1 = r.seq[:start]
        new_records.append(SeqRecord(chunk1, id=r.id+'_chunk1'))
        chunk2 = r.seq[start+chunk+1:]
        new_records.append(SeqRecord(chunk2, id=r.id+'_chunk2'))
    if to_file:
        # write to file
        with open(outfile, 'w') as handle:
            SeqIO.write(new_records, handle, format='fasta')
    else:
        return new_records
        
# outfile MUST end in .fa
def ContaminateGenome(genome, contaminant, outfile):
    # add 25% contaminant reads to fasta
    # first decomplete the genome
    genome_records = DecompleteGenome(genome, 0.75, outfile, to_file=False)
    contam_records = DecompleteGenome(contaminant, 0.25, outfile, to_file=False)
    with open(outfile, 'w') as handle:
        SeqIO.write(genome_records + contam_records, handle, format='fasta')
         
# decompletion code
def DecompleteGenome(genome, completeness, outfile, to_file=False):
    # for each sequence
    # get length of sequence to remove
    # pull it out and make two new records
    # save records to file
    with open(genome, 'r') as handle:
        records = [r for r in SeqIO.parse(handle, 'fasta')]
        
    new_records = []
    for r in records:
        # get length of record
        r_len = len(r.seq)
        # get a chunk to remove
        chunk = int(np.floor(r_len * (1 - completeness)))
        # determine start location
        start = np.random.randint(1, r_len - chunk - 1)
        # get chunks
        chunk1 = r.seq[:start]
        new_records.append(SeqRecord(chunk1, id=r.id+'_chunk1'))
        chunk2 = r.seq[start+chunk+1:]
        new_records.append(SeqRecord(chunk2, id=r.id+'_chunk2'))
    if to_file:
        # write to file
        with open(outfile, 'w') as handle:
            SeqIO.write(new_records, handle, format='fasta')
    else:
        return new_records
        
    
def create_sparse_series(index):
    length = len(index)
    # Start with a series of all zeroes
    series = pd.Series(0, index=index)
    
    # Determine the number of ones based on the length of the index
    if length > 1:
        num_ones = max(1, int(length * 0.1)) if length >= 10 else 1
        ones_indices = np.random.choice(index, size=num_ones, replace=False)
        series.loc[ones_indices] = 1
        
    return series

def download_contaminants(logger):
    current_dir = os.path.dirname(__file__)
    contam_path = os.path.join(current_dir, '../files', 'contaminants.txt')
    os.makedirs(os.path.join(cf.OUTPUT, 'contaminants'))
    command1 = [f'datasets download genome accession --inputfile', 
                contam_path,
                '--filename', os.path.join(cf.OUTPUT, 'contaminants', 'genomes_dataset.zip')]

    command2 = ['unzip', '-q', '-o',
                os.path.join(cf.OUTPUT, 'contaminants', 'genomes_dataset.zip'),
                '-d',
                os.path.join(cf.OUTPUT, 'contaminants')]

    command1 = ' '.join(command1)
    logger.info('getting contaminants:')
    logger.info(command1)
    p1 = subprocess.Popen(command1.split(' '), 
                        shell = False,)
    
    p1.wait()

    logger.info('unzipping contaminants:')
    logger.info(' '.join(command2))
    p2 = subprocess.Popen(command2, 
                        shell = False,)
    p2.wait()
        
def generate_alt_databases(simulation_data, logger, genome2file, acc2genbank):
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
                        subprocess.run(f'cp {genome_file} {outfile}', shell=True)
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


def run_sylph(simdata, n_threads, genome2file, acc2genbank, alt_dbs = False, containment=100):

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
        generate_alt_databases(simdata, logger, genome2file, acc2genbank)
        dbs = dbs + alt_dbs
    else:
        dbs = ['sylph_db']

    for sim in simdata['simid'].unique():

        for _db in dbs:
            print(f'Sketching {_db}...')
            logger.info(f'Sketching {_db}...')
            command = f'sylph sketch {BASE_PATH}/sim{sim}/genomes/{_db}/*.fna -o {BASE_PATH}/sim{sim}/dbs/{_db} -t {n_threads} -c {containment}'
            subprocess.run(command, shell=True)
            
            logger.info(f'Profiling with {_db} to output {BASE_PATH}/sim{sim}/sim{sim}_{_db}_profile.tsv')
            command = f'sylph profile {BASE_PATH}/sim{sim}/dbs/{_db}.syldb -1 {BASE_PATH}/sim{sim}/sim{sim}_R1.fq.gz -2 {BASE_PATH}/sim{sim}/sim{sim}_R2.fq.gz -c {containment} -o {BASE_PATH}/sim{sim}/sim{sim}_{_db}_profile.tsv -t {n_threads}'
            subprocess.run(command.split(' '), shell=False)
            print('\n\n')

