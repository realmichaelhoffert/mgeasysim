# simulate module
import warnings
import subprocess
import os
import glob

import pandas as pd
from rapidfuzz import process, fuzz
import numpy as np
from tqdm import tqdm

from mgeasysim import config as cf
from concurrent.futures import ThreadPoolExecutor


# Function to run a command
def run_command(command):
    try:
        result = subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
        print(f"Success: {command}")
    except subprocess.CalledProcessError as e:
        print(f"Error: {command}")
        print(e.stderr)

def simulate(simulation_data, N_READS, n_threads, acc2genbank, genome_lengths, genome2file):
    BASE_PATH = os.path.join(cf.OUTPUT, 'simulations')
    os.makedirs(BASE_PATH, exist_ok=True)
    
    # for each simulation
    total_reads = 0
    for sim in simulation_data['simid'].unique():
        print('SIMULATION NUMBER:', sim)
        # make a directory for the simulation
        if not os.path.exists(f'{BASE_PATH}/sim{sim}'):
            os.mkdir(f'{BASE_PATH}/sim{sim}')

        # make subdirs
        if not os.path.exists(f'{BASE_PATH}/sim{sim}/fqs'):
            os.mkdir(f'{BASE_PATH}/sim{sim}/fqs')
            os.mkdir(f'{BASE_PATH}/sim{sim}/dbs')
            os.mkdir(f'{BASE_PATH}/sim{sim}/genomes')

        #
        read_simulations = []
        for i, (index, row) in enumerate(simulation_data[simulation_data['simid'].eq(sim)].iterrows()):
            print(row)
            n_reads = int(np.ceil(row['abun'] * N_READS))
            fold_coverage = f'{(n_reads * 150 * 2) / genome_lengths.loc[row["index"]]:0.20f}'
            print(i, 'NREADS: ', n_reads, row['index'], ' ABUN:', row['abun'], 'FOLD COV', fold_coverage)
            
            total_reads += n_reads
            file = genome2file.loc[acc2genbank.loc[row['index']]]
            # print()
            print(i, 'NREADS: ', n_reads, row['index'], ' ABUN:', row['abun'])
            # copy genome to genomes file
            command = f'cp {file} {BASE_PATH}/sim{sim}/genomes/{row["index"]}_genomic.fna'
            result = subprocess.run(command.split(' '), shell=False, capture_output=True)
            
            # simulate reads
            outfile = f'{BASE_PATH}/sim{sim}/fqs/{row['index']}'
            if not os.path.exists(f'{outfile}_R1.fq.gz'):
                print('added')
                command1 = f'art_illumina -ss HS25 -i {file} -l 150 -f {fold_coverage} -d {row['index']} -m 300 -s 5 -o {outfile}_R -p'
                command2 = f'gzip {outfile}_R*'
                # print(command)
                # result = subprocess.run(command.split(' '), shell=False, capture_output=True)
                read_simulations.append(f'{command1} && {command2}')
                
        print('Number of samples to simulate: ', len(read_simulations))
        # run in parallel 
        with ThreadPoolExecutor(max_workers=n_threads) as executor:
            executor.map(run_command, read_simulations)
        
        print('Dealing with files...')
        command = f'cat {BASE_PATH}/sim{sim}/fqs/*_R1.fq.gz > {BASE_PATH}/sim{sim}/sim{sim}_R1.fq.gz'
        os.system(command)
        
        command = f'cat {BASE_PATH}/sim{sim}/fqs/*_R2.fq.gz > {BASE_PATH}/sim{sim}/sim{sim}_R2.fq.gz'
        os.system(command)
        
        
        print('Running sylph...')
        ### RUN SYLPH
        # ----------------------------------------------------------------------
        command = f'sylph sketch {BASE_PATH}/sim{sim}/genomes/*.fna -o {BASE_PATH}/sim{sim}/dbs/sylph_db_c100 -t 8 -c 100'
        subprocess.run(command, shell=True)
        
        command = f'sylph profile {BASE_PATH}/sim{sim}/dbs/sylph_db_c100.syldb -1 {BASE_PATH}/sim{sim}/sim{sim}_R1.fq.gz -2 {BASE_PATH}/sim{sim}/sim{sim}_R2.fq.gz -c 100 -o {BASE_PATH}/sim{sim}/sim{sim}_sylph_prof.tsv -t 8'
        subprocess.run(command.split(' '), shell=False)
            
        # clear_output(wait=True)

