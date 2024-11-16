"""
Module containing utilties for mgeasysim
"""
# simulate module
import subprocess
import os
import logging
import glob

import numpy as np
import pandas as pd

from mgeasysim import config as cf

from rapidfuzz import process, fuzz
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Function to find top match from slist for each item in qlist
def find_top_matches(qlist, slist):
    matches = []
    for query in qlist:
        # Extract the best match from slist for the current query
        match, score, _ = process.extractOne(query, slist)
        matches.append(match)
    return pd.Series(matches, index=qlist, name='top_match')

def get_genome2file():
    # save distribution of genome files
    genome2file = pd.Series()
    genome_folders = glob.glob(os.path.join(cf.OUTPUT, 'ncbi_dataset/data/*/*_genomic.fna'))
    for g_file in genome_folders:
        genome = g_file.split('/')[-2]
        genome2file.loc[genome] = g_file

    return genome2file

def get_genome_lengths():
    # get genome lengths
    gtdb_md = pd.read_csv(cf.GTDB_MD, sep='\t', index_col='accession')
    acc2genbank = gtdb_md['ncbi_genbank_assembly_accession']
    
    genome_lengths = pd.Series()
    genome2file = get_genome2file()
    for gtdb_acc, genbank_acc in acc2genbank.items():
        if genbank_acc in genome2file.index:
            file = genome2file[genbank_acc]
            with open(file, 'r') as handle:
                total_length = np.sum([len(r.seq) for r in SeqIO.parse(handle, 'fasta')])
            genome_lengths.loc[gtdb_acc] = total_length

    return genome_lengths, acc2genbank
    
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

class LastLogHandler(logging.Handler):
    """Custom handler to capture the last logged message."""
    def __init__(self):
        super().__init__()
        self.last_log = None  # Store the last log message

    def emit(self, record):
        self.last_log = self.format(record)  # Format and store the log record

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
    file_handler = logging.FileHandler(log_filename, mode='w')
    file_handler.setLevel(logging.INFO)
    
    # Set up custom handler to capture the last logged message
    last_log_handler = LastLogHandler()
    last_log_handler.setLevel(logging.INFO)
    
    # Define log format
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    last_log_handler.setFormatter(formatter)
    
    # Add handlers to logger
    logger.addHandler(file_handler)
    logger.addHandler(last_log_handler)
    
    # Attach the custom handler to the logger for easy access
    logger.last_log_handler = last_log_handler  # Allow retrieval of last log message
    
    return logger


# Function to run a command
def run_command(command, logger, verbose=False, error_message=None):
    """
    Runs a command using subprocess, logs outputs, and optionally prints them.
    
    Args:
        command (str): The command to run as a string.
        logger (logging.Logger): A logger object to write outputs.
        verbose (bool): If True, prints outputs to the console in addition to logging.
        error_message (str, optional): Custom error message to return on failure.
        
    Returns:
        str: The standard output of the command if successful.
        
    Raises:
        RuntimeError: If the command exits with a non-zero code.
    """
    try:
        # Run the command
        process = subprocess.Popen(
            command,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        stdout, stderr = process.communicate()
        
        # Log outputs
        if stdout:
            logger.info(stdout.strip())
            if verbose:
                print(stdout.strip())
        if stderr:
            logger.error(stderr.strip())
            if verbose:
                print(stderr.strip())
        
        # Check for non-zero exit code
        if process.returncode != 0:
            error_msg = error_message or f"Command failed with exit code {process.returncode}"
            logger.error(error_msg)
            raise RuntimeError(error_msg)
        
        return stdout.strip()
    
    except Exception as e:
        logger.exception("An error occurred while running the command.")
        raise RuntimeError(str(e)) from e
    
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