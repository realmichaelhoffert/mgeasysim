import pytest

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import numpy as np

import tempfile
from collections import defaultdict

from mgeasysim.utils import *


# Sample FASTA data
fasta_data = """>seq1
AGCTAGCTAGCTAGCTAG
>seq2
TGCATGCATGCATGCATGC
"""

genome_file = "test_genome.fasta"
output_file = "test_output.fasta"

# Mocking np.random.randint using pytest fixture
# @pytest.fixture
# def mock_randint(monkeypatch):
#     # Control the behavior of np.random.randint
#     monkeypatch.setattr(np.random, 'randint', lambda low, high: 5)  # Always return 5
#     return 5

@pytest.fixture
def test_fasta():
    """
    Pytest fixture that provides a temporary FASTA file for testing.
    Automatically handles cleanup after tests.
    """
    # Create sample sequence with mixed cases and multiple lines
    header = ">test_genome"
    sequence = get_random_sequence(10000)
    
    # Create temporary file
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as temp_file:
        # Write header and sequence to file
        temp_file.write(f"{header}\n")
        # Write sequence in chunks of 60 characters per line
        for i in range(0, len(sequence), 60):
            temp_file.write(f"{sequence[i:i+60]}\n")
    
    # Provide the file path to the test
    yield temp_file.name
    
    # Cleanup after test completes
    os.unlink(temp_file.name)
    
@pytest.fixture
def output_fasta():
    """
    Pytest fixture that provides a temporary path for output FASTA.
    Automatically handles cleanup after tests.
    """
    # Create temporary file path
    temp_file = tempfile.NamedTemporaryFile(delete=False, suffix='_reduced.fasta')
    temp_file.close()
    
    # Provide the file path to the test
    yield temp_file.name
    
    # Cleanup after test completes
    if os.path.exists(temp_file.name):
        os.unlink(temp_file.name)
    

@pytest.fixture
def contaminant_fasta():
    """Fixture to create a test contaminant FASTA"""
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as temp_file:
        temp_file.write(">contaminant\n")
        # Create a sequence with different composition
        # Using a repeated pattern of GGCC
        sequence = get_random_sequence(1000) + "GGCC" * 25 + get_random_sequence(1000)  # 100 bases total
        temp_file.write(f"{sequence}\n")
    
    yield temp_file.name
    os.unlink(temp_file.name)

def seq_line_contents(_str):
    return all(i.lower() in ["a","t","g","c","\n"] for i in _str)

def header_line_contents(_str):
    return '>' in _str

def get_random_sequence(l):
    return ''.join(np.random.choice(["A", "T", "G", "C"], l))

def find_sequence_origin(sequence, genome_seq, contaminant_seq, kmer_size=10):
    """
    Determines the origin of each region in the contaminated sequence
    using k-mer matching.
    
    Returns:
        dict: Counts of bases attributed to each source
    """
    origins = defaultdict(int)
    
    # Convert sequences to sets of kmers for each source
    genome_kmers = {genome_seq[i:i+kmer_size] 
                   for i in range(len(genome_seq)-kmer_size+1)}
    contaminant_kmers = {contaminant_seq[i:i+kmer_size] 
                        for i in range(len(contaminant_seq)-kmer_size+1)}
    
    # Scan the sequence with a sliding window
    for i in range(len(sequence)-kmer_size+1):
        kmer = sequence[i:i+kmer_size]
        
        if kmer in genome_kmers and kmer not in contaminant_kmers:
            origins['genome'] += 1
        elif kmer in contaminant_kmers and kmer not in genome_kmers:
            origins['contaminant'] += 1
        else:
            origins['ambiguous'] += 1
    
    # Convert from kmer counts to base counts
    # The first base of each k-mer represents one unique base
    return dict(origins)

class TestDecompleteGenome:

    def test_files(self, test_fasta, output_fasta):
        """
        Unit test function to verify file is a fasta file.
        """
        print('test')
        # try:
        # Reduce genome completeness by 50%
        DecompleteGenome(test_fasta, 0.50, output_fasta, to_file=True )
        
        # Verify output file exists and contents
        assert os.path.exists(output_fasta), "Output file was not created"
        
        # Read output file
        with open(output_fasta, 'r') as f:
            output_lines = f.readlines()
        
        # Basic validation
        assert all(seq_line_contents(s) or header_line_contents(s) for s in output_lines)

    def test_line_lengths(self, test_fasta, output_fasta):
        completeness = 0.50
        DecompleteGenome(test_fasta, completeness, output_fasta, to_file=True )

        with open(output_fasta, 'r') as handle:
            output_records = [r for r in SeqIO.parse(handle, 'fasta')]

        with open(test_fasta, 'r') as handle:
            input_records = [r for r in SeqIO.parse(handle, 'fasta')]

        assert ('_chunk' in r.id for r in output_records)

        total_len_out = np.sum([len(r.seq) for r in output_records])
        total_len_in = np.sum([len(r.seq) for r in input_records])

        assert total_len_in * completeness == pytest.approx(total_len_out, abs=1)

    def test_genome_reduction_invalid_percentage(test_fasta, output_fasta):
        """Test that invalid reduction percentages raise ValueError"""
        with pytest.raises(ValueError):
            DecompleteGenome(test_fasta, -1, output_fasta, to_file=True )
        
        with pytest.raises(ValueError):
            DecompleteGenome(test_fasta, 1, output_fasta, to_file=True )




    

class TestContaminateGenome:

    def test_contaminate_file(self, test_fasta, contaminant_fasta, output_fasta):
        """
        Test to ensure contaminated file is correct length
        """
        ContaminateGenome(test_fasta, contaminant_fasta, output_fasta, 
                     contamination_percentage=0.25)

        # Read all sequences
        genome_seq = ''.join([str(r.seq) for r in SeqIO.parse(test_fasta, "fasta")])
        contaminant_seq =''.join([str(r.seq) for r in SeqIO.parse(contaminant_fasta, "fasta")])
        contaminated_seq = ''.join([str(r.seq) for r in SeqIO.parse(output_fasta, "fasta")])

        # troubleshooting for test failure
        print([(len(r.seq), r.id) for r in SeqIO.parse(output_fasta, "fasta")])
        print(len(contaminated_seq))
        
        # Verify output length matches input length
        expected_length =((len(genome_seq) * 0.75) + (len(contaminant_seq) * 0.25))
        assert len(contaminated_seq) == pytest.approx(expected_length, abs=5), \
            "Output sequence length should match input genome length"


    def test_contaminate_genome(self, test_fasta, contaminant_fasta, output_fasta):
        """
        Test that ContaminateGenome produces correct proportion of contamination
        """
        cp = 0.25
        # Call the function being tested
        ContaminateGenome(test_fasta, contaminant_fasta, output_fasta, 
                        contamination_percentage=cp)
        
        # Read all sequences
        genome_seq = ''.join([str(r.seq) for r in SeqIO.parse(test_fasta, "fasta")])
        contaminant_seq =''.join([str(r.seq) for r in SeqIO.parse(contaminant_fasta, "fasta")])
        contaminated_seq = ''.join([str(r.seq) for r in SeqIO.parse(output_fasta, "fasta")])
        
        # Analyze sequence composition
        origins = find_sequence_origin(contaminated_seq, 
                                       genome_seq, 
                                       contaminant_seq, kmer_size=21)
        
        
        total_assigned_bases = origins['genome'] + origins['contaminant']
        contaminant_percentage = (origins['contaminant'] / total_assigned_bases)

        # Optional: Print detailed statistics
        print(f"\nSequence origin analysis:")
        print(f"Genome bases: {origins['genome']}")
        print(f"Contaminant bases: {origins['contaminant']}")
        print(f"Ambiguous bases: {origins.get('ambiguous', 0)}")
        print(f"Contaminant percentage: {contaminant_percentage:.3f}%")

        # expected amount of DNA from each
        bases_eg = (len(genome_seq) * (1 - cp))
        bases_ec = len(contaminant_seq) * cp

        expected_genome = bases_eg / len(contaminated_seq)
        expected_contaminant = bases_ec / len(contaminated_seq)
        
        if total_assigned_bases > 0:  # avoid division by zero
            contaminant_percentage = (origins['contaminant'] / total_assigned_bases)
            # Check if contaminant percentage is close to 25%
            assert contaminant_percentage == pytest.approx(expected_contaminant, abs=.005), \
                f"Contaminant percentage ({contaminant_percentage:.3f}%) " \
                f"is not close to target ({expected_contaminant}%)"
            
            og_percentage = (origins['genome'] / total_assigned_bases)
            # Check if contaminant percentage is close to 25%
            assert og_percentage == pytest.approx(expected_genome, abs=.005), \
                f"Contaminant percentage ({contaminant_percentage:.3f}%) " \
                f"is not close to target ({expected_genome}%)"
            
