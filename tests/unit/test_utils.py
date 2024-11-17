import pytest
from io import StringIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import numpy as np
import tempfile

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
    
def seq_line_contents(_str):
    return all(i.lower() in ["a","t","g","c","\n"] for i in _str)

def header_line_contents(_str):
    return '>' in _str

def get_random_sequence(l):
    return ''.join(np.random.choice(["A", "T", "G", "C"], l))

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


    # @pytest.mark.parametrize("reduction_percentage", [0, 25, 75, 100])
    # def test_various_perc(self, ):


    # # Test for DecompleteGenome with to_file=False (no file output)
    # def test_decomplete_genome(self, create_temp_fasta, create_temp_out):
    #     # Simulate file handling with StringIO

        
    #     # Run the function with 80% completeness
    #     new_records = DecompleteGenome(create_temp_fasta(), completeness=0.8, outfile=output_file, to_file=False)

    #     # Check that the records returned have the correct size
    #     for record in new_records:
    #         chunk_len = int(np.floor(len(record.seq) * (1 - 0.8)))  # 20% removal
    #         assert len(record.seq) == chunk_len

    # # Test for DecompleteGenome with to_file=True (writing to file)
    # def test_decomplete_genome_write_to_file(self, mock_randint):
    #     # Simulate file handling with StringIO
    #     mock_input = StringIO(fasta_data)
    #     mock_output = StringIO()
        
    #     # Run the function with 80% completeness and to_file=True
    #     DecompleteGenome(genome_file, completeness=0.8, outfile=output_file, to_file=True)

    #     # Check that the output contains the expected data
    #     mock_output.seek(0)  # Rewind to start of the StringIO buffer
    #     output_records = list(SeqIO.parse(mock_output, "fasta"))
        
    #     # Ensure 4 records are created (2 chunks per original record)
    #     assert len(output_records) == 4

    #     # Check that each chunk has the expected length
    #     for record in output_records:
    #         chunk_len = int(np.floor(len(record.seq) * (1 - 0.8)))
    #         assert len(record.seq) == chunk_len

    # # Test for DecompleteGenome with invalid completeness value
    # def test_decomplete_genome_invalid_completeness(self):
    #     # Test for invalid completeness value (should raise an error)
    #     with pytest.raises(ValueError):
    #         DecompleteGenome(genome_file, completeness=1.5, outfile=output_file, to_file=False)

    #     with pytest.raises(ValueError):
    #         DecompleteGenome(genome_file, completeness=0, outfile=output_file, to_file=False)

    # # Test for DecompleteGenome with to_file=False and verifying the returned records
    # def test_decomplete_genome_no_file_output(self, mock_randint):
    #     # Simulate file handling with StringIO
    #     mock_input = StringIO(fasta_data)
    #     mock_output = StringIO()

    #     # Run the function with to_file=False
    #     new_records = DecompleteGenome(genome_file, completeness=0.8, outfile=output_file, to_file=False)

    #     # Check that we get the expected number of records (2 chunks per original record)
    #     assert len(new_records) == 4

    #     # Check the length of each chunk
    #     for record in new_records:
    #         chunk_len = int(np.floor(len(record.seq) * (1 - 0.8)))
    #         assert len(record.seq) == chunk_len
