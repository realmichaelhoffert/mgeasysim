# mgeasysim

Cause simulating MG data should be easy

## Outline

MGEasySim is a pile of python code to streamline the process of generating a simulated metagenome from a desired list of GTDB genomes.

### Example

You're developing a tool to analyze the metagenomic sequencing of bacteria in hydrothermal vents. You've gotten a sample and inferred its composition using [sylph](https://github.com/bluenote-1577/sylph) but you haven't gotten the rest of the samples. Using the list of GTDB genomes present in your sample from sylph, you could use MGEasySim to generate as many simulated communities as you need to develop your metagenomics workflow.

### How it works

We're on version 0.2.0: not much functionality has been added.

1. You provide a newline-delimited list of GTDB species names.
2. A match-list is constructed: each of your accessions is matched to a GTDB accession. 
3. A power-distribution is used to generate abundances for each genome in a community of your desired size. You can also choose to adds strains, where a random set of genomes are included alongside a member of their GTDB species cluster at half or twice the abundance: perhaps helpful to test tools that must distinguish between >95% ANI bacteria.
4. You pick a number of simulations.
5. The set of required genomes for all simulations is downloaded using [NCBI datasets](https://github.com/ncbi/datasets).
6. [ART-illumina](https://www.niehs.nih.gov/research/resources/software/biostatistics/art) is used to generate reads for each genome, which are combined into each "simulated community".
7. A database of the genomes used for each community, along with relevant metadata, are saved.

## Installation

1. Install with mamba:

```bash
mamba create env -f environment.yml
```

2. Create a GTDB folder with the GTDB metadata file in it in a directory of your choosing. Right now, you only need metadata.

```bash
wget https://data.gtdb.ecogenomic.org/releases/latest/bac120_metadata.tsv.gz
```

3. Pip install the package:

```bash
# from repo base
pip install -e .
```

## Usage

### 1. Configuration (`config` module)
The package needs to know where your version of the GTDB database is and where you want to write output files. To run the configuration module:
```bash
usage: mgeasysim config [-h] --gtdb GTDB --output OUTPUT [--threads THREADS] [--verbose VERBOSE]

options:
  -h, --help            show this help message and exit
  --gtdb GTDB, -g GTDB  Location of GTDB database
  --output OUTPUT, -o OUTPUT
                        Location of outfile
  --threads THREADS, -@ THREADS
                        Number of threads
  --verbose VERBOSE, -v VERBOSE
                        Verbosity
```

To set the output path and path to the GTDB database, run the command:

```bash
mgeasysim config -g [path to GTDB folder] -o [output path] -@ [threads]
```
This writes a `config.yaml` file to the `mgeasysim/` folder which is readable by all modules, enabling global access to important variables.
Each time the `config` module is imported, it loads the current configuration file from the package directory `mgeasysim/config.yaml` and prints the contents:

```bash
Using config:
{'database': {'gtdb_loc': '/Users/michaelhoffert/Documents/mgeasysim/GTDB_r220'},
 'info': {'software': 'mgeasysim'},
 'locations': {'matches_path': '/Users/michaelhoffert/Documents/mgeasysim/output_test/matches.tsv.gz',
               'output': '/Users/michaelhoffert/Documents/mgeasysim/output_test'},
 'parameters': {'n_sims': 3,
                'n_species': 15,
                'n_strains': 3,
                'power_a': 0.5,
                'verbose': True},
 'runtime': {'threads': 2}}
```

Unfortunately, this means right now you can only have one active configuration at a time.

### 2. Setting up a community (`community` module)

```bash
usage: mgeasysim community [-h] --taxlist TAXLIST [--n_sims N_SIMS] [--n_species N_SPECIES] [--power_a POWER_A] [--n_strains N_STRAINS]

options:
  -h, --help            show this help message and exit
  --taxlist TAXLIST, -t TAXLIST
                        File to define community
  --n_sims N_SIMS, -n N_SIMS
                        Number of simulations
  --n_species N_SPECIES, -s N_SPECIES
                        Number of species in each simulation
  --power_a POWER_A, -a POWER_A
                        Power distribuion A parameter
  --n_strains N_STRAINS, -x N_STRAINS
                        Number of same-species strains to include in simulation
```

Arguments:

1. Taxlist: a newline-delimited list of the taxa you want in your community. Currently support exactly matching GTDB species.
2. n_sims: the number of unique metagenomes you want to simulate
3. n_species: the number of species you want in each simulated community
4. power_a: the A parameter of the power distribution
5. n_strains: the number of co-occurring strains (sub-species) among the n_species species in the community
6. jupyter: whether the code is being run in a jupyter notebook (not relevant) 

### Simulating communities (`simulate` module)

```bash
usage: mgeasysim simulate [-h] [--n_reads N_READS] [--alt_dbs ALT_DBS]

options:
  -h, --help            show this help message and exit
  --n_reads N_READS, -n N_READS
                        Number of reads per simulation
  --alt_dbs ALT_DBS, -a ALT_DBS
                        whether to simulate alternate genome databases
```

The n_reads item controls the number of reads per simulation (all other qualities are controlled in the 'community' step). `alt_dbs` creates several alternate databases: this package began its life as a method of testing the sylph metagenomic profiler. The alt dbs include:

* full: all strains
* incomplete: some strains have lowered completeness
* primary: only primary strains included
* secondary: only secondary strains included
* neither: neither genome for species with strains included
* contam: some genomes are contaminated with E. coli, P. acnes, or Cutibacterium.