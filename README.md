# mgeasysim

Cause simulating MG data should be easy

## Installation

1. Install with mamba:

```bash
mamba create env -f environment.yaml
```

2. Create a GTDB folder with the GTDB metadata file in it in a directory of your choosing.

## Usage

### 1. Configuration (`config` module)
The package needs to know where your version of the GTDB database is and where you want to write output files. To run the configuration module:

```bash
mgeasysim config -g [path to GTDB folder] -o [output path] -@ [threads]
```

This writes a `config.yaml` file to the `mgeasysim/` folder which is readable by all modules, enabling global access to important variables.

### 2. Setting up a community (`community` module)

```bash
usage: mgeasysim community [-h] --taxlist TAXLIST [--n_sims N_SIMS]
                           [--n_species N_SPECIES] [--power_a POWER_A]
                           [--n_strains N_STRAINS] [--jupyter JUPYTER]

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
                        Number of same-species strains to include in
                        simulation
  --jupyter JUPYTER, -j JUPYTER

```

Arguments:

1. Taxlist: a newline-delimited list of the taxa you want in your community. Currently support exactly matching GTDB species.
2. n_sims: the number of unique metagenomes you want to simulate
3. n_species: the number of species you want in each simulated community
4. power_a: the A parameter of the power distribution
5. n_strains: the number of co-occurring strains (sub-species) among the n_species species in the community
6. jupyter: whether the code is being run in a jupyter notebook (not relevant) 

### Simulating communities (`simulate` module)

