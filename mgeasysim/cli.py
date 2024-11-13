import argparse
import sys
import os

# Add the package's root directory to the sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import pandas as pd

from mgeasysim import config as cf
from mgeasysim import community, simulate

def main():
    parser = argparse.ArgumentParser(prog="mgeasysim", description="MGEasySim: Microbial Genome Easy Simulator")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # config subcommand
    config_parser = subparsers.add_parser("config", help="Configure database locations")
    config_parser.add_argument("--gtdb", '-g', type=str, help="Location of GTDB database", 
                               required=True,
                               )
    config_parser.add_argument("--output", '-o', type=str, help="Location of outfile",
                               required=True)
    config_parser.add_argument("--threads", '-@', type=int, help="Number of threads",
                               default=1, )
    
    # community subcommand
    community_parser = subparsers.add_parser("community", help="Set up a community")
    community_parser.add_argument('--taxlist', '-t', required=True, type=str,
                                  help='File to define community')
    community_parser.add_argument('--n_sims', '-n', default=5, type=int,
                                  help='Number of simulations')
    community_parser.add_argument('--n_species', '-s', default=50, type=int,
                                  help='Number of species in each simulation')
    community_parser.add_argument('--power_a', '-a', default=0.5, type=float,
                                  help='Power distribuion A parameter')
    community_parser.add_argument('--n_strains', '-x', default=2, type=int,
                                  help='Number of same-species strains to include in simulation')

    # Simulate subcommand
    simulate_parser = subparsers.add_parser("simulate", help="Run a genome simulation")
    
    simulate_parser.add_argument("--n_reads", '-n', type=int, help="Number of reads per simulation")
    simulate_parser.add_argument('--alt_dbs', '-a', type=bool, help='whether to simulate alternate genome databases')

    args = parser.parse_args()
    
    if args.command == 'config':
        cf.config._clear()
        cf.config.set('database', 'gtdb_loc', os.path.abspath(args.gtdb))
        cf.config.set('locations', 'output', os.path.abspath(args.output))
        cf.config.set('runtime', 'threads', args.threads)

    elif args.command == "community":
        
        os.makedirs(os.path.abspath(cf.OUTPUT), exist_ok=True)
        matches = community.get_matching_gtdb(args.taxlist)
        genbanks = list(matches['top_match_genbank'].dropna().unique()) + list(matches['alt_genbank'].dropna().unique())
        community.download_genomes(genbanks)
        community.rename_files(genbanks)
        matches = community.add_mashdist(matches)

        cf.config.set('locations', 'matches_path', os.path.join(cf.OUTPUT, 'matches.tsv.gz'))
        matches.to_csv(os.path.join(cf.OUTPUT, 'matches.tsv.gz'), sep='\t', compression='gzip')

        cf.config.set('parameters', 'n_sims', args.n_sims)
        cf.config.set('parameters', 'n_species', args.n_species)
        cf.config.set('parameters', 'power_a', args.power_a)
        cf.config.set('parameters', 'n_strains', args.n_strains)
        
        # get genome lengths and mapping of GTDB accessions to genbanks
        genome_lengths, acc2genbank = community.get_genome_lengths()
        # get genome (genbank) to file mapping from output dir
        genome2file = community.get_genome2file()

        simdata = community.generate_simulations(matches, 
                                       n_sims=args.n_sims, 
                                       n_species=args.n_species, 
                                       power_a=args.power_a, 
                                       n_strains=args.n_strains
                                       )
        
        simdata.to_csv(os.path.join(cf.OUTPUT, 'simulation_data.tsv.gz'), sep='\t')
    
    elif args.command == "simulate":

        # load simdata
        simdata = pd.read_csv(os.path.join(cf.OUTPUT, 'simulation_data.tsv.gz'), sep='\t', index_col=0)

        if os.path.exists(os.path.join(cf.OUTPUT, 'genome_lengths.pkl')):

            genome_lengths = pd.read_pickle(os.path.join(cf.OUTPUT, 'genome_lengths.pkl'))
            acc2genbank = pd.read_pickle(os.path.join(cf.OUTPUT, 'acc2genbank.pkl'))
            genome2file = pd.read_pickle(os.path.join(cf.OUTPUT, 'genome2file.pkl'))

        else:
            # get genome lengths and mapping of GTDB accessions to genbanks
            genome_lengths, acc2genbank = community.get_genome_lengths()
            # get genome (genbank) to file mapping from output dir
            genome2file = community.get_genome2file()

            genome_lengths.to_pickle(os.path.join(cf.OUTPUT, 'genome_lengths.pkl'))
            acc2genbank.to_pickle(os.path.join(cf.OUTPUT, 'acc2genbank.pkl'))
            genome2file.to_pickle(os.path.join(cf.OUTPUT, 'genome2file.pkl'))

        # construct simulated communites
        simulate.simulate(simdata, 
                 N_READS=args.n_reads, 
                 n_threads=cf.config.get('runtime', 'threads'),
                 acc2genbank=acc2genbank, 
                 genome_lengths=genome_lengths, 
                 genome2file=genome2file)
        
        # run sylph with dbs
        simulate.run_sylph(simdata, 
                           n_threads=cf.config.get('runtime', 'threads'),
                           genome2file=genome2file, 
                           acc2genbank=acc2genbank,
                           alt_dbs=args.alt_dbs
                        )

if __name__ == "__main__":
    main()
