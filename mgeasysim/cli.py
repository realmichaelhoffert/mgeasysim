import argparse
import sys
import os

# Add the package's root directory to the sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import pandas as pd

from mgeasysim import config as cf
from mgeasysim import community, simulate, utils

def main():
    parser = argparse.ArgumentParser(prog="mgeasysim", description="MGEasySim: Microbial Genome Easy Simulator")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # config subcommand
    config_parser = subparsers.add_parser("config", help="Configure database locations")
    config_parser.add_argument("--gtdb", '-g', type=str, help="Location of GTDB database", 
                               required=True,
                               )
    config_parser.add_argument("--config", '-c', type=str, help="Location of configuration",
                               required=True)
    config_parser.add_argument("--threads", '-@', type=int, help="Number of threads",
                               default=1, )
    config_parser.add_argument("--verbose", '-v', type=bool, help="Verbosity (T/F)",
                               default=True, )
    
    # community subcommand
    community_parser = subparsers.add_parser("community", help="Set up a community")
    
    # file args: not saved in config
    community_parser.add_argument('--config', '-c', required=True, type=str,
                                  help='Path of configuration')
    
    # locations
    community_parser.add_argument('--taxlist', '-t', required=True, type=str,
                                  help='newline-delimited file of taxa in community')
    community_parser.add_argument('--output', '-o', required=True, type=str,
                                  help='Where to write files')
    
    # community params: saved in config
    community_parser.add_argument('--n_comms', '-n', default=5, type=int,
                                  help='Number of communities built with this profile')
    community_parser.add_argument('--n_species', '-s', default=50, type=int,
                                  help='Number of species in each of N communities')
    community_parser.add_argument('--power_a', '-a', default=0.5, type=float,
                                  help='Power distribuion A parameter')
    community_parser.add_argument('--n_strains', '-x', default=2, type=int,
                                  help='Number of same-species strains to include in communities')

    # Simulate subcommand
    simulate_parser = subparsers.add_parser("simulate", help="Run a genome simulation")
    
    
    simulate_parser.add_argument("--config", '-c', type=str, help="Location of configuration",
                               required=True)
    simulate_parser.add_argument('--output', '-o', required=True, type=str,
                                  help='Where to write files')
    simulate_parser.add_argument('--use_sims', '-s', required=False, type=str,
                                 default = '1,2,3,4,5',
                                  help='Which simulations to use (1-indexed)')
    simulate_parser.add_argument("--n_reads", '-n', type=int, help="Number of reads per simulation")
    simulate_parser.add_argument('--alt_dbs', '-a', type=bool, help='whether to simulate alternate genome databases')

    args = parser.parse_args()
    
    if args.command == 'config':

        config = cf.load_configuration(os.path.abspath(args.config), write_new=True)
        
        config.set('database', 'gtdb_loc', os.path.abspath(args.gtdb))
        config.set('locations', 'config', os.path.abspath(args.config))
        config.set('runtime', 'threads', args.threads)
        config.set('parameters', 'verbose', args.verbose)

        config._save_config(config.config_path)

    elif args.command == "community":

        # check tax list exists
        assert os.path.exists(os.path.abspath(args.taxlist)), 'Taxon list not found'
        # check output loc is accessible
        assert os.access(os.path.abspath(args.output), os.W_OK), 'Output loc not found'
        
        # load config and setup logging
        config = cf.load_configuration(os.path.abspath(args.config), write_new=False)
        utils.configure_output(os.path.abspath(args.output))
        logger = utils.setup_logging_for_function('community')
        lvargs = {'logger':logger, 'verbose':config.get('parameters', 'verbose')}
        
        # set params in config
        config.set('parameters', 'n_comms', args.n_comms)
        config.set('parameters', 'n_species', args.n_species)
        config.set('parameters', 'power_a', args.power_a)
        config.set('parameters', 'n_strains', args.n_strains)

        # configure outputs for each module
        config.set('locations', 'outputs', os.path.abspath(args.output))
        community.configure_output(os.path.abspath(args.output))

        # set path to file of taxa -> genome matches
        matches_path = os.path.join(os.path.abspath(args.output), 'matches.tsv.gz')
        config.set('locations', 'matches_path', matches_path)

        # get matches of taxa -> genomes
        matches = community.get_matching_gtdb(config.get('database', 'gtdb_loc'),
                                              os.path.abspath(args.taxlist), 
                                              **lvargs)
        
        # get list of genomes to download
        genbanks = list(matches['top_match_genbank'].dropna().unique()) + list(matches['alt_genbank'].dropna().unique())

        # download genomes
        gate = community.download_genomes(genbanks, **lvargs)
        if gate and os.path.exists(matches_path):
            # rename and add mash distances
            community.rename_files(genbanks, **lvargs)
            matches = community.add_mashdist(matches, **lvargs)

        # save match file
        matches.to_csv(matches_path, sep='\t', compression='gzip')
        
        # get genome lengths and mapping of GTDB accessions to genbanks
        genome_lengths, acc2genbank = community.get_genome_lengths(config.get('database', 'gtdb_loc'))
        # get genome (genbank) to file mapping from output dir
        genome2file = community.get_genome2file()

        # generate file containing actual abundances
        simdata = community.generate_simulations(logger, matches, 
                                       n_sims=args.n_comms, 
                                       n_species=args.n_species, 
                                       power_a=args.power_a, 
                                       n_strains=args.n_strains, 
                                       )
        # save simulation data
        sim_path = os.path.join(os.path.abspath(args.output), 'simulation_data.tsv.gz')
        config.set('locations', 'simulations_path', sim_path)
        simdata.to_csv(sim_path, sep='\t', compression='gzip')

        config._save_config(config.config_path)
    
    elif args.command == "simulate":

        # check output loc is accessible
        assert os.access(os.path.abspath(args.output), os.W_OK), 'Output loc not found'
        
        # load config and setup logging
        config = cf.load_configuration(os.path.abspath(args.config), write_new=False)

        # load simdata
        simdata = pd.read_csv(config.get('locations', 'simulations_path'), sep='\t', index_col=0)

        output_loc = config.get('locations', 'outputs')
        if os.path.exists(os.path.join(output_loc, 'genome_lengths.pkl')):

            genome_lengths = pd.read_pickle(os.path.join(output_loc, 'genome_lengths.pkl'))
            acc2genbank = pd.read_pickle(os.path.join(output_loc, 'acc2genbank.pkl'))
            genome2file = pd.read_pickle(os.path.join(output_loc, 'genome2file.pkl'))

        else:
            # get genome lengths and mapping of GTDB accessions to genbanks
            genome_lengths, acc2genbank = community.get_genome_lengths(config.get('database', 'gtdb_loc'))
            # get genome (genbank) to file mapping from output dir
            genome2file = community.get_genome2file()

            genome_lengths.to_pickle(os.path.join(output_loc, 'genome_lengths.pkl'))
            acc2genbank.to_pickle(os.path.join(output_loc 'acc2genbank.pkl'))
            genome2file.to_pickle(os.path.join(output_loc, 'genome2file.pkl'))

        # construct simulated communites
        simulate.simulate(simdata, 
                 N_READS=args.n_reads, 
                 n_threads=config.get('runtime', 'threads'),
                 acc2genbank=acc2genbank, 
                 genome_lengths=genome_lengths, 
                 genome2file=genome2file,
                 verbose=config.get('parameters', 'verbose'))
        # if args.alt_dbs:
        #     PRINT("RUNNING EVEN THOUGH I SHOULDN'T")
        #     # run sylph with dbs
        #     simulate.run_sylph(simdata, 
        #                        n_threads=cf.config.get('runtime', 'threads'),
        #                        genome2file=genome2file, 
        #                        acc2genbank=acc2genbank,
        #                        alt_dbs=args.alt_dbs,
        #                        verbose=cf.config.get('parameters', 'verbose')
        #                     )

if __name__ == "__main__":
    main()
