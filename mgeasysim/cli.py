import argparse
import sys
import os

# Add the package's root directory to the sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


from mgeasysim.config import config
config._clear()
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
    community_parser.add_argument('--jupyter', '-j', default=False, type=bool,
                                  help='Whether functions are being run in jupyter')

    # Simulate subcommand
    simulate_parser = subparsers.add_parser("simulate", help="Run a genome simulation")
    
    simulate_parser.add_argument("--arg1", type=str, help="First argument for genome simulation")
    simulate_parser.add_argument("--arg2", type=int, help="Second argument for genome simulation")

    args = parser.parse_args()
    
    if args.command == 'config':
        config._clear()
        os.makedirs(os.path.abspath(args.output), exist_ok=True)
        config.set('database', 'gtdb_loc', os.path.abspath(args.gtdb))
        config.set('locations', 'output', os.path.abspath(args.output))
        config.set('runtime', 'threads', args.threads)

    elif args.command == "community":

        matches = matches.get_matching_gtdb(args.taxlist)
        genbanks = list(matches['top_match_genbank'].dropna().unique()) + list(matches['alt_genbank'].dropna().unique())
        download_genomes(genbanks, jupyter=args.jupyter)
        rename_files(genbanks)
        matches = add_mashdist(matches)
        simdata = generate_simulations(matches, 
                                       n_sims=args.n_sims, 
                                       n_species=args.species, 
                                       power_a=args.a, 
                                       n_strains=args.strains
                                       )
    elif args.command == "simulate":
        simulate.run(args.arg1, args.arg2)

if __name__ == "__main__":
    main()
