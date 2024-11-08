import argparse
import sys
import os

# Add the package's root directory to the sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from mgeasysim import community, simulate
from mgeasysim.config import config

def main():
    parser = argparse.ArgumentParser(prog="mgeasysim", description="MGEasySim: Microbial Genome Easy Simulator")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Community subcommand
    config_parser = subparsers.add_parser("config", help="Configure database locations")
    config_parser.add_argument("--gtdb", '-g', type=str, help="Location of GTDB database", 
                               required=True,
                               )
    config_parser.add_argument("--output", '-o', type=str, help="Location of outfile",
                               required=True)
    config_parser.add_argument("--threads", '-t', type=int, help="Number of threads",
                               default=1, )

    # Simulate subcommand
    simulate_parser = subparsers.add_parser("simulate", help="Run a genome simulation")
    simulate_parser.add_argument("--arg1", type=str, help="First argument for genome simulation")
    simulate_parser.add_argument("--arg2", type=int, help="Second argument for genome simulation")

    args = parser.parse_args()
    
    if args.command == 'config':
        
        # config._clear()

        config.set('database', 'gtdb_loc', args.gtdb)
        config.set('locations', 'output', args.output)
        config.set('runtime', 'threads', args.threads)

    elif args.command == "community":
        community.run(args.arg1, args.arg2)
    elif args.command == "simulate":
        simulate.run(args.arg1, args.arg2)

if __name__ == "__main__":
    main()
