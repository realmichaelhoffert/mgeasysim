from mgeasysim.config import config

def run(arg1, arg2):
    db_path = config.get("database", "path")
    print(f"Using database at {db_path}")
    # Rest of the code