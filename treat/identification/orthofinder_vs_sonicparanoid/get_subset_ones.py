import os
import shutil
import concurrent.futures
from pathlib import Path

"""Given orphan orthogroups from both tools, identify subset orthogroups"""

of_dir = Path("of_ogs")
sp_dir = Path("sp_ogs")
subset_dir = Path("subset_ones")
subset_dir.mkdir(exist_ok=True)

def is_subset(of_file, sp_file):
    with open(of_file, 'r') as f:
        of_lines = set(f.read().splitlines())
    
    with open(sp_file, 'r') as f:
        sp_lines = set(f.read().splitlines())
    
    if of_lines == sp_lines:
        return None  
    
    if of_lines.issubset(sp_lines):
        return sp_file  
    elif sp_lines.issubset(of_lines):
        return of_file  
    else:
        return None

def copy_superset_file(pair):
    of_file, sp_file = pair
    superset_file = is_subset(of_file, sp_file)
    
    if superset_file:
        of_filename = of_file.stem
        sp_filename = sp_file.stem
        dest_file = subset_dir / f"{of_filename}_{sp_filename}_sub.txt"
        
        shutil.copy(superset_file, dest_file)
        print(f"Copied {superset_file.name} to {dest_file.name}")

def main():
    of_files = list(of_dir.glob("*.txt"))
    sp_files = list(sp_dir.glob("*.txt"))
    
    pairs = [(of_file, sp_file) for of_file in of_files for sp_file in sp_files]
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=32) as executor:
        executor.map(copy_superset_file, pairs)

if __name__ == "__main__":
    main()

