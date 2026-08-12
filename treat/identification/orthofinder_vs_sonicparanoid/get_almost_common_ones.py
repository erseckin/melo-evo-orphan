import os
import concurrent.futures
from pathlib import Path

"""Given orphan orthogroups from both tools, identify almost common orthogroups given a threshold"""

of_dir = Path("of_ogs")
sp_dir = Path("sp_ogs")
output_dir = Path("almost_common")
output_dir.mkdir(exist_ok=True)

lower_threshold = 75

def calculate_similarity(of_file, sp_file):
    with open(of_file, 'r') as f:
        of_lines = set(f.read().splitlines())
    
    with open(sp_file, 'r') as f:
        sp_lines = set(f.read().splitlines())
    
    intersection = of_lines.intersection(sp_lines)
    union = of_lines.union(sp_lines)
    
    similarity = (len(intersection) / len(union)) * 100
    return similarity, of_lines, sp_lines

def merge_files(of_file, sp_file):
    similarity, of_lines, sp_lines = calculate_similarity(of_file, sp_file)
    
    if lower_threshold <= similarity < 100:
        merged_lines = sorted(of_lines.union(sp_lines))
        output_file = output_dir / f"{of_file.stem}_{sp_file.stem}_almost.txt"
        
        with open(output_file, 'w') as f:
            f.write('\n'.join(merged_lines))
        
        print(f"Merged {of_file.name} and {sp_file.name} into {output_file.name}")

def process_file_pair(pair):
    of_file, sp_file = pair
    merge_files(of_file, sp_file)

def main():
    of_files = list(of_dir.glob("*.txt"))
    sp_files = list(sp_dir.glob("*.txt"))
    
    pairs = [(of_file, sp_file) for of_file in of_files for sp_file in sp_files]
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=32) as executor:
        executor.map(process_file_pair, pairs)

if __name__ == "__main__":
    main()

