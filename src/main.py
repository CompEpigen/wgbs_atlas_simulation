import argparse, re, os, random, gc, pathlib, glob,json
import pandas as pd
import numpy as np

from read import simulate_reads
from region import merge_close_regions, find_cpg_overlaps
from database import get_reference_file
from pathlib import Path


pd.options.mode.chained_assignment = None 

def arg_parser():
	# Parse input arguments
	parser = argparse.ArgumentParser()

	parser.add_argument("-o", "--output_dir", type=str, default="./", help="Directory to save the results (default: ./)")
	parser.add_argument("-c", "--cores", type=int, default=1, help="Number of cores for multiprocessing. A larger number increases the computation speed. (default: 1)")
	parser.add_argument("-r", "--f_region", type=str, default=None, help="Selected regions for training MethylBERT.")
	parser.add_argument("-g", "--genome", type=str, default="hg19", help="Reference genome (either hg19 or hg38). Currently only hg19 is available. (default: hg19)")
	parser.add_argument("-f", "--f_input", required=True, help="Text file containing a list of .pat files OR path to a .pat file")
	parser.add_argument("-wd", "--wgbstools_dir", help="The directory where wgbstools is installed. Only needed to be provided once per genome.")
	parser.add_argument("-ov", "--overwrite_processed_files", action=argparse.BooleanOptionalAction, help="Force overwrite of already processed files.")
	

	return parser.parse_args()

def get_input_files(f_input: str) -> pd.DataFrame:
    """
    Parse input argument to get list of .pat files
    
    Args:
        f_input: Can be:
            - A single .pat or .pat.gz file
            - A text file containing list of .pat/.pat.gz files
            - A directory containing .pat/.pat.gz files
    
    Returns:
        DataFrame with columns: files, cell_type
    """
    input_path = Path(f_input)
    
    # Case 1: Single .pat or .pat.gz file
    if input_path.is_file() and (input_path.suffix == '.pat' or input_path.name.endswith('.pat.gz')):
        ctype = filename2ctype(str(input_path))
        return pd.DataFrame({
            "files": [str(input_path)],
            "cell_type": [ctype]
        })
    
    # Case 2: Directory containing .pat files
    elif input_path.is_dir():
        pat_files = list(input_path.glob('*.pat')) + list(input_path.glob('*.pat.gz'))
        if len(pat_files) == 0:
            raise ValueError(f"No .pat or .pat.gz files found in directory: {f_input}")
        
        df_files = pd.DataFrame({
            "files": [str(f) for f in pat_files]
        })
        df_files["cell_type"] = df_files["files"].apply(lambda x: filename2ctype(x))
        print(f"Found {len(pat_files)} .pat/.pat.gz files in directory")
        return df_files
    
    # Case 3: Text file with list of .pat files
    elif input_path.is_file():
        df_files = pd.read_csv(f_input, header=None)
        if df_files.shape[1] != 1:
            raise ValueError("The input file must have only one column containing .pat file names.")
        df_files.columns = ["files"]
        df_files["cell_type"] = df_files["files"].apply(lambda x: filename2ctype(x))
        return df_files
    
    else:
        raise ValueError(f"Invalid input: {f_input}. Must be a .pat file, a directory, or a text file with .pat file paths.")


def filename2ctype(f_name: str) -> str:
	
	'''
	
	Extract cell types from the file name
	It returns the cell type matching to the cell types written in the region information 
	e.g.,) Blood-Monocytes -> Blood-Mono+Macro
	
	'''

	f_name = os.path.basename(f_name)
	file_ctype = "_".join(f_name.split("_")[1:]).split("-Z")[0].strip()
	try:
		return df_ctype_match[file_ctype]
	except KeyError:
		print(f"{f_name} cannot be matched to the 39 cell types. NA is assigned for the cell type")
		# print(f_name, "---->", file_ctype)
		return "NA"

if __name__=="__main__":
	
	args = arg_parser()

	### Set-ups
		# output directory
	if not os.path.exists(args.output_dir):
		os.mkdir(args.output_dir)

	# load cell tyep match dictionary 
	global df_ctype_match
	f_cell_type_match = "data/cell_type_match.json"
	with open(f_cell_type_match, "r") as fp:
		df_ctype_match = json.load(fp)

	# Get input files (handles .pat, .pat.gz, directory, or list file)
	df_files = get_input_files(args.f_input)
	print(f"Processing {len(df_files)} file(s)")

	# Region selection
	if args.f_region is not None:
		if args.f_region.split(".")[-1] == "csv":
			sep=","
		elif args.f_region.split(".")[-1] == "tsv":
			sep="\t"
		else: 
			raise ValueError("The target regions must be either in .csv or .tsv format")
		df_region = pd.read_csv(args.f_region, sep=sep)
	else:
		raise ValueError("The target regions must be provided")

	# Renaming columns
	if ("target" in df_region.columns) and ("Type" not in  df_region.columns):
		df_region.rename(columns={"target":"ctype"}, inplace=True)
	elif "Type" in df_region.columns:
		df_region.rename(columns={"Type":"ctype"}, inplace=True)
	else:
		ValueError("Either 'target' or 'Type' column must be provided in the reference regions")

	unique_ctypes = df_files["cell_type"].unique()
	df_region = df_region.loc[df_region["ctype"].apply(lambda x: x in unique_ctypes), :]
	df_region["dmr_id"] = list(range(df_region.shape[0]))
	print(df_region)
	unique_ctypes=df_region["ctype"].unique()
	print(f"Regions for {unique_ctypes} are selected")

	# Read cpg file 

	print(f"Read {args.genome} cpg files")
	f_cpgs = get_reference_file(genome=args.genome, file_type='cpgs', wgbstools_ref_dir=args.wgbstools_dir)
	df_cpg = pd.read_csv(f_cpgs, sep=",")
	df_cpg = df_cpg.set_index("index")

	# To match 0-based and 1-based
	df_cpg["start"] -= 1
	df_cpg["end"] -= 1

	df_cpg["prev_cpg"] = [-1] + list(df_cpg["start"][:-1]) # cytosine pos of previous CpG
	df_cpg["next_cpg"] = list(df_cpg["start"][1:]) + [-1] # cytosine pos of next CpG

	# Reset columns for finding overlaps
	df_subject = df_cpg.loc[:, ["seqnames", "start", "end", "prev_cpg", "next_cpg"]]
	df_subject.columns = ["chr", "start", "end", "prev_cpg", "next_cpg"]
	df_region.index = df_region["dmr_id"]

	# Find CpGs overlapping with given regions
	print("Find CpGs in the selected regions")
	cpg_overlaps = find_cpg_overlaps(df_region.loc[:, ["chr", "start", "end", "ctype"]], 
									 df_subject, 
									 n_cores = args.cores)
	f_out = os.path.join(args.output_dir, "selected_cpgs.csv")
	cpg_overlaps.to_csv(f_out, header= True, sep="\t", index=True)
	print(cpg_overlaps)
	del df_cpg, df_region, df_subject
	gc.collect()

	for idx in df_files.index:

		# Read input .pat file
		input_file = df_files.loc[idx, "files"]
		try:
			print(f"{input_file} is being processed...")
			# Generate output filename (handle both .pat and .pat.gz)
			base_name = os.path.basename(input_file)
			if base_name.endswith('.pat.gz'):
				out_name = base_name.replace(".pat.gz", "_reads.csv")
			else:
				out_name = base_name.replace(".pat", "_reads.csv")
			
			f_out = os.path.join(args.output_dir, out_name)
			if not args.overwrite_processed_files:
				if(os.path.exists(f_out)):
					continue
			
			# Read .pat file (pandas automatically handles .gz compression)
			df_reads = pd.read_csv(input_file, sep="\t", header=None, compression='infer')
			df_reads.columns = ["chr", "index", "methyl", "n_reads"] 
			# Read simulation 
			res = simulate_reads(df_reads, cpg_overlaps, wgbstools_ref_dir=args.wgbstools_dir, genome=args.genome, n_cores=args.cores)
			res["ctype"] = df_files.loc[idx, "cell_type"]
			res.to_csv(f_out, header= True, sep="\t", index=False)
		except:
			print(f"Error procesing {input_file}. Moving to the next one")