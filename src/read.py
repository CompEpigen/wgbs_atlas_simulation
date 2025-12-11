import pandas as pd
import pickle as pkl
import multiprocessing as mp
import os, random, re, warnings, gc
from database import get_reference_file

def kmers(seq: str, methyl_seq: list, k=3):
	converted_seq = list()
	for seq_idx in range(len(seq)-k+1):
		token = seq[seq_idx:seq_idx+k]
		converted_seq.append(token)

	n_rm_methyl = int((k-1)/2)
	methyl_seq = methyl_seq[n_rm_methyl:-n_rm_methyl]

	if len(methyl_seq) != len(converted_seq):
		raise ValueError(f"Length of methylation and k-mer DNA seqs are different (%{len(methyl_seq)} vs %{len(converted_seq)})")

	return " ".join(converted_seq), methyl_seq

def get_processed_sequences(start: int, end: int, methyl_seq: str, dna_seq: str, first_cpg_idx: int, last_cpg_idx: int, k: int, cpg_overlaps):
	
	# methylation pattern conversion match in the data 
	methyl_patterns={"T":"0","C":"1", ".":"2"}

	cpg_idx = [t.start() for t in re.finditer("CG",dna_seq)]

	if len(cpg_idx) != len(methyl_seq):
		print(start, end, first_cpg_idx, last_cpg_idx)
		raise ValueError(f"Number of CpGs doesn't match: {len(cpg_idx)} in {dna_seq}\nGiven methyl pattern {methyl_seq}")

	res_methyl_seq = ["2" for i in range(len(dna_seq))]
	for c, m in zip(cpg_idx, methyl_seq):
		res_methyl_seq[c] = methyl_patterns[m]

	if len(res_methyl_seq) != len(dna_seq):
		raise ValueError(f"methyl len : {len(res_methyl_seq)} / DNA len: {len(dna_seq)}")

	# K-mer sequences
	if k>0:
		kmer_seq, _ = kmers(dna_seq, res_methyl_seq, k=k) # Fix: do not overwrite methyl_seq with k_mers variant so that all models that use original sequence can learn correctly. For MethylBERT it means that additional postprocessing must be done. 

	return kmer_seq, res_methyl_seq

def simulate_read_start_end(cpg_overlaps: pd.DataFrame, 
							first_idx: int, 
							last_idx: int, 
							seq_len: int=150):
	"""
		Find the start and end position of simulated read 

		cpg_overlaps: pd.DataFrame
			DataFrame of CpGs overlapping with the given region
		first_idx, last_idx: int, int
			Index of the first and the last index 
	"""
	start, end = cpg_overlaps.loc[first_idx, "start"]-1, cpg_overlaps.loc[last_idx, "start"] + 1 # -1 and + 1 to avoid methyl patterns at the beginning and the end of the read

	additional_bps = seq_len - (end - start + 1)
	allowed_bps_start = start - cpg_overlaps.loc[first_idx, "prev_cpg"] - 1
	allowed_bps_end = cpg_overlaps.loc[last_idx, "next_cpg"] - end

	if additional_bps > 10: # if the read is too short
		additional_start = random.randint(1, additional_bps)
		if additional_start > allowed_bps_start:
			additional_start = allowed_bps_start

		additional_end = additional_bps - additional_start
		if additional_end > allowed_bps_end:
			additional_end = allowed_bps_end

		start -= additional_start
		end += additional_end

	return start, end



def simulate_reads_chr(df_reads: pd.DataFrame, 
					   cpg_overlaps: pd.DataFrame, 
					   ref_string: str, 
					   chromosome: str, 
					   k: int=3):
	'''
	Simulate reads in the given chromosome 
	'''

	def single_read_simulation(read: pd.DataFrame, first_idx: int, last_idx: int, methyl_seq: str):

		# The second read is not fully overlapping with the region
		n_count = 0
		while first_idx not in cpg_overlaps.index:
			#print(methyl_seq, first_idx, cpg_overlaps.index)
			first_idx += 1
			n_count+=1
			methyl_seq = methyl_seq[1:]
			if len(methyl_seq) == 0:
				print(read, first_idx, last_idx, methyl_seq)
				raise ValueError(f"first_idx disappeared: {first_idx}, {last_idx}, {n_count}")
				return None

		n_count = 0
		while last_idx not in cpg_overlaps.index:
			#print(methyl_seq, first_idx, cpg_overlaps.index)
			last_idx -= 1
			methyl_seq = methyl_seq[:-1]
			n_count+=1
			
			if len(methyl_seq) == 0:
				print(read, first_idx, last_idx, methyl_seq)
				raise ValueError(f"last_idx disappeared: {first_idx}, {last_idx}, {n_count}")
				return None

		start,end =  simulate_read_start_end(cpg_overlaps, first_idx, last_idx)
		original_dna_seq = ref_string[start:end+1]
		kmer_seq, methyl_seq = get_processed_sequences(start, end, 
													   methyl_seq=methyl_seq,
													   dna_seq=original_dna_seq,
													   first_cpg_idx=first_idx,
													   last_cpg_idx=last_idx,
													   k=k, cpg_overlaps=cpg_overlaps)
		
		return pd.DataFrame({"ref_name": [chromosome],
								"ref_pos": [start],
								"original_seq" : [original_dna_seq],
								"dna_seq": [kmer_seq],
								"original_methyl": [read["methyl"]],
								"methyl_seq": ["".join(methyl_seq)],
								"dmr_label": [cpg_overlaps.loc[read["index"], "dmr_label"]],
								"dmr_ctype": [cpg_overlaps.loc[read["index"], "dmr_ctype"]],
								"dmr_coordinates": [cpg_overlaps.loc[read["index"], "dmr_coordinates"]]})

	df_res = list()
		
	for read_idx in range(df_reads.shape[0]):
		read = df_reads.iloc[read_idx,:]

		# Remove missing patterns at the beginning and the end
		read["original_pattern"] = read["methyl"]
		methyl_pattern = read["methyl"]
		n_count = 0
		for meth_idx, i in enumerate(methyl_pattern):
			if (i != ".") and (read["index"] + meth_idx in cpg_overlaps.index):
				break
			else:
				n_count += 1

		if n_count >= len(methyl_pattern):
			continue
		elif n_count > 0:
			methyl_pattern = methyl_pattern[n_count:]
			read["index"] += n_count
			read["nCG"] -= n_count

		n_count = 0
		for meth_idx, i in enumerate(methyl_pattern[::-1]):
			if (i != ".") and \
			   (read["index"] + len(methyl_pattern) - (meth_idx+1) in cpg_overlaps.index):
				break
			else:
				methyl_pattern = methyl_pattern[:-1]
				read["nCG"] -= 1

		if n_count >= len(methyl_pattern):
			continue
		elif n_count > 0:
			methyl_pattern = methyl_pattern[:-n_count]
			read["nCG"] -= n_count

		read["methyl"] = methyl_pattern
		last_idx = read["index"] + read["nCG"] -1

		if (last_idx not in cpg_overlaps.index):
			continue # the read is not fully overlapping 

		min_seq_len = cpg_overlaps.loc[last_idx, "start"] - cpg_overlaps.loc[read["index"], "start"] + 1
		#min_pos, max_pos = cpg_overlaps.loc[read["index"], "prev_cpg"]+1, cpg_overlaps.loc[last_idx, "next_cpg"]
		
		for r in range(read["n_reads"]):
			# Read length is <=150bp - single-read
			if ( min_seq_len <= 150 ) and \
			   ("." not in read["methyl"]) and \
			   (len(read["methyl"]) > 0):

				res = single_read_simulation(read=read, 
											 first_idx=read["index"], last_idx=last_idx, 
											 methyl_seq=read["methyl"])
				df_res.append(res)
				
			else:
				# Read len > 150 bp - divide the read into two reads

				missing_cpg_idces = [i.start() for i in re.finditer("\.", read["methyl"])]
				is_distant_read_pairs = (len(missing_cpg_idces) != 0)
				#if (missing_cpg_idces[0] == 0) or (missing_cpg_idces[-1] == len(read["methyl"]))

				# When the missing CpGs are not consecutive 
				# Assuming that paired-end reads have unprofiled CpGs beween two reads
				for a, b in zip(missing_cpg_idces[1:], missing_cpg_idces[:-1]):
					if a-b != 1:
						is_distant_read_pairs=False

				if is_distant_read_pairs: # Paried-end reads
					# First read
					first_idx = read["index"]
					last_idx = read["index"] + missing_cpg_idces[0] -1 # change the last idx of the read to the last idx of the last CpG in the first read

					res = single_read_simulation(read=read, 
												 first_idx=first_idx, last_idx=last_idx, 
												 methyl_seq=read["methyl"][:missing_cpg_idces[0]])
					if res is not None: df_res.append(res)

					# Second read
					first_idx = read["index"] + missing_cpg_idces[-1] + 1
					last_idx = read["index"] + read["nCG"] -1

					res = single_read_simulation(read=read, 
												 first_idx=first_idx, last_idx=last_idx, 
												 methyl_seq=read["methyl"][missing_cpg_idces[-1]+1:])
					if res is not None: df_res.append(res)

				else: # There is no missing part. The fragment is divded into two parts at the middle
					half_idx = int(len(read["methyl"])/2)
					# Divide the fragment into two reads
					# First read
					last_idx = read["index"] + half_idx -1
					#if last_idx in cpg_overlaps.index:

					res = single_read_simulation(read=read, 
											 first_idx=read["index"], last_idx=last_idx, 
											 methyl_seq=read["methyl"][:half_idx])
					if res is not None: df_res.append(res)

					# Second read
					first_idx = read["index"] + half_idx
					last_idx = read["index"] + read["nCG"] -1
					#if first_idx in cpg_overlaps.index:
					res = single_read_simulation(read=read, 
												 first_idx=first_idx, last_idx=last_idx, 
												 methyl_seq=read["methyl"][half_idx:])
					if res is not None: df_res.append(res)
				
	return pd.concat(df_res) if len(df_res) > 0 else None
	


def simulate_reads(df_reads: pd.DataFrame, 
				   cpg_overlaps: pd.DataFrame,
				   wgbstools_ref_dir: str, 
				   genome: str = "hg19", 
				   n_cores: int = 10):
	# Filter reads not fully overlapping with the regions
	df_reads = df_reads[df_reads["index"].isin(cpg_overlaps.index)]
	df_reads.loc[:, "nCG"] = df_reads['methyl'].apply(lambda x: len(x))
	DATA_DIR = os.path.join(os.getcwd(), "data") 

	# Read reference genome saved in a dictionary 
	f_genome = get_reference_file(genome=genome, file_type='genome', wgbstools_ref_dir=wgbstools_ref_dir)
	with open(f_genome, "rb") as fp:
			dict_ref = pkl.load(fp)

	# Arguments for multiprocessing - divide reads into each chromosome
	list_args = list()
	for chr_idx in range(1, 23):
		sub_reads = df_reads[df_reads["chr"]=="chr%d"%chr_idx]
		sub_cpgs = cpg_overlaps[cpg_overlaps["chr"]=="chr%d"%chr_idx]
		if (sub_reads.shape[0] == 0) or (len(sub_cpgs)==0):
			print(chr_idx, sub_reads.shape, sub_cpgs.shape)
			continue
		else:
			list_args.append((sub_reads, 
							  sub_cpgs, 
							  dict_ref["chr%d"%chr_idx], 
							  "chr%d"%chr_idx, 
							  3)) # order according to the "simulate_reads_chr" function
		del sub_reads, sub_cpgs 	

	# Multiprocessing - simulate reads in each chromosome	
	with mp.Pool(n_cores) as pool:
		processed_reads = pool.starmap(simulate_reads_chr, list_args)

	# Merge simulated reads 
	processed_reads = [p for p in processed_reads if p is not None]
	del dict_ref
	gc.collect()

	processed_reads = pd.concat(processed_reads).reset_index(drop=True)
	
	return processed_reads
