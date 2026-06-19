import pandas as pd
import pickle as pkl
import multiprocessing as mp
import os, random, re, warnings, gc
import numpy as np
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
		kmer_seq, _ = kmers(dna_seq, res_methyl_seq, k=k)

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


def _build_region_bounds(cpg_overlaps: pd.DataFrame) -> pd.DataFrame:
	"""
	Derive per-region CpG-index boundaries from the per-CpG cpg_overlaps table.
	Returns one row per dmr_label with region_startCpG and region_endCpG.
	"""
	reset = cpg_overlaps.reset_index()  # brings CpG integer index into 'index' column
	return (
		reset.groupby("dmr_label", as_index=False)
		.agg(
			chr=("chr", "first"),
			region_startCpG=("index", "min"),
			region_endCpG=("index", "max"),
			dmr_ctype=("dmr_ctype", "first"),
			dmr_coordinates=("dmr_coordinates", "first"),
		)
	)


def _assign_atlas_regions_chr(
	reads_chr: pd.DataFrame,
	regions_chr: pd.DataFrame,
) -> pd.DataFrame:
	"""
	Assign each read to the first overlapping atlas region using UXM's any-overlap
	condition in CpG-index space:
		read_startCpG <= region_endCpG  AND  read_endCpG >= region_startCpG

	Also computes cpgs_in_region = size of the overlap in CpG positions (including
	dots), matching how wgbstools counts for its --rlen filter.

	Reads with no matching region are dropped.
	Adds columns: region_startCpG, region_endCpG, dmr_label, dmr_ctype,
	              dmr_coordinates, cpgs_in_region.
	"""
	if reads_chr.empty or regions_chr.empty:
		return reads_chr.iloc[:0]

	reads = reads_chr.copy().sort_values("index").reset_index(drop=True)
	reads["_read_end"] = reads["index"] + reads["methyl"].str.len() - 1

	regions = regions_chr.sort_values("region_startCpG").reset_index(drop=True)
	r_starts = regions["region_startCpG"].to_numpy()
	r_ends   = regions["region_endCpG"].to_numpy()

	matched_idx = np.full(len(reads), -1, dtype=int)
	j = 0
	for i in range(len(reads)):
		rs  = int(reads.at[i, "index"])
		re_ = int(reads.at[i, "_read_end"])
		# advance past regions that end before this read starts
		while j < len(r_starts) and r_ends[j] < rs:
			j += 1
		# first remaining region: overlaps if it starts before read ends
		if j < len(r_starts) and r_starts[j] <= re_:
			matched_idx[i] = j

	mask    = matched_idx >= 0
	reads   = reads[mask].reset_index(drop=True)
	matched = regions.iloc[matched_idx[mask]].reset_index(drop=True)

	reads["region_startCpG"] = matched["region_startCpG"].values
	reads["region_endCpG"]   = matched["region_endCpG"].values
	reads["dmr_label"]       = matched["dmr_label"].values
	reads["dmr_ctype"]       = matched["dmr_ctype"].values
	reads["dmr_coordinates"] = matched["dmr_coordinates"].values
	reads["cpgs_in_region"]  = (
		np.minimum(reads["_read_end"], reads["region_endCpG"]) -
		np.maximum(reads["index"],     reads["region_startCpG"]) + 1
	).astype(int)
	return reads.drop(columns=["_read_end"])


def simulate_reads_chr(df_reads: pd.DataFrame,
					   cpg_overlaps: pd.DataFrame,
					   ref_string: str,
					   chromosome: str,
					   k: int = 3,
					   split_paired_end_reads: bool = True,
					   split_long_reads: bool = True):
	"""Simulate reads in the given chromosome."""

	def single_read_simulation(read, first_idx, last_idx, methyl_seq):
		# Safety: walk boundary indices into cpg_overlaps (for paired-end splits that
		# may land on a dot position not present in the atlas).
		while first_idx not in cpg_overlaps.index:
			first_idx += 1
			methyl_seq = methyl_seq[1:]
			if not methyl_seq:
				return None
		while last_idx not in cpg_overlaps.index:
			last_idx -= 1
			methyl_seq = methyl_seq[:-1]
			if not methyl_seq:
				return None

		start, end = simulate_read_start_end(cpg_overlaps, first_idx, last_idx)
		original_dna_seq = ref_string[start : end + 1]
		kmer_seq, methyl_out = get_processed_sequences(
			start, end,
			methyl_seq=methyl_seq,
			dna_seq=original_dna_seq,
			first_cpg_idx=first_idx,
			last_cpg_idx=last_idx,
			k=k,
			cpg_overlaps=cpg_overlaps,
		)
		return pd.DataFrame({
			"ref_name":        [chromosome],
			"ref_pos":         [start],
			"original_seq":    [original_dna_seq],
			"dna_seq":         [kmer_seq],
			"original_methyl": [read["original_methyl"]],
			"methyl_seq":      ["".join(methyl_out)],
			"dmr_label":       [read["dmr_label"]],
			"dmr_ctype":       [read["dmr_ctype"]],
			"dmr_coordinates": [read["dmr_coordinates"]],
			"cpgs_in_region":  [read["cpgs_in_region"]],
		})

	def paired_end_as_single(read, first_idx, last_idx, methyl_seq, missing_idces):
		"""
		Process a paired-end read as one unit without splitting.
		The unsequenced gap between the two mates is filled with 'N' in original_seq.
		Dots in the methylation pattern are kept and encoded as '2' in methyl_seq.
		"""
		while first_idx not in cpg_overlaps.index:
			first_idx += 1
			methyl_seq = methyl_seq[1:]
			if not methyl_seq:
				return None
		while last_idx not in cpg_overlaps.index:
			last_idx -= 1
			methyl_seq = methyl_seq[:-1]
			if not methyl_seq:
				return None

		start, end = simulate_read_start_end(cpg_overlaps, first_idx, last_idx)
		original_dna_seq = ref_string[start : end + 1]

		# N-fill the unsequenced region between the two mates.
		# The gap in genomic space runs from just after the G of the last mate-1
		# CpG up to just before the C of the first mate-2 CpG.
		mate1_last_cpg_idx  = first_idx + missing_idces[0] - 1
		mate2_first_cpg_idx = first_idx + missing_idces[-1] + 1
		if (mate1_last_cpg_idx  in cpg_overlaps.index and
				mate2_first_cpg_idx in cpg_overlaps.index):
			gap_seq_start = cpg_overlaps.loc[mate1_last_cpg_idx,  "start"] + 2 - start
			gap_seq_end   = cpg_overlaps.loc[mate2_first_cpg_idx, "start"]     - start
			n_filled_seq  = (
				original_dna_seq[:gap_seq_start]
				+ "N" * max(0, gap_seq_end - gap_seq_start)
				+ original_dna_seq[gap_seq_end:]
			)
		else:
			n_filled_seq = original_dna_seq

		# k-mers and methyl encoding use the actual reference (not N-filled) so
		# that CpG detection via re.finditer("CG", ...) remains intact.
		kmer_seq, methyl_out = get_processed_sequences(
			start, end,
			methyl_seq=methyl_seq,
			dna_seq=original_dna_seq,
			first_cpg_idx=first_idx,
			last_cpg_idx=last_idx,
			k=k,
			cpg_overlaps=cpg_overlaps,
		)
		return pd.DataFrame({
			"ref_name":        [chromosome],
			"ref_pos":         [start],
			"original_seq":    [n_filled_seq],
			"dna_seq":         [kmer_seq],
			"original_methyl": [read["original_methyl"]],
			"methyl_seq":      ["".join(methyl_out)],
			"dmr_label":       [read["dmr_label"]],
			"dmr_ctype":       [read["dmr_ctype"]],
			"dmr_coordinates": [read["dmr_coordinates"]],
			"cpgs_in_region":  [read["cpgs_in_region"]],
		})

	df_res = []

	for read_idx in range(df_reads.shape[0]):
		read = df_reads.iloc[read_idx, :].copy()
		read["original_methyl"] = read["methyl"]   # preserve full original pattern

		region_start = int(read["region_startCpG"])
		region_end   = int(read["region_endCpG"])
		read_start   = int(read["index"])
		read_end     = read_start + len(read["methyl"]) - 1

		# ── 1. Clip methyl pattern to atlas region boundaries ────────────────
		clip_left  = max(0, region_start - read_start)
		clip_right = max(0, read_end - region_end)

		methyl_pattern = read["methyl"]
		if clip_left > 0:
			methyl_pattern = methyl_pattern[clip_left:]
			read["index"] = read_start + clip_left
		if clip_right > 0:
			methyl_pattern = methyl_pattern[: len(methyl_pattern) - clip_right]

		if not methyl_pattern:
			continue

		# ── 2. Trim leading dots ─────────────────────────────────────────────
		n_lead = len(methyl_pattern) - len(methyl_pattern.lstrip("."))
		if n_lead >= len(methyl_pattern):
			continue
		if n_lead:
			methyl_pattern = methyl_pattern[n_lead:]
			read["index"]  = int(read["index"]) + n_lead

		# ── 3. Trim trailing dots ────────────────────────────────────────────
		n_trail = len(methyl_pattern) - len(methyl_pattern.rstrip("."))
		if n_trail >= len(methyl_pattern):
			continue
		if n_trail:
			methyl_pattern = methyl_pattern[:-n_trail]

		if not methyl_pattern:
			continue

		read["methyl"] = methyl_pattern
		read["nCG"]    = len(methyl_pattern)
		first_idx = int(read["index"])
		last_idx  = first_idx + len(methyl_pattern) - 1

		# Safety: both boundary CpGs must be in cpg_overlaps for sequence lookup
		if first_idx not in cpg_overlaps.index or last_idx not in cpg_overlaps.index:
			continue

		min_seq_len = (
			cpg_overlaps.loc[last_idx, "start"] - cpg_overlaps.loc[first_idx, "start"] + 1
		)

		for _ in range(read["n_reads"]):
			if min_seq_len <= 150 and "." not in methyl_pattern:
				# Short single-end read — no splitting needed regardless of flags
				res = single_read_simulation(read, first_idx, last_idx, methyl_pattern)
				if res is not None:
					df_res.append(res)
			else:
				missing_idces = [m.start() for m in re.finditer(r"\.", methyl_pattern)]
				# Paired-end: all dots form one consecutive block
				is_paired = bool(missing_idces) and all(
					b - a == 1
					for a, b in zip(missing_idces, missing_idces[1:])
				)

				if is_paired:
					if split_paired_end_reads:
						gap_start = missing_idces[0]
						gap_end   = missing_idces[-1]
						res = single_read_simulation(
							read, first_idx, first_idx + gap_start - 1,
							methyl_pattern[:gap_start],
						)
						if res is not None:
							df_res.append(res)
						res = single_read_simulation(
							read, first_idx + gap_end + 1, last_idx,
							methyl_pattern[gap_end + 1:],
						)
						if res is not None:
							df_res.append(res)
					else:
						res = paired_end_as_single(
							read, first_idx, last_idx, methyl_pattern, missing_idces
						)
						if res is not None:
							df_res.append(res)
				else:
					# Long read (no dot gap or non-consecutive dots)
					if split_long_reads:
						half = len(methyl_pattern) // 2
						res = single_read_simulation(
							read, first_idx, first_idx + half - 1, methyl_pattern[:half]
						)
						if res is not None:
							df_res.append(res)
						res = single_read_simulation(
							read, first_idx + half, last_idx, methyl_pattern[half:]
						)
						if res is not None:
							df_res.append(res)
					else:
						res = single_read_simulation(
							read, first_idx, last_idx, methyl_pattern
						)
						if res is not None:
							df_res.append(res)

	return pd.concat(df_res) if df_res else None


def simulate_reads(df_reads: pd.DataFrame,
				   cpg_overlaps: pd.DataFrame,
				   wgbstools_ref_dir: str,
				   genome: str = "hg19",
				   n_cores: int = 10,
				   split_paired_end_reads: bool = True,
				   split_long_reads: bool = True) -> pd.DataFrame:

	df_reads = df_reads.copy()
	df_reads["nCG"] = df_reads["methyl"].str.len()

	# Build one row per atlas region with its CpG-index boundaries
	region_bounds = _build_region_bounds(cpg_overlaps)

	# Load reference genome
	f_genome = get_reference_file(
		genome=genome, file_type="genome", wgbstools_ref_dir=wgbstools_ref_dir
	)
	with open(f_genome, "rb") as fp:
		dict_ref = pkl.load(fp)

	# Per-chromosome: assign atlas regions via any-overlap, then queue for simulation.
	# chrX and chrY are included to match UXM's chromosome coverage.
	list_args = []
	chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
	for chrom in chroms:
		sub_reads = df_reads[df_reads["chr"] == chrom]
		sub_cpgs  = cpg_overlaps[cpg_overlaps["chr"] == chrom]
		sub_regs  = region_bounds[region_bounds["chr"] == chrom]

		if sub_reads.empty or sub_cpgs.empty:
			print(chrom, sub_reads.shape, sub_cpgs.shape)
			continue

		# Overlap join: assigns region info and cpgs_in_region to each read
		sub_reads = _assign_atlas_regions_chr(sub_reads, sub_regs)
		if sub_reads.empty:
			continue

		if chrom not in dict_ref:
			print(f"Reference genome missing for {chrom}, skipping")
			continue

		list_args.append((sub_reads, sub_cpgs, dict_ref[chrom], chrom, 3,
						  split_paired_end_reads, split_long_reads))

	del dict_ref
	gc.collect()

	with mp.Pool(n_cores) as pool:
		processed_reads = pool.starmap(simulate_reads_chr, list_args)

	processed_reads = [p for p in processed_reads if p is not None]
	if not processed_reads:
		return pd.DataFrame()
	return pd.concat(processed_reads).reset_index(drop=True)
