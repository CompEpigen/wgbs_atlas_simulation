import os
import pickle
import pandas as pd
import gzip
from pathlib import Path
import subprocess

class ReferenceAdapter:
    """Adapter to convert wgbstools references to simulation-ready formats"""
    
    def __init__(self, genome_name='hg19', wgbstools_ref_dir=None):
        self.genome = genome_name
        if wgbstools_ref_dir is not None:
             self.wgbs_dir = Path(wgbstools_ref_dir) / "references" / genome_name
        else:
             self.wgbs_dir = None
        self.cache_dir = Path('data') / genome_name
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        
    def _validate_wgbs_setup(self):
        """Validate that wgbstools reference exists before generation"""
        if self.wgbs_dir is None:
            raise ValueError(
				f"wgbstools_ref_dir must be provided for first-time setup of {self.genome}. "
				f"Run: wgbstools init_genome {self.genome}"
			)
        if not self.wgbs_dir.exists():
            raise ValueError(
				f"Reference {self.genome} not found at {self.wgbs_dir}. "
				f"Run: wgbstools init_genome {self.genome}"
			)
             
    
    def get_cpg_locations(self):
        """Convert CpG.bed.gz to the CSV format expected by simulation"""
        cache_file = self.cache_dir / f"{self.genome}_cpgs.csv"
        if cache_file.exists():
            return str(cache_file)
        
        self._validate_wgbs_setup()
        print(f"Converting CpG locations for {self.genome}...")
        cpg_bed = self.wgbs_dir / "CpG.bed.gz"
        
        # Read bed.gz file (columns: chr, loc, site)
        df = pd.read_csv(cpg_bed, sep='\t', header=None, 
                        names=['chr', 'loc', 'site'])
        
        # Convert to expected format
        df_out = pd.DataFrame({
            'chr': df['chr'],
            'start': df['loc'],  # C position (1-based)
            'end': df['loc'] + 1,  # G position
            'index': df['site']
        })
        
        df_out.to_csv(cache_file, index=False)
        return str(cache_file)
    
    def get_genome_dict(self):
        """Convert FASTA to pickled dictionary"""
        cache_file = self.cache_dir / f"{self.genome}_genome.pk"
        
        if cache_file.exists():
            return str(cache_file)
        
        self._validate_wgbs_setup()
        
        print(f"Creating genome dictionary for {self.genome}...")
        fasta = self.wgbs_dir / "genome.fa.gz"
        if not fasta.exists():
            fasta = self.wgbs_dir / "genome.fa"
        
        # Read FASTA into dictionary
        genome_dict = {}
        current_chr = None
        current_seq = []
        
        opener = gzip.open if str(fasta).endswith('.gz') else open
        with opener(fasta, 'rt') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # Save previous chromosome
                    if current_chr is not None:
                        genome_dict[current_chr] = ''.join(current_seq)
                    # Start new chromosome
                    current_chr = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line.upper())
            
            # Save last chromosome
            if current_chr is not None:
                genome_dict[current_chr] = ''.join(current_seq)
        
        # Pickle it
        with open(cache_file, 'wb') as f:
            pickle.dump(genome_dict, f)
        
        return str(cache_file)
    
    def initialize_genome(self, force=False):
        """Initialize genome if not already done"""
        if not self.wgbs_dir.exists() or force:
            print(f"Initializing {self.genome} with wgbstools...")
            cmd = f"wgbstools init_genome {self.genome}"
            if force:
                cmd += " -f"
            subprocess.check_call(cmd, shell=True)
        
        return self

def get_reference_file(genome='hg19', file_type='cpgs', wgbstools_ref_dir=None):
    """
    Get reference file path, initializing if needed.
    
    Args:
        genome: Genome name (hg19, hg38, mm10, etc.)
        file_type: 'cpgs' or 'genome'
    """
    adapter = ReferenceAdapter(genome,wgbstools_ref_dir)
    
    if file_type == 'cpgs':
        return adapter.get_cpg_locations()
    elif file_type == 'genome':
        return adapter.get_genome_dict()
    else:
        raise ValueError(f"Unknown file_type: {file_type}")
