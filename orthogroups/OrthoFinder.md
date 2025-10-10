## OrthoFinder orthology assignment

### Retrieve the longest isoform per gene
Using [AGAT](https://agat.readthedocs.io/en/latest/) v0.8.0.

Running first on sloth bear:
```bash
mkdir long_iso

# Run the following in longIso_sloth.sh:
#!/bin/bash
#SBATCH --partition=kamiak
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --job-name=agat_sloth          ### Job Name
#SBATCH --time=0-06:00:00       ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1               ### Node count required for the job
#SBATCH --out=agat_sloth.out

module load anaconda3
source activate bwp_agat_env

# Filter GFF to retain longest isoform per gene
agat_sp_keep_longest_isoform.pl -gff sloth_out/complete.genomic.gff  -o long_iso/sloth.egapx.longestIsoform.gff

# Extract protein sequences for longest isoforms
agat_sp_extract_sequences.pl -g long_iso/sloth.egapx.longestIsoform.gff -f sloth_out/complete.genomic.fna -p -o long_iso/sloth.egapx.longestIso.prots.faa

# Extract CDS sequences for longest isoforms
agat_sp_extract_sequences.pl -g long_iso/sloth.egapx.longestIsoform.gff -f sloth_out/complete.genomic.fna -t cds -o long_iso/sloth.egapx.longestIso.cds.fna

```

The script above was then modified and run for the remaining seven EGAPx annotations. 

### Run OrthoFinder
We ran OrthoFinder using the eight EGAPx assemblies as well as the existing brown bear annotation from NCBI (UrsArc2.0). That will give us NCBI brown bear gene IDs for each orthogroup and allow us to focus on orthogroups that are supported by the inclusion of a gene from the high quality NCBI annotation. After OrthoFinder, the existing brown bear annotation will be ignored and analyses will focus only on the EGAPx annotations. 

```bash
mkdir orthofinder_newApr21
cd orthofinder_newApr21
mkdir input

cp ../ncbi_annot_pipeline/egapx/long_iso/*faa ./input/

# Run the following in orthofinder_oldBrownAllNew_apr22.sh
#!/bin/bash
#SBATCH --partition=kamiak
#SBATCH --cpus-per-task=38
#SBATCH --job-name=orthofinder        ### Job Name
#SBATCH --time=0-14:00:00       ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1               ### Node count required for the job
#SBATCH --out=orthofinder_out
#SBATCH --mem=64G

module load anaconda3
source activate bwp_orthofinder_env

orthofinder -t 38 -a 8 -f ./input -o ./output
```

OrthoFinder output was parsed using Rscript  `parseOrthogroups_04.07.25.R` and custom python script `combine_fastas_8spec.py` to produce unaligned fastas containing the longest isoform CDS sequence for all single-copy orthologs.

