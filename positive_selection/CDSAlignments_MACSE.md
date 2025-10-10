## Run MACSE alignment
The following code uses [MACSE v2.07](https://www.agap-ge2pop.org/macse/) to generate CDS alignments for use in HyPhy analyses. 

### Run initial alignment
Running this on longest isoform CDS sequences for all single-copy orthologs detected with OrthoFinder.
```bash
mkdir cds_aligned_macse
mkdir prot_aligned_macse

# Run the following in 1_macseAlign.sh:

#!/bin/bash
#SBATCH --job-name=macse
#SBATCH --time=4-00:00:00
#SBATCH --partition=kamiak
#SBATCH --cpus-per-task=1
#SBATCH --mem=12GB
#SBATCH --mail-type=ALL
#SBATCH --array=1-14656:1 
#SBATCH --output=./logs/JobArray_%A_%a.out
#SBATCH --error=./logs/JobArray_%A_%a.err

file=$(sed -n ''$SLURM_ARRAY_TASK_ID'p' cdsfiles.txt | awk '{print$1}')

module load anaconda3
source activate bwp_macse_env

macse -prog alignSequences -seq $file -out_NT ./cds_aligned_macse/$(basename ${file} .unaligned.8spec.fa).macse.aln.fna -out_AA ./prot_aligned_macse/$(basename ${file} .unaligned.8spec.fa).macse.aln.faa
```

### Remove frameshifts and terminal stop codons
```bash

ls cds_aligned_macse/* > macsealigns.txt
wc -l macsealigns.txt # 14666

# Run the following in 2_macseClean.sh:

#!/bin/bash
#SBATCH --job-name=macse
#SBATCH --time=1-00:00:00
#SBATCH --partition=kamiak
#SBATCH --cpus-per-task=1
#SBATCH --mem=12GB
#SBATCH --mail-type=ALL
#SBATCH --output=./logs/JobArray_%A_%a.out
#SBATCH --error=./logs/JobArray_%A_%a.err

file=$(sed -n ''$SLURM_ARRAY_TASK_ID'p' macsealigns.txt | awk '{print$1}')

module load anaconda3
source activate bwp_macse_env
macse -prog exportAlignment -align $file -codonForInternalStop NNN -codonForFinalStop --- -codonForExternalFS --- -charForRemainingFS - -codonForInternalFS --- -out_NT ./cds_aligned_macse/$(basename ${file} .macse.aln.fna).refined.macse.aln.fna -out_AA ./prot_aligned_macse/$(basename ${file} .macse.aln.fna).refined.macse.aln.faa

```

### Trim and clean MACSE alignments using the Alfix pipeline (HMMcleaner and MACSE)

```bash
# Download singularity image: 
https://www.agap-ge2pop.org/wp-content/uploads/macse/releases/MACSE_ALFIX_v01.sif

# Run the following in 3_hmmClean.sh:

#!/bin/bash
#SBATCH --job-name=hmmclean
#SBATCH --time=1-00:00:00
#SBATCH --partition=kamiak
#SBATCH --cpus-per-task=1
#SBATCH --mem=12GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=blperry@ucsc.edu
#SBATCH --array=1-14666:1 # Running on 2 to start
#SBATCH --output=./logs/JobArray_%A_%a.out
#SBATCH --error=./logs/JobArray_%A_%a.err

module load singularity
export SINGULARITY_BINDPATH="/path/to/scratch/dir/ortho_alignments/"

file=$(sed -n ''$SLURM_ARRAY_TASK_ID'p' macsealigns_refined.txt | awk '{print$1}')

singularity run ./MACSE_ALFIX_v01.sif --out_dir hmmcleaned_alignments --out_file_prefix $(basename ${file} .trim.macse.aln.fna)_hmmcleaned --in_seq_file $file --no_prefiltering --no_FS_detection --java_mem 10000m

```


