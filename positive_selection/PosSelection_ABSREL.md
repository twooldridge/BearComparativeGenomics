## HyPhy - Analysis of Positive Selection
Working on Kamiak scratch at: `hyphy_bears`

Running [aBSREL](https://stevenweaver.github.io/hyphy-site/methods/selection-methods/#absrel), which tests if positive selection has occurred on a proportion of branches. It does not test for selection at specific sites. 

HyPhy v2.5.8 installed in Conda environment: `bwp_hyphy_env`

### Run ABSREL on all branches
```bash

# Run the following in 1_absrelTestAll.sh
1_absrelAll.sh:

#!/bin/bash
#SBATCH --job-name=hyphy_abs
#SBATCH --time=0-01:00:00
#SBATCH --partition=kamiak
#SBATCH --cpus-per-task=1
#SBATCH --mem=12GB
#SBATCH --array=1-14663:1
#SBATCH --output=./logs/JobArray_%A_%a.out
#SBATCH --error=./logs/JobArray_%A_%a.err

module load hyphy/2.5.42

file=$(sed -n ''$SLURM_ARRAY_TASK_ID'p' clean_aligns.txt | awk '{print$1}')

hyphy aBSREL --alignment $file --tree custom_BL_rooted.unroot.nwk --srv Yes --output ./all_absrel/$(basename ${file} .refined.macse.aln.fna_hmmcleaned_final_align_NT.aln).absrel.all.json > all_absrel_logs/$(basename ${file} .refined.macse.aln.fna_hmmcleaned_final_align_NT.aln).absrel.all.out.txt

# Count number significant
grep 'positive test results' ./all_absrel/* | awk -F ':|,' '$3 > 0 {print $0}'

# Pull out all significant genes/branches from outfiles
grep -A 20 'Holm-Bonferroni corrected _p' all_absrel_logs/* | grep -v 'Holm-Bonferroni c' | grep 'p-value' > allsig_absrel_ap26.txt

```

### Parse aBSREL results into CSVs
Using [phyphy](https://github.com/sjspielman/phyphy/tree/master?tab=readme-ov-file#extracting-csvs-from-hyphy-output-json) v0.4.3 python library. 

```bash
ls all_absrel > absrelResultsJsons.txt
mkdir all_absrel_csv

# Create toCSV.py:
import phyphy
import sys

infile = './all_absrel/'+sys.argv[1]
outfile = './all_absrel_csv/'+sys.argv[1]+'.csv'
e = phyphy.Extractor(infile) 
e.extract_csv(outfile)  


# Run the following in 4_resToCsv.sh
#!/bin/bash
#SBATCH --job-name=res_to_csv
#SBATCH --time=0-00:10:00
#SBATCH --partition=kamiak
#SBATCH --cpus-per-task=1
#SBATCH --mem=12GB
#SBATCH --array=1-14645:1
#SBATCH --output=./logs/JobArray_%A_%a.out
#SBATCH --error=./logs/JobArray_%A_%a.err

module load anaconda3
source activate bwp_phyphy_env

file=$(sed -n ''$SLURM_ARRAY_TASK_ID'p' absrelResultsJsons.txt | awk '{print$1}')

python toCSV.py ${file}

###

# Use awk to pull out all significant genes

awk -F "," -v OFS="," '$8 < 0.05 {print FILENAME,$0}' > absrel_sigResults_05.25.25.csv
```

