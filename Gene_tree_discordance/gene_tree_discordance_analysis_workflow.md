## Generate Gene windows and extract window alignments 
Bed file of windows of the maco chromosomes was made by:

```bash
halStats --chromSizes brown (or sloth) bears.hal > brown_bears.genome 
awk '$2 >= 20000000 || $1 == "NC_003427.1"' brown_bears.genome > brown_bears_filtered.genome # filters out chromosmes with <20M bp and keeps the mitochondria chromosome
bedtools makewindows -g brown_bears_filtered.genome -w 50000 (or 10000 or 10kb windows) > brown_bears_50kb_windows.bed # use bedtools to make a bed file with 50kb windows of the genome
```
Then for each window, using a custom python script, we identify overlapping reference regions from a bed file consisting of all regions that are covered by 6 or 7 species and calculate the total overlap length. 
If the overlap exceeds 90% of the window size, the window is retained. The final filtered set of high-coverage windows is saved as a new BED file for further analysis. 
This ensures that only windows with sufficient species representation are used in downstream analyses.
Script linked here

Filtered windows were then used to get alignments across the different species

```bash
# Load modules
module load cactus
module load phast
module load seqkit

# Define paths
HAL_FILE="/path/to/file/bears.hal"
GENOME_FILE="/path/to/file/brown_bears_filtered.genome"
BED_FILE="/path/to/file/brown_bears_50kb_windows.bed"
ALIGNMENT_DIR="/path/to/output/directory/1_alignments_50kb_brown"

# Read the specific line from BED file
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${GENOME_FILE})

chrom=$(echo "${LINE}" | awk '{print $1}')
size=$(echo "${LINE}" | awk '{print $2}')

mkdir -p "$ALIGNMENT_DIR/mafs"
mkdir -p "$ALIGNMENT_DIR/${chrom}"

awk -v chr="$chrom" '$1 == chr' "$BED_FILE" > "${ALIGNMENT_DIR}/${chrom}_regions.bed"

# Loop through BED file and extract MAF per region
while IFS=$'\t' read -r chrom start end; do
    region="${chrom}_${start}_${end}"
    echo "Processing $region..."

    # Create a temporary BED file for the current region
    TEMP_BED=$(mktemp)
    echo -e "$chrom\t$start\t$end" > "$TEMP_BED"

    # Run hal2maf using the temporary BED file
    hal2maf --refGenome brown --refTargets "$TEMP_BED" --noAncestors --noDupes "$HAL_FILE" "$ALIGNMENT_DIR/mafs/$region.maf"

    # Remove the temporary BED file
    rm -f "$TEMP_BED"

done < "${ALIGNMENT_DIR}/${chrom}_regions.bed"

# Convert each MAF file to FASTA
for maf in "$ALIGNMENT_DIR"/mafs/*.maf; do
    fasta_out="${maf%.maf}.fasta"
    msa_view --in-format MAF --out-format FASTA --unmask "$maf" | sed 's/\*/N/g' | seqkit grep -n -i -r -v -p 'Anc' | seqkit grep -n -i -r -v -p ':' > $fasta_out
done

mv "$ALIGNMENT_DIR"/mafs/*.fasta "$ALIGNMENT_DIR"/"$chrom"

```
## Infer Maximum likelihood gene trees and compute concordance factors
ML estimates best tree topology by finding tree that maximizes the likelihood of the observed sequence data given a substitution model
Concordance factors (CFs) measure proportion of gene trees that support a given branch (or split) in species tree.
- Provide insight into topological heterogeneity across genome.
CFs reflect biological variation due to processes like incomplete lineage sorting (ILS) or introgression.

```bash
module load iq-tree

ALN_DIR="/path/to/input/directory/1_alignments_50kb_brown"
TREE_DIR="/path/to/MLtree/output/directory/2_MLtrees_50kb_brown"
CF_DIR="/path/to/CF/output/directory/3_CF_50kb_brown"
REF_TREE="/path/to/file/cladogram.nwk"
genome_names="/path/to/file/brown_bears_filtered.genome"

# Read the specific line from genome file for this array task
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${genome_names})
chrom=$(echo "${LINE}" | awk '{print $1}')
size=$(echo "${LINE}" | awk '{print $2}')

# Find FASTA files for this chromosome
fasta_files=("$ALN_DIR"/"${chrom}"/"${chrom}"_*.fasta)


# Run IQ-TREE on FASTA files
for fasta in "${fasta_files[@]}"; do
    base=$(basename "$fasta" .fasta)

    if [[ -f "$TREE_DIR/${base}.fasta.treefile" ]]; then
        echo "Treefile ${base}.fasta.treefile already exists. Skipping."
        continue
    fi

    echo "Running IQ-TREE on: $fasta"
    iqtree2 -s "$fasta" -o panda -bb 1000 -m GTR+I+G --prefix "$TREE_DIR/${base}.fasta"
done

# Concatenate all treefiles for this chrom
if [[ -f "$CF_DIR/${chrom}.treefile" ]]; then
    echo "Concatenated treefile $CF_DIR/${chrom}.treefile already exists. Skipping."
else
    cat "$TREE_DIR"/"${chrom}"_*.treefile > "$CF_DIR/${chrom}.treefile"
fi

# Compute gene concordance factors
iqtree2 -t "$REF_TREE" --gcf "$CF_DIR/${chrom}.treefile" --prefix "$CF_DIR/${chrom}.concord"

# Compute site concordance factor
iqtree2 -te "$REF_TREE" -p "$ALN_DIR/$chrom" --scfl 100 -m GTR+I+G --prefix "$CF_DIR/${chrom}.concord2"

echo "Pipeline completed for ${chrom}!"
```
Gene concordance factor concored.cf.stat files and site concordance factor concord2.cf.stat files were concatenated and visualized in R.

## Get topology weighting
Assess genome-wide support for alternative topologies to the species tree.

First, we concatenated the tree files from maximum likelihood gene trees and filtered out trees that had less than 8 taxa
```bash
> 4_TWISST_50kb_brown/all_trees.iqtree.treefiles  # start with a clean file

for file in 2_MLtrees_50kb_brown/*.treefile; do
    cat "$file" >> 4_TWISST_50kb_brown/all_trees.iqtree.treefiles
done

awk '/panda/ && /asiatic/ && /sloth/ && /sun/ && /black/ && /brown/ && /polar/ && /spectacled/' 4_TWISST_50kb_brown/all_trees.iqtree.treefiles > 4_TWISST_50kb_brown/all_trees.8taxa.treefile

```

Using the filtered tree file, we ran TWISST
groups_brown.tsv is the species grouping:
asiatic	A
sloth	B
sun	B
black	C
brown	C
polar	C
spectacled	D
panda	Outgroup

```bash
# Load environment
module load miniconda3

# Create output directories if they don't exist
mkdir -p /path/to/output/4_TWISST_50kb_brown/output-weights
mkdir -p /path/to/output/4_TWISST_50kb_brown/output-topologies

python twisst-0.2/twisst.py -t 4_TWISST_50kb_sloth/all_UFBtrees_50kb_sloth.8taxa.treefile -w 4_TWISST_50kb_sloth/output-weights/50kb_sloth.trees.weights.csv.gz \
  --outputTopos 4_TWISST_50kb_brown/output-topologies/50kb_brown.all_trees.output.topologies.trees \
	--method complete \
	-g Outgroup -g A -g B -g C -g D \
	--outgroup Outgroup \
	--groupsFile 4_TWISST_50kb_brown/groups_brown.tsv
```
Output weights and topologies were visualized in R.
