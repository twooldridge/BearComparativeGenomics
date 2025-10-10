## EGAPx annotation
The following code was run using the Washington State University Kamiak HPC. Specific paths have been generalized. 

### Download and set up EGAPx 
```bash
git clone https://github.com/ncbi/egapx.git # Version 0.3.2-alpha
cd egapx

idev # start interactive session
module load anaconda3

# set up virtual env
python -m venv ./egap_env
source ./egap_env/bin/activate
pip install -r requirements.txt

# Run example to generate config files that can be edited 
python3 ui/egapx.py ./examples/input_D_farinae_small.yaml -o example_out 
```

The `slurm.config` file generated from the example run was edited as follows to run on the WSU Kamiak HPC:
```bash
// Config for typical SLURM cluster
// Adjust the parameters for your cluster - queue name, temp dir, etc.
// Use your temp space instead of /data/$USER in the following config

singularity {
    enabled = true
    autoMounts = true
    cacheDir = "path/to/scratch/dir/ncbi_annot_pipeline/egapx/tmp/scratch/"
    envWhitelist='https_proxy,http_proxy,ftp_proxy,DISPLAY,SLURM_JOB_ID,SINGULARITY_BINDPATH'
}

env {
    SINGULARITY_CACHEDIR="path/to/scratch/dir/ncbi_annot_pipeline/egapx/tmp/singularity/"
    SINGULARITY_TMPDIR="path/to/scratch/dir/ncbi_annot_pipeline/egapx/tmp/tmp/"
}

process {
    resourceLimits = [
        memory: 120.GB,
        cpus: 16,
        time: 7.d
    ]
    executor = 'slurm'
    // Set other parameters specific to your cluser here
    maxRetries = 1
    queue = 'kamiak'
    queueSize = 200
    pollInterval = '2 min'
    queueStatInterval = '5 min'
    submitRateLimit = '6/1min'
    retry.maxAttempts = 1

    clusterOptions = ' --mem-per-cpu=8GB '

    scratch = 'path/to/scratch/dir/ncbi_annot_pipeline/egapx/tmp/'
    // with the default stageIn and stageOut settings using scratch can
    // result in humungous work folders
    // see https://github.com/nextflow-io/nextflow/issues/961 and
    //     https://www.nextflow.io/docs/latest/process.html?highlight=stageinmode
    stageInMode = 'symlink'
    stageOutMode = 'rsync'
}
```

## Initial run on sloth bear
Running a test using sloth bear to verify that the installation and configuration works. 

The parameters under `tasks` are based on [suggested modification to default parameters](https://github.com/ncbi/egapx?tab=readme-ov-file#modifying-default-parameters) when using RNA from other species (no RNA currently available for sloth bear). 

For `reads`, if we supply a SRA Study ID, then all of the data included under that study ID will be included. Using the following for this test:
- `SRP512493` - black bear transcriptome
- `SRP063482` - panda transcriptome
- `SRP373672` - brown bear transcriptome

### Input files are specified with a YAML file (`sloth_test.yaml`), shown below:
```bash
genome: path/to/scratch/dir/input_files/sloth.fasta.masked
taxid: 9636
reads:
  - SRP512493
  - SRP063482
  - SRP373672


tasks:
  rnaseq_collapse:
    rnaseq_collapse: -high-identity 0.8
  convert_from_bam:
    sam2asn: -filter 'pct_identity_gap >= 85'
  star_wnode:
    star_wnode: -pct-identity 85
```

### Run the following with `run_sloth.sh`:
```bash
#!/bin/bash
#SBATCH --partition=kamiak
#SBATCH --cpus-per-task=1
#SBATCH --mem=200G
#SBATCH --job-name=sloth_test        ### Job Name
#SBATCH --time=2-00:00:00       ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1               ### Node count required for the job
#SBATCH --out=sloth_test.out

module load anaconda3
source path/to/scratch/dir/ncbi_annot_pipeline/egapx/egap_env/bin/activate

module load nextflow
module load singularity

python3 ./ui/egapx.py ./sloth_test.yaml -e slurm -w sloth_working -o sloth_out

```

This ran successfully. The above `sloth.yaml` and `run_sloth.sh` files were modified to annotate each of the remaining seven species by replacing the input genome file and NCBI taxonomic ID. The same RNA reads and task parameters were used for all species. 

## Decontamination of the Asiatic black bear assembly

During the initial annotation of the Asiatic black bear assembly, EGAPx produced a warning about high levels of contamination. We used the NCBI Foreign Contamination Screen GX pipeline (https://github.com/ncbi/fcs) to screen and remove contamination using the code provided below. The cleaned assembly was then annotated using EGAPx. 

### Download and set up FCS-GX
```bash
# Download FCS-GX runner script and Singularity image
curl -LO https://github.com/ncbi/fcs/raw/main/dist/fcs.py

curl https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/releases/latest/fcs-gx.sif -Lo fcs-gx.sif
export FCS_DEFAULT_IMAGE=fcs-gx.sif

# Download FCS-GX database using the following slurm script:
#!/bin/bash
#SBATCH --job-name=fcs_db
#SBATCH --time=3-00:00:00
#SBATCH --partition=kamiak
#SBATCH --cpus-per-task=1
#SBATCH --mem=12GB
#SBATCH --out=slurm_out

module load singularity
module load anaconda3

export FCS_DEFAULT_IMAGE=fcs-gx.sif
SOURCE_DB_MANIFEST="https://ncbi-fcs-gx.s3.amazonaws.com/gxdb/latest/all.manifest"
LOCAL_DB="/path/to/scratch/dir/contamination_screen/db"
python3 fcs.py db get --mft "$SOURCE_DB_MANIFEST" --dir "$LOCAL_DB/gxdb"

```

### Screen the genome for contamination
```bash
#!/bin/bash
#SBATCH --job-name=fcs_screen
#SBATCH --time=3-00:00:00
#SBATCH --partition=kamiak
#SBATCH --cpus-per-task=1
#SBATCH --mem=600GB
#SBATCH --out=slurm_out

module load singularity
module load anaconda3

export FCS_DEFAULT_IMAGE=fcs-gx.sif

python3 ./fcs.py screen genome --fasta ../input_files/asiatic.fasta.masked --out-dir ./gx_out/ --gx-db ./db/gxdb --tax-id 9642
```

### Clean the genome
```bash
#!/bin/bash
#SBATCH --job-name=fcs_screen
#SBATCH --time=3-00:00:00
#SBATCH --partition=kamiak
#SBATCH --cpus-per-task=1
#SBATCH --mem=120GB
#SBATCH --out=slurm_out

module load singularity
module load anaconda3

#txid=9642

export FCS_DEFAULT_IMAGE=fcs-gx.sif

cat ../input_files/asiatic.fasta.masked | python3 ./fcs.py clean genome --action-report ./gx_out/asiatic.fasta.9642.fcs_gx_report.txt --output asiatic.fasta.masked.clean --contam-fasta-out asiatic.contam.fasta
```


