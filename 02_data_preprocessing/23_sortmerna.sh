#!/bin/bash
#SBATCH --job-name=23_sortmerna
#SBATCH --partition=fast
#SBATCH --array=1-16
#SBATCH --cpus-per-task=8
#SBATCH --mem=20GB
#SBATCH --output=%x.%A_%a.stdout
#SBATCH --error=%x.%A_%a.stderr


# modules/dependencies loading
module load sortmerna

# Main script
READIN1=$(ls *.cutadapt.trim_1P.fq.gz | sed -n ${SLURM_ARRAY_TASK_ID}p)
READIN2=${READIN1/_1P/_2P}

READOUT=${READIN1/.cutadapt.trim_1P.fq.gz/.cutadapt.trim.non_rrna}
WORKDIR=${READIN1/cutadapt.trim_1P.fq.gz/.cutadapt.trim.sortmerna}

rm -rf $WORKDIR/kvdb

srun sortmerna --ref smr_v4.3_default_db.fasta --reads $READIN1 --reads $READIN2 --workdir $WORKDIR --fastx --paired_out --other ./$READOUT --threads ${SLURM_CPUS_PER_TASK}

srun gzip ${READOUT}.fq

mv ${WORKDIR}/out/aligned.log ${WORKDIR}.log
