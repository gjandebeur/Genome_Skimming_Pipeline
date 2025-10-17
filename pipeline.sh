#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=nanopore_pipeline
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --gres=gpu:1

# ==================== USER INPUTS ====================
# Sample info
SAMPLE_NAME="sample1"

# Input/Output
INPUT_POD5_DIR="/path/to/input/pod5/"
OUTPUT_DIR="/path/to/output/${SAMPLE_NAME}/"

# Reference
REFERENCE_GENOME="/path/to/reference/GRCh38.fa"

# Software paths
DORADO_PATH="/path/to/dorado-1.0.1-linux-x64"
MINIMAP2_PATH="/path/to/minimap2"
MODKIT_PATH="/path/to/modkit"
APPTAINER_IMAGE="/scratch/$USER/apptainer_images/ubuntu20.sif"

# Basecalling parameters
MODEL="rna004_130bps_sup@v5.2.0"
MODIFICATIONS="sup,inosine_m6A_2OmeA,2OmeG,m5C_2OmeC,pseU_2OmeU"
MIN_QSCORE=10

# Alignment parameters
ALIGNMENT_TYPE="rna"  # Options: "rna" or "dna"
THREADS=8
# =====================================================

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Define output files
BASECALLED_BAM="${OUTPUT_DIR}/${SAMPLE_NAME}_basecalled.bam"
FASTQ_FILE="${OUTPUT_DIR}/${SAMPLE_NAME}_basecalled.fastq"
SAM_FILE="${OUTPUT_DIR}/${SAMPLE_NAME}_aligned.sam"
SORTED_BAM="${OUTPUT_DIR}/${SAMPLE_NAME}_aligned_sorted.bam"

echo "=========================================="
echo "Nanopore RNA-seq Pipeline"
echo "Sample: ${SAMPLE_NAME}"
echo "Output: ${OUTPUT_DIR}"
echo "=========================================="

# ==================== BASECALLING ====================
echo "[$(date)] Step 1: Basecalling..."

module load GCC PyTorch FlexiBLAS/3.3.1-GCC-12.3.0 FFmpeg/4.4.2-GCCcore-11.3.0 HTSlib protobuf

export APPTAINER_CACHEDIR=/scratch/$USER/apptainer_cache
export APPTAINER_TMPDIR=/scratch/$USER/apptainer_tmp
mkdir -p $APPTAINER_CACHEDIR $APPTAINER_TMPDIR

INPUT_POD5_MNT=$(echo ${INPUT_POD5_DIR} | sed 's|^/[^/]*|/mnt|')

apptainer exec --nv --bind /ourdisk:/mnt ${APPTAINER_IMAGE} \
    ${DORADO_PATH}/bin/dorado basecaller \
    --models-directory ${DORADO_PATH}/${MODEL} \
    ${MODIFICATIONS} \
    --min-qscore ${MIN_QSCORE} \
    --emit-moves \
    --estimate-poly-a \
    -r "${INPUT_POD5_MNT}" \
    > "${BASECALLED_BAM}"

echo "[$(date)] Basecalling complete!"

# ==================== ALIGNMENT ====================
echo "[$(date)] Step 2: Converting to FASTQ..."

module load SAMtools/1.16.1-GCC-11.3.0

samtools fastq -@ ${THREADS} -T "*" "${BASECALLED_BAM}" > "${FASTQ_FILE}"

echo "[$(date)] Step 3: Aligning..."

if [ "${ALIGNMENT_TYPE}" = "rna" ]; then
    MINIMAP_OPTS="-ax splice -y --secondary=no"
else
    MINIMAP_OPTS="-ax map-ont"
fi

${MINIMAP2_PATH} ${MINIMAP_OPTS} "${REFERENCE_GENOME}" "${FASTQ_FILE}" > "${SAM_FILE}"

echo "[$(date)] Step 4: Sorting..."

samtools sort -@ ${THREADS} -o "${SORTED_BAM}" "${SAM_FILE}"
samtools index "${SORTED_BAM}"

rm "${SAM_FILE}"

echo "[$(date)] Alignment complete!"

# ==================== STATS ====================
echo "[$(date)] Generating statistics..."

samtools flagstat "${SORTED_BAM}" > "${OUTPUT_DIR}/${SAMPLE_NAME}_stats.txt"

echo "=========================================="
echo "Pipeline Complete!"
echo "Basecalled: ${BASECALLED_BAM}"
echo "Aligned: ${SORTED_BAM}"
echo "Stats: ${OUTPUT_DIR}/${SAMPLE_NAME}_stats.txt"
echo "=========================================="
```
