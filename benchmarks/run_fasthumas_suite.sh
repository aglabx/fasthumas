#!/bin/bash
set -euo pipefail

BASE_DIR="/mnt/data/claude/benchmark_fasthumas"
RESULTS_DIR="${BASE_DIR}/results/fasthumas"
LOGS_DIR="${BASE_DIR}/logs"
BIN="/mnt/data/claude/2026-09-18_humas_hmmer/target/release/fasthumas"
HOR_HMM="/mnt/data/claude/2026-09-18_humas_hmmer/data/AS-HORs-hmmer3.3.2-120124.hmm"
SF_HMM="/mnt/data/claude/2026-09-18_humas_hmmer/data/AS-SFs-hmmer3.0.290621.hmm"
SUMMARY_FILE="${BASE_DIR}/fasthumas_summary.tsv"
THREADS=96

mkdir -p "${RESULTS_DIR}" "${LOGS_DIR}"

if [ ! -f "${SUMMARY_FILE}" ]; then
    echo -e "genome\tspecies\tassembly_path\twall_time\tuser_sec\tsys_sec\tmax_rss_kb\thor_sf_records\thor_records\tsf_records" > "${SUMMARY_FILE}"
fi

declare -A GENOMES
GENOMES["T2T-CHM13_v2.0"]="Homo sapiens|/mnt/data/t2t/CHM13/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna"
GENOMES["T2T-HG002_mat"]="Homo sapiens|/mnt/data/t2t/HG002/hg002v1.1.mat.fasta"
GENOMES["T2T-RPE1_hap1"]="Homo sapiens|/mnt/data/t2t/RPE1/RPE1v1.1.hap1.fasta"
GENOMES["Chimpanzee_mPanTro3"]="Pan troglodytes|/mnt/data/t2t/primates_t2t/chimpanzee/mPanTro3.cur.20231122.fasta"
GENOMES["Bonobo_mPanPan1"]="Pan paniscus|/mnt/data/t2t/primates_t2t/bonobo/mPanPan1.mat.cur.20231122.fasta"
GENOMES["Gorilla_mGorGor1"]="Gorilla gorilla|/mnt/data/t2t/primates_t2t/gorilla/mGorGor1.mat.cur.20231122.fasta"
GENOMES["Orangutan_mPonPyg2"]="Pongo pygmaeus|/mnt/data/t2t/primates_t2t/bornean_orangutan/mPonPyg2.hap1.cur.20231122.fasta"

ORDER=(
    "T2T-CHM13_v2.0"
    "T2T-HG002_mat"
    "T2T-RPE1_hap1"
    "Chimpanzee_mPanTro3"
    "Bonobo_mPanPan1"
    "Gorilla_mGorGor1"
    "Orangutan_mPonPyg2"
)

for NAME in "${ORDER[@]}"; do
    IFS="|" read -r SPECIES FASTA <<< "${GENOMES[$NAME]}"
    echo "=========================================================="
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting FastHumAS benchmark for ${NAME} (${SPECIES})"
    echo "FASTA: ${FASTA}"
    echo "=========================================================="

    PREFIX="${RESULTS_DIR}/${NAME}"
    TIMELOG="${LOGS_DIR}/${NAME}.time.log"
    PIPELOG="${LOGS_DIR}/${NAME}.pipeline.log"

    /usr/bin/time -v -o "${TIMELOG}" "${BIN}" \
        -i "${FASTA}" \
        -t "${THREADS}" \
        --hor "${HOR_HMM}" \
        --sf "${SF_HMM}" \
        -o "${PREFIX}" \
        > "${PIPELOG}" 2>&1 || true

    WALL_TIME=$(grep "Elapsed (wall clock) time" "${TIMELOG}" | awk -F': ' '{print $2}' || echo "N/A")
    USER_TIME=$(grep "User time (seconds)" "${TIMELOG}" | awk -F': ' '{print $2}' || echo "0")
    SYS_TIME=$(grep "System time (seconds)" "${TIMELOG}" | awk -F': ' '{print $2}' || echo "0")
    MAX_RSS=$(grep "Maximum resident set size" "${TIMELOG}" | awk -F': ' '{print $2}' || echo "0")

    HOR_SF_COUNT=$(wc -l < "${PREFIX}.AS-HOR+SF.bed" 2>/dev/null || echo "0")
    HOR_COUNT=$(wc -l < "${PREFIX}.AS-HOR.bed" 2>/dev/null || echo "0")
    SF_COUNT=$(wc -l < "${PREFIX}.AS-SF.bed" 2>/dev/null || echo "0")

    echo -e "${NAME}\t${SPECIES}\t${FASTA}\t${WALL_TIME}\t${USER_TIME}\t${SYS_TIME}\t${MAX_RSS}\t${HOR_SF_COUNT}\t${HOR_COUNT}\t${SF_COUNT}" >> "${SUMMARY_FILE}"

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished ${NAME}: Wall=${WALL_TIME}, Peak RSS=${MAX_RSS} KB, Records: HOR+SF=${HOR_SF_COUNT}, HOR=${HOR_COUNT}, SF=${SF_COUNT}"
done

echo "FastHumAS benchmark suite completed!"
