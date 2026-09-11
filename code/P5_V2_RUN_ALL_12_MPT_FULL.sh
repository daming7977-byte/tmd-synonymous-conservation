#!/usr/bin/env bash

set -euo pipefail

ROOT="${P5_PROJECT_ROOT:?Set P5_PROJECT_ROOT to the external project tree}"

MAPROOT="$ROOT/10_FORMAL_INDEX_BUILD"

FASTQ_DUMP="${P5_FASTQ_DUMP:?Set P5_FASTQ_DUMP to the installed executable}"

MAP_RUNTIME="$ROOT/08_EXECUTION_PREP/P5_V2_FORMAL_MAP210_RUNTIME.sh"

TOPHAT_RUNTIME="$MAPROOT/P5_V2_FORMAL_TOPHAT210_RUNTIME.sh"

MAPENV="${P5_MAP_ENV:?Set P5_MAP_ENV to the mapping environment}"

SAMTOOLS="$MAPENV/bin/samtools"

RRNA_INDEX="$MAPROOT/rRNA_index/formal_4rRNA"

GENOME_INDEX="$MAPROOT/hg38_index/hg38"

TRANSCRIPTOME_INDEX="$MAPROOT/transcriptome_index/known"

CUTENV="p5-rfoot-cutadapt41"

WORK="$PWD"

mkdir -p runs

exec > >(tee -a P5_V2_FULL_12_MPT_MASTER.log) 2>&1


fail () {

    echo
    echo "============================================================"
    echo "FULL_MPT_PROCESSING_GATE = RED"
    echo "FAILED_RUN = $1"
    echo "FAILED_STAGE = $2"
    echo "MPT_OUTCOME_REVEALED = FALSE"
    echo "============================================================"

    exit 1
}


echo "============================================================"
echo "P5-V2 FORMAL FULL MPT PROCESSING"
echo "============================================================"

echo "FULL_MPT_RAW_PROCESSING_AUTHORIZED = TRUE"
echo "MPT_OUTCOME_REVEALED = FALSE"

echo
echo "===== FORMAL TOOLCHAIN ====="

"$MAP_RUNTIME" bowtie2 --version 2>&1 | head -3
"$TOPHAT_RUNTIME" --version 2>&1
"$SAMTOOLS" --version 2>&1 | head -4

conda run \
  -n "$CUTENV" \
  cutadapt --version

"$FASTQ_DUMP" --version 2>&1 | head -5


# ============================================================
# VERIFIED EXTERNAL MANIFEST
# ============================================================

[ -s P5_V2_MPT_12RUN_MANIFEST.tsv ] || {
    echo "MANIFEST_GATE = RED"
    echo "REASON = MANIFEST_MISSING"
    exit 1
}

BAD_MANIFEST=$(
    awk -F '\t' '
        NF != 5 {
            n++
        }
        END {
            print n+0
        }
    ' P5_V2_MPT_12RUN_MANIFEST.tsv
)

if [ "$BAD_MANIFEST" -ne 0 ]; then
    echo "MANIFEST_GATE = RED"
    echo "BAD_MANIFEST_LINES = $BAD_MANIFEST"
    exit 1
fi

ROWS=$(
    awk 'END{print NR-1}' \
    P5_V2_MPT_12RUN_MANIFEST.tsv
)

if [ "$ROWS" -ne 12 ]; then
    echo "MANIFEST_GATE = RED"
    echo "MANIFEST_DATA_ROWS = $ROWS"
    exit 1
fi

echo "MANIFEST_GATE = GREEN"
echo "MANIFEST_DATA_ROWS = 12"



process_run () {

    RUN="$1"
    GSM="$2"
    FACTOR="$3"
    REP="$4"
    TYPE="$5"

    DIR="$WORK/runs/$RUN"

    RAW="$DIR/raw"
    TRIM="$DIR/trimmed"
    RRNA="$DIR/rrna"
    TOP="$DIR/tophat_out"
    FINAL="$DIR/final"
    LOG="$DIR/logs"

    mkdir -p \
      "$RAW" \
      "$TRIM" \
      "$RRNA" \
      "$FINAL" \
      "$LOG"

    if [ -e "$DIR/.DONE" ]; then

        echo
        echo "============================================================"
        echo "SKIP COMPLETED RUN = $RUN"
        echo "============================================================"

        return
    fi


    echo
    echo "============================================================"
    echo "RUN = $RUN"
    echo "GSM = $GSM"
    echo "FACTOR = $FACTOR"
    echo "REPLICATE = $REP"
    echo "TYPE = $TYPE"
    echo "============================================================"


    # ========================================================
    # 1. ENA DIRECT R1 ACQUISITION WITH TRUE RESUME
    # ========================================================

    echo
    echo "===== $RUN : ENA DIRECT R1 DOWNLOAD ====="

    # Never delete R1 .part here.
    # It is the resume checkpoint.
    rm -f \
      "$RAW/${RUN}_1.fastq" \
      "$RAW/${RUN}_2.fastq" \
      "$RAW/${RUN}_2.fastq.gz"

    ENA_REPORT="$LOG/01_ena_report.tsv"
    ENA_ENV="$LOG/01_ena_r1.env"
    ENA_LOG="$LOG/01_ena_download.log"

    if ! curl -fsSL \
      "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${RUN}&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes&format=tsv" \
      -o "$ENA_REPORT"
    then
        fail "$RUN" "ENA_METADATA_QUERY"
    fi

    python3 - \
      "$ENA_REPORT" \
      "$ENA_ENV" \
      <<'PYENA'
import sys
from pathlib import Path

report = Path(sys.argv[1])
env = Path(sys.argv[2])

lines = report.read_text().strip().splitlines()

if len(lines) != 2:
    raise SystemExit("ENA_REPORT_EXPECTED_ONE_DATA_ROW")

header = lines[0].split("\t")
values = lines[1].split("\t")

d = dict(zip(header, values))

urls = d["fastq_ftp"].split(";")
md5s = d["fastq_md5"].split(";")
sizes = d["fastq_bytes"].split(";")

if not (
    len(urls) ==
    len(md5s) ==
    len(sizes)
):
    raise SystemExit(
        "ENA_FASTQ_METADATA_LENGTH_MISMATCH"
    )

idx = None

for i, url in enumerate(urls):

    if url.endswith("_1.fastq.gz"):
        idx = i
        break

if idx is None:

    if len(urls) == 2:
        idx = 0
    else:
        raise SystemExit(
            "ENA_R1_COULD_NOT_BE_IDENTIFIED"
        )

url = urls[idx]

if not url.startswith(
    ("http://", "https://")
):
    url = "https://" + url

with env.open("w") as f:

    f.write(
        "ENA_R1_URL='%s'\n" % url
    )

    f.write(
        "ENA_R1_MD5='%s'\n" % md5s[idx]
    )

    f.write(
        "ENA_R1_BYTES='%s'\n" % sizes[idx]
    )

print("ENA_R1_URL =", url)
print("ENA_R1_MD5 =", md5s[idx])
print("ENA_R1_BYTES =", sizes[idx])
PYENA

    source "$ENA_ENV"

    R1="$RAW/${RUN}_1.fastq.gz"
    PART="$RAW/${RUN}_1.fastq.gz.part"

    echo "ENA_EXPECTED_R1_BYTES = $ENA_R1_BYTES"
    echo "ENA_EXPECTED_R1_MD5 = $ENA_R1_MD5"

    DOWNLOAD_REQUIRED=TRUE

    # --------------------------------------------------------
    # Already-complete R1?
    # --------------------------------------------------------

    if [ -s "$R1" ]; then

        EXISTING_BYTES=$(
          stat -f%z "$R1"
        )

        if [ "$EXISTING_BYTES" = "$ENA_R1_BYTES" ]; then

            EXISTING_MD5=$(
              md5 -q "$R1"
            )

            if [ "$EXISTING_MD5" = "$ENA_R1_MD5" ]; then

                echo "ENA_R1_ALREADY_COMPLETE = TRUE"

                DOWNLOAD_REQUIRED=FALSE

            else

                echo "EXISTING_COMPLETE_SIZE_BUT_MD5_MISMATCH = TRUE"

                # Preserve R1 until complete RUN_PROCESSING_GATE = GREEN.
            fi

        elif [ "$EXISTING_BYTES" -lt "$ENA_R1_BYTES" ]; then

            echo "MOVE_INCOMPLETE_R1_TO_RESUME_PART = TRUE"

            mv "$R1" "$PART"

        else

            echo "EXISTING_R1_TOO_LARGE = TRUE"

            rm -f "$R1"
        fi
    fi


    # --------------------------------------------------------
    # Download/resume
    # --------------------------------------------------------

    if [ "$DOWNLOAD_REQUIRED" = TRUE ]; then

        ATTEMPT=0

        while true
        do
            ATTEMPT=$((ATTEMPT + 1))

            if [ -e "$PART" ]; then
                CURRENT=$(stat -f%z "$PART")
            else
                CURRENT=0
            fi

            echo "ENA_DOWNLOAD_ATTEMPT = $ATTEMPT"
            echo "ENA_RESUME_FROM_BYTES = $CURRENT"

            if [ "$CURRENT" -eq "$ENA_R1_BYTES" ]; then
                echo "ENA_DOWNLOAD_BYTE_TARGET_REACHED = TRUE"
                break
            fi

            if [ "$CURRENT" -gt "$ENA_R1_BYTES" ]; then
                echo "ENA_PART_TOO_LARGE_RESET = TRUE"
                rm -f "$PART"
                CURRENT=0
            fi

            if [ "$CURRENT" -eq 0 ]; then
                curl \
                  --http1.1 \
                  -L \
                  --fail \
                  --connect-timeout 60 \
                  "$ENA_R1_URL" \
                  -o "$PART" \
                  >> "$ENA_LOG" 2>&1 \
                || true
            else
                curl \
                  --http1.1 \
                  -L \
                  --fail \
                  --connect-timeout 60 \
                  --continue-at - \
                  "$ENA_R1_URL" \
                  -o "$PART" \
                  >> "$ENA_LOG" 2>&1 \
                || true
            fi

            NEW=$(
              stat -f%z "$PART" 2>/dev/null \
              || echo 0
            )

            echo "ENA_CURRENT_BYTES_AFTER_ATTEMPT = $NEW"

            if [ "$NEW" -le "$CURRENT" ]; then
                sleep 20
            fi

            if [ "$ATTEMPT" -ge 30 ]; then
                fail "$RUN" "ENA_R1_DOWNLOAD_RETRY_LIMIT"
            fi
        done

        mv "$PART" "$R1"
    fi

    # --------------------------------------------------------
    # Hard archive integrity gate
    # --------------------------------------------------------

    [ -s "$R1" ] \
        || fail "$RUN" "ENA_R1_FILE_ABSENT"

    ACTUAL_BYTES=$(
      stat -f%z "$R1"
    )

    ACTUAL_MD5=$(
      md5 -q "$R1"
    )

    echo "ENA_ACTUAL_R1_BYTES = $ACTUAL_BYTES"
    echo "ENA_ACTUAL_R1_MD5 = $ACTUAL_MD5"

    [ "$ACTUAL_BYTES" = "$ENA_R1_BYTES" ] \
        || fail "$RUN" "ENA_R1_BYTE_COUNT"

    [ "$ACTUAL_MD5" = "$ENA_R1_MD5" ] \
        || fail "$RUN" "ENA_R1_MD5"

    echo "ENA_R1_ARCHIVE_INTEGRITY_GATE = GREEN"

    ls -lh "$R1"

    # ========================================================
    # 2. CUTADAPT
    # ========================================================

    echo
    echo "===== $RUN : CUTADAPT ====="

    TRIMMED="$TRIM/${RUN}.trimmed.fastq.gz"

    if ! conda run \
        -n "$CUTENV" \
        cutadapt \
        -a AAAAAAAA \
        -e 0.2 \
        -u 7 \
        -m 18 \
        -M 35 \
        -o "$TRIMMED" \
        "$R1" \
        > "$LOG/02_cutadapt.log" 2>&1
    then

        tail -100 "$LOG/02_cutadapt.log" || true

        fail "$RUN" "CUTADAPT"
    fi

    [ -s "$TRIMMED" ] \
        || fail "$RUN" "TRIMMED_FASTQ_ABSENT"

    gzip -t "$TRIMMED" \
        || fail "$RUN" "TRIMMED_GZIP_INTEGRITY"

    RAW_READS=$(
        grep -E \
          '^Total reads processed:' \
          "$LOG/02_cutadapt.log" \
        | tail -1 \
        | sed 's/.*:[[:space:]]*//' \
        | tr -d ', '
    )

    TRIMMED_READS=$(
        grep -E \
          '^Reads written \(passing filters\):' \
          "$LOG/02_cutadapt.log" \
        | tail -1 \
        | sed 's/.*:[[:space:]]*//' \
        | awk '{print $1}' \
        | tr -d ', '
    )

    echo "RAW_R1_READS = ${RAW_READS:-NA}"
    echo "CUTADAPT_RETAINED_READS = ${TRIMMED_READS:-NA}"

    rm -f "$R1"


    # ========================================================
    # 3. rRNA FILTER
    # ========================================================

    echo
    echo "===== $RUN : FOUR-rRNA FILTER ====="

    NONRRNA="$RRNA/${RUN}.non_rRNA.fastq"

    if ! "$MAP_RUNTIME" bowtie2 \
        --end-to-end \
        --very-sensitive \
        -a \
        --ignore-quals \
        --mp 6,6 \
        --np 6 \
        --rdg 100,100 \
        --rfg 100,100 \
        --score-min L,-12,0 \
        -x "$RRNA_INDEX" \
        -U "$TRIMMED" \
        --un "$NONRRNA" \
        -S /dev/null \
        2> "$LOG/03_rrna_bowtie2.log"
    then

        tail -100 "$LOG/03_rrna_bowtie2.log" || true

        fail "$RUN" "RRNA_BOWTIE2"
    fi

    [ -s "$NONRRNA" ] \
        || fail "$RUN" "NON_RRNA_FASTQ_ABSENT"

    NONRRNA_READS=$(
        awk 'END{print NR/4}' \
        "$NONRRNA"
    )

    echo "NON_RRNA_READS = $NONRRNA_READS"

    [ "$NONRRNA_READS" -gt 0 ] \
        || fail "$RUN" "NON_RRNA_ZERO"

    rm -f "$TRIMMED"


    # ========================================================
    # 4. TOPHAT
    # ========================================================

    echo
    echo "===== $RUN : TOPHAT ====="

    rm -rf "$TOP"

    if ! "$TOPHAT_RUNTIME" \
        -p 8 \
        --transcriptome-index="$TRANSCRIPTOME_INDEX" \
        -o "$TOP" \
        "$GENOME_INDEX" \
        "$NONRRNA" \
        > "$LOG/04_tophat.log" 2>&1
    then

        tail -150 "$LOG/04_tophat.log" || true

        fail "$RUN" "TOPHAT"
    fi

    BAM="$TOP/accepted_hits.bam"

    [ -s "$BAM" ] \
        || fail "$RUN" "ACCEPTED_HITS_BAM_ABSENT"

    ACCEPTED=$(
        "$SAMTOOLS" view \
          -c \
          "$BAM"
    )

    echo "TOPHAT_ACCEPTED_ALIGNMENTS = $ACCEPTED"

    [ "$ACCEPTED" -gt 0 ] \
        || fail "$RUN" "TOPHAT_ZERO_ALIGNMENTS"

    rm -f "$NONRRNA"


    # ========================================================
    # 5. UNIQUE NH:i:1
    # ========================================================

    echo
    echo "===== $RUN : UNIQUE NH:i:1 ====="

    UNIQUE_BAM="$FINAL/${RUN}.unique.bam"

    "$SAMTOOLS" view \
        -h \
        "$BAM" \
    | awk '
        /^@/ {
            print
            next
        }
        {
            keep=0
            for(i=12;i<=NF;i++){
                if($i=="NH:i:1"){
                    keep=1
                    break
                }
            }
            if(keep)
                print
        }
      ' \
    | "$SAMTOOLS" view \
        -S \
        -b \
        - \
    | "$SAMTOOLS" sort \
        -o "$UNIQUE_BAM"

    "$SAMTOOLS" index \
        "$UNIQUE_BAM"

    UNIQUE=$(
        "$SAMTOOLS" view \
          -c \
          "$UNIQUE_BAM"
    )

    echo "UNIQUE_MAPPED_READS = $UNIQUE"

    [ "$UNIQUE" -gt 0 ] \
        || fail "$RUN" "UNIQUE_ZERO"


    # ========================================================
    # 6. SELECT FOOTPRINT LENGTHS
    # ========================================================

    echo
    echo "===== $RUN : SELECT 18-20 / 26-29 ====="

    SELECTED_BAM="$FINAL/${RUN}.selected_lengths.bam"

    "$SAMTOOLS" view \
        -h \
        "$UNIQUE_BAM" \
    | awk '
        /^@/ {
            print
            next
        }

        length($10)==18 ||
        length($10)==19 ||
        length($10)==20 ||
        length($10)==26 ||
        length($10)==27 ||
        length($10)==28 ||
        length($10)==29 {
            print
        }
      ' \
    | "$SAMTOOLS" view \
        -S \
        -b \
        - \
    | "$SAMTOOLS" sort \
        -o "$SELECTED_BAM"

    "$SAMTOOLS" index \
        "$SELECTED_BAM"

    SELECTED=$(
        "$SAMTOOLS" view \
          -c \
          "$SELECTED_BAM"
    )

    echo "SELECTED_18_20_26_29_READS = $SELECTED"

    [ "$SELECTED" -gt 0 ] \
        || fail "$RUN" "SELECTED_ZERO"


    "$SAMTOOLS" view \
        "$UNIQUE_BAM" \
    | awk '
        {
            n=length($10)
            c[n]++
        }
        END {
            for(n in c)
                print n "\t" c[n]
        }
      ' \
    | sort -n \
    > "$FINAL/${RUN}.unique_length_distribution.tsv"


    # ========================================================
    # 7. +15 A-SITE
    # ========================================================

    echo
    echo "===== $RUN : +15 A-SITE ====="

    python3 - \
      "$SELECTED_BAM" \
      "$SAMTOOLS" \
      "$FINAL/${RUN}.asite_plus15.tsv" \
      "$FINAL/${RUN}.asite_summary.txt" \
      <<'PY'

import re
import subprocess
import sys

bam = sys.argv[1]
samtools = sys.argv[2]
out_path = sys.argv[3]
summary_path = sys.argv[4]

cigar_re = re.compile(
    r"(\d+)([MIDNSHP=X])"
)

p = subprocess.Popen(
    [samtools, "view", bam],
    stdout=subprocess.PIPE,
    universal_newlines=True,
)

total = 0
assigned = 0
failed = 0

with open(out_path, "w") as out:

    out.write(
        "qname\tchrom\tstrand\tread_length\tasite_0based\n"
    )

    for line in p.stdout:

        f = line.rstrip("\n").split("\t")

        if len(f) < 11:
            continue

        total += 1

        qname = f[0]
        flag = int(f[1])
        chrom = f[2]
        pos0 = int(f[3]) - 1
        cigar = f[5]
        seq = f[9]

        read_len = len(seq)

        reverse = bool(
            flag & 16
        )

        if reverse:

            target_q = (
                read_len - 1 - 15
            )

            strand = "-"

        else:

            target_q = 15
            strand = "+"

        q = 0
        r = pos0
        asite = None

        for ns, op in cigar_re.findall(
            cigar
        ):

            n = int(ns)

            if op in (
                "M",
                "=",
                "X",
            ):

                if (
                    q
                    <= target_q
                    < q + n
                ):

                    asite = (
                        r
                        + target_q
                        - q
                    )

                    break

                q += n
                r += n

            elif op in (
                "I",
                "S",
            ):

                if (
                    q
                    <= target_q
                    < q + n
                ):

                    break

                q += n

            elif op in (
                "D",
                "N",
            ):

                r += n

        if asite is None:

            failed += 1
            continue

        assigned += 1

        out.write(
            "%s\t%s\t%s\t%d\t%d\n"
            % (
                qname,
                chrom,
                strand,
                read_len,
                asite,
            )
        )

rc = p.wait()

if rc != 0:
    raise SystemExit(
        "samtools view failed"
    )

fraction = (
    float(assigned) / total
    if total
    else 0.0
)

with open(
    summary_path,
    "w",
) as fh:

    fh.write(
        "ASITE_INPUT_READS = %d\n"
        % total
    )

    fh.write(
        "ASITE_ASSIGNED_READS = %d\n"
        % assigned
    )

    fh.write(
        "ASITE_FAILED_READS = %d\n"
        % failed
    )

    fh.write(
        "ASITE_ASSIGNED_FRACTION = %.6f\n"
        % fraction
    )

print(
    "ASITE_INPUT_READS =",
    total
)

print(
    "ASITE_ASSIGNED_READS =",
    assigned
)

print(
    "ASITE_FAILED_READS =",
    failed
)

print(
    "ASITE_ASSIGNED_FRACTION = %.6f"
    % fraction
)

PY

    ASITE=$(
        awk \
          '$1=="ASITE_ASSIGNED_READS"{print $3}' \
          "$FINAL/${RUN}.asite_summary.txt"
    )

    echo "ASITE_ASSIGNED_READS = $ASITE"

    [ "${ASITE:-0}" -gt 0 ] \
        || fail "$RUN" "ASITE_ZERO"


    # ========================================================
    # 8. RUN QC SUMMARY
    # ========================================================

    cat > "$FINAL/${RUN}.QC.txt" <<EOF
RUN = $RUN
GSM = $GSM
FACTOR = $FACTOR
REPLICATE = $REP
TYPE = $TYPE

RAW_R1_READS = ${RAW_READS:-NA}
CUTADAPT_RETAINED_READS = ${TRIMMED_READS:-NA}
NON_RRNA_READS = $NONRRNA_READS
TOPHAT_ACCEPTED_ALIGNMENTS = $ACCEPTED
UNIQUE_MAPPED_READS = $UNIQUE
SELECTED_18_20_26_29_READS = $SELECTED
ASITE_ASSIGNED_READS = $ASITE

RUN_PROCESSING_GATE = GREEN
MPT_OUTCOME_REVEALED = FALSE
EOF


    # ========================================================
    # 9. CLEAN BULKY INTERMEDIATES
    # ========================================================

    if [ -s "$TOP/align_summary.txt" ]; then

        cp \
          "$TOP/align_summary.txt" \
          "$FINAL/${RUN}.tophat_align_summary.txt"
    fi

    rm -rf \
      "$RAW" \
      "$TRIM" \
      "$RRNA" \
      "$TOP"

    # selected_lengths BAM + A-site table are sufficient
    # for formal downstream outcome reconstruction.
    rm -f \
      "$UNIQUE_BAM" \
      "$UNIQUE_BAM.bai"

    touch "$DIR/.DONE"

    echo
    cat "$FINAL/${RUN}.QC.txt"

    echo
    echo "RUN_COMPLETE = $RUN"
}


tail -n +2 \
  P5_V2_MPT_12RUN_MANIFEST.tsv \
| while IFS=$'\t' read -r \
    RUN GSM FACTOR REP TYPE
do

    process_run \
      "$RUN" \
      "$GSM" \
      "$FACTOR" \
      "$REP" \
      "$TYPE"

done


echo
echo "============================================================"
echo "ALL 12 RUNS PROCESSED"
echo "============================================================"


echo
echo "===== FINAL 12-RUN QC TABLE ====="

{
    echo -e \
      "run\tfactor\treplicate\ttype\traw\ttrimmed\tnon_rrna\taccepted\tunique\tselected\tasite"

    tail -n +2 \
      P5_V2_MPT_12RUN_MANIFEST.tsv \
    | while IFS=$'\t' read -r \
        RUN GSM FACTOR REP TYPE
      do

        QC="runs/$RUN/final/${RUN}.QC.txt"

        RAW=$(
          awk '$1=="RAW_R1_READS"{print $3}' "$QC"
        )

        TRIM=$(
          awk '$1=="CUTADAPT_RETAINED_READS"{print $3}' "$QC"
        )

        NR=$(
          awk '$1=="NON_RRNA_READS"{print $3}' "$QC"
        )

        ACC=$(
          awk '$1=="TOPHAT_ACCEPTED_ALIGNMENTS"{print $3}' "$QC"
        )

        UNI=$(
          awk '$1=="UNIQUE_MAPPED_READS"{print $3}' "$QC"
        )

        SEL=$(
          awk '$1=="SELECTED_18_20_26_29_READS"{print $3}' "$QC"
        )

        AS=$(
          awk '$1=="ASITE_ASSIGNED_READS"{print $3}' "$QC"
        )

        echo -e \
          "$RUN\t$FACTOR\t$REP\t$TYPE\t$RAW\t$TRIM\t$NR\t$ACC\t$UNI\t$SEL\t$AS"

      done

} \
| tee P5_V2_FULL_12_MPT_QC.tsv


DONE=$(
  find runs \
    -name .DONE \
    | wc -l \
    | tr -d ' '
)

echo
echo "COMPLETED_RUNS = $DONE"

if [ "$DONE" -ne 12 ]; then

    fail \
      "MULTIPLE" \
      "NOT_ALL_12_COMPLETE"

fi


echo
echo "FULL_12_MPT_PROCESSING_GATE = GREEN"
echo "COMPLETED_RUNS = 12"
echo "MPT_OUTCOME_REVEALED = FALSE"
echo "NEXT_STEP = BUILD_CODON_LEVEL_MPT_OUTCOME_AND_REVEAL"
