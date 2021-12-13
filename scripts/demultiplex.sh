INDIR=$1
OUTDIR=$2

TAG_NAME=$3
TAG_SEQ=$4
RTAG_SEQ=$5
BATCHFILE=$6

mkdir -p $OUTDIR/DADA2_AS
mkdir -p $OUTDIR/DADA2_SS

while read INPUT_R1 INPUT_R2 PRIMER_F PRIMER_R MIN_LENGTH ; do

# Define binaries, temporary files and output files
CUTADAPT="$(which cutadapt) --discard-untrimmed --minimum-length ${MIN_LENGTH} -e 0"
CUTADAPT2="$(which cutadapt) -e 0"
VSEARCH=$(which vsearch)
TMP_FASTQ=$(mktemp)
MIN_F=$((${#PRIMER_F}))
MIN_R=$((${#PRIMER_R}))

REV_PRIMER_F="$(echo $PRIMER_F | rev | tr ATUGCYRKMBDHVN TAACGRYMKVHDBN)"
REV_PRIMER_R="$(echo $PRIMER_R | rev | tr ATUGCYRKMBDHVN TAACGRYMKVHDBN)"

rev="$(echo $primer | rev | tr ATUGCYRKMBDHVN TAACGRYMKVHDBN)"

FINAL_FASTQ="$OUTDIR/DADA2_SS/${TAG_NAME}_R1.fastq"

FTFP="$TAG_SEQ$PRIMER_F"
RTRP="$RTAG_SEQ$PRIMER_R"

# Trim tags, forward & reverse primers (search sense and antisense)
    cat "$INDIR/${INPUT_R1}" | \
    ${CUTADAPT} -g "${FTFP}" -e 0 -O "${#FTFP}" - | \
    ${CUTADAPT2} -a "${REV_PRIMER_R}" - > "${FINAL_FASTQ}"

FINAL_FASTQ="$OUTDIR/DADA2_AS/${TAG_NAME}_R2.fastq"

# Trim tags, forward & reverse primers (search sense and antisense)
    cat "$INDIR/${INPUT_R2}" | \
    ${CUTADAPT} -g "${FTFP}" -e 0 -O "${#FTFP}" - | \
    ${CUTADAPT2} -a "${REV_PRIMER_R}" - > "${FINAL_FASTQ}"

FINAL_FASTQ="$OUTDIR/DADA2_SS/${TAG_NAME}_R2.fastq"

# Trim tags, forward & reverse primers (search sense and antisense)
    cat "$INDIR/${INPUT_R2}" | \
    ${CUTADAPT} -g "${RTRP}" -e 0 -O "${#RTRP}" - | \
    ${CUTADAPT2} -a "${REV_PRIMER_F}" - > "${FINAL_FASTQ}"

FINAL_FASTQ="$OUTDIR/DADA2_AS/${TAG_NAME}_R1.fastq"

# Trim tags, forward & reverse primers (search sense and antisense)
    cat "$INDIR/${INPUT_R1}" | \
    ${CUTADAPT} -g "${RTRP}" -e 0 -O "${#RTRP}" - | \
    ${CUTADAPT2} -a "${REV_PRIMER_F}" - > "${FINAL_FASTQ}"

done < $INDIR/$BATCHFILE
