#!/usr/bin/bash
source /cephfs/covid/software/eagle-owl/scripts/hootstrap.sh
source "$EAGLEOWL_CONF/common.sh"
source "$EAGLEOWL_CONF/asklepian/conf.sh"

# Activate env
eval "$(conda shell.bash hook)"
conda activate $CONDA_ASKLEPIAN

while read var; do
      [ -z "${!var}" ] && { echo 'Global Asklepian variable '$var' is empty or not set. Environment likely uninitialised. Aborting go_genome.'; exit 64; }
done << EOF
AZURE_SAS
AZURE_END
ASKLEPIAN_DIR
WUHAN_FP
EOF

set -euo pipefail

WORKDIR=$1
OUTDIR=$2
TABLE_BASENAME=$3
SECONDS=0

# Make and push genome table
if [ ! -f "$OUTDIR/genome_table2.ok" ]; then
    python $ASKLEPIAN_DIR/make_genomes_table_v2.py --fasta $WORKDIR/naive_msa.fasta --meta $WORKDIR/consensus.metrics.tsv --best-ls $WORKDIR/best_refs.paired.ls | gzip --fast > $OUTDIR/${TABLE_BASENAME}.csv.gz
    touch $OUTDIR/genome_table2.ok
else
    echo "[NOTE] Skipping make_genomes_table (v2)"
fi
python -c "import datetime; print('make-genome', str(datetime.timedelta(seconds=$SECONDS)))"
SECONDS=0

if [ ! -f "$OUTDIR/genome_upload2.ok" ]; then
    python $ASKLEPIAN_DIR/upload_azure.py -c genomics -f $OUTDIR/${TABLE_BASENAME}.csv.gz
    touch $OUTDIR/genome_upload2.ok
else
    echo "[NOTE] Skipping genome upload (v2)"
fi
python -c "import datetime; print('push-genome', str(datetime.timedelta(seconds=$SECONDS)))"
SECONDS=0

