#!/usr/bin/bash
source /cephfs/covid/software/eagle-owl/scripts/hootstrap.sh
source "$EAGLEOWL_CONF/common.sh"
source "$EAGLEOWL_CONF/asklepian/conf.sh"

# Activate env
eval "$(conda shell.bash hook)"
conda activate $CONDA_ASKLEPIAN

while read var; do
      [ -z "${!var}" ] && { echo 'Global Asklepian variable '$var' is empty or not set. Environment likely uninitialised. Aborting go_variant.'; exit 64; }
done << EOF
AZURE_SAS
AZURE_END
ASKLEPIAN_DIR
WUHAN_FP
ASKLEPIAN_VARIANT_THREADS
EOF

set -euo pipefail

WORKDIR=$1
OUTDIR=$2
TABLE_BASENAME=$3
SECONDS=0

NUM_SEQUENCES=$(wc -l $WORKDIR/best_refs.paired.ls | cut -f1 -d' ')

# Make and push variant table
if [ ! -f "$OUTDIR/variant_table2.ok" ]; then
    python $ASKLEPIAN_DIR/make_variants_table_v2.py --ref $WUHAN_FP --msa $WORKDIR/naive_msa.fasta -n $NUM_SEQUENCES -t $ASKLEPIAN_VARIANT_THREADS > $OUTDIR/${TABLE_BASENAME}.csv
    touch $OUTDIR/variant_table2.ok
else
    echo "[NOTE] Skipping make_variants_table (v2)"
fi
python -c "import datetime; print('make-variant2', str(datetime.timedelta(seconds=$SECONDS)))"
SECONDS=0

if [ ! -f "$OUTDIR/variant_upload2.ok" ]; then
    python $ASKLEPIAN_DIR/upload_azure.py -c genomics -f $OUTDIR/${TABLE_BASENAME}.csv
    touch $OUTDIR/variant_upload2.ok
else
    echo "[NOTE] Skipping variant upload (v2)"
fi
python -c "import datetime; print('push-variant2', str(datetime.timedelta(seconds=$SECONDS)))"
SECONDS=0

