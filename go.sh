#!/usr/bin/bash
source /cephfs/covid/software/eagle-owl/scripts/hootstrap.sh
source "$EAGLEOWL_CONF/common.sh"
source "$EAGLEOWL_CONF/asklepian/conf.sh"
source "$EAGLEOWL_CONF/ocarina/service_elan.sh"

# Activate env
eval "$(conda shell.bash hook)"
conda activate $CONDA_ASKLEPIAN

## Global environment variable override
## You should provide these somewhere else but you can override them here if you feel like it
### upload_azure.py
# AZURE_SAS= # key
# AZURE_END= # endpoint URL
### Asklepian
# ASKLEPIAN_DIR= # git root dir
# ASKLEPIAN_OUTDIR= # private working dir
# ASKLEPIAN_PUBDIR= # asklepian root dir for consortium consumption
### Runtime
# ELAN_DATE="YYYYMMDD" # elan dir date
# WUHAN_FP= # path to wuhan ref
# ARTIFACTS_ROOT= # consortium readable root elan publish dir

while read var; do
      [ -z "${!var}" ] && { echo 'Global Asklepian variable '$var' is empty or not set. Environment likely uninitialised. Aborting.'; exit 64; }
done << EOF
AZURE_SAS
AZURE_END
ASKLEPIAN_DIR
ASKLEPIAN_OUTDIR
ASKLEPIAN_PUBDIR
ELAN_DATE
WUHAN_FP
ARTIFACTS_ROOT
MAJORA_DOMAIN
EOF

set -euo pipefail

# Export variables needed by go_genome and go_variant
export DATESTAMP=${ELAN_DATE} # set by mqtt wrapper
export AZURE_SAS=$AZURE_SAS
export AZURE_END=$AZURE_END
export ASKLEPIAN_DIR=$ASKLEPIAN_DIR
export WUHAN_FP=$WUHAN_FP

# For testing, update OUTDIR, PUBDIR, PUBROOT and ensure to override the TABLE_BASENAMEs
OUTDIR="$ASKLEPIAN_OUTDIR/$DATESTAMP"
PUBDIR="$ASKLEPIAN_PUBDIR/$DATESTAMP"
PUBROOT="$ASKLEPIAN_PUBDIR"
GENOME_TABLE_BASENAME="v2_genome_table_$DATESTAMP"
VARIANT_TABLE_BASENAME="variant_table_$DATESTAMP"
LAST_BEST_REFS="$ASKLEPIAN_PUBDIR/latest/best_refs.paired.ls"

# Init outdir
mkdir -p $OUTDIR
cd $OUTDIR
SECONDS=0

# Get best ref for each central_sample_id
# TODO There is no need to do this de novo but we do so here for simplicity. We would do well to replace this as it will likely become quite slow, especially with IO perf on login node.
#      We can work toward some future system where Majora keeps a "symlink" to the best QC pag for a PAG Group.
# NOTE samstudio8 2021-01-29
#      Indeed this has become a little slow. In the interim, I have now made I
#      possible to query for consensus metrics directly from Ocarina rather than
#      indivudally hitting each QC file on the FS. We also use the consensus
#      matched FASTA to write out sequences, saving considerable time.
#        We also mark sequences in the paired.ls file to indicate if the best
#      sequence has changed today (1) or not (0). This will allow later functions
#      to decide whether or not to re-process the artifacts (e.g. for variants).

if [ ! -f "$OUTDIR/ocarina.ok" ]; then
    ocarina --oauth --env get pag --test-name 'cog-uk-elan-minimal-qc' --pass --task-wait --task-wait-attempts 60 \
        --ofield consensus.pc_masked pc_masked 'XXX' \
        --ofield consensus.pc_acgt pc_acgt 'XXX' \
        --ofield consensus.current_path fasta_path 'XXX' \
        --ofield consensus.num_bases num_bases 0 \
        --ofield central_sample_id central_sample_id 'XXX' \
        --ofield run_name run_name 'XXX' \
        --ofield published_name published_name 'XXX' \
        --ofield published_date published_date 'XXX' \
        --ofield adm1 adm1 'XXX' \
        --ofield collection_pillar collection_pillar '' \
        --ofield collection_date collection_date '' \
        --ofield received_date received_date '' > $OUTDIR/consensus.metrics.tsv
    touch $OUTDIR/ocarina.ok
else
    echo "[NOTE] Skipping ocarina"
fi
python -c "import datetime; print('ocarina', str(datetime.timedelta(seconds=$SECONDS)))"
SECONDS=0

if [ ! -f "$OUTDIR/best.ok" ]; then
    python $ASKLEPIAN_DIR/get_best_ref.py --fasta $ARTIFACTS_ROOT/elan/latest/elan.consensus.fasta --metrics $OUTDIR/consensus.metrics.tsv --latest $LAST_BEST_REFS --out-ls $OUTDIR/best_refs.paired.ls > $OUTDIR/best_refs.paired.fasta 2> $OUTDIR/best_refs.log
    touch $OUTDIR/best.ok
else
    echo "[NOTE] Skipping get_best_ref.py"
fi
python -c "import datetime; print('best', str(datetime.timedelta(seconds=$SECONDS)))"
SECONDS=0

################################################################################
# The 'naive' MSA here is not filtered on person-level identifiers, neither does
# it consider any filtering or masking. We use it in the knowledge that some
# of the variants may be dirty. Outside of its use for the full variant table here,
# it should not be considered an equivalent drop-in for the datapipe MSA.
# NOTE 20210729 datafunk was replaced with gofasta to significantly improve performance

if [ ! -f "$OUTDIR/msa.ok" ]; then
    minimap2 -t 24 -a -x asm5 --score-N=0 $WUHAN_FP $OUTDIR/best_refs.paired.fasta 2> $OUTDIR/mm2.log | gofasta sam tomultialign -t 24 --reference $WUHAN_FP -o $OUTDIR/naive_msa.fasta 2> $OUTDIR/gofasta.log
    touch $OUTDIR/msa.ok
else
    echo "[NOTE] Skipping MSA"
fi
python -c "import datetime; print('MSA', str(datetime.timedelta(seconds=$SECONDS)))"
SECONDS=0
################################################################################

# Start the genome and variant tables in parallel, parameters are: $1=WORKDIR $2=OUTDIR $3=TABLE_BASENAME
# Here we pass OUTDIR as both WORKDIR and OUTDIR for Asklep prod, but alternative OUTDIR can be used for testing
# Note flags are written to OUTDIR
$ASKLEPIAN_DIR/go_genome.sh $OUTDIR $OUTDIR $GENOME_TABLE_BASENAME &
$ASKLEPIAN_DIR/go_variant_v2.sh $OUTDIR $OUTDIR ${VARIANT_TABLE_BASENAME} &

# Wait for jobs
wait
python -c "import datetime; print('(wait)', str(datetime.timedelta(seconds=$SECONDS)))"
SECONDS=0

if [ ! -f "$OUTDIR/genome_upload2.ok" ] || [ ! -f "$OUTDIR/variant_upload2.ok" ]; then
    # EX_SOFTWARE
    exit 70
fi

# Make today's pubdir
mkdir -p $PUBDIR

# Clean up and push new artifacts
if [ ! -f "$OUTDIR/latest.ok" ]; then
    # Clean
    rm -f $OUTDIR/best_refs.paired.fasta
    rm -f $OUTDIR/${GENOME_TABLE_BASENAME}.csv.gz
    rm -f $OUTDIR/consensus.metrics.tsv

    # Push
    mv $OUTDIR/naive_msa.fasta $PUBDIR
    mv $OUTDIR/${VARIANT_TABLE_BASENAME}.csv $PUBDIR/naive_variant_table.csv
    mv $OUTDIR/best_refs.paired.ls $PUBDIR
    ln -fn -s $PUBDIR $PUBROOT/latest
    touch $OUTDIR/latest.ok
fi

# Remove yesterdays resources and repoint
if [ ! -f "$OUTDIR/head.ok" ]; then
    rm -f $PUBROOT/head/best_refs.paired.ls
    rm -f $PUBROOT/head/naive_msa.fasta
    rm -f $PUBROOT/head/naive_variant_table.csv
    ln -fn -s $PUBDIR $PUBROOT/head
    touch $OUTDIR/head.ok
fi

python -c "import datetime; print('finish', str(datetime.timedelta(seconds=$SECONDS)))"
SECONDS=0

