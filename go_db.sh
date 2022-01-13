#!/usr/bin/bash
source ~/.bootstrap.sh
source "$EAGLEOWL_CONF/envs.env"
source "$EAGLEOWL_CONF/asklepian.env"

# Activate env
eval "$(conda shell.bash hook)"
conda activate $CONDA_ASKLEPIAN

set -euo pipefail

cd $ASKLEPIAN_PUBDIR

# Get latest date
LAST_DIR_NAME=`readlink latest`
LAST_DIR_DATE=`basename $LAST_DIR_NAME`

# Init new database
touch asklepian.${LAST_DIR_DATE}.db
chmod 600 asklepian.${LAST_DIR_DATE}.db
sqlite3 asklepian.${LAST_DIR_DATE}.db < $ASKLEPIAN_DIR/sqlite3.cmd

# Get line count (probably a little time consuming)
LINE_COUNT=$(wc -l < /cephfs/covid/bham/results/variants/latest/naive_variant_table.csv)
LINE_COUNT=$((LINE_COUNT-1)) # Subtract header

# Get new db row count
ROW_COUNT=`echo 'SELECT COUNT(*) from variants;' | sqlite3 asklepian.${LAST_DIR_DATE}.db`

if [ "$LINE_COUNT" -eq "$ROW_COUNT" ]; then
    OLD_DB=`readlink -f asklepian.latest.db`

    # Update symlink
    chmod 644 $ASKLEPIAN_PUBDIR/asklepian.${LAST_DIR_DATE}.db
    ln -fn -s $ASKLEPIAN_PUBDIR/asklepian.${LAST_DIR_DATE}.db $ASKLEPIAN_PUBDIR/asklepian.latest.db

    # Remove previous database
    rm ${OLD_DB}
else
    exit 2
fi

