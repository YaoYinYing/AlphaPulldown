#!/bin/bash
# a demo run script of AlphaPulldown example 1 (shorter)

source activate AlphaPulldown

set -e

export ALPHAFOLD_HOME='/repo/alphafold'

# inspect run script location
AFP_REPO=$(readlink -f $(dirname $0));
# database dir
DATABASE_DIR="/mnt/db/"

# Path and user config (change me if required)
mgnify_database_path="$DATABASE_DIR/mgnify/mgy_clusters.fa"


uniclust30_database_path="$DATABASE_DIR/uniref30_uc30/UniRef30_2020_06/UniRef30_2020_06"



# alphafold pretrained data {sorted in subdirs named as release date}
AF_PRETRAINED="/mnt/db/alphafold/"
AF_PRETRAINED_DATE='2022-03-02'

NPROC=32

TEST_DATA_BAITS=/mnt/data/yinying/tests/afpulldonw/baits.fasta
TEST_DATA_BAITS_INFO=/mnt/data/yinying/tests/afpulldonw/baits.txt

TEST_DATA_CANDIDATES=/mnt/data/yinying/tests/afpulldonw/example_1_sequences_shorter.fasta
TEST_DATA_CANDIDATES_INFO=/mnt/data/yinying/tests/afpulldonw/candidates_shorter.txt

# run parallel MSA building step for features
cmd="python $AFP_REPO/alphapulldown/create_individual_features.py \
  --data_dir=${DATABASE_DIR} \
  --uniclust30_database_path=${uniclust30_database_path} \
  --fasta_paths=${TEST_DATA_BAITS},${TEST_DATA_CANDIDATES} \
  --max_template_date=2020-10-30 \
  --db_preset=full_dbs \
  --output_dir=${PWD}/test \
  --mgnify_database_path=${mgnify_database_path} \
  --num_threads=${NPROC}"
echo "$cmd"
eval "$cmd"

# run alphapulldown
cmd="python $AFP_REPO/alphapulldown/run_multimer_jobs.py --mode=pulldown \
--num_cycle=3 \
--num_predictions_per_model=1 \
--output_path=${PWD}/test2 \
--data_dir=${DATABASE_DIR} \
--protein_lists=${TEST_DATA_BAITS_INFO},${TEST_DATA_CANDIDATES_INFO} \
--monomer_objects_dir=${PWD}/test \
--job_index=1"


echo "$cmd"
eval "$cmd"