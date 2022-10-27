#!/bin/bash
# a demo run script of AlphaPulldown example 1 (shorter)
# Usage:
# bash /repo/AlphaPulldown/run_pipeline.sh -b CYP736A167.fasta -B P450_baits.txt -c CPR_candidates.fasta  -C CPR_candidates.txt

# inspect run script location
AFP_REPO=$(readlink -f $(dirname $0));
# database dir
DATABASE_DIR="/mnt/db/"

# Path and user config
# change here if you wish to use a alternate version of any database.

mgnify_database_path="$DATABASE_DIR/mgnify/mgy_clusters.fa"
uniclust30_database_path="$DATABASE_DIR/uniref30_uc30/UniRef30_2020_06/UniRef30_2020_06"


usage() {
        echo ""
        echo "Usage: $0 <OPTIONS>"
        echo "Required Parameters:"
        # edited by Yinying
        echo "-b <baits>              Baits sequences"
        echo "-B <baits_info>         Baits info, one sequence description per line."
        echo "-c <candidates>         Candidates sequences"
        echo "-C <candidates_info>    Candidates info, one sequence description per line."
        echo "Optional Parameters:"
        echo "-o <save_dir>           Where to save the results"
        echo "-j <nproc>              number of parallel worker for MSA building"
        echo "-m <run_mode>           Run mode. pulldown<default>, all_vs_all, homo-oligomer, custom."
        echo ""
        exit 1
}

while getopts ":b:B:c:C:o:j:m:" i; do
        case "${i}" in

        b)
                baits=$OPTARG
        ;;
        B)
                baits_info=$OPTARG
        ;;
        c)
                candidates=$OPTARG
        ;;
        C)
                candidates_info=$OPTARG
        ;;
        j)
                nproc=$OPTARG
        ;;
        o)
                save_dir=$OPTARG
        ;;
        m)
                run_mode=$OPTARG
        ;;
        *)
                echo Unknown argument!
                usage
        ;;

        esac
done

if [[ ! -f "$baits" &&  ! -f "$baits_info" ]];then
  usage
fi


if [[ "$run_mode" != "" \
    &&"$run_mode" != "pulldown" \
    && "$run_mode" != "all_vs_all" \
    && "$run_mode" != "homo-oligomer" \
    && "$run_mode" != "custom"  ]];then
 usage
fi

if [[ "$save_dir" == "" ]];then
  save_dir=$PWD/save
fi

if [[ ! -p "$save_dir" ]];then
  mkdir -p $save_dir
fi

if [[ "$nproc" == "" ]];then
  nproc=32
fi

if [[ "$run_mode" == "" ]];then
  run_mode="pulldown"
fi

if [[ "$run_mode" == "pulldown" ]];then
  if [[ ! -f "$candidates" && ! -f "$candidates_info" ]];then
    echo "Your choose ${run_mode} as run_mode, yet the -c and -C is not available! "
    usage
  fi
fi

source activate AlphaPulldown

set -e

# alphafold pretrained data {sorted in subdirs named as release date}

echo  "run parallel MSA building step for features ..... "
cmd="python $AFP_REPO/alphapulldown/create_individual_features.py \
  --data_dir=${DATABASE_DIR} \
  --uniclust30_database_path=${uniclust30_database_path} \
  --fasta_paths=${baits},${candidates} \
  --max_template_date=2020-10-30 \
  --db_preset=full_dbs \
  --output_dir=${save_dir}/features \
  --mgnify_database_path=${mgnify_database_path} \
  --num_threads=${nproc} \
  --skip_existing "
echo "$cmd"
eval "$cmd"

echo 'run alphapulldown ...'
cmd="python $AFP_REPO/alphapulldown/run_multimer_jobs.py \
  --mode=${run_mode} \
  --num_cycle=3 \
  --num_predictions_per_model=1 \
  --output_path=${save_dir}/predict \
  --data_dir=${DATABASE_DIR} \
  --protein_lists=${baits_info},${candidates_info} \
  --monomer_objects_dir=${save_dir}/features "


echo "$cmd"
eval "$cmd"


echo "run analysis ..."
cmd="python $AFP_REPO/alphapulldown/analysis_pipeline/create_notebook.py \
  --cutoff=5.0 \
  --output_dir=${save_dir}/predict"

echo "$cmd"
eval "$cmd"

echo 'Done!'

echo "launch \`source activate AlphaPulldown; jupyter-lab ${save_dir}/output.ipynb --ServerApp.ip 0.0.0.0\` to enable remote analysis!"