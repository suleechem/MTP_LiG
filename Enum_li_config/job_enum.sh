#!/bin/bash
#SBATCH -J Enum              # jobname
#SBATCH -p 128core             # partition name

#SBATCH -n 24                 # n of cores (= --mtp-workers 자동 감지값)
#SBATCH -N 1                 # n of nodes

#SBATCH -o %x.o%j
#SBATCH -e %x.e%j

export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

cd "$SLURM_SUBMIT_DIR"

# ── Conda 환경 확인 (mlip env 필요) ─────────────────────────────────
if [[ "${CONDA_DEFAULT_ENV}" != "Stage_Enum" ]]; then
    echo "ERROR: Please activate the Stage_Enum conda environment!"
    echo "Current env    : ${CONDA_DEFAULT_ENV:-<none>}"
    echo "Current python : $(which python)"
    exit 1
fi

# ── 실행 환경 정보 출력 ───────────────────────────────────────────────
echo "========================================"
echo "Job ID       : $SLURM_JOB_ID"
echo "Job name     : $SLURM_JOB_NAME"
echo "Node list    : $SLURM_JOB_NODELIST"
echo "Num nodes    : $SLURM_JOB_NUM_NODES"
echo "Num tasks    : $SLURM_NTASKS"
echo "Running host : $(hostname)"
echo "Allocated hosts:"
srun hostname | sort | uniq -c
echo "Python       : $(which python)"
python --version
echo "MLP path     : $(which mlp)"
timeout 3 mlp --version 2>/dev/null || true
echo "========================================"


python Enum_li_conf.py --structure LiC6.cif \
   --ga \
   --ga-min-li-li-same-layer 4.00 \
   --ga-pop 150 \
   --ga-gens 100 \
   --ga-expand-li-sites \
   --stacking-sequences ALL \
   --soc-stacking-policy conservative \
   --supercell 2 2 4 \
   --energy-backend mtp-relax \
   --mlip-ini mlip.ini \
   --mtp-launcher srun \
   --n-workers $SLURM_NTASKS \
   --mtp-workers 1 \
   --output OPTB88_AB

