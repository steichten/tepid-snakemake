#!/bin/bash
logdir=raijin/log
if [ -d $logdir ]
then
  pushd $logdir >/dev/null
  if [ -n "$(ls *.ER *.OU 2>/dev/null)" ]
  then
    tar cf `date +%y%m%d_%H%M%S`.tar *.OU *.ER 2>/dev/null
    rm *.OU *.ER
  fi
  popd >/dev/null
else
  mkdir -p $logdir
fi

# if [ -d ".kdmwrappers/.git" ]
# then
#     cd .kdmwrappers && git pull && cd ..
# else
#     git clone 'https://github.com/kdmurray91/snakemake-wrappers.git' .kdmwrappers
# fi

QSUB="qsub -q {cluster.queue} -l ncpus={threads} -l jobfs={cluster.jobfs}"
QSUB="$QSUB -l walltime={cluster.time} -l mem={cluster.mem}"
QSUB="$QSUB -l wd -o $logdir -e $logdir -P xe2"

snakemake                                \
    -j 3000                              \
    --cluster-config raijin/cluster.yaml \
    --js raijin/jobscript.sh             \
    --rerun-incomplete                   \
    --keep-going                         \
    --cluster "$QSUB" $@
