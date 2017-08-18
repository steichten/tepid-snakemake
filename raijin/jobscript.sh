#!/bin/bash
# properties = {properties}

source /g/data1/xe2/.profile

set -ueo pipefail

export TMPDIR=$PBS_JOBFS

module load snakemake adapterremoval nextgenmap samtools bcftools freebayes bwa

{exec_job}
