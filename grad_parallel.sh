#!/bin/bash

#all arguments passed to the script
export args="$@"

scratchdir="scr"
#cmlfiles=$( ls $scratchdir/*/*.cml )
cmlfiles=$( ls $scratchdir/cmls/*.cml )

#echo "args are:"
#echo "$args"
#echo "cml files are:"
#echo "$cmlfiles"

#could also customize this to be job-specific
rundir=$( pwd )

#Load the gnu parallel module
module load parallel

#Create a file listing the nodes on which the job is running                              
uniq -c $PBS_NODEFILE | sed "s/$( hostname )/:/" | awk '{print $1"/"$2}' > $rundir/$PBS_JOBID.hosts

echo "$cmlfiles" | parallel --no-notice \
  --sshloginfile $rundir/$PBS_JOBID.hosts \
  --joblog $rundir/$PBS_JOBID.log \
  --env PATH \
  --env LD_LIBRARY_PATH \
  --env INCLUDE \
  --env PYTHONPATH \
  --env PYTHONUSERBASE \
  --env TMPFS \
  --env args \
  "cd $rundir;
  module load gcc cmake openblas/0.2.20 python/3.5.0;
  export OPENBLAS_NUM_THREADS=1;
  export MKL_NUM_THREADS=1;
  export PSI_SCRATCH=\"\$TMPFS/\$( printf %03d \$PARALLEL_SEQ )\";
  mkdir -p \$PSI_SCRATCH;
  [[ ! -z \$PSI_SCRATCH ]] && rm -rf \$PSI_SCRATCH/*;
  echo \"Process \$PARALLEL_SEQ on \$( hostname ) doing {}\";
  env > $rundir/env_\$PARALLEL_SEQ.txt;
  python Grad_Standalone.py {} \$args;"

