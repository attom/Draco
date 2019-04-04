#!/bin/bash
##---------------------------------------------------------------------------##
## File  : regression/darwin-job-launch.sh
## Date  : Tuesday, May 31, 2016, 14:48 pm
## Author: Kelly Thompson
## Note  : Copyright (C) 2016, Triad National Security, LLC.
##         All rights are reserved.
##---------------------------------------------------------------------------##

# called from regression-master.sh
# assumes the following variables are defined in regression-master.sh:
#    $regdir     - /scratch/regress
#    $rscriptdir - /scratch/regress/draco/regression (actually, the location
#                  where the active regression_master.sh is located)
#    $subproj    - 'draco', 'jaynne', 'core', etc.
#    $build_type - 'Debug', 'Release'
#    $extra_params - '', 'intel13', 'vtest', 'coverage', 'arm'

# command line arguments
args=( "$@" )
nargs=${#args[@]}
scriptname=${0##*/}
host=`uname -n`

#export MOABHOMEDIR=/opt/MOAB
export SHOWQ=`which squeue`
export MSUB=`which sbatch`

# Dependencies: wait for these jobs to finish
dep_jobids=""
for (( i=0; i < $nargs ; ++i )); do
   dep_jobids="${dep_jobids} ${args[$i]} "
done

# load some common bash functions
export rscriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
if [[ -f $rscriptdir/scripts/common.sh ]]; then
  source $rscriptdir/scripts/common.sh
else
  echo " "
  echo "FATAL ERROR: Unable to locate Draco's bash functions: "
  echo "   looking for .../regression/scripts/common.sh"
  echo "   searched rscriptdir = $rscriptdir"
  exit 1
fi

# sanity check
job_launch_sanity_checks

#available_queues=`sacctmgr -np list assoc user=$LOGNAME | sed -e 's/.*|\(.*dev.*\)|.*/\1/' | sed -e 's/|.*//'`
#                = {darwin}
#available_account=`sacctmgr -np list assoc user=$LOGNAME format=Account | sed -e 's/|//' | sed -e 's/,/ /g' | xargs -n 1 | sort -u | xargs`
#                = {all_users asc-priority}
#available_partition=`sacctmgr -np list assoc user=$LOGNAME format=Partition | sed -e 's/|//' | sed -e 's/,/ /g' | xargs -n 1 | sort -u | xargs`
#                = { }
#available_qos=`sacctmgr -np list assoc user=$LOGNAME format=Qos | sed -e 's/|//' | sed -e 's/,/ /g' | xargs -n 1 | sort -u | xargs`
#                = {debug long normal scaling}

# case $available_queues in
#   *access*) access_queue="-A access --qos=access" ;;
#   *dev*)    access_queue="--qos=dev" ;;
# esac

# case $machine_name_short in
#   # Must use interactive qos on grizzly because standrad has a 70 node minimum
#   gr*) access_queue="--qos=interactive";;
# esac

# Banner
print_job_launch_banner

# Prerequisits:
# Wait for all dependencies to be met before creating a new job

for jobid in ${dep_jobids}; do
    while [ `ps --no-headers -u ${USER} -o pid | grep -c ${jobid}` -gt 0 ]; do
       echo "   ${subproj}: waiting for jobid = $jobid to finish (sleeping 5 minutes)."
       sleep 5m
    done
done

if ! test -d $logdir; then
  mkdir -p $logdir
  chgrp draco $logdir
  chmod g+rwX $logdir
  chmod g+s $logdir
fi

# choose partition type (knl, power9, volta, arm, etc)
case $extra_params_sort_safe in
  *knl*)       partition_options="-p knl-quad_cache" ;;
  *arm*)       partition_options="-p arm" ;;
  *power9*)    partition_options="-p power9-asc -A asc-priority" ;;
  *gpu-volta*) partition_options="-p volta-x86" ;;
esac

# Tell ctest to checkout jayenne from a local directory instead of gitlab. Also,
# update the local repository to include any feature branches that might be
# requested.
# export gitroot=/usr/projects/draco/regress/git

darwin_regress_script="${rscriptdir}/darwin-regress.msub"

# How long should we reserve the allocation for?
howlong="-t 1:00:00"
case $subproj in
  jayenne | core | trt | npt)  howlong="-t 8:00:00"
esac

##---------------------------------------------------------------------------##
# Proposed: Init, Configure, Build, Test and Submit
# (1) Configure, build and test from the backend
# (2) Submit results from the front end.
# Submit from the front end

# Configure, Build, Test on back end
echo " "
echo "Configure, Build, Test:"
export REGRESSION_PHASE=cbt
logfile=${logdir}/darwin-${subproj}-${build_type}${epdash}${extra_params}${prdash}${featurebranch}-${REGRESSION_PHASE}.log
cmd="${MSUB} -v -o ${logfile} -J${subproj:0:5}-${featurebranch} -N 1 ${howlong} ${partition_options} ${darwin_regress_script}"
echo "${cmd}"
jobid=`eval ${cmd}`
# trim extra whitespace from number
jobid=`echo ${jobid//[^0-9]/}`

# Wait for CBT (Config, build, test) to finish
sleep 1m
let elapsed_min=0
while [[ "`$SHOWQ | grep $jobid`" != "" ]]; do
   $SHOWQ | grep $jobid
   echo "   ${subproj}: waiting for jobid = $jobid to finish"
   echo "               we have waited $elapsed_min min. Sleeping another 5 min."
   sleep 5m
   let "elapsed_min += 5"
   # if we wait 12 hours, cancel the job.
   if [[ "${elapsed_min}" == "720" ]]; then
     scancel $jobid
   fi
done

# Submit from the front end
echo " "
echo "Submit:"
export REGRESSION_PHASE=s
echo "Jobs done, now submitting ${build_type} results from darwin."
cmd="${darwin_regress_script} >& ${logdir}/darwin-${subproj}-${build_type}${epdash}${extra_params}${prdash}${featurebranch}-${REGRESSION_PHASE}.log"
echo "${cmd}"
eval "${cmd}"

echo "All done."

##---------------------------------------------------------------------------##
## End of darwin-job-launch.sh
##---------------------------------------------------------------------------##
