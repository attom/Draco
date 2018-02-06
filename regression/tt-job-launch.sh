#!/bin/bash
##---------------------------------------------------------------------------##
## File  : regression/tt-job-launch.sh
## Date  : Tuesday, May 31, 2016, 14:48 pm
## Author: Kelly Thompson
## Note  : Copyright (C) 2016-2018, Los Alamos National Security, LLC.
##         All rights are reserved.
##---------------------------------------------------------------------------##

# called from regression-master.sh
# assumes the following variables are defined in regression-master.sh:
#    $regdir     - /scratch/regress
#    $rscriptdir - /scratch/regress/draco/regression (actually, the location
#                  where the active regression_master.sh is located)
#    $subproj    - 'draco', 'jaynne', 'capsaicin', etc.
#    $build_type - 'Debug', 'Release'
#    $extra_params - '', 'intel13', 'pgi', 'coverage'

# Under cron, a basic environment might not be loaded yet.
if [[ `which sbatch 2>/dev/null | grep -c sbatch` == 0 ]]; then
  source /etc/bash.bashrc.local
fi

# command line arguments
args=( "$@" )
nargs=${#args[@]}
scriptname=${0##*/}
host=`uname -n`

# import some bash functions
export rscriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $rscriptdir/scripts/common.sh

export SHOWQ=`which squeue`
export MSUB=`which sbatch`

 # Dependencies: wait for these jobs to finish
dep_jobids=""
for (( i=0; i < $nargs ; ++i )); do
   dep_jobids="${dep_jobids} ${args[$i]} "
done

# sanity checks
job_launch_sanity_checks

# Banner
print_job_launch_banner

# Prerequisits:
# Wait for all dependencies to be met before creating a new job

for jobid in ${dep_jobids}; do
    while [ `ps --no-headers -u ${USER} -o pid | grep -c ${jobid}` -gt 0 ]; do
       echo "   ${subproj}: waiting for jobid = $jobid to finish (sleeping 5 min)."
       sleep 5m
    done
done

  # srun options: (also see config/setupMPI.cmake)
  # -------------------------------------------------
  # -N        limit job to a single node.
  # --gres=craynetwork:0 This option allows more than one srun to be running at
  #           the same time on the Cray. There are 4 gres “tokens” available. If
  #           unspecified, each srun invocation will consume all of
  #           them. Setting the value to 0 means consume none and allow the user
  #           to run as many concurrent jobs as there are cores available on the
  #           node. This should only be specified on the salloc/sbatch command.
  #           Gabe doesn't recommend this option for regression testing.
  # --vm-overcommit=disable|enable Do not allow overcommit of heap resources.
  # -p knl    Limit allocation to KNL nodes.
  # -t N:NN:NN Length of allocation (e.g.: 8:00:00 hours).
# Select haswell or knl partition
# Optional: Use -C quad,flat to select KNL mode
# sinfo -o "%45n %30b %65f" | cut -b 47-120 | sort | uniq -c

# Note that we build on the haswell back-end (even when building code for knl):
build_partition_options="-N 1 -t 8:00:00 --gres=craynetwork:0"
case $extra_params in
knl) partition_options="-N 1 -t 8:00:00 --gres=craynetwork:0 -p knl" ;;
*)   partition_options="-N 1 -t 8:00:00 --gres=craynetwork:0" ;;
esac

# Configure on front end
# Only the front-end can see the github and gitlab repositories.
echo "Configure on the front end..."
export REGRESSION_PHASE=c
echo " "
cmd="${rscriptdir}/tt-regress.msub >& ${logdir}/${machine_name_short}-${subproj}-${build_type}${epdash}${extra_params}${prdash}${featurebranch}-${REGRESSION_PHASE}.log"
echo "${cmd}"
eval "${cmd}"


# Wait for Configure to finish before starting the build on the worker node:
# Build on the back-end to avoid loading up the shared front end.
echo " "
export REGRESSION_PHASE=b
echo "Build from the back end..."
echo " "
logfile=${logdir}/${machine_name_short}-${subproj}-${build_type}${epdash}${extra_params}${prdash}${featurebranch}-${REGRESSION_PHASE}.log
cmd="$MSUB -o ${logfile} -J ${subproj:0:5}-${featurebranch} ${build_partition_options} ${rscriptdir}/tt-regress.msub"
echo "${cmd}"
jobid=`eval ${cmd}`
jobid=`echo $jobid | sed -e 's/.*[ ]//'`
echo "jobid = ${jobid}"

# Wait for testing to finish
sleep 1m
while test "`${SHOWQ} | grep $jobid`" != ""; do
   ${SHOWQ} | grep $jobid
   sleep 5m
done

# Wait for the Build to finish before starting the testing from the worker node:
# Test from a different back end.  This step is separate from the build because
# we build KNL binaries from a Haswell back-end, but the tests must be run from
# a KNL node.
echo " "
export REGRESSION_PHASE=t
echo "Test from the back end..."
echo " "
logfile=${logdir}/${machine_name_short}-${subproj}-${build_type}${epdash}${extra_params}${prdash}${featurebranch}-${REGRESSION_PHASE}.log
cmd="$MSUB -o ${logfile} -J ${subproj:0:5}-${featurebranch} ${partition_options} ${rscriptdir}/tt-regress.msub"
echo "${cmd}"
jobid=`eval ${cmd}`
# delete blank lines
#jobid=`echo $jobid | sed '/^$/d'`
# only keep the job number
jobid=`echo $jobid | sed -e 's/.*[ ]//'`
echo "jobid = ${jobid}"

# Wait for testing to finish
sleep 1m
while test "`${SHOWQ} | grep $jobid`" != ""; do
   ${SHOWQ} | grep $jobid
   sleep 5m
done

# Submit from the front end.  Only the front end can see our cdash server.
echo " "
echo "Submit:"
export REGRESSION_PHASE=s
echo "- jobs done, now submitting ${build_type} results from tt-fey."
cmd="${rscriptdir}/tt-regress.msub >& ${logdir}/${machine_name_short}-${subproj}-${build_type}${epdash}${extra_params}${prdash}${featurebranch}-${REGRESSION_PHASE}.log"
echo "${cmd}"
eval "${cmd}"

# Submit from the front end
echo "Jobs done."

##---------------------------------------------------------------------------##
## End of tt-job-launch.sh
##---------------------------------------------------------------------------##
