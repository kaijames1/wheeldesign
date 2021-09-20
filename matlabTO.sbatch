#!/bin/bash
###############################################################################
##                                                                           ##
##                           Campus Cluster                                  ##
##             Topology Optimization MATLAB Job Batch Script                 ##
##                                                                           ##
## SLURM Options (To view, run the following command below)                  ##
##                                                                           ##
##     man sbatch                                                            ##
##                                                                           ##
###############################################################################
#                                                                            
#SBATCH --time=11:30:00                  # Job run time (hh:mm:ss)
#SBATCH --nodes=1                        # Number of nodes
#SBATCH --ntasks-per-node=12             # Number of task (cores/ppn) per node
#SBATCH --job-name=A_TO                  # Name of batch job
#SBATCH --partition=eng-research         # Partition (queue)
#SBATCH --output=Zmatlab.o%j             # Name of batch job output file
##SBATCH --error=matlab.e%j              # Name of batch job error file
#SBATCH --mail-user=USERNAME@illinois.edu  # Send email notifications
#SBATCH --mail-type=BEGIN,END           # Type of email notifications to send
##SBATCH --dependency=afterany:2027905
#                                                                            
###############################################################################

# Change to the directory from which the batch job was submitted
# Note: SLURM defaults to running jobs in the directory where
# they are submitted, no need for cd'ing to $SLURM_SUBMIT_DIR

#cd ${SLURM_SUBMIT_DIR}

# Create the scratch directory
SCRATCH_DIR=/scratch/users/USERNAME
mkdir -p ${SCRATCH_DIR}
cd ${SCRATCH_DIR}


# Clear the value set in the DISPLAY environment variable
# to run the CLI version of MATLAB
unset DISPLAY

# Load MATLAB module (Enable MATLAB in user environment)
module load matlab/9.7

# Run matlab code and redirect output into a file whose name includes the JobID.
# matlab -nodisplay -r m-file-name >& ${SLURM_JOB_NAME}.${SLURM_JOB_ID}.oe
matlab -nodisplay -r MAIN_TO >& ${SLURM_JOB_NAME}.${SLURM_JOB_ID}.oe