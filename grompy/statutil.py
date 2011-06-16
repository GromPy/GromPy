PCA_CAN_VIEW       = (1<<5)
PCA_CAN_BEGIN      = (1<<6)
PCA_CAN_END        = (1<<7)
PCA_CAN_DT         = (1<<14)
PCA_CAN_TIME       = (PCA_CAN_BEGIN | PCA_CAN_END | PCA_CAN_DT)
PCA_TIME_UNIT      = (1<<15)
PCA_KEEP_ARGS      = (1<<8)
PCA_SILENT         = (1<<9)
PCA_CAN_SET_DEFFNM = (1<<10)
PCA_NOEXIT_ON_ARGS = (1<<11)
PCA_QUIET          = (1<<12)
PCA_BE_NICE        = (1<<13)
#from statutil.h:
#/*****************************************************
# *         Some command line parsing routines
# *****************************************************/
#
##define PCA_CAN_VIEW       (1<<5)
#/* add option -w to view output files (must be implemented in program) */
##define PCA_CAN_BEGIN      (1<<6)
##define PCA_CAN_END        (1<<7)
##define PCA_CAN_DT         (1<<14)
##define PCA_CAN_TIME       (PCA_CAN_BEGIN | PCA_CAN_END | PCA_CAN_DT)
#/* adds options -b and -e for begin and end time for reading trajectories */
##define PCA_TIME_UNIT      (1<<15)
#/* set time unit for output */
##define PCA_KEEP_ARGS      (1<<8)
#/* keep parsed args in argv (doesn't make sense without NOEXIT_ON_ARGS) */
##define PCA_SILENT         (1<<9)
#/* don't print options by default */
##define PCA_CAN_SET_DEFFNM (1<<10)
#/* does something for non-master mdrun nodes */
##define PCA_NOEXIT_ON_ARGS (1<<11)
#/* no fatal_error when invalid options are encountered */
##define PCA_QUIET          (1<<12)
#/* does something for non-master mdrun nodes */
##define PCA_BE_NICE        (1<<13)
#/* Default to low priority, unless configured with --disable-nice */