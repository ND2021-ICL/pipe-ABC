cd $PBS_O_WORKDIR

OUTPUT_FILE="DEseq_max_"`date +"%Y%m%d-%H:%M"`".out"

module load python/3.8.5


if [ -z $PBS_NUM_NODES ]
  then
  echo "For threading correct usage: qsub -lnodes=1:ppn=x:msc script_name"
  exit 1
fi
if [ $PBS_NUM_NODES != "1" ]
  then
  echo "Number of Nodes "$PBS_NUM_NODES
  echo "For threading you must only specify 1 node e.g. qsub -lnodes=1:ppn=x:msc script_name"
  exit 2
fi

NP=$(( $PBS_NUM_NODES * $PBS_NUM_PPN ))
NUMBER_OF_THREADS=$PBS_NUM_PPN

DEBUG=0
if [ $DEBUG == 1 ]; then
  MYLOGFILE=debug.out
  hostname -f                                                            > $MYLOGFILE
  echo "PBS_NUM_NODES,PBS_NUM_PPN,NP : "$PBS_NUM_NODES,$PBS_NUM_PPN,$NP >> $MYLOGFILE
  echo "PBS_O_WORKDIR                : "$PBS_O_WORKDIR                  >> $MYLOGFILE
  echo "PBS_NODEFILE                 : "$PBS_NODEFILE                   >> $MYLOGFILE
  echo "Here is PBS_NODEFILE:"     >> $MYLOGFILE
  cat $PBS_NODEFILE                >> $MYLOGFILE
  echo "End of PBS_NODEFILE file"  >> $MYLOGFILE
fi

# --- The '2>&1' string redirects any errors into the output file

echo "Starting Program at : "`date`  > $OUTPUT_FILE 2>&1
echo "On node "`hostname`           >> $OUTPUT_FILE 2>&1
echo "Using Python: "`which python` >> $OUTPUT_FILE 2>&1
echo "Modules loaded:"              >> $OUTPUT_FILE 2>&1
module list                         >> $OUTPUT_FILE 2>&1
echo " "                            >> $OUTPUT_FILE 2>&1

python SMCABC_DESeq_ruben_4ODEs_maxscore.py $NUMBER_OF_THREADS $1 >> $OUTPUT_FILE 2>&1

echo " "                            >> $OUTPUT_FILE 2>&1
echo "Finishing Program at: "`date` >> $OUTPUT_FILE 2>&1

