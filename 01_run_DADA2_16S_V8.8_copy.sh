#! /bin/bash
#$ -S /bin/bash
#$ -q bioinfo.q
#$ -V
#$ -cwd
#$ -N Run_DADA2
#$ -pe shared 20
#$ -e /gpfs0/biores/users/taliagab/SIP/err/) 
source activate /gpfs0/biores/users/taliagab/envs/
HOMEFOLDER=/gpfs0/biores/users/taliagab/SIP/
SCRIPTS=${HOMEFOLDER}
POOLING=${1:-FALSE} # pseudo, TRUE, FALSE
DATAFOLDER=${2:-"/gpfs0/biores/users/taliagab/SIP/noPrimers/"} # directory containing the zipped fastq files after untaring
RESOURCES=${3:-"/gpfs0/biores/users/taliagab/SIP/DADA2_FALSE/Resources/"} # RESOURCES should be referred to from $DATAFOLDER/DADA2_pseudo!
GENE=${4:-"16S"}
MOCKID=${5:-"Mock-com"} # change to "" if no mock sample was used
VERSION=V8.8
if [ -d $DATAFOLDER/DADA2_${POOLING} ]
then
 rm -rf $DATAFOLDER/DADA2_${POOLING}
fi
mkdir $DATAFOLDER/DADA2_${POOLING}

LOG=01_run_DADA2_${GENE}_${POOLING}_${VERSION}.log
touch ${HOMEFOLDER}/${LOG}

echo -e "Run DADA2 ${VERSION} on ${GENE} dataset with ${POOLING} option \n" >${HOMEFOLDER}/${LOG}
echo -e "Will process reads from the following folder:"  >> ${HOMEFOLDER}/${LOG}
echo -e $DATAFOLDER >> ${HOMEFOLDER}/${LOG}

echo -e "Starting R script \n" >> ${HOMEFOLDER}/${LOG}
cd $DATAFOLDER/DADA2_${POOLING}
/usr/bin/time -v Rscript --vanilla ${SCRIPTS}/01_DADA2_16S_merge_${VERSION}.R ../ $POOLING $RESOURCES $MOCKID >> ${HOMEFOLDER}/${LOG} 2>&1
cd ${HOMEFOLDER}
mv ${DATAFOLDER}/DADA2_${POOLING} ${HOMEFOLDER}/
rm -r ${DATAFOLDER}/filtered

