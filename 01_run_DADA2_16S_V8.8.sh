#!/bin/bash -   
#title          :01_run_DADA2_16S_V8.8.sh
#description    :Run DADA2 pipeline on a 16S MiSeq or MiniSeq run data
#author         :Roey Angel
#date           :20210317
#version        :8.8    
#usage          :go to the Analysis folder and run: 01_run_DADA2_16S_V8.8.sh $POOLING $DATAFOLDER $RESOURCES $GENE $MOCKID
#notes          :       
#bash_version   :5.0.17(1)-release
#============================================================================

HOMEFOLDER=$(pwd)
SCRIPTS=${HOMEFOLDER}
POOLING=${1:-pseudo} # pseudo, TRUE, FALSE
DATAFOLDER=${2:-"../Data/16S/noPrimers"} # directory containing the zipped fastq files after untaring
RESOURCES=${3:-"/proj/Resources/"} # RESOURCES should be referred to from $DATAFOLDER/DADA2_pseudo!
GENE=${4:-"16S"}
MOCKID=${5:-"mock"} # change to "" if no mock sample was used
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
mv -f ${DATAFOLDER}/DADA2_${POOLING} ${HOMEFOLDER}/
rm -r ${DATAFOLDER}/filtered

