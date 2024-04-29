#!/bin/bash

# This is a bash script of CA-align ver 1.0
# the following programs are required to be in the appropriate directories:
# CAMTX : Contact Area Caluculation
# CAalign : generate residue-residue contact area based alignment

#Please specify the directory of the CAB_ALIGN2
BIN_DIR='./'

USAGE(){
 echo -e "$0 [PDBfile1] [PDBfile2]"
 echo -e "Simple CAB-align"
 exit
}

CHK_FILE() {
 if [ ! -e "$1" ];then
  echo -e "No Such File or Directory. >> $1"
  exit
 fi
}

EXIST_CHK_FILE(){
 if [ -e "$1" ];then
  echo -e "$1 already exists"
  exit
 fi
}


if [ -z $2 ];then
 USAGE
fi

#Setting
exe1="$BIN_DIR/CAMTX_src"
exe2="$BIN_DIR/CAalign_src"
pdb1=$1
pdb2=$2

CHK_FILE "$pdb1"
CHK_FILE "$pdb2"
CHK_FILE "$exe1/CACDMTX"
CHK_FILE "$exe2/CAalign"

id1=`basename $pdb1 .pdb`
id1=`basename $id1 .ent`

id2=`basename $pdb2 .pdb`
id2=`basename $id2 .ent`

wk_dir=./tmp_cab_${id1}_${id2}

EXIST_CHK_FILE "$wk_dir"

mkdir -p $wk_dir

mtx1=$wk_dir/p1_$id1.mtx
mtx2=$wk_dir/p2_$id2.mtx

EXIST_CHK_FILE "$mtx1"
EXIST_CHK_FILE "$mtx2"
outfile=$wk_dir/${id1}_${id2}.align
EXIST_CHK_FILE "$outfile"


$exe1/CACDMTX -i $pdb1 > $mtx1
$exe1/CACDMTX -i $pdb2 > $mtx2

CHK_FILE "$mtx1"
CHK_FILE "$mtx2"

$exe2/CAalign -i $mtx1 -I $mtx2 -p $pdb1 -P $pdb2 > $wk_dir/${id1}_${id2}.align
cat $wk_dir/${id1}_${id2}.align


