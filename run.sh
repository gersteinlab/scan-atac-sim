#!/bin/bash

# Copyright 2020 Flynn Chen
# Contact: flynn.chen at yale.edu


#test if there is intel compiler
if ! command -v icc &> /dev/null
then
    echo "intel c++ compiler could not be found"
    exit
fi


#test if there is python
if ! command -v python &> /dev/null
then
    echo "python could not be found"
    exit
fi

# Constants
variance=0.5 #default variance
cell_number=10000 #default cell number
signal_to_noise=0.7 #default signal to noise ratio
frag_num=3000
min_frag=1000 #default minimum fragment per cell
max_frag=20000 #default maxmum fragment per cell
temp_dir="./temp/"
out_dir="./output/"
cell_types=""
extend_peak_size=1000
bin_size=1000
peak_dir=""
bam_dir=""
preprocess_bam="yes"

usage="$(basename "$0") [-h] [-n -v -c -s -f -i -a -t -o -l -e -b -p -r] -- a program to simulate single-cell atac-seq from bulk-tissue atac-seq"

#parse more complex flags
while test $# -gt 0; do
  case "$1" in
    -h|--help)
      echo $usage
      exit
      ;;
    -n|--nopreprocess)
      shift
      preprocess_bam="no"
      shift
      ;;
    -v|--variance)
      shift
      variance="${1}"
      shift
      ;;
    -c|--cell_number)
      shift
      cell_number="${1}"
      shift
      ;;
    -s|--signal_to_noise)
      shift
      signal_to_noise="${1}"
      shift
      ;;
    -f|--frag_num)
      shift
      frag_num="${1}"
      shift
      ;;
    -i|--min_frag)
      shift
      min_frag="${1}"
      shift
      ;;
    -a|--max_frag)
      shift
      max_frag="${1}"
      shift
      ;;
    -t|--temp_dir)
      shift
      temp_dir="$1/"
      shift
      ;;
    -o|--out_dir)
      shift
      out_dir="${1}/"
      shift
      ;;
    -l|--cell_types)
      shift
      if test $# -gt 0; then
        cell_types=$1
      else
        echo "no cell types specified"
        exit 1
      fi
      shift
      ;;
    -e|--extend_peak_size)
      shift
      extend_peak_size=$1
      shift
      ;;
    -b|--bin_size)
      shift
      bin_size=$1
      shift
      ;;
    -p|--peak_dir)
      shift
      if test $# -gt 0; then
        peak_dir="$1/"
      else
        echo "no dir for peaks specified"
        exit 1
      fi
      shift
      ;;
    -r|--bam_dir)
      shift
      if test $# -gt 0; then
        bam_dir="$1/"
      else
        echo "no dir for bams specified"
        exit 1
      fi
      shift
      ;;
  esac
done

echo "Input Parameters:"
echo variance=$variance 
echo cell_number=$cell_number
echo signal_to_noise=$signal_to_noise
echo frag_num=$frag_num
echo min_frag=$min_frag
echo max_frag=$max_frag
echo temp_dir=$temp_dir
echo out_dir=$out_dir
echo cell_types=$cell_types
echo extend_peak_size=$extend_peak_size
echo bin_size=$bin_size
echo peak_dir=$peak_dir
echo bam_dir=$bam_dir

#make and clear temp_dor
mkdir -p $temp_dir
mkdir -p $out_dir

#compile under multicore environment
echo
echo "compiling simulation script in multi-cpu env"
make clean
make weighted_sampling
make uniform_sampling

if [ $preprocess_bam == "yes" ]; then
	# preprocessing
	echo
	echo "start preprocessing"
	python preprocess.py -c $cell_types -e $extend_peak_size -b $bin_size -i $peak_dir -j $bam_dir -o $temp_dir
	if [ $? != 0 ];
	then
	    echo "preprocessing failed"
	    exit 1
	fi
else
	echo
	echo "skip preprocessing"
fi

#iterate through each cell type
for cell_name in $(echo $cell_types | sed "s/,/ /g")
do

	# execute weighted sampling without replacement
	echo
	echo "start weighted sampling of peaks"
	./weighted_sampling -f ${temp_dir}/${cell_name}.peak_counts.bed -b ${temp_dir}/bg_counts.bed -of ${temp_dir}/${cell_name}.foreground.sampled.bed -ob ${temp_dir}/${cell_name}.background.sampled.bed -n $frag_num -nv $variance -c $cell_number -u -s $signal_to_noise -min $min_frag -max $max_frag
	if [ $? != 0 ];
	then
	    echo "weighted_sampling failed"
	    exit 1
	fi
	
	
	#perform uniform sampling of read from peaks
	rm -f ${temp_dir}/${cell_name}.bed
	echo
	echo "start uniform sampling of foreground reads"
	./uniform_sampling ${temp_dir}/${cell_name}.peak_intersect.bed ${temp_dir}/${cell_name}.foreground.sampled.bed ${temp_dir}/${cell_name}.foreground.bed
	if [ $? != 0 ];
	then
	    echo "uniform_sampling of foreground failed"
	    exit 1
	fi
	
	echo
	echo "start uniform sampling of background reads"
	./uniform_sampling ${temp_dir}/bg_intersect.expanded.bed ${temp_dir}/${cell_name}.background.sampled.bed ${temp_dir}/${cell_name}.background.bed
	if [ $? != 0 ];
	then
	    echo "uniform_sampling of background failed"
	    exit 1
	fi
	
	echo
	echo "combining generated foreground and background reads"
	cat ${temp_dir}/${cell_name}.foreground.bed ${temp_dir}/${cell_name}.background.bed > $out_dir/${cell_name}.SCAN-ATAC-Sim.bed

done

echo
echo "simulation finished"



