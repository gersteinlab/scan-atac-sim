#!/usr/bin/env python
# coding: utf-8

# Assuming the user is starting from bulk cell ATAC-seq, processed via standard ENCODE pipeline

# In[1]:


#-----import packages-----#

#common python packages
import numpy as np
import os
import pickle
import csv
import argparse
import wget
import multiprocessing
import pysam
import pandas as pd

#biological packages
import pybedtools
from pybedtools import featurefuncs
import pyBigWig


# In[2]:


#parsing command line arguments
# -----parsing command line arguments-----#
parser = argparse.ArgumentParser(description='preprocessing for SCAN-ATAC simulator')
parser.add_argument('-c', '--cell_types', type=str, help='comma separated string of cell_types')
parser.add_argument('-e', '--extend_peak_size', type=str, help='prevent background from within e-bp of peak')
parser.add_argument('-b', '--bin_size', type=str, help='resolution of background bins')
parser.add_argument('-i', '--in_dir', type=str, help='directory containing peak files')
parser.add_argument('-j', '--bam_dir', type=str, help='directory containing bam files')
parser.add_argument('-o', '--out_dir', type=str, help='output directory')


#simulate command line input
# cmdline_str='-c ' + " CLP,Ery,MEP,Mono,NKcell " + \
#     ' -i ' + "/gpfs/ysm/scratch60/gerstein/zc264/linkseq/scan-atac-sim/data/peak/" +\
#     ' -j ' + "/gpfs/ysm/scratch60/gerstein/zc264/linkseq/scan-atac-sim/data/bam/" + \
#     ' -o ' + "/gpfs/ysm/scratch60/gerstein/zc264/linkseq/scan-atac-sim/temp/" + \
#     ' -e ' + "1000" + \
#     ' -b ' + "1000"

#check if the files are there

#args = parser.parse_args(cmdline_str.split())
args = parser.parse_args()
args.cell_types = args.cell_types.split(",")

peak_file_list = []
bam_file_list = []
out_bam_list = []
for cell in args.cell_types:
    
    peak_file = args.in_dir + cell + ".narrowPeak"
    if not os.path.exists(peak_file):
        print(peak_file + " file does not exist")
        #exit(1)
    peak_file_list.append(peak_file)
    
    bam_file = args.bam_dir + cell + ".bam"
    if not os.path.exists(bam_file):
        print(bam_file + " file does not exist")
        #exit(1)
    bam_file_list.append(bam_file)
    
    out_bam_list.append(args.out_dir + cell + ".bam")
    
print("all files found!")


# In[3]:


#construct a set of autosome + X chromosome names
chromosomes = []
for i in range(1,23):
    chromosomes.append("chr"+str(i))
chromosomes.append("chrX")


# In[12]:


#step 1: read in, sort, and merge peak bed files
#step 2.1: merge peaks with -d 250 (within 250bp of each other)
for i in range(len(peak_file_list)):
    if i == 0:
        merged_peaks = pybedtools.BedTool(peak_file_list[i]).sort().merge()
    else:
        merged_peaks = merged_peaks.cat(pybedtools.BedTool(peak_file_list[i]).sort().merge(), 
                                        postmerge=False, force_truncate=False)

merged_peaks = merged_peaks.sort().merge(d=250).filter(lambda x: x.chrom in chromosomes)

#restore iterator for pybedtools as a datastream via outputting and re-reading
merged_peaks.saveas(args.out_dir + "merged_peaks.bed")
merged_peaks = pd.read_csv(args.out_dir + "merged_peaks.bed", sep="\t",header=None)
merged_peaks[3] = merged_peaks.index
merged_peaks[3] = 'p_' + merged_peaks[3].astype(str)
merged_peaks.to_csv(args.out_dir + "merged_peaks.bed", sep="\t",header=None, index=False)
merged_peaks = pybedtools.BedTool(args.out_dir + "merged_peaks.bed")

#step 2.2 make extended merged peak file
extended_merged_peaks = merged_peaks.slop(b=int(args.extend_peak_size), genome="hg38")
print("foreground peaks merged")


# In[6]:


#step 3: define background regions

#generate 2000bp of the entire genome
hg38_windows = pybedtools.BedTool().window_maker(genome="hg38", w=int(args.bin_size)).filter(pybedtools.featurefuncs.greater_than, int(args.bin_size)-1).filter(lambda x: x.chrom in chromosomes)

#remove ENCODE blacklist regions
if not os.path.exists('./hg38.blacklist.bed.gz'):
    url = 'http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg38-human/hg38.blacklist.bed.gz'
    wget.download(url, './hg38.blacklist.bed.gz')
blacklist = pybedtools.BedTool('./hg38.blacklist.bed.gz')
hg38_windows = (hg38_windows - blacklist).sort()

background = hg38_windows - extended_merged_peaks
background.saveas(args.out_dir + "background.bed")
background = pd.read_csv(args.out_dir + "background.bed", sep="\t",header=None)
background[3] = background.index
background[3] = 'b_' + background[3].astype(str)
background.to_csv(args.out_dir + "background.bed", sep="\t",header=None, index=False)
background = pybedtools.BedTool(args.out_dir + "background.bed")
print("background made")


# In[7]:


#step 4: get deduplicated and paired reads

bed_name_lst = []

for i in range(len(bam_file_list)):
    
    bam_file_name = bam_file_list[i]
    bam_sorted_name = out_bam_list[i].replace(".bam", ".sorted.bam")
    bed_out_name = out_bam_list[i].replace(".bam", ".bed")
    print(bam_file_name)
    
    #sort and index bam file
    pysam.sort("-@", str(multiprocessing.cpu_count()), "-o", bam_sorted_name, bam_file_name)
    pysam.index(bam_sorted_name)
    
    #fitler for paired and deduplciated reads
    samfile = pysam.AlignmentFile(bam_sorted_name, "rb")
    tsvfile = open(bed_out_name, 'w')
    writer = csv.writer(tsvfile, delimiter='\t')
    read_num = 0
    for read in samfile.fetch():
        if read.is_paired and not read.is_duplicate:

            read_name = "r_" + str(read_num) #read.query_name
            read_num = read_num + 1
            #if (read_name == None):
            #    read_name = "None"

            strand = '+'
            if(read.is_reverse) :
                strand = '-'
            writer.writerow([read.reference_name, read.reference_start, read.reference_end, read_name, 0, strand])

            
    samfile.close()
    tsvfile.close()
    
    bed_name_lst.append(bed_out_name)
    
print("bam filtered")


# In[8]:


#step 5: get peak stats
#step 6.1: merge background reads
for i in range(len(bed_name_lst)):
    merged_peaks.intersect(pybedtools.BedTool(bed_name_lst[i]), c=True).saveas(args.out_dir + args.cell_types[i]+".peak_counts.bed")
    merged_peaks.intersect(pybedtools.BedTool(bed_name_lst[i]), wa=True, wb=True).saveas(args.out_dir + args.cell_types[i]+".peak_intersect.bed")
    if i == 0:
        merged_reads = pybedtools.BedTool(bed_name_lst[i])
    else:
        merged_reads = merged_reads.cat(pybedtools.BedTool(bed_name_lst[i]), postmerge=False, force_truncate=False)
#merged_reads = merged_reads.sort()

print("merged reads")


# In[9]:


#step 6.1: merge background reads

# merge reads
# merged_reads = reads_bedOjb_list[0]
# for bed in reads_bedOjb_list[1:]:
#     merged_reads = merged_reads.cat(bed, postmerge=False, force_truncate=False)
# merged_reads = merged_reads.sort()


# In[10]:


# remove reads that overlap extended reads
bg_reads = merged_reads - extended_merged_peaks
bg_reads.saveas(args.out_dir + "bg_reads.bed")
bg_reads = pd.read_csv(args.out_dir + "bg_reads.bed", sep="\t",header=None)
bg_reads[3] = bg_reads[0] + ":" + bg_reads[1].astype(str) + ":" + bg_reads[2].astype(str)
bg_reads.to_csv(args.out_dir + "bg_reads.bed", sep="\t",header=None, index=False)

print("got background reads")


# In[11]:


#step 6.2: get background stats
mid_bg_reads = pybedtools.BedTool(args.out_dir + "bg_reads.bed").each(pybedtools.featurefuncs.midpoint)
background.intersect(mid_bg_reads, c=True).saveas(args.out_dir + "bg_counts.bed")

mid_bg_reads = pybedtools.BedTool(args.out_dir + "bg_reads.bed").each(pybedtools.featurefuncs.midpoint)
background.intersect(mid_bg_reads, wa=True, wb=True).saveas(args.out_dir + "bg_intersect.bed")

# return the intersect mid points to original position according to name
bg_reads = pd.read_csv(args.out_dir + "bg_intersect.bed", sep="\t",header=None)

bg_reads = pd.concat([bg_reads.iloc[:, 0:4], bg_reads[7].str.split(":",expand=True), pd.Series(bg_reads.index), bg_reads.iloc[:, 8:10]], axis=1)
bg_reads.to_csv(args.out_dir + "bg_intersect.expanded.bed", sep="\t",header=None, index=False)


# In[ ]:




