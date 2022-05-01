#!/bin/bash
#  REQUIRED TOOLS:
#    bwa 0.7.17-r1188
#    ivar 1.3.1
#    samtools 1.10

filename="list_samples_ids.txt"

while read line; do

echo $line

#====B0====#

echo "performing variant calling"

sampleName=$line
projectName="SARS-CoV-2"
dataDir="/data/SARS-CoV-2/raw_data/${projectName}/"
genomeFa="/data/SARS-CoV-2/reference/SARS-CoV-2-ANC.fasta"
primersBed="/data/SARS-CoV-2/primers/nCoV-2019_v3.bed"
resultsDir="/data/SARS-CoV-2/results/${projectName}/"
bamDir=${resultsDir}bamDir/
vcfDir=${resultsDir}vcfDir/
coverageDir=${resultsDir}coverage/
jobs=1

mkdir -p $resultsDir
if [ ! -d "$resultsDir" ]; then
    echo "Error mkdir"
    exit 1
fi

#====E0====#

#====B1====#

echo "bwa mem -- mapping reads to SARS-CoV-2 reference"

mkdir -p $bamDir
if [ ! -d "$bamDir" ]; then
    echo "Error mkdir"
    exit 1
fi

/path_to_conda/conda/envs/SARS-CoV-2/bin/bwa mem -t $jobs $genomeFa ${dataDir}${sampleName}_1.fastq.gz ${dataDir}${sampleName}_2.fastq.gz > ${bamDir}${sampleName}_aln.sam

#====E1====#

#====B2====#

echo "samtools -- building sorted bam"

/path_to_conda/conda/envs/SARS-CoV-2/bin/samtools view -b -F 4 -F 2048 -T $genomeFa ${bamDir}${sampleName}_aln.sam > ${bamDir}${sampleName}_aln.bam
/path_to_conda/conda/envs/SARS-CoV-2/bin/samtools sort ${bamDir}${sampleName}_aln.bam > ${bamDir}${sampleName}_aln.sorted.bam
rm ${bamDir}${sampleName}_aln.sam
rm ${bamDir}${sampleName}_aln.bam

#====E2====#

#====B3====#

echo "ivar trim -- trimming off the primer sequences"

/path_to_conda/conda/envs/SARS-CoV-2/bin/ivar trim -i ${bamDir}${sampleName}_aln.sorted.bam -b $primersBed -e -m 30 -q 20 -s 4 -p ${bamDir}${sampleName}_trimmed
rm ${bamDir}${sampleName}_aln.sorted.bam
rm ${bamDir}${sampleName}_aln.sorted.bam.bai

#====E3====#

#====B4====#

echo "samtools -- building and indexing trimmed sorted bam"

/path_to_conda/conda/envs/SARS-CoV-2/bin/samtools sort ${bamDir}${sampleName}_trimmed.bam > ${bamDir}${sampleName}.bam
/path_to_conda/conda/envs/SARS-CoV-2/bin/samtools index ${bamDir}${sampleName}.bam
rm ${bamDir}${sampleName}_trimmed.bam

#====E4====#

#====B5====#

echo "ivar variants -- calling variants"

mkdir -p $vcfDir
if [ ! -d "$vcfDir" ]; then
    echo "Error mkdir"
    exit 1
fi

/path_to_conda/conda/envs/SARS-CoV-2/bin/samtools mpileup -A -d 0 --reference $genomeFa -Q 0 ${bamDir}${sampleName}.bam | /path_to_conda/conda/envs/SARS-CoV-2/bin/ivar variants -p ${vcfDir}${sampleName} -q 20 -t 0.03

#====E5====#

#====B6====#

mkdir -p $coverageDir
if [ ! -d "$coverageDir" ]; then
    echo "Error mkdir"
    exit 1
fi

echo "samtools depth -- extracting coverage information"

/path_to_conda/conda/envs/SARS-CoV-2/bin/samtools depth -a ${bamDir}${sampleName}.bam > ${coverageDir}${sampleName}.txt

#====E6====#

done < $filename
