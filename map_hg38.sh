#!/bin/bash

#SYSTEM REQS FOR PYTHON SCRIPTS
#sudo apt install bedtools
#sudo apt install python3
#sudo apt install python3-dev
#sudo apt install python3-pip
#sudo apt install pyvcf
#sudo apt install python3-pysam
#sudo apt install python3-pyvcf
#sudo apt install default-jre (for haplogrep)

# required temporary directory for samtools sort:
SORTDIR='/usr/local/geospiza/var/tmp/'



START=$(date +%s.%N)
clear
# setup parameters

YSEQID=${PWD##*/}
# YSEQID="1234" # (the above command simply gets the name of the last segment of the current working directory)

my_link="fastq"
if [ -L ${my_link} ] ; then
	if [ -e ${my_link} ] ; then
		echo "Good symlink ${my_link} already exists. Not searching for data."
	else
		echo "Broken link ${my_link}"
		exit 1
	fi
elif [ -e ${my_link} ] ; then
	echo "Note: ${my_link} is not a symlink! Using actual directory ${my_link}"
else
	echo "Missing symlink ${my_link}. Searching for fastq data of YSEQ ID ${YSEQID}..."
	ARCHDIR=$(find /genomes/arch* -type d -name "${YSEQID}" 2>/dev/null)
	FQDIR="${ARCHDIR}/fastq"
	if [ -d "$FQDIR" ];
	then
    	echo "Found fastq directory for ${YSEQID} at ${FQDIR}"
    	ln -s ${FQDIR} fastq
	else
		echo "/genomes/arch*/${YSEQID}/fastq directory does not exist."
		exit 1
	fi
fi

# Check file size in bytes
need=$(ls -la fastq/ | awk '{sum+=$5;} END {printf "%i", sum;}')

# Check for disk space in bytes
avail=$(df --block-size=1 $(pwd -P) | awk 'NR==2 {print $4}')

# Calculate available space after accounting for overhead in kilobytes
free=$(($avail - $need ))

# Convert free space to gigabytes
free_gb=$(awk "BEGIN {printf \"%.2f\", $free / 1024 / 1024 /1000 }")

if [ "$need" -gt "$avail" ]; then
	echo "Insufficient space on /work: ${free_gb} GB missing. Stopping script."
	exit 1
else
	echo "After mapping, approximately ${free_gb} GB left on /$(pwd -P | cut -d'/' -f2)."
fi


NUM_THREADS=$(getconf _NPROCESSORS_ONLN)
echo "We can use ${NUM_THREADS} threads."

PYTHON="python3"

BWA='bwa mem'
#if test -x "/usr/local/bin/bwa-mem2"; then
#    echo "We can use bwa-mem2"
#    BWA='bwa-mem2'
#fi


# Create a random password
#password creation omitted for securety reasons
PASSWD="12345"
echo "Random password created: ${PASSWD}"


# Generate summary file
cat /genomes/0/script-templates/x_result_summary.txt | sed "s/<kitnumber>/${YSEQID}/g" | sed "s/<passwd>/${PASSWD}/g" >${YSEQID}_result_summary.txt
echo "${YSEQID}_result_summary.txt created"

REF="/genomes/0/refseq/hg38/hg38.fa"
REFURL="http://ybrowse.org/WGS/ref/hg38/hg38.fa"

READS=""
BWA_PARAM="-M -t $NUM_THREADS"

READS_1="fastq/${YSEQID}_R1.fastq.gz"
READS_2="fastq/${YSEQID}_R2.fastq.gz"
READS_NANOPORE="fastq/${YSEQID}_nanopore.fastq.gz"
READS_400SE="fastq/${YSEQID}_400SE.fastq.gz"
READS_SANGER="fastq/${YSEQID}_Sanger.fastq"

if test -f "$READS_NANOPORE"; then
    echo "Found Nanopore FASTQ file. Using standard minimap2."
    READS=$READS_NANOPORE
    BWA="minimap2"
    BWA_PARAM="-ax map-ont -t ${NUM_THREADS}"
    REF="/genomes/0/refseq/hg38/hg38.fa.mmi"
fi

if test -f "$READS_SANGER"; then
    echo "Found Sanger FASTQ file. Using standard bwa mem."
    READS=$READS_SANGER
    BWA="bwa mem"
fi

if test -f "$READS_400SE"; then
    echo "Found 400SE FASTQ file"
    READS=$READS_400SE
    BWA="bwa mem"
fi

if test -f "$READS_1"; then
    echo "Found R1 FASTQ file"
    READS=$READS_1
fi

if test -f "$READS_2"; then
    echo "Found R2 FASTQ file"
    READS="$READS_1 $READS_2"
fi

# Place a file in the fastq directory to indicate that the mapping is being done
touch fastq/$USER.$HOSTNAME.work

BAMFILE="${YSEQID}_bwa-mem_hg38.bam"

BAMFILE_SORTED="${YSEQID}_bwa-mem_hg38_sorted.bam"
VCF_FILE="${YSEQID}_hg38.vcf"

cp /genomes/0/script-templates/addRows.pl .
cp /genomes/0/script-templates/multiply.pl .

cp /genomes/0/tree/trees/yfull/latest_YFull_YTree.json .
YFULLTREE="latest_YFull_YTree.json"

cp /genomes/0/script-templates/cladeFinder.py .
CLADEFINDER="cladeFinder.py"

#ANN_VCF_FILE="${YSEQID}_hg38_dbSNP150_ann.vcf"
#ANN_CLINVAR_VCF_FILE="${YSEQID}_hg38_clinvar_ann.vcf"

#DBSNP="/genomes/0/refseq/hg19/All_chr_20161121.vcf.gz"

# Read Group(s)
# READ_GROUPS="\"@RG\tID:1\tSM:${YSEQID}\tPL:ILLUMINA\tPU:UNIT1\tLB:${YSEQID}\""

if [[ $1 != "-resume" ]]
then


	# Pipeline commands:
	${BWA} ${BWA_PARAM} $REF $READS | \
	samtools view -@ $NUM_THREADS -b -t $REF -o $BAMFILE -
	samtools sort -@ $NUM_THREADS -T ${SORTDIR}sorted -o $BAMFILE_SORTED $BAMFILE
	samtools index -@ $NUM_THREADS $BAMFILE_SORTED
	samtools idxstats $BAMFILE_SORTED > ${BAMFILE_SORTED}.idxstats.tsv
fi

# Update summary file with BAM filesize:
BAM_FILESIZE=`du -kh "${BAMFILE_SORTED}" | cut -f1`
echo "Size of ${BAMFILE_SORTED} = ${BAM_FILESIZE}"
cat ${YSEQID}_result_summary.txt | sed "s/<bamsize>/${BAM_FILESIZE}byte/g" >${YSEQID}_result_summary.txt.tmp
rm -f ${YSEQID}_result_summary.txt
mv ${YSEQID}_result_summary.txt.tmp ${YSEQID}_result_summary.txt


# We're reading the first 100k sequences from the first fastq file to determine the average read length
READLENGTH=$(samtools bam2fq -@ ${NUM_THREADS} ${BAMFILE_SORTED} | head -n 400000 | awk '{if(NR%4==2) {count++; bases += length} } END{print(bases/count)}')
echo "The average read length is ${READLENGTH}"

# Calculating the basecounts from the idxstats output
MAPPED_READS_HG38=$(perl addRows.pl "--column=2" "--filename=${BAMFILE_SORTED}.idxstats.tsv")
#MAPPED_READS_HG38=`perl addRows.pl --column=2 --filename=${BAMFILE_SORTED}.idxstats.tsv`
MAPPED_BASES_HG38=$(perl multiply.pl "--multiplicands=$MAPPED_READS_HG38,$READLENGTH" "-round")

#MAPPED_BASES_HG38=$(awk -v rl="${READLENGTH}" '{m+=$3}END{printf("%i",m*rl)}' < ${BAMFILE_SORTED}.idxstats.tsv)
echo "The hg38 mapped basecount is ${MAPPED_BASES_HG38}" >>basecounts.txt
echo "The hg38 mapped basecount is ${MAPPED_BASES_HG38}"

UNMAPPED_READS=$(perl addRows.pl "--column=3" "--filename=${BAMFILE_SORTED}.idxstats.tsv")
#UNMAPPED_READS=`perl addRows.pl --column=3 --filename=${BAMFILE_SORTED}.idxstats.tsv`
UNMAPPED_BASES=$(perl multiply.pl "--multiplicands=$UNMAPPED_READS,$READLENGTH" "-round")


#UNMAPPED_BASES=$(awk -v rl="${READLENGTH}" '{u+=$4}END{printf("%i",u*rl)}' < ${BAMFILE_SORTED}.idxstats.tsv)
echo "The unmapped basecount is ${UNMAPPED_BASES}" >>basecounts.txt
echo "The unmapped basecount is ${UNMAPPED_BASES}"

SEQUENCED_BASES="$( printf '%s + %s\n' "$MAPPED_BASES_HG38" "$UNMAPPED_BASES" | bc )"
echo "Total sequenced bases: ${SEQUENCED_BASES}" >>basecounts.txt

echo "Total sequenced bases: ${SEQUENCED_BASES}"

# Calculating the coverage from the idxstats output
COVERAGE_HG38=$(awk -v rl="${READLENGTH}" '{x+=$2;m+=$3}END{print m*rl/x "x"}' < ${YSEQID}_bwa-mem_hg38_sorted.bam.idxstats.tsv)
echo "The hg38 coverage is ${COVERAGE_HG38}"

if [[ $1 != "-resume" ]]
then
	# Delete no longer needed large files
	rm -f $BAMFILE # keep $BAMFILE_SORTED
fi

# Easy tview file
echo "#!/bin/bash" > tview_${YSEQID}.sh
echo "samtools tview ${BAMFILE_SORTED} ${REF}" >> tview_${YSEQID}.sh
chmod a+x tview_${YSEQID}.sh



if [[ $1 != "-resume" ]]
then

	# Separate chrY & mtDNA BAM files
	samtools view -@ $NUM_THREADS -b -o "${YSEQID}_bwa-mem_hg38_chrY.bam" $BAMFILE_SORTED chrY &
	samtools view -@ $NUM_THREADS -b -o "${YSEQID}_bwa-mem_rCRS_chrM.bam" $BAMFILE_SORTED chrM &
	wait
	samtools index "${YSEQID}_bwa-mem_hg38_chrY.bam" &
	samtools index "${YSEQID}_bwa-mem_rCRS_chrM.bam" &
	wait


	#rsync -e "ssh" --progress -v ${YSEQID}_bwa-mem_rCRS_chrM.bam ${T4VPS}:/genomes/${YSEQID}/
	#rsync -e "ssh" --progress -v ${YSEQID}_bwa-mem_rCRS_chrM.bam.bai ${T4VPS}:/genomes/${YSEQID}/

	#rsync -e "ssh" --progress -v ${YSEQID}_bwa-mem_hg*_chrY.bam ${T4VPS}:/genomes/${YSEQID}/
	#rsync -e "ssh" --progress -v ${YSEQID}_bwa-mem_hg*_chrY.bam.bai ${T4VPS}:/genomes/${YSEQID}/

	# mtDNA allele calling $ FASTA file generation
	bcftools mpileup -r chrM -Ou -C 50 -f ${REF} ${BAMFILE_SORTED} | bcftools call --threads $NUM_THREADS -O z -v -m -P 0  > chrM_${VCF_FILE}.gz
	tabix chrM_${VCF_FILE}.gz
	samtools faidx $REF chrM | bcftools consensus chrM_${VCF_FILE}.gz -o ${YSEQID}_mtDNA.fasta


	# Annotate all known Y chromosome SNPs from Ybrowse
	wget http://ybrowse.org/gbrowse2/gff/snps_hg38.vcf.gz
	wget http://ybrowse.org/gbrowse2/gff/snps_hg38.vcf.gz.tbi

	bcftools mpileup -r chrY -C 0 --threads $NUM_THREADS -O z -f $REF $BAMFILE_SORTED > chrY_raw_${VCF_FILE}.gz
	tabix chrY_raw_${VCF_FILE}.gz

	bcftools merge -m both -O z snps_hg38.vcf.gz chrY_raw_${VCF_FILE}.gz > chrY_merged_${VCF_FILE}.gz
	tabix chrY_merged_${VCF_FILE}.gz
	echo "merging completed"

	bcftools call -O z -m -P 0 chrY_merged_${VCF_FILE}.gz > chrY_called_${VCF_FILE}.gz
	tabix chrY_called_${VCF_FILE}.gz
	echo "calling completed"

	bcftools view -O z -k -s ^HARRY,ALIEN chrY_called_${VCF_FILE}.gz >chrY_cleaned_${VCF_FILE}.gz
	tabix chrY_cleaned_${VCF_FILE}.gz
	echo "cleanup completed"

	bcftools filter -O z -i '(GT=="1/1" && AA==REF) || (GT=="0/0" && AA==ALT)' chrY_cleaned_${VCF_FILE}.gz >chrY_derived_${VCF_FILE}.gz &
	bcftools filter -O z -i '(GT=="0/0" && AA==REF) || (GT=="1/1" && AA==ALT)' chrY_cleaned_${VCF_FILE}.gz >chrY_ancestral_${VCF_FILE}.gz &
	wait
	tabix chrY_derived_${VCF_FILE}.gz &
	tabix chrY_ancestral_${VCF_FILE}.gz &
	wait
	echo "ancestral & derived filters completed"

	# Find the Y haplogroup
	bcftools query -f '%ID,' chrY_derived_${VCF_FILE}.gz | sed ':a;N;$!ba;s/\n//g' > ${YSEQID}_positives.txt &
	bcftools query -f '%ID,' chrY_ancestral_${VCF_FILE}.gz | sed ':a;N;$!ba;s/\n//g' > ${YSEQID}_negatives.txt &
	wait
fi


${PYTHON} ${CLADEFINDER} ${YFULLTREE} ${YSEQID}_positives.txt ${YSEQID}_negatives.txt ${YSEQID}_cladeFinderOutput.csv
echo "Y haplogroup checked"
echo "The 5 best hits are:"
head -n 5 ${YSEQID}_cladeFinderOutput.csv


YFULLHG="unknown"
YFULLPATH="unknown"
OLDIFS=$IFS
IFS=$'\t'
LINECOUNTER=$((0))

[ ! -f ${YSEQID}_cladeFinderOutput.csv ] && { echo "${YSEQID}_cladeFinderOutput.csv file not found"; }
while read -r haplogroup path purity depth score conflicting_negatives
do
	LINECOUNTER=$((LINECOUNTER + 1))
	if [ $LINECOUNTER -eq $((2)) ]
	then
		YFULLHG=$haplogroup
		YFULLPATH=$path
	fi		

done < ${YSEQID}_cladeFinderOutput.csv
IFS=$OLDIFS


echo "YFULLHG: ${YFULLHG}"
echo "YFULLPATH: ${YFULLPATH}"

cp /genomes/0/script-templates/create-mtDNA.py .

cp /genomes/0/script-templates/getEquivalentAndDownstreamSNPs.py .
PHYLOEQ_SNPS_FILE=${YSEQID}_PHYLOEQ_SNPS.tsv
DOWNSTR_SNPS_FILE=${YSEQID}_DOWNSTR_SNPS.tsv
${PYTHON} getEquivalentAndDownstreamSNPs.py ${YFULLTREE} "${YFULLHG}" chrY_cleaned_${VCF_FILE}.gz ${YSEQID}_positives.txt ${YSEQID}_negatives.txt ${PHYLOEQ_SNPS_FILE} ${DOWNSTR_SNPS_FILE} ${YSEQID}_bwa-mem_hg38_sorted.bam.idxstats.tsv

cp /genomes/0/script-templates/getMTDNADifferences.py .
MTDNA_SNPS_FILE=${YSEQID}_MTDNA_SNPS.tsv
HAPLOGREP_JAR=haplogrep-2.1.25.jar
cp /genomes/0/haplogrep/${HAPLOGREP_JAR} .

MTHAPLOGROUP=$(${PYTHON} getMTDNADifferences.py -process chrM_${VCF_FILE}.gz ${HAPLOGREP_JAR} ${MTDNA_SNPS_FILE})


# Update summary file with Y haplogroup:

# Note: Since the YFull tree uses forward slashes for identical SNPs with different names, we replace the sed separator "/" with "|"
# This keeps us from escaping every forward slash. But don't confuse it with the pipe character in the same command!

bcftools filter -i 'TYPE="snp" && QUAL>=30 && GT=="1/1" && DP4[2] + DP4[3] >=2 && (DP4[0] + DP4[1]) / DP < 0.1' chrY_called_${VCF_FILE}.gz | bcftools view -n -O z -s ^HARRY,ALIEN >chrY_novel_SNPs_${VCF_FILE}.gz
tabix chrY_novel_SNPs_${VCF_FILE}.gz
zcat chrY_novel_SNPs_${VCF_FILE}.gz | grep -v "##" >chrY_novel_SNPs_${VCF_FILE}.tsv
echo "novel SNP calling completed"

# Create Novel SNP template after checking for unreliable ranges and identity within Y Chromosome
cp /genomes/0/script-templates/identityResolutionTemplateCreator.py .
${PYTHON} identityResolutionTemplateCreator.py -batch chrY_novel_SNPs_${VCF_FILE}.tsv chrY_novel_SNPs_${YSEQID}_hg38.tsv ${REF}

bcftools filter -i 'TYPE="indel" && QUAL>=30 && GT=="1/1" && DP4[2] + DP4[3] >=2 && (DP4[0] + DP4[1]) / DP < 0.1' chrY_called_${VCF_FILE}.gz | bcftools view -n -O z -s ^HARRY,ALIEN >chrY_INDELs_${VCF_FILE}.gz
tabix chrY_INDELs_${VCF_FILE}.gz
echo "INDEL calling completed"



if [[ $1 != "-noautosomal" ]]
then

	# Parallel SNP calling by chromosome

	PARAMC=0

	bcftools mpileup -r chr1 -Ou -C ${PARAMC} -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0 > chr1_${VCF_FILE}.gz &
	bcftools mpileup -r chr2 -Ou -C ${PARAMC} -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0 > chr2_${VCF_FILE}.gz &
	bcftools mpileup -r chr3 -Ou -C ${PARAMC} -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0 > chr3_${VCF_FILE}.gz &
	bcftools mpileup -r chr4 -Ou -C ${PARAMC} -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0 > chr4_${VCF_FILE}.gz &
	bcftools mpileup -r chr5 -Ou -C ${PARAMC} -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0 > chr5_${VCF_FILE}.gz &
	bcftools mpileup -r chr6 -Ou -C ${PARAMC} -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0 > chr6_${VCF_FILE}.gz &
	bcftools mpileup -r chr7 -Ou -C ${PARAMC} -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0 > chr7_${VCF_FILE}.gz &
	bcftools mpileup -r chr8 -Ou -C ${PARAMC} -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0 > chr8_${VCF_FILE}.gz &
	bcftools mpileup -r chr9 -Ou -C ${PARAMC} -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0 > chr9_${VCF_FILE}.gz &
	bcftools mpileup -r chr10 -Ou -C ${PARAMC} -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0 > chr10_${VCF_FILE}.gz &
	bcftools mpileup -r chr11 -Ou -C ${PARAMC} -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0 > chr11_${VCF_FILE}.gz &
	bcftools mpileup -r chr12 -Ou -C ${PARAMC} -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0 > chr12_${VCF_FILE}.gz &
	bcftools mpileup -r chr13 -Ou -C ${PARAMC} -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0 > chr13_${VCF_FILE}.gz &
	bcftools mpileup -r chr14 -Ou -C ${PARAMC} -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0 > chr14_${VCF_FILE}.gz &
	bcftools mpileup -r chr15 -Ou -C ${PARAMC} -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0 > chr15_${VCF_FILE}.gz &
	bcftools mpileup -r chr16 -Ou -C ${PARAMC} -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0 > chr16_${VCF_FILE}.gz &
	bcftools mpileup -r chr17 -Ou -C ${PARAMC} -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0 > chr17_${VCF_FILE}.gz &
	bcftools mpileup -r chr18 -Ou -C ${PARAMC} -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0 > chr18_${VCF_FILE}.gz &
	bcftools mpileup -r chr19 -Ou -C ${PARAMC} -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0 > chr19_${VCF_FILE}.gz &
	bcftools mpileup -r chr20 -Ou -C ${PARAMC} -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0 > chr20_${VCF_FILE}.gz &
	bcftools mpileup -r chr21 -Ou -C ${PARAMC} -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0 > chr21_${VCF_FILE}.gz &
	bcftools mpileup -r chr22 -Ou -C ${PARAMC} -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0 > chr22_${VCF_FILE}.gz &
	bcftools mpileup -r chrX -Ou -C ${PARAMC} -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0  > chrX_${VCF_FILE}.gz &
	bcftools mpileup -r chrY -Ou -C ${PARAMC} -f $REF $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -v -V indels -m -P 0  > chrY_${VCF_FILE}.gz &
	wait

	# Concatenate all chromosome VCFs to one big file
	bcftools concat -O z chr[1-9]_${VCF_FILE}.gz chr[1-2][0-9]_${VCF_FILE}.gz chr[M,X-Y]_${VCF_FILE}.gz > ${VCF_FILE}.gz
	tabix ${VCF_FILE}.gz

	# Delete no longer needed VCFs
	rm -f chr[0-9]_${VCF_FILE}.gz chr[1-2][0-9]_${VCF_FILE}.gz chrX_${VCF_FILE}.gz

	rm -f fastq.fastq
	rm -f bed.bed
	rm -f alignment.paf

	#samtools stats $BAMFILE_SORTED > ${BAMFILE_SORTED}.stats &

	chgrp -R users ./*
	chmod -R g+rw ./*





fi

cat ${YSEQID}_result_summary.txt

END=$(date +%s.%N)

DIFF=$(echo "$END - $START" | bc)

echo ${DIFF}




