#
# lanes have been concatenated for each pool

# go straight into process_radtags
# do not use PEAR. TSV2BAM can handle overlap in ddRAD PE reads directly

# BUT: process_radtags cannot detect R2 reads with readthrough into barcodes
# need to pre-process with CUTADAPT
# inspection of reads reveals the following adapters (including restriction sites and barcodes)
# R1: AATTAGATCGGAAGAGCACACGTCTGAACTCCAGTCA
# R2: CTGCANNNNNNNNAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

module load seqtk-1.3
module load cutadapt-2.3
# make FASTA of R2 adapters
cut -f 1 raw/*.barcode | sort | uniq | awk '{print ">revbarcode"NR"\n"$0}' | seqtk seq -r | awk '!/^>/{print "CTGCA"$0"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"}/^>/{print}' > raw/cutadapt.revadapters.fa

# discard read pairs with evidence of adapters
for pool in pool0{4..9}
do
	mkdir -p clean_cutadapt/$pool
	sbatch -c 8 --mem 16G --wrap "cutadapt -a AATTAGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
				-A file:raw/cutadapt.revadapters.fa \
				--discard -O 4 -j 8 \
				-o clean_cutadapt/$pool/$pool.R1.fq.gz \
				-p clean_cutadapt/$pool/$pool.R2.fq.gz \
				raw/$pool/concatenated_1.gz \
				raw/$pool/concatenated_2.gz"
done

stacksdir=~/sharedscratch/apps/stacks-2.58

mkdir -p clean_stacks2
for pool in pool0{4..9}
do
	sbatch --mem 16G --wrap "$stacksdir/process_radtags \
		-1 clean_cutadapt/$pool/$pool.R1.fq.gz \
		-2 clean_cutadapt/$pool/$pool.R2.fq.gz \
		-b raw/*$pool.barcode --inline-null --barcode-dist-1 1 \
		-o clean_stacks2 \
		-c -q -r -y gzfastq \
		--renz_1 pstI --renz_2 ecoRI \
		--adapter_1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
		--adapter_2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
		--adapter_mm 2 \
		--filter-illumina"
done

# check for degree of read overlap using PEAR or FLASH
# this can help identify libraries with different size selection
for pool in pool0{4..9}
do
	sbatch -c 8 --mem 8G --wrap "~/sharedscratch/apps/FLASH-1.2.11-Linux-x86_64/flash -t 8 -c \
		clean_cutadapt/$pool/$pool.R1.fq.gz \
		clean_cutadapt/$pool/$pool.R2.fq.gz \
		> /dev/null 2> clean_cutadapt/$pool/FLASH.log"
done
