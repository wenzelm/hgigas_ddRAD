stacksdir=~/sharedscratch/apps/stacks-2.58

# 1) build loci in USTACKS from R1 files
# do a sweep from M0 to M15
# to aid blacklisting of potential paralogs:
# - no --deleverage
# - maximum of two alleles per locus
# - no gappy alignments

for M in {0..15}
do
	outdir=stacks_m3M${M}
	mkdir -p $outdir
	sbatch -c 4 --mem 16G -a 1-$(wc -l < popmaps/final_group_Mariana.txt) --wrap "$stacksdir/ustacks \
		-f clean_stacks/\$(cut -f 1 popmaps/final_group_Mariana.txt | head -n \$SLURM_ARRAY_TASK_ID | tail -n 1).1.fq.gz \
		-i \$SLURM_ARRAY_TASK_ID -o $outdir \
		-m 3 -M $M -H --disable-gapped --max-gaps 0 --min-aln-len 1.0 --max-locus-stacks 2 \
		--model-type bounded --bound-high 0.05 -p 4"
done

# remove the ".1." suffix from output files to avoid errors in cstacks
rename ".1." "." stacks_m3M*/*.tsv.gz

# to decide which M is best, look at numbers of 1-allelic, 2-allelic, and 3+-allelic loci
# 2-allelic loci should peak early on
# 1-allelic loci should reach a valley early on

for a in stacks_m3M!(*_*)/[0-9]*.tags.tsv.gz
do
	echo "$a"
	srun zless $a | awk '$3=="consensus" && $5=="0" && $6=="0" && $7=="0"{print A; L=$2} 
			$3=="primary"{A=1+$4; if(A>2) A="3+"}END{print A}' \
			| sort | uniq -c | sort -g -k 2 | grep -v "^\s*1\s" > $a.locusploidy.txt
done
# summarise
function ploidysum {
	grep -h "3+$" stacks_m3M$1/*ploidy* | awk '{print $1}'
}
paste <(ploidysum 0) <(ploidysum 1) <(ploidysum 2) <(ploidysum 3) <(ploidysum 4) <(ploidysum 5) <(ploidysum 6) <(ploidysum 7) <(ploidysum 8) <(ploidysum 9) <(ploidysum 10) <(ploidysum 11) <(ploidysum 12) <(ploidysum 13) <(ploidysum 14) <(ploidysum 15)

# total loci
for M in {0..15}; do echo "$M"; zgrep -B1 "completed" stacks_m3M$M/*tags.tsv.gz | grep -v "#" | cut -f 2 > stacks_m3M$M/ustacks.totalloci.txt; done

# also count how many loci were blacklisted
# .tags flags are: $5 = deleveraged (always 0 because switched off); $6 = blacklisted by deleverage; $7 = lumberjack stack
for a in stacks_m3M!(*_*)
do
	for b in $a/[^c]*.tags.tsv.gz
	do
		echo "$b"
		paste <(zless $b | awk '$3=="consensus" && $6=="1"' | wc -l) <(zless $b | awk '$3=="consensus" && $7=="1"' | wc -l) >> $a/ustacksflags.txt
	done
done

# 2) build locus catalog
# - usual way is n=M
# - for paralog detection it's better to do n=0 or n=1
for M in {0..15}
do
	outdir=stacks_m3M${M}
	mkdir -p $outdir
	rename ".1." "." $outdir/*
	sbatch -c 16 --mem 96G --wrap "$stacksdir/cstacks \
		-P $outdir \
		-M popmaps/final_group_Mariana.txt \
		-n $M -p 16 --disable-gapped --max-gaps 0 --min-aln-len 1.0"
done
# n=0
for M in {0..15}
do
	outdir=stacks_m3M${M}_n0
	mkdir -p $outdir
	cd $outdir
	ln -s ../sweep_m3M${M}/[0-9]*.[sta]*.tsv.gz .
	cd -
	sbatch -c 32 --mem 96G --wrap "$stacksdir/cstacks \
		-P $outdir \
		-M popmaps/final_group_Mariana.txt \
		-n 0 -p 32 --disable-gapped --max-gaps 0 --min-aln-len 1.0"
done
# n=1
for M in {0..15}
do
	outdir=stacks_m3M${M}_n1
	mkdir -p $outdir
	cd $outdir
	ln -s ../sweep_m3M${M}/[0-9]*.[sta]*.tsv.gz .
	cd -
	sbatch -c 32 --mem 96G --wrap "$stacksdir/cstacks \
		-P $outdir \
		-M popmaps/final_group_Mariana.txt \
		-n 1 -p 32 --disable-gapped --max-gaps 0 --min-aln-len 1.0"
done

# 3) Paralog check (PMERGE-style)
# - cluster catalogue loci at varying %
# - anything that clusters is potential paralogs
for ctl in stacks_m3M*_n*/catalog.tags.tsv.gz
do
	zless $ctl | awk '$3=="consensus"{print ">"$2"\n"$6}' > $ctl.fa
done
for ctl in stacks_m3M*_n*/catalog.tags.tsv.gz.fa
do
	for id in {1..9}
	do
		sbatch -c 8 --mem 16G --wrap "./vsearch --cluster_fast $ctl --id 0.$id --iddef 1 --qmask none --leftjust --maxgaps 0 --maxrejects 0 --threads 8 --uc $ctl.vsearch${id}0"
	done
done
# print summaries on screen
for a in stacks_m3M*/catalog.tags.tsv.*fa*.vsearch*0; do echo "$a: $(awk '$1=="C"' $a | wc -l)"; done
# write lists of paralogs
for a in stacks_m3M*/catalog.tags.tsv.*fa*.vsearch*0; do awk '$1=="H"{print $9"\n"$10}' $a | sort | uniq > $a.paralogs; done

# 4) complete STACKS genotyping pipeline (sstacks, tsv2bam, gstacks)
# sstacks
for outdir in stacks_m3M{4,5,6,7,12}_n*
do
	sbatch -c 32 --mem 96G --wrap "$stacksdir/sstacks \
		-P $outdir \
		-M popmaps/final_group_Mariana.txt \
		-p 32 --disable-gapped"
done
# tsv2bam
for outdir in stacks_m3M{4,5,6,7,12}_n*
do
	sbatch -c 32 --mem 96G --wrap "$stacksdir/tsv2bam \
		-P $outdir \
		-M popmaps/final_group_Mariana.txt \
		-R clean_stacks \
		-t 32"
done
# gstacks
for outdir in stacks_m3M{4,5,6,7,12}_n*
do
	sbatch -c 32 --mem 96G --wrap "$stacksdir/gstacks \
		-P $outdir \
		-M popmaps/final_group_Mariana.txt \
		-t 32 \
		--model marukilow"
done

# 5) identify paralogues within individuals from catalogue matches
# this is equivalent to "haplotyping": if multiple loci match to the same catalogue locus, it's a paralog
for ctl in stacks_m3M*_n*/catalog.fa.gz
do
	sbatch --wrap "zless $(dirname $ctl)/*.matches.tsv.gz | cut -f 1,2,3 | sort | uniq \
		| cut -f 1,2 | sort | uniq -c \
		| awk '\$1>1{print \$2}' | sort | uniq > $ctl.matches.paralogs"
done

# 6) identify contaminant sequences among the catalogue
# use BLAST and KRAKEN2
# PlusPF-16 database from https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_16gb_20210517.tar.gz
# https://benlangmead.github.io/aws-indexes/k2
module load blast-2.9.0
export BLASTDB=~/sharedscratch/scops_RAD_stacks2.58/blast
module load kraken2-2.0.8_beta
DBNAME=~/sharedscratch/scops_RAD_stacks2.58/kraken_db

for ctl in stacks_m3M*_n*/catalog.fa.gz
do
	# RefSeq BLAST
	sbatch -c 12 --mem 24G --wrap "zcat $ctl | blastn -db blastalias -query - \
					-evalue 1e-3 \
					-outfmt 6 -num_threads 12 -out $ctl.contaminants_blastn"
	# KRAKEN2
	sbatch -c 12 --mem 96G --wrap "kraken2 --db $DBNAME --use-names \
		--threads 12 --output $ctl.contaminants_kraken2 --gzip-compressed $ctl"
done

# make contaminant blacklist
for ctl in stacks_m3M*_n*/catalog.fa.gz
do
	cat <(awk '$1=="C"' $ctl.contaminants_kraken2 | cut -f 2) \
	<(cut -f 1 $ctl.contaminants_blastn) | sort | uniq > $ctl.contaminants_blacklist
done


# 7) final paralog check with paralog-finder (HDPlot)
# extract VCF for paralog-finder (HDPlot)
# ensure that SNPs are filtered according to final intentions (maf, het), 
# otherwise the whitelist is way too strict
for dset in stacks_m3M{4,5,6,7,12}_n*
do
	for r in 0 
	do
		popout=$dset/paralogfinder_p5r${r}
		mkdir -p $popout
		sbatch -c 16 --mem 96G --wrap "$stacksdir/populations \
			-P $dset \
			-B $dset/catalog.fa.gz.contaminants_blacklist \
			-O $popout \
			-M popmaps/final_group_Mariana.txt \
			-p 5 -r $r -t 16 --min-maf 0.01 --max-obs-het 0.6 \
			--vcf
			rm $popout/*.tsv $popout/*haps*"
	done
done
# run paralog-finder
# z-score filter of +-2 is strict enough. +-1 is too strict.
conda activate paralog-finder
module load r-3.6.1
for vcf in stacks_m3M*_n*/paralogfinder_p5r0/populations.snps.vcf
do
	sbatch --wrap "~/sharedscratch/apps/paralog-finder/HDplot_process_vcf.py -i $vcf
		Rscript ~/sharedscratch/apps/paralog-finder/HDplot_graphs.R -i $vcf.HDPlot.depthsBias
		~/sharedscratch/apps/paralog-finder/blacklist_paralogs.py -i $vcf.HDPlot.depthsBias --maxH 0.6 --minD -2 --maxD 2
		rename 'Bias' 'BiasD2' $vcf.HDPlot.depthsBias_*list
		~/sharedscratch/apps/paralog-finder/blacklist_paralogs.py -i $vcf.HDPlot.depthsBias --maxH 0.6 --minD -1 --maxD 1
		rename 'Bias' 'BiasD1' $vcf.HDPlot.depthsBias_*list"
done

# 8) construct final singleton/paralog locus lists
# - singletons
#    blacklist contaminants + paralogs from all methods
# - paralogs
#    whitelist paralogs from all methods and remove contaminants
# - do the same for the intersection of methods instead of union
# 
for a in stacks_m3M{4,5,6,7,12}_n*
do
	for D in D1 D2
	do
	# no VSEARCH
	cat $a/*contaminants_blacklist \
		$a/*matches.paralogs \
		$a/paralogfinder_p5r0/*Bias${D}_paralogs.blacklist \
			| sort | uniq > $a/contams+paralogs100_${D}_blacklist.txt
	cat $a/*matches.paralogs \
		$a/paralogfinder_p5r0/*Bias${D}_paralogs.blacklist \
		| sort | uniq | grep -v -Fw -f $a/*contaminants_blacklist > $a/combined_paralogs100_${D}_whitelist.txt
	# VSEARCH 40% or 80%
	for vs in 40 80
	do
		cat $a/*contaminants_blacklist \
			$a/*vsearch$vs.paralogs \
			$a/*matches.paralogs \
			$a/paralogfinder_p5r0/*Bias${D}_paralogs.blacklist \
				| sort | uniq > $a/contams+paralogs${vs}_${D}_blacklist.txt
		cat $a/*vsearch$vs.paralogs \
			$a/*matches.paralogs \
			$a/paralogfinder_p5r0/*Bias${D}_paralogs.blacklist \
				| sort | uniq | grep -v -Fw -f $a/*contaminants_blacklist > $a/combined_paralogs${vs}_${D}_whitelist.txt
	done
	# (HDplot AND VSEARCH40) OR .matches
	cat $a/paralogfinder_p5r0/*Bias${D}_paralogs.blacklist $a/*vsearch40.paralogs | sort | uniq -c | awk '$1==2{print $2}' | cat - $a/*matches.paralogs $a/*contaminants_blacklist | sort | uniq > $a/contams+paralogsHD+VS_${D}_blacklist.txt
	cat $a/paralogfinder_p5r0/*Bias${D}_paralogs.blacklist $a/*vsearch40.paralogs | sort | uniq -c | awk '$1==2{print $2}' | cat - $a/*matches.paralogs | sort | uniq | grep -v -Fw -f $a/*contaminants_blacklist > $a/combined_paralogsHD+VS_${D}_whitelist.txt
	done
done

# 9) Extract final singleton and paralog datasets
for dset in stacks_m3M{4,5,6,7,12}_n*
do
	for vs in 40 80 100 HD+VS
	do
		for r in 0.5 1.0
		do
			for D in D1 D2
			do
	# singletons
	popout=$dset/singletons${vs}_${D}_p5r${r}
	mkdir -p $popout
	sbatch -c 16 --mem 64G --wrap "$stacksdir/populations \
		-P $dset \
		-O $popout \
		-B $dset/contams+paralogs${vs}_${D}_blacklist.txt \
		-M popmaps/final_group_Mariana.txt \
		-p 5 -r $r -t 16 --min-maf 0.01 --max-obs-het 0.6 \
		--fstats --hwe --fst-correction --vcf --plink --phylip-var-all --phylip-var --radpainter --fasta-samples"
	# paralogs
	popout=$dset/paralogs${vs}_${D}_p5r${r}
	mkdir -p $popout
	sbatch -c 16 --mem 64G --wrap "$stacksdir/populations \
		-P $dset \
		-O $popout \
		-W $dset/combined_paralogs${vs}_${D}_whitelist.txt \
		-M popmaps/final_group_Mariana.txt \
		-p 5 -r $r -t 16 --min-maf 0.01 --max-obs-het 0.6 \
		--fstats --hwe --fst-correction --vcf --plink --phylip-var-all --phylip-var --radpainter --fasta-samples"
			done
		done
	done
done
# extract data for pool7 only
popout=stacks_m3M7_n0/singletonsHD+VS_D2_POOL7_p3r0.8/
	mkdir -p $popout
	sbatch -c 16 --mem 64G --wrap "$stacksdir/populations \
		-P $dset \
		-O $popout \
		-B $dset/contams+paralogs${vs}_${D}_blacklist.txt \
		-M popmaps/final_group_Mariana.pool7.txt \
		-p 3 -r $r -t 16 --min-maf 0.01 --max-obs-het 0.6 \
		--fstats --hwe --fst-correction --vcf --plink --phylip-var-all --phylip-var --radpainter --fasta-samples --genepop"

# 10) annotate locus catalogue
# repeatmasker
module load repeatmasker-4.0.9_p2
RepeatMasker -pa 32 -species arthropoda -dir ~/localscratch stacks_m3M7_n0/catalog.fa.gz
mv ~/localscratch/catalog.fa.out stacks_m3M7_n0/catalog.fa.gz.RepeatMasker.out
mv ~/localscratch/catalog.fa.tbl stacks_m3M7_n0/catalog.fa.gz.RepeatMasker.tbl
awk '{print $5"\t"$11}' stacks_m3M7_n0/catalog.fa.gz.RepeatMasker.out | tail -n +4 | less > stacks_m3M7_n0/catalog.fa.gz.RepeatMasker.out.summary

# hirondellea transcriptome
module load blast-2.9.0
export BLASTDB=~/sharedscratch/scops_RAD_stacks2.58_new/blast
for ctl in stacks_m3M7_n0/catalog.fa.gz
do
	sbatch -c 12 --mem 24G --wrap "zcat $ctl | blastn -db GEZX01.1.fsa_nt -query - \
					-evalue 1e-3 \
					-outfmt 6 -num_threads 12 -out $ctl.GEZX01.1_blastn
		cut -f 1 stacks_m3M7_n0/catalog.fa.gz.GEZX01.1_blastn | uniq > stacks_m3M7_n0/catalog.fa.gz.GEZX01.1_blastn.loci"
done

# uniprot and GO
module load blast-2.9.0
# makedb
zless uniprot-taxonomy__Arthropoda6656_.fasta.gz | makeblastdb -in - -dbtype prot -out uniprot-taxonomy__Arthropoda6656 -title uniprot-taxonomy__Arthropoda6656

# do in parallel because blastx takes a long time
module load bbmap-38.51
partition.sh in=stacks_m3M7_n0/catalog.fa.gz out=uniprot_split/split%.fa ways=50
for ctl in stacks_m3M7_n0/catalog.fa.gz
do
        sbatch -a 0-49 -c 4 --mem 16G --wrap "blastx -db uniprot-taxonomy__Arthropoda6656 -query uniprot_split/split\${SLURM_ARRAY_TASK_ID}.fa \
                                        -evalue 1e-3 \
                                        -outfmt 6 -num_threads 4 -out uniprot_split/split\${SLURM_ARRAY_TASK_ID}.fa.blastx"
done


# organise GO table
zless uniprot-taxonomy__Arthropoda6656_GO.tab.gz | cut -f 1,2,4 | tr ';' '\n' | awk -F'\t' '/reviewed/{p=$1; $1=$3}{print p, $1}' OFS='\t' | grep "GO:" | gzip > uniprot-taxonomy__Arthropoda6656_GO.tsv.gz

