#!/bin/bash
#SBATCH --mem 16G
#
# carry out RADPainter and fineRADstructure analysis
#
module load gsl-2.5
module load r-3.6.1

infile=$1
outdir=$(dirname $infile)/RADpainter
mkdir $outdir

# extract number of pops from directory name
npops=$(echo "$infile" | grep -o "_p[0-9]" | grep -o [0-9])

# make population file
for M in M3 M4 M7 M8 M9
do
	#echo "$M($(grep $M popmaps/final_group_Mariana.txt | cut -f 1 | paste -s -d ','))"
	echo "($(grep $M popmaps/final_group_Mariana.txt | cut -f 1 | grep -Fw -f <(head -n 1 $infile | tr '\t' '\n') | paste -s -d ','))"
done | grep -v "()" | tr -d '\n' | awk '{print "<Pop>"$0"</Pop>"}' > $outdir/pops.txt

# infer coancestry matrix
~/sharedscratch/apps/fineRADstructure/RADpainter paint $infile

infile=${infile/.radpainter/}
outfile=$(basename $infile)_chunks.out
mv $infile*out $outdir

# 2) infer populations and tree
~/sharedscratch/apps/fineRADstructure/finestructure \
	-X -Y -x 1000000 -y 1000000 -z 1000 \
	$outdir/$outfile \
	$outdir/$outfile.mcmc.xml
~/sharedscratch/apps/fineRADstructure/finestructure \
	-e max -X -Y \
	$outdir/$outfile \
	$outdir/$outfile.mcmc.xml \
	$outdir/$outfile.mcmc.max
~/sharedscratch/apps/fineRADstructure/finestructure \
	-m T -X -Y -x 100000 \
	$outdir/$outfile \
	$outdir/$outfile.mcmc.xml \
	$outdir/$outfile.tree.xml
~/sharedscratch/apps/fineRADstructure/finestructure \
	-e TREE -X -Y \
	$outdir/$outfile \
	$outdir/$outfile.mcmc.xml \
	$outdir/$outfile.tree.xml \
	$outdir/$outfile.tree.nw
# now for K=npops
# 2) infer populations and tree
~/sharedscratch/apps/fineRADstructure/finestructure \
	-I $npops -K -X -Y -x 1000000 -y 1000000 -z 1000 \
	$outdir/$outfile \
	$outdir/$outfile.K.mcmc.xml
~/sharedscratch/apps/fineRADstructure/finestructure \
	-e max -X -Y \
	$outdir/$outfile \
	$outdir/$outfile.K.mcmc.xml \
	$outdir/$outfile.K.mcmc.max
~/sharedscratch/apps/fineRADstructure/finestructure \
	-m T -X -Y -x 100000 \
	$outdir/$outfile \
	$outdir/$outfile.K.mcmc.xml \
	$outdir/$outfile.K.tree.xml
~/sharedscratch/apps/fineRADstructure/finestructure \
	-e TREE -X -Y \
	$outdir/$outfile \
	$outdir/$outfile.K.mcmc.xml \
	$outdir/$outfile.K.tree.xml \
	$outdir/$outfile.K.tree.nw	

# 3) infer admixture from best populations or from pre-defined populations
~/sharedscratch/apps/fineRADstructure/finestructure \
	-m admixture -X -Y -x 1000000 -y 1000000 -z 1000 \
	$outdir/$outfile \
	$outdir/$outfile.mcmc.max \
	$outdir/$outfile.admix
~/sharedscratch/apps/fineRADstructure/finestructure \
	-m admixture -X -Y -x 1000000 -y 1000000 -z 1000 \
	$outdir/$outfile \
	$outdir/pops.txt \
	$outdir/$outfile.admix.K
~/sharedscratch/apps/fineRADstructure/finestructure \
	-e admixture -X -Y -o $outdir/$outfile.admix \
	$outdir/$outfile \
	$outdir/$outfile.admix \
	$outdir/$outfile.admix.Q
~/sharedscratch/apps/fineRADstructure/finestructure \
	-e admixture -X -Y -o $outdir/$outfile.admix.K \
	$outdir/$outfile \
	$outdir/$outfile.admix.K \
	$outdir/$outfile.admix.K.Q

# extract Q matrix from best iteration
for ad in $outdir/$outfile.admix $outdir/$outfile.admix.K
do
	grep "Likelihood" $ad | sort | tail -n 1 | grep -f - -m 1 -A 53 $ad | tail -n +5 > $ad.Qbest
done

# visualise
Rscript ~/sharedscratch/apps/fineRADstructure/fineRADstructurePlot.R $outdir/$outfile
