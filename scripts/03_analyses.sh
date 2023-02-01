# 1) remove batch effect, identify outliers and generate clean input files
Rscript 04_filter_SNPs.R stacks_m3M7_n0/singletonsHD+VS_D2_p5r0.5

# 2) RADpainter on clean data and outliers
sbatch RADpainter.sh stacks_m3M7_n0/singletonsHD+VS_D2_p5r0.5/populations.haps.nonoutl.radpainter
sbatch RADpainter.sh stacks_m3M7_n0/singletonsHD+VS_D2_p5r0.5/populations.haps.outliers.radpainter

sbatch RADpainter.sh stacks_m3M7_n0/singletonsHD+VS_D2_POOL7_p3r0.5/populations.haps.nonoutl.radpainter
sbatch RADpainter.sh stacks_m3M7_n0/singletonsHD+VS_D2_POOL7_p3r0.5/populations.haps.outliers.radpainter

sbatch RADpainter.sh stacks_m3M7_n0/singletonsHD+VS_D2_noPOOL7_p4r0.5/populations.haps.nonoutl.radpainter
sbatch RADpainter.sh stacks_m3M7_n0/singletonsHD+VS_D2_noPOOL7_p4r0.5/populations.haps.outliers.radpainter


# 3) fastSTRUCTURE on clean data and outliers
for indir in stacks_m3M7_n0/singletonsHD+VS_D2_p5r0.5 stacks_m3M7_n0/singletonsHD+VS_D2_*POOL*r0.5
do
	srun bash fastSTRUCTURE.sh $indir/populations.plink.nonoutl
	srun bash fastSTRUCTURE.sh $indir/populations.plink.outliers
done


#
# For Figure 1, we need:
# RADPainter run on full data
#
srun 03_RADpainter.sh stacks_m3M7_n0/singletonsHD+VS_D2_p5r0.5/populations.haps.clean.radpainter


# 1) run PCA and DAPC to identify important SNPs for population structure
# 2) run fastSTRUCTURE
# 3) make multipanel figure

# 1) run BAYESCAN

for indir in stacks_m3M*_n*/singletons*0.8
do
	outdir=$indir/BayeScan
	mkdir -p $outdir
	# convert VCF to BAYESCAN format
	sbatch -c 8 --mem 32G --wrap "java -Xmx20G -jar PGDSpider_2.1.1.5/PGDSpider2-cli.jar -spid vcf2bayescan.spid \
		-inputfile $indir/populations.snps.vcf -inputformat VCF \
		-outputfile $outdir/input.snps.txt -outputformat GESTE_BAYE_SCAN
	~/sharedscratch/apps/BayeScan2.1/binaries/BayeScan2.1_linux64bits -od $outdir -o output -threads 8 $outdir/input.snps.txt"
done

# 1) PCA/DAPC in R
module load adegenet-2.1.3
for a in stacks_m3M*_n*/singletons*0.8
do
	sbatch --mem 8G --wrap "Rscript 04_DAPC_new.R $a"
done

# 2) fastSTRUCTURE
module load plink-1.90b6.10
module load python-2.7.18-gcc-10.1.0
module load gsl-2.5
export PYTHONPATH=$PYTHONPATH:/uoa/home/s02mw5/sharedscratch/apps/fastStructure/vars

for indir in stacks_m3M4_n0/singletons*
do
	plink --file $indir/populations.plink --out $indir/populations.plink --allow-extra-chr
	plink --file $indir/populations.plink.snpzip --out $indir/populations.plink.snpzip --allow-extra-chr
	mkdir -p $indir/fastStructure
	sbatch -a 2-4 --mem 8G --wrap "python ~/sharedscratch/apps/fastStructure/structure.py -K \$SLURM_ARRAY_TASK_ID --input=$indir/populations.plink.snpzip --output=$indir/fastStructure/snpzip_logistic --prior=logistic"
done
#sbatch -a 2-4 --mem 64G --wrap "python ~/sharedscratch/apps/fastStructure/structure.py -K \$SLURM_ARRAY_TASK_ID --input=$indir/populations.plink --output=$indir/fastStructure/allsnps_logistic --prior=logistic"
#sbatch -a 2-4 --mem 8G --wrap "python ~/sharedscratch/apps/fastStructure/structure.py -K \$SLURM_ARRAY_TASK_ID --input=$indir/populations.plink.snpzip --output=$indir/fastStructure/snpzip_simple --prior=simple"

### make plots
for indir in stacks_m3M[47]_n0/singletons*
do
	Rscript 04_structureplots.R $indir
done



# 3) RADpainter
module load gsl-2.5
module load r-3.6.1
indir=sweep_m3M8_n0/singletonsHD+VS_p5r0.5
mkdir $indir/RADpainter
for M in M3 M4 M7 M8 M9; do echo "$M($(grep $M popmaps/final_group_Mariana.txt | cut -f 1 | paste -s -d ','))"; done > $indir/RADpainter/pops.txt
srun ~/sharedscratch/apps/fineRADstructure/RADpainter paint $indir/populations.haps.radpainter.snpzip
mv $indir/populations.haps.radpainter_chunks.out $indir/RADpainter/snpzip_chunks.out
srun ~/sharedscratch/apps/fineRADstructure/finestructure -I 5 -K -M 4 -X -Y -x 100000 -y 100000 -z 1000 $indir/RADpainter/snpzip_chunks.out $indir/RADpainter/snpzip_chunks.out.mcmc.xml
srun ~/sharedscratch/apps/fineRADstructure/finestructure -m T -x 100000 -X -Y $indir/RADpainter/snpzip_chunks.out $indir/RADpainter/snpzip_chunks.out.mcmc.xml $indir/RADpainter/snpzip_chunks.out.tree.xml
grep "Tree" $indir/RADpainter/snpzip_chunks.out.tree.xml | sed 's/<\/*Tree>//g' > $indir/RADpainter/snpzip_chunks.out.tree
# visualise
srun Rscript ~/sharedscratch/apps/fineRADstructure/fineRADstructurePlot.R $indir/RADpainter/snpzip_chunks.out

mat <- read.table("sweep_m3M8_n0/singletonsHD+VS_p5r0.5/RADpainter/haps_chunks.out")
mat <- as.dist(mat[,-1])
cs <- cmdscale(mat)

