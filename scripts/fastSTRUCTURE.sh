module load plink-1.90b6.10
module load python-2.7.18-gcc-10.1.0
module load gsl-2.5
export PYTHONPATH=$PYTHONPATH:/uoa/home/s02mw5/sharedscratch/apps/fastStructure/vars

infile=$1
indir=$(dirname $infile)
mkdir -p $indir/fastStructure

plink --file $infile --out $infile --allow-extra-chr --make-bed
plink --file $infile --out $infile.ld --allow-extra-chr --indep-pairwise 1000000 1 0.9 --make-bed

sbatch -a 2-5 --mem 64G --wrap "python ~/sharedscratch/apps/fastStructure/structure.py \
		-K \$SLURM_ARRAY_TASK_ID --input=$infile --output=$indir/fastStructure/$(basename $infile).logistic --prior=logistic"
sbatch -a 2-5 --mem 64G --wrap "python ~/sharedscratch/apps/fastStructure/structure.py \
		-K \$SLURM_ARRAY_TASK_ID --input=$infile --output=$indir/fastStructure/$(basename $infile).simple --prior=simple"

sbatch -a 2-5 --mem 64G --wrap "python ~/sharedscratch/apps/fastStructure/structure.py \
		-K \$SLURM_ARRAY_TASK_ID --input=$infile.ld --output=$indir/fastStructure/$(basename $infile.ld).logistic --prior=logistic"
sbatch -a 2-5 --mem 64G --wrap "python ~/sharedscratch/apps/fastStructure/structure.py \
		-K \$SLURM_ARRAY_TASK_ID --input=$infile.ld --output=$indir/fastStructure/$(basename $infile.ld).simple --prior=simple"
