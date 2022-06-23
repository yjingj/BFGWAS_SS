#########################################################
############## Test on 1KG Example Data
############## Run BFGWAS_SS with 4 genome block
#########################################################

# Specify BFGWAS tool directory and working directory with writing access
BFGWAS_SS_dir="/home/jyang/GIT/BFGWAS_SS"
wkdir="/home/jyang/GIT/BFGWAS_SS/1KG_example/Test_Wkdir"
# file with all genome block file name heads
filehead=${BFGWAS_SS_dir}/1KG_example/ExData/fileheads_4region.txt
# genotype directory
geno_dir=${BFGWAS_SS_dir}/1KG_example/ExData/vcfs
# phenotype directory
pheno=${BFGWAS_SS_dir}/1KG_example/ExData/phenoAMD_1KG.txt

############ Step 1 ###############
####### First generate single variant GWAS Zscore statistic and LD correlation files
# LD coefficient directory
LDdir=${wkdir}/LDcorr
# Zscore statistics directory
Zscore_dir=${wkdir}/Zscore
LDwindow=1000000

cd $wkdir
${BFGWAS_SS_dir}/bin/GetRefLD.sh --wkdir ${wkdir} \
--toolE ${BFGWAS_SS_dir}/bin/Estep_mcmc \
--geno_dir ${geno_dir} --filehead ${filehead} --pheno ${pheno} \
--GTfield GT --genofile_type vcf --maf 0.001 --LDwindow ${LDwindow}  \
--Zscore_dir ${Zscore_dir} --LDdir ${LDdir} &

############ Step 2 ###############
########## Generate make file
Nsample=2540 # Sample size
maf=0.001 # MAF threshold
Nburnin=10000; # Burn-in iterations in MCMC
Nmcmc=10000; # MCMC iteration number
mkfile=${wkdir}/BFGWAS.mk

anno_dir=${BFGWAS_SS_dir}/1KG_example/ExData/annos
anno_code=${BFGWAS_SS_dir}/1KG_example/ExData/AnnoCode6.txt
hfile=${BFGWAS_SS_dir}/1KG_example/ExData/InitValues6.txt

em=3; # EM steps
cd $wkdir
${BFGWAS_SS_dir}/bin/gen_mkf.pl \
--wkdir ${wkdir} --tool_dir ${BFGWAS_SS_dir} \
--LDdir ${LDdir} --Zscore_dir ${Zscore_dir} \
--filehead ${filehead} --anno_dir ${anno_dir} --anno_code ${anno_code} \
--hfile ${hfile} --maf ${maf} --Nsample ${Nsample} \
--Nburnin ${Nburnin} --Nmcmc ${Nmcmc} --NmcmcLast ${Nmcmc} \
--em ${em} --mf ${mkfile}

############ Step 3 ###############
## Clean all jobs in the make file when you need rerun everything
make -f ${wkdir}/BFGWAS.mk clean

######### Submit the job for running the makefile
j=4 # Number of cores to request
qsub -q b.q -j y -pe smp ${j} -wd ${wkdir} -N BFGWAS ${BFGWAS_SS_dir}/bin/run_make.sh --wkdir ${wkdir} --mkfile ${mkfile} --njob ${j}

#########################################################
#########  See example Rscript for Analyzing BFGWAS results
#########################################################
# See "/1KG_example/AnalyzeResults/Analysis.r" for details of loading data and make plots
# /home/jyang/ResearchProjects/ROSMAP_GWAS/BFGWAS/Scripts/BFGWAS.r


#########################################################
############## Output under ${wkdir} ###################
#########################################################

######### /OUT/ : saves all screen outputs

######### /output/ : saves all MCMC output from E-step / will be overridden by the next E-step
# filehead.log.txt contains log file for MCMC

# filehead.paramtemp contains the following columns. function LoadEMdata() in bin/R_funcs.r can be used to read this paramtemp file into R, requiring library "data.table" and "ggplot2"
	# "CHR" : chromosome number
	# "POS" : base pair position
	# "ID" : variant ID
	# "REF" : reference allel
	# "ALT" : alternative allel
	# "MAF" : MAF of the variant
	# "AnnoFunc_code" : annotation code used in --ac FuncAnno4.txt
	# "pi" : causal probability for each variant
	# "Beta" : effect size estimate by MCMC
	# "mBeta" : Marginal effect size based on single variant test
	# "Chisq" : Chisq test statistic
	# "Pval" : Pvalue by single variant test
	# "Rank" : rank within block based on p-values, 0:top significant variant by pvalue

# filehead.hyptemp contains estimates required for M-step

# filehead.mcmc contains all included variants (id:chr:pos:ref:alt, seperated by ";") in the M-step, one row per MCMC iteration (can be used for calculating regional posterior inclusion probabilities)


######### /slurm_err/ : saves all error file for slurm jobs


###### under folder /Eoutput/ : main results ########
# let i denote the ith EM iteration
# log${i}.txt contains all log files for all blocks
# hyptemp${i}.txt contains all hyptemp files for all blocks
# paramtemp${i}.txt contains all paramtemp files for all blocks
# EM_result.txt contains all hyper parameter estimates with columns "EM_iter_#", "PVE/heritability", "likelihood", every 4 of the following columns denotes the "causal probability" "causal probability SE" "effect size variance" "effect size variance SE" for group 0, 1, 2, ...
# R function LoadEMhyp() in bin/R_funcs.r can be used to read EM_result.txt file
# R function CItable() in bin/R_funcs.r can be used to convert one row of EM_result.txt to a table of annotations


############ Additional Example Commands ###############
## Clean all jobs in the make file when you need rerun everything
make -f ${wkdir}/BFGWAS.mk clean

## Test run Mstep for one EM iteration
Rscript --vanilla ${BFGWAS_SS_dir}/bin/Mstep.r ${BFGWAS_SS_dir}/1KG_exampe/Test_Wkdir/Eoutput/hyptemp0.txt 0 1e-6 0.1 ${BFGWAS_SS_dir}/1KG_example/Test_Wkdir/Eoutput/EM_result.txt ${BFGWAS_SS_dir}/1KG_example/Test_Wkdir/hypval.current


#########################################################
######### Test Run `/bin/Estep_mcmc` with one genome block of `CFH_REGION_1KG.vcf.gz`
#########################################################
filehead=CFH_REGION_1KG

############# Command to Generate single variant GWAS Summary Stat and LD files
${BFGWAS_SS_dir}/bin/Estep_mcmc \
-vcf ${BFGWAS_SS_dir}/1KG_example/ExData/vcfs/${filehead}.vcf.gz \
-p ${BFGWAS_SS_dir}/1KG_example/ExData/phenoAMD_1KG.txt \
-GTfield GT -maf 0.001 -o ${filehead} \
-LDwindow 1000000 -saveSS -zipSS

### Run `/bin/Estep_mcmc` with Summary data
${BFGWAS_SS_dir}/bin/Estep_mcmc \
-Zscore ${BFGWAS_SS_dir}/1KG_example/Test_Wkdir/output/${filehead}.Zscore.txt.gz \
-LDcorr ${BFGWAS_SS_dir}/1KG_example/Test_Wkdir/output/${filehead}.LDcorr.txt.gz \
-a ${BFGWAS_SS_dir}/1KG_example/ExData/annos/Anno_${filehead}.txt.gz \
-fcode ${BFGWAS_SS_dir}/1KG_example/ExData/AnnoCode6.txt \
-hfile ${BFGWAS_SS_dir}/1KG_example/ExData/InitValues6.txt \
-inputSS -bvsrm -maf 0.005 -win 100 -smin 0 -smax 10 -w 10000 -s 100000 \
-o ${filehead}_SS -seed 2017

## Test run Estep for one genome block
${BFGWAS_SS_dir}/bin/Estep_mcmc -inputSS \
-Zscore ${BFGWAS_SS_dir}/1KG_example/Test_Wkdir/Zscore/${filehead}.Zscore.txt.gz \
-LDcorr ${BFGWAS_SS_dir}/1KG_example/Test_Wkdir/LDcorr/${filehead}.LDcorr.txt.gz \
-a ${BFGWAS_SS_dir}/1KG_example/ExData/annos/Anno_${filehead}.txt.gz \
-fcode ${BFGWAS_SS_dir}/1KG_example/ExData/AnnoCode6.txt \
-hfile ${BFGWAS_SS_dir}/1KG_example/Test_Wkdir/hypval.current \
-maf 0.001 -bvsrm -smin 0 -smax 5 -win 100 -o ${filehead} \
-w 10000 -s 10000 -seed 2022


################## Command to Save Genotype
BFGWAS_SS_dir="/home/jyang/GIT/BFGWAS_SS"
${BFGWAS_SS_dir}/bin/Estep_mcmc \
-vcf ${BFGWAS_SS_dir}/1KG_example/ExData/vcfs/All_4REGION_1KG.vcf.gz \
-p ${BFGWAS_SS_dir}/1KG_example/ExData/phenoAMD_1KG.txt \
-GTfield GT -maf 0.001 -o All_4REGION_1KG \
-saveGeno
mv -f output/*.geno /home/jyang/GIT/BFGWAS_SS/1KG_example/ExData/genos/
