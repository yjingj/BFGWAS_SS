######################### Test on 1KG Example Data
############## Run BFGWAS_SS with multiple genome block
# Specify BFGWAS tool directory and working directory with writing access
BFGWAS_SS_dir="/home/jyang/GIT/BFGWAS_SS"
wkdir="/home/jyang/GIT/BFGWAS_SS/1KG_example/Test_Wkdir"

# file with all genome block file name heads
filehead=${BFGWAS_SS_dir}/1KG_example/ExData/fileheads_4region.txt

# genotype directory
geno_dir=${BFGWAS_SS_dir}/1KG_example/ExData/vcfs
# phenotype directory
pheno=${BFGWAS_SS_dir}/1KG_example/ExData/phenoAMD_1KG.txt

####### First generate single variant GWAS score statistic and LD files
# LD coefficient directory
LDdir=${wkdir}/LDcorr
# Score statistics directory
Score_dir=${wkdir}/ScoreStat
LDwindow=500000

cd $wkdir
${BFGWAS_SS_dir}/bin/GetRefLD.sh --wkdir ${wkdir} \
--toolE ${BFGWAS_SS_dir}/bin/Estep_mcmc \
--geno_dir ${geno_dir} --filehead ${filehead} --pheno ${pheno} \
--GTfield GT --genofile_type vcf --maf 0 --LDwindow ${LDwindow}  \
--Score_dir ${Score_dir} --LDdir ${LDdir} &


########## Generate make file
N=2504 ; # sample size
pheno_var=0.250064; # phenotype variance
maf=0.005
em=3; # EM steps
Nburnin=100; # Burn-in iterations in MCMC
Nmcmc=100; # MCMC iteration number
mkfile=${wkdir}/BFGWAS.mk

anno_dir=${BFGWAS_SS_dir}/1KG_example/ExData/annos
anno_code=${BFGWAS_SS_dir}/1KG_example/ExData/AnnoCode6.txt
hfile=${BFGWAS_SS_dir}/1KG_example/ExData/InitValues6.txt

cd $wkdir
${BFGWAS_SS_dir}/bin/gen_mkf.pl \
--wkdir ${wkdir} --tool_dir ${BFGWAS_SS_dir} \
--LDdir ${LDdir} --Score_dir ${Score_dir} \
--filehead ${filehead} --anno_dir ${anno_dir} --anno_code ${anno_code} \
--hfile ${hfile} --pv ${pheno_var} --Nsample ${N} --maf ${maf} \
--Nburnin ${Nburnin} --Nmcmc ${Nmcmc} --NmcmcLast ${Nmcmc} \
--em ${em} --mf ${mkfile}

######### Submit the job for running the makefile
j=4 # Number of cores to request
qsub -q b.q -j y -pe smp ${j} -wd ${wkdir} -N BFGWAS ${BFGWAS_SS_dir}/bin/run_make.sh --wkdir ${wkdir} --mkfile ${mkfile} --njob ${j}

######### update the following 09/10/2021 ####
#########  See example Rscript for Analyzing BFGWAS results
/home/jyang/ResearchProjects/ROSMAP_GWAS/BFGWAS/Scripts/BFGWAS.r

##### perl script is used to generate makefile #############

${BFGWAS_dir}/bin/gen_mkf.pl \
-w ${wkdir} \
--Estep ${BFGWAS_dir}/bin/Estep_mcmc \
--ad ${BFGWAS_SS_dir}/1KG_example/ExData/annos \
--geno vcf --ac ${BFGWAS_SS_dir}/1KG_example/ExData/AnnoCode6.txt  \
--gd ${BFGWAS_SS_dir}/1KG_example/ExData/vcfs \
--pheno ${BFGWAS_SS_dir}/1KG_example/ExData/phenoAMD_1KG.txt \
--hyp ${BFGWAS_SS_dir}/1KG_example/ExData/InitValues6.txt \
-f ${BFGWAS_SS_dir}/1KG_example/ExData/fileheads_4region.txt \
--rs ${BFGWAS_dir}/bin/Mstep.r \
-G GT --maf 0.005 --smax 10 -b 10000 -N 10000 --NL 10000 \
--em 3 -j testjob -l local \
--mf ${wkdir}/Test_BFGWAS.mk

## Run Makefile with 4 parallel jobs
make -k -C ${wkdir} -f ${wkdir}/Test_BFGWAS.mk -j 10 > ${wkdir}/make.output 2> ${wkdir}/make.err  &

## Clean all jobs when you need rerun everything
# make -f ${wkdir}/Test_BFGWAS.mk clean


#################
############## Output under ${wkdir} ###################

######### /OUT/ : saves all screen outputs

######### /output/ : saves all MCMC output from E-step / will be overridden by the next E-step
# filehead.log.txt contains log file for MCMC
# filehead.paramtemp contains estimates for each variant with columns "ID", "chr", "bp", "ref", "alt", "maf", "func", "beta", "pi", "Zscore", "SE_beta", "LRT", "pval_LRT", "rank";

# "ID" : variant ID
# "chr" : chromosome number
# "bp" : base pair position
# "REF" : reference allel
# "ALT" : alternative allel
# "maf" : MAF of the variant
# "func" : annotation code used in --ac FuncAnno4.txt
# "beta" : effect size estimate
# "pi" : causal probability for each variant
# "Zscore" : Zscore by single likelihood ratio test
# "SE_beta" : from single likelihood ratio test
# "LRT" : test statistic by single likelihood ratio test
# "pval_LRT" : pvalue by single likelihood ratio test
# "rank" : rank within block based on p-values, 0:top significant variant by pvalue

# function LoadEMdata() in bin/R_funcs.r can be used to read this paramtemp file into R, requiring library "data.table" and "ggplot2"

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


############ example R code for analysis
# see "/1KG_example/AnalyzeResults/Analysis.r" for details of loading data and make plots



### Test Run `/bin/Estep_mcmc` with individual level data and ANNO file for the genome block of `CFH_REGION_1KG.vcf.gz`
filehead=CFH_REGION_1KG

#### change `-GTfield` argument to the dosage field name if imputed genotype data are used
#### Save single variant test GWAS summary data and LD  with flag `-saveSS`
${BFGWAS_SS_dir}/bin/Estep_mcmc \
-vcf ${BFGWAS_SS_dir}/1KG_example/ExData/vcfs/${filehead}.vcf.gz \
-p ${BFGWAS_SS_dir}/1KG_example/ExData/phenoAMD_1KG.txt \
-a ${BFGWAS_SS_dir}/1KG_example/ExData/annos/Anno_${filehead}.txt.gz \
-fcode ${BFGWAS_SS_dir}/1KG_example/ExData/AnnoCode6.txt \
-hfile ${BFGWAS_SS_dir}/1KG_example/ExData/InitValues6.txt \
-GTfield GT -bvsrm -maf 0.005 -win 100 -smin 0 -smax 10 -w 10000 -s 100000 \
-o ${filehead} -initype 3 -LDwindow 500000 -seed 2017 \
-saveSS -zipSS > ${wkdir}/${filehead}_indv.output.txt &

############# Command to Generate single variant GWAS Summary Stat and LD files
filehead=C2_REGION_1KG
${BFGWAS_SS_dir}/bin/Estep_mcmc \
-vcf ${BFGWAS_SS_dir}/1KG_example/ExData/vcfs/${filehead}.vcf.gz \
-p ${BFGWAS_SS_dir}/1KG_example/ExData/phenoAMD_1KG.txt \
-GTfield GT -maf 0.005 -o ${filehead} \
-LDwindow 500000 -saveSS -zipSS > ${wkdir}/${filehead}_generateSS.output.txt &


### Run `/bin/Estep_mcmc` with Summary data
${BFGWAS_SS_dir}/bin/Estep_mcmc \
-score ${BFGWAS_SS_dir}/1KG_example/Test_Wkdir/output/${filehead}.score.txt.gz \
-LDcorr ${BFGWAS_SS_dir}/1KG_example/Test_Wkdir/output/${filehead}.LDcorr.txt.gz \
-a ${BFGWAS_SS_dir}/1KG_example/ExData/annos/Anno_${filehead}.txt.gz \
-fcode ${BFGWAS_SS_dir}/1KG_example/ExData/AnnoCode6.txt \
-hfile ${BFGWAS_SS_dir}/1KG_example/ExData/InitValues6.txt \
-n 2504 -pv 0.250064 -inputSS \
-bvsrm -maf 0.005 -win 100 -smin 0 -smax 10 -w 10000 -s 100000 \
-o ${filehead}_SS -initype 3 -seed 2017 > ${filehead}_SS.output.txt &


################## Command to Save Genotype
${BFGWAS_SS_dir}/bin/Estep_mcmc \
-vcf ${BFGWAS_SS_dir}/1KG_example/ExData/vcfs/${filehead}.vcf.gz \
-p ${BFGWAS_SS_dir}/1KG_example/ExData/phenoAMD_1KG.txt \
-GTfield GT -maf 0 -o ${filehead} \
-saveGeno > ${wkdir}/${filehead}_saveGeno.output.txt &