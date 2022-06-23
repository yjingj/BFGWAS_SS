# rm(list=ls(all=TRUE))
library(tidyverse)

####### Source Util Functions and set data directories
source("/home/jyang/GIT/BFGWAS_SS/bin/R_funcs.r")

DataDir = "/home/jyang/GIT/BFGWAS_SS/1KG_example/ExData/" # example data directory
OutDir = "/home/jyang/GIT/BFGWAS_SS/1KG_example/AnalyzeResults/" # directory to save plots
ResultDir = "/home/jyang/GIT/BFGWAS_SS/1KG_example/Test_Wkdir" # result directory

setwd("/home/jyang/GIT/BFGWAS_SS/1KG_example/AnalyzeResults/")

## True phenotype : "/net/fantasia/home/yjingj/GIT/SFBA_example/ExData/phenoAMD_1KG.txt"
OR = read.table(file = "/home/jyang/GIT/BFGWAS_SS/1KG_example/ExData/vcfs/causalSNP_OR.txt", header = TRUE)
paramdata_bfgwas[OR$rsID, ]

######## Compare results
paramdata_bfgwas = LoadEMdata(filename="/home/jyang/GIT/BFGWAS_SS/1KG_example/Test_Wkdir/Eoutput/paramtemp3.txt", header = TRUE)
#head(paramdata_bfgwas)
sum(paramdata_bfgwas$Pi)

## Manhantton plot
paramdata_bfgwas_sig <- filter(paramdata_bfgwas, Pi > 0.1)
dim(paramdata_bfgwas_sig )
paramdata_bfgwas_sig

ggplot(paramdata_bfgwas, aes(x=POS, y = -log10(Pval), color = Pi)) +
	geom_point() +scale_color_gradient(low="blue", high="red") +
	facet_grid(cols = vars(CHR), scales = "free") +
	# geom_point(paramdata_bfgwas_sig, aes(x=POS, y = -log10(Pval), color = Pi)) +
	geom_hline(yintercept=-log10(5e-8))
ggsave("/home/jyang/GIT/BFGWAS_SS/1KG_example/AnalyzeResults/mp_06_2022.pdf")

## Load Phenotype
library(data.table)
library(tidyverse)

pheno = read.table("/home/jyang/GIT/BFGWAS_SS/1KG_example/ExData/phenoAMD_1KG.txt", header = FALSE)
y = pheno$V2
y = scale(y)
names(y) = pheno$V1

### Read all genotype Data
geno_data <- fread("/home/jyang/GIT/BFGWAS_SS/1KG_example/ExData/genos/All_4REGION_1KG.geno", sep = "\t", header = TRUE)
setkey(geno_data, "ID")

Zscore <- fread("/home/jyang/GIT/BFGWAS_SS/1KG_example/Test_Wkdir/Zscore/All_region.Zscore.txt", sep = "\t", header = TRUE)
setkey(Zscore, "ID")

######
paramdata_bfgwas = LoadEMdata(filename="/home/jyang/GIT/BFGWAS_SS/1KG_example/Test_Wkdir/Eoutput/paramtemp3.txt", header = TRUE)
#head(paramdata_bfgwas)
sum(paramdata_bfgwas$Pi)

SNP_vec <- paramdata_bfgwas[paramdata_bfgwas$Pi > 0.001, ]$ID
length(SNP_vec)
X <- geno_data[SNP_vec, -c(1:5)] %>% t()
colnames(X) <- SNP_vec

fit = lm(y ~ ., data = data.frame(y=y, X))
summary(fit)
summary(fit)$adj.r.squared # Adjusted R-squared: 0.2536693

beta <- paramdata_bfgwas[SNP_vec, ]$Beta
pred_y <- as.matrix(X) %*% beta
cor(pred_y[, 1], y)^2 #  0.2966076

ggplot(paramdata_bfgwas[SNP_vec, ], aes(x = mBeta, y = Beta, col = Pi)) +
	geom_point() + geom_abline(intercept=0, slope = 1) +
	scale_color_gradient(low="blue", high="red")
ggsave("/home/jyang/GIT/BFGWAS_SS/1KG_example/AnalyzeResults/beta_06_2022.pdf")

SNP_vec <- paramdata_bfgwas[paramdata_bfgwas$Pval < 1e-5, ]$ID
X <- geno_data[SNP_vec, -c(1:5)] %>% as.matrix()
summary(lm(y ~ t(X)))$adj.r.squared # Adjusted R-squared:  0.1966276
beta <- paramdata_bfgwas[SNP_vec, ]$mBeta
pred_y <- t(X) %*% beta
cor(pred_y[, 1], y)^2 # 0.1318686

SNP_vec <- OR$rsID[c(1, 3, 4, 6, 7, 9, 10, 11, 12)]
summary(lm(y ~ t(X)))$adj.r.squared # Adjusted R-squared:  0.1966276
X <- geno_data[SNP_vec, -c(1:5)] %>% as.matrix()
beta <- paramdata_bfgwas[SNP_vec, ]$Beta
pred_y <- t(X) %*% beta
cor(pred_y[, 1], y)^2 #  0.07832969

################################################
###### REMARK: you might want to check if the estimated variants with high association probability are in high LD with the true causal ones

###### Although we are not able to identify all causal variants, we can estimate the probability that there exist at least one causal variant per genome block
system("cat \`ls /home/jyang/GIT/BFGWAS_SS/1KG_example/Test_Wkdir/output/** | grep log\` | grep Region_PIP | cut -d\" \" -f4 > /home/jyang/GIT/BFGWAS_SS/1KG_example/Test_Wkdir/Eoutput/Region_PP.txt")
Region_PP = scan("/home/jyang/GIT/BFGWAS_SS/1KG_example/Test_Wkdir/Eoutput/Region_PP.txt", what = double())
sum(Region_PP > 0.95)

######## Load results of hyper parameters ####### 
test_hyp <- LoadEMhyp(filename = "/home/jyang/GIT/BFGWAS_SS/1KG_example/Test_Wkdir/Eoutput/EM_result.txt", header = TRUE)

# the labels should be consistant to 0,1,2... defined in the annotation code file
group_labels <- as.factor(c("Coding", "UTR", "Promoter", "DHS", "Intronic", "Intergenic"))

test_CI_table <- CItable(test_hyp[nrow(test_hyp), ], n_type = 6, alpha = 0.95,
		funcgroup = group_labels)

# Set the values of these categories without associations at NAs, change prior_pp value accordingly
prior_pp = 1e-6
test_CI_table [test_CI_table$pi == prior_pp, 1:6] <- NA 
 
# plot causal proportion estimates, requiring library "ggplot2"
PlotCI_groupPP(hyp_table=test_CI_table, pdfname="/home/jyang/GIT/BFGWAS_SS/1KG_example/AnalyzeResults/test_groupPP.pdf", size = 18, tit = "")


######## Compare hyper estimates to the genome-wide averages
n_group = 6
pp_cols = (1:n_group) * 4
group_pp = as.vector( test_hyp[6, (pp_cols)])
no_asso_groups = (group_pp == prior_pp)
group_pp[no_asso_groups] = 0
pp_se <- as.vector( test_hyp[6, (pp_cols + 1)])
pp_se[no_asso_groups] = 0

# as.vector(table(paramdata$func))

n_vec = c(9221, 203, 0, 1354, 0, 5627) # number of variants per annotation, 
# should be consistant to 0,1,2... defined in the annotation code file

ncausal_vec = n_vec * group_pp # number of associations
group_sigma2 <- as.vector( test_hyp[6, (pp_cols + 2)])
group_sigma2[no_asso_groups] = 0
sigma2_se <- as.vector( test_hyp[6, (pp_cols + 3)])
sigma2_se[no_asso_groups] = 0

### Construct comparison results
comp_group_sigma2 <-comp_group_pp <- data.frame(log_ratio_lcl = rep(NA, n_group), 
							log_ratio=rep(NA, n_group), 
							log_ratio_ucl = rep(NA, n_group))

for(i in 1:n_group){
		comp_group_pp[i, ] = comp_groupEst(est = group_pp, est_se = pp_se, n_vec = n_vec, i = i, conf = 0.95)
		comp_group_sigma2[i, ] = comp_groupEst(est = group_sigma2, est_se = sigma2_se, n_vec = ncausal_vec, i = i, conf = 0.95)
}

exp(comp_group_pp[, 2])
exp(comp_group_sigma2[, 2])


### Plot comparison results
PlotRatio(comp_dat = data.frame(exp(comp_group_pp), funcgroup=group_labels), 
	tit ="Regulatory Annotation", 
	pdfname = paste(OutDir, "comp_test_pp.pdf", sep =""), 
	ymode = 1, size = 28, wid = 8) 


PlotRatio(comp_dat = data.frame(exp(comp_group_sigma2), funcgroup=group_labels), 
	tit ="Regulatory Annotation", 
	pdfname = paste(OutDir, "comp_test_sigma2.pdf", sep =""), 
	ymode = 0, size = 28, wid = 8) 

######## END ################










