Sys.setlocale('LC_ALL', 'C')
options(stringsAsFactors=F)
options(warn=-1)
# source("/home/jyang/GIT/bfGWAS_SS/bin/R_funcs.r")

# library(tidyverse)
library(data.table)

######## Need to pass args(hypfile, paramfile, k, hypcurrent_file) from bash
args <- commandArgs(TRUE)
# print(args)

hypfile=args[[1]] # hyper parameter values
k=as.numeric(args[[2]]) # EM iteration
pp = as.numeric(args[[3]]) # prior causal probability
a_gamma = as.numeric(args[[4]]) # prior shape parameter in Gamma
b_gamma = as.numeric(args[[5]]) # prior rate parameter in Gamma
gwas_n = as.numeric(args[[6]]) # GWAS sample size
EM_result_file = args[[7]]
hypcurrent_file=args[[8]]

###### Define functions to be used
CI_fish_pi <- function(gamma, m, a, b){
	# gamma: PIP vector
	# m : number of SNPs per category
	# a, b : hyper parameters in prior beta distribution
	pi_hat = (sum(gamma) + a - 1.0) / (m + a + b - 2.0)
	if(pi_hat <= 0 || pi_hat > 1){
		pi_hat = a/(a+b)
		se_pi = 0
	}else{
		se_pi = sqrt(pi_hat * (1-pi_hat) / (m + a + b - 2.0))
		if(se_pi < 0){
			se_pi=0
		}
	}
	return(c(pi_hat, se_pi))
}

# Posterior estimates of tau
Est_tau <- function(beta2, gamma, a, b, gwas_n){
	# beta2: squared beta estimates
	# gamma: PIP vector
	# a, b: hyper parameters in the prior inverse-gamma distribution
	tau_hat = (sum(gamma) + 2 * (a - 1)) / (gwas_n * sum(beta2 * gamma) + 2 * b)
	return(tau_hat)
}

CI_fish_tau <- function(beta2, gamma, a, b, gwas_n){
	tau_hat = Est_tau(beta2, gamma, a, b, gwas_n)
	temp = 0.5 * sum(gamma)   + (a - 1)
	if( temp > 0)
	{
		se_tau = tau_hat / sqrt( temp )
	}else{
		se_tau = 0
	}
	return(c(tau_hat, se_tau))
}

## log prior functions
logprior_tau <- function(a, b, x){ return((a - 1) * log(x) - b * x) }
logprior_pi <- function(a, b, x){ return( (a - 1) * log(x) + (b - 1) * log(1 - x)) }


ptm <- proc.time()
########### Load hypfile ....

# hypfile="/mnt/YangFSS/data2/jyang/test_BFGWAS/test_wkdir/Eoutput/hyptemp0.txt"

hypdata = read.table(hypfile, sep="\t", header=FALSE)
n_type = (dim(hypdata)[2] - 3)/3
print(paste(" Total Annotation categories : ", n_type))
temp_col_names <- c("block", "loglike", "r2")
em_out_names <- c("EM_iteration", "r2", "loglike")
for(i in 1:n_type){
	temp_col_names <- c(temp_col_names, 
			paste(c("m", "sum_gamma", "sum_Ebeta2"), (i-1), sep = "_"))
	em_out_names <- c(em_out_names, paste(c("pi_est", "pi_se", "tau_est", "tau_se"), (i-1), sep = "_"))
}
colnames(hypdata) <-  temp_col_names

########### Update hyper parameter values
# hypcurrent_file="/mnt/YangFSS/data2/jyang/test_BFGWAS/test_wkdir/hypval.current"
prehyp <- fread(hypcurrent_file, header=TRUE)
print("hyper parameter values before MCMC: ")
print(prehyp)

######### Set hierarchical parameter values
m_vec = rep(0, n_type)
for(i in 1:n_type){
	 m_vec[i] <- sum(hypdata[, paste("m", (i-1), sep="_")])
}

#### updating hyper pi and tau_beta values for each group
hypcurrent <- NULL
hypmat <- NULL
sum_gamma = 0
for(i in 1:n_type){
	# print(i)
	if(m_vec[i] > 0){
		a_beta = m_vec[i] * pp; b_beta = m_vec[i] - a_beta;
		}else{a_beta=1; b_beta = (1/pp) - 1;}
	sum_gamma_temp = hypdata[, paste("sum_gamma", (i-1), sep="_")]
	sum_gamma = sum_gamma + sum(sum_gamma_temp)
	pi_temp = CI_fish_pi(sum_gamma_temp, m_vec[i], a_beta, b_beta)
	tau_temp = c(prehyp[i, 2], 0) # do not update tau_beta
	#CI_fish_tau(hypdata[, paste("sum_Ebeta2", (i-1), sep="_")], sum_gamma_temp,  a_gamma, b_gamma, gwas_n)
	hypcurrent <- c(hypcurrent, pi_temp, tau_temp)
	# print(cbind(pi_temp, tau_temp))
	hypmat <- rbind(hypmat, c(pi_temp[1], tau_temp[1]))
}

########## Write out updated hyper parameter values
colnames(hypmat) <- c("pi", "tau")
print("hyper parameter values updates after MCMC: ")
print(hypmat)
write.table(format(data.frame(hypmat), scientific=TRUE, digits = 3),
	file=hypcurrent_file, 
	quote = FALSE, sep = "\t", row.names=FALSE, col.names=TRUE)

#### Summarize log-likelihood
loglike_total = sum(hypdata$loglike, na.rm = TRUE)
	for(i in 1:n_type){
		if(m_vec[i] > 0){
			a_beta = m_vec[i] * pp; b_beta = m_vec[i] - a_beta;
		}else{a_beta=1; b_beta = (1/pp) - 1;}
		loglike_total = loglike_total +
				logprior_pi(a_beta, b_beta, prehyp[i, 1])
}

########## Write out updated hyper parameter values and se to EM_result_file
# EM_result_file="/mnt/YangFSS/data2/jyang/test_BFGWAS/test_wkdir/Eoutput/EM_result.txt"
pve = sum(hypdata[, "r2"], na.rm = TRUE)
print(paste("Sum PIP = ", sum_gamma))
print(paste("Regression R2 = ", pve))
print(paste("Posterior log likelihood = ", loglike_total))

hypcurrent = c(pve, loglike_total, hypcurrent)
hypcurrent <- format(hypcurrent, scientific = TRUE)
if(k==0){
	write.table(matrix(c(k, hypcurrent), nrow=1, dimnames = list(NULL, em_out_names)), file = EM_result_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, append=FALSE)
}else{
	write.table(matrix(c(k, hypcurrent), nrow=1), file = EM_result_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append=TRUE)
}


print("EM step time cost (in minutes) : ")
print((proc.time() - ptm)/60)








