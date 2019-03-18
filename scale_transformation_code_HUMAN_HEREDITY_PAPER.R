#############################################################
#############################################################
### Transforming summary statistics from logistic ###########
### regression to the liability scale: Application to #######
### genetic and environmental risk scores ###################
#############################################################
### Human Heredity. DOI: 10.1159/000495697 ##################
#############################################################
### Contact details: alex.gillett@icloud.com ################
### Please contact me with any problems or if you want to ###
### use this method for your work and require assistance ####
#############################################################
### ORDER OF SCRIPT #########################################
#############################################################
### 1. FIGURE 1: Normal and standardised Logistic CDFs
### 2.1. FIGURE 2: Log OR versus liability effect size for a 
# binary risk factor, for a range of disease prevalences and 
# risk factor PDFs
### 2.2. FIGURES A1 and A2: Difference between equations 4
# and 5, and 4 and 6, respectively, for a binary risk factor, 
# for a range of disease prevalences and risk factor PDFs
### 3. SCHIZOPHRENIA RISK EXAMPLE
### 3.1. Liability threshold models (LTMs)
### 3.1.1. General code and functions
### 3.1.2. Parameterising PRS only, ERS only and joint LTMs
### 3.1.3. Risk estimation; all LTMs
### 3.1.4. FIGURE 3: Schizophrenia risk; all LTMs
### 3.1.5. Specific risk examples from text in paper 
### 3.2. Logistic ERS only model
### 3.2.1. General code and functions
### 3.2.2. Parameterising model
### 3.1.3. Risk estimation ERS only logistic model
### 3.2.4. FIGURE 4: ERS only logistic model expected risk 
# versus ERS only liability model expected risk
### 3.3. Additional possible approximate logistic models 
# (presented in the Appendix of the paper)
### Sincere apologies that the code for this section is 
# not easy to follow, and irregularly named. I aim to change 
# this after my maternity leave (by start of 2020)
### 3.3.1. Possible approximate PRS only logistic model
### 3.3.2. Possible approximate joint PRS and ERS logistic 
# model
### 3.3.3. FIGURE A3. Estimated liability thresmodel model
# risks against estimated logistic model risks for ERS only, 
# PRS only and, PRS and ERS models
### NOTE!!!
### PLEASE DOWNLOAD AND SAVE FILES:
### A. scz_risk_ers_prs.txt 
# Contains the expected schizophrenia risks for the joint PRS and
# ERS liability threshold model (downloaded the pre-calculated risks
# will save the reader time)
### B. scz_prs_snps.txt
# Originally downloaded from the PGC website, contained the ORs
# for the SNPs used in the PRS analysis (therefore contains 
# more SNPs than the final PRS)
### C. SCZ_1000G_AF.txt
# Originally obtained from the 1000G project website, contains
# the MAF for the PRS SNPs in the European 1000G population
# This is available from:
# https://github.com/alexgillett/scale_transformation
#############################################################
### Libraries needed and source other scripts ###############
#############################################################
### INSTALL REQUIRED LIBRARIES:
library(BB) # dfsane function for solving simultaneous eqns
library(ggplot2)
library(data.table)
#############################################################
### 1. FIGURE 1 #############################################
#############################################################
### Normal and standardised Logistic CDFs
xrange <- seq(from = -3, to = 3, by = 0.01)
CDF.norm <- pnorm(xrange)
CDF.logis <- plogis(xrange, s=sqrt(3)/pi)
plot(x=xrange, y=CDF.norm, bty="l", xlab="x", ylab="Probability", type="l", lwd=2)
points(x=xrange, y=CDF.logis, col="dark grey", type="l", lty=2, lwd=2)
leg.txt <- c("Standard normal CDF", "Standardised logistic CDF")
legend("topleft", legend=leg.txt, col=c("black", "dark grey"), lty=c(1,2), lwd=c(2,2), text.col=c("black", "dark grey"), bty="n")
#############################################################
### 2.1. FIGURE 2 ###########################################
#############################################################
### Log OR versus liability effect size for a binary risk 
# factor, for a range of disease prevalences and risk factor 
# PDFs
### Input values:
K.vec <- c(1/10000, 1/1000, 1/100, 1/10, 1/5)
OR.vec <- c(1.001, 1.005, 1.01, 1.05, 1.1, 1.5, 2, 5, 10)
logOR.vec <- log(OR.vec)
pdfRF.vec <- c(1/10000, 1/100, 1/5)
### Scale transformation functions
### Equation 3 in paper is equivalent to equation 4 for a binary RF
### Eqns 4, 5 and 6 function:
logit2liab.f <- function(OR, K, prob){
	functionb0.e4 <- function(p){
		r <- rep(NA, length(p))
		r[1] <- ((1-prob)*(1/(1 + exp(-p[1])))) + (prob*(1/(1 + exp(-p[1])*exp(-log(OR))))) - K
		r
	}
	root.out <- uniroot(functionb0.e4, c(-100,100))
	b0.e4 <- root.out$root
	pDlogit.RF1.e4 <- 1/(1 + (exp(-b0.e4)*exp(-log(OR))))
	pDlogit.RF0.e4 <- 1/(1 + (exp(-b0.e4)))
	b0.e6 <- log(K/(1-K))
	pDlogit.RF1.e6 <- 1/(1 + (exp(-b0.e6)*exp(-log(OR))))
	pDlogit.RF0.e6 <- 1/(1 + (exp(-b0.e6)))
	varRF <- prob*(1-prob)
	a <- qnorm(pDlogit.RF1.e4) - qnorm(pDlogit.RF0.e4)
	tau.e4 <- a/sqrt(1 + (a^2)*varRF)
	tau.e5 <- a
	tau.e6 <- qnorm(plogis(log(K/(1-K)) + log(OR))) - qnorm(K)
	out <- c(b0.e4, b0.e6, pDlogit.RF0.e4, pDlogit.RF1.e4, pDlogit.RF0.e6, pDlogit.RF1.e6, tau.e4, tau.e5, tau.e6)
	out
}
### Generation of output data
logit2liab.out <- NULL
names.logit2liab.out <- c("K", "OR", "logOR", "prob", "b0.e4", "b0.e6", "pDlogit.RF0.e4", "pDlogit.RF1.e4", "pDlogit.RF0.e6", "pDlogit.RF1.e6", "tau.e4", "tau.e5", "tau.e6")
for(i in 1:length(K.vec)){
	for(j in 1:length(OR.vec)){
		for(k in 1:length(pdfRF.vec)){
			out.ijk <- c(K.vec[i], OR.vec[j], logOR.vec[j], pdfRF.vec[k])
			work.out.ijk <- logit2liab.f(OR=OR.vec[j], K=K.vec[i], prob= pdfRF.vec[k])
			out.ijk <- c(out.ijk, work.out.ijk)
			out.ijk <- t(matrix(out.ijk))
			logit2liab.out <- rbind(logit2liab.out, out.ijk)
			rm(out.ijk, work.out.ijk)
		}
	}
}
colnames(logit2liab.out) <- names.logit2liab.out
logit2liab.out <- data.frame(logit2liab.out)
K_prob <- rep("", dim(logit2liab.out)[1])
for(i in 1:dim(logit2liab.out)[1]){
	K_prob[i] <- paste(as.character(logit2liab.out$K[i]), as.character(logit2liab.out$prob[i]), sep="_")
}
logit2liab.out$K_prob <- K_prob
logit2liab.out$d.e4e5 <- logit2liab.out$tau.e4 - logit2liab.out$tau.e5
logit2liab.out$d.e4e6 <- logit2liab.out$tau.e4 - logit2liab.out$tau.e6
# Change inputs (K and PDF values) into factors for ease of plotting
K.f <- rep("", dim(logit2liab.out)[1])
prob.f <- rep("", dim(logit2liab.out)[1])
K.f[logit2liab.out$K == 1/1000] <- "1/1000"
K.f[logit2liab.out$K == 1/100] <- "1/100"
K.f[logit2liab.out$K == 1/5] <- "1/5"
K.f[logit2liab.out$K == 1/10000] <- "1/10000"
K.f[logit2liab.out$K == 1/10] <- "1/10"
K.f <- factor(K.f, levels=c("1/5","1/10", "1/100","1/1000", "1/10000"), ordered=T)
logit2liab.out$Prevalence <- K.f
prob.f[logit2liab.out$prob==1/10000] <- "p(X=1)=1/10000"
prob.f[logit2liab.out$prob==1/100] <- "p(X=1)=1/100"
prob.f[logit2liab.out$prob==1/5] <- "p(X=1)=1/5"
prob.f <- factor(prob.f, levels = c("p(X=1)=1/10000","p(X=1)=1/100","p(X=1)=1/5"), ordered=T)
logit2liab.out$prob.f <- prob.f
### Figure 2
p2 <- ggplot(logit2liab.out, aes(x=logOR, y=tau.e4, colour=Prevalence)) + geom_point() + geom_line(lwd=1) + xlab("Log odds ratio") + ylab("Effect size estimate on the liability scale") + scale_color_grey() +theme_bw()
p2 + facet_grid(prob.f~.)
#############################################################
### 2.2. FIGURES A1 and A2 ##################################
#############################################################
# FIGURES A1 and A2: Difference between equations 4
# and 5, and 4 and 6, respectively, for a binary risk factor, 
# for a range of disease prevalences and risk factor PDFs
pA1 <- ggplot(logit2liab.out, aes(x=logOR, y=d.e4e5, colour=Prevalence)) + geom_point() + geom_line(lwd=1) + xlab("Log odds ratio") + ylab("Difference in liability effect size: Equations (4) - (5)") + scale_color_grey() +theme_bw()
pA1 + facet_grid(prob.f~.)
pA2 <- ggplot(logit2liab.out, aes(x=logOR, y=d.e4e6, colour=Prevalence)) + geom_point() + geom_line(lwd=1) + xlab("Log odds ratio") + ylab("Difference in liability effect size: Equations (4) - (6)") + scale_color_grey() +theme_bw()
pA2 + facet_grid(prob.f~.)
#############################################################
#############################################################
### 3. SCHIZOPHRENIA RISK EXAMPLE ###########################
#############################################################
#############################################################
### 3.1. Liability threshold models (LTMs) ##################
#############################################################
### 3.1.1. General code and functions to be used ############
#############################################################
### Joint risk function: ERS and PRS in model
risk_function <- function(ers, prs, E_ERS, V_ERS, E_PRS=0, V_PRS, K){
	E_L <- E_PRS + E_ERS
	L_T <- E_L - qnorm(K)
	denom <- sqrt(1 - V_ERS - V_PRS)
	risk <- pnorm(-(L_T - ers - prs)/denom)
	risk
}
### Risk function: only ERS:
risk_single_ERS_function <- function(ers, E_ERS, E_PRS = 0, V_ERS, K){
	E_L <- E_PRS + E_ERS
	L_T <- E_L - qnorm(K)
	denom <- sqrt(1 - V_ERS)
	risk <- pnorm(-(L_T - ers)/denom)
	risk
}
### Risk function: only PRS:
risk_single_PRS_function <- function(prs, E_ERS, E_PRS = 0, V_PRS, K){
	E_L <- E_PRS + E_ERS
	L_T <- qnorm(1-K)
	denom <- sqrt(1 - V_PRS)
	risk <- pnorm(-(L_T - prs)/denom)
	risk
}
#############################################################
### 3.1.2. Parameterising PRS only, ERS only and joint LTMs #
#############################################################
### General known parameters:
###################################################
K <- 0.01 ### prevalence
H2 <- 0.81 ### Broad-sense heritability
h2 <- 0.81 ### Narrow-sense heritability
VML <- 0.07 ### Var[PRS] on L-scale
VD <- H2-h2 ### (Quasi) Dominant variance component
###################################################
### Environmental risk variables:
###################################################
### For each of the 5 environmental risk factors 
# calculate the effect size estimates on the liability
# scale (denoted tau)
###################################################
### E1 = cannabis usage 
### Number of categories contained within E1
levels_E1 <- 3
### probability distribution function (known)
pvecE1 <- c(0.70, 0.15, 0.15)
### OR estimates (known)
OR_E1_0 <- 1
OR_E1_1 <- 1.41
OR_E1_2 <- 2.78
lOR_E1_1 <- log(OR_E1_1)
lOR_E1_2 <- log(OR_E1_2)
### Penetrance function calculated using K and OR
functionb0logit_E1 <- function(p){
	r <- rep(NA, length(p))
	r[1] <- (pvecE1[1]*(1/(1 + exp(-p[1])))) + (pvecE1[2]*(1/(1 + exp(-p[1])*exp(-lOR_E1_1)))) + (pvecE1[3]*(1/(1 + exp(-p[1])*exp(-lOR_E1_2)))) - K
	r
}
functionb0logit_E1_2 <- function(p){
	(pvecE1[1]*(1/(1 + exp(-p)))) + (pvecE1[2]*(1/(1 + exp(-p)*exp(-lOR_E1_1)))) + (pvecE1[3]*(1/(1 + exp(-p)*exp(-lOR_E1_2)))) - K
}
p0 <- rep(0.001, 1)
outE1 <- dfsane(par=p0, fn=functionb0logit_E1, control=list(trace=F))
outE1
#spg(par=p0, fn=functionb0logit_E1_2, control=list(ftol = 1.e-16, gtol=1.e-16))
uniroot(functionb0logit_E1_2, c(-10,10))
b0logitE1 <- outE1$par[1]
#pD_E1_0 <- 0.007562538
pD_E1_0 <- 1/(1 + exp(-1*b0logitE1))
1/(1 + exp(-1*-12.00201))
#pD_E1_1 <- 0.010630218
pD_E1_1 <- 1/(1 + exp(-1*b0logitE1)*exp(-1* lOR_E1_1))
#pD_E1_2 <- 0.020744605
pD_E1_2 <- 1/(1 + exp(-1*b0logitE1)*exp(-1* lOR_E1_2))
b1logitE1 <- lOR_E1_1
b2logitE1 <- lOR_E1_2
### 
tau_work <- c(b0logitE1, b0logitE1+b1logitE1, b0logitE1+b2logitE1)
tau_work <- qnorm(plogis(tau_work)) - qnorm(plogis(b0logitE1))
C <- sum((tau_work^2)*pvecE1) - sum(as.matrix(tau_work*pvecE1)%*%t(as.matrix(tau_work*pvecE1)))
tauE1vec <- sqrt(1/(1+C))*tau_work
CE1 <- C
VE1 <- CE1/(1+CE1)
tau0E1 <- tauE1vec[1]
tau1E1 <- tauE1vec[2]
tau2E1 <- tauE1vec[3]
E_E1 <- sum(tauE1vec*pvecE1)
###################################################
### E2 = Migration
levels_E2 <- 2
pvecE2 <- c(0.924, 1-0.924)
OR_E2_0 <- 1
OR_E2_1 <- 2.3
lOR_E2_1 <- log(OR_E2_1)
### Penetrance function calculated using K and OR
functionb0logit_E2 <- function(p){
	r <- rep(NA, length(p))
	r[1] <- (pvecE2[1]*(1/(1 + exp(-p[1])))) + (pvecE2[2]*(1/(1 + exp(-p[1])*exp(-lOR_E2_1)))) - K
	r
}
p0 <- rep(0.001, 1)
outE2 <- dfsane(par=p0, fn=functionb0logit_E2, control=list(trace=F))
outE2
b0logitE2 <- outE2$par[1]
pD_E2_0 <- 1/(1 + exp(-1*b0logitE2))
pD_E2_1 <- 1/(1 + exp(-1*b0logitE2)*exp(-1*lOR_E2_1))
#b0logitE2 <- log(pD_E2_0/(1-pD_E2_0))
b1logitE2 <- lOR_E2_1
b0probitE2 <- qnorm(plogis(b0logitE2))
b1probitE2 <- qnorm(plogis(b0logitE2 + b1logitE2)) - b0probitE2
tau0E2 <- 0
### Calculate effect size estimate on the L-scale. Note 
# different method used compared to E1. Both give 
# the same answer
tau1E2 <- b1probitE2/sqrt(1 + ((b1probitE2^2)*pvecE2[2]*(1-pvecE2[2])))
tauE2vec <- c(tau0E2, tau1E2)
VE2 <- ((tau1E2^2)*pvecE2[2]*(1-pvecE2[2]))
E_E2 <- sum(tauE2vec*pvecE2)
###################################################
### E3 = Urbanicity
levels_E3 <- 3
pvecE3 <- c(0.25, 0.25, 0.5)
OR_E3_0 <- 1
OR_E3_1 <- 1.5
OR_E3_2 <- 2
lOR_E3_1 <- log(OR_E3_1)
lOR_E3_2 <- log(OR_E3_2)
### Penetrance function calculated using K and OR
functionb0logit_E3 <- function(p){
	r <- rep(NA, length(p))
	r[1] <- (pvecE3[1]*(1/(1 + exp(-p[1])))) + (pvecE3[2]*(1/(1 + exp(-p[1])*exp(-lOR_E3_1)))) + (pvecE3[3]*(1/(1 + exp(-p[1])*exp(-lOR_E3_2)))) - K
	r
}
p0 <- rep(0.001, 1)
outE3 <- dfsane(par=p0, fn=functionb0logit_E3, control=list(trace=F))
outE3
b0logitE3 <- outE3$par[1]
pD_E3_0 <- 1/(1 + exp(-1*b0logitE3))
pD_E3_1 <- 1/(1 + exp(-1*b0logitE3)*exp(-1*lOR_E3_1))
pD_E3_2 <- 1/(1 + exp(-1*b0logitE3)*exp(-1*lOR_E3_2))
###
b1logitE3 <- lOR_E3_1
b2logitE3 <- lOR_E3_2
b0probitE3 <- qnorm(plogis(b0logitE3))
b1probitE3 <- qnorm(plogis(b0logitE3 + b1logitE3)) - b0probitE3
b2probitE3 <- qnorm(plogis(b0logitE3 + b2logitE3)) - b0probitE3
###
tau_work <- c(b0logitE3, b0logitE3+b1logitE3, b0logitE3+b2logitE3)
tau_work <- qnorm(plogis(tau_work)) - qnorm(plogis(b0logitE3))
C <- sum((tau_work^2)*pvecE3) - sum(as.matrix(tau_work*pvecE3)%*%t(as.matrix(tau_work*pvecE3)))
tauE3vec <- sqrt(1/(1+C))*tau_work
CE3 <- C
VE3 <- CE3/(1+CE3)
tau0E3 <- tauE3vec[1]
tau1E3 <- tauE3vec[2]
tau2E3 <- tauE3vec[3]
E_E3 <- sum(tauE3vec*pvecE3)
###################################################
### E4 = Paternal age
levels_E4 <- 7
pvecE4 <- c(0.342, 0.204, 0.252, 0.123, 0.052, 0.019, 0.008)
OR_E4_0 <- 1
OR_E4_1 <- 1.06
OR_E4_2 <- 1.06
OR_E4_3 <- 1.13
OR_E4_4 <- 1.22
OR_E4_5 <- 1.21
OR_E4_6 <- 1.66
lOR_E4_1 <- log(OR_E4_1)
lOR_E4_2 <- log(OR_E4_2)
lOR_E4_3 <- log(OR_E4_3)
lOR_E4_4 <- log(OR_E4_4)
lOR_E4_5 <- log(OR_E4_5)
lOR_E4_6 <- log(OR_E4_6)
### Penetrance function calculated using K and OR
functionb0logit_E4 <- function(p){
	r <- rep(NA, length(p))
	r[1] <- (pvecE4[1]*(1/(1 + exp(-p[1])))) + (pvecE4[2]*(1/(1 + exp(-p[1])*exp(-lOR_E4_1)))) + (pvecE4[3]*(1/(1 + exp(-p[1])*exp(-lOR_E4_2)))) + (pvecE4[4]*(1/(1 + exp(-p[1])*exp(-lOR_E4_3)))) + (pvecE4[5]*(1/(1 + exp(-p[1])*exp(-lOR_E4_4)))) + (pvecE4[6]*(1/(1 + exp(-p[1])*exp(-lOR_E4_5)))) + (pvecE4[7]*(1/(1 + exp(-p[1])*exp(-lOR_E4_6)))) - K
	r
}
p0 <- rep(0.001, 1)
outE4 <- dfsane(par=p0, fn=functionb0logit_E4, control=list(trace=F))
outE4
b0logitE4 <- outE4$par[1]
pD_E4_0 <- 1/(1 + exp(-1*b0logitE4))
pD_E4_1 <- 1/(1 + exp(-1*b0logitE4)*exp(-1*lOR_E4_1))
pD_E4_2 <- 1/(1 + exp(-1*b0logitE4)*exp(-1*lOR_E4_2))
pD_E4_3 <- 1/(1 + exp(-1*b0logitE4)*exp(-1*lOR_E4_3))
pD_E4_4 <- 1/(1 + exp(-1*b0logitE4)*exp(-1*lOR_E4_4))
pD_E4_5 <- 1/(1 + exp(-1*b0logitE4)*exp(-1*lOR_E4_5))
pD_E4_6 <- 1/(1 + exp(-1*b0logitE4)*exp(-1*lOR_E4_6))
b1logitE4 <- lOR_E4_1
b2logitE4 <- lOR_E4_2
b3logitE4 <- lOR_E4_3
b4logitE4 <- lOR_E4_4
b5logitE4 <- lOR_E4_5
b6logitE4 <- lOR_E4_6
###
tau_work <- c(b0logitE4, b0logitE4+b1logitE4, b0logitE4+b2logitE4, b0logitE4+b3logitE4, b0logitE4+b4logitE4,b0logitE4+b5logitE4, b0logitE4+b6logitE4)
tau_work <- qnorm(plogis(tau_work)) - qnorm(plogis(b0logitE4))
C <- sum((tau_work^2)*pvecE4) - sum(as.matrix(tau_work*pvecE4)%*%t(as.matrix(tau_work*pvecE4)))
tauE4vec <- sqrt(1/(1+C))*tau_work
CE4 <- C
VE4 <- CE4/(1+CE4)
tau0E4 <- tauE4vec[1]
tau1E4 <- tauE4vec[2]
tau2E4 <- tauE4vec[3]
tau3E4 <- tauE4vec[4]
tau4E4 <- tauE4vec[5]
tau5E4 <- tauE4vec[6]
tau6E4 <- tauE4vec[7]
E_E4 <- sum(tauE4vec*pvecE4)
###################################################
### E5 = Childhood adversity
levels_E5 <- 2
pvecE5 <- c(0.73, 1-0.73)
OR_E5_0 <- 1
OR_E5_1 <- 2.78
lOR_E5_1 <- log(OR_E5_1)
### Penetrance function calculated using K and OR
functionb0logit_E5 <- function(p){
	r <- rep(NA, length(p))
	r[1] <- (pvecE5[1]*(1/(1 + exp(-p[1])))) + (pvecE5[2]*(1/(1 + exp(-p[1])*exp(-lOR_E5_1)))) - K
	r
}
p0 <- rep(0.001, 1)
outE5 <- dfsane(par=p0, fn=functionb0logit_E5, control=list(trace=F))
outE5
b0logitE5 <- outE5$par[1]
pD_E5_0 <- 1/(1 + exp(-1*b0logitE5))
pD_E5_1 <- 1/(1 + exp(-1*b0logitE5)*exp(-1*lOR_E5_1))
b1logitE5 <- lOR_E5_1
b0probitE5 <- qnorm(plogis(b0logitE5))
b1probitE5 <- qnorm(plogis(b0logitE5 + b1logitE5)) - b0probitE5
tau0E5 <- 0
tau1E5 <- b1probitE5/sqrt(1 + ((b1probitE5^2)*pvecE5[2]*(1-pvecE5[2])))
tauE5vec <- c(tau0E5, tau1E5)
VE5 <- ((tau1E5^2)*pvecE5[2]*(1-pvecE5[2]))
E_E5 <- sum(tauE5vec*pvecE5)
###################################################
### Extra parameters relating to the ERS:
# Mean of the ERS (on the liaility scale)
E_ERS <- E_E1 + E_E2 + E_E3 + E_E4 + E_E5
# Variance of the ERS (on the liaility scale)
V_ERS <- VE1 + VE2 + VE3 + VE4 + VE5
# Disease threshold approximation
D_threshold <- qnorm(1-K) + E_ERS
# All possible combinations of the environmental risk factors
combos_E <- expand.grid(1:levels_E1, 1:levels_E2, 1:levels_E3, 1:levels_E4, 1:levels_E5)
# ers_vec: vector containing all possible (unique) ERS values
# prob_ers_vec: vector containing the probability for each entry in ers_vec
ers_vec <- rep(0, dim(combos_E)[1])
prob_ers_vec <- rep(0, dim(combos_E)[1])
for(i in 1:dim(combos_E)[1]){
	EI_i <- unlist(combos_E[i,])
	prob_ers_vec[i] <- pvecE1[EI_i[1]]*pvecE2[EI_i[2]]*pvecE3[EI_i[3]]*pvecE4[EI_i[4]]*pvecE5[EI_i[5]]
	ers_vec[i] <- tauE1vec[EI_i[1]] + tauE2vec[EI_i[2]] + tauE3vec[EI_i[3]] + tauE4vec[EI_i[4]] + tauE5vec[EI_i[5]]
}
range(ers_vec)
length(unique(ers_vec)) ### 216 versus 252 combinations
### multiple combinations of the environmental risk variables
# can therefore give the same ERS value: we want unique values
# only in the ERS vector- combine duplicates and update 
# ERS probability to reflect this:
ers_vec_old <- ers_vec
ers_vec <- unique(ers_vec)
prob_ers_vec_old <- prob_ers_vec
prob_ers_vec <- rep(0, length(ers_vec))
for(i in 1:length(ers_vec)){
	ers_i <- ers_vec[i]
	match_i <- rep(0, length(ers_vec_old))
	for(j in 1:length(ers_vec_old)){
		match_i[j] <- as.numeric(ers_vec_old[j] == ers_i)
	}
	probs_ij <- prob_ers_vec_old[match_i == 1]
	prob_ers_vec[i] <- sum(probs_ij)
}
### Obtain ERS descriptive statistics
# Median ERS:
med_ers <- quantile(ers_vec)[3]
med_ers ## 0.7156674
# Modal ERS:
max(prob_ers_vec) ## 0.1076534
modeERS <- ers_vec[prob_ers_vec == max(prob_ers_vec)] ## 0.2737725
modeERS
###################################################
### Range of ERS and PRS to use in risk estimation:
# ERS values to use:
ersI <- c(ers_vec, 0.72, 0.27)
# PRS values to use:
# Standardised/ normalised PRS:
prsI <- c(0, -1.96, 1.96, -1.64, 1.64, 0.67, -0.67, seq(from=-3.5, to = 3, length.out=1000))
# Corresponding  PRS values on the liability scale:
prsI_L <- sqrt(VML)*prsI
#############################################################
### 3.1.3. Risk estimation; all LTMs ########################
#############################################################
### Joint risk: PRS, ERS
###################################################
### Schizophrenia risk generation for the above 
# defined ranges of ERS and PRS
### Blanked out as data can be downloaded form GitHub
### Code provided so that user can generate their own
# dataset if required
### Please update FILEPATH in code!!!
###################################################
# # names.out <- c("ERS", "PRS", "specialE", "E1", "E2", "E3", "E4", "E5", "risk")
# write.table(t(matrix(names.out)), file="FILEPATH/scz_risk_ers_prs.txt", row.names=F, col.names=F)
# ### NOTE specialE == {Min ERS, Median ERS, Max ERS, Modal ERS} = {1, 2, 3, 4}
# for(i in 1:length(ersI)){
	# for(j in 1:length(prsI)){
		# risk_ij <- risk_function(ers=ersI[i], prs=prsI_L[j], E_ERS= E_ERS, V_ERS= V_ERS, V_PRS=VML, K=K)
		# out_ij <- rep(0, length(names.out))
		# out_ij[1] <- ersI[i]
		# out_ij[2] <- prsI[j]
		# if(out_ij[1] == 0){
			# out_ij[3] <- 1
		# }
		# if(out_ij[1] == 0.72){
			# out_ij[3] <- 2
		# }
		# if(out_ij[1] == max(ersI)){
			# out_ij[3] <- 3
		# }
		# if(out_ij[1] == 0.27){
			# out_ij[3] <- 4
		# }
		# if(i <= dim(combos_E)[1]){
			# out_ij[4:8] <- unlist(combos_E[i,])
		# }else{
			# out_ij[4:8] <- NA
		# }
		# out_ij[9] <- risk_ij
		# out_ij <- t(matrix(out_ij))
		# write.table(out_ij, file="FILEPATH/scz_risk_ers_prs.txt", row.names=F, col.names=F, append=T, quote=F)
		# rm(risk_ij)
	# }
# }
###################################################
### READ IN GENERATED SCHIZOPHRENIA RISK DATA
###################################################
### Please update FILEPATH!!!
scz <- read.table(file="FILEPATH/scz_risk_ers_prs.txt", header=T)
dim(scz) # [1] 219526      9
scz <- data.frame(scz)
###################################################
### Single risk: PRS
###################################################
risk_prs <- risk_single_PRS_function(prs=sqrt(VML)*unique(scz$PRS), E_ERS=E_ERS, E_PRS = 0, V_PRS=VML, K=K)
scz_noERS <- cbind(unique(scz$PRS), risk_prs)
colnames(scz_noERS) <- c("PRS", "risk")
scz_noERS <- data.frame(scz_noERS)
scz_noERS <- scz_noERS[order(scz_noERS$PRS), ]
###################################################
### Single risk: ERS
###################################################
risk_ers <- risk_single_ERS_function(ers=unique(scz$ERS), E_ERS=E_ERS, E_PRS = 0, V_ERS= V_ERS, K=K)
scz_noPRS <- cbind(unique(scz$ERS), risk_ers)
colnames(scz_noPRS) <- c("ERS", "risk")
scz_noPRS <- data.frame(scz_noPRS)
scz_noPRS <- scz_noPRS[order(scz_noPRS$ERS), ]
#############################################################
### 3.1.4. FIGURE 3: Schizophrenia risk; all LTMs ###########
#############################################################
### NOTE: axis labels added externally (in Word)
###################################################
### Create Figure 2(a)
###################################################
scz1 <- scz[scz$specialE == 1, ]
scz2 <- scz[scz$specialE == 2, ]
scz3 <- scz[scz$specialE == 3, ]
scz4 <- scz[scz$specialE == 4, ]
scz1 <- scz1[order(scz1$risk), ]
scz2 <- scz2[order(scz2$risk), ]
scz3 <- scz3[order(scz3$risk), ]
scz4 <- scz4[order(scz4$risk), ]
### Group = ERS, PRS ranging betwen +/- 3:
nf <- layout(matrix(c(1,2),2,1,byrow = TRUE), c(4,3), c(3,1), TRUE)
layout.show(nf)
par(mar = c(3,3,1,1))
plot(x=scz_noERS$PRS, y=scz_noERS$risk, type="l", ylim = c(0, 0.3), xlab = "Polygenic risk score", ylab="Risk of schizophrenia", bty="n", axes=F, lwd=3, xlim=c(-3, 3), col="dark grey")
points(x=scz1$PRS, y=scz1$risk, type="l", lwd=2)
points(x=scz2$PRS, y=scz2$risk, type="l", lty=2, lwd=2)
points(x=scz3$PRS, y=scz3$risk, type="l", lty=3, lwd=2)
points(x=scz4$PRS, y=scz4$risk, type="l", lty=4, lwd=2)
leg.txt <- c("No ERS", "Maximum (ERS = 1.53)", "Median (ERS = 0.72)", "Mode (ERS = 0.27)","Minimum (ERS = 0)")
legend("topleft", legend=leg.txt, bty="n", lty=c(1, 3,2,4,1), col = c("dark grey", rep("black", 4)), lwd=c(3,rep(2,4)), text.col = c("dark grey", rep("black", 4)))
box(bty="l")
axis(2)
par(mar = c(3,3,1,1))
curve(dnorm(x, sd=1), from=-3, to=3, axes=F)
box(bty="l")
axis(1)
###################################################
### Create Figure 2(b)
###################################################
PRS1 <- scz[scz$PRS == 0, ]
### 95% CI
PRS2 <- scz[scz$PRS == 1.96, ]
PRS1$High95 <- PRS2$risk
PRS3 <- scz[scz$PRS == -1.96, ]
PRS1$Low95 <- PRS3$risk
PRS1 <- PRS1[order(PRS1$ERS), ]
nf <- layout(matrix(c(1,2),2,1,byrow = TRUE), c(4,3), c(3,1), TRUE)
pos.gray <- gray(seq(0.1,0.9,length=10))
layout.show(nf)
par(mar = c(3,3,1,1))
plot(x = scz_noPRS$ERS, y= scz_noPRS$risk, type="l", lwd = 3, ylim=c(0, 0.3), ylab="Risk of schizophrenia", bty="n", axes=F, col="dark grey")
points(x = PRS1$ERS, y=PRS1$risk, type="l", lwd = 2)
points(x = PRS1$ERS, y=PRS1$Low95, type="l", lwd = 2, lty=2, col=pos.gray[2])
points(x = PRS1$ERS, y=PRS1$High95, type="l", lwd = 2, lty=2, col=pos.gray[2])
leg.txt <- c("No PRS", "Mean PRS", "PRS = +/- 1.96")
legend("topleft", legend=leg.txt, lwd=c(3,2,2), lty=c(1,1,2), col = c("dark grey", "black", pos.gray[2]), text.col=c("dark grey", "black", pos.gray[2]), bty="n")
box(bty="l")
axis(2)
par(mar = c(3,3,1,1))
plot(ers_vec[order(ers_vec)], prob_ers_vec[order(ers_vec)], type="h", lwd=2, axes=F)
box(bty="l")
axis(1)
axis(2, labels = c("0.00", "0.11"), at=range(prob_ers_vec))
#############################################################
### 3.1.5. Specific risk examples from text in paper ########
#############################################################
### PRS only LTM:
# 95% of individuals have estimated risk between:
PRS_risk95L <- risk_single_PRS_function(prs = -1.96*sqrt(VML), E_ERS = E_ERS, V_PRS= VML, K=K)
PRS_risk95U <- risk_single_PRS_function(prs = 1.96*sqrt(VML), E_ERS = E_ERS, V_PRS= VML, K=K)
PRS_risk95L ### 0.0016
PRS_risk95U ### 0.0304
### an individual with PRS > 95th percentile risk and relative risk
risk_single_PRS_function(prs = qnorm(0.95)*sqrt(VML), E_ERS = E_ERS, V_PRS= VML, K=K) ### 0.025
risk_single_PRS_function(prs = qnorm(0.95)*sqrt(VML), E_ERS = E_ERS, V_PRS= VML, K=K)/K
### Risk PRS == 0
risk_single_PRS_function(prs = 0*sqrt(VML), E_ERS = E_ERS, V_PRS= VML, K=K) ### 0.0079
### ERS only LTM:
# Minimum risk (will be for minimum ERS here due to reference
# category definitions)
minERS_risk <- risk_single_ERS_function(ers = 0, E_ERS = E_ERS, V_ERS=V_ERS, K=K)
minERS_risk
# minERS_risk/K
# K/minERS_risk
# Maximum risk (will be for maximum ERS here due to reference
# category definitions)
maxERS_risk <- risk_single_ERS_function(ers = max(ers_vec), E_ERS = E_ERS, V_ERS=V_ERS, K=K)
maxERS_risk
maxERS_risk/K
# 95% of individuals have estimated risk between [X,Y]:
# Because ERS doesn't follow a standard distribution we opt
# to calculate intervals using bootstrap approach as follows...
# [NOTE: we also provide code to calculate the bootstrap 95% CI
# for risks from all LTMs (ERS only, PRS only and joint)]
###################################################
### A. GENERATE POPULATION SAMPLE:
###################################################
### This is to approximate 95% CI and other summary measures
# for joint risk estimates in the population under the 
# liability model
### N = size of population to generate
N = 1500000
#E_ers <- sum(ers_vec*prob_ers_vec)
# PRS value (on the liability scale) for simulated population
prs_sample <- rnorm(N, sd=sqrt(VML))
# Set up vectors to contain ERS values, joint risk estimates,
# PRS-only risk estimates and ERS-only risk estimates:
ers_sample <- rep(0, N)
risk_sample <- rep(0, N)
ers_risk_sample <- rep(0, N)
prs_risk_sample <- rep(0, N)
# Set up vector to contain joint risk estimates when PRS = 0, ERS varying:
risk_prs0 <- rep(0, N)
# Set up vector to contain joint risk estimates when PRS = 1.96, ERS varying:
risk_prsU95 <- rep(0, N)
# Set up vector to contain joint risk estimates when PRS = -1.96, ERS varying:
risk_prsL95 <- rep(0, N)
# Set up vector to contain joint risk estimates when PRS varying, ERS = min = 0:
risk_ers0 <- rep(0, N)
# Set up vector to contain joint risk estimates when PRS varying, ERS = median = 0.27:
risk_ers0.27 <- rep(0, N)
# Set up vector to contain joint risk estimates when PRS varying, ERS = max = 1.53:
risk_ers1.53 <- rep(0, N)
# Set up vector to contain joint risk estimates when PRS varying, ERS = mean = 0.3884:
risk_mean_ers <- rep(0, N)
### Loop to generate risk ERS values and risk in sample:
for(i in 1:N){
	sample_i <- unique((1:length(ers_vec))*rmultinom(1,1,prob=prob_ers_vec))
	sample_i <- sample_i[sample_i>0]
	ers_sample[i] <- ers_vec[sample_i]
	risk_sample[i] <- risk_function(ers=ers_sample[i], prs=prs_sample[i],E_ERS= E_ERS, V_ERS= V_ERS, V_PRS=VML, K=K)
	ers_risk_sample[i] <- risk_single_ERS_function(ers = ers_sample[i], E_ERS = E_ERS, V_ERS=V_ERS, K=K)
	prs_risk_sample[i] <- risk_single_PRS_function(prs = prs_sample[i], E_ERS = E_ERS, V_PRS=VML, K=K)
	risk_prs0[i] <- risk_function(ers=ers_sample[i], prs=0,E_ERS= E_ERS, V_ERS= V_ERS, V_PRS=VML, K=K)
	risk_prsU95[i] <- risk_function(ers=ers_sample[i], prs=1.96*sqrt(VML),E_ERS= E_ERS, V_ERS= V_ERS, V_PRS=VML, K=K)
	risk_prsL95[i] <- risk_function(ers=ers_sample[i], prs=-1.96*sqrt(VML),E_ERS= E_ERS, V_ERS= V_ERS, V_PRS=VML, K=K)
	risk_ers0[i] <- risk_function(ers=0, prs=prs_sample[i],E_ERS= E_ERS, V_ERS= V_ERS, V_PRS=VML, K=K)
	risk_ers0.27[i] <- risk_function(ers=0.27, prs=prs_sample[i],E_ERS= E_ERS, V_ERS= V_ERS, V_PRS=VML, K=K)
	risk_ers1.53[i] <- risk_function(ers=max(ers_vec), prs=prs_sample[i],E_ERS= E_ERS, V_ERS= V_ERS, V_PRS=VML, K=K)
	risk_mean_ers[i] <- risk_function(ers=E_ERS, prs=prs_sample[i],E_ERS= E_ERS, V_ERS= V_ERS, V_PRS=VML, K=K)
}
pop_sample <- cbind(prs_sample, ers_sample, risk_sample, ers_risk_sample, prs_risk_sample, risk_prs0, risk_prsL95, risk_prsU95, risk_ers0, risk_ers0.27, risk_ers1.53, risk_mean_ers)
pop_sample <- data.frame(pop_sample)
colnames(pop_sample)[1:5] <- c("prs", "ers", "risk", "ers_risk", "prs_risk")
###################################################
### B. Approx 95% of sampled individuals have a risk
# between...
###################################################
### PRS:
quantile(pop_sample$prs_risk, probs=c(0.025, 0.975))
### ERS:
#quantile(pop_sample$ers_risk, probs=c(0.025, 0.975)) 
### Use below due to ERS distribution properties:
quantile(pop_sample$ers_risk, probs=c(0.0, 0.95)) 
### PRS and ERS:
quantile(pop_sample$risk, probs=c(0.025, 0.975))
###################################################
### PRS and ERS LTM:
# Minimum and maximum risk when PRS == 0
risk_function(ers=0, prs=0, E_ERS= E_ERS, V_ERS= V_ERS, V_PRS=VML, K=K) ### 0.0017
risk_function(ers=max(ers_vec), prs=0, E_ERS= E_ERS, V_ERS= V_ERS, V_PRS=VML, K=K) ### 0.1018
# See B. above for approximate 95% CI for risk
#############################################################
### 3.2. Logistic ERS only model ############################
#############################################################
### 3.2.1. General code and functions #######################
#############################################################
### Finding intercept (b0) for univariate logistic models 
# Uses the penetrance function calculated using K and OR as
# inputs
### Function for variable E1 (Cannabis usage)
functionb0logit_E1 <- function(p){
	r <- rep(NA, length(p))
	r[1] <- (pvecE1[1]*(1/(1 + exp(-p[1])))) + (pvecE1[2]*(1/(1 + exp(-p[1])*exp(-lOR_E1_1)))) + (pvecE1[3]*(1/(1 + exp(-p[1])*exp(-lOR_E1_2)))) - K
	r
}
### Function for variable E2 (Migration)
functionb0logit_E2 <- function(p){
	r <- rep(NA, length(p))
	r[1] <- (pvecE2[1]*(1/(1 + exp(-p[1])))) + (pvecE2[2]*(1/(1 + exp(-p[1])*exp(-lOR_E2_1)))) - K
	r
}
### Function for variable E3 (Urbanicity)
functionb0logit_E3 <- function(p){
	r <- rep(NA, length(p))
	r[1] <- (pvecE3[1]*(1/(1 + exp(-p[1])))) + (pvecE3[2]*(1/(1 + exp(-p[1])*exp(-lOR_E3_1)))) + (pvecE3[3]*(1/(1 + exp(-p[1])*exp(-lOR_E3_2)))) - K
	r
}
### Function for variable E4 (Paternal age)
functionb0logit_E4 <- function(p){
	r <- rep(NA, length(p))
	r[1] <- (pvecE4[1]*(1/(1 + exp(-p[1])))) + (pvecE4[2]*(1/(1 + exp(-p[1])*exp(-lOR_E4_1)))) + (pvecE4[3]*(1/(1 + exp(-p[1])*exp(-lOR_E4_2)))) + (pvecE4[4]*(1/(1 + exp(-p[1])*exp(-lOR_E4_3)))) + (pvecE4[5]*(1/(1 + exp(-p[1])*exp(-lOR_E4_4)))) + (pvecE4[6]*(1/(1 + exp(-p[1])*exp(-lOR_E4_5)))) + (pvecE4[7]*(1/(1 + exp(-p[1])*exp(-lOR_E4_6)))) - K
	r
}
### Function for variable E5 (Childhood adversity)
functionb0logit_E5 <- function(p){
	r <- rep(NA, length(p))
	r[1] <- (pvecE5[1]*(1/(1 + exp(-p[1])))) + (pvecE5[2]*(1/(1 + exp(-p[1])*exp(-lOR_E5_1)))) - K
	r
}
#############################################################
### 3.2.2. Parameterising model #############################
#############################################################
# To parameterise the logistic ERS only model we need to 
# fully parameterise each univariate logistic model for 
# variables E1 - E5 (specifically find each intercept). We 
# then use numerical optimization with the univariate model 
# parameters as inputs (this is function joint_logit_function 
# given below in subheading 'Parameterising ERS only logistic 
# model').
###################################################
### Variable E1 (Cannabis usage)
###################################################
# Find the intercept using functionb0logit_E1
p0 <- rep(0.001, 1)
outE1 <- dfsane(par=p0, fn=functionb0logit_E1, control=list(trace=F))
outE1
b0logitE1 <- outE1$par[1]
# E1 penetrance values for the logistic model
pD_E1_0 <- 1/(1 + exp(-1*b0logitE1))
pD_E1_1 <- 1/(1 + exp(-1*b0logitE1)*exp(-1*lOR_E1_1))
pD_E1_2 <- 1/(1 + exp(-1*b0logitE1)*exp(-1*lOR_E1_2))
# Additional univariate beta estimates
b1logitE1 <- lOR_E1_1
b2logitE1 <- lOR_E1_2
###################################################
### Variable E2 (Migration)
###################################################
# Find the intercept using functionb0logit_E2
p0 <- rep(0.001, 1)
outE2 <- dfsane(par=p0, fn=functionb0logit_E2, control=list(trace=F))
b0logitE2 <- outE2$par[1]
# E2 penetrance values for the logistic model
pD_E2_0 <- 1/(1 + exp(-1*b0logitE2))
pD_E2_1 <- 1/(1 + exp(-1*b0logitE2)*exp(-1*lOR_E2_1))
# Additional univariate beta estimates
b1logitE2 <- lOR_E2_1
###################################################
### Variable E3 (Urbanicity)
###################################################
# Find the intercept using functionb0logit_E3
p0 <- rep(0.001, 1)
outE3 <- dfsane(par=p0, fn=functionb0logit_E3, control=list(trace=F))
b0logitE3 <- outE3$par[1]
# E3 penetrance values for the logistic model
pD_E3_0 <- 1/(1 + exp(-1*b0logitE3))
pD_E3_1 <- 1/(1 + exp(-1*b0logitE3)*exp(-1*lOR_E3_1))
pD_E3_2 <- 1/(1 + exp(-1*b0logitE3)*exp(-1*lOR_E3_2))
# Additional univariate beta estimates
b1logitE3 <- lOR_E3_1
b2logitE3 <- lOR_E3_2
###################################################
### Variable E4 (Paternal age)
###################################################
# Find the intercept using functionb0logit_E4
p0 <- rep(0.001, 1)
outE4 <- dfsane(par=p0, fn=functionb0logit_E4, control=list(trace=F))
b0logitE4 <- outE4$par[1]
# E4 penetrance values for the logistic model
pD_E4_0 <- 1/(1 + exp(-1*b0logitE4))
pD_E4_1 <- 1/(1 + exp(-1*b0logitE4)*exp(-1*lOR_E4_1))
pD_E4_2 <- 1/(1 + exp(-1*b0logitE4)*exp(-1*lOR_E4_2))
pD_E4_3 <- 1/(1 + exp(-1*b0logitE4)*exp(-1*lOR_E4_3))
pD_E4_4 <- 1/(1 + exp(-1*b0logitE4)*exp(-1*lOR_E4_4))
pD_E4_5 <- 1/(1 + exp(-1*b0logitE4)*exp(-1*lOR_E4_5))
pD_E4_6 <- 1/(1 + exp(-1*b0logitE4)*exp(-1*lOR_E4_6))
# Additional univariate beta estimates
b1logitE4 <- lOR_E4_1
b2logitE4 <- lOR_E4_2
b3logitE4 <- lOR_E4_3
b4logitE4 <- lOR_E4_4
b5logitE4 <- lOR_E4_5
b6logitE4 <- lOR_E4_6
###################################################
### Variable E5 (Childhood adversity)
###################################################
# Find the intercept using functionb0logit_E5
p0 <- rep(0.001, 1)
outE5 <- dfsane(par=p0, fn=functionb0logit_E5, control=list(trace=F))
b0logitE5 <- outE5$par[1]
# E5 penetrance values for the logistic model
pD_E5_0 <- 1/(1 + exp(-1*b0logitE5))
pD_E5_1 <- 1/(1 + exp(-1*b0logitE5)*exp(-1*lOR_E5_1))
# Additional univariate beta estimates
b1logitE5 <- lOR_E5_1
###################################################
### Parameterising ERS only logistic model
###################################################
# Create a design matrix (called combos_p) which would
# be multipled by ERS only model parameters to calculate
# ERS
combos_p <- matrix(0, nrow=dim(combos_E)[1], ncol=(1 + (levels_E1-1) + (levels_E2-1) + (levels_E3-1) + (levels_E4-1) + (levels_E5-1)))
combos_p[,1] <- 1
for(i in 1:dim(combos_E)[1]){
	out.work <- combos_E[i,]
	if(out.work[1] == 2){combos_p[i,2] <- 1}
	if(out.work[1] == 3){combos_p[i,3] <- 1}
	if(out.work[2] == 2){combos_p[i,4] <- 1}
	if(out.work[3] == 2){combos_p[i,5] <- 1}
	if(out.work[3] == 3){combos_p[i,6] <- 1}
	if(out.work[4] == 2){combos_p[i,7] <- 1}
	if(out.work[4] == 3){combos_p[i,8] <- 1}
	if(out.work[4] == 4){combos_p[i,9] <- 1}
	if(out.work[4] == 5){combos_p[i,10] <- 1}
	if(out.work[4] == 6){combos_p[i,11] <- 1}
	if(out.work[4] == 7){combos_p[i,12] <- 1}
	if(out.work[5] == 2){combos_p[i,13] <- 1}
}
### Probability vector for risk factor combinations in 
# combos_p
prob_ers_vec_logistic_old <- rep(0, dim(combos_E)[1])
for(i in 1:dim(combos_E)[1]){
	EI_i <- unlist(combos_E[i,])
	prob_ers_vec_logistic_old[i] <- pvecE1[EI_i[1]]*pvecE2[EI_i[2]]*pvecE3[EI_i[3]]*pvecE4[EI_i[4]]*pvecE5[EI_i[5]]
}
### Function to use in numerical optimization to find ERS only 
### logistic model parameters
# Each model parameter requires a separate line of code...
joint_logit_function <- function(p){
	r <- rep(NA, length(p))
	r[1] <- t(matrix(prob_ers_vec_logistic_old))%*%plogis(combos_p%*%matrix(p)) - K
	## E1 Level 2 and 3
	r[2] <- t(matrix(prob_ers_vec_logistic_old[combos_E[,1] == 2]/pvecE1[2]))%*%plogis(combos_p[combos_E[,1] == 2,]%*%matrix(p)) - plogis(b0logitE1 + b1logitE1)
	r[3] <- t(matrix(prob_ers_vec_logistic_old[combos_E[,1] == 3]/pvecE1[3]))%*%plogis(combos_p[combos_E[,1] == 3,]%*%matrix(p)) - plogis(b0logitE1 + b2logitE1)
	## E2 Level 2
	r[4] <- t(matrix(prob_ers_vec_logistic_old[combos_E[,2] == 2]/pvecE2[2]))%*%plogis(combos_p[combos_E[,2] == 2,]%*%matrix(p)) - plogis(b0logitE2 + b1logitE2)
	## E3 Level 2 and 3
	r[5] <- t(matrix(prob_ers_vec_logistic_old[combos_E[,3] == 2]/pvecE3[2]))%*%plogis(combos_p[combos_E[,3] == 2,]%*%matrix(p)) - plogis(b0logitE3 + b1logitE3)
	r[6] <- t(matrix(prob_ers_vec_logistic_old[combos_E[,3] == 3]/pvecE3[3]))%*%plogis(combos_p[combos_E[,3] == 3,]%*%matrix(p)) - plogis(b0logitE3 + b2logitE3)
	## E4 Level 2-7
	r[7] <- t(matrix(prob_ers_vec_logistic_old[combos_E[,4] == 2]/pvecE4[2]))%*%plogis(combos_p[combos_E[,4] == 2,]%*%matrix(p)) - plogis(b0logitE4 + b1logitE4)
	r[8] <- t(matrix(prob_ers_vec_logistic_old[combos_E[,4] == 3]/pvecE4[3]))%*%plogis(combos_p[combos_E[,4] == 3,]%*%matrix(p)) - plogis(b0logitE4 + b2logitE4)
	r[9] <- t(matrix(prob_ers_vec_logistic_old[combos_E[,4] == 4]/pvecE4[4]))%*%plogis(combos_p[combos_E[,4] == 4,]%*%matrix(p)) - plogis(b0logitE4 + b3logitE4)
	r[10] <- t(matrix(prob_ers_vec_logistic_old[combos_E[,4] == 5]/pvecE4[5]))%*%plogis(combos_p[combos_E[,4] == 5,]%*%matrix(p)) - plogis(b0logitE4 + b4logitE4)
	r[11] <- t(matrix(prob_ers_vec_logistic_old[combos_E[,4] == 6]/pvecE4[6]))%*%plogis(combos_p[combos_E[,4] == 6,]%*%matrix(p)) - plogis(b0logitE4 + b5logitE4)
	r[12] <- t(matrix(prob_ers_vec_logistic_old[combos_E[,4] == 7]/pvecE4[7]))%*%plogis(combos_p[combos_E[,4] == 7,]%*%matrix(p)) - plogis(b0logitE4 + b6logitE4)
	## E5 Level 2
	r[13] <- t(matrix(prob_ers_vec_logistic_old[combos_E[,5] == 2]/pvecE5[2]))%*%plogis(combos_p[combos_E[,5] == 2,]%*%matrix(p)) - plogis(b0logitE5 + b1logitE5)
	r
}
### Use function to obtain model parameters
p0 <- c(-4, b1logitE1, b2logitE1, b1logitE2, b1logitE3, b2logitE3, b1logitE4, b2logitE4, b3logitE4, b4logitE4, b5logitE4, b6logitE4, b1logitE5)
out_logistic_joint <- dfsane(par=p0, fn= joint_logit_function, control=list(trace=F))
out_logistic_joint
logOR_joint <- out_logistic_joint$par
### Parameters for each risk factor
jlOR_E1 <- c(0, logOR_joint[2:3])
jlOR_E2 <- c(0, logOR_joint[4])
jlOR_E3 <- c(0, logOR_joint[5:6])
jlOR_E4 <- c(0, logOR_joint[7:12])
jlOR_E5 <- c(0, logOR_joint[13])
### Variance component, E1
vE1_jli <- sum((jlOR_E1^2)*pvecE1) - sum(matrix(jlOR_E1)%*%t(matrix(jlOR_E1))*matrix(pvecE1)%*%t(matrix(pvecE1)))
### Variance component, E2
vE2_jli <- sum((jlOR_E2^2)*pvecE2) - sum(matrix(jlOR_E2)%*%t(matrix(jlOR_E2))*matrix(pvecE2)%*%t(matrix(pvecE2)))
### Variance component, E3
vE3_jli <- sum((jlOR_E3^2)*pvecE3) - sum(matrix(jlOR_E3)%*%t(matrix(jlOR_E3))*matrix(pvecE3)%*%t(matrix(pvecE3)))
### Variance component, E4
vE4_jli <- sum((jlOR_E4^2)*pvecE4) - sum(matrix(jlOR_E4)%*%t(matrix(jlOR_E4))*matrix(pvecE4)%*%t(matrix(pvecE4)))
### Variance component, E5
vE5_jli <- sum((jlOR_E5^2)*pvecE5) - sum(matrix(jlOR_E5)%*%t(matrix(jlOR_E5))*matrix(pvecE5)%*%t(matrix(pvecE5)))
### Total ERS variance
vERSi <- vE1_jli + vE2_jli + vE3_jli + vE4_jli + vE5_jli
### Proportion of variability in the logistic-liability of the ERS 
# only logistic model attributable to each risk factor
prop_v_RF <- (c(vE1_jli, vE2_jli, vE3_jli, vE4_jli, vE5_jli)/(vERSi + ((pi^2)/3)))
prop_v_RF
### Proportion of variability in the logistic-liability of the ERS 
# only logistic model attributable to the ERS
prop_vERS <- sum(c(vE1_jli, vE2_jli, vE3_jli, vE4_jli, vE5_jli)/(vERSi + ((pi^2)/3)))
prop_vERS
#############################################################
### 3.2.3. Risk estimation (ERS only logistic model) ########
#############################################################
### Find ERS for every risk factor combination under the ERS
# only logistic model
score_logistic_joint <- combos_p%*%matrix(logOR_joint)
logistic_risk_joint <- plogis(combos_p%*%matrix(logOR_joint))
#############################################################
### 3.2.4. FIGURE 4 #########################################
#############################################################
### ERS only logistic model expected risk versus ERS only 
# liability model expected risk
### Need vector of expected risks from the ERS only LTM for
# all risk factor combinations
liability_risk2 <- rep(0, length(ers_vec_old))
for(i in 1:length(ers_vec_old)){
	liability_risk2[i] <- risk_single_ERS_function(ers = ers_vec_old[i], E_ERS = E_ERS, V_ERS=V_ERS, K=K)
}
### Fig. 4
plot(x=liability_risk2, y=logistic_risk_joint, xlab="Liability risk", ylab="Logistic risk")
abline(a=0, b=1)
#############################################################
#############################################################
### 3.3. Additional possible approximate logistic models ####
#############################################################
### Presented in the appendix of the paper ##################
#############################################################
### Logistic risk functions:
### Joint risk function: ERS and PRS in model
logis_risk_function <- function(ers, prs, E_ERS, V_ERS, E_PRS=0, V_PRS, K){
	E_L <- E_PRS + E_ERS
	L_T <- E_L - qlogis(K, s = sqrt(3)/pi)
	denom <- sqrt(1 - V_ERS - V_PRS)
	risk <- plogis(-(L_T - ers - prs)/denom,  s = sqrt(3)/pi)
	risk
}
### Risk function: only ERS:
logis_risk_single_ERS_function <- function(ers, E_ERS, E_PRS = 0, V_ERS, K){
	E_L <- E_PRS + E_ERS
	L_T <- E_L - qlogis(K, s = sqrt(3)/pi)
	denom <- sqrt(1 - V_ERS)
	risk <- plogis(-(L_T - ers)/denom,  s = sqrt(3)/pi)
	risk
}
### Risk function: only PRS:
logis_risk_single_PRS_function <- function(prs, E_ERS, E_PRS = 0, V_PRS, K){
	E_L <- E_PRS + E_ERS
	L_T <- - qlogis(K, s = sqrt(3)/pi)
	denom <- sqrt(1 - V_PRS)
	risk <- plogis(-(L_T - prs)/denom,  s = sqrt(3)/pi)
	risk
}
#########################################################
### 3.3.1. Possible approximate PRS only logistic model #
#########################################################
### To find this approximate model we first need the ORs from the
# PGC Schizophrenia working groups PRS SNPs. 
### We also need the corresponding MAF. We are working with OR 
# estimates from a European population and therefore require
# MAFs from a European population- we use the 1000G European
# population.
### Both datasets are provided onthe GitHub
### scz_prs_snps.txt == PGC schizophrenia SNP info
### SCZ_1000G_AF.txt == 1000G European MAFs
### CHANGE FILEPATH in below code
###################################################
### Read in PGC data
PGCdata <- fread(file="/FILEPATH/scz_prs_snps.txt")
### Select final PRS SNPs
PRSdata <- PGCdata[p <= 0.05,]
dim(PRSdata)
length(unique(PRSdata[,snpid]))
### Isolate rs names for future use
PRS_rs <- data.frame(PRSdata[,1])
PRS_rs <- as.character(unlist(PRS_rs))
colnames(PRS_rs) <- NULL
### Read in 1000G MAF data
scz_maf <- read.table(file="/FILEPATH/SCZ_1000G_AF.txt")
scz_maf$SNP1 <- as.character(unlist(scz_maf$SNP))
scz_maf <- data.table(scz_maf)
### Match SNPs between 2 datasets
setkey(scz_maf, SNP1)
PRSdata[, SNP1:= as.character(unlist(snpid))]
RRdata1 <- merge(scz_maf, PRSdata, by="SNP1")
dim(RRdata1)
### Create risk allele frequencies (RAF) for PRS SNPs
# and find ORs corresponding to these RAF (ORprs)
RAF <- rep(0, dim(RRdata1)[1])
ORprs <- RRdata1$or
for(i in 1:dim(RRdata1)[1]){
	mA1i <- as.character(unlist(RRdata1$A1[i]))
	if(RRdata1$or[i] < 1){
		ORprs[i] <- 1/RRdata1$or[i]
		refA1i <- as.character(unlist(RRdata1$a2[i]))
		if(mA1i == refA1i){
			RAF[i] <- 1 - RRdata1$MAF[i]
		}else{
			RAF[i] <- RRdata1$MAF[i]
		}
	}else{
		refA1i <- as.character(unlist(RRdata1$a1[i]))
		if(mA1i == refA1i){
			RAF[i] <- 1 - RRdata1$MAF[i]
		}else{
			RAF[i] <- RRdata1$MAF[i]
		}
	}
}
### Rename RAF due to some previous code I didn't want to change!
fcontrol_prs <- RAF
### Obtain the effect size estimates (betas) from the ORs
beta_logit <- log(ORprs)
### Variance components
# SNP variances
varG <- 2*fcontrol_prs*(1-fcontrol_prs)
# PRS variance component using original log ORs:
vPRS_logit <- sum((beta_logit^2)*varG) ### 10.68
### Work to tranform log ORs for the PRS to the liability scale
# and so calculate the variability on the liability scale from
# these transformed effect size estimates...
beta_logit_joint <- beta_logit/(sqrt(1 + ((beta_logit^2)*varG)))
beta_logit_var_std <- beta_logit/(sqrt((pi^2)/3 + ((beta_logit^2)*varG)))
sum((beta_logit_var_std^2)*2*fcontrol_prs*(1-fcontrol_prs))
# Approximating joint intercept...
b0work <- log(K/(1-K))
# Calculating the corresponding LTM effect size estimates
tauPRS_work <- qnorm(plogis(b0work + beta_logit)) - qnorm(K)
# The variance component on the liability scale using these 
# LTM effect size estimates
VML_work <- sum((tauPRS_work^2)*varG) ### 1.48
### The ratio of the PGC present PRS heritability (VML == 0.07)
# and the estimated variance component using the transformed 
# log OR:
C_var_work <- VML/VML_work
### We use this ratio C_var_work to approximate the proportion 
# of variability on the logistic-liability scale attributable to
# the PRS when all PRS SNPs are included simultaneously in the 
# logistc model
VM_logistic <- C_var_work*sum((beta_logit_var_std^2)*2*fcontrol_prs*(1-fcontrol_prs)) ### 0.15
############################################################################
### PRS only Risk estimates ################################################
############################################################################
### Expected value of joint ERS values
# Needed because even though ERS is unobserved for PRS only model, it is
# collasped into the mean liability in this example...
E_ERS_all2 <- sum((jlOR_E1/sqrt((vERSi + ((pi^2)/3))))*pvecE1) + 
sum((jlOR_E2/sqrt((vERSi + ((pi^2)/3))))*pvecE2) +
sum((jlOR_E3/sqrt((vERSi + ((pi^2)/3))))*pvecE3) +
sum((jlOR_E4/sqrt((vERSi + ((pi^2)/3))))*pvecE4) +
sum((jlOR_E5/sqrt((vERSi + ((pi^2)/3))))*pvecE5)
### Setting up prs values to calculate risk for
prs_vec_std <- c(-4, -1.96, -1.64, -1.28, -0.67, -0.12, 0, 0.12, 0.67, 1.28, 1.64, 1.96, 4)
### Setting up output matrix
out_df_PRS <- matrix(0, nrow=length(prs_vec_std), ncol=3)
out_df_PRS[,1] <- prs_vec_std
colnames(out_df_PRS) <- c("PRS", "Liability", "Logistic")
### Calculating risk for prs values for both the PRS only liability model and the 
# approximate PRS only logistic model
for(i in 1:length(prs_vec_std)){
	out_df_PRS[i,2] <- risk_single_PRS_function(prs=prs_vec_std[i]*sqrt(VML), E_ERS= E_ERS, E_PRS = 0, V_PRS=VML, K=K)
	out_df_PRS[i,3] <- logis_risk_single_PRS_function(prs=prs_vec_std[i]*sqrt(VM_logistic), E_ERS= E_ERS_all2, E_PRS = 0, V_PRS= VM_logistic, K=K)
}
out_df_PRS <- data.frame(out_df_PRS)
### Quick plot to see data
plot(out_df_PRS$Liability, out_df_PRS$Logistic)
abline(a=0, b=1)
############################################################
### 3.3.2. Possible approximate joint PRS and ERS logistic #
########## model ###########################################
############################################################
### Re-standardise previously calculated joint ERS effect size 
# estimates
updated_ERS_params2 <- logOR_joint/sqrt((vERSi + ((pi^2)/3)))
### Set intercept here equal to zero as it will now be captured
# the mean
updated_ERS_params2[1] <- 0
### Calculate updated ERS
ers_values2 <- combos_p%*%matrix(updated_ERS_params2)
ers_values2 <- as.vector(ers_values2)
### Re-name the approximated logistic joint PRS variance due to
# another name being used in following code
v_PRS_all2 <- VM_logistic
### Re-name the estimated logistic joint ERS variance due to
# another name being used in following code
v_ERS_all2 <- prop_vERS
### Create space for risk output
out.joint <- NULL
### Risk calculation for liability model and approximate logistic
# model
for(i in 1:length(prs_vec_std)){
	for(j in 1:length(ers_vec_old)){
		score_liab_ij <- ers_vec_old[j] + sqrt(VML)*prs_vec_std[i]
		score_log_ij <- ers_values2[j] + sqrt(v_PRS_all2)*prs_vec_std[i]
		risk_liab_ij <- risk_function(ers=ers_vec_old[j], prs=prs_vec_std[i]*sqrt(VML), E_ERS= E_ERS, V_ERS= V_ERS, E_PRS=0, V_PRS=VML, K=K)
		risk_log_ij <- logis_risk_function(ers= ers_values2[j], prs = prs_vec_std[i]*sqrt(v_PRS_all2), E_ERS= E_ERS_all2, V_ERS= v_ERS_all2, E_PRS=0, V_PRS= v_PRS_all2, K=K)
		out_ij <- t(matrix(c(score_liab_ij, score_log_ij, risk_liab_ij, risk_log_ij)))
		out.joint <- rbind(out.joint, out_ij)
	}
}
colnames(out.joint) <- c("Liab_score", "Logit_score", "Liab_risk", "Logit_risk")
out.joint <- data.frame(out.joint)
# Quick plot to check data
plot(out.joint$Liab_risk, out.joint$Logit_risk)
abline(a = 0, b=1)
#############################################################
### 3.3.3. FIGURE A3. #######################################
#############################################################
### Estimated liability thresmodel model risks against 
# estimated logistic model risks for ERS only, PRS only and, 
# PRS and ERS models
### Combined all risk estimates into single dataset
df.all <- rbind(cbind(out.joint$Liab_risk, out.joint$Logit_risk, rep("ERS & PRS", length(out.joint$Logit_risk))),cbind(out_df_PRS$Liability, out_df_PRS$Logistic, rep("PRS only", length(out_df_PRS$Logistic))),
cbind(liability_risk2, logistic_risk_joint, rep("ERS only", length(logistic_risk_joint))))
colnames(df.all) <- c("Liability", "Logistic", "Model")
df.all <- data.frame(df.all)
df.all$Liability <- as.numeric(as.character(df.all$Liability))
df.all$Logistic <- as.numeric(as.character(df.all$Logistic))
pA3 <- ggplot(df.all, aes(x=Liability, y=Logistic)) + geom_point() + geom_abline(intercept=0, slope=1) + scale_color_grey() +theme_bw() + xlab("Joint liability model risk") + ylab("Joint logistic model risk")
pA3 + facet_grid(Model~.)