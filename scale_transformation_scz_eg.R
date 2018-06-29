###################################################
### Schizophrenia Example #########################
###################################################
### Libraries needed and source other scripts #####
###################################################
### INSTALL REQUIRED LIBRARIES:
library(BB) ## uniroot- when one parameter, or dfsane when multiple parameters
library(ggplot2)
library(gridExtra)
library(ggExtra)
### General known parameters:
K <- 0.01 ### prevalence
H2 <- 0.81 ### Broad-sense heritability
h2 <- 0.81 ### Narrow-sense heritability
VML <- 0.07 ### Var[PRS] on L-scale
VD <- H2-h2 ### (Quasi) Dominant variance component
###################################################
### Environmental variables:
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
E_ERS <- E_E1 + E_E2 + E_E3 + E_E4 + E_E5
V_ERS <- VE1 + VE2 + VE3 + VE4 + VE5
D_threshold <- qnorm(1-K) + E_ERS
pnorm(-D_threshold/(sqrt(1 - V_ERS - VML)))
###################################################
risk_function <- function(ers, prs, E_ERS, V_ERS, E_PRS=0, V_PRS, K){
	E_L <- E_PRS + E_ERS
	L_T <- E_L - qnorm(K)
	denom <- sqrt(1 - V_ERS - V_PRS)
	risk <- pnorm(-(L_T - ers - prs)/denom)
	risk
}

risk_single_ERS_function <- function(ers, E_ERS, E_PRS = 0, V_ERS, K){
	E_L <- E_PRS + E_ERS
	L_T <- E_L - qnorm(K)
	denom <- sqrt(1 - V_ERS)
	risk <- pnorm(-(L_T - ers)/denom)
	risk
}
risk_single_PRS_function <- function(prs, E_ERS, E_PRS = 0, V_PRS, K){
	E_L <- E_PRS + E_ERS
	L_T <- qnorm(1-K)
	denom <- sqrt(1 - V_PRS)
	risk <- pnorm(-(L_T - prs)/denom)
	risk
}
###################################################
combos_E <- expand.grid(1:levels_E1, 1:levels_E2, 1:levels_E3, 1:levels_E4, 1:levels_E5)
ers_vec <- rep(0, dim(combos_E)[1])
prob_ers_vec <- rep(0, dim(combos_E)[1])
for(i in 1:dim(combos_E)[1]){
	EI_i <- unlist(combos_E[i,])
	prob_ers_vec[i] <- pvecE1[EI_i[1]]*pvecE2[EI_i[2]]*pvecE3[EI_i[3]]*pvecE4[EI_i[4]]*pvecE5[EI_i[5]]
	ers_vec[i] <- tauE1vec[EI_i[1]] + tauE2vec[EI_i[2]] + tauE3vec[EI_i[3]] + tauE4vec[EI_i[4]] + tauE5vec[EI_i[5]]
}
range(ers_vec)
length(unique(ers_vec)) ### 216 versus 252 combinations
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

med_ers <- quantile(ers_vec)[3]
table(round(ers_vec, digits =2) == round(med_ers, digits = 2))
table(round(ers_vec_old, digits =2) == round(med_ers, digits = 2))
ers_vec_old[round(ers_vec_old, digits =2) == round(med_ers, digits = 2)]
sum(prob_ers_vec[round(ers_vec, digits =2) == round(med_ers, digits = 2)])
combos_E[round(ers_vec_old, digits =2) == round(med_ers, digits = 2), ]
### most frequent value:
max(prob_ers_vec) ## 0.1076534
ers_vec[prob_ers_vec == max(prob_ers_vec)] ## 0.2737725
### Round most frequent value to 0.27:
sum(prob_ers_vec[round(ers_vec, digits =2) == 0.27]) ## 0.1163041
table(round(ers_vec_old, digits =2) == 0.27)
combos_E[round(ers_vec_old, digits =2) == 0.27, ]

ersI <- c(ers_vec, 0.72, 0.27)
prsI <- c(0, -1.96, 1.96, -1.64, 1.64, 0.67, -0.67, seq(from=-3.5, to = 3, length.out=1000))
prsI_L <- sqrt(VML)*prsI
###################################################
### Single risk: ERS
###################################################
minERS_risk <- risk_single_ERS_function(ers = 0, E_ERS = E_ERS, V_ERS=V_ERS, K=K)
medERS_risk <- risk_single_ERS_function(ers = 0.72, E_ERS = E_ERS, V_ERS=V_ERS, K=K)
modeERS_risk <- risk_single_ERS_function(ers = 0.27, E_ERS = E_ERS, V_ERS=V_ERS, K=K)
maxERS_risk <- risk_single_ERS_function(ers = max(ers_vec), E_ERS = E_ERS, V_ERS=V_ERS, K=K)
plot(ers_vec, prob_ers_vec)
plot(ers_vec[order(ers_vec)], prob_ers_vec[order(ers_vec)], type="l")
plot(ers_vec[order(ers_vec)], prob_ers_vec[order(ers_vec)], type="h", lwd=2)
#plot(prob_ers_vec, ers_vec)
ers_vec_2dps <- round(ers_vec, digits=2)
prob_ers_vec2dp <- prob_ers_vec[order(ers_vec_2dps)]
ers_vec_2dps <- ers_vec_2dps[order(ers_vec_2dps)]
plot(ers_vec_2dps, prob_ers_vec2dp, type="l")

# N = 1000000
# ers_sample <- rep(0, N)
# for(i in 1:N){
	# sample_i <- unique((1:length(ers_vec))*rmultinom(1,1,prob=prob_ers_vec))
	# sample_i <- sample_i[sample_i>0]
	# ers_sample[i] <- ers_vec[sample_i]
# }
# hist(ers_sample)
###################################################
### Single risk: PRS
###################################################
minPRS_risk <- risk_single_PRS_function(prs = -4*sqrt(VML), E_ERS = E_ERS, V_PRS=VML, K=K)
medPRS_risk <- risk_single_PRS_function(prs = 0, E_ERS = E_ERS, V_PRS= VML, K=K)
maxPRS_risk <- risk_single_PRS_function(prs = 4, E_ERS = E_ERS, V_PRS= VML, K=K)
PRS_risk95U <- risk_single_PRS_function(prs = 1.96*sqrt(VML), E_ERS = E_ERS, V_PRS= VML, K=K)
PRS_risk95L <- risk_single_PRS_function(prs = -1.96*sqrt(VML), E_ERS = E_ERS, V_PRS= VML, K=K)
PRS_risk99U <- risk_single_PRS_function(prs = qnorm(0.99)*sqrt(VML), E_ERS = E_ERS, V_PRS= VML, K=K)
PRS_risk99L <- risk_single_PRS_function(prs = -1*qnorm(0.99)*sqrt(VML), E_ERS = E_ERS, V_PRS= VML, K=K)
###################################################
### Joint risk: PRS, ERS
###################################################
# # names.out <- c("ERS", "PRS", "specialE", "E1", "E2", "E3", "E4", "E5", "risk")
# write.table(t(matrix(names.out)), file="/Users/alexg/Documents/PAPERS_THESIS/Transforming_summ_stats/R_code/scz_risk_ers_prs20180521.txt", row.names=F, col.names=F)

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
		# write.table(out_ij, file="/Users/alexg/Documents/PAPERS_THESIS/Transforming_summ_stats/R_code/scz_risk_ers_prs20180521.txt", row.names=F, col.names=F, append=T, quote=F)
		# rm(risk_ij)
	# }
# }

modeERS <- ers_vec[prob_ers_vec == max(prob_ers_vec)]
modeERS2dp <- 0.27

scz <- read.table(file="/Users/alexg/Documents/PAPERS_THESIS/Transforming_summ_stats/R_code/scz_risk_ers_prs20180521.txt", header=T)
dim(scz)
scz <- data.frame(scz)
ERS_restricted <- scz[scz$specialE != 0, ]
table(ERS_restricted[,3])
ERS_restricted <- data.frame(ERS_restricted)
ERS_restricted$ERS_profiles <- "Minimum (ERS = 0)"
ERS_restricted$ERS_profiles[ERS_restricted$specialE == 2] <- "Median (ERS = 0.72)"
ERS_restricted$ERS_profiles[ERS_restricted$specialE == 3] <- "Maximum (ERS = 1.53)"
ERS_restricted$ERS_profiles[ERS_restricted$specialE == 4] <- "Mode (ERS = 0.27)"

### x = PRS, groupings = ERS
p1 <- ggplot(ERS_restricted) + geom_line(aes(x = PRS, y = risk, linetype=ERS_profiles)) + geom_hline(yintercept=K, colour="blue") + xlab("Polygenic risk score") + ylab("Estimated risk of Schizophrenia") 
p1 + scale_linetype_discrete(name  ="ERS profile")
# ### multi ggplot
# hist_top <- ggplot()+geom_histogram(aes(rnorm(1000000)))
# empty <- ggplot()+geom_point(aes(1,1), colour="white")+
         # theme(axis.ticks=element_blank(), 
               # panel.background=element_blank(), 
               # axis.text.x=element_blank(), axis.text.y=element_blank(),           
               # axis.title.x=element_blank(), axis.title.y=element_blank())

# grid.arrange(hist_top, p1, ncol=1, nrow=2)
###
#ggMarginal(p1, margins="x")
### multi:
ers_vec1 <- ers_vec[order(ers_vec)]
## FYI 95th pc is approx cell number 142 of ers_vec1
risk_range <- range(scz$risk)
scz1 <- scz[scz$specialE == 1, ]
scz2 <- scz[scz$specialE == 2, ]
scz3 <- scz[scz$specialE == 3, ]
scz4 <- scz[scz$specialE == 4, ]
scz1 <- scz1[order(scz1$risk), ]
scz2 <- scz2[order(scz2$risk), ]
scz3 <- scz3[order(scz3$risk), ]
scz4 <- scz4[order(scz4$risk), ]
scz5 <- scz[scz$ERS == ers_vec1[142], ]
scz5 <- scz5[order(scz5$risk), ]
# par(mfrow=c(2,1))
# plot(x=scz1$PRS, y=scz1$risk, type="l", ylim = c(0, 0.4), xlab = "Polygenic risk score", ylab="Risk of Schizophrenia", bty="l", lwd=2)
# points(x=scz2$PRS, y=scz2$risk, type="l", lty=2, lwd=2)
# points(x=scz3$PRS, y=scz3$risk, type="l", lty=3, lwd=2)
# points(x=scz4$PRS, y=scz4$risk, type="l", lty=4, lwd=2)
# curve(dnorm(x, sd=1), from=min(scz$PRS), to=max(scz$PRS))

#nf <- layout(matrix(c(2,0,1,3),2,2,byrow = TRUE), c(3,1), c(1,3), TRUE)
nf <- layout(matrix(c(1,2),2,1,byrow = TRUE), c(4,3), c(3,1), TRUE)
layout.show(nf)
par(mar = c(3,3,1,1))
plot(x=scz1$PRS, y=scz1$risk, type="l", ylim = c(0, 0.4), xlab = "Polygenic risk score", ylab="Risk of Schizophrenia", bty="l", lwd=2)
#points(x=scz2$PRS, y=scz2$risk, type="l", lty=2, lwd=2)
points(x=scz3$PRS, y=scz3$risk, type="l", lty=3, lwd=2)
points(x=scz4$PRS, y=scz4$risk, type="l", lty=4, lwd=2)
leg.txt <- c("Maximum (ERS = 1.53)", "Median (ERS = 0.72)", "Mode (ERS = 0.27)","Minimum (ERS = 0)")
legend("topleft", legend=leg.txt, bty="n", lty=c(3,2,4,1), lwd=2)
par(mar = c(3,3,1,1))
# curve(dnorm(x, sd=sqrt(VML)), from=min(scz$PRS), to=max(scz$PRS), bty="l", lwd=2, yaxt="n", xaxt="n")
curve(dnorm(x, sd=1), from=min(scz$PRS), to=max(scz$PRS), axes=F)
box(bty="l")

risk_prs <- risk_single_PRS_function(prs=sqrt(VML)*unique(scz$PRS), E_ERS=E_ERS, E_PRS = 0, V_PRS=VML, K=K)
scz_noERS <- cbind(unique(scz$PRS), risk_prs)
colnames(scz_noERS) <- c("PRS", "risk")
scz_noERS <- data.frame(scz_noERS)
scz_noERS <- scz_noERS[order(scz_noERS$PRS), ]

risk_ers <- risk_single_ERS_function(ers=unique(scz$ERS), E_ERS=E_ERS, E_PRS = 0, V_ERS= V_ERS, K=K)
scz_noPRS <- cbind(unique(scz$ERS), risk_ers)
colnames(scz_noPRS) <- c("ERS", "risk")
scz_noPRS <- data.frame(scz_noPRS)
scz_noPRS <- scz_noPRS[order(scz_noPRS$ERS), ]

# ### PLOT 1:
# ### Group = ERS, PRS ranging betwen +/- 2:
# scz1res <- scz1[abs(scz1$PRS) <= 2, ]
# nf <- layout(matrix(c(1,2),2,1,byrow = TRUE), c(4,3), c(3,1), TRUE)
# layout.show(nf)
# par(mar = c(3,3,1,1))
# plot(x=scz_noERS$PRS, y=scz_noERS$risk, type="l", ylim = c(0, 0.3), xlab = "Polygenic risk score", ylab="Risk of Schizophrenia", bty="n", axes=F, lwd=3, xlim=c(-2, 2), col="dark grey")
# points(x=scz1$PRS, y=scz1$risk, type="l", lwd=2)
# points(x=scz2$PRS, y=scz2$risk, type="l", lty=2, lwd=2)
# points(x=scz3$PRS, y=scz3$risk, type="l", lty=3, lwd=2)
# points(x=scz4$PRS, y=scz4$risk, type="l", lty=4, lwd=2)
# leg.txt <- c("No ERS", "Maximum (ERS = 1.53)", "Median (ERS = 0.72)", "Mode (ERS = 0.27)","Minimum (ERS = 0)")
# legend("topleft", legend=leg.txt, bty="n", lty=c(1, 3,2,4,1), col = c("dark grey", rep("black", 4)), lwd=c(3,rep(2,4)), text.col = c("dark grey", rep("black", 4)))
# box(bty="l")
# axis(1, labels=F)
# axis(2)
# par(mar = c(3,3,1,1))
# # curve(dnorm(x, sd=1), from=min(scz1res$PRS), to=max(scz1res$PRS), bty="l", lwd=2, yaxt="n", xaxt="n")
# curve(dnorm(x, sd=1), from=min(scz1res$PRS), to=max(scz1res$PRS), axes=F)
# box(bty="l")
# axis(1)

### PLOT 1:
### Group = ERS, PRS ranging betwen +/- 3:
nf <- layout(matrix(c(1,2),2,1,byrow = TRUE), c(4,3), c(3,1), TRUE)
layout.show(nf)
par(mar = c(3,3,1,1))
plot(x=scz_noERS$PRS, y=scz_noERS$risk, type="l", ylim = c(0, 0.3), xlab = "Polygenic risk score", ylab="Risk of Schizophrenia", bty="n", axes=F, lwd=3, xlim=c(-3, 3), col="dark grey")
points(x=scz1$PRS, y=scz1$risk, type="l", lwd=2)
points(x=scz2$PRS, y=scz2$risk, type="l", lty=2, lwd=2)
points(x=scz3$PRS, y=scz3$risk, type="l", lty=3, lwd=2)
points(x=scz4$PRS, y=scz4$risk, type="l", lty=4, lwd=2)
leg.txt <- c("No ERS", "Maximum (ERS = 1.53)", "Median (ERS = 0.72)", "Mode (ERS = 0.27)","Minimum (ERS = 0)")
legend("topleft", legend=leg.txt, bty="n", lty=c(1, 3,2,4,1), col = c("dark grey", rep("black", 4)), lwd=c(3,rep(2,4)), text.col = c("dark grey", rep("black", 4)))
box(bty="l")
axis(1, labels=F)
axis(2)
par(mar = c(3,3,1,1))
# curve(dnorm(x, sd=1), from=min(scz1res$PRS), to=max(scz1res$PRS), bty="l", lwd=2, yaxt="n", xaxt="n")
curve(dnorm(x, sd=1), from=-3, to=3, axes=F)
box(bty="l")
axis(1)

# ### PLOT 1 OPTION 2 (95% CI)
# ### Group = ERS, PRS ranging betwen +/- 2:
# scz1res <- scz1[abs(scz1$PRS) <= 2, ]
# nf <- layout(matrix(c(1,2),2,1,byrow = TRUE), c(4,3), c(3,1), TRUE)
# layout.show(nf)
# par(mar = c(3,3,1,1))
# plot(x=scz_noERS$PRS, y=scz_noERS$risk, type="l", ylim = c(0, 0.2), xlab = "Polygenic risk score", ylab="Risk of Schizophrenia", bty="n", axes=F, lwd=3, xlim=c(-2, 2), col="dark grey")
# points(x=scz1$PRS, y=scz1$risk, type="l", lwd=2, lty=2)
# points(x=scz5$PRS, y=scz5$risk, type="l", lty=2, lwd=2)
# #points(x=scz2$PRS, y=scz2$risk, type="l", lty=2, lwd=2)
# #points(x=scz3$PRS, y=scz3$risk, type="l", lty=3, lwd=2)
# points(x=scz4$PRS, y=scz4$risk, type="l", lwd=2)
# leg.txt <- c("No ERS", "Mode (ERS = 0.27)","ERS between [0, 0.85]")
# legend("topleft", legend=leg.txt, bty="n", lty=c(1, 3,2,4,1), col = c("dark grey", rep("black", 4)), lwd=c(3,rep(2,4)), text.col = c("dark grey", rep("black", 4)))
# box(bty="l")
# axis(1, labels=F)
# axis(2)
# par(mar = c(3,3,1,1))
# # curve(dnorm(x, sd=1), from=min(scz1res$PRS), to=max(scz1res$PRS), bty="l", lwd=2, yaxt="n", xaxt="n")
# curve(dnorm(x, sd=1), from=min(scz1res$PRS), to=max(scz1res$PRS), axes=F)
# box(bty="l")
# axis(1)

### Plots with x = ERS...
####
PRS_restricted <- rbind(scz[scz$PRS == 0, ], scz[abs(scz$PRS) == 1.96, ], scz[abs(scz$PRS) == 1.64, ], scz[abs(scz$PRS) == 0.67, ])
PRS_restricted <- data.frame(PRS_restricted)
ggplot(PRS_restricted) + geom_line(aes(x = ERS, y = risk, linetype=factor(PRS))) + geom_hline(yintercept=K, colour="blue")
### The above is probably better not in ggplot...

### ribbon ### have to reshape data:
PRS1 <- scz[scz$PRS == 0, ]
#colnames(PRS1)[9] <- "Expected"
### 95% CI
PRS2 <- scz[scz$PRS == 1.96, ]
PRS1$High95 <- PRS2$risk
PRS3 <- scz[scz$PRS == -1.96, ]
PRS1$Low95 <- PRS3$risk
### 90% CI
PRS4 <- scz[scz$PRS == 1.64, ]
PRS1$High90 <- PRS4$risk
PRS5 <- scz[scz$PRS == -1.64, ]
PRS1$Low90 <- PRS5$risk
### 50% CI
PRS6 <- scz[scz$PRS == 0.67, ]
PRS1$High50 <- PRS6$risk
PRS7 <- scz[scz$PRS == -0.67, ]
PRS1$Low50 <- PRS7$risk

ggplot(data=PRS1) +  geom_line(aes(x = ERS, y = risk)) + geom_ribbon(aes(x=ERS, ymin=Low90, ymax=High90), alpha = .25) + geom_ribbon(aes(x=ERS, ymin=Low95, ymax=High95), alpha = .25) + geom_ribbon(aes(x=ERS, ymin=Low50, ymax=High50), alpha = .25) + xlab("Environmental risk score") + ylab("Estimated risk of Schizophrenia")

PRS1 <- PRS1[order(PRS1$ERS), ]

pos.gray <- gray(seq(0.1,0.9,length=10))
plot(x = PRS1$ERS, y=PRS1$risk, type="l", lwd = 2, ylim=c(0, 0.4), ylab="Risk of Schizophrenia")
points(x = PRS1$ERS, y=PRS1$Low50, type="l", lwd = 2, lty=2, col=pos.gray[2])
points(x = PRS1$ERS, y=PRS1$High50, type="l", lwd = 2, lty=2, col=pos.gray[2])
points(x = PRS1$ERS, y=PRS1$Low90, type="l", lwd = 2, lty=3, col=pos.gray[5])
points(x = PRS1$ERS, y=PRS1$High90, type="l", lwd = 2, lty=3, col=pos.gray[5])
points(x = PRS1$ERS, y=PRS1$Low95, type="l", lwd = 2, lty=4, col=pos.gray[8])
points(x = PRS1$ERS, y=PRS1$High95, type="l", lwd = 2, lty=4, col=pos.gray[8])

### PLOT 2: Risk versus ERS.
nf <- layout(matrix(c(1,2),2,1,byrow = TRUE), c(4,3), c(3,1), TRUE)
layout.show(nf)
par(mar = c(3,3,1,1))
#plot(x = PRS1$ERS, y=PRS1$risk, type="l", lwd = 2, ylim=c(0, 0.3), ylab="Risk of Schizophrenia", bty="n", axes=F)
plot(x = scz_noPRS$ERS, y= scz_noPRS$risk, type="l", lwd = 3, ylim=c(0, 0.3), ylab="Risk of Schizophrenia", bty="n", axes=F, col="dark grey")
points(x = PRS1$ERS, y=PRS1$risk, type="l", lwd = 2)
#points(x = PRS1$ERS, y=PRS1$Low50, type="l", lwd = 2, lty=2, col=pos.gray[2])
#points(x = PRS1$ERS, y=PRS1$High50, type="l", lwd = 2, lty=2, col=pos.gray[2])
#points(x = PRS1$ERS, y=PRS1$Low90, type="l", lwd = 2, lty=3, col=pos.gray[5])
#points(x = PRS1$ERS, y=PRS1$High90, type="l", lwd = 2, lty=3, col=pos.gray[5])
#points(x = PRS1$ERS, y=PRS1$Low95, type="l", lwd = 2, lty=4, col=pos.gray[8])
#points(x = PRS1$ERS, y=PRS1$High95, type="l", lwd = 2, lty=4, col=pos.gray[8])
points(x = PRS1$ERS, y=PRS1$Low95, type="l", lwd = 2, lty=2, col=pos.gray[2])
points(x = PRS1$ERS, y=PRS1$High95, type="l", lwd = 2, lty=2, col=pos.gray[2])
leg.txt <- c("No PRS", "Mean PRS", "PRS = +/- 1.96")
legend("topleft", legend=leg.txt, lwd=c(3,2,2), lty=c(1,1,2), col = c("dark grey", "black", pos.gray[2]), text.col=c("dark grey", "black", pos.gray[2]), bty="n")
box(bty="l")
axis(1, labels=F)
axis(2)
par(mar = c(3,3,1,1))
#plot(ers_vec[order(ers_vec)], prob_ers_vec[order(ers_vec)], type="h", lwd=2, bty="l", xaxt="n", yaxt = "n")
plot(ers_vec[order(ers_vec)], prob_ers_vec[order(ers_vec)], type="h", lwd=2, axes=F)
box(bty="l")
axis(1)
axis(2, labels = c("0.00", "0.11"), at=range(prob_ers_vec))

### Joint risk ERS == min, PRS == 0.99L:
risk_function(ers=0, prs=-1*qnorm(0.99)*sqrt(VML), E_ERS= E_ERS, V_ERS= V_ERS, V_PRS=VML, K=K)
### Joint risk ERS == max, PRS == 0.99U:
risk_function(ers=max(ers_vec), prs=1*qnorm(0.99)*sqrt(VML), E_ERS= E_ERS, V_ERS= V_ERS, V_PRS=VML, K=K)
### Joint risk ERS == med = 0.71, PRS == 0:
risk_function(ers=0.71, prs=0, E_ERS= E_ERS, V_ERS= V_ERS, V_PRS=VML, K=K)
### Joint risk ERS == mode = 0.27, PRS == 0:
risk_function(ers=0.27, prs=0, E_ERS= E_ERS, V_ERS= V_ERS, V_PRS=VML, K=K)

### Joint risk plot...
### 1. we need a population sample of PRS, ERS and risk
N = 1500000
E_ers <- sum(ers_vec*prob_ers_vec)
prs_sample <- rnorm(N, sd=sqrt(VML))
ers_sample <- rep(0, N)
risk_sample <- rep(0, N)
ers_risk_sample <- rep(0, N)
prs_risk_sample <- rep(0, N)
risk_prs0 <- rep(0, N)
risk_prsU95 <- rep(0, N)
risk_prsL95 <- rep(0, N)
risk_ers0 <- rep(0, N)
risk_ers0.27 <- rep(0, N)
risk_ers1.53 <- rep(0, N)
risk_mean_ers <- rep(0, N)
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
	risk_mean_ers[i] <- risk_function(ers=E_ers, prs=prs_sample[i],E_ERS= E_ERS, V_ERS= V_ERS, V_PRS=VML, K=K)
}
pop_sample <- cbind(prs_sample, ers_sample, risk_sample, ers_risk_sample, prs_risk_sample, risk_prs0, risk_prsL95, risk_prsU95, risk_ers0, risk_ers0.27, risk_ers1.53, risk_mean_ers)
pop_sample <- data.frame(pop_sample)
##rm(prs_sample, ers_sample, risk_sample)
colnames(pop_sample)[1:5] <- c("prs", "ers", "risk", "ers_risk", "prs_risk")
pop_ers <- pop_sample
pop_prs <- pop_sample
pop_sample$rs <- pop_sample$prs + pop_sample$ers
pop_sample <- pop_sample[order(pop_sample$rs), ]
Fn_pop_rs <- ecdf(pop_sample$rs)
pop_sample$pc_rs <- Fn_pop_rs(pop_sample$rs)

#plot(x=pop_sample$pc_rs, y=pop_sample$risk/K, type="l", ylab= "Population relative risk", xlab="Joint risk score percentile")
#abline(h = 1)
#plot(x=pop_sample$rs, y=pop_sample$risk/K, type="l", ylab= "Population relative risk", xlab="Joint risk score")
#abline(h = 1)
quantile(pop_sample$risk)
quantile(scz$risk)
#approxK <- pop_sample[round(pop_sample$risk, digits = 8) == 0.01000000, ]
#range(approxK$rs)
### 0.5523135
risk_function(ers=0.5523135/2, prs=0.5523135/2, E_ERS= E_ERS, V_ERS= V_ERS, V_PRS=VML, K=K)
risk_function(ers=0.5523135, prs=0, E_ERS= E_ERS, V_ERS= V_ERS, V_PRS=VML, K=K)
Fn_pop_rs(0.5523135) ## 0.681665

pop_ers <- pop_ers[order(pop_ers$ers), ]
Fn_pop_ers <- ecdf(pop_ers$ers)
pop_ers$pc_ers <- Fn_pop_ers(pop_ers$ers)

pop_prs <- pop_prs[order(pop_prs$prs), ]
Fn_pop_prs <- ecdf(pop_prs$prs)
pop_prs$pc_prs <- Fn_pop_prs(pop_prs$prs)

plot(x=pop_sample$pc_rs, y=pop_sample$risk/K, type="l", ylab= "Population relative risk", xlab="Risk score percentile")
abline(h = 1)
points(x=pop_ers$pc_ers, y = pop_ers$ers_risk/K, type="l", col="red")
points(x=pop_prs$pc_prs, y = pop_prs$prs_risk/K, type="l", col="blue")
points(y = max(pop_prs$prs_risk/K), x=1, col = "Blue")
points(y = max(pop_sample$risk/K), x=1)
points(y = max(pop_ers$ers_risk/K), x=1, col = "red")

plot(x=pop_sample$pc_rs[pop_sample$pc_rs > 0.9], y=pop_sample$risk[pop_sample$pc_rs > 0.9]/K, type="l", ylab= "Population relative risk", xlab="Risk score percentile", ylim = c(0, 30))
#abline(h = 1)
points(x=pop_ers$pc_ers[pop_ers$pc_ers > 0.9], y = pop_ers$ers_risk[pop_ers$pc_ers > 0.9]/K, type="l", col="red")
points(x=pop_prs$pc_prs[pop_prs$pc_prs > 0.9], y = pop_prs$prs_risk[pop_prs$pc_prs > 0.9]/K, type="l", col="blue")
points(y = max(pop_prs$prs_risk/K), x=1, col = "Blue")
points(y = max(pop_sample$risk/K), x=1)
points(y = max(pop_ers$ers_risk/K), x=1, col = "red")

# plot(x=pop_ers$pc_ers, y = pop_ers$ers_risk/K, type="l")
# abline(h = 1)
# plot(x=pop_prs$pc_prs, y = pop_prs$prs_risk/K, type="l")
# abline(h = 1)

risk.df <- data.frame(c(rep("Joint", N), rep("PRS only", N), rep("ERS only", N)), c(pop_sample$risk, pop_prs$prs_risk, pop_ers$ers_risk))
colnames(risk.df) <- c("Model", "Risk")

jointrisk0.9 <- quantile(pop_sample$risk, probs=c(0.9))
prsrisk0.9 <- quantile(pop_prs$prs_risk, probs=c(0.9))
ersrisk0.9 <- quantile(pop_ers$ers_risk, probs=c(0.9))
ggplot(data=risk.df) + geom_density(aes(x = Risk, colour=Model)) + geom_vline(xintercept = jointrisk0.9)

df.sample90 <- pop_sample[pop_sample$risk >= jointrisk0.9, ]
df.prs90 <- pop_prs[pop_prs$prs_risk >= prsrisk0.9, ]
df.ers90 <- pop_ers[pop_ers$ers_risk >= ersrisk0.9, ]

colnames(df.sample90)
dim(df.sample90)
dim(df.prs90)
dim(df.ers90)

range(df.sample90$ers)
range(df.prs90$ers)
range(df.ers90$ers)
hist(df.sample90$ers)
hist(df.prs90$ers)
hist(df.ers90$ers, xlim=c(0, 1.54))

range(df.sample90$prs)
range(df.prs90$prs)
range(df.ers90$prs)
hist(df.sample90$prs)
hist(df.prs90$prs)
hist(df.ers90$prs)

df.sample90$Model <- "Joint"
df.prs90$Model <- "PRS only"
df.ers90$Model <- "ERS only"

df.sample90r <- cbind(df.sample90[,1:3], df.sample90$Model) 
df.prs90r <- cbind(df.prs90[,1:2], df.prs90$prs_risk, df.prs90$Model) 
df.ers90r <- cbind(df.ers90[,1:2], df.ers90$ers_risk, df.ers90$Model) 
colnames(df.sample90r) <- c("prs", "ers", "risk", "Model")
colnames(df.prs90r) <- c("prs", "ers", "risk", "Model")
colnames(df.ers90r) <- c("prs", "ers", "risk", "Model")

df.all90 <- rbind(df.sample90r, df.prs90r, df.ers90r)
#ggplot(data=df.all90) + geom_density(aes(x = prs, colour=Model, fill=Model), alpha = 0.25)
#ggplot(data=df.all90) + stat_density(aes(x = prs, colour=Model, fill=Model), alpha = 0.25)
#ggplot(data=df.all90) + geom_histogram(aes(x = prs, colour=Model, fill=Model), alpha = 0.25, bins=100)
ggplot(data=df.all90) + geom_freqpoly(aes(x = prs, colour=Model, fill=Model), bins=100)

ggplot() + geom_histogram(aes(x=df.sample90r$prs), bins=200)+ xlim(-3,3) + geom_histogram(aes(x=df.prs90r$prs), col = "Red", bins=200) + geom_histogram(aes(x=df.ers90r$prs), col = "Blue", bins=200)

ggplot() + geom_freqpoly(aes(x=df.sample90r$prs), bins=200)+ xlim(-3,3) + geom_freqpoly(aes(x=df.prs90r$prs), col = "Red", bins=200) + geom_freqpoly(aes(x=df.ers90r$prs), col = "Blue", bins=200)


risk_single_ERS_function(ers = 0.27, E_ERS = E_ERS, V_ERS=V_ERS, K=K)
risk_single_PRS_function(prs = 0, E_ERS = E_ERS, V_PRS=VML, K=K)
risk_function(ers= 0.27, prs=0, E_ERS= E_ERS, V_ERS= V_ERS, V_PRS=VML, K=K)
risk_function(ers= 0.2737725, prs=0, E_ERS= E_ERS, V_ERS= V_ERS, V_PRS=VML, K=K)
risk_function(ers= E_ERS, prs=0, E_ERS= E_ERS, V_ERS= V_ERS, V_PRS=VML, K=K)
risk_function(ers=0, prs=0, E_ERS= E_ERS, V_ERS= V_ERS, V_PRS=VML, K=K)
(risk_function(ers=0, prs=0.1, E_ERS= E_ERS, V_ERS= V_ERS, V_PRS=VML, K=K) - risk_function(ers=0, prs=0, E_ERS= E_ERS, V_ERS= V_ERS, V_PRS=VML, K=K))/risk_function(ers=0, prs=0, E_ERS= E_ERS, V_ERS= V_ERS, V_PRS=VML, K=K)
(risk_function(ers= 0.2737725, prs=0.1, E_ERS= E_ERS, V_ERS= V_ERS, V_PRS=VML, K=K) - risk_function(ers= 0.2737725, prs=0, E_ERS= E_ERS, V_ERS= V_ERS, V_PRS=VML, K=K))/risk_function(ers= 0.2737725, prs=0, E_ERS= E_ERS, V_ERS= V_ERS, V_PRS=VML, K=K)
(risk_function(ers= 1.532693, prs=0.1, E_ERS= E_ERS, V_ERS= V_ERS, V_PRS=VML, K=K) - risk_function(ers= 1.532693, prs=0, E_ERS= E_ERS, V_ERS= V_ERS, V_PRS=VML, K=K))/risk_function(ers= 1.532693, prs=0, E_ERS= E_ERS, V_ERS= V_ERS, V_PRS=VML, K=K)

risk_function(ers= 1.532693, prs=0, E_ERS= E_ERS, V_ERS= V_ERS, V_PRS=VML, K=K)
risk_function(ers=0.27, prs=-1.96*sqrt(VML), E_ERS= E_ERS, V_ERS= V_ERS, V_PRS=VML, K=K)
risk_function(ers=0.27, prs=1.96*sqrt(VML), E_ERS= E_ERS, V_ERS= V_ERS, V_PRS=VML, K=K)
risk_function(ers= 0.2737725, prs=-1.96*sqrt(VML), E_ERS= E_ERS, V_ERS= V_ERS, V_PRS=VML, K=K)
risk_function(ers= 0.2737725, prs=1.96*sqrt(VML), E_ERS= E_ERS, V_ERS= V_ERS, V_PRS=VML, K=K)
## RR:
risk_single_PRS_function(prs = 0, E_ERS = E_ERS, V_PRS=VML, K=K)/K
risk_function(ers=0, prs=0, E_ERS= E_ERS, V_ERS= V_ERS, V_PRS=VML, K=K)/K
risk_function(ers= 1.532693, prs=0, E_ERS= E_ERS, V_ERS= V_ERS, V_PRS=VML, K=K)/K

risk_single_ERS_function(ers = 0.27, E_ERS = E_ERS, V_ERS=V_ERS, K=K)/K
risk_function(ers=0.27, prs=-1.96*sqrt(VML), E_ERS= E_ERS, V_ERS= V_ERS, V_PRS=VML, K=K)/K
risk_function(ers=0.27, prs=1.96*sqrt(VML), E_ERS= E_ERS, V_ERS= V_ERS, V_PRS=VML, K=K)/K
# plot(x=pop_sample$pc_rs, y=pop_sample$ers_risk/K, type="l", col="red")
# plot(x=pop_sample$pc_rs, y=pop_sample$prs_risk/K, type="l", col="blue")

# plot(x=pop_sample$pc_rs, y=pop_sample$risk/K, type="l", ylab= "Population relative risk", xlab="Risk score percentile")
# points(x=pop_sample$pc_rs, y=pop_sample$prs_risk/K, col="blue")
# points(x=pop_sample$pc_rs, y=pop_sample$ers_risk/K, col="red")
# abline(h=1)
# points(x=pop_sample$pc_rs, y=pop_sample$risk/K, type="l", lwd=2)

### PLOT 3: ERS percentile
plot(x=pop_ers$pc_ers, y = pop_ers$ers_risk/K, pch=20, ylim=c(0,25), ylab="Population relative risk", xlab="ERS percentile", lwd=2, bty="l")
#abline(h = 1)
points(x=pop_ers$pc_ers, y = pop_ers$risk_prs0/K, type="l", lwd=2)
points(x=pop_ers$pc_ers, y = pop_ers$risk_prsL95/K, type="l", col = pos.gray[7], lwd=2)
points(x=pop_ers$pc_ers, y = pop_ers$risk_prsU95/K, type="l", col = pos.gray[4], lwd=2)
leg.txt <- c("No PRS", "PRS = 0", "PRS = 1.96", "PRS = -1.96")
legend("topleft", legend=leg.txt, col=c("black", "black", pos.gray[4], pos.gray[7]), bty="n", text.col=c("black","black", pos.gray[4], pos.gray[7]), lty=c(NA, rep(1,3)), pch=c(20, rep(NA, 3)), lwd=c(2,2,2,2))

### 
#max(pop_prs$risk_ers1.53/K)
plot(x=pop_prs$pc_prs, y = pop_prs$prs_risk/K, pch=20, ylim=c(0,70), ylab="Population relative risk", xlab="PRS percentile", lwd=2, bty="l")
#abline(h = 1)
points(x=pop_prs$pc_prs, y = pop_prs$risk_ers0.27/K, type="l", col = "black", lwd=2)
#points(x=pop_prs$pc_prs, y = pop_prs$risk_mean_ers/K, type="l", col = pos.gray[3])
points(x=pop_prs$pc_prs, y = pop_prs$risk_ers0/K, type="l", col = pos.gray[7], lwd=2)
points(x=pop_prs$pc_prs, y = pop_prs$risk_ers1.53/K, type="l", col = pos.gray[4], lwd=2)
#leg.txt <- c("no ERS", "Maximum ERS", "Mean ERS", "Minimum ERS")
#legend("topleft", legend=leg.txt, col=c("black", pos.gray[3], pos.gray[5], pos.gray[7]), bty="n", text.col=c("black", pos.gray[3], pos.gray[5], pos.gray[7]), , lwd=c(3,1,1,1))
leg.txt <- c("No ERS", "Most frequent ERS", "Maximum ERS", "Minimum ERS")
legend("topleft", legend=leg.txt, col=c("black","black", pos.gray[4], pos.gray[7]), bty="n", text.col=c("black","black", pos.gray[4], pos.gray[7]), lty=c(NA,1,1,1),pch=c(20, NA,NA,NA), lwd=rep(2,4))
# scz$jointRS <- (scz$PRS*sqrt(VML)) + scz$ERS
# scz_pc <- scz[order(scz$jointRS), ]
# Fn_jointRS <- ecdf(scz_pc$jointRS)
# scz_pc$pcjointRS <- Fn_jointRS(scz_pc$jointRS)
# plot(x=scz_pc$jointRS, y=scz_pc$risk/K, type="l")
# plot(x=scz_pc$pcjointRS, y=scz_pc$risk/K, type="l")
# ### Plot option 3:
# risk_ers_only <- rep(0, length(ers_vec))
# for(i in 1:length(ers_vec)){
	# risk_ers_only[i] <- risk_single_ERS_function(ers = ers_vec[i], E_ERS = E_ERS, V_ERS=V_ERS, K=K)
# }
# risk_prs_only <- rep(0, length(prsI_L))
# for(i in 1:length(prsI_L)){
	# risk_prs_only[i] <- risk_single_PRS_function(prs = prsI_L[i], E_ERS = E_ERS, V_PRS=VML, K=K)
# }
# ### percentiles:
# ers_vec1 <- ers_vec[order(ers_vec)]
# risk_ers_only1<- risk_ers_only[order(ers_vec)]
# ers_pc <- 


### Approx 95% of sampled individuals: risk
colnames(pop_sample)
### PRS:
quantile(pop_sample$prs_risk, probs=c(0.025, 0.975))
### ERS:
quantile(pop_sample$ers_risk, probs=c(0.025, 0.975)) ### 0.0025 to 0.0385 (0.0379)
quantile(pop_sample$ers_risk, probs=c(0.0, 0.95)) ### 0.0025 to 0.0267
### PRS+ERS:
quantile(pop_sample$risk, probs=c(0.025, 0.975))
### ERS value- not risk
quantile(pop_sample$ers, probs=c(0.025, 0.975))
quantile(pop_sample$ers, probs=c(0.94, 0.95, 0.96))
quantile(pop_sample$ers, probs=c(0))
ers_vec1 <- ers_vec[order(ers_vec)]
prob_vec1 <- prob_ers_vec[order(ers_vec)]
cum_prob <- cumsum(prob_vec1)
plot(ers_vec1, cum_prob)
abline(h = 0.025, col="red")
abline(h = 0.975, col="red")
abline(h = 0.95, col="blue")
hist(pop_sample$rs)
hist(pop_sample$ers_risk)
### CDF normal logistic plot...
#curve(pnorm, from=-4, to=4)
#curve(plogis, from=-4, to=4)
s1 = 1/1.702
s_wilson <- 1/1.86
s_JK <- 1
s_t <- 1/1.595769
xrange <- seq(from = -3, to = 3, by = 0.01)
CDF.norm <- pnorm(xrange)
CDF.logis <- plogis(xrange, s=sqrt(3)/pi)
CDF.logis1 <- plogis(xrange, s=s1)
CDF.logis.wilson <- plogis(xrange, s=s_wilson)
CDF.logis.JK <- plogis(xrange, s=s_JK)
CDF.logis.t <- plogis(xrange, s=s_t)

plot(x=xrange, y=CDF.norm, bty="l", xlab="x", ylab="Probability", type="l", lwd=2)
points(x=xrange, y=CDF.logis, col="dark grey", type="l", lty=2, lwd=2)
points(x=xrange, y=CDF.logis1, col="red", type="l", lty=2, lwd=2)
#points(x=xrange, y=CDF.logis.wilson, col="green", type="l", lty=2, lwd=2)
points(x=xrange, y=CDF.logis.JK, col="blue", type="l", lty=2, lwd=2)
points(x=xrange, y=CDF.logis.t, col="purple", type="l", lty=2, lwd=2)

leg.txt <- c("Standard normal CDF", "Standardised logistic CDF")
legend("topleft", legend=leg.txt, col=c("black", "dark grey"), lty=c(1,2), lwd=c(2,2), text.col=c("black", "dark grey"), bty="n")

