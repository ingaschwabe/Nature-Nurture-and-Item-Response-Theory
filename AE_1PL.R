#y_dz = Item responses of DZ twins (matrix)
#y_mz = Item responses of MZ twins (matrix)
#n_mz = Number of MZ twin pairs, 
#n_dz = Number of DZ twin pairs,
#n_items = Number of items administered
#b = item parameters, assumed known in the analysis

#Required structure of the y_dz/y_mz data matrix:
#y_dz[i,k] = kth datapoint from the ith DZ twin pair
#y_mz[i,k] = kth datapoint from the ith MZ twin pair

#This results in a matrix of n_mz (or, in 
#case of y_dz, n_dz) rows and 2*n_items columns
#(e.g., y_mz[1,22] is the response of 
#MZ twin 1 from family 1 to item 22 

#When item parameters are unknown in the analysis,
#following code can be integrated in the script:

#for (i in 1:n.items){
#  b[i] ~ dnorm(0,.1)
#  }
#Mu then has to be set to zero to identify the scale

#JAGS uses precision parameters for the variance
#parameters. Therefore, after running the script, 
#these precision parameters should be inverted. 
#e.g.: var_c = 1/outputAnalysis$tau_c[,,1] 
#with the rjags package

model{

##MZ twins
for (fam in 1:n_mz){
  c_mz[fam] ~ dnorm(mu, tau_c)
  f_mz[fam] ~ dnorm(c_mz[fam], tau_a) 
  a_mz[fam] <- f_mz[fam] - c_mz[fam] 
  tau_e[fam] <- 1/(exp(beta0 + (beta1*a_mz[fam])))
  
 for (twin in 1:2){
   pheno_mz[fam,twin]~ dnorm(f_mz[fam],tau_e[fam])   
  }
   		
#1pl model twin1
 for (k in 1:n_items){
   logit(p[fam,k]) <- pheno_mz[fam,1] - b[k]
   y_mz[fam,k] ~ dbern(p[fam,k])
  }   		

#1pl model twin2
 for (k in (n_items+1):(2*n_items)){
	 logit(p[fam,k]) <- pheno_mz[fam,2] - b[k-n_items]
	 y_mz[fam,k] ~ dbern(p[fam,k])
 }
}

##DZ twins
for (fam in 1:ndz){
  c_dz[fam] ~ dnorm(mu, tau_c)
  f0_dz[fam] ~ dnorm(c_dz[fam], doubletau_a)
				
 for (twin in 1:2){										
	 f_dz[fam,twin] ~ dnorm(f0_dz[fam], doubletau_a)
	 a_dz[fam,twin] <- f_dz[fam,twin] - c_dz[fam]
	 tau_e_dz[fam,twin] <- 1/(exp(beta0 + 
                               (beta1*a_dz[fam,twin])))
	 pheno_dz[fam,twin] ~ dnorm(f_dz[fam,twin], 
                              tau_e_dz[fam,twin])
 }

#1pl model twin1 (DZ)
 for (k in 1:n.items){
  logit(p2[fam,k]) <- pheno_dz[fam,1] - b[k]
   Ydz[fam,k] ~ dbern(p2[fam,k])
 }

#1pl model twin2 (DZ)
 for (k in (n_items+1):(2*n_items)){
	 logit(p2[fam,k]) <- pheno_dz[fam,2] - b[k-n_items]
	 Ydz[fam,k] ~ dbern(p2[fam,k])
	}

}

#Prior distributions
mu ~ dnorm(0,.1)
beta1 ~ dnorm(0,.1) 
beta0 ~ dnorm(0,1)

doubletau_a <- 2*tau_a
tau_a ~ dgamma(1,1)   
tau_c ~ dgamma(1,1)
}