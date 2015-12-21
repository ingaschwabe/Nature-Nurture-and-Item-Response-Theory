#y_dz = Item responses of DZ twins (matrix)
#y_mz = Item responses of MZ twins (matrix)
#n_mz = Number of MZ twin pairs
#n_dz = Number of DZ twin pairs
#n_items = Number of phenotypic items administered
#b = Vector with item difficulty parameters, assumed known here 
#x_MZ = moderator variable values for every MZ twin (matrix)
#x_DZ = moderator variable values for every DZ twin  (matrix)

#Required structure of the y_mz/y_dz data matrix:
#y_dz[i,k] = kth datapoint from the ith DZ twin pair
#y_mz[i,k] = kth datapoint from the ith MZ twin pair
#This results in a matrix of n_mz (or, in case of y_dz, n_dz)
#rows and 2*n_items columns. e.g. y_mz[1,22] is the response 
#of MZ twin 1 from family 1 to item 22 if n_items = 22

#Required structure of the x_MZ/x_DZ data matrix: 
#x_MZ[i,1:2] = Vector with moderator variable values 
#for the first (x_MZ[i,1]) and second (x_MZ[i,2])
#twin of every ith MZ twin pair. 
#x_DZ[i,1:2] = Vector with moderator variable values 
#for the first (x_MZ[i,1]) and second (x_MZ[i,2])
#twin of every ith MZ twin pair. 

#JAGS uses precision parameters for the variance parameters. 
#Therefore, after running the script, these precision parameters 
#have to be inverted. For example:
#var_a <-  1/outputAnalysis$tau_a[,,1] with the rjags package

model{
##MZ twins
for (fam in 1:n_mz){	
  c_mz[fam] ~ dnorm(0, 1) 
  a_mz[fam] ~ dnorm(0, 1) 
    
  tau_c_mz[fam,1] <- exp(beta_0c + beta_1c*x_MZ[fam,1]) 
  tau_c_mz[fam,2] <- exp(beta_0c + beta_1c*x_MZ[fam,2]) 
        
  tau_a_mz[fam,1] <- exp(beta_0a + beta_1a*x_MZ[fam,1])
  tau_a_mz[fam,2] <- exp(beta_0a + beta_1a*x_MZ[fam,2])
        
  tau_e_mz[fam,1] <- 1/(exp(beta_0e + beta_1e*x_MZ[fam,1] )) 
  tau_e_mz[fam,2] <- 1/(exp(beta_0e + beta_1e*x_MZ[fam,2] ))
        
  c_mz_twin1[fam] <- c_mz[fam] * sqrt(tau_c_mz[fam,1])
  c_mz_twin2[fam] <- c_mz[fam] * sqrt(tau_c_mz[fam,2])
        
  a_mz_twin1[fam] <- a_mz[fam] * sqrt(tau_a_mz[fam,1])
  a_mz_twin2[fam] <- a_mz[fam] * sqrt(tau_a_mz[fam,2])
               
  pheno_mz[fam,1] ~ dnorm(mu_mz + beta_1m_mz * x_MZ[fam,1] + 
                          beta_1mc_mz * x_MZ[fam,2] + 
                          c_mz_twin1[fam] + a_mz_twin1[fam],
                          tau_e_mz[fam,1])
  pheno_mz[fam,2] ~ dnorm(mu_mz + beta_1m_mz * x_MZ[fam,2] + 
                          beta_1mc_mz * x_MZ[fam,1] + 
                          c_mz_twin2[fam] + a_mz_twin2[fam],
                          tau_e_mz[fam,2])
                
  #### Measurement model for phenotypic item data: 
  for (k in 1:n_items){
    logit(p[fam,k]) <- pheno_mz[fam,1] - b[k]
    y_mz[fam,k] ~ dbern(p[fam,k])
  }   		
        
  for (k in (n_items+1):(2*n_items)){
    logit(p[fam,k]) <- pheno_mz[fam,2] - b[k-n_items]
    y_mz[fam,k] ~ dbern(p[fam,k])
  }
} #end MZ twins 									  
    
##DZ twins
for (fam in 1:n_dz){
  c_dz[fam] ~ dnorm(0,1) 
  a1_dz[fam] ~ dnorm(0,2) 
  a2_dz[fam,1] ~ dnorm(a1_dz[fam], 2)
  a2_dz[fam,2] ~ dnorm(a1_dz[fam], 2)
        
  tau_c_dz[fam,1] <- exp(beta_0c + beta_1c*x_DZ[fam,1])
  tau_c_dz[fam,2] <- exp(beta_0c + beta_1c*x_DZ[fam,2])
        
  tau_a_dz[fam,1] <- exp(beta_0a + beta_1a*x_DZ[fam,1])
  tau_a_dz[fam,2] <- exp(beta_0a + beta_1a*x_DZ[fam,2])
        
  tau_e_dz[fam,1] <- 1/(exp(beta_0e + beta_1e*x_DZ[fam,1] ))
  tau_e_dz[fam,2] <- 1/(exp(beta_0e + beta_1e*x_DZ[fam,2] ))
        
  c_dz_twin1[fam] <- c_dz[fam] * sqrt(tau_c_dz[fam,1])
  c_dz_twin2[fam] <- c_dz[fam] * sqrt(tau_c_dz[fam,2])
        
  a_dz_twin1[fam] <- a2_dz[fam,1] * sqrt(tau_a_dz[fam,1])
  a_dz_twin2[fam] <- a2_dz[fam,2] * sqrt(tau_a_dz[fam,2])
        
  pheno_dz[fam,1] ~ dnorm(mu_dz + beta_1m_dz * x_DZ[fam,1] + 
                          beta_1mc_dz * x_DZ[fam,2] + 
                          c_dz_twin1[fam] + a_dz_twin1[fam],
                          tau_e_dz[fam,1])
  pheno_dz[fam,2] ~ dnorm(mu_dz + beta_1m_dz * x_DZ[fam,2] + 
                          beta_1mc_dz * x_DZ[fam,1] + 
                          c_dz_twin2[fam] + a_dz_twin2[fam],
                          tau_e_dz[fam,2])
						  
  #### Measurement model for phenotypic item data: 
  for (k in 1:n_items){
    logit(p_dz[fam,k]) <- pheno_dz[fam,1] - b[k]
    y_dz[fam,k] ~ dbern(p_dz[fam,k])
  }   		
        
  for (k in (n_items+1):(2*n_items)){
    logit(p_dz[fam,k]) <- pheno_dz[fam,2] - b[k-n_items]
    y_dz[fam,k] ~ dbern(p_dz[fam,k])
  }
        
} #end DZ twins 
    
#Priors
mu_mz ~ dnorm(0, .1) 
mu_dz ~ dnorm(0, .1)
    
beta_1m_mz ~ dnorm(0, .1)
beta_1m_dz ~ dnorm(0, .1)
beta_1mc_mz ~ dnorm(0, .1)
beta_1mc_dz ~ dnorm(0, .1)
    
beta_1a ~ dnorm(0, .1)
beta_1c ~ dnorm(0, .1)
beta_1e ~ dnorm(0, .1)
    
beta_0a ~ dnorm(-1, .5)
beta_0c ~ dnorm(-1, .5)
beta_0e ~ dnorm(-1, .5)
}