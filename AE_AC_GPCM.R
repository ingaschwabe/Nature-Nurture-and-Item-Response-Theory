#y_dz = Item responses of DZ twins (matrix)
#y_mz = Item responses of MZ twins (matrix)
#n_mz = Number of MZ twin pairs
#n_dz = Number of DZ twin pairs
#n_items = Number of items administered

#Required structure of the y_dz/y_mz data matrix:
#y_dz[i,k] = kth datapoint from the ith DZ twin pair
#y_mz[i,k] = kth datapoint from the ith MZ twin pair

#This results in a matrix of n_mz (or, in case of y_dz, n_dz)
#rows and 2*n_items columns. e.g. y_mz[1,22] is the response 
#of MZ twin 1 from family 1 to item 22 if n_items = 22

#JAGS uses precision parameters for the variance parameters. 
#Therefore, after running the script, these precision parameters 
#have to be inverted. For example:
#var_a <-  1/outputAnalysis$tau_a[,,1] with the rjags package

model{
#MZ twins
for (i in 1:n_mz){
  c_mz[i] ~ dnorm(0, tau_c_mz[i])
  a_mz[i] ~ dnorm(0, tau_a)
  
  tau_c_mz[i] <- 1/(exp(gamma0 + gamma1*a_mz[i]))
  tau_e_mz[i] <- 1/(exp(beta0 + beta1*a_mz[i]))
        
  #Phenotypic values: 
  mz[i,1] ~ dnorm(a_mz[i] + c_mz[i], tau_e_mz[i]) 
  mz[i,2] ~ dnorm(a_mz[i] + c_mz[i], tau_e_mz[i]) 
   
  for (j in 1:n_items){
    for (k in 1:3){
	  eta[i,j,k] <- alpha[j] * 
                    (mz[i,1]-beta[j,k])
	  psum[i,j,k] <- sum(eta[i,j,1: k])
	  exp_psum[i,j,k] <- exp(psum [i,j,k])
	  prob[i,j,k] <- exp_psum[i,j,k]/sum(exp_psum [i,j,1:3])
     }	 
   }
                
  for (j in (n_items+1):(2*n_items)){
    for (k in 1:3){
	  eta[i,j,k] <- alpha[j-n_items] * 
                    (mz[i,2]-beta[j-n_items,k])
	  psum [i,j,k] <- sum(eta[i,j,1: k])
	  exp_psum [i,j,k] <- exp( psum [i , j , k])
	  prob [i,j,k] <- exp_psum[i,j,k]/sum(exp_psum [i,j,1:3])
    } 
  }
                
  for (j in 1:(2*n_items)){
    y_mz[i,j] ~ dcat (prob[i,j,1:3])
  }

} #end MZ twins									  

#DZ twins
for (i in 1:n_dz){					
  c_dz[i] ~ dnorm(0, 1)

  a1_dz[i] ~ dnorm(0, doubletau_a) 
  a2_dz[i,1] ~ dnorm(a1_dz[i], doubletau_a)
  a2_dz[i,2] ~ dnorm(a1_dz[i], doubletau_a)

  tau_c_dz[i,1] <- exp(gamma0 + gamma1*a2_dz[i,1])
  tau_c_dz[i,2] <- exp(gamma0 + gamma1*a2_dz[i,2])
  
  c_dz_twin1[i] <- c_dz[i] * sqrt(tau_c_dz[i,1])
  c_dz_twin2[i] <- c_dz[i] * sqrt(tau_c_dz[i,2])
 
  tau_e_dz[i,1] <- 1/(exp(beta0 + beta1*a2_dz[i,1]))
  tau_e_dz[i,2] <- 1/(exp(beta0 + beta1*a2_dz[i,2]))
	
  dz[i,1] ~ dnorm(a2_dz[i,1] + c_dz_twin1[i], tau_e_dz[i,1])
  dz[i,2] ~ dnorm(a2_dz[i,2] + c_dz_twin2[i], tau_e_dz[i,2])
	
  for (j in 1:n_items){
    for (k in 1:3){
	  etadz[i,j,k] <- alpha [j] * (dz[i,1]-beta [j,k])
	  psumdz [i,j,k] <- sum(etadz[i,j,1: k])
	  exp_psumdz[i,j,k] <- exp(psumdz [i,j,k])
	  probdz [i,j,k] <- exp_psumdz [i,j,k]/
                        sum(exp_psumdz[i,j,1:3])
    } 
  }
                
   for (j in (n_items+1):(2*n_items)){
     for (k in 1:3){
	  etadz[i,j,k] <- alpha[j-n_items] * 
                     (dz[i,2]-beta [j-n_items,k])
	  psumdz [i,j,k] <- sum(etadz[i,j,1:k])
	  exp_psumdz [i,j,k] <- exp(psumdz [i,j,k])
	  probdz [i,j,k] <- exp_psumdz [i,j,k]/
                        sum(exp_psumdz [i,j,1:3])
    }
  }
		
  for (j in 1:(2*n_items)){
    y_dz[i,j] ~ dcat (probdz[i,j,1:3])
   }
} #end DZ twins

#DZ twins genetic correlation 0.5:
doubletau_a <- 2*tau_a    

#Set alpha of item 3 to 1 to identify the scale
alpha[3] <- 1

#for the rest of the alpha parameters: 
#lognormal prior with expectation of 0 
#and variance of 10
alpha[1] ~ dlnorm(0, .1)
alpha[2] ~ dlnorm(0, .1)

for (j in 4:n_items){
  alpha[j] ~ dlnorm(0, .1)
}

#Beta IRT parameters: 
for (j in 1:n_items){
  beta [j , 1] <- 0.0
  for (k in 2:3){
    beta [j , k] ~ dnorm (0, .1)
  }
}

#Priors distributions: 
tau_a ~ dgamma(1,1)
beta0 ~ dnorm(-1,.5)
beta1 ~ dnorm(0,.1)
gamma0 ~ dnorm(-1,.5)
gamma1 ~ dnorm(0,.1)
}