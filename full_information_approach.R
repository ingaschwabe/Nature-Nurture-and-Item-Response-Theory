#Install OpenMx
source('http://openmx.psyc.virginia.edu/getOpenMx.R')
install.packages("OpenMx") #Install OpenMx package

#Load package
library(OpenMx)

#OpenMx analysis
twin_ACE_cov <- mxModel("twinACE",
 #Matrices X,Y,Z to store a,c,e path coefficients
 mxMatrix(type="Full", nrow=1, ncol=1, free=TRUE, 
          values=.6, label="a", name="X" ),
 mxMatrix(type="Full", nrow=1, ncol=1, free=TRUE, 
          values=.6, label="c", name="Y" ),
 mxMatrix(type="Full", nrow=1, ncol=1, free=TRUE, 
          values=.6, label="e", name="Z" ),
 mxMatrix(type="Full", nrow=1, ncol=(2+Nvar*2), 
          free=TRUE, values= 0, 
          label=c('meanfeno','meanfeno', 
                rep('mean',Nvar*2)), name="expMean"),
 mxMatrix(type ='Lower', nrow=Nvar , ncol=Nvar, 
          values=0.5, free=TRUE, name="CholCovW"),
 mxMatrix(type ='Lower', nrow=Nvar , ncol=Nvar,
          values=0.5, free=TRUE, name="CholCovB"),
    
 mxAlgebra(expression=CholCovW %*% t(CholCovW), 
           name="CovW" ),
 mxAlgebra(expression=CholCovB %*% t(CholCovB), 
           name="CovB" ),
 mxAlgebra(expression=CovB + CovW, 
           name="CovWplusB"),
 
 #Matrices A, C,E + compute variance components
 mxAlgebra(expression=X %*% t(X), name="A"),
 mxAlgebra(expression=Y %*% t(Y), name="C"),
 mxAlgebra(expression=Z %*% t(Z), name="E"),
                                                  
 #Declare a vector for the regression parameters
 mxMatrix(type="Full", nrow=Nvar, ncol=1, free=TRUE, 
          values= 0, 
          label=c("beta1","beta2","beta3",
                  "beta4","beta5"), 
          name="beta"),
                             
 #Algebra for expected variance/covariance matrix
 #in MZ twins
 mxAlgebra(
 expression=rbind(cbind(
                    A+C+E+t(beta)%*%CovWplusB%*%beta, 
                    A+C+t(beta)%*%CovB%*%beta, 
                    t(beta)%*%CovWplusB,
                    t(beta)%*% CovB),
                  cbind(
                    A+C+t(beta)%*%CovB%*%beta, 
                    A+C+E+t(beta)%*%CovWplusB%*%beta,
                    t(beta)%*%CovB,
                    t(beta)%*%CovWplusB),
                  cbind(
                    CovWplusB%*%beta,
                    CovB%*%beta, 
                    CovWplusB, CovB),
                  cbind(
                    CovB%*%beta, 
                    CovWplusB%*%beta, 
                    CovB, CovWplusB)), 
                  name="expCovMZ"),
                             
 #Algebra for expected variance/covariance matrix
 #in DZ twins
 mxAlgebra(expression=rbind(
                        cbind(
                            A+C+E+t(beta)%*%CovWplusB%*%beta, 
                            0.5%x%A+C+t(beta)%*%CovB%*%beta, 
                            t(beta)%*%CovWplusB,
                            t(beta)%*%CovB),
                        cbind(
                            0.5%x%A+C+t(beta)%*%CovB%*%beta, 
                            A+C+E+ t(beta)%*%CovWplusB%*%beta, 
                            t(beta)%*%CovB,
                            t(beta)%*%CovWplusB),
                        cbind(
                            CovWplusB%*%beta,
                            CovB%*%beta,
                            CovWplusB, CovB),
                        cbind(
                            CovB%*%beta,
                            CovWplusB%*%beta,
                            CovB, CovWplusB)), 
                        name="expCovDZ"),
    
 mxModel("MZ", 
 mxData( observed=mzData, type="raw" ),
 
 #Algebra for making the means a function 
 #of the definition variables
 mxFIMLObjective(covariance="twinACE.expCovMZ", 
                 means="twinACE.expMean", 
                 dimnames=names(mzData))),
                             
 mxModel("DZ", 
    mxData( observed=dzData, type="raw"),
    mxFIMLObjective(covariance="twinACE.expCovDZ", 
                    means="twinACE.expMean", 
                    dimnames=names(dzData))),
    mxAlgebra(expression=MZ.objective + DZ.objective, 
              name="twin" ),
    mxAlgebraObjective("twin")
)
