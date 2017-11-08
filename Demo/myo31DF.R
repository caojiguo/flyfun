# please install all these packages first 
install.packages('DiceKriging')
install.packages('dplyr')
install.packages('tidyr')
install.packages('matrixcalc')
install.packages('devtools')
install.packages('fda')
install.packages('broom') 

# this command will install our package flyfuns from github
install_github('YunlongNie/flyfuns')


# load all the R packages 
library(fda)
library(DiceKriging)
library(matrixcalc)
library(plyr)
library(dplyr)
library(flyfuns)
library(broom)



load('demo_data.Rdata') 
# the data file contains two list:  fdlist_20genes and Lchol 
# fdlist_20genes is a list with 20 elements and each element correponds to a gene's smoothing trajectory 
# Lchol is the precondioning matrix for each gene 

GeneName= matched_id("Myo31DF")$CG_ID
yindex = which(names(fdlist_20genes)==GeneName)

res_true= regfun_slos_fun(lambda=10, gamma=1e-3, xfdlist=fdlist_20genes, yfd=fdlist_20genes[[yindex]]  , Lcholy = Lchol[[GeneName]][1:24,1:24], time_obs=0:23,  yname= GeneName,xnames =names(fdlist_20genes),lambdaI = 1e4,maxiteration=500,
    d =4, K = 5,maxabs = 1e-6,verbose=FALSE,type="derivatives"
)


par(mfrow=c(1,nrow(res_true$estimated_fd) ))
for (j in 1:nrow(res_true$estimated_fd) )
{
  
  plot(res_true$estimated_fd$regfd[[j]],main=flyfuns::matched_id(res_true$estimated_fd$xname[j])$genesymbol[1],col=j,lwd=2)

}


