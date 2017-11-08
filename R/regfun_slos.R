#' This function compute the regulation functions via slos penalty
#' @param xfdlist a list with each element being a fd object, for all the predictors
#' @param yfd a fd object for the response 
#' @param time_obs  evaluated time points 
#' @param xnames  matches the xfdlist
#' @param yname  reponse name 
#' @param d the order 
#' @param K the number of inside knots 
#' @param lambda the fscad penalty parameter 
#' @param gamma the smoothing parameter
#' @export
#' @import tidyr
#' @import DiceKriging
#' @import matrixcalc
#' @import dplyr 
#' @importFrom plyr alply
#' @examples
#' \dontrun{
#' yfd = fdlist[[14]]
#' xfdlist  = fdlist
#' time_obs = (eval_times = timepoints$realtime)
#' xnames = names(xfdlist)
#' res = regfun_slos(xfdlist,yfd,time_obs,yname=xnames[14],xnames)
#' }

     

regfun_slos= function(xfdlist, yfd,   Lcholy , time_obs=seq(0,23,len=24), yname, xnames,
  lambda = 1e-3,gamma=1e-3,maxiteration=500,lambdaI=1e5,
  d =5, K = 5,maxabs = 1e-6,verbose=TRUE,type="derivatives"
     )
{

if (type=="derivatives") { y = eval.fd(time_obs,yfd,Lfdobj=1)%>%as.vector} else {y = eval.fd(time_obs,yfd,Lfdobj)}
intery = mean(y)
y = y - intery


#### compute the trajactories for the x list 
if (missing(xnames)) xnames=1:length(xfdlist)

xbasis = data_frame(xname =xnames,xfd = xfdlist)%>%group_by(xname)%>%do(xfd=.$xfd[[1]],xbasis=create.bspline.basis(rangeval = range(c(eval.fd(c(time_obs,seq(0,max(time_obs),len=1001)),.$xfd[[1]]))),nbasis=d+K,norder=d))

(xnames=xbasis$xname)

x = xbasis%>%do(x = eval.fd(time_obs,.$xfd))%>%first%>%do.call(cbind,.)
colnames(x) = xnames

Phi_matrix_big0 =  (xbasis%>%group_by(xname)%>%
      do(Phi_matrix= t(eval.basis(as.vector(eval.fd(time_obs,.$xfd[[1]])),.$xbasis[[1]],Lfdobj=0)))%>%broom::tidy(Phi_matrix)%>%data.frame%>%dplyr::select(-matches('name'))%>%as.matrix)
Phi_matrix_big = Phi_matrix_big0%*%t(Lcholy)


# wtide_fun = function(xfd,xbasis,time_obs){
# rangeX = range(c(eval.fd(c(time_obs,seq(0,max(time_obs),len=1001)),xfd)))
# weightm = eval.basis(seq(min(rangeX),max(rangeX),length=10000),xbasis)%>%colSums%>%as.numeric
# wtilde = matrix(rep(weightm, each=length(weightm)),ncol=length(weightm))*matrix(rep(weightm, each=length(weightm)),ncol=length(weightm),byrow=TRUE)
# return(wtilde)  
# }

wtide_fun = function(xfd,xbasis,time_obs){
Xtra = eval.fd(time_obs,xfd)%>%as.numeric
weightm = eval.basis(Xtra,xbasis)%>%colSums%>%as.numeric
wtilde = matrix(rep(weightm, each=length(weightm)),ncol=length(weightm))*matrix(rep(weightm, each=length(weightm)),ncol=length(weightm),byrow=TRUE)
return(wtilde)  
}



Wtilde = xbasis%>%group_by(xname)%>%do(wtide = {wtide_fun(.$xfd[[1]],.$xbasis[[1]],time_obs)})%>%dplyr::select(wtide)%>%ungroup%>%first%>%do.call(bdiag,.)

############################################
########   roughness penalty on time    ############
############################################

# V_all

# # part 1  use fdlist each element in fdlist the smoothed fd for each gene 
inside_function1 = function(t,xfd0,xbasis0)
     {
      weight = c(1,rep(c(4,2),(length(t)-3)/2),4,1)*(length(time_obs)-1)/(length(t)-1)/3
      
      (x_t = as.numeric(eval.fd(t,xfd0))) # the gene expression at time t
      first_de   = xfd0%>%eval.fd(t,.)%>%as.vector%>%eval.basis(.,xbasis0,Lfdobj=1)
      second_de = xfd0%>%eval.fd(t,.)%>%as.vector%>%eval.basis(.,xbasis0,Lfdobj=2)
    
      first_x = eval.fd(t,xfd0,Lfdobj=1)^2
      second_x = eval.fd(t,xfd0,Lfdobj=2)

      mat1 = second_de*matrix(rep(first_x,ncol(second_de)),ncol=ncol(second_de))

      mat2 =  first_de*matrix(rep(second_x,ncol(first_de)),ncol=ncol(first_de))

      mat3  = mat1 + mat2
      temp =  alply(mat3,1,function(x)
      {
        tmat = matrix(rep(x,length(x)),ncol=length(x))
        tmat*t(tmat)
        })

      tempt = mapply(function(x,y) x*y,as.list(weight),temp,SIMPLIFY = FALSE)
      return(Reduce('+',tempt))
      }

V_all = xbasis%>%do(V_mat={inside_function1(t=seq(0,max(time_obs),length=1001),xfd0 = .$xfd,xbasis0=.$xbasis)})%>%first%>%do.call(bdiag,.)


wtemp  = function(xbasis0){
xbasis0%>%knots(,interior=FALSE)%>%unique%>%
data.frame(knot =.)%>%mutate(knotlead= lead(knot))%>%dplyr::filter(!is.na(knotlead))%>%
rowwise()%>%do(temp = eval.penalty(xbasis0,int2Lfd(0),rng=c(.$knot,.$knotlead)))
}





Wtemp = xbasis%>%group_by(xname)%>%do(w_temp  = {wtemp(.$xbasis[[1]])})
W_temp = Wtemp$w_temp

(nspline = d+K)
set.seed(2015)
beta0=rnorm(nspline*length(xfdlist))
break_while=FALSE
beta = 2*beta0
beta_new = beta0

if (verbose)  cat('\n################\nComputing when fscad penalty parameter is :', lambda ,'and smoothing parameter is: ', gamma)

nreg = length(xfdlist)
i=1
while(max(abs(beta_new-beta)) >maxabs)
 {

index_nonzero = rep(TRUE,length(beta0))
beta= beta_new
#library(DiceKriging) 
W0=list()
for (p in 1:nreg)
{
betap = beta[1:nspline+(p-1)*nspline]
res = matrix(0,nspline,nspline);
for (j in 2:(length(W_temp[[p]]$temp)+1))
{
          
          rangetemp = range(x[,p])
          (betaNorm = sqrt(as.numeric(betap^T%*%W_temp[[p]]$temp[[j-1]]%*%betap)))
          (temp1 =  SCAD.derivative(sqrt((K+1)/(rangetemp[2]-rangetemp[1]))*betaNorm,lambda=lambda))
          if (temp1!=0) # if temp1 is too small makes it to be zero
          {

              if (betaNorm < 1e-6)  # absTol cannot be too small then temp will be way too large absTol cannot be too small compared with lambda and (rangetemp[2]-rangetemp[1])/(K+1)
          {
         
          index_nonzero[(j-1):(j-1+d-1) + (p-1)*nspline]=FALSE; 
          }  else 
          {
             
             temp2 = betaNorm*sqrt((rangetemp[2]-rangetemp[1])/(K+1))
             tr = as.numeric((temp1/temp2))
             temp = tr*W_temp[[p]]$temp[[j-1]]
             res = res + temp
               
          }
        }
}
W0[[p]] =0.5*res
}
W0_all = do.call(bdiag,W0)
nonezero= index_nonzero;
if (sum(nonezero)==0) {beta_new = rep(0,length(beta));cat('lambda is too large! \n');break}
Phi_nonzero_bigg = Phi_matrix_big[nonezero,]
G_nonzero = crossprod(t(Phi_nonzero_bigg),t(Phi_nonzero_bigg))
V_nonzero = V_all[nonezero,nonezero]
W0_nonzero = W0_all[nonezero,nonezero]
right_nonzero =  crossprod(t(Phi_nonzero_bigg),Lcholy%*%matrix(y))
mat= (G_nonzero + length(y)*gamma*V_nonzero + length(y)*W0_nonzero+length(y)*lambdaI*Wtilde[nonezero,nonezero])
beta_new_nonzero =  try(solve(mat,right_nonzero))

if (class(beta_new_nonzero) == "try-error"|any(is.na(beta_new_nonzero))) {beta_new=NA;break_while=TRUE;cat('\nmat singular!loop break!\n');break}
beta_new = rep(0,length(beta))
beta_new[nonezero] = as.numeric(beta_new_nonzero)
i=i+1
if (max(abs(beta_new-beta))==Inf) {cat('beta is inf!\n'); break_while=TRUE;break}
if (i%%100==0&verbose==TRUE) sprintf("\n#\n%sth iteration: the difference between two consective interations is %s; nonezero terms left: %s out of %s", i,format(max(abs(beta_new-beta)),digits=3),sum(nonezero),length(beta))%>%cat
if (i > maxiteration) 
  {
    sprintf('\nWarning: alogrithm did not coverage within %s iteration\n', maxiteration)%>%cat ;break
  }
}

if (break_while) break

beta_new_nonzero[abs(beta_new_nonzero)<=maxabs]=0
haty_transform = as.numeric(t(beta_new_nonzero)%*%as.matrix(Phi_nonzero_bigg))
y_transform = Lcholy%*%matrix(y)
mse_transform = mean((y_transform -haty_transform)^2)
(trace= matrix.trace(as.matrix(crossprod(as.matrix(Phi_nonzero_bigg),solve(mat,as.matrix(Phi_nonzero_bigg))))))

AIC = length(y)*log(mse_transform) + 2*trace
AICc = AIC + 2*trace*(trace+1)/(length(y) - trace- 1)
BIC = length(y)*log(mse_transform) + trace*log(length(y))
gcv = mse_transform*length(y)/((length(y)-  trace)/length(y))^2

haty = as.numeric(beta_new%*%as.matrix(Phi_matrix_big0))
mse = mean((haty-y)^2)


if (verbose==TRUE)  
  {
    sprintf('\n#\nThe different between the final two interations are %s\nThe total number of nonzero parameters are %s out of %s \nThe total number of interation is %s\n', format(max(abs(beta_new-beta)),digits=3),sum(nonezero),length(beta_new),i)%>%cat

}

res_base  = data_frame(xname = xnames,beta = matrix(beta_new,ncol=length(xnames))%>%data.frame%>%as.list)%>%left_join(xbasis,by="xname")

res_final = res_base%>%group_by(xname)%>%do(regfd = fd(coef=.$beta[[1]], basisobj=.$xbasis[[1]]))%>%left_join(res_base,by="xname")
names(res_final) = c('xname',"regfd","coefreg","xtrafd","regbasisfd")

coefm = matrix(beta_new,ncol=length(xnames));colnames(coefm) = xnames

return(list(coef = coefm, estimated_fd =res_final, x=x, xname= xnames, transformy =y_transform , hatytransform = haty_transform, y=y,haty = haty, intery=intery, time=time_obs,trace=trace, AIC=AIC, AICc=AICc,BIC=BIC,gcv = gcv, mse = mse ,mse_transform = mse_transform, targent_gene = yname,type=type,lambda=lambda,gamma=gamma,Lcholy = Lcholy))

}



