MixtureNormInfoM_2comp_verify=function(x,mean,sd,lambda)
{
  normDu=function(mean,sd,x)
  {
    z=(x-mean)/sd
    f=dnorm(x,mean,sd)
    return(z*f/sd)
  }
  normDuDu=function(mean,sd,x)
  {
    z=(x-mean)/sd
    f=dnorm(x,mean,sd)
    f*(z/sd)^2-1/sd^2*f
  }
  normDsigma=function(mean,sd,x)
  {
    z=(x-mean)/sd
    f=dnorm(x,mean,sd)
    -1/sd*f+f*z^2/sd
  }
  normDsigmaDsigma=function(mean,sd,x)
  {
    z=(x-mean)/sd
    f=dnorm(x,mean,sd)
    1/sd^2*f-normDsigma(mean,sd,x)/sd+normDsigma(mean,sd,x)*z^2/sd-3*f*z^2/sd^2
  }
  normDuDsiama=function(mean,sd,x)
  {
    z=(x-mean)/sd
    f=dnorm(x,mean,sd)
    -1/sd*normDu(mean,sd,x)+z^2/sd*normDu(mean,sd,x)-2*f*z/sd^2
  }
  
  #calculat the likelihood
  normMix_2comp_lik=function(mean,sd,lambda,x)
  {
    val=lambda*dnorm(x,mean = mean[1],sd = sd[1])+
      (1-lambda)*dnorm(x,mean = mean[2],sd = sd[2])
    val
  }
  
  h=matrix(ncol = 5, nrow = 5)
  f1=dnorm(x,mean[1],sd[1])
  f2=dnorm(x,mean[2],sd[2])
  Likeli=f1*lambda+(1-lambda)*f2
  #LogLikeli=log(Likeli)
  
  Ivec=rep(0,5)
  Ivec[1]=1/Likeli*(f1-f2)
  Ivec[2]=1/Likeli*lambda*normDu(mean[1],sd[1],x)
  Ivec[3]=1/Likeli*lambda*normDsigma(mean[1],sd[1],x)
  Ivec[4]=1/Likeli*(1-lambda)*normDu(mean[2],sd[2],x)
  Ivec[5]=1/Likeli*(1-lambda)*normDsigma(mean[2],sd[2],x)
  
  h=Ivec %*% t(Ivec)
  h=h*Likeli
  h[is.na(h)]=0
  h[is.infinite(h)]=10000
  h
}

#calculate the expected infoM
MixtureNormInfoM_2comp_exp_verify=function(mean,sd,lambda)
{
  
  Index=function(x,mean,sd,lambda,idx=c(1,1))
  {
    h=MixtureNormInfoM_2comp_verify(x,mean,sd,lambda)
    h[idx[1],idx[2]]
  }
  Index_vec=Vectorize(Index,vectorize.args = "x")
  Eh=matrix(ncol = 5,nrow = 5)
  for(i in 1:5)
  {
    for(j in i:5)
    {
      Eh[i,j]=integrate(Index_vec,-Inf,Inf,mean=mean,sd=sd,lambda=lambda,
                        idx=c(i,j))$value
      Eh[j,i]=Eh[i,j]
    }
  }
  Eh
}


mixnormQuantile=function(mean,sd,lambda,q=0.9)
{
  cdf_mnorm=function(x,mean,sd,lambda)
  {
    lambda*pnorm(x,mean[1],sd[1])+(1-lambda)*pnorm(x,mean[2],sd[2])
  }
  pdf_mnorm=function(x,mean,sd,lambda)
  {
    lambda*dnorm(x,mean[1],sd[1])+(1-lambda)*dnorm(x,mean[2],sd[2])
  }
  x1=sum(mean)-sum(sd)*10;x2=sum(mean)+sum(sd)*10
  f1=cdf_mnorm(x1,mean,sd,lambda)-q;f2=cdf_mnorm(x2,mean,sd,lambda)-q
  xx=(x1+x2)/2;fxx=cdf_mnorm(xx,mean,sd,lambda)-q
  while(abs(cdf_mnorm(xx,mean,sd,lambda)-q)>1e-6)
  {
    if(fxx==0)
    {
      return(xx)
    }
    if(fxx*f1<=0)
    {
      x2=xx
      xx=(xx+x1)/2
    } else if(fxx*f2<=0)
    {
      x1=xx
      xx=(xx+x2)/2
    
    } else#all failed, e
    {
      print("pick new star value")
      x1=x1*2;x2=x2*2;xx=(x1+x2)/2
    }
    f1=cdf_mnorm(x1,mean,sd,lambda)-q
    f2=cdf_mnorm(x2,mean,sd,lambda)-q
    fxx=cdf_mnorm(xx,mean,sd,lambda)-q
    
    #print(c(x1,f1,x2,f2,xx,fxx))
    #Sys.sleep(1)
  }
  return(xx)
}
#f=function(x,mean,sd,lambda,q)
#cdf_mnorm(mixnormQuantile(c(0,1),c(1,2),0.5,0.1),c(0,1),c(1,2),0.5)
###pass the unit text!
#calculate the partial dev
mixnormParitaldev=function(mean,sd,lambda,q=0.9)
{
  xq=mixnormQuantile(mean,sd,lambda,q)
  cdf_mnorm=function(x,mean,sd,lambda)
  {
    lambda*pnorm(x,mean[1],sd[1])+(1-lambda)*pnorm(x,mean[2],sd[2])
  }
  pdf_mnorm=function(x,mean,sd,lambda)
  {
    lambda*dnorm(x,mean[1],sd[1])+(1-lambda)*dnorm(x,mean[2],sd[2])
  }
  pardev=rep(0,5)
  pardev[1]=(pnorm(xq,mean[1],sd[1])-pnorm(xq,mean[2],sd[2]))/
    pdf_mnorm(xq,mean,sd,lambda)
  pardev[2]=lambda*dnorm(xq,mean[1],sd[1])/pdf_mnorm(xq,mean,sd,lambda)
  pardev[3]=lambda*dnorm(xq,mean[1],sd[1])/pdf_mnorm(xq,mean,sd,lambda)*
    (xq-mean[1])/sd[1]
  pardev[4]=(1-lambda)*dnorm(xq,mean[2],sd[2])/pdf_mnorm(xq,mean,sd,lambda)
  pardev[5]=(1-lambda)*dnorm(xq,mean[2],sd[2])/pdf_mnorm(xq,mean,sd,lambda)*
    (xq-mean[2])/sd[2]
  pardev
}

emTrynormalkem=function(x,k,itera=100)
{
  normalkem=function(x,k,itera=100)
  {
    n=length(x)
    Tm=matrix(ncol=k,nrow=n)
    Mu=rep(NA,k)
    Sig=rep(NA,k)
    Pie=rep(NA,k)
    #choosing starting value
    kmod=kmeans(x,centers = k)
    for(i in 1:k)
    {
      Mu[i]=mean(x[kmod$cluster==i])
      Sig[i]=sd(x[kmod$cluster==i])
      Pie[i]=sum(kmod$cluster==i)/length(x)
    }
    ss=1
    embool=F
    while(ss<=itera)
    {
      #construct T matrix
      for(i in 1:n)
      {
        #update each column of Tm
        Tm[i,]=Pie*dnorm(x[i],mean=Mu,sd=Sig)
        tmps=sum(Tm[i,])
        Tm[i,]=Tm[i,]/tmps
      }
      #update Pie
      for(j in 1:k)
      {
        #pie
        Pie[j]=sum(Tm[,j])/n
        Mu[j]=sum(x*Tm[,j])/sum(Tm[,j])
        Sig[j]=sqrt(sum(Tm[,j]*(x-Mu[j])^2)/sum(Tm[,j]))
      }
      ss=ss+1
    }
    embool=!all(is.na(Sig))
    list(Mu=Mu,Sig=Sig,Pie=Pie,embool=embool)
  }
  
  emlokikeli=function(x,emmod)
  {
    #calculate the incomplete likelihood
    n=length(x)
    k=length(emmod$Pie)
    tmp=rep(NA,n)
    for(i in 1:n)
    {
      tmp[i]=sum(emmod$Pie*dnorm(x[i],mean=emmod$Mu,sd=emmod$Sig))
    }
    sum(log(tmp))
  }
  
  emAIC=function(x,emmod)
  {
    bic=emlokikeli(x,emmod)
    k=length(emmod$Mu)
    n=length(x)
    2*k-2*bic
  }
  
  emBIC=function(x,emmod)
  {
    bic=emlokikeli(x,emmod)
    k=length(emmod$Mu)
    n=length(x)
    log(n)*(2*k+k-1)-2*bic
  }
  
  succeed=F
  ll=1
  while(!succeed)
  {
    emmod=normalkem(x,k,itera)
    if(emmod$embool)
    {
      succeed=T
    }
    else
    {
      Sys.sleep(2)
      set.seed(Sys.time())
      print("NA occurs, try once again!")
    }
    ll=ll+1
    if(ll>100)
    {
      print("too many tials, check the x")
      return(emmod)
      break
    }
  }
  emmod$loglik=emlokikeli(x,emmod)
  emmod$AIC=emAIC(x,emmod)
  emmod$BIC=emBIC(x,emmod)
  emmod
}