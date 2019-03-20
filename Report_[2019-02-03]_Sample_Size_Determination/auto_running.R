passed_arg = commandArgs(trailingOnly=TRUE)

cat("loading the data in  ",passed_arg,"\n")

x=scan(passed_arg)
x=(x-mean(x))/sd(x)
source("CI_function_verify.R")
#em_mod_mixtool=normalmixEM(x)
em_mod=emTrynormalkem(x,2)

# I_mixtool=MixtureNormInfoM_2comp_exp_verify(mean=em_mod_mixtool$mu,
#                                             sd=em_mod_mixtool$sigma,
#                                             lambda = em_mod_mixtool$lambda[1])
I=MixtureNormInfoM_2comp_exp_verify(mean=em_mod$Mu,
                                            sd=em_mod$Sig,
                                            lambda = em_mod$Pie[1])

pardev1=mixnormParitaldev(mean=em_mod$Mu,sd=em_mod$Sig,lambda = em_mod$Pie[1],q=0.1)
pardev9=mixnormParitaldev(mean=em_mod$Mu,sd=em_mod$Sig,lambda = em_mod$Pie[1],q=0.9)

varll=sqrt(t(pardev1) %*% solve(I) %*% (pardev1))
varhh=sqrt(t(pardev9) %*% solve(I) %*% (pardev9))
x1=quantile(x,0.1)
x9=quantile(x,0.9)

lowvar=highvar=rep(0,900)
for(i in 1:900)
{
  lowvar[i]=varll/x1/sqrt(i)
  highvar[i]=varhh/x9/sqrt(i)
}
png("a.png")
plot(1:900,highvar,type = "l",col="red",ylim=c(-1.5,3),xlab="sample size")
points(1:900,lowvar,type = "l",col="blue")
legend("topright",c("hign","low"),col=c("red","blue"),lty=c(1,1))
dev.off()

write.csv(data.frame(lowvar,highvar),file="out.csv")
