library(zoo)
library(scales)
library(dplyr)
library(reshape2)
library(lubridate)
library(fda)
library(ggpubr)
library(objfunction)
library(objectfunctionsbooster)
library(parallel)
library(ggsci)

#search: h = 
#----------- Boundary Kernel and Smooth function----------------
get_Boundary_Weight <- function(t, N, h){
  get_muk <- function(l,p){
    f <- function(u) 0.75*u^l*(1-u^2)
    muk = integrate(f,p,1)
    return(muk$value)
  }
  
  get_K <- function(t){
    if(0.75*(1-t^2)>0){
      return(0.75*(1-t^2))
    }else{return(0)}
  }
  
  get_B <- function(pt,t){
    K = get_K(t)
    if(pt <= -1){
      return(K)
    }else if(pt <= 0){
      muk0 = get_muk(0,pt)
      muk1 = get_muk(1,pt)
      muk2 = get_muk(2,pt)
      B = (muk2-muk1*t) / (muk0*muk2-muk1^2)*K
      return(B)
    }else{return(NA)}
  }
  
  Weight = c()
  
  weight_cal <- function(x, pt, choice)
  {
    if(choice=='positive')
      return(get_B(pt,(t-x)/h))
    else
      return(get_B(pt,-(t-x)/h))
  }
  if(t > (N+1)/2){
    pt = (t-N) / h
    Weight = sapply(1:N, weight_cal, pt, choice='positive')
  }else{
    pt = (1-t) / h
    Weight = sapply(1:N, weight_cal, pt, choice='negative')
  }
  
  return(Weight)
}

## One-point Smooth
Weighted.Ave <- function(x,h){
  N = length(x)
  
  Weighted.m <- function(ii, x)
  {
    weighted.mean(x,get_Boundary_Weight(ii, N, h), na.rm=T)
  }
  return(sapply(1:N, Weighted.m, x))
}
## Smooth function ; 
mav <- function(x,h){
  dx=diff(x)
  dy = Weighted.Ave(dx,h)
  # dy = Local.Linear(dx,h)
  y = cumsum(c(x[1],dy))
  y = as.double(y)
  return(y)
}

mav_increa <- function(x,h){
  dx=diff(x)
  dy = sapply(Weighted.Ave(dx,h), max, 0)
  # dy = Local.Linear(dx,h)
  y = cumsum(c(x[1],dy))
  y = as.double(y)
  return(y)
}

vSEIDR_Data_reinfect <- function(opendir, h, Omicronst){
  A1 = read.csv(opendir, stringsAsFactors = F) 
  L_dat=nrow(A1)
  
  A1$mu_r=0
  Omi_st=which(A1$date==Omicronst)
  A1$mu_r[(Omi_st-1):L_dat]=mu_r
  for(ss in Omi_st:L_dat)
  {
    A1$I[ss]=max(A1$I[ss-1]+A1$N[ss]-A1$N[ss-1]-A1$I[ss-1]/14-(A1$R.d[ss]-A1$R.d[ss-1]), 0)
    A1$R.r[ss]=max(A1$R.r[ss-1]+A1$I[ss-1]/14-A1$mu_r[ss-1]*A1$R.r[ss-1], 0)
  }
  
  # A1$R.r[Omi_st:L_dat]=sapply(mav(A1$R.r[Omi_st:L_dat], h), max, 0)
  A1$R = A1$R.d + A1$R.r
  # A1$I[Omi_st:L_dat] = sapply(mav(A1$I[Omi_st:L_dat], h), max, 0)
  # A1$N[A1$N<A1$I+A1$R]=(A1$I+A1$R)[A1$N<A1$I+A1$R]
  A1$date=as.character(A1$date)
  
  print(all((A1$N-A1$I-A1$R)[Omi_st:L_dat]>=0))
  # apply(A1[,c('infected','dead','recovered','people_vaccinated', 'people_partly_vaccinated', 'people_fully_vaccinated',
  #              'N', 'R.d', 'R.r', 'G1', 'G2')], 2, range, na.rm=T)
  return(A1)
}

provinces=read.csv("~/Documents/COVID-19/data/finalnames.csv",stringsAsFactors = F)
# read.csv(paste0("~/Documents/COVID-19/data/population.csv"),stringsAsFactors = F)
countries=setdiff(sort(unique(provinces$provinces.en)),'NonHB')
Omicron_period=read.csv('~/Documents/COVID-19/data/Omicron_bd.csv', stringsAsFactors = F)

mu_r=1/(16*30)
for(j1 in c(4, 8, 13, 18, 25, 26, 27))
{
  set.seed(1234)
  
  country=countries[j1] #country='US'
  print(country)
  opendir = paste0('~/Documents/COVID-19/data/real_data_0410/',country,'.csv')
  # st=as.Date(provinces$start[provinces$provinces.en==country])
  
  dat = vSEIDR_Data_reinfect(opendir, h = 15, as.Date(Omicron_period$Omicron_b[Omicron_period$province==country]))
  dat=dat[,c('province', 'population', 'date','infected','dead','recovered', 
             'people_vaccinated', 'people_partly_vaccinated', 'people_fully_vaccinated',
             'N', 'R.d', 'R.r', 'R', 'I', 'G1', 'G2', 'G3', 'people_boosters', 't', 'mu_r')]

  write.csv(dat, paste0('~/Documents/COVID-19/data/real_data_0410/',country,'_reinfect.csv'))
}

#----------------------
kernel_reg <- function(y, A4, nrows, w, yname, xname)
{
  W = c(get_Boundary_Weight(y, nrows, w/2),rep(0,nrow(A4)-nrows))
  return(sum(A4[,yname] *A4[,xname]*W,na.rm=T)/sum(W*(A4[,xname])^2,na.rm=T))
}

Ea0.dist = function(delta, dat, theta, alpha, r.beta, M, w, B = 300, t_nV, nrow_dat){
  Ea = rep(NA, nrow_dat)
  
  Ea[1] = delta*dat$Ep.hat[1]; Ra=rep(0, nrow_dat)
  for(t in 2 : (nrow_dat-1)){
    Ea[t] = (1 - theta) * dat$Gt[t-1] + (1 - dat$gamma.r[t-1]) * Ea[t-1]
    Ra[t] = Ra[t-1] + dat$gamma.r[t-1]*Ea[t-1] - dat$mu_r[t-1]*Ra[t-1]
  }
  
  t=t+1; Ra[t]=Ra[t-1] + dat$gamma.r[t-1]*Ea[t-1] - dat$mu_r[t-1]*Ra[t-1]
  
  dat$Ea = Ea; dat$Ra=Ra
  dat$S = M - dat$I - dat$R - dat$Ep.hat - dat$Ea - dat$Ra
  dat$S[!is.na(dat$S) & dat$S<0]=0
  
  dat$y.temp = dat$Gt
  dat$x.temp = ((dat$I + dat$Ea)/r.beta + dat$Ep.hat)*dat$S/M
  beta.temp = rep(NA, nrow_dat)
  
  beta.temp[2:(nrow_dat-1)]=sapply(1:(nrow_dat-2), kernel_reg, dat, nrow_dat-2, w, 'y.temp', 'x.temp')
  
  dat$beta=beta.temp
  dat$beta[dat$beta<0] = 0.001
  dat.2 = dat[!is.na(dat$beta) & dat$t <= t_nV,]
  
  cal_N_boot <- function(x, t, y, alpha, beta, gamma.d, gamma.r, theta, M, r.beta)
  {
    dat.3 = sim_vsveipdr(t, y, alpha, beta, gamma.d, gamma.r,
                         theta, M, 0, 0, 0, 0, 0, 0, r.beta)
    return(dat.3[, 14])
  }
  
  N.matr = sapply(1:B, cal_N_boot, t=nrow(dat.2)-1, y = c(0,0,as.numeric(dat.2[1, c('Ea', 'Ep.hat', 'I', 'R.d', 'R.r', 'Ra')])),
                  alpha = rep(alpha, nrow(dat.2)), beta = dat.2$beta, gamma.d = dat.2$gamma.d, 
                  gamma.r = dat.2$gamma.r, theta = theta, M = M, r.beta)
  N.boot = apply(N.matr, 1, mean)
  
  distan=abs(1-N.boot/dat.2$N)
  distan=distan[-union(1, which(dat.2$N==0))]
  return(mean(distan,na.rm=T))
}

vSEPIDR_Coef_mmd_Gs = function(A4, w, r, D, M, alpha, t_nV){
  nrow_dat=nrow(A4)
  gamma.d.temp = rep(NA,nrow_dat)
  gamma.r.temp = rep(NA,nrow_dat)
  beta.temp = rep(NA,nrow_dat)
  
  A4$R.d.diff2 = c(diff(A4$R.d),NA)
  A4$R.r.diff2 = c(diff(A4$R.r),NA)+A4$R_S
  A4$I[A4$I<=0]=0.1
  
  gamma.d.temp[2:(nrow_dat)]=sapply(1:(nrow_dat-1), kernel_reg, A4, nrow_dat-1, w, 'R.d.diff2', 'I')
  gamma.r.temp[2:(nrow_dat)]=sapply(1:(nrow_dat-1), kernel_reg, A4, nrow_dat-1, w, 'R.r.diff2', 'I')
  
  gamma.r.temp[is.nan(gamma.r.temp) | gamma.r.temp<0] = 0.01; gamma.d.temp[is.nan(gamma.d.temp) | gamma.d.temp<0] = 0.001
  A4$gamma.r = gamma.r.temp; A4$gamma.d = gamma.d.temp
  A4$gamma.r[1] = A4$gamma.r[2]; A4$gamma.d[1] = A4$gamma.d[2]
  
  A4$Ep.hat = c(mav(c(diff(A4$N) / alpha), h = 15),NA)
  A4$Ep.hat[A4$Ep.hat<0]=0
  
  A4$Gt=c((A4$Ep.hat[-c(1, nrow_dat)] - (1 - alpha) * A4$Ep.hat[-c(nrow_dat-1, nrow_dat)])/A4$theta[-c(nrow_dat-1, nrow_dat)], NA, NA)
  A4$Gt[A4$Gt <= 0] = 0.1
  
  Ea.hat = Ra.hat = rep(NA, nrow_dat)
  Ra.hat[1]=0  
  
  # de.ops=c(seq(1, 9.5, 0.5)/10, seq(1,10,0.5))
  de.ops=c(seq(1,10,0.5))
  if(A4$Ep.hat[1]<=0)
  {
    delta = NULL
    Ea.hat[1]=0
  }
  else
  {
    de.dist=sapply(de.ops, Ea0.dist, dat = A4, theta = A4$theta[1], alpha = alpha, r.beta= r, M = M, w = w, t_nV=t_nV, nrow_dat=nrow_dat)
    delta = de.ops[which.min(de.dist)]
    Ea.hat[1] = delta*A4$Ep.hat[1]
  }
  
  for(t in 2 : (nrow_dat-1)){
    Ea.hat[t] = (1 - A4$theta[t-1]) * A4$Gt[t-1] + (1 - A4$gamma.r[t-1]) * Ea.hat[t-1]
    Ra.hat[t] = min(Ra.hat[t-1] + A4$gamma.r[t-1]*Ea.hat[t-1]-A4$mu_r[t-1]*Ra.hat[t-1], 
                    M - A4$I[t] - A4$R[t] - A4$Ep.hat[t] - Ea.hat[t])
  }
  t=t+1; Ra.hat[t]=min(Ra.hat[t-1] + A4$gamma.r[t-1]*Ea.hat[t-1]-A4$mu_r[t-1]*Ra.hat[t-1], 
                       M - A4$I[t] - A4$R[t] - A4$Ep.hat[t] - Ea.hat[t]); A4$Ra.hat = Ra.hat
  A4$Ea.hat = Ea.hat
  A4$S.hat = M - A4$I - A4$R - A4$Ep.hat - A4$Ea.hat-A4$Ra.hat
  A4$S.hat[!is.na(A4$S.hat)&A4$S.hat<0]=0
  A4$S.tilde=A4$S.hat/M
  
  A4$y.temp = A4$Gt
  A4$allinfected=(A4$I + A4$Ea.hat)/r + A4$Ep.hat
  A4$x.temp = A4$allinfected*A4$S.hat/M
  
  beta.temp[2:(nrow_dat-1)]=sapply(1:(nrow_dat-2), kernel_reg, A4, nrow_dat-2, w, 'y.temp', 'x.temp')
  beta.temp[beta.temp<=0]=0.001
  
  A4$alpha = alpha
  
  # result
  coef.one = data.frame(gamma = gamma.r.temp+gamma.d.temp,
                        gamma.d = gamma.d.temp,
                        gamma.r = gamma.r.temp,
                        beta = beta.temp,
                        Rt = beta.temp *((1/alpha+D/r)*A4$theta+(1-A4$theta)* D/r)*A4$S.tilde)
  
  if('province' %in%colnames(A4)) A4 = A4[,c('t','province', 'population', 'date', 'infected', 'dead', 'recovered',
                                             'people_vaccinated', 'people_partly_vaccinated', 'people_fully_vaccinated',
                                             'N', 'R.d', 'R.r', 'R', 'I', 'G1', 'G2', 'G3', 'people_boosters', 'R_S', 'mu_r',  
                                             'alpha','theta','S.hat','Ra.hat','Ep.hat','Ea.hat','x.temp','y.temp','allinfected')]
  else  A4 = A4[,c('t','N', 'R.d', 'R.r', 'R', 'I', 'G1', 'G2', 'G3', 'R_S', 'mu_r', 'alpha', 'theta',
                   'S.hat','Ra.hat','Ep.hat','Ea.hat','x.temp','y.temp','allinfected')]
  coef.one = cbind(A4,coef.one)
  return(list(coef = coef.one, op.par = delta))
}

path_table='~/Documents/COVID-19/data/real_data_V1_3/'
path_plot='~/Documents/COVID-19/data/plot_real_data_V1_3/'

#---------------set parameters----------------
iter=300
D=14
r.beta = r = 5
t_b=15
t_a=15
t_V=50
t_choose=30
criteria=c('Ep', 'I', 'N')

country_knots=read.csv('~/Documents/COVID-19/data/nknots_0410.csv',stringsAsFactors = F)
time_variants=read.csv('~/Documents/COVID-19/data/COVID-19-Delta.csv', stringsAsFactors = F)
time_variants$Earliest.report=as.character(as.Date(time_variants$Earliest.report))
time_variants$Date.local.transmission=as.character(as.Date(time_variants$Date.local.transmission))
time_variants$province[time_variants$province=="United Kingdom"]='UK'
time_variants$province[time_variants$province=="United States"]='US'
Delta_period=read.csv('~/Documents/COVID-19/data/Delta_dominate_new.csv', stringsAsFactors = F)

interval_criteria<- function(y, inter_bs, inter_l, t_V)
{
  domain=rep(inter_bs, each=inter_l+1)+rep(0:(inter_l), length(inter_bs))
  domains=c(domain, domain+t_V, domain+t_V*2)
  # c('dist.Ep.2', 'dist.Ep.8', 'dist.Ep.14', 'dist.Ep.20',
  #   'dist.I.2', 'dist.I.8', 'dist.I.14', 'dist.I.20',
  #   'dist.N.2', 'dist.N.8', 'dist.N.14', 'dist.N.20')
  return(sqrt(apply(matrix(as.numeric(y[domains]), nrow=inter_l+1), 2, mean, na.rm=T)))
}

balance_weight <- function(country)
{
  if(country == c('Germany'))
    return(0.966)
  else if(country == c('Italy'))
    return(0.963)
  else if(country == c('Turkey'))
    return(0.977)
  else if(country == c('UK'))
    return(0.97)
  else if(country == c('US'))
    return(0.962)
  else
    return(0.965)
}

coefs_pre_Omi=c()
for(j1 in c(4, 8, 13, 18, 25, 26, 27))
{
  country=countries[j1] #country='US'
  print(country)
  dat=read.csv(paste0('~/Documents/COVID-19/data/real_data_0410/',country,'_reinfect.csv'), stringsAsFactors = F)
  dat=dat[,c('province', 'population', 'date','infected','dead','recovered', 
             'people_vaccinated', 'people_partly_vaccinated', 'people_fully_vaccinated',
             'N', 'R.d', 'R.r', 'R', 'I', 'G1', 'G2', 't', 'G3', 'people_boosters', 'mu_r')]
  nrow_dat=nrow(dat)
  t_nV=min(which(dat$G1>0))-1
  t_booster=min(which(dat$G3>0))-1
  t_delta=which(as.character(as.Date(dat$date))==as.character(
    time_variants$Earliest.report[time_variants$province==country&time_variants$Variant=='B.1.617.2']))-1
  t_delta_d=which(as.character(as.Date(dat$date))==as.character(as.Date(
    Delta_period$dominate[Delta_period$province==country])))-1
  t_omicron=which(as.character(as.Date(dat$date))==as.character(as.Date(
    Omicron_period$Omicron_b[Omicron_period$province==country])))-1
  t_omicron_d=which(as.character(as.Date(dat$date))==as.character(as.Date(
    Omicron_period$Omicron_d[Omicron_period$province==country])))-1
  M=dat$population[1] 
  w0=provinces$w[provinces$provinces.en==country]
  
  dat[1:t_nV, c('G1', 'G2', 'G3')]=0
  date.c=dat$date[t_nV]
  print(date.c)
  # dat$theta=c(rep(0.8, t_delta), rev(seq(from=0.6, to=0.8, length.out=t_omicron-t_delta+1)[-1]),
  #             rev(seq(from=0.15, to=0.6, length.out=t_omicron_d-t_omicron+1)[-1]),
  #             rep(0.15, nrow_dat-t_omicron_d))
  dat$theta=c(rep(0.8, t_delta), rev(seq(from=0.6, to=0.8, length.out=t_omicron-t_delta+1)[-1]),
              rev(seq(from=0.1, to=0.6, length.out=nrow_dat-t_omicron)))
  dat$R_S=c(diff(dat$N-dat$I-dat$R), NA)
  
  estimates_or=read.csv(paste0(path_table,country,'_estimates_Omicron_dominate.csv'), stringsAsFactors = F)
  als=estimates_or$alpha[1]
  
  ## result of vSIADR model, local linear
  result=vSEPIDR_Coef_mmd_Gs(dat, w0, r, D, M, als, t_nV)
  coef_SEPIDR=result$coef
  colnames(coef_SEPIDR)[which(colnames(coef_SEPIDR)=="Ep.hat")]="E.hat"
  
  coef.es=coef_SEPIDR
  ## estimates of S(t)+V1(t)+V2(t)+V3(t) (The estimate of S(t) from vSIADR model is actually the estimate of S(t)+V1(t)+V2(t)+V3(t) in vSVIADR model)
  coef.es$SV.hat=coef.es$S.hat
  ## estimates of V1(t): the cumulative number of people partly vaccinated at day t, G1(t)-G2(t). True V1(t) is overestimated
  coef.es$V1.hat=coef.es$G1-coef.es$G2
  ## estimates of V2(t): the cumulative number of people fully vaccinated at day t, G2(t)-G3(t). True V2(t) is overestimated
  coef.es$V2.hat=coef.es$G2-coef.es$G3
  ## estimates of V3(t): the cumulative number of people receiving boosters at day t, G3(t). True V3(t) is overestimated
  coef.es$V3.hat=coef.es$G3
  ## estimates of S(t): \hat{S(t)+V1(t)+V2(t)+V3(t)}-G1(t), true S(t) is underestimated
  coef.es$S.hat=coef.es$SV.hat-coef.es$G1
  coef.es$dN=c(diff(coef.es$N), NA)
  coef.es$dG1=c(diff(coef.es$G1), NA)
  coef.es$dG2=c(diff(coef.es$G2), NA)
  coef.es$dG3=c(diff(coef.es$G3), NA)
  coef.es[, paste0('scale_',c('dN', 'dG1', 'dG2', 'dG3', 'G1', 'G2', 'G3', 'N', 'I'))]=
    coef.es[, c('dN', 'dG1', 'dG2', 'dG3', 'G1', 'G2', 'G3', 'N', 'I')]/M
  
  gamma.d.temp = rep(NA,nrow_dat)
  gamma.r.temp = rep(NA,nrow_dat)
  R.d.diff2 = c(diff(coef.es$R.d),NA)
  R.r.diff2 = c(diff(coef.es$R.r),NA)+coef.es$R_S
  posi_I = coef.es$I
  posi_I[posi_I<=0]=0.1
  A4=data.frame(R.d.diff2 = R.d.diff2, R.r.diff2 = R.r.diff2, I=posi_I)
  gamma.d.temp[2:(nrow_dat)]=sapply(1:(nrow_dat-1), kernel_reg, A4, nrow_dat-1, 11, 'R.d.diff2', 'I')
  gamma.r.temp[2:(nrow_dat)]=sapply(1:(nrow_dat-1), kernel_reg, A4, nrow_dat-1, 11, 'R.r.diff2', 'I')
  gamma.r.temp[is.nan(gamma.r.temp) | gamma.r.temp<0] = 0.01; gamma.d.temp[is.nan(gamma.d.temp) | gamma.d.temp<0] = 0.001
  coef.es$gamma.r.11 = gamma.r.temp; coef.es$gamma.d.11 = gamma.d.temp
  
  coef.es=left_join(coef.es, estimates_or[, c('date', 'kappa', 'varphikappa', 'omegakappa', 'rho', 'omega')])
  coef.es$kappa[(t_omicron+1):nrow_dat]=coef.es$kappa[t_omicron]
  coef.es$varphikappa[(t_omicron+1):nrow_dat]=coef.es$varphikappa[t_omicron]
  coef.es$omegakappa[(t_omicron+1):nrow_dat]=coef.es$omegakappa[t_omicron]
  coef.es$rho[(t_omicron+1):nrow_dat]=coef.es$rho[t_omicron]
  coef.es$omega[(t_omicron+1):nrow_dat]=coef.es$omega[t_omicron]
  coef.es$mu_1=c(rep(1/(30*2), t_omicron-1), rep(1/(8*7), nrow_dat-t_omicron+1))
  coef.es$mu_2=c(rep(1/(30*8), t_omicron-1), rep(1/(30*3), nrow_dat-t_omicron+1))
  coef.es$mu_3=c(rep(1/(30*10), t_omicron-1), rep(1/(20*7), nrow_dat-t_omicron+1))
  
  S.modify=coef.es$S.hat
  V1.modify=coef.es$V1.hat
  V2.modify=coef.es$V2.hat
  V3.modify=coef.es$V3.hat
  
  S_infected=V1_infected=V2_infected=V3_infected=rep(0, nrow_dat)
  for(j in t_nV:(nrow_dat-1))
  {
    #transfer in [j, j+1]
    V3_infected[j]=max(min(min(coef.es$y.temp[j]/(S.modify[j]+coef.es$varphikappa[j]*V1.modify[j]+coef.es$kappa[j]*V2.modify[j]+
                                                    coef.es$omegakappa[j]*V3.modify[j])*coef.es$omegakappa[j]*V3.modify[j],
                               coef.es$SV.hat[j]-coef.es$SV.hat[j+1]+coef.es$mu_r[j]*coef.es$Ra.hat[j]+coef.es$R_S[j]),
                           V3.modify[j]+(coef.es$G3[j+1]-coef.es$G3[j]) - coef.es$mu_3[j]*V3.modify[j]), 0)
    V2_infected[j]=max(min(min(coef.es$y.temp[j]/(S.modify[j]+coef.es$varphikappa[j]*V1.modify[j]+coef.es$kappa[j]*V2.modify[j]+
                                                    coef.es$omegakappa[j]*V3.modify[j])*coef.es$kappa[j]*V2.modify[j], 
                               coef.es$SV.hat[j]-coef.es$SV.hat[j+1]+coef.es$mu_r[j]*coef.es$Ra.hat[j]+coef.es$R_S[j]-V3_infected[j]),
                           V2.modify[j]+(coef.es$G2[j+1]-coef.es$G2[j]) - (coef.es$G3[j+1]-coef.es$G3[j]) - coef.es$mu_2[j]*V2.modify[j]), 0)
    V1_infected[j]=max(min(min(coef.es$y.temp[j]/(S.modify[j]+coef.es$varphikappa[j]*V1.modify[j]+coef.es$kappa[j]*V2.modify[j]+
                                                    coef.es$omegakappa[j]*V3.modify[j])*coef.es$varphikappa[j]*V1.modify[j],
                               coef.es$SV.hat[j]-coef.es$SV.hat[j+1]+coef.es$mu_r[j]*coef.es$Ra.hat[j]+coef.es$R_S[j]-V3_infected[j]-V2_infected[j]),
                           V1.modify[j]+(coef.es$G1[j+1]-coef.es$G1[j])-(coef.es$G2[j+1]-coef.es$G2[j]) - coef.es$mu_1[j]*V1.modify[j]), 0)
    # S_infected[j]=max(min(coef.es$y.temp[j]/(S.modify[j]+coef.es$varphikappa[j]*V1.modify[j]+coef.es$kappa[j]*V2.modify[j]+
    #                                            coef.es$omegakappa[j]*V3.modify[j])*S.modify[j], 
    #                       coef.es$SV.hat[j]-coef.es$SV.hat[j+1]+coef.es$mu_r[j]*coef.es$Ra.hat[j]+coef.es$R_S[j]-
    #                       V3_infected[j]-V2_infected[j]-V1_infected[j]), 0)
    S_infected[j]=coef.es$y.temp[j]-V3_infected[j]-V2_infected[j]-V1_infected[j]
    # S.modify[j+1]=max(coef.es$SV.hat[j+1]-V1.modify[j+1]-V2.modify[j+1]-V3.modify[j+1], 0)
    S.modify[j+1]=S.modify[j]-S_infected[j]-(coef.es$G1[j+1]-coef.es$G1[j]) + 
      coef.es$mu_1[j]*V1.modify[j] + coef.es$mu_2[j]*V2.modify[j] + coef.es$mu_3[j]*V3.modify[j]+
      coef.es$mu_r[j]*coef.es$Ra.hat[j]+coef.es$R_S[j]
    V3.modify[j+1]=max(V3.modify[j]+(coef.es$G3[j+1]-coef.es$G3[j]) - coef.es$mu_3[j]*V3.modify[j]-V3_infected[j], 0)
    V2.modify[j+1]=max(V2.modify[j]+(coef.es$G2[j+1]-coef.es$G2[j]) - (coef.es$G3[j+1]-coef.es$G3[j]) - 
                         coef.es$mu_2[j]*V2.modify[j]-V2_infected[j], 0)
    # V1.modify[j+1]=max(min(V1.modify[j]+(coef.es$G1[j+1]-coef.es$G1[j])-(coef.es$G2[j+1]-coef.es$G2[j]) - 
    #                          coef.es$mu_1[j]*V1.modify[j]-V1_infected[j],
    #                        coef.es$SV.hat[j+1]-V2.modify[j+1]-V3.modify[j+1]), 0)
    V1.modify[j+1]=max(coef.es$SV.hat[j+1]-S.modify[j+1]-V2.modify[j+1]-V3.modify[j+1], 0)
    w_balan=balance_weight(country)
    if(V2.modify[j]+(coef.es$G2[j+1]-coef.es$G2[j]) - (coef.es$G3[j+1]-coef.es$G3[j]) - coef.es$mu_2[j]*V2.modify[j]<0 |
       V1.modify[j]+(coef.es$G1[j+1]-coef.es$G1[j])-(coef.es$G2[j+1]-coef.es$G2[j]) - coef.es$mu_1[j]*V1.modify[j]<0)
    {
      S.modify[j+1]=w_balan*max(S.modify[j+1], 0) + (1-w_balan)*max(coef.es$SV.hat[j+1]-V1.modify[j+1]-V2.modify[j+1]-V3.modify[j+1], 0)
      S_infected[j]=max(S.modify[j]-S.modify[j+1]-(coef.es$G1[j+1]-coef.es$G1[j]) + 
                          coef.es$mu_1[j]*V1.modify[j] + coef.es$mu_2[j]*V2.modify[j] + coef.es$mu_3[j]*V3.modify[j]+
                          coef.es$mu_r[j]*coef.es$Ra.hat[j]+coef.es$R_S[j], 0)
    }
  }
  coef.es$S_infected=S_infected
  coef.es$V1_infected=V1_infected
  coef.es$V2_infected=V2_infected
  coef.es$V3_infected=V3_infected
  coef.es$S.modify=S.modify
  coef.es$V1.modify=V1.modify
  coef.es$V2.modify=V2.modify
  coef.es$V3.modify=V3.modify
  
  coef.es$x1.temp = coef.es$allinfected*(coef.es$S.modify+coef.es$varphikappa*coef.es$V1.modify+
                                           coef.es$kappa*coef.es$V2.modify+coef.es$omegakappa*coef.es$V3.modify)/M
  coef.es$beta.mo=c(NA, sapply(1:(nrow_dat-2), kernel_reg, coef.es, nrow_dat-2, w0, 'y.temp', 'x1.temp'), NA)
  coef.es$beta.mo[coef.es$beta.mo <= 0]=0.001
  coef.es$beta.mo.11=c(NA, sapply(1:(nrow_dat-2), kernel_reg, coef.es, nrow_dat-2, 11, 'y.temp', 'x1.temp'), NA)
  coef.es$beta.mo.11[coef.es$beta.mo.11 <= 0]=0.001
  coef.es$Rt.mo=coef.es$beta.mo.11 *((1/coef.es$alpha[1]+D/r)*coef.es$theta+(1-coef.es$theta)* D/r)*
    (coef.es$S.modify+coef.es$varphikappa*coef.es$V1.modify+coef.es$kappa*coef.es$V2.modify+coef.es$omegakappa*coef.es$V3.modify)/M
  coef.es$Rt.nova=coef.es$beta.mo.11 *((1/coef.es$alpha[1]+D/r)*coef.es$theta+(1-coef.es$theta)* D/r)*coef.es$SV.hat/M
  coef.es$t1=1:nrow_dat
  coef.es$t1=coef.es$t1-t_nV
  coef.es$t2=1:nrow_dat
  coef.es$t2=coef.es$t2-t_booster
  
  write.csv(coef.es, paste0(path_table,country,'_estimates_preOmicron_reinfect.csv'))
  
  coefs_pre_Omi=rbind(coefs_pre_Omi, coef.es)
}

#------Omicron period--------
ranges=read.csv('~/Documents/COVID-19/data/range_Omicron.csv', stringsAsFactors = F)
for(j1 in c(27, 26, 25, 18, 13, 8, 4))
{
  set.seed(12)
  country=countries[j1] #country='US'
  print(country)
  estimates_final=read.csv(paste0(path_table,country,'_estimates_preOmicron_reinfect.csv'), stringsAsFactors = F)
  
  nknots5=country_knots$nknots5[country_knots$j1==j1]
  bspline6=bsplineS(x= (0:t_V)/t_V, breaks=seq(from=0, to=1, length.out = nknots5+2), norder=4)
  M=estimates_final$population[1]
  w0=provinces$w[provinces$provinces.en==country]
  t_omicron=which(as.character(as.Date(estimates_final$date))==as.character(as.Date(
    Omicron_period$Omicron_b[Omicron_period$province==country])))-1
  t_begin=t_omicron
  
  beta.after=estimates_final$beta.mo[estimates_final$t >= t_begin & estimates_final$t <= (t_begin+t_V)]
  beta.before=estimates_final$beta.mo[estimates_final$t >= (t_begin-t_V) & estimates_final$t <= t_begin]
  
  sp.after=as.data.frame(cbind(beta.after, bspline6))
  colnames(sp.after)=c('beta.after', paste0('base_', 0:(nknots5+3)))
  lm.spline=lm(beta.after ~ 0 + . , data = sp.after)
  mse.spline=mean((as.numeric(lm.spline$fitted.values)/beta.after-1)^2)
  summary(lm.spline)
  print(paste0('B-spline: ', sqrt(mse.spline)))
  coef.spline=as.numeric(coef(lm.spline))
  
  sp.before=as.data.frame(cbind(beta.before, bspline6))
  colnames(sp.before)=c('beta.before', paste0('base_', 0:(nknots5+3)))
  lm.spline.b=lm(beta.before ~ 0 + . , data = sp.before)
  mse.spline.b=mean((as.numeric(lm.spline.b$fitted.values)/beta.before-1)^2)
  summary(lm.spline.b)
  print(paste0('B-spline: ', sqrt(mse.spline.b)))
  coef.spline.b=as.numeric(coef(lm.spline.b))
  
  beta.mod=data.frame(t=c((-t_V):t_V),fit=c(beta.before, beta.after[-1]))
  beta.fit=data.frame(t=c((-t_V):0, 0:t_V), fit=c(lm.spline.b$fitted.values, lm.spline$fitted.values), type=rep(c('before','after'),each=t_V+1))
  ggplot(data=beta.fit)+geom_line(aes(x=t, y=fit, group=type, color=type))+
    geom_line(data=beta.mod, aes(x=t, y=fit),color='black')
  
  dat_sim=estimates_final[t_begin:(t_begin+t_V),]
  
  beta.range = round(apply(rbind(coef.spline, coef.spline.b), 2, range), 2)
  #varphi, kappa, omega, B-spline coefficients
  lower=c(1.1, ranges$kappa_l[ranges$j1==j1], 0.2, max(0.01, beta.range[1,1]-0.01), beta.range[1,-1]-0.01)
  upper=c(3.5, ranges$kappa_u[ranges$j1==j1], 0.8, min(1, beta.range[2,1]), beta.range[2,-1])
  
  rho.ops=data.frame(rho.ops=c(seq(from=lower[1], to=upper[1], by=0.1)))
  kappa.ops=data.frame(kappa.ops=seq(from=lower[2], to=upper[2], by=0.01))
  omega.ops=data.frame(omega.ops=seq(from=lower[3], to=upper[3], by=0.1))
  rho.k=merge(merge(rho.ops, kappa.ops), omega.ops)
  rho.k=rho.k[rho.k[,1]*rho.k[,2] < ranges$rhokappa_u[ranges$j1==j1],]
  rho.k=rho.k[rho.k[,1]*rho.k[,2] > ranges$rhokappa_l[ranges$j1==j1],]
  rho.k=rho.k[rho.k[,3]*rho.k[,2] < ranges$omegakappa_u[ranges$j1==j1],]
  rho.k=rho.k[rho.k[,3]*rho.k[,2] > ranges$omegakappa_l[ranges$j1==j1],]
  
  print(nrow(rho.k)*prod((upper-lower)[-(1:3)]/0.01+1))
  
  beta.o0.ops=data.frame(beta.o0.ops=seq(from=lower[4], to=upper[4],by=0.01))
  beta.o1.ops=data.frame(beta.o1.ops=seq(from=lower[5], to=upper[5],by=0.01))
  beta.o2.ops=data.frame(beta.o2.ops=seq(from=lower[6], to=upper[6],by=0.01))
  kb.ops=merge(merge(merge(rho.k, beta.o0.ops), beta.o1.ops), beta.o2.ops)
  colnames(kb.ops)=c('rho','kappa','omega',paste0('beta.o',0:2))
  
  for(cc in 0:nknots5)
  {
    beta.ops=data.frame(beta.ops=seq(from=lower[7+cc], to=upper[7+cc],by=0.01))
    kb.ops=merge(kb.ops, beta.ops)
    colnames(kb.ops)[7+cc]=paste0('beta.o', cc+3)
  }
  
  #-------grid search----------
  dats_sim=as.matrix(dat_sim[, c('t', 'N', 'R.d', 'R.r', 'R', 'I', 'G1', 'G2', 'alpha', 'Ra.hat',
                                 'E.hat', 'Ea.hat', 'gamma.d', 'gamma.r', 'V1.modify', 'V2.modify', 'G3', 'V3.modify')])
  colnames(dats_sim)=NULL
  
  inter_b_min = 11; inter_b_max = 20; inter_l = 30
  ptm = proc.time()
  cl <- makeCluster(64)
  object.fun = t(parApply(cl, kb.ops, 1, object_kappa_booster_reinfect, dats_sim, 
                          t_V, iter, dat_sim$theta, M, dat_sim$mu_1[1], dat_sim$mu_2[1], dat_sim$mu_3[1], dat_sim$mu_r[1], 
                          r.beta, bspline6, inter_b_min, inter_b_max, inter_l))
  stopCluster(cl)
  print(paste0('Grid search: ',round(as.numeric((proc.time() - ptm)[3]), 3)))
  
  object.f.ess=cbind(kb.ops[apply(object.fun, 2, which.min),], apply(object.fun, 2, min))
  object.f.ess=as.data.frame(object.f.ess)
  colnames(object.f.ess)=c('rho','kappa','omega',paste0('beta.o',0:(nknots5+3)),'object') 
  rownames(object.f.ess)=NULL
  object.f.ess$inter_b=inter_b_min:inter_b_max
  object.f.ess$inter_l=inter_l
  object.f.ess$alpha=estimates_final$alpha[1]
  
  write.csv(object.f.ess,paste0(path_table,country,'_criteria_grid_final_inter_Omicron_reinfect.csv'))
  print(object.f.ess)
  
  object.fun=as.data.frame(object.fun)
  colnames(object.fun)=as.character(inter_b_min:inter_b_max)
  object.fun = cbind(kb.ops, object.fun)
  save(object.fun, file = paste0(path_table,country,'_kappa_objfunction_Omicron_reinfect.RData'))
  
}

country_kappas=c()
country_data=c()
kappa_estis=c()
for(j1 in c(27, 26, 25, 18, 13, 8, 4))
{
  country=countries[j1] #country='US'
  print(country)
  estimates_final=read.csv(paste0(path_table,country,'_estimates_preOmicron_reinfect.csv'), stringsAsFactors = F)
  
  M=estimates_final$population[1]
  w0=provinces$w[provinces$provinces.en==country]
  t_nV=min(which(estimates_final$G1>0))-1
  t_omicron=which(as.character(as.Date(estimates_final$date))==as.character(as.Date(
    Omicron_period$Omicron_b[Omicron_period$province==country])))-1
  
  object.f.ess = read.csv(paste0(path_table,country,'_criteria_grid_final_inter_Omicron_reinfect.csv'), stringsAsFactors = F)
  print(object.f.ess)
  
  choose_inter_b=18
  if(j1 %in% c(26)) choose_inter_b=17
  if(j1 %in% c(4, 13, 25)) choose_inter_b=20
  
  country_kappas=rbind(country_kappas, data.frame(country=country, rho=object.f.ess$rho, kappa=object.f.ess$kappa,
                                                  omega=object.f.ess$omega, inter_b=object.f.ess$inter_b))
  coef.es=estimates_final
  nrow_dat=nrow(coef.es)
  rho.es=object.f.ess$rho[object.f.ess$inter_b==choose_inter_b]
  kappa.es=object.f.ess$kappa[object.f.ess$inter_b==choose_inter_b]
  omega.es=object.f.ess$omega[object.f.ess$inter_b==choose_inter_b]
  
  kappa_estis=rbind(kappa_estis, c(country, rho.es*kappa.es, kappa.es, omega.es*kappa.es))
  
  loop_b=t_omicron
  coef.es$kappa[(loop_b+1):nrow_dat]=kappa.es
  coef.es$varphikappa[(loop_b+1):nrow_dat]=rho.es*kappa.es
  coef.es$omegakappa[(loop_b+1):nrow_dat]=omega.es*kappa.es
  
  S.modify=coef.es$S.modify
  V1.modify=coef.es$V1.modify
  V2.modify=coef.es$V2.modify
  V3.modify=coef.es$V3.modify
  
  S_infected=coef.es$S_infected
  V1_infected=coef.es$V1_infected
  V2_infected=coef.es$V2_infected
  V3_infected=coef.es$V3_infected
  for(j in loop_b:(nrow_dat-1))
  {
    #transfer in [j, j+1]
    V3_infected[j]=max(min(min(coef.es$y.temp[j]/(S.modify[j]+coef.es$varphikappa[j]*V1.modify[j]+coef.es$kappa[j]*V2.modify[j]+
                                                    coef.es$omegakappa[j]*V3.modify[j])*coef.es$omegakappa[j]*V3.modify[j],
                               coef.es$SV.hat[j]-coef.es$SV.hat[j+1]+coef.es$mu_r[j]*coef.es$Ra.hat[j]+coef.es$R_S[j]),
                           V3.modify[j]+(coef.es$G3[j+1]-coef.es$G3[j]) - coef.es$mu_3[j]*V3.modify[j]), 0)
    V2_infected[j]=max(min(min(coef.es$y.temp[j]/(S.modify[j]+coef.es$varphikappa[j]*V1.modify[j]+coef.es$kappa[j]*V2.modify[j]+
                                                    coef.es$omegakappa[j]*V3.modify[j])*coef.es$kappa[j]*V2.modify[j], 
                               coef.es$SV.hat[j]-coef.es$SV.hat[j+1]+coef.es$mu_r[j]*coef.es$Ra.hat[j]+coef.es$R_S[j]-V3_infected[j]),
                           V2.modify[j]+(coef.es$G2[j+1]-coef.es$G2[j]) - (coef.es$G3[j+1]-coef.es$G3[j]) - coef.es$mu_2[j]*V2.modify[j]), 0)
    V1_infected[j]=max(min(min(coef.es$y.temp[j]/(S.modify[j]+coef.es$varphikappa[j]*V1.modify[j]+coef.es$kappa[j]*V2.modify[j]+
                                                    coef.es$omegakappa[j]*V3.modify[j])*coef.es$varphikappa[j]*V1.modify[j],
                               coef.es$SV.hat[j]-coef.es$SV.hat[j+1]+coef.es$mu_r[j]*coef.es$Ra.hat[j]+coef.es$R_S[j]-V3_infected[j]-V2_infected[j]),
                           V1.modify[j]+(coef.es$G1[j+1]-coef.es$G1[j])-(coef.es$G2[j+1]-coef.es$G2[j]) - coef.es$mu_1[j]*V1.modify[j]), 0)
    # S_infected[j]=max(min(coef.es$y.temp[j]/(S.modify[j]+coef.es$varphikappa[j]*V1.modify[j]+coef.es$kappa[j]*V2.modify[j]+
    #                                            coef.es$omegakappa[j]*V3.modify[j])*S.modify[j], 
    #                       coef.es$SV.hat[j]-coef.es$SV.hat[j+1]+coef.es$mu_r[j]*coef.es$Ra.hat[j]+coef.es$R_S[j]-
    #                       V3_infected[j]-V2_infected[j]-V1_infected[j]), 0)
    S_infected[j]=coef.es$y.temp[j]-V3_infected[j]-V2_infected[j]-V1_infected[j]
    # S.modify[j+1]=max(coef.es$SV.hat[j+1]-V1.modify[j+1]-V2.modify[j+1]-V3.modify[j+1], 0)
    S.modify[j+1]=S.modify[j]-S_infected[j]-(coef.es$G1[j+1]-coef.es$G1[j]) + 
      coef.es$mu_1[j]*V1.modify[j] + coef.es$mu_2[j]*V2.modify[j] + coef.es$mu_3[j]*V3.modify[j]+
      coef.es$mu_r[j]*coef.es$Ra.hat[j]+coef.es$R_S[j]
    V3.modify[j+1]=max(V3.modify[j]+(coef.es$G3[j+1]-coef.es$G3[j]) - coef.es$mu_3[j]*V3.modify[j]-V3_infected[j], 0)
    V2.modify[j+1]=max(V2.modify[j]+(coef.es$G2[j+1]-coef.es$G2[j]) - (coef.es$G3[j+1]-coef.es$G3[j]) - 
                         coef.es$mu_2[j]*V2.modify[j]-V2_infected[j], 0)
    # V1.modify[j+1]=max(min(V1.modify[j]+(coef.es$G1[j+1]-coef.es$G1[j])-(coef.es$G2[j+1]-coef.es$G2[j]) - 
    #                          coef.es$mu_1[j]*V1.modify[j]-V1_infected[j],
    #                        coef.es$SV.hat[j+1]-V2.modify[j+1]-V3.modify[j+1]), 0)
    V1.modify[j+1]=max(coef.es$SV.hat[j+1]-S.modify[j+1]-V2.modify[j+1]-V3.modify[j+1], 0)
    w_balan=balance_weight(country)
    if(V2.modify[j]+(coef.es$G2[j+1]-coef.es$G2[j]) - (coef.es$G3[j+1]-coef.es$G3[j]) - coef.es$mu_2[j]*V2.modify[j]<0 |
       V1.modify[j]+(coef.es$G1[j+1]-coef.es$G1[j])-(coef.es$G2[j+1]-coef.es$G2[j]) - coef.es$mu_1[j]*V1.modify[j]<0)
    {
      S.modify[j+1]=w_balan*max(S.modify[j+1], 0) + (1-w_balan)*max(coef.es$SV.hat[j+1]-V1.modify[j+1]-V2.modify[j+1]-V3.modify[j+1], 0)
      S_infected[j]=max(S.modify[j]-S.modify[j+1]-(coef.es$G1[j+1]-coef.es$G1[j]) + 
                          coef.es$mu_1[j]*V1.modify[j] + coef.es$mu_2[j]*V2.modify[j] + coef.es$mu_3[j]*V3.modify[j]+
                          coef.es$mu_r[j]*coef.es$Ra.hat[j]+coef.es$R_S[j], 0)
    }
  }
  coef.es$S_infected=S_infected
  coef.es$V1_infected=V1_infected
  coef.es$V2_infected=V2_infected
  coef.es$V3_infected=V3_infected
  coef.es$S.modify=S.modify
  coef.es$V1.modify=V1.modify
  coef.es$V2.modify=V2.modify
  coef.es$V3.modify=V3.modify
  
  coef.es$x1.temp = coef.es$allinfected*(coef.es$S.modify+coef.es$varphikappa*coef.es$V1.modify+
                                           coef.es$kappa*coef.es$V2.modify+coef.es$omegakappa*coef.es$V3.modify)/M
  coef.es$beta.mo=c(NA, sapply(1:(nrow_dat-2), kernel_reg, coef.es, nrow_dat-2, w0, 'y.temp', 'x1.temp'), NA)
  coef.es$beta.mo[coef.es$beta.mo <= 0]=0.001
  coef.es$beta.mo.11=c(NA, sapply(1:(nrow_dat-2), kernel_reg, coef.es, nrow_dat-2, 11, 'y.temp', 'x1.temp'), NA)
  coef.es$beta.mo.11[coef.es$beta.mo.11 <= 0]=0.001
  coef.es$Rt.mo=coef.es$beta.mo.11 *((1/coef.es$alpha[1]+D/r)*coef.es$theta+(1-coef.es$theta)* D/r)*
    (coef.es$S.modify+coef.es$varphikappa*coef.es$V1.modify+coef.es$kappa*coef.es$V2.modify+coef.es$omegakappa*coef.es$V3.modify)/M
  coef.es$Rt.nova=coef.es$beta.mo.11 *((1/coef.es$alpha[1]+D/r)*coef.es$theta+(1-coef.es$theta)* D/r)*coef.es$SV.hat/M
  
  write.csv(coef.es, paste0(path_table,country,'_estimates_Omicron_reinfect.csv'))
  
  country_data=rbind(country_data, coef.es[, c("t", "province", "population", "date", "infected", "dead", "recovered", "R.d", "R.r",
                                               "N", "R", "I", "G1", "G2", "G3", 'R_S', 'mu_r', "alpha", "theta", "S.hat", "Ra.hat", "E.hat", "Ea.hat",
                                               "x.temp", "y.temp", "allinfected", "gamma", "gamma.d", "gamma.r", "beta", "Rt", 
                                               "SV.hat", "V1.hat", "V2.hat", "V3.hat", "dN", "dG1", "dG2", "dG3", 
                                               "scale_dN", "scale_dG1", "scale_dG2", "scale_dG3", 
                                               "scale_G1", "scale_G2", "scale_G3", "scale_N", "scale_I", 
                                               "S_infected", "V1_infected", "V2_infected", "V3_infected",
                                               "S.modify", "V1.modify", "V2.modify", "V3.modify", "x1.temp", 
                                               "kappa", "varphikappa", "omegakappa", 'mu_1', 'mu_2', 'mu_3', 
                                               "beta.mo", "beta.mo.11", "Rt.mo", "Rt.nova", 
                                               "gamma.r.11", "gamma.d.11", "t1", "t2")])
}
country_kappas$kappa1=country_kappas$rho*country_kappas$kappa
country_kappas$kappa3=country_kappas$omega*country_kappas$kappa
plot_k12=melt(country_kappas, id=c('country', 'inter_b'), measure=c('kappa1', 'kappa', 'kappa3'))
plot_k12$variable=factor(plot_k12$variable, levels=c('kappa1', 'kappa', 'kappa3'), labels=c('kappa1', 'kappa2', 'kappa3'))
ggplot(plot_k12, aes(x=inter_b, y=value, color=factor(variable)))+
  geom_point()+geom_line()+geom_text(data=plot_k12, aes(x=inter_b, y=value, label=value), color='black',
                                     size=3.5, hjust = 0, nudge_x = -0.8, check_overlap = TRUE)+
  facet_wrap(.~ country, nrow=2, ncol=4, scales = "fixed")+
  labs(x='The left end point of the interval', y='kappa', title='kappa', color='')+
  theme(legend.position = 'bottom')

VEs2=as.data.frame(kappa_estis)
VEs2[,-1]=1-apply(kappa_estis[,-1], 2, as.numeric)
colnames(VEs2)=c('province', 'partial', 'full', 'booster')
VEs2$period='Omicron-reported'

#------Omicron dominated period--------
ranges=read.csv('~/Documents/COVID-19/data/range_Omicron_dominate.csv', stringsAsFactors = F)
for(j1 in c(27, 26, 25, 18, 13, 8, 4))
{
  set.seed(12)
  country=countries[j1] #country='US'
  print(country)
  estimates_final=read.csv(paste0(path_table,country,'_estimates_Omicron_reinfect.csv'), stringsAsFactors = F)
  
  nknots6=country_knots$nknots6[country_knots$j1==j1]
  bspline7=bsplineS(x= (0:t_V)/t_V, breaks=seq(from=0, to=1, length.out = nknots6+2), norder=4)
  M=estimates_final$population[1]
  w0=provinces$w[provinces$provinces.en==country]
  t_omicron_d=which(as.character(as.Date(estimates_final$date))==as.character(as.Date(
    Omicron_period$Omicron_d[Omicron_period$province==country])))-1
  t_begin=t_omicron_d
  
  beta.after=estimates_final$beta.mo[estimates_final$t >= t_begin & estimates_final$t <= (t_begin+t_V)]
  beta.before=estimates_final$beta.mo[estimates_final$t >= (t_begin-t_V) & estimates_final$t <= t_begin]
  
  sp.after=as.data.frame(cbind(beta.after, bspline7))
  colnames(sp.after)=c('beta.after', paste0('base_', 0:(nknots6+3)))
  lm.spline=lm(beta.after ~ 0 + . , data = sp.after)
  mse.spline=mean((as.numeric(lm.spline$fitted.values)/beta.after-1)^2)
  summary(lm.spline)
  print(paste0('B-spline: ', sqrt(mse.spline)))
  coef.spline=as.numeric(coef(lm.spline))
  
  sp.before=as.data.frame(cbind(beta.before, bspline7))
  colnames(sp.before)=c('beta.before', paste0('base_', 0:(nknots6+3)))
  lm.spline.b=lm(beta.before ~ 0 + . , data = sp.before)
  mse.spline.b=mean((as.numeric(lm.spline.b$fitted.values)/beta.before-1)^2)
  summary(lm.spline.b)
  print(paste0('B-spline: ', sqrt(mse.spline.b)))
  coef.spline.b=as.numeric(coef(lm.spline.b))
  
  beta.mod=data.frame(t=c((-t_V):t_V),fit=c(beta.before, beta.after[-1]))
  beta.fit=data.frame(t=c((-t_V):0, 0:t_V), fit=c(lm.spline.b$fitted.values, lm.spline$fitted.values), type=rep(c('before','after'),each=t_V+1))
  ggplot(data=beta.fit)+geom_line(aes(x=t, y=fit, group=type, color=type))+
    geom_line(data=beta.mod, aes(x=t, y=fit),color='black')
  
  dat_sim=estimates_final[t_begin:(t_begin+t_V),]
  
  beta.range = round(apply(rbind(coef.spline, coef.spline.b), 2, range), 2)
  #varphi, kappa, omega, B-spline coefficients
  lower=c(1.1, ranges$kappa_l[ranges$j1==j1], 0.2, max(0.01, beta.range[1,1]-0.01), beta.range[1,-1]-0.01)
  upper=c(3.5, ranges$kappa_u[ranges$j1==j1], 0.8, min(1, beta.range[2,1]), beta.range[2,-1])
  
  rho.ops=data.frame(rho.ops=c(seq(from=lower[1], to=upper[1], by=0.2)))
  kappa.ops=data.frame(kappa.ops=seq(from=lower[2], to=upper[2], by=0.01))
  omega.ops=data.frame(omega.ops=seq(from=lower[3], to=upper[3], by=0.2))
  rho.k=merge(merge(rho.ops, kappa.ops), omega.ops)
  rho.k=rho.k[rho.k[,1]*rho.k[,2] < ranges$rhokappa_u[ranges$j1==j1],]
  rho.k=rho.k[rho.k[,1]*rho.k[,2] > ranges$rhokappa_l[ranges$j1==j1],]
  rho.k=rho.k[rho.k[,3]*rho.k[,2] < ranges$omegakappa_u[ranges$j1==j1],]
  rho.k=rho.k[rho.k[,3]*rho.k[,2] > ranges$omegakappa_l[ranges$j1==j1],]
  
  print(nrow(rho.k)*prod((upper-lower)[-(1:3)]/0.01+1))
  
  beta.o0.ops=data.frame(beta.o0.ops=seq(from=lower[4], to=upper[4],by=0.01))
  beta.o1.ops=data.frame(beta.o1.ops=seq(from=lower[5], to=upper[5],by=0.01))
  beta.o2.ops=data.frame(beta.o2.ops=seq(from=lower[6], to=upper[6],by=0.01))
  kb.ops=merge(merge(merge(rho.k, beta.o0.ops), beta.o1.ops), beta.o2.ops)
  colnames(kb.ops)=c('rho','kappa','omega',paste0('beta.o',0:2))
  
  for(cc in 0:nknots6)
  {
    beta.ops=data.frame(beta.ops=seq(from=lower[7+cc], to=upper[7+cc],by=0.01))
    kb.ops=merge(kb.ops, beta.ops)
    colnames(kb.ops)[7+cc]=paste0('beta.o', cc+3)
  }
  
  #-------grid search----------
  dats_sim=as.matrix(dat_sim[, c('t', 'N', 'R.d', 'R.r', 'R', 'I', 'G1', 'G2', 'alpha', 'Ra.hat',
                                 'E.hat', 'Ea.hat', 'gamma.d', 'gamma.r', 'V1.modify', 'V2.modify', 'G3', 'V3.modify')])
  colnames(dats_sim)=NULL
  
  inter_b_min = 11; inter_b_max = 20; inter_l = 30
  ptm = proc.time()
  cl <- makeCluster(64)
  object.fun = t(parApply(cl, kb.ops, 1, object_kappa_booster_reinfect, dats_sim, 
                          t_V, iter, dat_sim$theta, M, dat_sim$mu_1[1], dat_sim$mu_2[1], dat_sim$mu_3[1], dat_sim$mu_r[1], 
                          r.beta, bspline7, inter_b_min, inter_b_max, inter_l))
  stopCluster(cl)
  print(paste0('Grid search: ',round(as.numeric((proc.time() - ptm)[3]), 3)))
  
  object.f.ess=cbind(kb.ops[apply(object.fun, 2, which.min),], apply(object.fun, 2, min))
  object.f.ess=as.data.frame(object.f.ess)
  colnames(object.f.ess)=c('rho','kappa','omega',paste0('beta.o',0:(nknots6+3)),'object') 
  rownames(object.f.ess)=NULL
  object.f.ess$inter_b=inter_b_min:inter_b_max
  object.f.ess$inter_l=inter_l
  object.f.ess$alpha=estimates_final$alpha[1]
  
  write.csv(object.f.ess,paste0(path_table,country,'_criteria_grid_final_inter_Omicron_dominate_reinfect.csv'))
  print(object.f.ess)
  
  object.fun=as.data.frame(object.fun)
  colnames(object.fun)=as.character(inter_b_min:inter_b_max)
  object.fun = cbind(kb.ops, object.fun)
  save(object.fun, file = paste0(path_table,country,'_kappa_objfunction_Omicron_dominate_reinfect.RData'))
}

country_kappas=c()
country_data=c()
kappa_estis=c()
for(j1 in c(27, 26, 25, 18, 13, 8, 4))
{
  country=countries[j1] #country='US'
  print(country)
  estimates_final=read.csv(paste0(path_table,country,'_estimates_Omicron_reinfect.csv'), stringsAsFactors = F)
  
  M=estimates_final$population[1]
  w0=provinces$w[provinces$provinces.en==country]
  t_nV=min(which(estimates_final$G1>0))-1
  t_omicron_d=which(as.character(as.Date(estimates_final$date))==as.character(as.Date(
    Omicron_period$Omicron_d[Omicron_period$province==country])))-1
  
  object.f.ess = read.csv(paste0(path_table,country,'_criteria_grid_final_inter_Omicron_dominate_reinfect.csv'), stringsAsFactors = F)
  print(object.f.ess)
  
  choose_inter_b=18
  if(j1==8) choose_inter_b=20
  
  country_kappas=rbind(country_kappas, data.frame(country=country, rho=object.f.ess$rho, kappa=object.f.ess$kappa,
                                                  omega=object.f.ess$omega, inter_b=object.f.ess$inter_b))
  coef.es=estimates_final
  nrow_dat=nrow(coef.es)
  rho.es=object.f.ess$rho[object.f.ess$inter_b==choose_inter_b]
  kappa.es=object.f.ess$kappa[object.f.ess$inter_b==choose_inter_b]
  omega.es=object.f.ess$omega[object.f.ess$inter_b==choose_inter_b]
  
  kappa_estis=rbind(kappa_estis, c(country, rho.es*kappa.es, kappa.es, omega.es*kappa.es))
  
  loop_b=t_omicron_d
  coef.es$kappa[(loop_b+1):nrow_dat]=kappa.es
  coef.es$varphikappa[(loop_b+1):nrow_dat]=rho.es*kappa.es
  coef.es$omegakappa[(loop_b+1):nrow_dat]=omega.es*kappa.es
  
  S.modify=coef.es$S.modify
  V1.modify=coef.es$V1.modify
  V2.modify=coef.es$V2.modify
  V3.modify=coef.es$V3.modify
  
  S_infected=coef.es$S_infected
  V1_infected=coef.es$V1_infected
  V2_infected=coef.es$V2_infected
  V3_infected=coef.es$V3_infected
  for(j in loop_b:(nrow_dat-1))
  {
    #transfer in [j, j+1]
    V3_infected[j]=max(min(min(coef.es$y.temp[j]/(S.modify[j]+coef.es$varphikappa[j]*V1.modify[j]+coef.es$kappa[j]*V2.modify[j]+
                                                    coef.es$omegakappa[j]*V3.modify[j])*coef.es$omegakappa[j]*V3.modify[j],
                               coef.es$SV.hat[j]-coef.es$SV.hat[j+1]+coef.es$mu_r[j]*coef.es$Ra.hat[j]+coef.es$R_S[j]),
                           V3.modify[j]+(coef.es$G3[j+1]-coef.es$G3[j]) - coef.es$mu_3[j]*V3.modify[j]), 0)
    V2_infected[j]=max(min(min(coef.es$y.temp[j]/(S.modify[j]+coef.es$varphikappa[j]*V1.modify[j]+coef.es$kappa[j]*V2.modify[j]+
                                                    coef.es$omegakappa[j]*V3.modify[j])*coef.es$kappa[j]*V2.modify[j], 
                               coef.es$SV.hat[j]-coef.es$SV.hat[j+1]+coef.es$mu_r[j]*coef.es$Ra.hat[j]+coef.es$R_S[j]-V3_infected[j]),
                           V2.modify[j]+(coef.es$G2[j+1]-coef.es$G2[j]) - (coef.es$G3[j+1]-coef.es$G3[j]) - coef.es$mu_2[j]*V2.modify[j]), 0)
    V1_infected[j]=max(min(min(coef.es$y.temp[j]/(S.modify[j]+coef.es$varphikappa[j]*V1.modify[j]+coef.es$kappa[j]*V2.modify[j]+
                                                    coef.es$omegakappa[j]*V3.modify[j])*coef.es$varphikappa[j]*V1.modify[j],
                               coef.es$SV.hat[j]-coef.es$SV.hat[j+1]+coef.es$mu_r[j]*coef.es$Ra.hat[j]+coef.es$R_S[j]-V3_infected[j]-V2_infected[j]),
                           V1.modify[j]+(coef.es$G1[j+1]-coef.es$G1[j])-(coef.es$G2[j+1]-coef.es$G2[j]) - coef.es$mu_1[j]*V1.modify[j]), 0)
    # S_infected[j]=max(min(coef.es$y.temp[j]/(S.modify[j]+coef.es$varphikappa[j]*V1.modify[j]+coef.es$kappa[j]*V2.modify[j]+
    #                                            coef.es$omegakappa[j]*V3.modify[j])*S.modify[j], 
    #                       coef.es$SV.hat[j]-coef.es$SV.hat[j+1]+coef.es$mu_r[j]*coef.es$Ra.hat[j]+coef.es$R_S[j]-
    #                       V3_infected[j]-V2_infected[j]-V1_infected[j]), 0)
    S_infected[j]=coef.es$y.temp[j]-V3_infected[j]-V2_infected[j]-V1_infected[j]
    # S.modify[j+1]=max(coef.es$SV.hat[j+1]-V1.modify[j+1]-V2.modify[j+1]-V3.modify[j+1], 0)
    S.modify[j+1]=S.modify[j]-S_infected[j]-(coef.es$G1[j+1]-coef.es$G1[j]) + 
      coef.es$mu_1[j]*V1.modify[j] + coef.es$mu_2[j]*V2.modify[j] + coef.es$mu_3[j]*V3.modify[j]+
      coef.es$mu_r[j]*coef.es$Ra.hat[j]+coef.es$R_S[j]
    V3.modify[j+1]=max(V3.modify[j]+(coef.es$G3[j+1]-coef.es$G3[j]) - coef.es$mu_3[j]*V3.modify[j]-V3_infected[j], 0)
    V2.modify[j+1]=max(V2.modify[j]+(coef.es$G2[j+1]-coef.es$G2[j]) - (coef.es$G3[j+1]-coef.es$G3[j]) - 
                         coef.es$mu_2[j]*V2.modify[j]-V2_infected[j], 0)
    # V1.modify[j+1]=max(min(V1.modify[j]+(coef.es$G1[j+1]-coef.es$G1[j])-(coef.es$G2[j+1]-coef.es$G2[j]) - 
    #                          coef.es$mu_1[j]*V1.modify[j]-V1_infected[j],
    #                        coef.es$SV.hat[j+1]-V2.modify[j+1]-V3.modify[j+1]), 0)
    V1.modify[j+1]=max(coef.es$SV.hat[j+1]-S.modify[j+1]-V2.modify[j+1]-V3.modify[j+1], 0)
    w_balan=balance_weight(country)
    if(V2.modify[j]+(coef.es$G2[j+1]-coef.es$G2[j]) - (coef.es$G3[j+1]-coef.es$G3[j]) - coef.es$mu_2[j]*V2.modify[j]<0 |
       V1.modify[j]+(coef.es$G1[j+1]-coef.es$G1[j])-(coef.es$G2[j+1]-coef.es$G2[j]) - coef.es$mu_1[j]*V1.modify[j]<0)
    {
      S.modify[j+1]=w_balan*max(S.modify[j+1], 0) + (1-w_balan)*max(coef.es$SV.hat[j+1]-V1.modify[j+1]-V2.modify[j+1]-V3.modify[j+1], 0)
      S_infected[j]=max(S.modify[j]-S.modify[j+1]-(coef.es$G1[j+1]-coef.es$G1[j]) + 
                          coef.es$mu_1[j]*V1.modify[j] + coef.es$mu_2[j]*V2.modify[j] + coef.es$mu_3[j]*V3.modify[j]+
                          coef.es$mu_r[j]*coef.es$Ra.hat[j]+coef.es$R_S[j], 0)
    }
  }
  coef.es$S_infected=S_infected
  coef.es$V1_infected=V1_infected
  coef.es$V2_infected=V2_infected
  coef.es$V3_infected=V3_infected
  coef.es$S.modify=S.modify
  coef.es$V1.modify=V1.modify
  coef.es$V2.modify=V2.modify
  coef.es$V3.modify=V3.modify
  
  coef.es$x1.temp = coef.es$allinfected*(coef.es$S.modify+coef.es$varphikappa*coef.es$V1.modify+
                                           coef.es$kappa*coef.es$V2.modify+coef.es$omegakappa*coef.es$V3.modify)/M
  coef.es$beta.mo=c(NA, sapply(1:(nrow_dat-2), kernel_reg, coef.es, nrow_dat-2, w0, 'y.temp', 'x1.temp'), NA)
  coef.es$beta.mo[coef.es$beta.mo <= 0]=0.001
  coef.es$beta.mo.11=c(NA, sapply(1:(nrow_dat-2), kernel_reg, coef.es, nrow_dat-2, 11, 'y.temp', 'x1.temp'), NA)
  coef.es$beta.mo.11[coef.es$beta.mo.11 <= 0]=0.001
  # coef.es$beta.spline[loop_b:(loop_b+t_V)]=sapply(sapply(bspline5 %*% as.vector(as.numeric(
  #   object.f.ess[object.f.ess$inter_b==choose_inter_b, paste0('beta.o',0:(nknots4+3))])), max, 0.001), min, 1)
  coef.es$S.modify=coef.es$SV.hat-coef.es$V1.modify-coef.es$V2.modify-coef.es$V3.modify
  coef.es$Rt.mo=coef.es$beta.mo.11 *((1/coef.es$alpha[1]+D/r)*coef.es$theta+(1-coef.es$theta)* D/r)*
    (coef.es$S.modify+coef.es$varphikappa*coef.es$V1.modify+coef.es$kappa*coef.es$V2.modify+coef.es$omegakappa*coef.es$V3.modify)/M
  coef.es$Rt.one=coef.es$beta.mo.11 *((1/coef.es$alpha[1]+D/r)*coef.es$theta+(1-coef.es$theta)* D/r)*
    (coef.es$S.modify+coef.es$varphikappa*(coef.es$V1.modify+coef.es$V2.modify+coef.es$V3.modify))/M
  coef.es$Rt.noboost=coef.es$beta.mo.11 *((1/coef.es$alpha[1]+D/r)*coef.es$theta+(1-coef.es$theta)* D/r)*
    (coef.es$S.modify+coef.es$varphikappa*coef.es$V1.modify+coef.es$kappa*(coef.es$V2.modify+coef.es$V3.modify))/M
  coef.es$Rt.nova=coef.es$beta.mo.11 *((1/coef.es$alpha[1]+D/r)*coef.es$theta+(1-coef.es$theta)* D/r)*coef.es$SV.hat/M
  
  coef.es$rho=coef.es$varphikappa/coef.es$kappa
  coef.es$rho[is.nan(coef.es$rho)]=0
  coef.es$omega=coef.es$omegakappa/coef.es$kappa
  coef.es$omega[is.nan(coef.es$omega)]=0
  coef.es$dRd=c(diff(coef.es$R.d), NA)
  coef.es$dRr=c(diff(coef.es$R.r), NA)
  
  write.csv(coef.es, paste0(path_table,country,'_estimates_Omicron_dominate_reinfect.csv'))
  
  country_data=rbind(country_data, coef.es[, c("t", "province", "population", "date", "infected", "dead", "recovered", "R.d", "R.r",
                                               "N", "R", "I", "G1", "G2", "G3", 'R_S', 'mu_r', "alpha", "theta", "S.hat", "Ra.hat", "E.hat", "Ea.hat",
                                               "x.temp", "y.temp", "allinfected", "gamma", "gamma.d", "gamma.r", "beta", "Rt", 
                                               "SV.hat", "V1.hat", "V2.hat", "V3.hat", "dN", "dG1", "dG2", "dG3", 
                                               "scale_dN", "scale_dG1", "scale_dG2", "scale_dG3", 
                                               "scale_G1", "scale_G2", "scale_G3", "scale_N", "scale_I", 
                                               "S_infected", "V1_infected", "V2_infected", "V3_infected",
                                               "S.modify", "V1.modify", "V2.modify", "V3.modify", "x1.temp", 
                                               "kappa", "varphikappa", "omegakappa", 'mu_1', 'mu_2', 'mu_3', "beta.mo", "beta.mo.11", 
                                               "Rt.mo", "Rt.noboost", "Rt.one", "Rt.nova", 
                                               "gamma.r.11", "gamma.d.11", "t1", "t2", "rho", "omega", "dRd", "dRr")])
}
write.csv(country_data, '~/Documents/COVID-19/data/country_data_theta1_reinfect.csv')
country_kappas$kappa1=country_kappas$rho*country_kappas$kappa
country_kappas$kappa3=country_kappas$omega*country_kappas$kappa
plot_k12=melt(country_kappas, id=c('country', 'inter_b'), measure=c('kappa1', 'kappa', 'kappa3'))
plot_k12$variable=factor(plot_k12$variable, levels=c('kappa1', 'kappa', 'kappa3'), labels=c('kappa1', 'kappa2', 'kappa3'))
ggplot(plot_k12, aes(x=inter_b, y=value, color=factor(variable)))+
  geom_point()+geom_line()+geom_text(data=plot_k12, aes(x=inter_b, y=value, label=value), color='black',
                                     size=3.5, hjust = 0, nudge_x = -0.8, check_overlap = TRUE)+
  facet_wrap(.~ country, nrow=2, ncol=4, scales = "fixed")+
  labs(x='The left end point of the interval', y='kappa', title='kappa', color='')+
  theme(legend.position = 'bottom')

VEs3=as.data.frame(kappa_estis)
VEs3[,-1]=1-apply(kappa_estis[,-1], 2, as.numeric)
colnames(VEs3)=c('province', 'partial', 'full', 'booster')
VEs3$period='Omicron-dominated'

country_VEs=rbind(#VEso1, VEso2, VEso3, VEs1, 
  VEs2, VEs3)
plot_VE=melt(country_VEs, id=c('province', 'period'), measure=c('partial', 'full', 'booster'))
plot_VE$period=factor(plot_VE$period, levels=c(
  # 'Pre-Delta', 'Intervening', 'Delta-dominated', "Pre-Omicron", 
  "Omicron-reported", "Omicron-dominated"))
write.csv(plot_VE, '~/Documents/COVID-19/data/plot_VE_theta1_reinfect.csv')
VE_countries_noreinfect=read.csv(paste0(path_table, 'VEs_postbooster_ad.csv'), stringsAsFactors = F)
plot_VE$doses='one'
plot_VE$doses[plot_VE$variable=='full']='two'
plot_VE$doses[plot_VE$variable=='booster']='three'
plot_VE$periods=paste(plot_VE$period, 'period')
VE_countries_reinfect=left_join(VE_countries_noreinfect, plot_VE[, c('province', 'value', 'doses', 'periods')])
VE_countries_reinfect$q975[!is.na(VE_countries_reinfect$value)]=sapply(VE_countries_reinfect$value[!is.na(VE_countries_reinfect$value)]+
                                                                         VE_countries_reinfect$q975[!is.na(VE_countries_reinfect$value)]-
                                                                         VE_countries_reinfect$original[!is.na(VE_countries_reinfect$value)], min, 0.99)
VE_countries_reinfect$q025[!is.na(VE_countries_reinfect$value)]=sapply(VE_countries_reinfect$value[!is.na(VE_countries_reinfect$value)]-
                                                                         VE_countries_reinfect$original[!is.na(VE_countries_reinfect$value)]+
                                                                         VE_countries_reinfect$q025[!is.na(VE_countries_reinfect$value)], max, 0.01)
VE_countries_reinfect$original[!is.na(VE_countries_reinfect$value)]=VE_countries_reinfect$value[!is.na(VE_countries_reinfect$value)]
write.csv(VE_countries_reinfect, paste0(path_table, 'VEs_postbooster_theta1_reinfect_ad.csv'))

plot_VE_theta1_reinfect=read.csv('~/Documents/COVID-19/data/plot_VE_theta1_reinfect.csv', stringsAsFactors = F)
plot_VE_theta1_reinfect$type='Reinfection'
plot_VE_theta1_reinfect_mur6=read.csv('~/Documents/COVID-19/data/plot_VE_theta1_reinfect_mur6.csv', stringsAsFactors = F)
plot_VE_theta1_reinfect_mur6$type='6 months'
plot_VE_theta1=read.csv('~/Documents/COVID-19/data/plot_VE_theta1.csv', stringsAsFactors = F)
plot_VE_theta1=plot_VE_theta1[plot_VE_theta1$period %in% c("Omicron-reported", "Omicron-dominated"),]
plot_VE_theta1$type='Main analysis'
plot_VE_theta2=read.csv('~/Documents/COVID-19/data/plot_VE_theta2.csv', stringsAsFactors = F)
plot_VE_theta2=plot_VE_theta2[plot_VE_theta2$period %in% c("Omicron-reported", "Omicron-dominated"),]
plot_VE_theta2$type='Asymptomatic rate' 
summary(arrange(plot_VE_theta1_reinfect, variable, period, province)$value-
          arrange(plot_VE_theta1_reinfect_mur6, variable, period, province)$value)
summary(abs(arrange(plot_VE_theta1_reinfect, variable, period, province)$value-
              arrange(plot_VE_theta1_reinfect_mur6, variable, period, province)$value))
summary(arrange(plot_VE_theta2, variable, period, province)$value-
          arrange(plot_VE_theta1_reinfect, variable, period, province)$value)
summary(abs(arrange(plot_VE_theta2, variable, period, province)$value-
              arrange(plot_VE_theta1_reinfect, variable, period, province)$value))
VEs_sens=rbind(plot_VE_theta1_reinfect[, c("province","period","variable","value","type")],
               plot_VE_theta1[, c("province","period","variable","value","type")],
               plot_VE_theta2[, c("province","period","variable","value","type")])
VEs_sens$period=factor(VEs_sens$period, levels=c("Omicron-reported","Omicron-dominated"),
                       labels=c('Intervening II', 'Omicron-dominated'))
VEs_sens$variable=factor(VEs_sens$variable, levels=c("partial","full","booster"))
ggplot(VEs_sens[VEs_sens$type=='Main analysis',])+
  geom_line(aes(x=period, y=value, group=variable, color=variable), linetype = 1)+
  geom_line(data=VEs_sens[VEs_sens$type == 'Asymptomatic rate',], aes(x=period, y=value, group=variable, color=variable), linetype=2)+
  geom_line(data=VEs_sens[VEs_sens$type == "Reinfection",], aes(x=period, y=value, group=variable, color=variable), linetype=3)+
  #geom_hline(yintercept = 0.5, linetype='dashed', color='gray')+
  facet_wrap(.~ province, nrow=2, scales = "fixed")+labs(x='Period', y='Vaccine efficacy', color='Vaccination')+ylim(0, 0.8)+
  theme_bw() + theme(axis.title = element_text(size = 10), 
                     text = element_text(face = "bold"),
                     strip.text = element_text(size = 10,face = 'bold'),
                     strip.background = element_rect(color="black", fill="white", linetype="solid"),
                     axis.title.x = element_text(size=10, face = 'bold', hjust = 0.5),
                     axis.title.y = element_text(size=10, face = 'bold', hjust = 0.5),
                     axis.text.x = element_text(size=10, angle=0, hjust = 0.5, face = 'bold'),
                     axis.text.y = element_text(size=10, face = 'bold'),
                     #axis.title.y.right=element_text(color='black'),
                     plot.title = element_text(size=15, face = 'bold', hjust = 0.5),
                     legend.title = element_text(size=10, face = 'bold'),
                     legend.text = element_text(size=10, face = 'bold'),
                     legend.key.width  = unit(.3,"inches"),
                     legend.key.height = unit(.3,"inches"))+theme(legend.position = c(0.9, 0.2))
ggsave(paste0(path_plot,'VEsensitivity.png'), units="in",width=12.3, height=6.5, dpi=300)

country_data$date=as.Date(country_data$date)
data_plot=country_data[country_data$t>14 & country_data$date < as.Date("2022-04-07")-14,]

ggplot(data_plot)+
  geom_line(aes(x =date, y=gamma.d.11))+
  geom_vline(data=time_point[, c('province', 'date')],
             aes(xintercept=date),color = "black",linetype = "dashed",size=0.7)+
  geom_rect(data=shade,
            mapping=aes(xmin=x1, xmax=x2, ymin = -Inf, ymax = Inf, fill=periods), alpha=0.4)+
  scale_fill_simpsons()+
  scale_x_date(date_breaks = "2 months",date_labels=format("%B %Y"))+
  #scale_x_date(date_breaks = "2 months",date_labels=format("%m/%d"))+
  facet_wrap(.~ province, nrow=2, scales = "free_y") + 
  mytheme+theme(legend.position = 'bottom')+
  labs(x = "Date", y = expression(bold(gamma['d,t'])))
ggsave(paste0(path_plot,'gammad_omicrond_reinfect.png'), units="in",width=13, height=6, dpi=300)

ggplot(data_plot)+
  geom_line(aes(x =date, y=gamma.r.11))+
  geom_vline(data=time_point[, c('province', 'date')],
             aes(xintercept=date),color = "black",linetype = "dashed",size=0.7)+
  geom_rect(data=shade,
            mapping=aes(xmin=x1, xmax=x2, ymin = -Inf, ymax = Inf, fill=periods), alpha=0.4)+
  scale_fill_simpsons()+
  scale_x_date(date_breaks = "2 months",date_labels=format("%B %Y"))+
  #scale_x_date(date_breaks = "2 months",date_labels=format("%m/%d"))+
  facet_wrap(.~ province, nrow=2, scales = "free_y") + 
  mytheme+theme(legend.position = 'bottom')+
  labs(x = "Date", y = expression(bold(gamma['r,t'])))
ggsave(paste0(path_plot,'gammar_omicrond_reinfect.png'), units="in",width=13, height=6, dpi=300)

country_SV=melt(data_plot, id=c('t1', 'province', 'date'), measure=c('S_infected', 'V1_infected', 'V2_infected', 'V3_infected'))
country_SV$variable=factor(country_SV$variable, levels=c('S_infected', 'V1_infected', 'V2_infected', 'V3_infected'),
                           labels=c('S_infected', 'V1_infected', 'V2_infected', 'V3_infected'))
ggplot(country_SV[country_SV$t1>=0,])+
  geom_line(aes(x =date, y=value, color=variable))+
  geom_vline(data=time_point[, c('province', 'date')],
             aes(xintercept=date),color = "black",linetype = "dashed",size=0.7)+
  geom_rect(data=shade,
            mapping=aes(xmin=x1, xmax=x2, ymin = -Inf, ymax = Inf, fill=periods), alpha=0.4)+
  scale_fill_simpsons()+
  scale_x_date(date_breaks = "2 months",date_labels=format("%B %Y"))+
  scale_color_manual(values=c('#CC0000','#0000FF',"#009966", "#9933CC"),
                     labels=c(expression(bold(S*"(t)")),expression(bold(V[1]*"(t)")),
                              expression(bold(V[2]*"(t)")),expression(bold(V[3]*"(t)"))))+
  facet_wrap(.~ province, nrow=2, scales = "free_y") + 
  labs(x = "Date", y = 'The number of daily infected people', title=expression(bold('The number of daily infected people from S(t), 
                                                                                    '*V[1]*"(t), "*V[2]*"(t) or "*V[3]*"(t)")))+
  # scale_x_date(date_breaks = "2 months",date_labels=format("%m/%d"))+
  mytheme+theme(legend.position = 'bottom')
ggsave(paste0(path_plot,'country_SV_omicrond_reinfect.png'), units="in",width=13, height=6, dpi=300)

ggplot(data_plot)+
  geom_line(aes(x =date, y=beta.mo.11))+
  geom_vline(data=time_point[, c('province', 'date')],
             aes(xintercept=date),color = "black",linetype = "dashed",size=0.7)+
  geom_rect(data=shade,
            mapping=aes(xmin=x1, xmax=x2, ymin = -Inf, ymax = Inf, fill=periods), alpha=0.4)+
  scale_fill_simpsons()+
  scale_x_date(date_breaks = "2 months",date_labels=format("%B %Y"))+
  facet_wrap(.~ province, nrow=2, scales = "fixed") + 
  # ylim(0, 2)+
  labs(x = "Date", y = expression(bold(beta[t])), title=expression(bold(beta[t])))+
  #scale_x_date(date_breaks = "2 months",date_labels=format("%m/%d"))+
  mytheme+theme(legend.position = 'bottom')
ggsave(paste0(path_plot,'country_beta_omicrond_reinfect.png'), units="in",width=13, height=6, dpi=300)

#-------Fitting performance-----------
B=1000
RelativeD_countries=RelativeD_seg_countries=c()
for(j1 in c(4, 8, 13, 18, 25, 26, 27))
{
  country=countries[j1] #country='US'
  print(country)
  
  estimates_final=read.csv(paste0(path_table,country,'_estimates_Omicron_dominate_reinfect.csv'), stringsAsFactors = F)
  nrow_dat=nrow(estimates_final)
  t_nV=min(which(estimates_final$G1>0))-1
  t_delta=which(as.character(as.Date(estimates_final$date))==as.character(
    time_variants$Earliest.report[time_variants$province==country&time_variants$Variant=='B.1.617.2']))-1
  t_delta_d=which(as.character(as.Date(estimates_final$date))==as.character(as.Date(
    Delta_period$dominate[Delta_period$province==country])))-1
  t_booster=min(which(estimates_final$G3>0))-1
  t_omicron=which(as.character(as.Date(estimates_final$date))==as.character(as.Date(
    Omicron_period$Omicron_b[Omicron_period$province==country])))-1
  t_omicron_d=which(as.character(as.Date(estimates_final$date))==as.character(as.Date(
    Omicron_period$Omicron_d[Omicron_period$province==country])))-1 
  
  M=estimates_final$population[1] 
  w0=provinces$w[provinces$provinces.en==country]
  
  ind.b=t_nV-40
  if(j1==13)
    ind.b=t_nV-50
  if(j1==25)
    ind.b=t_nV-70
  if(j1==18)
    ind.b=t_nV-180
  
  dat_0=estimates_final[1:(ind.b-1), c('t', 'S.modify', 'V1.modify', 'V2.modify', 'Ea.hat', 'E.hat', 'I', 
                                       'Ra.hat', 'R.d', 'R.r', 'G1', 'G2', 'R', 'N', 'V3.modify', 'G3')]
  colnames(dat_0)=c('t', 'S', 'V1', 'V2', 'Ea', 'Ep', 'I', 'R.a', 
                    'R.d', 'R.r', 'G1', 'G2', 'R', 'N', 'V3', 'G3')
  y=as.numeric(estimates_final[ind.b, c('V1.modify', 'V2.modify', 'Ea.hat', 'E.hat', 'I', 'R.d', 'R.r', 'Ra.hat', 'V3.modify', 'N')])
  beta_sim=estimates_final$beta.mo.11
  beta_sim[1]=beta_sim[2]; beta_sim[nrow_dat]=beta_sim[nrow_dat-1]
  gamma_d_sim=estimates_final$gamma.d.11
  gamma_d_sim[1]=gamma_d_sim[2]
  gamma_r_sim=estimates_final$gamma.r.11
  gamma_r_sim[1]=gamma_r_sim[2]
  
  simulations=array(0, dim=c(nrow_dat,4, B))
  simulations_seg=array(0, dim=c(nrow_dat,4,B))
  set.seed(12)
  for(hh in 1:B)
  {
    # Fitting performance of the whole period
    dat_sim1=sim_vsveipdr_Gs_booster_reinfect(t_nV-ind.b, y,  
                                              estimates_final$alpha[ind.b:t_nV], beta_sim[ind.b:t_nV],
                                              gamma_d_sim[ind.b:t_nV], gamma_r_sim[ind.b:t_nV],
                                              estimates_final$theta[ind.b:t_nV], M, 
                                              rep(0, t_nV-ind.b), rep(0, t_nV-ind.b), rep(0, t_nV-ind.b),
                                              estimates_final$mu_1[ind.b], estimates_final$mu_2[ind.b], estimates_final$mu_3[ind.b], estimates_final$mu_r[ind.b],
                                              0, 0, 0, r.beta)
    dat_sim1[,1]=dat_sim1[,1]+ind.b-1
    
    dat_sim2=sim_vsveipdr_Gs_booster_reinfect(t_delta-t_nV, dat_sim1[nrow(dat_sim1), c(3,4,5,6,7,9,10,8,15,14)], 
                                              estimates_final$alpha[t_nV:t_delta], beta_sim[t_nV:t_delta], 
                                              gamma_d_sim[t_nV:t_delta], gamma_r_sim[t_nV:t_delta],
                                              estimates_final$theta[t_nV:t_delta], M, 
                                              diff(estimates_final$G1[t_nV:t_delta]), 
                                              diff(estimates_final$G2[t_nV:t_delta]), 
                                              diff(estimates_final$G3[t_nV:t_delta]), 
                                              estimates_final$mu_1[t_nV], estimates_final$mu_2[t_nV], estimates_final$mu_3[t_nV], estimates_final$mu_r[t_nV], 
                                              estimates_final$kappa[t_nV+1], estimates_final$rho[t_nV+1], estimates_final$omega[t_nV+1], r.beta)
    dat_sim2[,1]=dat_sim2[,1]+t_nV-1
    dat_sim2[,11]=dat_sim2[,11]+dat_sim1[nrow(dat_sim1),11]
    dat_sim2[,12]=dat_sim2[,12]+dat_sim1[nrow(dat_sim1),12]
    dat_sim2[,16]=dat_sim2[,16]+dat_sim1[nrow(dat_sim1),16]
    
    dat_sim3=sim_vsveipdr_Gs_booster_reinfect(t_delta_d-t_delta, dat_sim2[nrow(dat_sim2), c(3,4,5,6,7,9,10,8,15,14)], 
                                              estimates_final$alpha[t_delta:t_delta_d], beta_sim[t_delta:t_delta_d], 
                                              gamma_d_sim[t_delta:t_delta_d], gamma_r_sim[t_delta:t_delta_d],
                                              estimates_final$theta[t_delta:t_delta_d], M, diff(estimates_final$G1[t_delta:t_delta_d]), 
                                              diff(estimates_final$G2[t_delta:t_delta_d]), diff(estimates_final$G3[t_delta:t_delta_d]), 
                                              estimates_final$mu_1[t_delta], estimates_final$mu_2[t_delta], estimates_final$mu_3[t_delta], estimates_final$mu_r[t_delta], 
                                              estimates_final$kappa[t_delta+1], estimates_final$rho[t_delta+1], estimates_final$omega[t_delta+1], r.beta)
    dat_sim3[,1]=dat_sim3[,1]+t_delta-1
    dat_sim3[,11]=dat_sim3[,11]+dat_sim2[nrow(dat_sim2),11]
    dat_sim3[,12]=dat_sim3[,12]+dat_sim2[nrow(dat_sim2),12]
    dat_sim3[,16]=dat_sim3[,16]+dat_sim2[nrow(dat_sim2),16]
    
    dat_sim4=sim_vsveipdr_Gs_booster_reinfect(t_booster-t_delta_d, dat_sim3[nrow(dat_sim3), c(3,4,5,6,7,9,10,8,15,14)], 
                                              estimates_final$alpha[t_delta_d:t_booster], beta_sim[t_delta_d:t_booster], 
                                              gamma_d_sim[t_delta_d:t_booster], gamma_r_sim[t_delta_d:t_booster],
                                              estimates_final$theta[t_delta_d:t_booster], M, diff(estimates_final$G1[t_delta_d:t_booster]), 
                                              diff(estimates_final$G2[t_delta_d:t_booster]), diff(estimates_final$G3[t_delta_d:t_booster]), 
                                              estimates_final$mu_1[t_delta_d], estimates_final$mu_2[t_delta_d], estimates_final$mu_3[t_delta_d], estimates_final$mu_r[t_delta_d],  
                                              estimates_final$kappa[t_delta_d+1], estimates_final$rho[t_delta_d+1], estimates_final$omega[t_delta_d+1], r.beta)
    dat_sim4[,1]=dat_sim4[,1]+t_delta_d-1
    dat_sim4[,11]=dat_sim4[,11]+dat_sim3[nrow(dat_sim3),11]
    dat_sim4[,12]=dat_sim4[,12]+dat_sim3[nrow(dat_sim3),12]
    dat_sim4[,16]=dat_sim4[,16]+dat_sim3[nrow(dat_sim3),16]
    
    dat_sim5=sim_vsveipdr_Gs_booster_reinfect(t_omicron-t_booster, dat_sim4[nrow(dat_sim4), c(3,4,5,6,7,9,10,8,15,14)], 
                                              estimates_final$alpha[t_booster:t_omicron], beta_sim[t_booster:t_omicron], 
                                              gamma_d_sim[t_booster:t_omicron], gamma_r_sim[t_booster:t_omicron],
                                              estimates_final$theta[t_booster:t_omicron], M, diff(estimates_final$G1[t_booster:t_omicron]), 
                                              diff(estimates_final$G2[t_booster:t_omicron]), diff(estimates_final$G3[t_booster:t_omicron]), 
                                              estimates_final$mu_1[t_booster], estimates_final$mu_2[t_booster], estimates_final$mu_3[t_booster], estimates_final$mu_r[t_booster],   
                                              estimates_final$kappa[t_booster+1], estimates_final$rho[t_booster+1], estimates_final$omega[t_booster+1], r.beta)
    dat_sim5[,1]=dat_sim5[,1]+t_booster-1
    dat_sim5[,11]=dat_sim5[,11]+dat_sim4[nrow(dat_sim4),11]
    dat_sim5[,12]=dat_sim5[,12]+dat_sim4[nrow(dat_sim4),12]
    dat_sim5[,16]=dat_sim5[,16]+dat_sim4[nrow(dat_sim4),16]
    
    dat_sim6=sim_vsveipdr_Gs_booster_reinfect(t_omicron_d-t_omicron, dat_sim5[nrow(dat_sim5), c(3,4,5,6,7,9,10,8,15,14)], 
                                              estimates_final$alpha[t_omicron:t_omicron_d], beta_sim[t_omicron:t_omicron_d], 
                                              gamma_d_sim[t_omicron:t_omicron_d], gamma_r_sim[t_omicron:t_omicron_d],
                                              estimates_final$theta[t_omicron:t_omicron_d], M, diff(estimates_final$G1[t_omicron:t_omicron_d]), 
                                              diff(estimates_final$G2[t_omicron:t_omicron_d]), diff(estimates_final$G3[t_omicron:t_omicron_d]), 
                                              estimates_final$mu_1[t_omicron], estimates_final$mu_2[t_omicron], estimates_final$mu_3[t_omicron], estimates_final$mu_r[t_omicron],   
                                              estimates_final$kappa[t_omicron+1], estimates_final$rho[t_omicron+1], estimates_final$omega[t_omicron+1], r.beta)
    dat_sim6[,1]=dat_sim6[,1]+t_omicron-1
    dat_sim6[,11]=dat_sim6[,11]+dat_sim5[nrow(dat_sim5),11]
    dat_sim6[,12]=dat_sim6[,12]+dat_sim5[nrow(dat_sim5),12]
    dat_sim6[,16]=dat_sim6[,16]+dat_sim5[nrow(dat_sim5),16]
    
    dat_sim7=sim_vsveipdr_Gs_booster_reinfect(nrow_dat-t_omicron_d, dat_sim6[nrow(dat_sim6), c(3,4,5,6,7,9,10,8,15,14)], 
                                              estimates_final$alpha[t_omicron_d:nrow_dat], beta_sim[t_omicron_d:nrow_dat], 
                                              gamma_d_sim[t_omicron_d:nrow_dat], gamma_r_sim[t_omicron_d:nrow_dat],
                                              estimates_final$theta[t_omicron_d:nrow_dat], M, diff(estimates_final$G1[t_omicron_d:nrow_dat]), 
                                              diff(estimates_final$G2[t_omicron_d:nrow_dat]), diff(estimates_final$G3[t_omicron_d:nrow_dat]), 
                                              estimates_final$mu_1[t_omicron_d], estimates_final$mu_2[t_omicron_d], estimates_final$mu_3[t_omicron_d], estimates_final$mu_r[t_omicron_d],   
                                              estimates_final$kappa[t_omicron_d+1], estimates_final$rho[t_omicron_d+1], estimates_final$omega[t_omicron_d+1], r.beta)
    dat_sim7[,1]=dat_sim7[,1]+t_omicron_d-1
    dat_sim7[,11]=dat_sim7[,11]+dat_sim6[nrow(dat_sim6),11]
    dat_sim7[,12]=dat_sim7[,12]+dat_sim6[nrow(dat_sim6),12]
    dat_sim7[,16]=dat_sim7[,16]+dat_sim6[nrow(dat_sim6),16]
    
    datas_sim=rbind(dat_sim1, dat_sim2[-1,], dat_sim3[-1,], dat_sim4[-1,], dat_sim5[-1,], dat_sim6[-1,], dat_sim7[-1,])
    datas_sim=as.data.frame(datas_sim)
    colnames(datas_sim)=c('t', 'S', 'V1', 'V2', 'Ea', 'Ep', 'I', 'R.a', 
                          'R.d', 'R.r', 'G1', 'G2', 'R', 'N', 'V3', 'G3')
    
    datas_sim=rbind(dat_0, datas_sim)
    datas_sim$dN=c(diff(datas_sim$N), NA)
    simulations[,,hh]=as.matrix(datas_sim[,c('t', 'I', 'N', 'dN')])
    
    # Fitting performance of separate periods
    dat_sim_seg1=sim_vsveipdr_Gs_booster_reinfect(t_nV-ind.b, y,  
                                                  estimates_final$alpha[ind.b:t_nV], beta_sim[ind.b:t_nV],
                                                  gamma_d_sim[ind.b:t_nV], gamma_r_sim[ind.b:t_nV],
                                                  estimates_final$theta[ind.b:t_nV], M, 
                                                  rep(0, t_nV-ind.b), rep(0, t_nV-ind.b), rep(0, t_nV-ind.b),
                                                  estimates_final$mu_1[ind.b], estimates_final$mu_2[ind.b], estimates_final$mu_3[ind.b], estimates_final$mu_r[ind.b],
                                                  0, 0, 0, r.beta)
    dat_sim_seg1[,1]=dat_sim_seg1[,1]+ind.b-1
    
    dat_sim_seg2=sim_vsveipdr_Gs_booster_reinfect(t_delta-t_nV, 
                                                  as.numeric(estimates_final[t_nV, c('V1.modify', 'V2.modify', 'Ea.hat', 'E.hat', 'I', 'R.d', 'R.r', 'Ra.hat', 'V3.modify', 'N')]),
                                                  estimates_final$alpha[t_nV:t_delta], beta_sim[t_nV:t_delta], 
                                                  gamma_d_sim[t_nV:t_delta], gamma_r_sim[t_nV:t_delta],
                                                  estimates_final$theta[t_nV:t_delta], M, 
                                                  diff(estimates_final$G1[t_nV:t_delta]), 
                                                  diff(estimates_final$G2[t_nV:t_delta]), 
                                                  diff(estimates_final$G3[t_nV:t_delta]), 
                                                  estimates_final$mu_1[t_nV], estimates_final$mu_2[t_nV], estimates_final$mu_3[t_nV], estimates_final$mu_r[t_nV], 
                                                  estimates_final$kappa[t_nV+1], estimates_final$rho[t_nV+1], estimates_final$omega[t_nV+1], r.beta)
    dat_sim_seg2[,1]=dat_sim_seg2[,1]+t_nV-1
    dat_sim_seg2[,11]=dat_sim_seg2[,11]+estimates_final$G1[t_nV] 
    dat_sim_seg2[,12]=dat_sim_seg2[,12]+estimates_final$G2[t_nV] 
    dat_sim_seg2[,16]=dat_sim_seg2[,16]+estimates_final$G3[t_nV] 
    
    dat_sim_seg3=sim_vsveipdr_Gs_booster_reinfect(t_delta_d-t_delta, 
                                                  as.numeric(estimates_final[t_delta, c('V1.modify', 'V2.modify', 'Ea.hat', 'E.hat', 'I', 'R.d', 'R.r', 'Ra.hat', 'V3.modify', 'N')]),
                                                  estimates_final$alpha[t_delta:t_delta_d], beta_sim[t_delta:t_delta_d], 
                                                  gamma_d_sim[t_delta:t_delta_d], gamma_r_sim[t_delta:t_delta_d],
                                                  estimates_final$theta[t_delta:t_delta_d], M, diff(estimates_final$G1[t_delta:t_delta_d]), 
                                                  diff(estimates_final$G2[t_delta:t_delta_d]), diff(estimates_final$G3[t_delta:t_delta_d]), 
                                                  estimates_final$mu_1[t_delta], estimates_final$mu_2[t_delta], estimates_final$mu_3[t_delta], estimates_final$mu_r[t_delta], 
                                                  estimates_final$kappa[t_delta+1], estimates_final$rho[t_delta+1], estimates_final$omega[t_delta+1], r.beta)
    dat_sim_seg3[,1]=dat_sim_seg3[,1]+t_delta-1
    dat_sim_seg3[,11]=dat_sim_seg3[,11]+estimates_final$G1[t_delta] 
    dat_sim_seg3[,12]=dat_sim_seg3[,12]+estimates_final$G2[t_delta] 
    dat_sim_seg3[,16]=dat_sim_seg3[,16]+estimates_final$G3[t_delta] 
    
    dat_sim_seg4=sim_vsveipdr_Gs_booster_reinfect(t_booster-t_delta_d, 
                                                  as.numeric(estimates_final[t_delta_d, c('V1.modify', 'V2.modify', 'Ea.hat', 'E.hat', 'I', 'R.d', 'R.r', 'Ra.hat', 'V3.modify', 'N')]),
                                                  estimates_final$alpha[t_delta_d:t_booster], beta_sim[t_delta_d:t_booster], 
                                                  gamma_d_sim[t_delta_d:t_booster], gamma_r_sim[t_delta_d:t_booster],
                                                  estimates_final$theta[t_delta_d:t_booster], M, diff(estimates_final$G1[t_delta_d:t_booster]), 
                                                  diff(estimates_final$G2[t_delta_d:t_booster]), diff(estimates_final$G3[t_delta_d:t_booster]), 
                                                  estimates_final$mu_1[t_delta_d], estimates_final$mu_2[t_delta_d], estimates_final$mu_3[t_delta_d], estimates_final$mu_r[t_delta_d],  
                                                  estimates_final$kappa[t_delta_d+1], estimates_final$rho[t_delta_d+1], estimates_final$omega[t_delta_d+1], r.beta)
    dat_sim_seg4[,1]=dat_sim_seg4[,1]+t_delta_d-1
    dat_sim_seg4[,11]=dat_sim_seg4[,11]+estimates_final$G1[t_delta_d] 
    dat_sim_seg4[,12]=dat_sim_seg4[,12]+estimates_final$G2[t_delta_d] 
    dat_sim_seg4[,16]=dat_sim_seg4[,16]+estimates_final$G3[t_delta_d] 
    
    dat_sim_seg5=sim_vsveipdr_Gs_booster_reinfect(t_omicron-t_booster, 
                                                  as.numeric(estimates_final[t_booster, c('V1.modify', 'V2.modify', 'Ea.hat', 'E.hat', 'I', 'R.d', 'R.r', 'Ra.hat', 'V3.modify', 'N')]),
                                                  estimates_final$alpha[t_booster:t_omicron], beta_sim[t_booster:t_omicron], 
                                                  gamma_d_sim[t_booster:t_omicron], gamma_r_sim[t_booster:t_omicron],
                                                  estimates_final$theta[t_booster:t_omicron], M, diff(estimates_final$G1[t_booster:t_omicron]), 
                                                  diff(estimates_final$G2[t_booster:t_omicron]), diff(estimates_final$G3[t_booster:t_omicron]), 
                                                  estimates_final$mu_1[t_booster], estimates_final$mu_2[t_booster], estimates_final$mu_3[t_booster], estimates_final$mu_r[t_booster],   
                                                  estimates_final$kappa[t_booster+1], estimates_final$rho[t_booster+1], estimates_final$omega[t_booster+1], r.beta)
    dat_sim_seg5[,1]=dat_sim_seg5[,1]+t_booster-1
    dat_sim_seg5[,11]=dat_sim_seg5[,11]+estimates_final$G1[t_booster] 
    dat_sim_seg5[,12]=dat_sim_seg5[,12]+estimates_final$G2[t_booster] 
    dat_sim_seg5[,16]=dat_sim_seg5[,16]+estimates_final$G3[t_booster] 
    
    dat_sim_seg6=sim_vsveipdr_Gs_booster_reinfect(t_omicron_d-t_omicron, 
                                                  as.numeric(estimates_final[t_omicron, c('V1.modify', 'V2.modify', 'Ea.hat', 'E.hat', 'I', 'R.d', 'R.r', 'Ra.hat', 'V3.modify', 'N')]),
                                                  estimates_final$alpha[t_omicron:t_omicron_d], beta_sim[t_omicron:t_omicron_d], 
                                                  gamma_d_sim[t_omicron:t_omicron_d], gamma_r_sim[t_omicron:t_omicron_d],
                                                  estimates_final$theta[t_omicron:t_omicron_d], M, diff(estimates_final$G1[t_omicron:t_omicron_d]), 
                                                  diff(estimates_final$G2[t_omicron:t_omicron_d]), diff(estimates_final$G3[t_omicron:t_omicron_d]), 
                                                  estimates_final$mu_1[t_omicron], estimates_final$mu_2[t_omicron], estimates_final$mu_3[t_omicron], estimates_final$mu_r[t_omicron],   
                                                  estimates_final$kappa[t_omicron+1], estimates_final$rho[t_omicron+1], estimates_final$omega[t_omicron+1], r.beta)
    dat_sim_seg6[,1]=dat_sim_seg6[,1]+t_omicron-1
    dat_sim_seg6[,11]=dat_sim_seg6[,11]+estimates_final$G1[t_omicron] 
    dat_sim_seg6[,12]=dat_sim_seg6[,12]+estimates_final$G2[t_omicron] 
    dat_sim_seg6[,16]=dat_sim_seg6[,16]+estimates_final$G3[t_omicron] 
    
    dat_sim_seg7=sim_vsveipdr_Gs_booster_reinfect(nrow_dat-t_omicron_d, 
                                                  as.numeric(estimates_final[t_omicron_d, c('V1.modify', 'V2.modify', 'Ea.hat', 'E.hat', 'I', 'R.d', 'R.r', 'Ra.hat', 'V3.modify', 'N')]),
                                                  estimates_final$alpha[t_omicron_d:nrow_dat], beta_sim[t_omicron_d:nrow_dat], 
                                                  gamma_d_sim[t_omicron_d:nrow_dat], gamma_r_sim[t_omicron_d:nrow_dat],
                                                  estimates_final$theta[t_omicron_d:nrow_dat], M, diff(estimates_final$G1[t_omicron_d:nrow_dat]), 
                                                  diff(estimates_final$G2[t_omicron_d:nrow_dat]), diff(estimates_final$G3[t_omicron_d:nrow_dat]), 
                                                  estimates_final$mu_1[t_omicron_d], estimates_final$mu_2[t_omicron_d], estimates_final$mu_3[t_omicron_d], estimates_final$mu_r[t_omicron_d],   
                                                  estimates_final$kappa[t_omicron_d+1], estimates_final$rho[t_omicron_d+1], estimates_final$omega[t_omicron_d+1], r.beta)
    dat_sim_seg7[,1]=dat_sim_seg7[,1]+t_omicron_d-1
    dat_sim_seg7[,11]=dat_sim_seg7[,11]+estimates_final$G1[t_omicron_d] 
    dat_sim_seg7[,12]=dat_sim_seg7[,12]+estimates_final$G2[t_omicron_d] 
    dat_sim_seg7[,16]=dat_sim_seg7[,16]+estimates_final$G3[t_omicron_d] 
    
    datas_sim_seg=rbind(dat_sim_seg1[-nrow(dat_sim_seg1),], dat_sim_seg2[-nrow(dat_sim_seg2),], dat_sim_seg3[-nrow(dat_sim_seg3),], 
                        dat_sim_seg4[-nrow(dat_sim_seg4),], dat_sim_seg5[-nrow(dat_sim_seg5),], dat_sim_seg6[-nrow(dat_sim_seg6),], dat_sim_seg7)
    datas_sim_seg=as.data.frame(datas_sim_seg)
    colnames(datas_sim_seg)=c('t', 'S', 'V1', 'V2', 'Ea', 'Ep', 'I', 'R.a', 
                              'R.d', 'R.r', 'G1', 'G2', 'R', 'N', 'V3', 'G3')
    
    datas_sim_seg=rbind(dat_0, datas_sim_seg)
    datas_sim_seg$dN=c(diff(datas_sim_seg$N), NA)
    simulations_seg[,,hh]=as.matrix(datas_sim_seg[,c('t', 'I', 'N', 'dN')])
  }
  simu_mean=apply(simulations,c(1,2),mean)
  RelativeD=simu_mean
  RelativeD[,-1]=RelativeD[,-1]/as.matrix(estimates_final[, c('I', 'N', 'dN')])-1
  RelativeD[is.nan(RelativeD)]=0
  colnames(RelativeD)=c('t', 'I', 'N', 'dN')
  RelativeD=as.data.frame(RelativeD)
  RelativeD$t1=RelativeD$t-t_nV
  RelativeD$t2=RelativeD$t-t_booster
  RelativeD$t3=RelativeD$t-t_omicron
  RelativeD$province=country
  RelativeD$date=estimates_final$date
  RelativeD$beta.mo.11=estimates_final$beta.mo.11
  RelativeD$scale_dN=estimates_final$scale_dN
  RelativeD$Rt.mo=estimates_final$Rt.mo
  RelativeD_countries=rbind(RelativeD_countries, RelativeD)
  
  simu_mean_seg=apply(simulations_seg,c(1,2),mean)
  RelativeD_seg=simu_mean_seg
  RelativeD_seg[,-1]=RelativeD_seg[,-1]/as.matrix(estimates_final[, c('I', 'N', 'dN')])-1
  RelativeD_seg[is.nan(RelativeD_seg)]=0
  colnames(RelativeD_seg)=c('t', 'I', 'N', 'dN')
  RelativeD_seg=as.data.frame(RelativeD_seg)
  RelativeD_seg$t1=RelativeD_seg$t-t_nV
  RelativeD_seg$t2=RelativeD_seg$t-t_booster
  RelativeD_seg$t3=RelativeD_seg$t-t_omicron
  RelativeD_seg$province=country
  RelativeD_seg$date=estimates_final$date
  RelativeD_seg$beta.mo.11=estimates_final$beta.mo.11
  RelativeD_seg$scale_dN=estimates_final$scale_dN
  RelativeD_seg$Rt.mo=estimates_final$Rt.mo
  RelativeD_seg$dN[c(t_nV, t_delta, t_delta_d, t_booster, t_omicron, t_omicron_d)-1]=NA
  RelativeD_seg_countries=rbind(RelativeD_seg_countries, RelativeD_seg)
}
RelativeD_countries$date=as.Date(RelativeD_countries$date)
RelativeD_countries$province[RelativeD_countries$province=='UK']="United Kingdom"
RelativeD_countries$province[RelativeD_countries$province=='US']="United States"
RelativeD_seg_countries$date=as.Date(RelativeD_seg_countries$date)
RelativeD_seg_countries$province[RelativeD_seg_countries$province=='UK']="United Kingdom"
RelativeD_seg_countries$province[RelativeD_seg_countries$province=='US']="United States"
as.data.frame(RelativeD_countries%>%group_by(province)%>% 
                summarise(RD_N=max(abs(N), na.rm=T), RD_dN=max(abs(dN), na.rm=T)))

as.data.frame(RelativeD_seg_countries%>%group_by(province)%>% 
                summarise(RD_N=max(abs(N), na.rm=T), RD_dN=max(abs(dN), na.rm=T)))
write.csv(RelativeD_countries, paste0(path_table, 'RelativeD_countries_theta1_reinfect.csv'))
write.csv(RelativeD_seg_countries, paste0(path_table, 'RelativeD_seg_countries_theta1_reinfect.csv'))

RelativeD_countries_reinfect=read.csv("~/Documents/COVID-19/data/real_data_V1_3/RelativeD_countries_theta1_reinfect.csv", stringsAsFactors = F)
RelativeD_countries_reinfect$type='Sensitivity analysis on reinfection'
RelativeD_countries_theta1=read.csv("~/Documents/COVID-19/data/real_data_V1_3/RelativeD_countries_theta1.csv", stringsAsFactors = F)
RelativeD_countries_theta1$type='Main analysis'
RelativeD_countries_theta2=read.csv("~/Documents/COVID-19/data/real_data_V1_31/RelativeD_countries_theta2.csv", stringsAsFactors = F)
RelativeD_countries_theta2$type='Sensitivity analysis on asymptomatic rate'
RelativeD_countries_theta2_reinfect=read.csv("~/Documents/COVID-19/data/real_data_V1_3/RelativeD_countries_theta2_reinfect.csv", stringsAsFactors = F)
RelativeD_countries_theta2_reinfect$type='Sensitivity analysis on asymptomatic rate'
RelativeD_countries_theta1_reinfect_mur6=read.csv("~/Documents/COVID-19/data/real_data_V1_3/RelativeD_countries_theta1_reinfect_mur6.csv", stringsAsFactors = F)
RelativeD_countries_theta1_reinfect_mur6$type='Sensitivity analysis on duration of natural immunity'
RelativeDs=rbind(RelativeD_countries_reinfect[,c('t', 'I', 'N', 'dN', 't1', 't2', 'province', 'date', 'type')],
                 RelativeD_countries_theta2_reinfect[,c('t', 'I', 'N', 'dN', 't1', 't2', 'province', 'date', 'type')],
                 RelativeD_countries_theta1_reinfect_mur6[,c('t', 'I', 'N', 'dN', 't1', 't2', 'province', 'date', 'type')])
RelativeDs$date=as.Date(RelativeDs$date)
ggplot(RelativeDs[RelativeDs$date<as.Date('2022-03-16')&
                    RelativeDs$type %in% c('Sensitivity analysis on duration of natural immunity', 
                                           'Sensitivity analysis on asymptomatic rate'),])+
  scale_color_manual(values=c('#CC0000','#0000FF'),
                     labels=c('Asymptomatic rate', 'Duration of natural immunity'))+
  geom_line(aes(x=date, y=N, color=type))+
  geom_line(data=RelativeDs[RelativeDs$date<as.Date('2022-03-16')&
                              RelativeDs$type == 'Sensitivity analysis on reinfection',],aes(x=date, y=N, color=type), color='black')+
  geom_hline(yintercept = 0, color='gray', linetype='dashed')+
  # geom_vline(data=time_point[time_point$province%in%c("Brazil","Germany","Italy","Peru","Turkey","United Kingdom","United States"), c('province', 'date')],
  #            aes(xintercept=date),color = "black",linetype = "dashed",size=0.6)+
  geom_rect(data=shades[shades$province%in%c("Brazil","Germany","Italy","Peru","Turkey","United Kingdom","United States"),],
            mapping=aes(xmin=x1, xmax=x2, ymin = -Inf, ymax = Inf, fill=Periods), alpha=0.4)+
  scale_fill_simpsons()+
  scale_x_date(date_breaks = "2 months",date_labels=format("%B %Y"))+
  facet_wrap(.~ province, nrow=2, scales = "fixed")+
  labs(x='Date', y='Relative error of numbers of cumulative confirmed cases', color='Sensitivity analyses')+theme_bw() + theme(axis.title = element_text(size = 15), 
                                                                                                 text = element_text(face = "bold"),
                                                                                                 strip.text = element_text(size = 12,face = 'bold'),
                                                                                                 strip.background = element_rect(color="black", fill="white", linetype="solid"),
                                                                                                 axis.title.x = element_text(size=15, face = 'bold', hjust = 0.5),
                                                                                                 axis.title.y = element_text(size=15, face = 'bold', hjust = 0.5),
                                                                                                 axis.text.x = element_text(size=13, angle=45,hjust = 1, face = 'bold'),
                                                                                                 axis.text.y = element_text(size=13, face = 'bold'),
                                                                                                 #axis.title.y.right=element_text(color='black'),
                                                                                                 #plot.title = element_text(size=15, face = 'bold', hjust = 0.5),
                                                                                                 plot.title = element_blank(),
                                                                                                 legend.title = element_text(size=14, face = 'bold'),
                                                                                                 legend.text = element_text(size=13, face = 'bold'),
                                                                                                 legend.key.width  = unit(.3,"inches"),
                                                                                                 legend.key.height = unit(.3,"inches"),
                                                                                                 panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + 
  theme(legend.position = c(0.9, 0.13))
ggsave(paste0(path_plot, 'RDN_Omicron_sensitivity_reinfect.png'),units="in",width=16, height=9, dpi=300)

#-----Booster vaccination scenario--------------
mytheme1=theme_bw() + theme(axis.title = element_text(size = 10), 
                            text = element_text(face = "bold"),
                            strip.text = element_text(size = 10,face = 'bold'),
                            strip.background = element_rect(color="black", fill="white", linetype="solid"),
                            axis.title.x = element_text(size=10, face = 'bold', hjust = 0.5),
                            axis.title.y = element_text(size=10, face = 'bold', hjust = 0.5),
                            axis.text.x = element_text(size=10, angle=45,hjust = 1, face = 'bold'),
                            axis.text.y = element_text(size=10, face = 'bold'),
                            #axis.title.y.right=element_text(color='black'),
                            plot.title = element_text(size=15, face = 'bold', hjust = 0.5),
                            legend.title = element_blank(),
                            legend.text = element_text(size=10, face = 'bold'),
                            legend.key.width  = unit(.3,"inches"),
                            legend.key.height = unit(.3,"inches"),
                            panel.grid.major=element_blank(), panel.grid.minor=element_blank())

B=100
multis=c(0, 0.5, 2)
multins=c('0', '0_5', '2')
for(mm in 1:length(multis))
{
  multi=multis[mm]
  multin=multins[mm]
  countries_boosters=c()
  for(j1 in c(4, 8, 13, 18, 25, 26, 27))
  {
    country=countries[j1] #country='US'
    print(country)
    estimates_final=read.csv(paste0(path_table,country,'_estimates_Omicron_dominate_reinfect.csv'), stringsAsFactors = F)
    nrow_dat=nrow(estimates_final)
    
    t_nV=min(which(estimates_final$G1>0))-1
    t_delta=which(as.character(as.Date(estimates_final$date))==as.character(
      time_variants$Earliest.report[time_variants$province==country&time_variants$Variant=='B.1.617.2']))-1
    t_delta_d=which(as.character(as.Date(estimates_final$date))==as.character(as.Date(
      Delta_period$dominate[Delta_period$province==country])))-1
    t_booster=min(which(estimates_final$G3>0))-1
    t_omicron=which(as.character(as.Date(estimates_final$date))==as.character(as.Date(
      Omicron_period$Omicron_b[Omicron_period$province==country])))-1
    t_omicron_d=which(as.character(as.Date(estimates_final$date))==as.character(as.Date(
      Omicron_period$Omicron_d[Omicron_period$province==country])))-1
    
    dat=estimates_final[t_booster:nrow_dat,]
    y=as.numeric(dat[1, c('V1.modify', 'V2.modify', 'Ea.hat', 'E.hat', 'I', 'R.d', 'R.r', 'Ra.hat', 'V3.modify', 'N')])
    M=dat$population[1] 
    dat$infectious = dat$I + dat$E.hat + dat$Ea.hat
    
    estimates_final$G3_scen=estimates_final$G3
    estimates_final$G3_scen[t_booster:nrow_dat]=cumsum(c(estimates_final$G3[t_booster], 
                                                         multi*diff(estimates_final$G3[t_booster:nrow_dat])))
    # plot(estimates_final$G2-estimates_final$G3_scen)
    estimates_final$G3_scen=apply(estimates_final[, c('G2', 'G3_scen')], 1, min)
    
    set.seed(12)
    simulations=array(0, dim=c(nrow_dat-t_booster+1, 20, B))
    for(hh in 1:B)
    {
      dat_sim5=sim_vsveipdr_Gs_booster_reinfect(t_omicron-t_booster, y, 
                                                estimates_final$alpha[t_booster:t_omicron], estimates_final$beta.mo.11[t_booster:t_omicron], 
                                                estimates_final$gamma.d.11[t_booster:t_omicron], estimates_final$gamma.r.11[t_booster:t_omicron],
                                                estimates_final$theta[t_booster:t_omicron], M, diff(estimates_final$G1[t_booster:t_omicron]), 
                                                diff(estimates_final$G2[t_booster:t_omicron]), diff(estimates_final$G3_scen[t_booster:t_omicron]), 
                                                estimates_final$mu_1[t_booster], estimates_final$mu_2[t_booster], estimates_final$mu_3[t_booster], estimates_final$mu_r[t_booster],   
                                                estimates_final$kappa[t_booster+1], estimates_final$rho[t_booster+1], estimates_final$omega[t_booster+1], r.beta)
      dat_sim5[,1]=dat_sim5[,1]+t_booster-1
      dat_sim5[,11]=dat_sim5[,11]+estimates_final$G1[t_booster]
      dat_sim5[,12]=dat_sim5[,12]+estimates_final$G2[t_booster]
      dat_sim5[,16]=dat_sim5[,16]+estimates_final$G3_scen[t_booster]
      
      dat_sim6=sim_vsveipdr_Gs_booster_reinfect(t_omicron_d-t_omicron, dat_sim5[nrow(dat_sim5), c(3,4,5,6,7,9,10,8,15,14)], 
                                                estimates_final$alpha[t_omicron:t_omicron_d], estimates_final$beta.mo.11[t_omicron:t_omicron_d], 
                                                estimates_final$gamma.d.11[t_omicron:t_omicron_d], estimates_final$gamma.r.11[t_omicron:t_omicron_d],
                                                estimates_final$theta[t_omicron:t_omicron_d], M, diff(estimates_final$G1[t_omicron:t_omicron_d]), 
                                                diff(estimates_final$G2[t_omicron:t_omicron_d]), diff(estimates_final$G3_scen[t_omicron:t_omicron_d]), 
                                                estimates_final$mu_1[t_omicron], estimates_final$mu_2[t_omicron], estimates_final$mu_3[t_omicron], estimates_final$mu_r[t_omicron],   
                                                estimates_final$kappa[t_omicron+1], estimates_final$rho[t_omicron+1], estimates_final$omega[t_omicron+1], r.beta)
      dat_sim6[,1]=dat_sim6[,1]+t_omicron-1
      dat_sim6[,11]=dat_sim6[,11]+dat_sim5[nrow(dat_sim5),11]
      dat_sim6[,12]=dat_sim6[,12]+dat_sim5[nrow(dat_sim5),12]
      dat_sim6[,16]=dat_sim6[,16]+dat_sim5[nrow(dat_sim5),16]
      
      dat_sim7=sim_vsveipdr_Gs_booster_reinfect(nrow_dat-t_omicron_d, dat_sim6[nrow(dat_sim6), c(3,4,5,6,7,9,10,8,15,14)], 
                                                estimates_final$alpha[t_omicron_d:nrow_dat], estimates_final$beta.mo.11[t_omicron_d:nrow_dat], 
                                                estimates_final$gamma.d.11[t_omicron_d:nrow_dat], estimates_final$gamma.r.11[t_omicron_d:nrow_dat],
                                                estimates_final$theta[t_omicron_d:nrow_dat], M, diff(estimates_final$G1[t_omicron_d:nrow_dat]), 
                                                diff(estimates_final$G2[t_omicron_d:nrow_dat]), diff(estimates_final$G3_scen[t_omicron_d:nrow_dat]), 
                                                estimates_final$mu_1[t_omicron_d], estimates_final$mu_2[t_omicron_d], estimates_final$mu_3[t_omicron_d], estimates_final$mu_r[t_omicron_d],   
                                                estimates_final$kappa[t_omicron_d+1], estimates_final$rho[t_omicron_d+1], estimates_final$omega[t_omicron_d+1], r.beta)
      dat_sim7[,1]=dat_sim7[,1]+t_omicron_d-1
      dat_sim7[,11]=dat_sim7[,11]+dat_sim6[nrow(dat_sim6),11]
      dat_sim7[,12]=dat_sim7[,12]+dat_sim6[nrow(dat_sim6),12]
      dat_sim7[,16]=dat_sim7[,16]+dat_sim6[nrow(dat_sim6),16]
      
      dat_simss = rbind(dat_sim5, dat_sim6[-1,], dat_sim7[-1,])
      dat_simss=cbind(dat_simss, c(diff(dat_simss[,14]), NA))
      dat_simss=cbind(dat_simss, c(diff(dat_simss[,9]), NA))
      dat_simss=cbind(dat_simss, c(diff(dat_simss[,10]), NA))
      dat_simss=cbind(dat_simss, apply(dat_simss[,5:7], 1, sum))
      
      simulations[,,hh]=dat_simss
    }
    
    simu_mean=apply(simulations,c(1,2),mean)
    simu_q75=apply(simulations,c(1,2),quantile,0.975, na.rm=T)
    simu_q25=apply(simulations,c(1,2),quantile,0.025, na.rm=T)
    colnames(simu_mean)=colnames(simu_q75)=colnames(simu_q25)=c('t', 'S', 'V1', 'V2', 'Ea', 'Ep', 'I', 'R.a', 
                                                                'R.d', 'R.r', 'G1', 'G2', 'R', 'N', 'V3', 'G3', 
                                                                'dN', 'dRd', 'dRr', 'infectious')
    simu_mean=as.data.frame(simu_mean)
    simu_q75=as.data.frame(simu_q75)
    simu_q25=as.data.frame(simu_q25)
    simu_mean$type='mean'
    simu_q75$type='97.5%'
    simu_q25$type='2.5%'
    
    simu_q25_1=simu_q25
    simu_q75_1=simu_q75
    colnames(simu_q25_1)[2:20]= paste0(colnames(simu_mean)[2:20],'.q25')
    colnames(simu_q75_1)[2:20]= paste0(colnames(simu_mean)[2:20],'.q75')
    dat_plot=cbind(simu_mean[1:20],simu_q25_1[paste0(colnames(simu_mean)[2:20],'.q25')],
                   simu_q75_1[paste0(colnames(simu_mean)[2:20],'.q75')])
    colnames(dat_plot)[2:20]=paste0(colnames(simu_mean)[2:20],'.mean')
    dat_plot$date=as.Date(dat$date[1]) + 0:(nrow_dat-t_booster)
    
    country_boosters=cbind(dat_plot, dat[, c("province", "population", "infected", "dead", "recovered", "R.d", "R.r",
                                             "N", "R", "I", "G1", "G2", "Ra.hat", "E.hat", "Ea.hat", "SV.hat", 
                                             "S.modify", "V1.modify", "V2.modify", "V3.modify", "G3", 
                                             "dN", "dRd", 'dRr', 'infectious', 't1', 't2')])
    country_boosters$dN_postbooster_mean=country_boosters$N.mean-country_boosters$N.mean[country_boosters$t2==0]
    country_boosters$dN_postbooster_q75=country_boosters$N.q75-country_boosters$N.mean[country_boosters$t2==0]
    country_boosters$dN_postbooster_q25=country_boosters$N.q25-country_boosters$N.mean[country_boosters$t2==0]
    country_boosters$dN_postbooster=country_boosters$N-country_boosters$N[country_boosters$t2==0]
    country_boosters$dRd_postbooster_mean=country_boosters$R.d.mean-country_boosters$R.d.mean[country_boosters$t2==0]
    country_boosters$dRd_postbooster_q75=country_boosters$R.d.q75-country_boosters$R.d.mean[country_boosters$t2==0]
    country_boosters$dRd_postbooster_q25=country_boosters$R.d.q25-country_boosters$R.d.mean[country_boosters$t2==0]
    country_boosters$dRd_postbooster=country_boosters$R.d-country_boosters$R.d[country_boosters$t2==0]
    country_boosters$dN_postvaccine_mean=country_boosters$N.mean-estimates_final$N[estimates_final$t1==0]
    country_boosters$dN_postvaccine_q75=country_boosters$N.q75-estimates_final$N[estimates_final$t1==0]
    country_boosters$dN_postvaccine_q25=country_boosters$N.q25-estimates_final$N[estimates_final$t1==0]
    country_boosters$dN_postvaccine=country_boosters$N-estimates_final$N[estimates_final$t1==0]
    country_boosters$dRd_postvaccine_mean=country_boosters$R.d.mean-estimates_final$R.d[estimates_final$t1==0]
    country_boosters$dRd_postvaccine_q75=country_boosters$R.d.q75-estimates_final$R.d[estimates_final$t1==0]
    country_boosters$dRd_postvaccine_q25=country_boosters$R.d.q25-estimates_final$R.d[estimates_final$t1==0]
    country_boosters$dRd_postvaccine=country_boosters$R.d-estimates_final$R.d[estimates_final$t1==0]
    countries_boosters=rbind(countries_boosters, country_boosters)
  }
  write.csv(countries_boosters, paste0('~/Documents/COVID-19/data/ND_',multin,'_booster_Omicron_reinfect.csv'))
  
}

#----No vaccination scenario--------------
countries_nova=c()
for(j1 in c(4, 8, 13, 18, 25, 26, 27))
{
  country=countries[j1] #country='US'
  print(country)
  estimates_final=read.csv(paste0(path_table,country,'_estimates_Omicron_dominate_reinfect.csv'), stringsAsFactors = F)
  nrow_dat=nrow(estimates_final)
  t_nV=min(which(estimates_final$G1>0))-1
  t_omicron=which(as.character(as.Date(estimates_final$date))==as.character(as.Date(
    Omicron_period$Omicron_b[Omicron_period$province==country])))-1
  dat=estimates_final[t_nV:nrow_dat,]
  y=as.numeric(dat[1, c('Ea.hat', 'E.hat', 'I', 'R.d', 'R.r', 'Ra.hat')])
  M=dat$population[1] 
  dat$infectious = dat$I + dat$E.hat + dat$Ea.hat
  
  set.seed(12)
  simulations=array(0, dim=c(nrow_dat-t_nV+1,20, B))
  for(hh in 1:B)
  {
    # as.vector(as.numeric(dat[1, c('V1.modify', 'V2.modify', 'Ea.hat', 'E.hat', 'I', 'R.d', 'R.r', 'Ra.hat', 'V3.modify')]))  
    dat_sim1=sim_vsveipdr_Gs_booster(t_omicron-t_nV, c(0,0,y,0), 
                                     estimates_final$alpha[t_nV:t_omicron], estimates_final$beta.mo.11[t_nV:t_omicron], 
                                     estimates_final$gamma.d.11[t_nV:t_omicron], estimates_final$gamma.r.11[t_nV:t_omicron], estimates_final$theta[t_nV:t_omicron], M, 
                                     dG1s=rep(0, t_omicron-t_nV), dG2s=rep(0, t_omicron-t_nV), 
                                     dG3s=rep(0, t_omicron-t_nV),
                                     0, 0, 0, 0, 0, 0, r.beta)
    dat_sim1[,1]=dat_sim1[,1]+t_nV-1
    
    dat_sim2=sim_vsveipdr_Gs_booster_reinfect(nrow_dat-t_omicron, dat_sim1[nrow(dat_sim1), c(3,4,5,6,7,9,10,8,15,14)], 
                                              estimates_final$alpha[t_omicron:nrow_dat], estimates_final$beta.mo.11[t_omicron:nrow_dat], 
                                              estimates_final$gamma.d.11[t_omicron:nrow_dat], estimates_final$gamma.r.11[t_omicron:nrow_dat],
                                              estimates_final$theta[t_omicron:nrow_dat], M, rep(0, nrow_dat-t_omicron),
                                              rep(0, nrow_dat-t_omicron), rep(0, nrow_dat-t_omicron), 
                                              0, 0, 0, estimates_final$mu_r[t_omicron],   
                                              0, 0, 0, r.beta)
    dat_sim2[,1]=dat_sim2[,1]+t_omicron-1
    
    dat_simss = rbind(dat_sim1, dat_sim2[-1,])
    dat_simss=cbind(dat_simss, c(diff(dat_simss[,14]), NA))
    dat_simss=cbind(dat_simss, c(diff(dat_simss[,9]), NA))
    dat_simss=cbind(dat_simss, c(diff(dat_simss[,10]), NA))
    dat_simss=cbind(dat_simss, apply(dat_simss[,5:7], 1, sum))
    
    simulations[,,hh]=dat_simss
  }
  
  simu_mean=apply(simulations,c(1,2),mean)
  simu_q75=apply(simulations,c(1,2),quantile,0.975, na.rm=T)
  simu_q25=apply(simulations,c(1,2),quantile,0.025, na.rm=T)
  colnames(simu_mean)=colnames(simu_q75)=colnames(simu_q25)=c('t', 'S', 'V1', 'V2', 'Ea', 'Ep', 'I', 'R.a', 
                                                              'R.d', 'R.r', 'G1', 'G2', 'R', 'N', 'V3', 'G3', 
                                                              'dN', 'dRd', 'dRr', 'infectious')
  simu_mean=as.data.frame(simu_mean)
  simu_q75=as.data.frame(simu_q75)
  simu_q25=as.data.frame(simu_q25)
  simu_mean$type='mean'
  simu_q75$type='97.5%'
  simu_q25$type='2.5%'
  
  simu_q25_1=simu_q25
  simu_q75_1=simu_q75
  colnames(simu_q25_1)[2:20]= paste0(colnames(simu_mean)[2:20],'.q25')
  colnames(simu_q75_1)[2:20]= paste0(colnames(simu_mean)[2:20],'.q75')
  dat_plot=cbind(simu_mean[1:20],simu_q25_1[paste0(colnames(simu_mean)[2:20],'.q25')],
                 simu_q75_1[paste0(colnames(simu_mean)[2:20],'.q75')])
  colnames(dat_plot)[2:20]=paste0(colnames(simu_mean)[2:20],'.mean')
  dat_plot$date=as.Date(dat$date)
  
  country_nova=cbind(dat_plot, dat[, c("province", "population", "infected", "dead", "recovered", "R.d", "R.r",
                                       "N", "R", "I", "G1", "G2", "Ra.hat", "E.hat", "Ea.hat", "SV.hat", 
                                       "S.modify", "V1.modify", "V2.modify", "V3.modify", "G3", 
                                       "dN", "dRd", 'dRr', 'infectious', 't1', 't2')])
  country_nova$dN_postbooster_mean=country_nova$N.mean-country_nova$N.mean[country_nova$t2==0]
  country_nova$dN_postbooster_q75=country_nova$N.q75-country_nova$N.mean[country_nova$t2==0]
  country_nova$dN_postbooster_q25=country_nova$N.q25-country_nova$N.mean[country_nova$t2==0]
  country_nova$dN_postbooster=country_nova$N-country_nova$N[country_nova$t2==0]
  country_nova$dRd_postbooster_mean=country_nova$R.d.mean-country_nova$R.d.mean[country_nova$t2==0]
  country_nova$dRd_postbooster_q75=country_nova$R.d.q75-country_nova$R.d.mean[country_nova$t2==0]
  country_nova$dRd_postbooster_q25=country_nova$R.d.q25-country_nova$R.d.mean[country_nova$t2==0]
  country_nova$dRd_postbooster=country_nova$R.d-country_nova$R.d[country_nova$t2==0]
  country_nova$dN_postvaccine_mean=country_nova$N.mean-country_nova$N.mean[country_nova$t1==0]
  country_nova$dN_postvaccine_q75=country_nova$N.q75-country_nova$N.mean[country_nova$t1==0]
  country_nova$dN_postvaccine_q25=country_nova$N.q25-country_nova$N.mean[country_nova$t1==0]
  country_nova$dN_postvaccine=country_nova$N-country_nova$N[country_nova$t1==0]
  country_nova$dRd_postvaccine_mean=country_nova$R.d.mean-country_nova$R.d.mean[country_nova$t1==0]
  country_nova$dRd_postvaccine_q75=country_nova$R.d.q75-country_nova$R.d.mean[country_nova$t1==0]
  country_nova$dRd_postvaccine_q25=country_nova$R.d.q25-country_nova$R.d.mean[country_nova$t1==0]
  country_nova$dRd_postvaccine=country_nova$R.d-country_nova$R.d[country_nova$t1==0]
  countries_nova=rbind(countries_nova, country_nova)
}
write.csv(countries_nova, '~/Documents/COVID-19/data/ND_novaccination_Omicron_reinfect.csv')

#-----Partial vaccination scenario--------------
countries_partialva=c()
for(j1 in c(4, 8, 13, 18, 25, 26, 27))
{
  country=countries[j1] #country='US'
  print(country)
  estimates_final=read.csv(paste0(path_table,country,'_estimates_Omicron_dominate_reinfect.csv'), stringsAsFactors = F)
  nrow_dat=nrow(estimates_final)
  t_nV=min(which(estimates_final$G1>0))-1
  t_delta=which(as.character(as.Date(estimates_final$date))==as.character(
    time_variants$Earliest.report[time_variants$province==country&time_variants$Variant=='B.1.617.2']))-1
  t_delta_d=which(as.character(as.Date(estimates_final$date))==as.character(as.Date(
    Delta_period$dominate[Delta_period$province==country])))-1
  t_booster=min(which(estimates_final$G3>0))-1
  t_omicron=which(as.character(as.Date(estimates_final$date))==as.character(as.Date(
    Omicron_period$Omicron_b[Omicron_period$province==country])))-1
  t_omicron_d=which(as.character(as.Date(estimates_final$date))==as.character(as.Date(
    Omicron_period$Omicron_d[Omicron_period$province==country])))-1
  
  dat=estimates_final[t_nV:nrow_dat,]
  y=as.numeric(dat[1, c('V1.modify', 'V2.modify', 'Ea.hat', 'E.hat', 'I', 'R.d', 'R.r', 'Ra.hat', 'V3.modify', 'N')])
  M=dat$population[1] 
  dat$infectious = dat$I + dat$E.hat + dat$Ea.hat
  
  set.seed(12)
  simulations=array(0, dim=c(nrow_dat-t_nV+1, 20, B))
  for(hh in 1:B)
  {
    dat_sim2=sim_vsveipdr_Gs_booster_reinfect(t_delta-t_nV, y, 
                                              estimates_final$alpha[t_nV:t_delta], estimates_final$beta.mo.11[t_nV:t_delta], 
                                              estimates_final$gamma.d.11[t_nV:t_delta], estimates_final$gamma.r.11[t_nV:t_delta],
                                              estimates_final$theta[t_nV:t_delta], M, 
                                              diff(estimates_final$G1[t_nV:t_delta]), 
                                              rep(0, t_delta-t_nV), rep(0, t_delta-t_nV),
                                              estimates_final$mu_1[t_nV], estimates_final$mu_2[t_nV], estimates_final$mu_3[t_nV], estimates_final$mu_r[t_nV],
                                              estimates_final$kappa[t_nV+1], estimates_final$rho[t_nV+1], 
                                              estimates_final$omega[t_nV+1], r.beta)
    dat_sim2[,1]=dat_sim2[,1]+t_nV-1
    
    dat_sim3=sim_vsveipdr_Gs_booster_reinfect(t_delta_d-t_delta, dat_sim2[nrow(dat_sim2), c(3,4,5,6,7,9,10,8,15,14)], 
                                              estimates_final$alpha[t_delta:t_delta_d], estimates_final$beta.mo.11[t_delta:t_delta_d], 
                                              estimates_final$gamma.d.11[t_delta:t_delta_d], estimates_final$gamma.r.11[t_delta:t_delta_d],
                                              estimates_final$theta[t_delta:t_delta_d], M, diff(estimates_final$G1[t_delta:t_delta_d]), 
                                              rep(0, t_delta_d-t_delta), rep(0, t_delta_d-t_delta), 
                                              estimates_final$mu_1[t_delta], estimates_final$mu_2[t_delta], estimates_final$mu_3[t_delta], estimates_final$mu_r[t_delta], 
                                              estimates_final$kappa[t_delta+1], estimates_final$rho[t_delta+1], 
                                              estimates_final$omega[t_delta+1], r.beta)
    dat_sim3[,1]=dat_sim3[,1]+t_delta-1
    dat_sim3[,11]=dat_sim3[,11]+dat_sim2[nrow(dat_sim2),11]
    dat_sim3[,12]=dat_sim3[,12]+dat_sim2[nrow(dat_sim2),12]
    dat_sim3[,16]=dat_sim3[,16]+dat_sim2[nrow(dat_sim2),16]
    
    dat_sim4=sim_vsveipdr_Gs_booster_reinfect(t_booster-t_delta_d, dat_sim3[nrow(dat_sim3), c(3,4,5,6,7,9,10,8,15,14)], 
                                              estimates_final$alpha[t_delta_d:t_booster], estimates_final$beta.mo.11[t_delta_d:t_booster], 
                                              estimates_final$gamma.d.11[t_delta_d:t_booster], estimates_final$gamma.r.11[t_delta_d:t_booster],
                                              estimates_final$theta[t_delta_d:t_booster], M, diff(estimates_final$G1[t_delta_d:t_booster]), 
                                              rep(0, t_booster-t_delta_d), rep(0, t_booster-t_delta_d), 
                                              estimates_final$mu_1[t_delta_d], estimates_final$mu_2[t_delta_d], estimates_final$mu_3[t_delta_d], estimates_final$mu_r[t_delta_d],
                                              estimates_final$kappa[t_delta_d+1], estimates_final$rho[t_delta_d+1], 
                                              estimates_final$omega[t_delta_d+1], r.beta)
    dat_sim4[,1]=dat_sim4[,1]+t_delta_d-1
    dat_sim4[,11]=dat_sim4[,11]+dat_sim3[nrow(dat_sim3),11]
    dat_sim4[,12]=dat_sim4[,12]+dat_sim3[nrow(dat_sim3),12]
    dat_sim4[,16]=dat_sim4[,16]+dat_sim3[nrow(dat_sim3),16]
    
    dat_sim5=sim_vsveipdr_Gs_booster_reinfect(t_omicron-t_booster, dat_sim4[nrow(dat_sim4), c(3,4,5,6,7,9,10,8,15,14)], 
                                              estimates_final$alpha[t_booster:t_omicron], estimates_final$beta.mo.11[t_booster:t_omicron], 
                                              estimates_final$gamma.d.11[t_booster:t_omicron], estimates_final$gamma.r.11[t_booster:t_omicron],
                                              estimates_final$theta[t_booster:t_omicron], M, diff(estimates_final$G1[t_booster:t_omicron]), 
                                              rep(0, t_omicron-t_booster), rep(0, t_omicron-t_booster), 
                                              estimates_final$mu_1[t_booster], estimates_final$mu_2[t_booster], estimates_final$mu_3[t_booster], estimates_final$mu_r[t_booster],
                                              estimates_final$kappa[t_booster+1], estimates_final$rho[t_booster+1], 
                                              estimates_final$omega[t_booster+1], r.beta)
    dat_sim5[,1]=dat_sim5[,1]+t_booster-1
    dat_sim5[,11]=dat_sim5[,11]+dat_sim4[nrow(dat_sim4),11]
    dat_sim5[,12]=dat_sim5[,12]+dat_sim4[nrow(dat_sim4),12]
    dat_sim5[,16]=dat_sim5[,16]+dat_sim4[nrow(dat_sim4),16]
    
    dat_sim6=sim_vsveipdr_Gs_booster_reinfect(t_omicron_d-t_omicron, dat_sim5[nrow(dat_sim5), c(3,4,5,6,7,9,10,8,15,14)], 
                                              estimates_final$alpha[t_omicron:t_omicron_d], estimates_final$beta.mo.11[t_omicron:t_omicron_d], 
                                              estimates_final$gamma.d.11[t_omicron:t_omicron_d], estimates_final$gamma.r.11[t_omicron:t_omicron_d],
                                              estimates_final$theta[t_omicron:t_omicron_d], M, diff(estimates_final$G1[t_omicron:t_omicron_d]), 
                                              rep(0, t_omicron_d-t_omicron), rep(0, t_omicron_d-t_omicron),  
                                              estimates_final$mu_1[t_omicron], estimates_final$mu_2[t_omicron], estimates_final$mu_3[t_omicron], estimates_final$mu_r[t_omicron],   
                                              estimates_final$kappa[t_omicron+1], estimates_final$rho[t_omicron+1], estimates_final$omega[t_omicron+1], r.beta)
    dat_sim6[,1]=dat_sim6[,1]+t_omicron-1
    dat_sim6[,11]=dat_sim6[,11]+dat_sim5[nrow(dat_sim5),11]
    dat_sim6[,12]=dat_sim6[,12]+dat_sim5[nrow(dat_sim5),12]
    dat_sim6[,16]=dat_sim6[,16]+dat_sim5[nrow(dat_sim5),16]
    
    dat_sim7=sim_vsveipdr_Gs_booster_reinfect(nrow_dat-t_omicron_d, dat_sim6[nrow(dat_sim6), c(3,4,5,6,7,9,10,8,15,14)], 
                                              estimates_final$alpha[t_omicron_d:nrow_dat], estimates_final$beta.mo.11[t_omicron_d:nrow_dat], 
                                              estimates_final$gamma.d.11[t_omicron_d:nrow_dat], estimates_final$gamma.r.11[t_omicron_d:nrow_dat],
                                              estimates_final$theta[t_omicron_d:nrow_dat], M, diff(estimates_final$G1[t_omicron_d:nrow_dat]), 
                                              rep(0, nrow_dat-t_omicron_d), rep(0, nrow_dat-t_omicron_d),
                                              estimates_final$mu_1[t_omicron_d], estimates_final$mu_2[t_omicron_d], estimates_final$mu_3[t_omicron_d], estimates_final$mu_r[t_omicron_d],   
                                              estimates_final$kappa[t_omicron_d+1], estimates_final$rho[t_omicron_d+1], estimates_final$omega[t_omicron_d+1], r.beta)
    dat_sim7[,1]=dat_sim7[,1]+t_omicron_d-1
    dat_sim7[,11]=dat_sim7[,11]+dat_sim6[nrow(dat_sim6),11]
    dat_sim7[,12]=dat_sim7[,12]+dat_sim6[nrow(dat_sim6),12]
    dat_sim7[,16]=dat_sim7[,16]+dat_sim6[nrow(dat_sim6),16]
    
    dat_simss = rbind(dat_sim2, dat_sim3[-1,], dat_sim4[-1,], dat_sim5[-1,], dat_sim6[-1,], dat_sim7[-1,])
    dat_simss=cbind(dat_simss, c(diff(dat_simss[,14]), NA))
    dat_simss=cbind(dat_simss, c(diff(dat_simss[,9]), NA))
    dat_simss=cbind(dat_simss, c(diff(dat_simss[,10]), NA))
    dat_simss=cbind(dat_simss, apply(dat_simss[,5:7], 1, sum))
    
    simulations[,,hh]=dat_simss
  }
  simu_mean=apply(simulations,c(1,2),mean)
  simu_q75=apply(simulations,c(1,2),quantile,0.975, na.rm=T)
  simu_q25=apply(simulations,c(1,2),quantile,0.025, na.rm=T)
  colnames(simu_mean)=colnames(simu_q75)=colnames(simu_q25)=c('t', 'S', 'V1', 'V2', 'Ea', 'Ep', 'I', 'R.a', 
                                                              'R.d', 'R.r', 'G1', 'G2', 'R', 'N', 'V3', 'G3', 
                                                              'dN', 'dRd', 'dRr', 'infectious')
  simu_mean=as.data.frame(simu_mean)
  simu_q75=as.data.frame(simu_q75)
  simu_q25=as.data.frame(simu_q25)
  simu_mean$type='mean'
  simu_q75$type='97.5%'
  simu_q25$type='2.5%'
  
  simu_q25_1=simu_q25
  simu_q75_1=simu_q75
  colnames(simu_q25_1)[2:20]= paste0(colnames(simu_mean)[2:20],'.q25')
  colnames(simu_q75_1)[2:20]= paste0(colnames(simu_mean)[2:20],'.q75')
  dat_plot=cbind(simu_mean[1:20],simu_q25_1[paste0(colnames(simu_mean)[2:20],'.q25')],
                 simu_q75_1[paste0(colnames(simu_mean)[2:20],'.q75')])
  colnames(dat_plot)[2:20]=paste0(colnames(simu_mean)[2:20],'.mean')
  dat_plot$date=as.Date(dat$date[1]) + 0:(nrow_dat-t_nV)
  
  country_partialva=cbind(dat_plot, dat[, c("province", "population", "infected", "dead", "recovered", "R.d", "R.r",
                                            "N", "R", "I", "G1", "G2", "Ra.hat", "E.hat", "Ea.hat", "SV.hat", 
                                            "S.modify", "V1.modify", "V2.modify", "V3.modify", "G3", 
                                            "dN", "dRd", 'dRr', 'infectious', 't1', 't2')])
  country_partialva$dN_postbooster_mean=country_partialva$N.mean-country_partialva$N.mean[country_partialva$t2==0]
  country_partialva$dN_postbooster_q75=country_partialva$N.q75-country_partialva$N.mean[country_partialva$t2==0]
  country_partialva$dN_postbooster_q25=country_partialva$N.q25-country_partialva$N.mean[country_partialva$t2==0]
  country_partialva$dN_postbooster=country_partialva$N-country_partialva$N[country_partialva$t2==0]
  country_partialva$dRd_postbooster_mean=country_partialva$R.d.mean-country_partialva$R.d.mean[country_partialva$t2==0]
  country_partialva$dRd_postbooster_q75=country_partialva$R.d.q75-country_partialva$R.d.mean[country_partialva$t2==0]
  country_partialva$dRd_postbooster_q25=country_partialva$R.d.q25-country_partialva$R.d.mean[country_partialva$t2==0]
  country_partialva$dRd_postbooster=country_partialva$R.d-country_partialva$R.d[country_partialva$t2==0]
  country_partialva$dN_postvaccine_mean=country_partialva$N.mean-country_partialva$N.mean[country_partialva$t1==0]
  country_partialva$dN_postvaccine_q75=country_partialva$N.q75-country_partialva$N.mean[country_partialva$t1==0]
  country_partialva$dN_postvaccine_q25=country_partialva$N.q25-country_partialva$N.mean[country_partialva$t1==0]
  country_partialva$dN_postvaccine=country_partialva$N-country_partialva$N[country_partialva$t1==0]
  country_partialva$dRd_postvaccine_mean=country_partialva$R.d.mean-country_partialva$R.d.mean[country_partialva$t1==0]
  country_partialva$dRd_postvaccine_q75=country_partialva$R.d.q75-country_partialva$R.d.mean[country_partialva$t1==0]
  country_partialva$dRd_postvaccine_q25=country_partialva$R.d.q25-country_partialva$R.d.mean[country_partialva$t1==0]
  country_partialva$dRd_postvaccine=country_partialva$R.d-country_partialva$R.d[country_partialva$t1==0]
  countries_partialva=rbind(countries_partialva, country_partialva)
}
write.csv(countries_partialva, '~/Documents/COVID-19/data/ND_partialvaccination_Omicron_reinfect.csv')

countries_nova_noreinfect=read.csv('~/Documents/COVID-19/data/ND_novaccination_Omicron_ad.csv', stringsAsFactors = F)
countries_nova=read.csv('~/Documents/COVID-19/data/ND_novaccination_Omicron_reinfect.csv', stringsAsFactors = F)
countries_nova$Rd.nova_Rd=countries_nova$R.d.mean/countries_nova$R.d
countries_nova$N.nova_N=countries_nova$N.mean/countries_nova$N
countries_nova$scenario='No vaccination'
countries_nova=left_join(countries_nova, countries_nova_noreinfect[,c("province", "date", "bCI_N", "bCI_I", "bCI_Ea", "bCI_Ep",
                                                                      "bCI_infectious", "bCI_R.d", "bCI_dN", "bCI_dRd")])
countries_nova[, c('N.q25', 'I.q25', 'Ea.q25', 'Ep.q25', 'infectious.q25', 'R.d.q25', 'dN.q25', 'dRd.q25')]=
  countries_nova[, c('N.mean', 'I.mean', 'Ea.mean', 'Ep.mean', 'infectious.mean', 'R.d.mean', 'dN.mean', 'dRd.mean')]-
  countries_nova[, c("bCI_N", "bCI_I", "bCI_Ea", "bCI_Ep", "bCI_infectious", "bCI_R.d", "bCI_dN", "bCI_dRd")]
countries_nova[, c('N.q75', 'I.q75', 'Ea.q75', 'Ep.q75', 'infectious.q75', 'R.d.q75', 'dN.q75', 'dRd.q75')]=
  countries_nova[, c('N.mean', 'I.mean', 'Ea.mean', 'Ep.mean', 'infectious.mean', 'R.d.mean', 'dN.mean', 'dRd.mean')]+
  countries_nova[, c("bCI_N", "bCI_I", "bCI_Ea", "bCI_Ep", "bCI_infectious", "bCI_R.d", "bCI_dN", "bCI_dRd")]
NRd_t_V=countries_nova[countries_nova$t1==0, c('province', 'N', 'R.d')]
colnames(NRd_t_V)=c('province', 'N_t_V', 'R.d_t_V')
NRd_t_booster=countries_nova[countries_nova$t2==0, c('province', 'N.mean', 'R.d.mean')]
colnames(NRd_t_booster)=c('province', 'N_t_booster', 'R.d_t_booster')
countries_nova=left_join(left_join(countries_nova, NRd_t_V), NRd_t_booster)
countries_nova$dN_postbooster_q25=countries_nova$N.q25-countries_nova$N_t_booster
countries_nova$dN_postbooster_q75=countries_nova$N.q75-countries_nova$N_t_booster
countries_nova$dRd_postbooster_q25=countries_nova$R.d.q25-countries_nova$R.d_t_booster
countries_nova$dRd_postbooster_q75=countries_nova$R.d.q75-countries_nova$R.d_t_booster
countries_nova$dN_postvaccine_q25=countries_nova$N.q25-countries_nova$N_t_V
countries_nova$dN_postvaccine_q75=countries_nova$N.q75-countries_nova$N_t_V
countries_nova$dRd_postvaccine_q25=countries_nova$R.d.q25-countries_nova$R.d_t_V
countries_nova$dRd_postvaccine_q75=countries_nova$R.d.q75-countries_nova$R.d_t_V
write.csv(countries_nova, '~/Documents/COVID-19/data/ND_novaccination_Omicron_reinfect_ad.csv')

countries_partialva_noreinfect=read.csv(paste0('~/Documents/COVID-19/data/ND_partialvaccination_Omicron_ad.csv'), stringsAsFactors = F)
countries_partialva=read.csv('~/Documents/COVID-19/data/ND_partialvaccination_Omicron_reinfect.csv', stringsAsFactors = F)
countries_partialva$Rd.nova_Rd=countries_partialva$R.d.mean/countries_partialva$R.d
countries_partialva$N.nova_N=countries_partialva$N.mean/countries_partialva$N
countries_partialva$scenario='Only the first dose'
countries_partialva=left_join(countries_partialva, countries_partialva_noreinfect[,c("province", "date", "bCI_N", "bCI_I", "bCI_Ea", "bCI_Ep",
                                                                                     "bCI_infectious", "bCI_R.d", "bCI_dN", "bCI_dRd")])
countries_partialva[, c('N.q25', 'I.q25', 'Ea.q25', 'Ep.q25', 'infectious.q25', 'R.d.q25', 'dN.q25', 'dRd.q25')]=
  countries_partialva[, c('N.mean', 'I.mean', 'Ea.mean', 'Ep.mean', 'infectious.mean', 'R.d.mean', 'dN.mean', 'dRd.mean')]-
  countries_partialva[, c("bCI_N", "bCI_I", "bCI_Ea", "bCI_Ep", "bCI_infectious", "bCI_R.d", "bCI_dN", "bCI_dRd")]
countries_partialva[, c('N.q75', 'I.q75', 'Ea.q75', 'Ep.q75', 'infectious.q75', 'R.d.q75', 'dN.q75', 'dRd.q75')]=
  countries_partialva[, c('N.mean', 'I.mean', 'Ea.mean', 'Ep.mean', 'infectious.mean', 'R.d.mean', 'dN.mean', 'dRd.mean')]+
  countries_partialva[, c("bCI_N", "bCI_I", "bCI_Ea", "bCI_Ep", "bCI_infectious", "bCI_R.d", "bCI_dN", "bCI_dRd")]
NRd_t_booster=countries_partialva[countries_partialva$t2==0, c('province', 'N.mean', 'R.d.mean')]
colnames(NRd_t_booster)=c('province', 'N_t_booster', 'R.d_t_booster')
countries_partialva=left_join(left_join(countries_partialva, NRd_t_V), NRd_t_booster)
countries_partialva$dN_postbooster_q25=countries_partialva$N.q25-countries_partialva$N_t_booster
countries_partialva$dN_postbooster_q75=countries_partialva$N.q75-countries_partialva$N_t_booster
countries_partialva$dRd_postbooster_q25=countries_partialva$R.d.q25-countries_partialva$R.d_t_booster
countries_partialva$dRd_postbooster_q75=countries_partialva$R.d.q75-countries_partialva$R.d_t_booster
countries_partialva$dN_postvaccine_q25=countries_partialva$N.q25-countries_partialva$N_t_V
countries_partialva$dN_postvaccine_q75=countries_partialva$N.q75-countries_partialva$N_t_V
countries_partialva$dRd_postvaccine_q25=countries_partialva$R.d.q25-countries_partialva$R.d_t_V
countries_partialva$dRd_postvaccine_q75=countries_partialva$R.d.q75-countries_partialva$R.d_t_V
write.csv(countries_partialva, '~/Documents/COVID-19/data/ND_partialvaccination_Omicron_reinfect_ad.csv')

countries_2boosterva_noreinfect=read.csv(paste0('~/Documents/COVID-19/data/ND_2_booster_Omicron_ad.csv'), stringsAsFactors = F)
countries_2boosterva=read.csv('~/Documents/COVID-19/data/ND_2_booster_Omicron_reinfect.csv', stringsAsFactors = F)
countries_2boosterva$Rd.nova_Rd=countries_2boosterva$R.d.mean/countries_2boosterva$R.d
countries_2boosterva$N.nova_N=countries_2boosterva$N.mean/countries_2boosterva$N
countries_2boosterva$scenario='2 times the boosters'
countries_2boosterva=left_join(countries_2boosterva, countries_2boosterva_noreinfect[,c("province", "date", "bCI_N", "bCI_I", "bCI_Ea", "bCI_Ep",
                                                                                        "bCI_infectious", "bCI_R.d", "bCI_dN", "bCI_dRd")])
countries_2boosterva[, c('N.q25', 'I.q25', 'Ea.q25', 'Ep.q25', 'infectious.q25', 'R.d.q25', 'dN.q25', 'dRd.q25')]=
  countries_2boosterva[, c('N.mean', 'I.mean', 'Ea.mean', 'Ep.mean', 'infectious.mean', 'R.d.mean', 'dN.mean', 'dRd.mean')]-
  countries_2boosterva[, c("bCI_N", "bCI_I", "bCI_Ea", "bCI_Ep", "bCI_infectious", "bCI_R.d", "bCI_dN", "bCI_dRd")]
countries_2boosterva[, c('N.q75', 'I.q75', 'Ea.q75', 'Ep.q75', 'infectious.q75', 'R.d.q75', 'dN.q75', 'dRd.q75')]=
  countries_2boosterva[, c('N.mean', 'I.mean', 'Ea.mean', 'Ep.mean', 'infectious.mean', 'R.d.mean', 'dN.mean', 'dRd.mean')]+
  countries_2boosterva[, c("bCI_N", "bCI_I", "bCI_Ea", "bCI_Ep", "bCI_infectious", "bCI_R.d", "bCI_dN", "bCI_dRd")]
NRd_t_booster=countries_2boosterva[countries_2boosterva$t2==0, c('province', 'N.mean', 'R.d.mean')]
colnames(NRd_t_booster)=c('province', 'N_t_booster', 'R.d_t_booster')
countries_2boosterva=left_join(left_join(countries_2boosterva, NRd_t_V), NRd_t_booster)
countries_2boosterva$dN_postbooster_q25=countries_2boosterva$N.q25-countries_2boosterva$N_t_booster
countries_2boosterva$dN_postbooster_q75=countries_2boosterva$N.q75-countries_2boosterva$N_t_booster
countries_2boosterva$dRd_postbooster_q25=countries_2boosterva$R.d.q25-countries_2boosterva$R.d_t_booster
countries_2boosterva$dRd_postbooster_q75=countries_2boosterva$R.d.q75-countries_2boosterva$R.d_t_booster
countries_2boosterva$dN_postvaccine_q25=countries_2boosterva$N.q25-countries_2boosterva$N_t_V
countries_2boosterva$dN_postvaccine_q75=countries_2boosterva$N.q75-countries_2boosterva$N_t_V
countries_2boosterva$dRd_postvaccine_q25=countries_2boosterva$R.d.q25-countries_2boosterva$R.d_t_V
countries_2boosterva$dRd_postvaccine_q75=countries_2boosterva$R.d.q75-countries_2boosterva$R.d_t_V
write.csv(countries_2boosterva, '~/Documents/COVID-19/data/ND_2_booster_Omicron_reinfect_ad.csv')

countries_05boosterva_noreinfect=read.csv(paste0('~/Documents/COVID-19/data/ND_0_5_booster_Omicron_ad.csv'), stringsAsFactors = F)
countries_05boosterva=read.csv('~/Documents/COVID-19/data/ND_0_5_booster_Omicron_reinfect.csv', stringsAsFactors = F)
countries_05boosterva$Rd.nova_Rd=countries_05boosterva$R.d.mean/countries_05boosterva$R.d
countries_05boosterva$N.nova_N=countries_05boosterva$N.mean/countries_05boosterva$N
countries_05boosterva$scenario='0.5 times the boosters'
countries_05boosterva=left_join(countries_05boosterva, countries_05boosterva_noreinfect[,c("province", "date", "bCI_N", "bCI_I", "bCI_Ea", "bCI_Ep",
                                                                                           "bCI_infectious", "bCI_R.d", "bCI_dN", "bCI_dRd")])
countries_05boosterva[, c('N.q25', 'I.q25', 'Ea.q25', 'Ep.q25', 'infectious.q25', 'R.d.q25', 'dN.q25', 'dRd.q25')]=
  countries_05boosterva[, c('N.mean', 'I.mean', 'Ea.mean', 'Ep.mean', 'infectious.mean', 'R.d.mean', 'dN.mean', 'dRd.mean')]-
  countries_05boosterva[, c("bCI_N", "bCI_I", "bCI_Ea", "bCI_Ep", "bCI_infectious", "bCI_R.d", "bCI_dN", "bCI_dRd")]
countries_05boosterva[, c('N.q75', 'I.q75', 'Ea.q75', 'Ep.q75', 'infectious.q75', 'R.d.q75', 'dN.q75', 'dRd.q75')]=
  countries_05boosterva[, c('N.mean', 'I.mean', 'Ea.mean', 'Ep.mean', 'infectious.mean', 'R.d.mean', 'dN.mean', 'dRd.mean')]+
  countries_05boosterva[, c("bCI_N", "bCI_I", "bCI_Ea", "bCI_Ep", "bCI_infectious", "bCI_R.d", "bCI_dN", "bCI_dRd")]
NRd_t_booster=countries_05boosterva[countries_05boosterva$t2==0, c('province', 'N.mean', 'R.d.mean')]
colnames(NRd_t_booster)=c('province', 'N_t_booster', 'R.d_t_booster')
countries_05boosterva=left_join(left_join(countries_05boosterva, NRd_t_V), NRd_t_booster)
countries_05boosterva$dN_postbooster_q25=countries_05boosterva$N.q25-countries_05boosterva$N_t_booster
countries_05boosterva$dN_postbooster_q75=countries_05boosterva$N.q75-countries_05boosterva$N_t_booster
countries_05boosterva$dRd_postbooster_q25=countries_05boosterva$R.d.q25-countries_05boosterva$R.d_t_booster
countries_05boosterva$dRd_postbooster_q75=countries_05boosterva$R.d.q75-countries_05boosterva$R.d_t_booster
countries_05boosterva$dN_postvaccine_q25=countries_05boosterva$N.q25-countries_05boosterva$N_t_V
countries_05boosterva$dN_postvaccine_q75=countries_05boosterva$N.q75-countries_05boosterva$N_t_V
countries_05boosterva$dRd_postvaccine_q25=countries_05boosterva$R.d.q25-countries_05boosterva$R.d_t_V
countries_05boosterva$dRd_postvaccine_q75=countries_05boosterva$R.d.q75-countries_05boosterva$R.d_t_V
write.csv(countries_05boosterva, '~/Documents/COVID-19/data/ND_0_5_booster_Omicron_reinfect_ad.csv')

countries_nobooster_noreinfect=read.csv(paste0('~/Documents/COVID-19/data/ND_nobooster_Omicron_ad.csv'), stringsAsFactors = F)
countries_nobooster=read.csv('~/Documents/COVID-19/data/ND_0_booster_Omicron_reinfect.csv', stringsAsFactors = F)
countries_nobooster$Rd.nova_Rd=countries_nobooster$R.d.mean/countries_nobooster$R.d
countries_nobooster$N.nova_N=countries_nobooster$N.mean/countries_nobooster$N
countries_nobooster$scenario='No booster'
countries_nobooster=left_join(countries_nobooster, countries_nobooster_noreinfect[,c("province", "date", "bCI_N", "bCI_I", "bCI_Ea", "bCI_Ep",
                                                                                     "bCI_infectious", "bCI_R.d", "bCI_dN", "bCI_dRd")])
countries_nobooster[, c('N.q25', 'I.q25', 'Ea.q25', 'Ep.q25', 'infectious.q25', 'R.d.q25', 'dN.q25', 'dRd.q25')]=
  countries_nobooster[, c('N.mean', 'I.mean', 'Ea.mean', 'Ep.mean', 'infectious.mean', 'R.d.mean', 'dN.mean', 'dRd.mean')]-
  countries_nobooster[, c("bCI_N", "bCI_I", "bCI_Ea", "bCI_Ep", "bCI_infectious", "bCI_R.d", "bCI_dN", "bCI_dRd")]
countries_nobooster[, c('N.q75', 'I.q75', 'Ea.q75', 'Ep.q75', 'infectious.q75', 'R.d.q75', 'dN.q75', 'dRd.q75')]=
  countries_nobooster[, c('N.mean', 'I.mean', 'Ea.mean', 'Ep.mean', 'infectious.mean', 'R.d.mean', 'dN.mean', 'dRd.mean')]+
  countries_nobooster[, c("bCI_N", "bCI_I", "bCI_Ea", "bCI_Ep", "bCI_infectious", "bCI_R.d", "bCI_dN", "bCI_dRd")]
NRd_t_booster=countries_nobooster[countries_nobooster$t2==0, c('province', 'N.mean', 'R.d.mean')]
colnames(NRd_t_booster)=c('province', 'N_t_booster', 'R.d_t_booster')
countries_nobooster=left_join(left_join(countries_nobooster, NRd_t_V), NRd_t_booster)
countries_nobooster$dN_postbooster_q25=countries_nobooster$N.q25-countries_nobooster$N_t_booster
countries_nobooster$dN_postbooster_q75=countries_nobooster$N.q75-countries_nobooster$N_t_booster
countries_nobooster$dRd_postbooster_q25=countries_nobooster$R.d.q25-countries_nobooster$R.d_t_booster
countries_nobooster$dRd_postbooster_q75=countries_nobooster$R.d.q75-countries_nobooster$R.d_t_booster
countries_nobooster$dN_postvaccine_q25=countries_nobooster$N.q25-countries_nobooster$N_t_V
countries_nobooster$dN_postvaccine_q75=countries_nobooster$N.q75-countries_nobooster$N_t_V
countries_nobooster$dRd_postvaccine_q25=countries_nobooster$R.d.q25-countries_nobooster$R.d_t_V
countries_nobooster$dRd_postvaccine_q75=countries_nobooster$R.d.q75-countries_nobooster$R.d_t_V
write.csv(countries_nobooster, '~/Documents/COVID-19/data/ND_0_booster_Omicron_reinfect_ad.csv')

plot_scenarios=rbind(countries_nova[, c('province', 'date', 'population', 't', 'N.mean', 'N.q25', 'N.q75', 'I.mean', 'I.q25', 'I.q75', 
                                        'Ea.mean', 'Ea.q25', 'Ea.q75', 'Ep.mean', 'Ep.q25', 'Ep.q75', 'V1.mean', 'V2.mean', 'V3.mean',
                                        'infectious.mean', 'infectious.q25', 'infectious.q75', 'S.mean', 'S.q25', 'S.q75', 'R.mean',
                                        'R.d.mean', 'R.d.q25', 'R.d.q75', 'dN.mean', 'dN.q25', 'dN.q75', 'R.a.mean', 'dRd', 'infectious',
                                        'dRd.mean', 'dRd.q25', 'dRd.q75', "G1.mean", "G2.mean", "G3.mean", 'scenario')],
                     countries_partialva[, c('province', 'date', 'population', 't', 'N.mean', 'N.q25', 'N.q75', 'I.mean', 'I.q25', 'I.q75', 
                                             'Ea.mean', 'Ea.q25', 'Ea.q75', 'Ep.mean', 'Ep.q25', 'Ep.q75', 'V1.mean', 'V2.mean', 'V3.mean',
                                             'infectious.mean', 'infectious.q25', 'infectious.q75', 'S.mean', 'S.q25', 'S.q75', 'R.mean',
                                             'R.d.mean', 'R.d.q25', 'R.d.q75', 'dN.mean', 'dN.q25', 'dN.q75', 'R.a.mean', 'dRd', 'infectious',
                                             'dRd.mean', 'dRd.q25', 'dRd.q75', "G1.mean", "G2.mean", "G3.mean", 'scenario')],
                     countries_nobooster[, c('province', 'date', 'population', 't', 'N.mean', 'N.q25', 'N.q75', 'I.mean', 'I.q25', 'I.q75', 
                                             'Ea.mean', 'Ea.q25', 'Ea.q75', 'Ep.mean', 'Ep.q25', 'Ep.q75', 'V1.mean', 'V2.mean', 'V3.mean',
                                             'infectious.mean', 'infectious.q25', 'infectious.q75', 'S.mean', 'S.q25', 'S.q75', 'R.mean',
                                             'R.d.mean', 'R.d.q25', 'R.d.q75', 'dN.mean', 'dN.q25', 'dN.q75', 'R.a.mean', 'dRd', 'infectious',
                                             'dRd.mean', 'dRd.q25', 'dRd.q75', "G1.mean", "G2.mean", "G3.mean", 'scenario')],
                     countries_05boosterva[, c('province', 'date', 'population', 't', 'N.mean', 'N.q25', 'N.q75', 'I.mean', 'I.q25', 'I.q75', 
                                               'Ea.mean', 'Ea.q25', 'Ea.q75', 'Ep.mean', 'Ep.q25', 'Ep.q75', 'V1.mean', 'V2.mean', 'V3.mean',
                                               'infectious.mean', 'infectious.q25', 'infectious.q75', 'S.mean', 'S.q25', 'S.q75', 'R.mean',
                                               'R.d.mean', 'R.d.q25', 'R.d.q75', 'dN.mean', 'dN.q25', 'dN.q75', 'R.a.mean', 'dRd', 'infectious',
                                               'dRd.mean', 'dRd.q25', 'dRd.q75', "G1.mean", "G2.mean", "G3.mean", 'scenario')],
                     countries_2boosterva[, c('province', 'date', 'population', 't', 'N.mean', 'N.q25', 'N.q75', 'I.mean', 'I.q25', 'I.q75', 
                                              'Ea.mean', 'Ea.q25', 'Ea.q75', 'Ep.mean', 'Ep.q25', 'Ep.q75', 'V1.mean', 'V2.mean', 'V3.mean',
                                              'infectious.mean', 'infectious.q25', 'infectious.q75', 'S.mean', 'S.q25', 'S.q75', 'R.mean',
                                              'R.d.mean', 'R.d.q25', 'R.d.q75', 'dN.mean', 'dN.q25', 'dN.q75', 'R.a.mean', 'dRd', 'infectious',
                                              'dRd.mean', 'dRd.q25', 'dRd.q75', "G1.mean", "G2.mean", "G3.mean", 'scenario')])

plot_scenarios$date=as.Date(plot_scenarios$date)
plot_scenarios$scenario=factor(plot_scenarios$scenario, levels = c("No vaccination", "Only the first dose",
                                                                   "No booster", "0.5 times the boosters", 
                                                                   "2 times the boosters"))
countries_nova$date=as.Date(countries_nova$date)
time_point = read.csv('~/Documents/COVID-19/data/time_points_Omicron.csv', stringsAsFactors = F)
time_point$date=as.Date(time_point$date)
time_point$Earliest.report=as.Date(time_point$Earliest.report)
time_point$dominate=as.Date(time_point$dominate)
time_point$booster_b=as.Date(time_point$booster_b)
time_point$Omicron_b=as.Date(time_point$Omicron_b)
time_point$Omicron_d=as.Date(time_point$Omicron_d)
shade=data.frame(province=rep(time_point$province, 6), x1=c(time_point$date, time_point$Earliest.report, time_point$dominate, 
                                                            time_point$booster_b, time_point$Omicron_b, time_point$Omicron_d), 
                 x2=c(time_point$Earliest.report, time_point$dominate, time_point$booster_b, 
                      time_point$Omicron_b, time_point$Omicron_d, rep('2022-04-07', 9)), 
                 periods=rep(paste0(c('Pre-Delta', 'Intervening', 'Delta-dominated', 
                                      'Pre-Omicron', 'Omicron', 'Omicron-dominated'), ' period'), each=9))
shade$x1=as.Date(shade$x1); shade$x2=as.Date(shade$x2)
shade$periods=factor(shade$periods, levels=paste0(c('Pre-Delta', 'Intervening', 'Delta-dominated', 
                                                    'Pre-Omicron', 'Omicron', 'Omicron-dominated'), ' period'))
shade$x1[shade$province=='Turkey'&shade$periods=='Pre-Omicron period']=
  shade$x1[shade$province=='Turkey'&shade$periods=='Delta-dominated period']
shade[shade$periods=='Omicron-dominated period', 'x2']=as.Date("2022-03-15")
shade=shade %>% rowwise %>% mutate(mean.date = mean.Date(c(x1, x2))) %>% data.frame()
shade$labels=rep(c('Pre-Delta', NA, 'Delta-\ndominated', 'Pre-Omicron', NA, 'Omicron-\ndominated'), each=9)
shade=shade[shade$province!='Turkey' | shade$periods!='Delta-dominated period',]
shade=shade[shade$province%in%c("Brazil","Germany","Italy","Peru","Turkey","United Kingdom","United States"), ]
time_point=time_point[time_point$province%in%c("Brazil","Germany","Italy","Peru","Turkey","United Kingdom","United States"), ]
y_position=as.data.frame(plot_scenarios[plot_scenarios$date < as.Date("2022-03-16"),]%>%group_by(province)%>%
                           summarise(y_posi=max(I.q75/10^6, na.rm = T)))
shade=left_join(shade, y_position)
shade$y_posi[shade$province%in% c('Brazil', 'Peru', 'Italy', 'Germany')&shade$periods=='Delta-dominated period']=
  0.77*shade$y_posi[shade$province%in% c('Brazil', 'Peru', 'Italy', 'Germany')&shade$periods=='Delta-dominated period']
shade$y_posi[shade$province%in% c('Brazil', 'Italy', "United Kingdom", "United States", 'Peru')&shade$periods=='Pre-Omicron period']=
  0.91*shade$y_posi[shade$province%in% c('Brazil', 'Italy', "United Kingdom", "United States", 'Peru')&shade$periods=='Pre-Omicron period']
shades=shade
shades$Periods=factor(shades$periods, levels=c('Pre-Delta period', 'Intervening period', 'Delta-dominated period', 
                                               'Pre-Omicron period', 'Omicron period', 'Omicron-dominated period'),
                      labels=c('Pre-Delta', 'Intervening I', 'Delta-dominated', 
                               'Pre-Omicron', 'Intervening II', 'Omicron-dominated'))
plot_D=plot_scenarios[plot_scenarios$date < as.Date("2022-03-16"),]
plot_D$Scenarios=plot_D$scenario
plot_D$Scenarios=factor(plot_D$Scenarios, levels=c('No vaccination', 'Only the first dose', 'No booster', 
                                                   '0.5 times the boosters', '2 times the boosters'),
                        label=c('No vaccination', 'Only partial vaccination', 'No booster',
                                '0.5 times the boosters', '2 times the boosters'))
library(ggnewscale)
# p3=
ggplot(plot_D) +
  geom_line(aes(x = date, y = I.mean/10^6, color=Scenarios),size=0.6)+
  geom_ribbon(aes(x = date, ymin =I.q25/10^6, ymax =I.q75/10^6, fill=Scenarios),linetype=0,alpha = 0.6) +
  scale_color_manual(values=c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF","black"))+ 
  scale_fill_manual(values=c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF","black"))+ 
  geom_line(data=countries_nova[countries_nova$date < as.Date("2022-03-16"),], aes(x=date, y=I/10^6),
            color='black', linetype='dashed', size=0.6)+
  # geom_vline(data=time_point[, c('province', 'date')],
  #            aes(xintercept=date),color = "black",linetype = "dashed",size=0.6)+
  new_scale("fill") +
  geom_rect(data=shades,
            mapping=aes(xmin=x1, xmax=x2, ymin = -Inf, ymax = Inf, fill=Periods), alpha=0.2)+
  scale_fill_simpsons()+
  geom_text(data=shades, aes(x=mean.date, y=1.05*y_posi, label=labels),
            check_overlap = F, size=3.6, fontface='bold')+
  facet_wrap(.~ province, nrow=2, ncol=4, scales = "free_y") + 
  scale_x_date(date_breaks = "2 months",date_labels=format("%B %Y"))+
  # scale_x_date(date_breaks = "2 months",date_labels=format("%m/%d"))+
  theme_bw() + theme(axis.title = element_text(size = 15), 
                     text = element_text(face = "bold"),
                     strip.text = element_text(size = 12,face = 'bold'),
                     strip.background = element_rect(color="black", fill="white", linetype="solid"),
                     axis.title.x = element_text(size=15, face = 'bold', hjust = 0.5),
                     axis.title.y = element_text(size=15, face = 'bold', hjust = 0.5),
                     axis.text.x = element_text(size=10.2, angle=45,hjust = 1, face = 'bold'),
                     axis.text.y = element_text(size=13, face = 'bold'),
                     #axis.title.y.right=element_text(color='black'),
                     #plot.title = element_text(size=15, face = 'bold', hjust = 0.5),
                     plot.title = element_blank(),
                     legend.title = element_text(size=12, face = 'bold'),
                     legend.text = element_text(size=11, face = 'bold'),
                     legend.key.width  = unit(.3,"inches"),
                     legend.key.height = unit(.3,"inches"),
                     panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + 
  theme(legend.position = c(0.9, 0.116))+
  labs(x = "Date", y = 'Active confirmed cases (million)', title='(b) Active confirmed cases')
ggsave(paste0(path_plot, 'D_scenarios_reinfect.pdf'),units="in",width=16.5, height=8.9, dpi=300)

dat_S=plot_D
dat_S$S.mean=dat_S$population-
  dat_S$V1.mean-dat_S$V2.mean-
  dat_S$V3.mean-dat_S$Ea.mean-
  dat_S$Ep.mean-dat_S$I.mean-
  dat_S$R.a.mean-dat_S$R.mean
dat_S$S.q25=dat_S$S.mean-(plot_D$S.mean-plot_D$S.q25)
dat_S$S.q75=dat_S$S.mean+(plot_D$S.q75-plot_D$S.mean)
ggplot(dat_S[dat_S$scenario %in% c('No vaccination', 'Only the first dose', 'No booster'), ]) +
  geom_line(aes(x = date, y = S.mean/10^6, color=Scenarios),size=0.6)+
  geom_ribbon(aes(x = date, ymin =S.q25/10^6, ymax =S.q75/10^6, fill=Scenarios),linetype=0,alpha = 0.6) +
  scale_color_manual(values=c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF","black"))+ 
  scale_fill_manual(values=c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF","black"))+ 
  geom_line(data=countries_nova[countries_nova$date < as.Date("2022-03-16"),], aes(x=date, y=S.modify/10^6),
            color='black', linetype='dashed', size=0.6)+
  geom_hline(yintercept = 0, color='white')+
  # geom_vline(data=time_point[, c('province', 'date')],
  #            aes(xintercept=date),color = "black",linetype = "dashed",size=0.6)+
  new_scale("fill") +
  geom_rect(data=shades,
            mapping=aes(xmin=x1, xmax=x2, ymin = -Inf, ymax = Inf, fill=Periods), alpha=0.2)+
  scale_fill_simpsons()+
  facet_wrap(.~ province, nrow=2, ncol=4, scales = "free_y") + 
  scale_x_date(date_breaks = "2 months",date_labels=format("%B %Y"))+
  # scale_x_date(date_breaks = "2 months",date_labels=format("%m/%d"))+
  theme_bw() + theme(axis.title = element_text(size = 15), 
                     text = element_text(face = "bold"),
                     strip.text = element_text(size = 12,face = 'bold'),
                     strip.background = element_rect(color="black", fill="white", linetype="solid"),
                     axis.title.x = element_text(size=15, face = 'bold', hjust = 0.5),
                     axis.title.y = element_text(size=15, face = 'bold', hjust = 0.5),
                     axis.text.x = element_text(size=10.2, angle=45,hjust = 1, face = 'bold'),
                     axis.text.y = element_text(size=13, face = 'bold'),
                     #axis.title.y.right=element_text(color='black'),
                     #plot.title = element_text(size=15, face = 'bold', hjust = 0.5),
                     plot.title = element_blank(),
                     legend.title = element_text(size=14, face = 'bold'),
                     legend.text = element_text(size=13, face = 'bold'),
                     legend.key.width  = unit(.3,"inches"),
                     legend.key.height = unit(.3,"inches"),
                     panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + 
  theme(legend.position = c(0.9, 0.14))+
  labs(x = "Date", y = 'Numbers of susceptible people without the vaccine immunity (million)', title='')
ggsave(paste0(path_plot, 'S_scenarios_reinfect.pdf'),units="in",width=16.5, height=9.5, dpi=300)

ggplot(plot_D) +
  geom_line(aes(x = date, y = dRd.mean/1000, color=Scenarios),size=0.6)+
  geom_ribbon(aes(x = date, ymin =dRd.q25/1000, ymax =dRd.q75/1000, fill=Scenarios),linetype=0,alpha = 0.6) +
  scale_color_manual(values=c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF","black"))+ 
  scale_fill_manual(values=c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF","black"))+ 
  geom_line(data=countries_nova[countries_nova$date < as.Date("2022-03-16"),], aes(x=date, y=dRd/1000),
            color='black', linetype='dashed', size=0.6)+
  geom_hline(yintercept = 0, color='white')+
  # geom_vline(data=time_point[, c('province', 'date')],
  #            aes(xintercept=date),color = "black",linetype = "dashed",size=0.6)+
  new_scale("fill") +
  geom_rect(data=shades,
            mapping=aes(xmin=x1, xmax=x2, ymin = -Inf, ymax = Inf, fill=Periods), alpha=0.2)+
  scale_fill_simpsons()+
  facet_wrap(.~ province, nrow=2, ncol=4, scales = "free_y") + 
  scale_x_date(date_breaks = "2 months",date_labels=format("%B %Y"))+
  # scale_x_date(date_breaks = "2 months",date_labels=format("%m/%d"))+
  theme_bw() + theme(axis.title = element_text(size = 15), 
                     text = element_text(face = "bold"),
                     strip.text = element_text(size = 12,face = 'bold'),
                     strip.background = element_rect(color="black", fill="white", linetype="solid"),
                     axis.title.x = element_text(size=15, face = 'bold', hjust = 0.5),
                     axis.title.y = element_text(size=15, face = 'bold', hjust = 0.5),
                     axis.text.x = element_text(size=10.2, angle=45,hjust = 1, face = 'bold'),
                     axis.text.y = element_text(size=13, face = 'bold'),
                     #axis.title.y.right=element_text(color='black'),
                     #plot.title = element_text(size=15, face = 'bold', hjust = 0.5),
                     plot.title = element_blank(),
                     legend.title = element_text(size=14, face = 'bold'),
                     legend.text = element_text(size=13, face = 'bold'),
                     legend.key.width  = unit(.3,"inches"),
                     legend.key.height = unit(.3,"inches"),
                     panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + 
  theme(legend.position = c(0.9, 0.14))+
  labs(x = "Date", y = 'Numbers of Daily new deaths (thousand)', title='')
ggsave(paste0(path_plot, 'dRd_scenarios_reinfect.pdf'),units="in",width=16.5, height=9.5, dpi=300)

countries_scenarios=rbind(countries_nova[, c('province', 'date', 't2', 'R.d', 'R.d.mean', 'R.d.q25', 'R.d.q75', 
                                             'N', 'N.mean', 'N.q25', 'N.q75', 
                                             'dN_postbooster_mean', 'dN_postbooster_q75', 'dN_postbooster_q25', 'dN_postbooster',
                                             'dRd_postbooster_mean', 'dRd_postbooster_q75', 'dRd_postbooster_q25', 'dRd_postbooster',
                                             'dN_postvaccine_mean', 'dN_postvaccine_q75', 'dN_postvaccine_q25', 'dN_postvaccine',
                                             'dRd_postvaccine_mean', 'dRd_postvaccine_q75', 'dRd_postvaccine_q25', 'dRd_postvaccine',
                                             'Rd.nova_Rd', 'N.nova_N', 
                                             "G1.mean", "G2.mean", "G3.mean", 'scenario')],
                          countries_partialva[, c('province', 'date', 't2', 'R.d', 'R.d.mean', 'R.d.q25', 'R.d.q75', 
                                                  'N', 'N.mean', 'N.q25', 'N.q75', 
                                                  'dN_postbooster_mean', 'dN_postbooster_q75', 'dN_postbooster_q25', 'dN_postbooster',
                                                  'dRd_postbooster_mean', 'dRd_postbooster_q75', 'dRd_postbooster_q25', 'dRd_postbooster',
                                                  'dN_postvaccine_mean', 'dN_postvaccine_q75', 'dN_postvaccine_q25', 'dN_postvaccine',
                                                  'dRd_postvaccine_mean', 'dRd_postvaccine_q75', 'dRd_postvaccine_q25', 'dRd_postvaccine',
                                                  'Rd.nova_Rd', 'N.nova_N', 
                                                  "G1.mean", "G2.mean", "G3.mean", 'scenario')],
                          countries_nobooster[, c('province', 'date', 't2', 'R.d', 'R.d.mean', 'R.d.q25', 'R.d.q75', 
                                                  'N', 'N.mean', 'N.q25', 'N.q75', 
                                                  'dN_postbooster_mean', 'dN_postbooster_q75', 'dN_postbooster_q25', 'dN_postbooster',
                                                  'dRd_postbooster_mean', 'dRd_postbooster_q75', 'dRd_postbooster_q25', 'dRd_postbooster',
                                                  'dN_postvaccine_mean', 'dN_postvaccine_q75', 'dN_postvaccine_q25', 'dN_postvaccine',
                                                  'dRd_postvaccine_mean', 'dRd_postvaccine_q75', 'dRd_postvaccine_q25', 'dRd_postvaccine',
                                                  'Rd.nova_Rd', 'N.nova_N', 
                                                  "G1.mean", "G2.mean", "G3.mean", 'scenario')],
                          countries_05boosterva[, c('province', 'date', 't2', 'R.d', 'R.d.mean', 'R.d.q25', 'R.d.q75', 
                                                    'N', 'N.mean', 'N.q25', 'N.q75', 
                                                    'dN_postbooster_mean', 'dN_postbooster_q75', 'dN_postbooster_q25', 'dN_postbooster',
                                                    'dRd_postbooster_mean', 'dRd_postbooster_q75', 'dRd_postbooster_q25', 'dRd_postbooster',
                                                    'dN_postvaccine_mean', 'dN_postvaccine_q75', 'dN_postvaccine_q25', 'dN_postvaccine',
                                                    'dRd_postvaccine_mean', 'dRd_postvaccine_q75', 'dRd_postvaccine_q25', 'dRd_postvaccine',
                                                    'Rd.nova_Rd', 'N.nova_N', 
                                                    "G1.mean", "G2.mean", "G3.mean", 'scenario')],
                          countries_2boosterva[, c('province', 'date', 't2', 'R.d', 'R.d.mean', 'R.d.q25', 'R.d.q75', 
                                                   'N', 'N.mean', 'N.q25', 'N.q75', 
                                                   'dN_postbooster_mean', 'dN_postbooster_q75', 'dN_postbooster_q25', 'dN_postbooster',
                                                   'dRd_postbooster_mean', 'dRd_postbooster_q75', 'dRd_postbooster_q25', 'dRd_postbooster',
                                                   'dN_postvaccine_mean', 'dN_postvaccine_q75', 'dN_postvaccine_q25', 'dN_postvaccine',
                                                   'dRd_postvaccine_mean', 'dRd_postvaccine_q75', 'dRd_postvaccine_q25', 'dRd_postvaccine',
                                                   'Rd.nova_Rd', 'N.nova_N', 
                                                   "G1.mean", "G2.mean", "G3.mean", 'scenario')])
increase_ND_boot=countries_scenarios[countries_scenarios$date==as.Date("2022-03-15"),
                                     c('province', 'R.d', 'R.d.mean', 'R.d.q25', 'R.d.q75', 
                                       'N', 'N.mean', 'N.q25', 'N.q75', 
                                       'dN_postbooster_mean', 'dN_postbooster_q75', 'dN_postbooster_q25', 'dN_postbooster',
                                       'dRd_postbooster_mean', 'dRd_postbooster_q75', 'dRd_postbooster_q25', 'dRd_postbooster',
                                       'dN_postvaccine_mean', 'dN_postvaccine_q75', 'dN_postvaccine_q25', 'dN_postvaccine',
                                       'dRd_postvaccine_mean', 'dRd_postvaccine_q75', 'dRd_postvaccine_q25', 'dRd_postvaccine',
                                       'Rd.nova_Rd', 'N.nova_N', 'scenario')]
increase_ND_boot$Rd.increase=increase_ND_boot$R.d.mean-increase_ND_boot$R.d
increase_ND_boot$N.increase=increase_ND_boot$N.mean-increase_ND_boot$N
increase_ND_boot$dN_postbooster_percentage=increase_ND_boot$dN_postbooster_mean/increase_ND_boot$dN_postbooster
increase_ND_boot$dRd_postbooster_percentage=increase_ND_boot$dRd_postbooster_mean/increase_ND_boot$dRd_postbooster
write.csv(increase_ND_boot[increase_ND_boot$scenario%in% c("No booster", "0.5 times the boosters","2 times the boosters"),], 
          '~/Documents/COVID-19/data/scenario_results_boots_reinfect.csv')

tb_dNRd=increase_ND_boot[increase_ND_boot$scenario%in% c("No booster", "0.5 times the boosters","2 times the boosters"),
                         c('province', 'dN_postbooster', 'dN_postbooster_mean', 'dN_postbooster_q75', 'dN_postbooster_q25',  
                           'dRd_postbooster', 'dRd_postbooster_mean', 'dRd_postbooster_q75', 'dRd_postbooster_q25', 'scenario')]
tb_dNRd$pr_dN_postbooster=paste0(round(tb_dNRd$dN_postbooster_mean/1000), ' (', round(tb_dNRd$dN_postbooster_q25/1000), 
                                 ', ', round(tb_dNRd$dN_postbooster_q75/1000), ')')
tb_dNRd$pr_dN_postbooster_percentage=paste0(round(tb_dNRd$dN_postbooster_mean/tb_dNRd$dN_postbooster, 2)*100, 
                                            ' (', round(tb_dNRd$dN_postbooster_q25/tb_dNRd$dN_postbooster, 2)*100, 
                                            ', ', round(tb_dNRd$dN_postbooster_q75/tb_dNRd$dN_postbooster, 2)*100, ')')
tb_dNRd$pr_dRd_postbooster=paste0(round(tb_dNRd$dRd_postbooster_mean/1000), ' (', round(tb_dNRd$dRd_postbooster_q25/1000), 
                                  ', ', round(tb_dNRd$dRd_postbooster_q75/1000), ')')
tb_dNRd$pr_dRd_postbooster_percentage=paste0(round(tb_dNRd$dRd_postbooster_mean/tb_dNRd$dRd_postbooster, 2)*100, ' (', 
                                             round(tb_dNRd$dRd_postbooster_q25/tb_dNRd$dRd_postbooster, 2)*100, 
                                             ', ', round(tb_dNRd$dRd_postbooster_q75/tb_dNRd$dRd_postbooster, 2)*100, ')')
tb_dNRd$dN_postbooster=round(tb_dNRd$dN_postbooster/1000)
tb_dNRd$dRd_postbooster=round(tb_dNRd$dRd_postbooster/1000)
tb_0=tb_dNRd[tb_dNRd$scenario=="No booster", c('province', 'dN_postbooster', 'pr_dN_postbooster', 'pr_dN_postbooster_percentage',
                                               'dRd_postbooster', 'pr_dRd_postbooster', 'pr_dRd_postbooster_percentage')]
tb_0_5=tb_dNRd[tb_dNRd$scenario=="0.5 times the boosters", c('province', 'pr_dN_postbooster', 'pr_dN_postbooster_percentage',
                                                             'pr_dRd_postbooster', 'pr_dRd_postbooster_percentage')]
colnames(tb_0_5)[-1]=paste0(c('pr_dN_postbooster', 'pr_dN_postbooster_percentage',
                              'pr_dRd_postbooster', 'pr_dRd_postbooster_percentage'), '_05')
tb_2=tb_dNRd[tb_dNRd$scenario=="2 times the boosters", c('province', 'pr_dN_postbooster', 'pr_dN_postbooster_percentage',
                                                         'pr_dRd_postbooster', 'pr_dRd_postbooster_percentage')]
colnames(tb_2)[-1]=paste0(c('pr_dN_postbooster', 'pr_dN_postbooster_percentage',
                            'pr_dRd_postbooster', 'pr_dRd_postbooster_percentage'), '_2')
tb_pr=inner_join(inner_join(tb_0, tb_0_5), tb_2)
tb_pr=apply(tb_pr, 2, as.character)
print(xtable(tb_pr[, c("province", "dN_postbooster", "pr_dN_postbooster", "pr_dN_postbooster_percentage", 
                       "pr_dN_postbooster_05", "pr_dN_postbooster_percentage_05",
                       "pr_dN_postbooster_2", "pr_dN_postbooster_percentage_2")]), include.rownames=F)
print(xtable(tb_pr[, c("province", "dRd_postbooster", "pr_dRd_postbooster", "pr_dRd_postbooster_percentage", 
                       "pr_dRd_postbooster_05", "pr_dRd_postbooster_percentage_05",
                       "pr_dRd_postbooster_2", "pr_dRd_postbooster_percentage_2")]), include.rownames=F)

dNdRd_pb=as.data.frame(increase_ND_boot[increase_ND_boot$scenario%in% c("No booster", "0.5 times the boosters","2 times the boosters"),
                                        c('province', 'dN_postbooster', 'dN_postbooster_mean', 'dN_postbooster_q75', 'dN_postbooster_q25',  
                                          'dRd_postbooster', 'dRd_postbooster_mean', 'dRd_postbooster_q75', 'dRd_postbooster_q25', 'scenario')]
                       %>% group_by(scenario) %>% summarise(totaldN_pb=sum(dN_postbooster), totaldN_pbm=sum(dN_postbooster_mean), 
                                                            totaldN_pbq25=sum(dN_postbooster_q25), totaldN_pbq75=sum(dN_postbooster_q75),
                                                            totaldRd_pb=sum(dRd_postbooster), totaldRd_pbm=sum(dRd_postbooster_mean),
                                                            totaldRd_pbq25=sum(dRd_postbooster_q25), totaldRd_pbq75=sum(dRd_postbooster_q75)))
round(dNdRd_pb$totaldN_pb[1]/1000)
paste0(round(dNdRd_pb$totaldN_pbm/1000), ' (', round(dNdRd_pb$totaldN_pbq25/1000), 
       ',', round(dNdRd_pb$totaldN_pbq75/1000), ')')
paste0(round(dNdRd_pb$totaldN_pbm/dNdRd_pb$totaldN_pb, 2)*100, ' (', round(dNdRd_pb$totaldN_pbq25/dNdRd_pb$totaldN_pb, 2)*100,
       ',', round(dNdRd_pb$totaldN_pbq75/dNdRd_pb$totaldN_pb, 2)*100, ')')

round(dNdRd_pb$totaldRd_pb[1]/1000)
paste0(round(dNdRd_pb$totaldRd_pbm/1000), ' (', round(dNdRd_pb$totaldRd_pbq25/1000), 
       ',', round(dNdRd_pb$totaldRd_pbq75/1000), ')')
paste0(round(dNdRd_pb$totaldRd_pbm/dNdRd_pb$totaldRd_pb, 2)*100, ' (', round(dNdRd_pb$totaldRd_pbq25/dNdRd_pb$totaldRd_pb, 2)*100,
       ',', round(dNdRd_pb$totaldRd_pbq75/dNdRd_pb$totaldRd_pb, 2)*100, ')')

#-----------
increase_ND_boot=countries_scenarios[countries_scenarios$t2==0,
                                     c('province', 'R.d', 'R.d.mean', 'R.d.q25', 'R.d.q75', 
                                       'N', 'N.mean', 'N.q25', 'N.q75', 
                                       'dN_postbooster_mean', 'dN_postbooster_q75', 'dN_postbooster_q25', 'dN_postbooster',
                                       'dRd_postbooster_mean', 'dRd_postbooster_q75', 'dRd_postbooster_q25', 'dRd_postbooster',
                                       'dN_postvaccine_mean', 'dN_postvaccine_q75', 'dN_postvaccine_q25', 'dN_postvaccine',
                                       'dRd_postvaccine_mean', 'dRd_postvaccine_q75', 'dRd_postvaccine_q25', 'dRd_postvaccine',
                                       'Rd.nova_Rd', 'N.nova_N', 'scenario')]
increase_ND_boot$Rd.increase=increase_ND_boot$R.d.mean-increase_ND_boot$R.d
increase_ND_boot$N.increase=increase_ND_boot$N.mean-increase_ND_boot$N
increase_ND_boot$dN_postvaccine_percentage=increase_ND_boot$dN_postvaccine_mean/increase_ND_boot$dN_postvaccine
increase_ND_boot$dRd_postvaccine_percentage=increase_ND_boot$dRd_postvaccine_mean/increase_ND_boot$dRd_postvaccine
write.csv(increase_ND_boot[increase_ND_boot$scenario%in% c("No vaccination", "Only the first dose"),], 
          '~/Documents/COVID-19/data/scenario_results_before_booster_reinfect.csv')

tb_dNRd=increase_ND_boot[increase_ND_boot$scenario%in% c("No vaccination", "Only the first dose"),
                         c('province', 'dN_postvaccine', 'dN_postvaccine_mean', 'dN_postvaccine_q75', 'dN_postvaccine_q25',  
                           'dRd_postvaccine', 'dRd_postvaccine_mean', 'dRd_postvaccine_q75', 'dRd_postvaccine_q25', 'scenario')]
tb_dNRd$pr_dN_postvaccine=paste0(round(tb_dNRd$dN_postvaccine_mean/1000), ' (', round(tb_dNRd$dN_postvaccine_q25/1000), 
                                 ', ', round(tb_dNRd$dN_postvaccine_q75/1000), ')')
tb_dNRd$pr_dN_postvaccine_percentage=paste0(round(tb_dNRd$dN_postvaccine_mean/tb_dNRd$dN_postvaccine, 2)*100, 
                                            ' (', round(tb_dNRd$dN_postvaccine_q25/tb_dNRd$dN_postvaccine, 2)*100, 
                                            ', ', round(tb_dNRd$dN_postvaccine_q75/tb_dNRd$dN_postvaccine, 2)*100, ')')
tb_dNRd$pr_dRd_postvaccine=paste0(round(tb_dNRd$dRd_postvaccine_mean/1000), ' (', round(tb_dNRd$dRd_postvaccine_q25/1000), 
                                  ', ', round(tb_dNRd$dRd_postvaccine_q75/1000), ')')
tb_dNRd$pr_dRd_postvaccine_percentage=paste0(round(tb_dNRd$dRd_postvaccine_mean/tb_dNRd$dRd_postvaccine, 2)*100, ' (', 
                                             round(tb_dNRd$dRd_postvaccine_q25/tb_dNRd$dRd_postvaccine, 2)*100, 
                                             ', ', round(tb_dNRd$dRd_postvaccine_q75/tb_dNRd$dRd_postvaccine, 2)*100, ')')
tb_dNRd$dN_postvaccine=round(tb_dNRd$dN_postvaccine/1000)
tb_dNRd$dRd_postvaccine=round(tb_dNRd$dRd_postvaccine/1000)
tb_no=tb_dNRd[tb_dNRd$scenario=="No vaccination", c('province', 'dN_postvaccine', 'pr_dN_postvaccine', 'pr_dN_postvaccine_percentage',
                                                    'dRd_postvaccine', 'pr_dRd_postvaccine', 'pr_dRd_postvaccine_percentage')]
tb_partial=tb_dNRd[tb_dNRd$scenario=="Only the first dose", c('province', 'pr_dN_postvaccine', 'pr_dN_postvaccine_percentage',
                                                              'pr_dRd_postvaccine', 'pr_dRd_postvaccine_percentage')]
colnames(tb_partial)[-1]=paste0(c('pr_dN_postvaccine', 'pr_dN_postvaccine_percentage',
                                  'pr_dRd_postvaccine', 'pr_dRd_postvaccine_percentage'), '_partial')
tb_pr=inner_join(tb_no, tb_partial)
tb_pr=apply(tb_pr, 2, as.character)
print(xtable(tb_pr[, c("province", "dN_postvaccine", "pr_dN_postvaccine", "pr_dN_postvaccine_percentage", 
                       "pr_dN_postvaccine_partial", "pr_dN_postvaccine_percentage_partial")]), include.rownames=F)
print(xtable(tb_pr[, c("province", "dRd_postvaccine", "pr_dRd_postvaccine", "pr_dRd_postvaccine_percentage", 
                       "pr_dRd_postvaccine_partial", "pr_dRd_postvaccine_percentage_partial")]), include.rownames=F)

dNdRd_pv=as.data.frame(increase_ND_boot[increase_ND_boot$scenario%in% c("No vaccination", "Only the first dose"),
                                        c('province', 'dN_postvaccine', 'dN_postvaccine_mean', 'dN_postvaccine_q75', 'dN_postvaccine_q25',  
                                          'dRd_postvaccine', 'dRd_postvaccine_mean', 'dRd_postvaccine_q75', 'dRd_postvaccine_q25', 'scenario')]
                       %>% group_by(scenario) %>% summarise(totaldN_pv=sum(dN_postvaccine), totaldN_pvm=sum(dN_postvaccine_mean), 
                                                            totaldN_pvq25=sum(dN_postvaccine_q25), totaldN_pvq75=sum(dN_postvaccine_q75),
                                                            totaldRd_pv=sum(dRd_postvaccine), totaldRd_pvm=sum(dRd_postvaccine_mean),
                                                            totaldRd_pvq25=sum(dRd_postvaccine_q25), totaldRd_pvq75=sum(dRd_postvaccine_q75)))
round(dNdRd_pv$totaldN_pv[1]/1000)
paste0(round(dNdRd_pv$totaldN_pvm/1000), ' (', round(dNdRd_pv$totaldN_pvq25/1000), 
       ',', round(dNdRd_pv$totaldN_pvq75/1000), ')')
paste0(round(dNdRd_pv$totaldN_pvm/dNdRd_pv$totaldN_pv, 2)*100, ' (', round(dNdRd_pv$totaldN_pvq25/dNdRd_pv$totaldN_pv, 2)*100,
       ',', round(dNdRd_pv$totaldN_pvq75/dNdRd_pv$totaldN_pv, 2)*100, ')')

round(dNdRd_pv$totaldRd_pv[1]/1000)
paste0(round(dNdRd_pv$totaldRd_pvm/1000), ' (', round(dNdRd_pv$totaldRd_pvq25/1000), 
       ',', round(dNdRd_pv$totaldRd_pvq75/1000), ')')
paste0(round(dNdRd_pv$totaldRd_pvm/dNdRd_pv$totaldRd_pv, 2)*100, ' (', round(dNdRd_pv$totaldRd_pvq25/dNdRd_pv$totaldRd_pv, 2)*100,
       ',', round(dNdRd_pv$totaldRd_pvq75/dNdRd_pv$totaldRd_pv, 2)*100, ')')

#--------------
# countries_nova   countries_partialva   countries_nobooster   countries_05boosterva   countries_2boosterva
data_I=countries_nova[countries_nova$date < as.Date("2022-03-16"), 
                      c('province', 'date', 't1', 't2', 'I', 'scenario')]
I_max=as.data.frame(data_I %>% group_by(province) %>% 
                      summarise(I_peak=max(I, na.rm=T)))
data_I_scenarios=rbind(countries_nova[, c('province', 'date', 't1', 't2', 'I', 'I.mean', 'I.q25', 'I.q75', 'scenario')],
                       countries_partialva[, c('province', 'date', 't1', 't2', 'I', 'I.mean', 'I.q25', 'I.q75', 'scenario')],
                       countries_nobooster[, c('province', 'date', 't1', 't2', 'I', 'I.mean', 'I.q25', 'I.q75', 'scenario')],
                       countries_05boosterva[, c('province', 'date', 't1', 't2', 'I', 'I.mean', 'I.q25', 'I.q75', 'scenario')],
                       countries_2boosterva[, c('province', 'date', 't1', 't2', 'I', 'I.mean', 'I.q25', 'I.q75', 'scenario')])
data_I_scenarios=data_I_scenarios[as.Date(data_I_scenarios$date) < as.Date("2022-03-16"), ]
data_I_scenarios=left_join(data_I_scenarios, I_max)
Im_day_scenarios=as.data.frame(data_I_scenarios %>% group_by(scenario, province) %>% 
                                 summarise(Im_peak=max(I.mean, na.rm=T), day_more=sum(I.mean>I_peak)))
Im_day_05=Im_day_scenarios[Im_day_scenarios$scenario=='0.5 times the boosters', c('province', 'Im_peak', 'day_more')]
colnames(Im_day_05)[-1]=paste0(colnames(Im_day_05)[-1], '_05')
Im_day_2=Im_day_scenarios[Im_day_scenarios$scenario=='2 times the boosters', c('province', 'Im_peak', 'day_more')]
colnames(Im_day_2)[-1]=paste0(colnames(Im_day_2)[-1], '_2')
Im_day_nb=Im_day_scenarios[Im_day_scenarios$scenario=='No booster', c('province', 'Im_peak', 'day_more')]
colnames(Im_day_nb)[-1]=paste0(colnames(Im_day_nb)[-1], '_nb')
Im_day_nv=Im_day_scenarios[Im_day_scenarios$scenario=='No vaccination', c('province', 'Im_peak', 'day_more')]
colnames(Im_day_nv)[-1]=paste0(colnames(Im_day_nv)[-1], '_nv')
Im_day_partial=Im_day_scenarios[Im_day_scenarios$scenario=='Only the first dose', c('province', 'Im_peak', 'day_more')]
colnames(Im_day_partial)[-1]=paste0(colnames(Im_day_partial)[-1], '_partial')

peaks=inner_join(inner_join(inner_join(inner_join(inner_join(I_max, Im_day_nv), Im_day_partial), 
                                       Im_day_nb), Im_day_05), Im_day_2)
tb_peak=peaks
tb_peak$I_peak_pr=round(tb_peak$I_peak/10^6, 2)
tb_peak$Im_peak_nv_pr=paste0(round(tb_peak$Im_peak_nv/10^6, 2), ' (',
                             round(tb_peak$Im_peak_nv/tb_peak$I_peak, 2)*100, ')')
mean(tb_peak$Im_peak_nv/tb_peak$I_peak)
sd(tb_peak$Im_peak_nv/tb_peak$I_peak)/sqrt(7)
tb_peak$Im_peak_partial_pr=paste0(round(tb_peak$Im_peak_partial/10^6, 2), ' (',
                                  round(tb_peak$Im_peak_partial/tb_peak$I_peak, 2)*100, ')')
mean(tb_peak$Im_peak_partial/tb_peak$I_peak)
sd(tb_peak$Im_peak_partial/tb_peak$I_peak)/sqrt(7)
tb_peak$Im_peak_nb_pr=paste0(round(tb_peak$Im_peak_nb/10^6, 2), ' (',
                             round(tb_peak$Im_peak_nb/tb_peak$I_peak, 2)*100, ')')
sum(tb_peak$Im_peak_nb-tb_peak$I_peak)
mean(tb_peak$Im_peak_nb/tb_peak$I_peak)
sd(tb_peak$Im_peak_nb/tb_peak$I_peak)/sqrt(7)
tb_peak$Im_peak_05_pr=paste0(round(tb_peak$Im_peak_05/10^6, 2), ' (',
                             round(tb_peak$Im_peak_05/tb_peak$I_peak, 2)*100, ')')
mean(tb_peak$Im_peak_05/tb_peak$I_peak)
sd(tb_peak$Im_peak_05/tb_peak$I_peak)/sqrt(7)
tb_peak$Im_peak_2_pr=paste0(round(tb_peak$Im_peak_2/10^6, 2), ' (',
                            round(tb_peak$Im_peak_2/tb_peak$I_peak, 2)*100, ')')
sum(tb_peak$Im_peak_2-tb_peak$I_peak)
mean(tb_peak$Im_peak_2/tb_peak$I_peak)
sd(tb_peak$Im_peak_2/tb_peak$I_peak)/sqrt(7)

print(xtable(apply(tb_peak[, c('province', 'I_peak_pr', 'Im_peak_nv_pr', 'day_more_nv',
                               'Im_peak_partial_pr', 'day_more_partial',
                               'Im_peak_nb_pr', 'day_more_nb',
                               'Im_peak_05_pr', 'day_more_05',
                               'Im_peak_2_pr', 'day_more_2')], 2, as.character)), include.rownames=F)

#--------------
country_data1=read.csv('~/Documents/COVID-19/data/country_data_theta1.csv', stringsAsFactors = F)
country_data1$type='Main analysis'
country_data2=read.csv('~/Documents/COVID-19/data/country_data_theta2.csv', stringsAsFactors = F)
country_data2$type='Asymptomatic rate' 
country_data3=read.csv('~/Documents/COVID-19/data/country_data_theta1_reinfect.csv', stringsAsFactors = F)
country_data3$type='Reinfection'
country_datass=rbind(country_data1, country_data2[, colnames(country_data1)], country_data3[, colnames(country_data1)])
country_datass$date=as.Date(country_datass$date)
data_plot=country_datass[country_datass$t>14 & country_datass$date < as.Date("2022-04-07")-14,]

mytheme_sen=theme_bw() + theme(axis.title = element_text(size = 15), 
                   text = element_text(face = "bold"),
                   strip.text = element_text(size = 12,face = 'bold'),
                   strip.background = element_rect(color="black", fill="white", linetype="solid"),
                   axis.title.x = element_text(size=15, face = 'bold', hjust = 0.5),
                   axis.title.y = element_text(size=15, face = 'bold', hjust = 0.5),
                   axis.text.x = element_text(size=10.2, angle=45,hjust = 1, face = 'bold'),
                   axis.text.y = element_text(size=13, face = 'bold'),
                   #axis.title.y.right=element_text(color='black'),
                   #plot.title = element_text(size=15, face = 'bold', hjust = 0.5),
                   plot.title = element_blank(),
                   legend.title = element_text(size=12, face = 'bold'),
                   legend.text = element_text(size=11, face = 'bold'),
                   legend.key.width  = unit(.3,"inches"),
                   legend.key.height = unit(.3,"inches"),
                   panel.grid.major=element_blank(), panel.grid.minor=element_blank())

ggplot(data_plot[data_plot$date<as.Date('2022-03-16')&data_plot$type%in%c('Asymptomatic rate', 'Reinfection'),])+
  scale_color_manual(values=c('#CC0000','#0000FF'),
                     labels=c('Asymptomatic rate', 'Reinfection'))+
  geom_line(aes(x=date, y=gamma.d.11, color=type))+
  geom_line(data=data_plot[data_plot$date<as.Date('2022-03-16')&data_plot$type=='Main analysis',], aes(x =date, y=gamma.d.11), color='black')+
  # geom_vline(data=time_point[, c('province', 'date')],
  #            aes(xintercept=date),color = "black",linetype = "dashed",size=0.7)+
  geom_rect(data=shades[shades$province%in%c("Brazil","Germany","Italy","Peru","Turkey","United Kingdom","United States"),],
            mapping=aes(xmin=x1, xmax=x2, ymin = -Inf, ymax = Inf, fill=Periods), alpha=0.4)+
  scale_fill_simpsons()+
  scale_x_date(date_breaks = "2 months",date_labels=format("%B %Y"))+
  facet_wrap(.~ province, nrow=2, scales = "free_y") + 
  mytheme_sen+theme(legend.position = 'bottom')+
  labs(x = "Date", y = expression(bold(gamma['d,t'])), color='Sensitivity analyses')
ggsave(paste0(path_plot,'gammad_omicrond_sensitivity.png'), units="in",width=13, height=7, dpi=300)

ggplot(data_plot[data_plot$date<as.Date('2022-03-16')&data_plot$type%in%c('Asymptomatic rate', 'Reinfection'),])+
  scale_color_manual(values=c('#CC0000','#0000FF'),
                     labels=c('Asymptomatic rate', 'Reinfection'))+
  geom_line(aes(x =date, y=gamma.r.11, color=type))+
  geom_line(data=data_plot[data_plot$date<as.Date('2022-03-16')&data_plot$type=='Main analysis',], aes(x =date, y=gamma.r.11), color='black')+
  # geom_vline(data=time_point[, c('province', 'date')],
  #            aes(xintercept=date),color = "black",linetype = "dashed",size=0.7)+
  geom_rect(data=shades[shades$province%in%c("Brazil","Germany","Italy","Peru","Turkey","United Kingdom","United States"),],
            mapping=aes(xmin=x1, xmax=x2, ymin = -Inf, ymax = Inf, fill=Periods), alpha=0.4)+
  scale_fill_simpsons()+
  scale_x_date(date_breaks = "2 months",date_labels=format("%B %Y"))+
  facet_wrap(.~ province, nrow=2, scales = "free_y") + 
  mytheme_sen+theme(legend.position = 'bottom')+
  labs(x = "Date", y = expression(bold(gamma['r,t'])), color='Sensitivity analyses')
ggsave(paste0(path_plot,'gammar_omicrond_sensitivity.png'), units="in",width=13, height=7, dpi=300)

ggplot(data_plot[data_plot$date<as.Date('2022-03-16')&data_plot$type%in%c('Asymptomatic rate', 'Reinfection'),])+
  scale_color_manual(values=c('#CC0000','#0000FF'),
                     labels=c('Asymptomatic rate', 'Reinfection'))+
  geom_line(aes(x =date, y=beta.mo.11, color=type))+
  geom_line(data=data_plot[data_plot$date<as.Date('2022-03-16')&data_plot$type=='Main analysis',], aes(x =date, y=beta.mo.11), color='black')+
  # geom_vline(data=time_point[, c('province', 'date')],
  #            aes(xintercept=date),color = "black",linetype = "dashed",size=0.7)+
  geom_rect(data=shades[shades$province%in%c("Brazil","Germany","Italy","Peru","Turkey","United Kingdom","United States"),],
            mapping=aes(xmin=x1, xmax=x2, ymin = -Inf, ymax = Inf, fill=Periods), alpha=0.4)+
  scale_fill_simpsons()+
  scale_x_date(date_breaks = "2 months",date_labels=format("%B %Y"))+
  facet_wrap(.~ province, nrow=2, scales = "fixed") + 
  labs(x = "Date", y = 'Infection rates', color='Sensitivity analyses')+
  #labs(x = "Date", y = expression(bold(beta[t])), color='Sensitivity analyses')+
  mytheme_sen+theme(legend.position = 'bottom')
ggsave(paste0(path_plot,'country_beta_sensitivity.png'), units="in",width=13, height=7, dpi=300)

ggplot(data_plot[data_plot$date<as.Date('2022-03-16')&data_plot$type%in%c('Asymptomatic rate', 'Reinfection'),])+
  scale_color_manual(values=c('#CC0000','#0000FF'),
                     labels=c('Asymptomatic rate', 'Reinfection'))+
  geom_line(aes(x =date, y=Rt.mo, color=type))+
  geom_line(data=data_plot[data_plot$date<as.Date('2022-03-16')&data_plot$type=='Main analysis',], aes(x =date, y=Rt.mo), color='black')+
  geom_hline(yintercept=1,color = "gray",linetype = "dashed",size=0.7)+
  # geom_vline(data=time_point[, c('province', 'date')],
  #            aes(xintercept=date),color = "black",linetype = "dashed",size=0.7)+
  geom_rect(data=shades[shades$province%in%c("Brazil","Germany","Italy","Peru","Turkey","United Kingdom","United States"),],
            mapping=aes(xmin=x1, xmax=x2, ymin = -Inf, ymax = Inf, fill=Periods), alpha=0.4)+
  scale_fill_simpsons()+
  scale_x_date(date_breaks = "2 months",date_labels=format("%B %Y"))+
  facet_wrap(.~ province, nrow=2, scales = "fixed") + 
  labs(x = "Date", y = 'Rt', color='Sensitivity analyses')+
  mytheme_sen+theme(legend.position = 'bottom')
ggsave(paste0(path_plot,'country_Rt_sensitivity.png'), units="in",width=13, height=7, dpi=300)

