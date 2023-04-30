
##########################################
# Bayesian Model Averaging 结果代码
##########################################

bma=function(data,y,xvars){
  library(BMS)
  ndata=data[,c(y,xvars)]
  att = bms(ndata, mprior = "uniform", g = "UIP", user.int = F, mcmc="enumeration")
  coef=coef(att)
  pip=0.5
  bmaselect=names(coef[,1][coef[,1]>pip])
  sigx=paste(bmaselect, sep = " ",collapse = "+")
  formula=paste(y,sigx,sep = "~")
  model_bma <- lm(formula, data = data)
  model_bma
}

y='PM2.5'
xvars=c("ciswur", "ctrstp", "dgaca",  "gcrca",  "ieno",   "iesdio",
        "pss",    "psi",    "gpc",    "pd")
model=bma(data,y,xvars)
model1=bma(dat1,y,xvars)
model2=bma(dat2,y,xvars)
library(stargazer)
stargazer(model,model1,model2, column.sep.width="-8pt",
          align=TRUE,label = "bma",report = "vc*s",no.space=TRUE, omit.stat=c("LL","ser","f",'rsq','adj.rsq'))



##########################################
# Quantile regression结果代码
##########################################
taus=seq(0.1, 0.9, by = 0.1)

qr=function(data,y,xvars,taus){
  model=bma(data,y,xvars)
  coes=paste(xvars, sep = " ",collapse = "+")
  formula=paste(y,coes,sep = "~")
  library(quantreg)
  multi_rqfit <- rq(formula, data = data, tau = taus)
  multi_rqfit
}

qrr=function(data,y,xvars,taus){
  model=bma(data,y,xvars)
  coes_=xvars
  coes=paste(coes_, sep = " ",collapse = "+")
  formula=paste(y,coes,sep = "~")
  library(quantreg)
  ss=rq(formula, data = data, tau = taus)
  #https://www.rdocumentation.org/packages/quantreg/versions/5.85/topics/summary.rq
  ss=summary(ss,se='iid')
  
  l=list()
  for (i in 1:length(taus)){
    coefs=round(ss[[i]][[3]][,1],2)
    coefs_p=ss[[i]][[3]][,4]
    library(weights)
    sigs=starmaker(coefs_p, p.levels=c(.01, .05, .1), symbols=c("***", "**", "*"))
    r=paste0(coefs,sigs)
    l[[i]]=r
  }
  mat=simplify2array(l)
  colnames(mat)=taus
  #rownames(mat)=rownames(ss[[1]][[3]])[-1]
  rownames(mat)=c('Intercept',coes_)
  mat
}

library(xtable)
rdata=qrr(data,y,xvars,taus)
xtable(rdata)


# 画图
qrconf=function(data,y,xvars,taus){
  model=bma(data,y,xvars)
  coes_=xvars
  coes=paste(coes_, sep = " ",collapse = "+")
  formula=paste(y,coes,sep = "~")
  library(quantreg)
  ss=rq(formula, data = data, tau = taus)
  ss=summary(ss,se='nid')
  l=list()
  for (i in 1:length(taus)){
    coefs=round(ss[[i]][[3]][,1:2],3)
    l[[i]]=as.matrix(coefs)
  }
  l
}


plotf=function(j,data,y,xvars,taus){
  confl=qrconf(data,y,xvars,taus)
  
  cl=list()
  for (i in 1:length(taus)) {
    comat=confl[[i]]
    mean=comat[,1]
    ups=comat[,1]+1.96*comat[,2]
    lows=comat[,1]-1.96*comat[,2]
    cmat=cbind(mean,ups,lows)
    cl[[i]]=cmat
  }
  
  model=bma(data,y,xvars)
  s=summary(model)
  ss=s[[4]][,1:2]
  coef=model$coefficients
  coes_=names(model$coefficients)
  ll=list()
  for (i in 1:length(taus)) {
    ll[[i]]=cl[[i]][j,]
  }
  llmat=simplify2array(ll)
  
  plot(taus,llmat[1,],type = 'l',ylim = c(min(llmat),max(llmat)),xlab = '',ylab = coes_[j],cex.lab=1.8, cex.axis=1.8)
  polygon(c(rev(taus), taus), c(rev(llmat[2,]), llmat[3,]),col = 'grey',lty = 3)
  lines(taus,llmat[1,],col='black')
  abline(h=ss[j,1],col=rgb(0, 0, 0, 0.25))
  abline(h=ss[j,1]+1.96*ss[j,2],col=rgb(0, 0, 0, 0.25),lty=2)
  abline(h=ss[j,1]-1.96*ss[j,2],col=rgb(0, 0, 0, 0.25),lty=2)
}

pf=function(data,y,xvars,taus){
  model=bma(data,y,xvars)
  conames=names(model$coefficients)
  len=length(conames)
  for (i in 1:len) {
    nm=paste(as.character(i),'.jpeg',sep = '')
    #jpeg(nm,pointsize=6, width=2500, height=2000, res=600)
    plotf(i,data,y,xvars,taus)
    #dev.off()
  }
}

pf(data,y,xvars,taus)

##########################################
# LASSO 结果代码
##########################################

lasso=function(data,y,xvars){
  xvars=data[,xvars]
  yvar=data[,y]
  library(glmnet)
  set.seed(888)
  lambda_seq <- 10^seq(1, -1, by = -.01)
  cv_output=cv.glmnet(data.matrix(xvars), yvar, family = "gaussian", alpha = 1, lambda = lambda_seq)
  best_lam <- cv_output$lambda.min
  lasso_best <- glmnet(data.matrix(xvars), yvar, alpha = 1, lambda = best_lam)
  lasso_best$beta
  beta=as.numeric(lasso_best$beta)
  names(beta)=colnames(xvars)
  sigxx=beta[which(beta!=0)]
  sigx=paste(names(sigxx), sep = " ",collapse = "+")
  formula=paste(y,sigx,sep = "~")
  model <- lm(formula, data = data)
  model
}

model=lasso(data,y,xvars)
model1=lasso(dat1,y,xvars)
model2=lasso(dat2,y,xvars)

library(stargazer)
stargazer(model,model1,model2, column.sep.width="-8pt",
          align=TRUE,label = "lasso",report = "vc*s",no.space=TRUE, omit.stat=c("LL","ser","f",'rsq','adj.rsq'))

##########################################
# 变量重要性分析 结果代码
##########################################

library(caret)       # for general model fitting
library(rpart)       # for fitting decision trees
library(ipred)       # for fitting bagged decision trees
ames_bag0 <- train(
  PM2.5~.,
  data = data,
  method = "treebag",
  trControl = trainControl(method = "cv", number = 10),
  nbagg = 100,  
  control = rpart.control(minsplit = 2, cp = 0)
)

library(vip)
#jpeg('ames_bag01.jpeg',pointsize=6, width=2500, height=2000, res=600)
vip::vip(ames_bag0)
#vip::vip(ames_bag0,aesthetics = list(color='blue',fill='blue'))
#dev.off()

# 使用该代码，请引用文献：
# Cui X, Huang W, Deng W, Jia C. Spatial Patterns, Drivers and Heterogeneous Effects of PM2.5: Experience from China. Polish Journal of Environmental Studies. 2022;31(6):5633-5647. doi:10.15244/pjoes/152165.


