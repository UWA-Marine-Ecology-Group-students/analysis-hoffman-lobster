

library(dplyr)
library(magrittr)

dat <- read.csv('Lobster refined for R.csv',stringsAsFactors = F)

func1 <- function(x){
  x <- x[order(x)]
  x <- paste(x, collapse = "")
  return(x)}

dat %<>% group_by(rn=rownames(dat)) %>% mutate(scent=func1(c(A,B)))

tmp <- data.frame(scent=unique(dat$scent))
 ## BLANK - still waiting M's input
tmp$res <- 0

dat$res <- tmp$res[match(paste(dat$scent,dat$choice),paste(tmp$scent,tmp$choice))] 
dat$res[is.na(dat$res)] <- 1
dat$tankside <- ifelse(dat$A==dat$choice, paste('A',dat$tank,sep=''),paste('B',dat$tank,sep=''))

head(dat)
habitats <- c('Blank','Adult','Juvenile','Mix')


pin <- c('Bnk'=-3,'Ad'=-0.7,'Ju'=-3,'Mx'=-1.5,  #h- where do the numbers come from?
         A2=0.02,A3=0.03,A4=0.04,A5=0.05,A6=0.06,
         B1=0.0001,B2=0.002,B3=0.003,B4=0.004,B5=0.005,B6=0.006) 

(x <- pin)

dat <- dat[!is.na(dat$choice),] 

mod <- function(x,flag='solve'){
  pars <- x
  pref <- exp(pars[1:4])
  #pref <- pref+abs(min(pref))
  tanks <- c(A1=0, pars[5:length(x)]) #<- pars[5:length(x)]/mean(pars[5:length(x)])
  tanks <- tanks-mean(tanks)# + abs(min(tanks))
  dat$tA <- pref[match(dat$A,names(pref))]
  dat$tB <- pref[match(dat$B,names(pref))]
  dat$tankEff <- tanks[match(as.character(dat$tankside), names(tanks))]
  dat$res <- ifelse(dat$A==dat$choice, dat$tA/(dat$tA+dat$tB), dat$tB/(dat$tA+dat$tB))+dat$tankEff
  if(flag=='print') return(dat)
  LL <- -sum(log(dat$res))
  if(flag=='solve') return(LL)
}

mod(pin)
out <- nlminb(pin, mod)
out
x <- out$par



#### Observations vs Estimates
dout <- mod(out$par,flag='print')
tmp1 <- tapply(dout$trial, list(dout$choice, dout$scent),length)
tmp1 <- tmp1/colSums(tmp1,na.rm=TRUE) 

tmp <- tapply(dout$res, list(dout$choice, dout$scent),mean)
tmp <- tmp/colSums(tmp,na.rm=TRUE) 
par(mfrow=c(2,1),las=1,xpd=T, mar=c(5,5,4,2))
Col <- c(grey(0.95),'green4','wheat2','chartreuse')
barplot(tmp1,beside = T,col=Col,ylim=c(0,1),ylab='Proportion',main='Observed')
legend('bottom',ncol=4,fill=Col,legend=(rownames(tmp)),inset = -0.45)
barplot(tmp,beside = T,col=Col,ylim=c(0,1),ylab='Proportion',main='Estimated')
par(xpd=F)


## Habitat Preference
pref <- exp(x[1:4])
pref2 <- matrix(pref, ncol=4, nrow=4, dimnames = list(a=names(pref),b=names(pref)))
pref3 <- matrix(pref, ncol=4, nrow=4, dimnames = list(a=names(pref),b=names(pref)),byrow = T)
pref4 <- pref2/(pref3 + pref2)-0.5 
pref4[pref4==0] <- NA

par(mfrow=c(2,1),las=1,xpd=T)
barplot(pref4,beside = T, col=Col)
legend('bottom',ncol=4,fill=Col,legend=names(pref),inset = -0.45)


## Tank effect
tanks <- x[5:length(x)]-mean(x[5:length(x)])
barplot(tanks,beside = T,ylab='Tank effect',ylim=c(-.2,0.2))
par(xpd=F)

odat <- dat  # keep original data

























