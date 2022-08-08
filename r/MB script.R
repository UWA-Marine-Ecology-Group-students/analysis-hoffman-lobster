setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dplyr)
library(magrittr)
dat <- read.csv('Data/BROAD_ymaze_results.csv',stringsAsFactors = F)

func1 <- function(x){
  x <- x[order(x)]
  x <- paste(x, collapse = "")
  return(x)}

dat %<>% group_by(rn=rownames(dat)) %>% mutate(habs=func1(c(A,B)))

tmp <- data.frame(habs=unique(dat$habs))
tmp$choice <- c('Bnk', 'Snd', 'Sg', 'Snd', 'Bnk', 'Bnk')
tmp$res <- 0

dat$res <- tmp$res[match(paste(dat$habs,dat$choice),paste(tmp$habs,tmp$choice))]
dat$res[is.na(dat$res)] <- 1
dat$tankside <- ifelse(dat$A==dat$choice, paste('A',dat$tank,sep=''),paste('B',dat$tank,sep=''))

head(dat)
habitats <- c('Blank','SeaGrass','Sand','Turf')

pin <- c('Bnk'=-3,'Sg'=-0.7,'Snd'=-3,'Trf'=-1.5, 
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

## Observations vs Estimates
dout <- mod(out$par,flag='print')
tmp1 <- tapply(dout$trial, list(dout$choice, dout$habs),length)
tmp1 <- tmp1/colSums(tmp1,na.rm=T) 

tmp <- tapply(dout$res, list(dout$choice, dout$habs),mean)
tmp <- tmp/colSums(tmp,na.rm=T) 
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

##Bootstrap
bstrap <- 100
Prefout <- array(NA, dim=c(4,4,bstrap),dimnames = list(a=names(pin)[1:4],b=names(pin)[1:4],ob=1:bstrap))
Tankout <- matrix(NA, nrow=bstrap, ncol=length(pin)-4,dimnames = list(ob=1:bstrap, par=names(pin)[5:length(pin)]))

for(i in 1:bstrap){
  dat <- odat[sample(1:nrow(odat),nrow(odat),replace = T),]
  out <- nlminb(pin, mod)
  Prefout[,,i] <- matrix(exp(out$par[1:4]),ncol=4,nrow=4)/(matrix(exp(out$par[1:4]),ncol=4,nrow=4)+matrix(exp(out$par[1:4]),ncol=4,nrow=4,byrow = T))-0.5
  Tankout[i,] <- out$par[5:length(pin)]
}
#warnings()

##Bootstrapped preference
## Habitat Preference
Prefout2 <- apply(Prefout, c(1,2), quantile,probs=0.025)
Prefout50 <- apply(Prefout, c(1,2), quantile,probs=0.5)
Prefout97 <- apply(Prefout, c(1,2), quantile,probs=0.975)

par(mfrow=c(2,1),las=1, mar=c(5,5,2,2))
bp <- barplot(Prefout50,beside = T, col=Col, ylim=c(-0.5,0.5),names.arg = paste(habitats,'vs ..'), ylab='Relative preference', xlab='Habitat Options')
arrows(bp, Prefout2,y1=Prefout97,code=3,angle=90,length=0.05)
par(xpd=T)
legend('top',ncol=4,fill=Col,legend=habitats,inset = -0.1)
par(xpd=F)

## Bootstrapped tank effect
Tankout2 <- apply(Tankout-rowMeans(Tankout),2,quantile,probs=c(0.025,0.5,0.975))
mx <- ceiling(max(abs(Tankout))*10)/10
bp <- barplot(Tankout2['50%',],beside = T,ylab='Tank effect from A1',ylim=c(-mx,mx), xlab='Tank side and number')
arrows(bp, Tankout2['2.5%',],y1=Tankout2['97.5%',],code=3,angle=90,length=0.05)
lines(bp,rep(0,length(bp)),lty=1)

# Publication Plots:

# Full Plot:

# Set up the plot grid:
m = c(1.1,1.1,0.5,1.5)
par(mfrow=c(1,1), mar=c(5,5,2,2),mai=m,las=1)
bp <- barplot(Prefout50,beside = T, col=Col, ylim=c(-0.5,0.5),names.arg = paste(habitats,'vs ..'), ylab='Relative preference', xlab='Habitat Options')
arrows(bp, Prefout2,y1=Prefout97,code=3,angle=90,length=0.05)
par(xpd=T)
legend('right',ncol=1,fill=Col,legend=habitats,inset = -0.22)
par(xpd=F)

str(Prefout50)
as.data.frame(Prefout50)


#Re-work the data for publication plots:

# Change the coloumn and thus the plotting order:
col.order <- c("Trf","Snd","Bnk","Sg")
(Pref_data <- Prefout50[col.order,col.order])
(Pref_upper <- Prefout97[col.order,col.order])
(Pref_lower <- Prefout2[col.order,col.order])

# Remove the duplicate data values that dont need plotting: 
Pref_data[c(1),c(1:4)] <- NA
Pref_data[c(2),c(2:4)] <- NA
Pref_data[c(3),c(3:4)] <- NA
Pref_data[c(4),c(4)] <- NA

Pref_upper[c(1),c(1:4)] <- NA
Pref_upper[c(2),c(2:4)] <- NA
Pref_upper[c(3),c(3:4)] <- NA
Pref_upper[c(4),c(4)] <- NA

Pref_lower[c(1),c(1:4)] <- NA
Pref_lower[c(2),c(2:4)] <- NA
Pref_lower[c(3),c(3:4)] <- NA
Pref_lower[c(4),c(4)] <- NA


#Save comaprisons data to csv for plotting:

Mean <- Pref_data[lower.tri(Pref_data, diag = F)]
Comparison <- paste(row.names(Pref_data)[col(Pref_data)],
                    colnames(Pref_data)[row(Pref_data)], sep= '-')[lower.tri(Pref_data, diag = F)]
Broad_Means<- data.frame(Comparison, Mean)

Upper <- Pref_upper[lower.tri(Pref_upper, diag = F)]
Comparison <- paste(row.names(Pref_upper)[col(Pref_upper)],
                    colnames(Pref_upper)[row(Pref_upper)], sep= '-')[lower.tri(Pref_upper, diag = F)]
Broad_Upper <- data.frame(Comparison, Upper)

Lower <- Pref_lower[lower.tri(Pref_lower, diag = F)]
Comparison <- paste(row.names(Pref_lower)[col(Pref_lower)],
                    colnames(Pref_lower)[row(Pref_lower)], sep= '-')[lower.tri(Pref_lower, diag = F)]
Broad_Lower <- data.frame(Comparison, Lower)

Broad_Data <- Broad_Means %>%
  full_join(Broad_Upper, by = "Comparison")%>%
  full_join(Broad_Lower, by = "Comparison")%>%
  glimpse()


# write.csv(Broad_Data, paste("Data/PLOTTING_Broad_Hab_Analysis_",Sys.Date(),".csv",sep = ""))

# Create dummy data to plot the axis colours:
x_coord <- c(1,5,6,10,11,15,16,20)
y_coord <- c(0,0,0,0,0,0,0,0)
colours <- c('chartreuse',grey(0.95),'wheat2','green4')


# Re-order the habitats and the corresponding colours:
habitats <- c('Turf','Blank','Sand','Seagrass')
Col <- c('chartreuse',grey(0.95),'wheat2','green4')

# Plotting to check it all looks good:
m = c(1.1,1.1,0.5,1.5)
par(mfrow=c(1,1), mar=c(5,5,2,2),mai=m,las=1)
bp <- barplot(Pref_data[c(1:4),c(1:3)],beside = T, col=Col, ylim=c(-0.5,0.5),names.arg = paste(habitats[1:3],'vs ..'), ylab='Frequency chosen', xlab='Habitat Options')
lines(x_coord[1:2],y_coord[1:2],col= colours[1], type="l", lwd=2)
lines(x_coord[3:4],y_coord[3:4],col= colours[2], type="l", lwd=2)
lines(x_coord[5:6],y_coord[5:6],col= colours[3], type="l", lwd=2)
# lines(x_coord[7:8],y_coord[7:8],col= colours[4], type="l", lwd=2)
arrows(bp, Pref_lower[c(1:4),c(1:3)],y1=Pref_upper[c(1:4),c(1:3)],code=3,angle=90,length=0.05)
par(xpd=T)
# legend('right',ncol=1,fill=Col,legend=habitats,inset = -0.22)
par(xpd=F)


# Publication Plotting:
png(filename = paste("Plots/Broad_Hab_Cats_",Sys.Date(), "_V2.png",sep = ""),
    width = 200, height = 120, units = "mm", bg = "White", res = 600, pointsize = 8)

m = c(1.1,1.1,0.5,1.5)
par(mfrow=c(1,1), mar=c(5,5,2,2),mai=m,las=1)
bp <- barplot(Pref_data[c(1:4),c(1:3)],beside = T, col=Col, ylim=c(-0.5,0.5),names.arg = paste(habitats[1:3],'vs ..'), ylab='Frequency chosen', xlab='Habitat Options')
lines(x_coord[1:2],y_coord[1:2],col= colours[1], type="l", lwd=2)
lines(x_coord[3:4],y_coord[3:4],col= colours[2], type="l", lwd=2)
lines(x_coord[5:6],y_coord[5:6],col= colours[3], type="l", lwd=2)
# lines(x_coord[7:8],y_coord[7:8],col= colours[4], type="l", lwd=2)
arrows(bp, Pref_lower[c(1:4),c(1:3)],y1=Pref_upper[c(1:4),c(1:3)],code=3,angle=90,length=0.05)
par(xpd=T)
legend('right',ncol=1,fill=Col,legend=habitats,inset = -0.25)
par(xpd=F)

dev.off()