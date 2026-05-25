# The covariance between the intercept and slope was fixed to be 0.
# normal distribution
pp <- 1
setwd('/scratch/gux8df/ML2/C1/pp1')
#setwd(path)
### generate data: normal data
library(magrittr) # %>%
library(mvtnorm)
library(dplyr) # if_else
set.seed(20250428)
# the number of measurement occasions: one level
T	<-	4
# missing rate: Four levels
miss.rate0 <- 0
miss.rate1 <- 0.05
miss.rate2 <- 0.10
miss.rate3	<-	0.15
miss.rate4	<-	0.3

# model related parameter
beta		<-	c(6,2)
lambda	<-	matrix(c(1,1,1,1,0,1,2,3),T,2)
mu		<-	lambda%*%beta

# parameters can be changed
for(N in c(100, 200, 500, 1000, 5000, 10000, 50000)){
  NT <- N*T
  diaE <- 1
  varE <-	diag(diaE,T)
  covD <- 0
  D		<-	array(c(1,covD,covD,1), dim=c(2,2))


  for(rep1 in 1:300){
    id <- rep(rep1,N)
    e	<-	matrix(rnorm(NT,0,sqrt(diaE)),N,T)
    u 	<-	rmvnorm(N,c(0,0),D)
    y	<-	matrix(rep(mu,N),N,T,byrow=TRUE)+t(lambda%*%t(u)) + e
    yout <- cbind(id,y)
    name <- paste('y.','pp',pp,'.N',N,'.txt',sep='')
    write.table(yout,file=name,append=TRUE,row.names=F,col.names=F)
    ##MCAR
    yMAR1	<-	y
    yMAR2	<-	y
    yMAR3	<-	y
    yMAR4	<-	y
    miss1 <- round(2*N*miss.rate1/(T-1))
    miss2 <- round(2*N*miss.rate2/(T-1))
    miss3 <- round(2*N*miss.rate3/(T-1))
    miss4 <- round(2*N*miss.rate4/(T-1))
    for (t in 1:(T-1)){
      yMAR1	<-	yMAR1[order(yMAR1[,t], decreasing=T),]
      yMAR1[(N-t*miss1+1):N,(t+1)]	<-	-99
      yMAR2	<-	yMAR2[order(yMAR2[,t], decreasing=T),]
      yMAR2[(N-t*miss2+1):N,(t+1)]	<-	-99
      yMAR3	<-	yMAR3[order(yMAR3[,t], decreasing=T),]
      yMAR3[(N-t*miss3+1):N,(t+1)]	<-	-99
      yMAR4	<-	yMAR4[order(yMAR4[,t], decreasing=T),]
      yMAR4[(N-t*miss4+1):N,(t+1)]	<-	-99
    }
    for(i in 1:N){
      for(j in 2:T){
        if(yMAR1[i,j]==-99){yMAR1[i,j]<-NA}
        if(yMAR2[i,j]==-99){yMAR2[i,j]<-NA}
        if(yMAR3[i,j]==-99){yMAR3[i,j]<-NA}
        if(yMAR4[i,j]==-99){yMAR4[i,j]<-NA}
      }
    }
    yMAR1out <- cbind(id,yMAR1)
    yMAR2out <- cbind(id,yMAR2)
    yMAR3out <- cbind(id,yMAR3)
    yMAR4out <- cbind(id,yMAR4)

    name1	<-	paste('MAR1.','pp',pp,'.N',N,'.txt',sep='')
    write.table(yMAR1out,file=name1,row.names=F,col.names=F,append=TRUE)
    name2	<-	paste('MAR2.','pp',pp,'.N',N,'.txt',sep='')
    write.table(yMAR2out,file=name2,row.names=F,col.names=F,append=TRUE)
    name3	<-	paste('MAR3.','pp',pp,'.N',N,'.txt',sep='')
    write.table(yMAR3out,file=name3,row.names=F,col.names=F,append=TRUE)
    name4	<-	paste('MAR4.','pp',pp,'.N',N,'.txt',sep='')
    write.table(yMAR4out,file=name4,row.names=F,col.names=F,append=TRUE)

    ##MNAR
    b1	<-	u[,2]+beta[2]
    corAb <- 0.8
    a <- corAb/sqrt(1-corAb^2)
    # var(Aux) =1
    Aux <- a*b1+rnorm(N,0,1)
    # Aux~N(a*beta[2],sqrt(a^2+1))
    yMNAR1	<-	y
    yMNAR2	<-	y
    yMNAR3	<-	y
    yMNAR4	<-	y
    miss1 <- 2*miss.rate1/(T-1)
    miss2 <- 2*miss.rate2/(T-1)
    miss3 <- 2*miss.rate3/(T-1)
    miss4 <- 2*miss.rate4/(T-1)
    for(j in 2:T){
      crit1 <- qnorm((1-(j-1)*miss1),a*beta[2],sqrt(a^2+1))
      yMNAR1[which(Aux>crit1),j] <- NA
      crit2 <- qnorm((1-(j-1)*miss2),a*beta[2],sqrt(a^2+1))
      yMNAR2[which(Aux>crit2),j] <- NA
      crit3 <- qnorm((1-(j-1)*miss3),a*beta[2],sqrt(a^2+1))
      yMNAR3[which(Aux>crit3),j] <- NA
      crit4 <- qnorm((1-(j-1)*miss4),a*beta[2],sqrt(a^2+1))
      yMNAR4[which(Aux>crit4),j] <- NA
    }
    yMNAR1out <- cbind(id,yMNAR1)
    yMNAR2out <- cbind(id,yMNAR2)
    yMNAR3out <- cbind(id,yMNAR3)
    yMNAR4out <- cbind(id,yMNAR4)

    name5	<-	paste('MNAR1.','pp',pp,'.N',N,'.txt',sep='')
    write.table(yMNAR1out,file=name5,row.names=F,col.names=F,append=TRUE)
    name6	<-	paste('MNAR2.','pp',pp,'.N',N,'.txt',sep='')
    write.table(yMNAR2out,file=name6,row.names=F,col.names=F,append=TRUE)
    name7	<-	paste('MNAR3.','pp',pp,'.N',N,'.txt',sep='')
    write.table(yMNAR3out,file=name7,row.names=F,col.names=F,append=TRUE)
    name8	<-	paste('MNAR4.','pp',pp,'.N',N,'.txt',sep='')
    write.table(yMNAR4out,file=name8,row.names=F,col.names=F,append=TRUE)
  }# rep1 1 to 500
  # diaE
}# N


