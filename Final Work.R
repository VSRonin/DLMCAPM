setwd("E:\\Documenti Luca\\In Elaborazione\\Eclipse\\Assignment 2")
###############################################################################
# Uncomment to save output to file
# Testo<- file("FinalWork.log")
# sink(Testo, append=TRUE)
# sink(Testo, append=TRUE, type="message")
###############################################################################
# 
# Consider data on returns of two market assets and of a risk-free asset.
# First consider a univariate dynamic capital asset pricing model (CAPM).
# Provide estimates and forecasts by MLE and in a Bayesian approach.
# Now study the two market assets jointly, through a bivariate extension
# of the dynamic CAPM given by a seemingly unrelated regression
# model.
#
# Authors:
#  Luca Beldi - 1351030
#  Paolo Leonetti - 1302108
#
###############################################################################
rm(list=ls())

require(tseries)
require(dlm)
# We choose the biggest company (as for market capitalization) and the first in alphabetic order of the S&P500
ExxonPrice<-get.hist.quote("XOM",start="1990-1-1",end="2012-10-1",quote = "AdjClose")
ThreeMPrice<-get.hist.quote("MMM",start="1990-1-1",end="2012-10-1",quote = "AdjClose")
SePPrice<-get.hist.quote("^GSPC",start="1990-1-1",end="2012-10-1",quote = "AdjClose")
png(filename="Prices.png")
plot(ExxonPrice,ylab="Price",xlab="Time",main="Prices of Exxon Mobil and 3M in time")
lines(ThreeMPrice,col="gray",lty=2)
legend(7190,max(max(ExxonPrice),max(ThreeMPrice))-2,legend=c("Exxon Mobil","3M"),col=c("black","gray"),lty=c(1,2),cex=0.8)
# We use the 3 month TBILL as risk free rete. Data is from the Federal Reserve. Yields are in continuous Compunding
RiskFreeRate<-read.csv("TBILL.csv")
RiskFreeRate<-as.vector(RiskFreeRate$DailyYield)
RiskFreeRate<-RiskFreeRate[length(RiskFreeRate):2]
# We Build Daily returns
Exxon<-expm1(diff(log(ExxonPrice)))
ThreeM<-expm1(diff(log(ThreeMPrice)))
SeP<-expm1(diff(log(SePPrice)))

# Standard CAPM
x<-SeP-RiskFreeRate
Exxon<-Exxon-RiskFreeRate
ThreeM<-ThreeM-RiskFreeRate
ExxonCAPM<-lm(Exxon~x)
summary(ExxonCAPM)
ThreeMCAPM<-lm(ThreeM~x)
summary(ThreeMCAPM)

# CAPM as a dlm
# MLE
buildCAPM<-function(u){
	dlmModReg(x,dV=exp(u[1]),dW=exp(u[2:3]))
}
ExxonDlmCAPMp <- dlmMLE(Exxon, parm=rep(0,3),buildCAPM)
ExxonDlmCAPMp$convergence #The optimization converges
ExxonDlmCAPM<- buildCAPM(ExxonDlmCAPMp$par)
ExxonDlmCAPM["V"]
ExxonDlmCAPM["W"]
ThreeMDlmCAPMp <- dlmMLE(ThreeM, parm=rep(0,3),buildCAPM)
ThreeMDlmCAPMp$convergence #The optimization converges
ThreeMDlmCAPM<- buildCAPM(ThreeMDlmCAPMp$par)
ThreeMDlmCAPM["V"]
ThreeMDlmCAPM["W"]
ExxonSmooth <- dlmSmooth(Exxon,ExxonDlmCAPM)
ThreeMSmooth <- dlmSmooth(ThreeM,ThreeMDlmCAPM)

#Plot them
png(filename="Exxon CAPM Coefficients.png")
split.screen(c(2,1))
screen(1)
plot(as.Date(7305:(7304+length(dropFirst(ExxonSmooth$s)[,1]))),dropFirst(ExxonSmooth$s)[,1],type='l',ylab="Alpha",xlab="",main="Components of the CAPM for Exxon")
screen(2)
plot(as.Date(7305:(7304+length(dropFirst(ExxonSmooth$s)[,1]))),dropFirst(ExxonSmooth$s)[,2],type='l',ylab="Beta",xlab="Time")
png(filename="3M CAPM Coefficients.png")
split.screen(c(2,1))
screen(1)
plot(as.Date(7305:(7304+length(dropFirst(ThreeMSmooth$s)[,1]))),dropFirst(ThreeMSmooth$s)[,1],type='l',ylab="Alpha",xlab="",main="Components of the CAPM for 3M")
screen(2)
plot(as.Date(7305:(7304+length(dropFirst(ThreeMSmooth$s)[,1]))),dropFirst(ThreeMSmooth$s)[,2],type='l',ylab="Beta",xlab="Time")


#Bayesian Inference

Simulations<-25000
Burn<-5000
# Gibbs sampling takes several hours
################################################################
# Exxon.MC <- dlmGibbsDIG(Exxon,
#		dlmModReg(x, m0 = c(0, 1),C0 = diag(c(1000, 1000))),
#		a.y=1, b.y=1000, a.theta=10, b.theta = 1000,
#		thin=1, n.sample=Simulations, save.states=FALSE)
# save(Exxon.MC, file="Exxon.MC")
# ThreeM.MC <- dlmGibbsDIG(ThreeM,
#		dlmModReg(x, m0 = c(0, 1),C0 = diag(c(1000, 1000))),
#		a.y=1, b.y=1000, a.theta=10, b.theta = 1000,
#		thin=1, n.sample=Simulations, save.states=FALSE)
# save(ThreeM.MC, file="ThreeM.MC")
################################################################
#Load Pre Calculated Data
load("Exxon.MC")
load("ThreeM.MC")
################################################################

#Plot Simulations Diagnosis
png(filename="Exxon Simulation Diagnosis.png")
par(mfrow=c(4,3),mar=c(0,5,5,0), oma=c(3,3,0,3))
plot(1:(Simulations-Burn),Exxon.MC$dV[-(1:Burn)],type="p",main="V",xlab="", ylab="")
plot(1:(Simulations-Burn),(Exxon.MC$dW[,1])[-(1:Burn)],type="p",xlab="", ylab="",main=expression("W"[1]))
plot(1:(Simulations-Burn),(Exxon.MC$dW[,2])[-(1:Burn)],type="p",main=expression("W"[2]),xlab="", ylab="")
plot(ergMean(Exxon.MC$dV[-(1:Burn)],0.05*(Simulations-Burn)),type="l", xlab="", ylab="")
plot(ergMean((Exxon.MC$dW[,1])[-(1:Burn)],0.05*(Simulations-Burn)),type="l", xlab="", ylab="")
plot(ergMean((Exxon.MC$dW[,2])[-(1:Burn)],0.05*(Simulations-Burn)),type="l", xlab="", ylab="")
acf(Exxon.MC$dV[-(1:Burn)],ylab="",main="")
acf((Exxon.MC$dW[,1])[-(1:Burn)],ylab="",main="")
acf((Exxon.MC$dW[,2])[-(1:Burn)],ylab="",main="")
par(mar=c(5,5,3,0))
plot(Exxon.MC$dV[-(1:Burn)],(Exxon.MC$dW[,1])[-(1:Burn)],pch=19,cex=.5,xlab="V",ylab=expression("W"[1]))
plot(Exxon.MC$dV[-(1:Burn)],(Exxon.MC$dW[,2])[-(1:Burn)],pch=19,cex=.5,xlab="V",ylab=expression("W"[2]))
plot((Exxon.MC$dW[,1])[-(1:Burn)],(Exxon.MC$dW[,2])[-(1:Burn)],pch=19,cex=.5,ylab=expression("W"[2]),xlab=expression("W"[1]))

png(filename="3M Simulation Diagnosis.png")
par(mfrow=c(4,3),mar=c(0,5,5,0), oma=c(3,3,0,3))
plot(1:(Simulations-Burn),ThreeM.MC$dV[-(1:Burn)],type="p",main="V",xlab="", ylab="")
plot(1:(Simulations-Burn),(ThreeM.MC$dW[,1])[-(1:Burn)],type="p",xlab="", ylab="",main=expression("W"[1]))
plot(1:(Simulations-Burn),(ThreeM.MC$dW[,2])[-(1:Burn)],type="p",main=expression("W"[2]),xlab="", ylab="")
plot(ergMean(ThreeM.MC$dV[-(1:Burn)],0.05*(Simulations-Burn)),type="l", xlab="", ylab="")
plot(ergMean((ThreeM.MC$dW[,1])[-(1:Burn)],0.05*(Simulations-Burn)),type="l", xlab="", ylab="")
plot(ergMean((ThreeM.MC$dW[,2])[-(1:Burn)],0.05*(Simulations-Burn)),type="l", xlab="", ylab="")
acf(ThreeM.MC$dV[-(1:Burn)],ylab="",main="")
acf((ThreeM.MC$dW[,1])[-(1:Burn)],ylab="",main="")
acf((ThreeM.MC$dW[,2])[-(1:Burn)],ylab="",main="")
par(mar=c(5,5,3,0))
plot(ThreeM.MC$dV[-(1:Burn)],(ThreeM.MC$dW[,1])[-(1:Burn)],pch=19,cex=.5,xlab="V",ylab=expression("W"[1]))
plot(ThreeM.MC$dV[-(1:Burn)],(ThreeM.MC$dW[,2])[-(1:Burn)],pch=19,cex=.5,xlab="V",ylab=expression("W"[2]))
plot((ThreeM.MC$dW[,1])[-(1:Burn)],(ThreeM.MC$dW[,2])[-(1:Burn)],pch=19,cex=.5,ylab=expression("W"[2]),xlab=expression("W"[1]))

Exxon.Bay.DLM<- dlmModReg(x, m0 = c(0, 1), C0 = diag(c(1000, 1000)),
		dV=mean(Exxon.MC$dV[-(1:Burn)]), dW=c(mean((Exxon.MC$dW[,1])[-(1:Burn)]),mean((Exxon.MC$dW[,2])[-(1:Burn)])))
Exxon.Bay.DLM["V"]
Exxon.Bay.DLM["W"]
ThreeM.Bay.DLM<- dlmModReg(x, m0 = c(0, 1), C0 = diag(c(1000, 1000)),
		dV=mean(ThreeM.MC$dV[-(1:Burn)]), dW=c(mean((ThreeM.MC$dW[,1])[-(1:Burn)]),mean((ThreeM.MC$dW[,2])[-(1:Burn)])))
ThreeM.Bay.DLM["V"]
ThreeM.Bay.DLM["W"]
Exxon.Bay.CAPM<- dlmSmooth(Exxon,Exxon.Bay.DLM)
ThreeM.Bay.CAPM<- dlmSmooth(ThreeM,ThreeM.Bay.DLM)

png(filename="Exxon Bayesian CAPM Coefficients.png")
split.screen(c(2,1))
screen(1)
plot(as.Date(7305:(7304+length(dropFirst(Exxon.Bay.CAPM$s)[,1]))),dropFirst(Exxon.Bay.CAPM$s)[,1],type='l',ylab="Alpha",xlab="",main="Components of the CAPM for Exxon")
screen(2)
plot(as.Date(7305:(7304+length(dropFirst(Exxon.Bay.CAPM$s)[,1]))),dropFirst(Exxon.Bay.CAPM$s)[,2],type='l',ylab="Beta",xlab="Time")

png(filename="3M Bayesian CAPM Coefficients.png")
split.screen(c(2,1))
screen(1)
plot(as.Date(7305:(7304+length(dropFirst(ThreeM.Bay.CAPM$s)[,1]))),dropFirst(ThreeM.Bay.CAPM$s)[,1],type='l',ylab="Alpha",xlab="",main="Components of the CAPM for 3M")
screen(2)
plot(as.Date(7305:(7304+length(dropFirst(ThreeM.Bay.CAPM$s)[,1]))),dropFirst(ThreeM.Bay.CAPM$s)[,2],type='l',ylab="Beta",xlab="Time")

#Multivariate Extension
BuildSUR <- function(u){
	# u is a vector contaning, in order:
	#	Variance of the observation for Exxon
	#	Variance of the observation for 3M
	#	Variance of the state proces for the Exxon alpha
	#	Variance of the state proces for the 3M alpha
	#	Covariance between the alpha states
	#	Variance of the state proces for the Exxon beta
	#	Variance of the state proces for the 3M beta
	#	Covariance between the beta states
	SateV1 <- matrix(c(exp(u[3]),u[5],u[5],exp(u[4])),nrow=2)
	SateV2 <- matrix(c(exp(u[6]),u[8],u[8],exp(u[7])),nrow=2)
	Result <- dlmModReg(x)
	Result$FF <- Result$FF %x% diag(2)
	Result$GG <- Result$GG %x% diag(2)
	Result$JFF <- Result$JFF %x% diag(2)
	Result$V <- diag(c(exp(u[1]),exp(u[2])))
	Result$W <- bdiag(SateV1,SateV2)
	Result$m0 <- c(rep(0,2),rep(1,2))
	Result$C0 <- diag(1000,nr=4)
	dlm(Result)
}
Stocks<-matrix(c(as.vector(Exxon),as.vector(ThreeM)),ncol=2)
colnames(Stocks)<-c("Exxon","3M")
# MLE estimation may require about half an hour of computation
################################################################
# StocksCAPMp <- dlmMLE(Stocks,rep(1e-7,8),BuildSUR,
#	method="SANN",control=list(trace=6, maxit=10000))
# save(StocksCAPMp,file="StocksCAPMp")
################################################################
#Load Pre Calculated Data
load("StocksCAPMp")
################################################################
StocksCAPMp$convergence #The optimization converges
StocksCAPM<-BuildSUR(StocksCAPMp$par)
StocksCAPM["V"]
StocksCAPM["W"]
SmoothedSUR<-dlmSmooth(Stocks,StocksCAPM)

png(filename="SUR CAPM Coefficients.png")
split.screen(c(2,2))
screen(1)
plot(as.Date(7305:(7304+length(dropFirst(SmoothedSUR$s)[,1]))),dropFirst(SmoothedSUR$s)[,1],type='l',ylab="Alpha",xlab="",main="Components for Exxon")
screen(3)
plot(as.Date(7305:(7304+length(dropFirst(SmoothedSUR$s)[,1]))),dropFirst(SmoothedSUR$s)[,3],type='l',ylab="Beta",xlab="Time")
screen(2)
plot(as.Date(7305:(7304+length(dropFirst(SmoothedSUR$s)[,1]))),dropFirst(SmoothedSUR$s)[,2],type='l',ylab="Alpha",xlab="",main="Components for 3M")
screen(4)
plot(as.Date(7305:(7304+length(dropFirst(SmoothedSUR$s)[,1]))),dropFirst(SmoothedSUR$s)[,4],type='l',ylab="Beta",xlab="Time")


#Close all Painting Devices
while (dev.off()>1){}
# Restore Output to Console
# sink() 
# sink(type="message")


