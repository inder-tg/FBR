
library(raster)
source("myFBRfunctions.R")

##----------- Preparación de datos ----------##

listFiles <- list.files(path = paste0(getwd(),"/data/"),
                        pattern = ".tif", full.names = TRUE)

zona <- raster(listFiles[1])
fraction <- raster(listFiles[2])
indice <- raster(listFiles[3])

# --- shp file

rootDir <- getwd()

dsn <- system.file("vectors", package = "rgdal")

dirCalidad <- paste0(getwd(),"/data/puntos_MN")

setwd(dirCalidad)

shpCalidad <- shapefile('MN_points_2008082')

calidad_raster <- rasterize(shpCalidad, indice, field="Qual_num")

setwd(rootDir)

# ---

fraction_zona <- fraction
fraction_zona[zona != 1] <- NA
indice_zona <- indice
indice_zona[zona != 1] <- NA

# fraction_bestQuality_nonzero <- fraction_zona
# fraction_bestQuality_nonzero[calidad_raster!=0] <- NA
# fraction_bestQuality_nonzero[fraction_zona==0] <- NA

fraction_goodQuality_nonzero <- fraction_zona
fraction_goodQuality_nonzero[calidad_raster!=1] <- NA
fraction_goodQuality_nonzero[fraction_zona==0] <- NA

# indice_bestQuality_nonzero <- indice_zona
# indice_bestQuality_nonzero[calidad_raster!=0] <- NA
# indice_bestQuality_nonzero[fraction_zona==0] <- NA

indice_goodQuality_nonzero <- indice_zona
indice_goodQuality_nonzero[calidad_raster!=1] <- NA
indice_goodQuality_nonzero[fraction_zona==0] <- NA

# ---

# points_fractionBest <- rasterToPoints(fraction_bestQuality_nonzero)
# points_indiceBest <- rasterToPoints(indice_bestQuality_nonzero)

points_fractionGood <- rasterToPoints(fraction_goodQuality_nonzero)
points_indiceGood <- rasterToPoints(indice_goodQuality_nonzero)

# str(points_fractionBest)

# par(mar=c(4.5,4.5,2,2))
# plot(points_indiceBest[,3], points_fractionBest[,3] * 1e-2,
#      pch=20, col = "gray",cex=1.25,
#      xlab = "MNDWI6", ylab = "Fracción de agua")

par(mar=c(4.5,4.5,2,2))
plot(points_indiceGood[,3], points_fractionGood[,3] * 1e-2,
     pch=20, col = "gray",cex=1.5,cex.lab=1.5,
     xlab = "MNDWI6", ylab = "Fracción de agua")

# ---

#-------------- Valores iniciales ---------#

dataX <- points_indiceGood[,3]
dataY <- points_fractionGood[,3] * 1e-2
lm_orig <- lm(dataY~dataX)
lm_fit <- lm_orig$coefficients[1] + lm_orig$coefficients[2]*dataX
summary(lm_orig)

par(mar=c(4.5,4.5,2,2))
# yRan <- range(dataY,lm_fit)
plot(dataX, dataY,#ylim=yRan,
     pch=20, col = "gray",cex.lab=1.5,cex=1.5,cex.axis=1.5,
     xlab = "MNDWI6", ylab = "Fracción de agua")
# par(new=TRUE)
points(dataX,lm_fit,col="orange",cex=1.25,pch=20)
legend("bottomright", pch=rep(20,2), col=c("gray","orange"),
       pt.cex=rep(1.75,2),
       legend=c("Datos", "Ajuste (RL)"), bty="n",cex=1.25)


par(mar=c(4.5,4.5,2,2))
plot(dataY, lm_orig$fitted.values,pch=20,
     col = "gray",cex=1.5,cex.lab=1.5,cex.axis=1.5,
     xlab = "Fracción de agua",
     ylab = "Ajustado (Regresión Lineal)")
abline(h=0,lty=2,col="blue",lwd=2)
abline(h=1,lty=2,col="blue",lwd=2)

# --- plots for conference paper

par(mar=c(0,4.5,2,0),adj=0,mfrow=c(2,1))
plot(0,0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', ylab = "",
     main = "A", cex.main=2)
par(mar=c(4.5,4.5,2,2),new=TRUE,adj=0.5)
plot(dataX, dataY,
     pch=20, col = "gray",cex.lab=1.5,cex=1.5,cex.axis=1.5,
     xlab = "MNDWI6", ylab = "Water fraction")
points(dataX,lm_fit,col="orange",cex=1.25,pch=20)
legend("bottomright", pch=rep(20,2), col=c("gray","orange"),
       pt.cex=rep(1.75,2),
       legend=c("Observations", "Linear regression fit"), bty="n",cex=1.25)

par(mar=c(0,4.5,2,0),adj=0)
plot(0,0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', ylab = "",
     main = "B",cex.main=2)
par(mar=c(4.5,4.5,2,2),new=TRUE,adj=0.5)
plot(dataY, lm_orig$fitted.values,pch=20,
     col = "gray",cex=1.5,cex.lab=1.5,cex.axis=1.5,
     xlab = "Water fraction",
     ylab = "Linear regression fit")
abline(h=0,lty=2,col="blue",lwd=2)
abline(h=1,lty=2,col="blue",lwd=2)


## --- Prueba con posterior_phi ---##

b0 <- as.numeric(lm_orig$coefficients[1])
b1 <- as.numeric(lm_orig$coefficients[2])

# --- this was run in CONABIO system

# chains <- multipleGibbsFBR(x=dataX,y=dataY,p_ini=0.75,omega_ini=0.75,
#                            phi_ini=3,beta_ini=c(b0,b1),numChains=20,
#                            samplesByChain=8000,numCores=23,
#                            stDevPhi=0.125,g=0.1,k=35,
#                            stDevBetaJump=c(0.15,0.15),stDevBetaPrior=c(5,5))

# ---
# par(mfrow=c(5,4))
# for(i in 1:20){
#   plot(chains[[i]]$p,ylab="p")
# }
#
# par(mfrow=c(5,4))
# for(i in 1:20){
#   acf(chains[[i]]$p,ylab="p")
# }
#
# thin_p <- thinning(lag=15,delta=2,sampleSize=500,chainSize = 4000)
#
# par(mfrow=c(5,4))
# for(i in 1:20){
#   plot(chains[[i]]$p[thin_p$sample],ylab="p")
# }
#
# par(mfrow=c(5,4))
# for(i in 1:20){
#   acf(chains[[i]]$p[thin_p$sample],ylab="p")
# }
# # ---
# par(mfrow=c(5,4))
# for(i in 1:20){
#   plot(chains[[i]]$phi,ylab="p")
# }
#
# par(mfrow=c(5,4))
# for(i in 1:20){
#   acf(chains[[i]]$phi,ylab="p")
# }
#
# thin_phi <- thinning(lag=15,delta=2,sampleSize=500,chainSize = 4000)
#
# par(mfrow=c(5,4))
# for(i in 1:20){
#   plot(chains[[i]]$phi[thin_phi$sample],ylab="phi")
# }
#
# par(mfrow=c(5,4))
# for(i in 1:20){
#   acf(chains[[i]]$phi[thin_p$sample],ylab="p")
# }
# # ---
# for(i in 1:20){
#   cat(chains[[i]]$acceptancePhi,"\n")
# }
# # ---
# par(mfrow=c(5,4))
# for(i in 1:20){
#   plot(chains[[i]]$beta0,ylab="beta0")
# }
#
# par(mfrow=c(5,4))
# for(i in 1:20){
#   acf(chains[[i]]$beta0,ylab="beta0")
# }
#
# thin_beta0 <- thinning(lag=15,delta=2,sampleSize=500,chainSize = 4000)
#
# par(mfrow=c(5,4))
# for(i in 1:20){
#   plot(chains[[i]]$beta0[thin_beta0$sample],ylab="beta0")
# }
#
# par(mfrow=c(5,4))
# for(i in 1:20){
#   acf(chains[[i]]$beta0[thin_beta0$sample],ylab="beta0")
# }
# #---
# par(mfrow=c(5,4))
# for(i in 1:20){
#   plot(chains[[i]]$beta1,ylab="beta1")
# }
#
# par(mfrow=c(5,4))
# for(i in 1:20){
#   acf(chains[[i]]$beta1,ylab="beta1")
# }
#
# thin_beta1 <- thinning(lag=15,delta=2,sampleSize=500,chainSize = 4000)
#
# par(mfrow=c(5,4))
# for(i in 1:20){
#   plot(chains[[i]]$beta1[thin_beta1$sample],ylab="beta1")
# }
#
# par(mfrow=c(5,4))
# for(i in 1:20){
#   acf(chains[[i]]$beta1[thin_beta1$sample],ylab="beta1")
# }
# # ---
# for(i in 1:20){
#   cat(chains[[i]]$acceptanceBeta,"\n")
# }
# # ---
#
# p <- list(20)
# phi <- list(20)
# beta0 <- list(20)
# beta1 <- list(20)
# for(k in 1:20){
#   p[[k]] <- chains[[k]]$p[thin_p$sample]
#   phi[[k]] <- chains[[k]]$phi[thin_phi$sample]
#   beta0[[k]] <- chains[[k]]$beta0[thin_beta0$sample]
#   beta1[[k]] <- chains[[k]]$beta1[thin_beta1$sample]
# }
#
# save(chains,file = paste0(getwd(),"/RData/outputRaw_FBRMArismas_conabio.RData"))
# save(p,phi,beta0,beta1,file = paste0(getwd(),"/RData/output_FBRMarismas_conabio.RData"))


# ---

load(paste0(getwd(),"/RData/output_FBRMarismas_conabio.RData"))

p_output <- unlist(p)
plot(p_output)
acf(p_output)

phi_output <- unlist(phi)
plot(phi_output)
acf(phi_output)

beta0_output <- unlist(beta0)
beta1_output <- unlist(beta1)

par(mfrow=c(2,3))
plot(beta0_output,type="l",lty=2,
     # pch=20,
     col="gray", xlab="",
     ylab=expression(beta[0]), cex.lab=1.75)
acf(beta0_output, main="", xlab="")
hist(beta0_output,main=expression(beta[0]),breaks=30,
     xlab="",freq=FALSE, cex.main=2)
plot(beta1_output,type="l", lty=2,
     # pch=20,
     col="gray", xlab="Muestras",
     ylab=expression(beta[1]), cex.lab=1.75)
acf(beta1_output, main="", xlab="Rezago")
hist(beta1_output,main=expression(beta[1]),breaks=30,
     xlab="",freq=FALSE,cex.main=2)

# --- plots for conference paper

par(mar=c(0,4.5,2,0),mfrow=c(2,2),adj=0)
plot(0,0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', ylab = "",
     main = "A", cex.main=2)
par(mar=c(4.5,4.5,2,2),new=TRUE,adj=0.5)
plot(p_output,
     pch=20, col = "gray",cex=1.5,cex.lab=1.5,cex.axis=1.5,
     ylab=expression(p),xlab="nSamples")
par(mar=c(0,4.5,2,0),adj=0)
plot(0,0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', ylab = "",
     main = "B", cex.main=2)
par(mar=c(4.5,4.5,2,2),new=TRUE,adj=0.5)
acf(p_output,main="",cex=1.5,cex.lab=1.5,cex.axis=1.5)

# ---

par(mar=c(0,4.5,2,0),adj=0)
plot(0,0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', ylab = "",
     main = "C", cex.main=2)
par(mar=c(4.5,4.5,2,2),new=TRUE,adj=0.5)
plot(phi_output,
     pch=20,
     col = "gray",cex=1.5,cex.lab=1.5,cex.axis=1.5,
     ylab=expression(phi),xlab="nSamples")

par(mar=c(0,4.5,2,0),adj=0)
plot(0,0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', ylab = "",
     main = "D", cex.main=2)
par(mar=c(4.5,4.5,2,2),new=TRUE,adj=0.5)
acf(phi_output,main="",cex=1.5,cex.lab=1.5,cex.axis=1.5)

# ---

par(mar=c(0,4.5,2,0),mfrow=c(2,2),adj=0)
plot(0,0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', ylab = "",
     main = "A", cex.main=2)
par(mar=c(4.5,4.5,2,2),new=TRUE,adj=0.5)
hist(beta0_output,main=expression(beta[0]),breaks=30,
     xlab="", freq=FALSE, cex.main=2)
par(mar=c(0,4.5,2,0),adj=0)
plot(0,0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', ylab = "",
     main = "B", cex.main=2)
par(mar=c(4.5,4.5,2,2),new=TRUE,adj=0.5)
acf(beta0_output, main="",cex=1.5,cex.lab=1.5,cex.axis=1.5)

# ---

par(mar=c(0,4.5,2,0),adj=0)
plot(0,0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', ylab = "",
     main = "C", cex.main=2)
par(mar=c(4.5,4.5,2,2),new=TRUE,adj=0.5)
hist(beta1_output,main=expression(beta[1]),breaks=30,
     xlab="",freq=FALSE,cex.main=2)
par(mar=c(0,4.5,2,0),adj=0)
plot(0,0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', ylab = "",
     main = "D", cex.main=2)
par(mar=c(4.5,4.5,2,2),new=TRUE,adj=0.5)
acf(beta1_output, main="",cex=1.5,cex.lab=1.5,cex.axis=1.5)


# ---

b0FBR <- median(beta0_output) #
b1FBR <- median(beta1_output) #

muFBR <- exp( b0FBR + dataX * b1FBR ) / (1 + exp( b0FBR + dataX * b1FBR ))

pFBR <- median(p_output)

Omega <- matrix(nrow = length(dataX), ncol = 2)
Omega[,1] <- muFBR/pFBR
Omega[,2] <- (1-muFBR)/(1-pFBR)
omegaTilde <- apply(Omega,1,min)

lambda1 <- muFBR + (1-pFBR) * omegaTilde
lambda2 <- muFBR - pFBR * omegaTilde

par(mar=c(4.5,4.5,2,2),mfrow=c(1,1))
plot(points_indiceGood[,3], points_fractionGood[,3] * 1e-2,
     pch=20, col = "gray", cex=1.5, cex.lab=1.5,cex.axis=1.5,
     xlab = "MNDWI6", ylab = "Fracción de agua")
# points(points_indiceGood[,3], lambda1, col = "green")
points(points_indiceGood[,3], lambda2, col = "orange", pch=20,
       cex=1.5)
legend("bottomright", pch=rep(20,2), col=c("gray","orange"),
       pt.cex=rep(1.75,2),
       legend=c("Datos", "Ajuste (RBF)"), bty="n", cex=1.25)

yRan <- range(dataY,lambda2)
yRan[1] <- yRan[1]-0.05
yRan[2] <- yRan[2]+0.05
par(mar=c(4.5,4.5,2,2))
plot(dataY, lambda2, ylim=yRan,
     pch=20, col = "gray",cex=1.5,cex.lab=1.5,cex.axis=1.5,
     xlab = "Fracción de agua", ylab = "Ajustado (Beta Flexible)")
abline(h=0,lty=2,col="blue",lwd=2)
abline(h=1,lty=2,col="blue",lwd=2)


par(mfrow=c(2,1))
# par(mar=c(4.5,4.5,2,2))
par(mar=c(2.5,4.5,2,2))
plot(points_indiceGood[,3], type="l", col="gray",
     cex=1.5,cex.axis=1.5,cex.lab=1.5,
     xlab = "", ylab = "MNDWI6")
par(mar=c(4.5,4.5,0,2))
plot(lambda2, type="l", col="gray",
     cex=1.5,cex.axis=1.5,cex.lab=1.5,
     xlab = "Observaciones (pixeles)", ylab = "RBF")

par(mfrow=c(2,1))
# par(mar=c(4.5,4.5,2,2))
par(mar=c(2.5,4.5,2,2))
plot(dataY, type="l", col="gray",
     cex=1.5,cex.axis=1.5,cex.lab=1.5,
     xlab = "", ylab = "Fraccion")
par(mar=c(4.5,4.5,0,2))
plot(lambda2, type="l", col="gray",
     cex=1.5,cex.axis=1.5,cex.lab=1.5,
     xlab = "Observaciones (pixeles)", ylab = "RBF")

# --- plots conference paper

# (b0_hdi <- ci(beta0_output, method = "HDI", ci=0.95))
# (b1_hdi <- ci(beta1_output, method = "HDI", ci=0.95))
#
# muFBR_hdi_low <- exp( b0_hdi$CI_low + dataX * b1_hdi$CI_low ) / (1 + exp( b0_hdi$CI_low + dataX * b1_hdi$CI_low ))
# muFBR_hdi_high <- exp( b0_hdi$CI_high + dataX * b1_hdi$CI_high ) / (1 + exp( b0_hdi$CI_high + dataX * b1_hdi$CI_high ))
#
# # pFBR_hdi <- median(p_output)
# (p_hdi <- ci(p_output, method = "HDI", ci=0.95))
#
# Omega_low <- matrix(nrow = length(dataX), ncol = 2)
# Omega_low[,1] <- muFBR_hdi_low/p_hdi$CI_low
# Omega_low[,2] <- (1-muFBR_hdi_low)/(1-p_hdi$CI_low)
# omegaTilde_low <- apply(Omega_low,1,min)
#
# lambda1_low <- muFBR_hdi_low + (1-p_hdi$CI_low) * omegaTilde_low
# lambda2_low <- muFBR_hdi_low - p_hdi$CI_low * omegaTilde_low
#
# Omega_high <- matrix(nrow = length(dataX), ncol = 2)
# Omega_high[,1] <- muFBR_hdi_high/p_hdi$CI_high
# Omega_high[,2] <- (1-muFBR_hdi_high)/(1-p_hdi$CI_high)
# omegaTilde_high <- apply(Omega_high,1,min)
#
# lambda1_high <- muFBR_hdi_high + (1-p_hdi$CI_high) * omegaTilde_high
# lambda2_high <- muFBR_hdi_high - p_hdi$CI_high * omegaTilde_high

b0FBR <- median(beta0_output) #
b1FBR <- median(beta1_output) #

muFBR <- exp( b0FBR + dataX * b1FBR ) / (1 + exp( b0FBR + dataX * b1FBR ))

pFBR <- median(p_output)

Omega <- matrix(nrow = length(dataX), ncol = 2)
Omega[,1] <- muFBR/pFBR
Omega[,2] <- (1-muFBR)/(1-pFBR)
omegaTilde <- apply(Omega,1,min)

lambda1 <- muFBR + (1-pFBR) * omegaTilde
lambda2 <- muFBR - pFBR * omegaTilde


# ---

par(mar=c(0,4.5,2,0),adj=0,mfrow=c(2,1))
plot(0,0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', ylab = "",
     main = "A", cex.main=2)
par(mar=c(4.5,4.5,2,2),new=TRUE,adj=0.5)
plot(points_indiceGood[,3], points_fractionGood[,3] * 1e-2,
     pch=20, col = "gray", cex=1.5, cex.lab=1.5,cex.axis=1.5,
     xlab = "MNDWI6", ylab = "Water fraction")
points(points_indiceGood[,3], lambda2, col = "orange", pch=20,
       cex=1.5)
legend("bottomright", pch=rep(20,2), col=c("gray","orange"),
       pt.cex=rep(1.75,2),
       legend=c("Observations", "FBR fit"), bty="n", cex=1.25)

yRan <- range(dataY,lambda2)
yRan[1] <- yRan[1]-0.025
yRan[2] <- yRan[2]+0.025

par(mar=c(0,4.5,2,0),adj=0)
plot(0,0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', ylab = "",
     main = "B",cex.main=2)
par(mar=c(4.5,4.5,2,2),new=TRUE,adj=0.5)
plot(dataY, lambda2, ylim=yRan,
     pch=20, col = "gray",cex=1.5,cex.lab=1.5,cex.axis=1.5,
     xlab = "Water fraction", ylab = "FBR fitted")
abline(h=0,lty=2,col="blue",lwd=2)
abline(h=1,lty=2,col="blue",lwd=2)

# ---

par(mfrow=c(2,1))
par(mar=c(2.5,4.5,2,2))
plot(points_indiceGood[,3], type="l", col="gray",
     cex=1.5,cex.axis=1.5,cex.lab=1.5,
     xlab = "", ylab = "MNDWI6")
par(mar=c(4.5,4.5,0,2))
plot(lambda2, type="l", col="gray",
     cex=1.5,cex.axis=1.5,cex.lab=1.5,
     xlab = "Observaciones (pixeles)", ylab = "RBF")

par(mfrow=c(2,1))
# par(mar=c(4.5,4.5,2,2))
par(mar=c(2.5,4.5,2,2))
plot(dataY, type="l", col="gray",
     cex=1.5,cex.axis=1.5,cex.lab=1.5,
     xlab = "", ylab = "Fraccion")
par(mar=c(4.5,4.5,0,2))
plot(lambda2, type="l", col="gray",
     cex=1.5,cex.axis=1.5,cex.lab=1.5,
     xlab = "Observaciones (pixeles)", ylab = "RBF")




