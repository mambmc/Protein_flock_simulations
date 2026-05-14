# Coarse-grained Brownian model of protein flocks with 2 species
library(tidyr)
library(ijtiff)
library(EBImage)

# Initialize simulation and image parameters
nPixel <- 32 # pixels in the X and Y dimensions in the TIF image
xUCL <- 1000 # length of unit cell in X dimension in nm
yUCL <- 1000 # size of unit cell in Y dimension in nm
zUCL <- 1000 # size of unit cell in Z dimension in nm
nFrm <- 20 # frames
nMol <- 10000 # mean number of molecules
ppm <- 2 # number of photons per molecule in 10 us 
maxMolPrt <- 20 # max number of molecules per particle
oLap <- 1 # species overlap in particles
cDiff <- 20 # diffusion coefficient of an average protein monomer in um2/s or nm2/us
movT <- 10 # dwell time or snapshot time interval in us 
cf <- 0.84 # Correction factor for Gaussian approximation to photon sampling  
nPrt <- nMol/maxMolPrt  
nMolPrt <- ceiling(runif(nPrt,0,2*maxMolPrt-1))
nMolPrt1 <- rep(NA,nPrt)
nMolPrt2 <- rep(NA,nPrt)
hP <- floor(nPrt/2)
nMolPrt1[1:hP] <- round(nMolPrt[1:hP]*oLap/2)
nMolPrt2[1:hP] <- nMolPrt[1:hP]-nMolPrt1[1:hP]
nMolPrt2[(hP+1):nPrt] <- round(nMolPrt[(hP+1):nPrt]*oLap/2)
nMolPrt1[(hP+1):nPrt] <- nMolPrt[(hP+1):nPrt]-nMolPrt2[(hP+1):nPrt]
nSim <- nPixel^2
pixelSize <- xUCL/nPixel

# Create arrays and vectors
ACXGlobal1 <- array(0,dim=c(15,nFrm))
ACXGlobal2 <- array(0,dim=c(15,nFrm))
Bdataset1 <- rep(0,nFrm)
Bdataset2 <- rep(0,nFrm)
cRICS <- rep(0,nFrm)
newX <- rep(0,nPrt)
newY <- rep(0,nPrt)
newZ <- rep(0,nPrt)
movD <- rep(0,nPrt)
lastX <- rep(0,nPrt)
lastY <- rep(0,nPrt)
lastZ <- rep(0,nPrt)
imgStack <- array(0,dim=c(nPixel,nPixel,2,nSim))
imgMatrix1 <- matrix(0,nPixel,nPixel)
imgMatrix2 <- matrix(0,nPixel,nPixel)
rImg <- array(0,dim=c(nPixel,nPixel,2,nFrm))
wImg <- array(0,dim=c(nPixel,nPixel,2,nFrm))

# Assign initial random positions and calculate average movD
lastX <- runif(nPrt,0,xUCL-1)
lastY <- runif(nPrt,0,yUCL-1)
lastZ <- runif(nPrt,0,zUCL-1)
movD <- (2*cDiff*movT/(nMolPrt^0.3))^0.5

# Photon sampling function
rSampling <- function(x) rnorm(1,mean=x,sd=cf*(x^0.5))

# Perform simulations every movT (us)
for (iFrm in (1:nFrm)){
  print(paste("SIM iFrm= ",iFrm))
  imgStack[] <- 0
  for (iSim in (1:nSim)){
    # Brownian perturbations
    newX <- lastX+movD*rnorm(nPrt,mean=0,sd=1)
    newY <- lastY+movD*rnorm(nPrt,mean=0,sd=1)
    newZ <- lastZ+movD*rnorm(nPrt,mean=0,sd=1)
    # Reflect if out of bounds
    newX[which(newX<0)] <- -newX[which(newX<0)]
    newY[which(newY<0)] <- -newY[which(newY<0)]
    newZ[which(newZ<0)] <- -newZ[which(newZ<0)]
    newX[which(newX>xUCL)] <- 2*xUCL-newX[which(newX>xUCL)]
    newY[which(newY>yUCL)] <- 2*yUCL-newY[which(newY>yUCL)]
    newZ[which(newZ>zUCL)] <- 2*zUCL-newZ[which(newZ>zUCL)]
    # Save new positions
    lastX <- newX
    lastY <- newY
    lastZ <- newZ
    # Add counts to pixels in imgStack by position  
    cX <- ceiling(lastX/pixelSize)
    cY <- ceiling(lastY/pixelSize)
    imgMatrix1[] <- 0
    imgMatrix2[] <- 0
    for (iPos in (1:nPrt)){
      imgMatrix1[cY[iPos],cX[iPos]] <- imgMatrix1[cY[iPos],cX[iPos]] + ppm*nMolPrt1[iPos]
      imgMatrix2[cY[iPos],cX[iPos]] <- imgMatrix2[cY[iPos],cX[iPos]] + ppm*nMolPrt2[iPos]
    }
    # Apply a Gaussian filter to mimic PSF
    imgMatrix1 <- gblur(imgMatrix1,2)
    imgMatrix2 <- gblur(imgMatrix2,2)
    # Apply photon sampling variability
    imgMatrix1 <- matrix(vapply(imgMatrix1,rSampling,numeric(1)),nPixel,nPixel)
    imgMatrix2 <- matrix(vapply(imgMatrix2,rSampling,numeric(1)),nPixel,nPixel)
    imgMatrix1 <- round(imgMatrix1)
    imgMatrix2 <- round(imgMatrix2)
    imgMatrix1[which(imgMatrix1<0)] <- 0
    imgMatrix2[which(imgMatrix2<0)] <- 0
    # Add snapshot to stack
    imgStack[,,1,iSim] <- imgMatrix1[,]
    imgStack[,,2,iSim] <- imgMatrix2[,]
  }
  # Create TIFF image simulating raster scanning in confocal microscopy
  for (x in (1:nPixel)){
    for (y in (1:nPixel)){
      rImg[x,y,1,iFrm] <- imgStack[x,y,1,nPixel*(y-1)+x]
      rImg[x,y,2,iFrm] <- imgStack[x,y,2,nPixel*(y-1)+x]
      wImg[y,x,1,iFrm] <- imgStack[x,y,1,nPixel*(y-1)+x]
      wImg[y,x,2,iFrm] <- imgStack[x,y,2,nPixel*(y-1)+x]
    }
  }
}
# write_tif writes [y,x] matrix or stack as an ImageJ object
write_tif(wImg,"raster.tif",bits_per_sample=32,overwrite=TRUE,msg=FALSE)

# Perform RICS analysis
# Obtain mean and sd from pixel values per frame
meanImgPV1 <- rep(0,nFrm)
meanImgPV2 <- rep(0,nFrm)
sdImgPV1 <- rep(0,nFrm)
sdImgPV2 <- rep(0,nFrm)
for (iFrm in (1:nFrm)){
  meanImgPV1[iFrm] <- mean(rImg[,,1,iFrm])
  meanImgPV2[iFrm] <- mean(rImg[,,2,iFrm])
  sdImgPV1[iFrm] <- sd(rImg[,,1,iFrm])
  sdImgPV2[iFrm] <- sd(rImg[,,2,iFrm])
}
# Detrending procedures
detInt <- 2
dImg <- array(0,dim=c(nPixel,nPixel,2,nFrm))
for (x in (1:nPixel)){
  for (y in (1:nPixel)){
    for (iFrm in (1:nFrm)){
      pos1 <- iFrm-detInt
      if (pos1<1) pos1 <- 1
      pos2 <- iFrm+detInt
      if (pos2>nFrm) pos2 <- nFrm
      meanPV1 <- mean(rImg[x,y,1,(pos1:pos2)])
      meanPV2 <- mean(rImg[x,y,2,(pos1:pos2)])
      dImg[x,y,1,iFrm] <- rImg[x,y,1,iFrm]-meanPV1+meanImgPV1[iFrm]
      dImg[x,y,2,iFrm] <- rImg[x,y,2,iFrm]-meanPV2+meanImgPV2[iFrm]
    }
  }
}
thisIndex <- which(dImg<0)
dImg[thisIndex] <- 0  
# Autocorrelations and B maps
BImg <- array(0,dim=c(nPixel,nPixel,2,nFrm))
for (iFrm in (1:nFrm)){
  print(paste("ACX iFrm= ",iFrm))
  ACXdata1 <- array(0,dim=c(nPixel,nPixel,15))
  ACXdata2 <- array(0,dim=c(nPixel,nPixel,15))
  # ACXdataCC <- array(0,dim=c(nPixel,nPixel,15))
  ACXpoint1 <- array(0,dim=c(nPixel,nPixel,15))
  ACXpoint2 <- array(0,dim=c(nPixel,nPixel,15))
  # ACXpointCC <- array(0,dim=c(nPixel,nPixel,15))
  ACXGlobalCtr1 <- rep(0,15)
  ACXGlobalCtr2 <- rep(0,15)
  # ACXGlobalCtrCC <- rep(0,15)
  meanPV1 <- mean(dImg[,,1,iFrm])
  meanPV2 <- mean(dImg[,,2,iFrm])
  for (x in (1:nPixel)){
    for (y in (1:nPixel)){
      for (i in (1:15)){
        if (x+i<=nPixel){
          ACXdata1[x,y,i] <- ((dImg[x,y,1,iFrm]*dImg[(x+i),y,1,iFrm])/(meanPV1^2))-1
          ACXdata2[x,y,i] <- ((dImg[x,y,2,iFrm]*dImg[(x+i),y,2,iFrm])/(meanPV2^2))-1
          ACXGlobal1[i,iFrm] <- ACXGlobal1[i,iFrm] + ACXdata1[x,y,i]
          ACXGlobalCtr1[i] <- ACXGlobalCtr1[i] + 1
          ACXGlobal2[i,iFrm] <- ACXGlobal2[i,iFrm] + ACXdata2[x,y,i]
          ACXGlobalCtr2[i] <- ACXGlobalCtr2[i] + 1
        }
      }
    }
  }
  ACXGlobal1[,iFrm] <- ACXGlobal1[,iFrm]/ACXGlobalCtr1
  ACXGlobal2[,iFrm] <- ACXGlobal2[,iFrm]/ACXGlobalCtr2
  for (x1 in (1:nPixel)){
    for (y1 in (1:nPixel)){
      for (i in (1:15)){
        thisCtr <- 0
        for (x2 in (-7:7)){
          for (y2 in (-7:7)){
            if (x1+x2>0 && x1+x2<=nPixel && y1+y2>0 && y1+y2<=nPixel){
              ACXpoint1[x1,y1,i] <- ACXpoint1[x1,y1,i] + ACXdata1[x1+x2,y1+y2,i]
              ACXpoint2[x1,y1,i] <- ACXpoint2[x1,y1,i] + ACXdata2[x1+x2,y1+y2,i]
              thisCtr <- thisCtr + 1
            }
          }
        }
        ACXpoint1[x1,y1,i] <- ACXpoint1[x1,y1,i]/thisCtr;
        ACXpoint2[x1,y1,i] <- ACXpoint2[x1,y1,i]/thisCtr;
      }
    }
  }
  Bdata1 <- matrix(0,nPixel,nPixel)
  Bdata2 <- matrix(0,nPixel,nPixel)
  pn <- seq(1,15)
  for (x in (1:nPixel)){
    for (y in (1:nPixel)){
      acf <- ACXpoint1[x,y,]
      yIntc <- 0
      try(yIntc <- lm(acf ~ pn)$coefficients[1],silent=TRUE)
      if (yIntc<0 | !is.numeric(yIntc)){
        Bdata1[x,y] <- 0
      } else {
        Bdata1[x,y] <- yIntc*dImg[x,y,1,iFrm]
      }
      acf <- ACXpoint2[x,y,]
      yIntc <- 0
      try(yIntc <- lm(acf ~ pn)$coefficients[1],silent=TRUE)
      if (yIntc<0 | !is.numeric(yIntc)){
        Bdata2[x,y] <- 0
      } else {
        Bdata2[x,y] <- yIntc*dImg[x,y,2,iFrm]
      }
    }
  }
  BImg[,,1,iFrm] <- t(Bdata1)
  BImg[,,2,iFrm] <- t(Bdata2)
  Bdataset1[iFrm] <- mean(Bdata1)
  Bdataset2[iFrm] <- mean(Bdata2)
}
write_tif(BImg,"Bmaps.tif",overwrite=TRUE,bits_per_sample=32,msg=FALSE)

# Coincidence analysis (Pearson's correlation of B map values above mean)
for (iFrm in (1:nFrm)){
  x <- as.vector(BImg[,,1,iFrm])
  xMean <- mean(x)
  y <- as.vector(BImg[,,2,iFrm])
  yMean <- mean(y)
  thisIndex <- which(x>xMean & y>yMean)
  x <- x[thisIndex]
  y <- y[thisIndex]
  cRICS[iFrm] <- 0
  try(cRICS[iFrm] <- summary(lm(y~x))$r.squared,silent=TRUE)
}
save(Bdataset1,file="Bdataset1")
save(Bdataset2,file="Bdataset2")
save(cRICS,file="cRICS")
