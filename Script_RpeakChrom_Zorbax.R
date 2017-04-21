
library(RpeakChrom)



### link to the R data for RpeakChrom package

#    https://github.com/maialba3/RpeakChromData


############ STARTS THE SCRIPT

## I Reading and Processing Peaks

##column


setwd("H:/Projects/RpeakChrom/DATA/ZORBAX/KBR")



peak01kbr <- readChrom("KBRF01A.txt", do.plot = TRUE, t1=12, t2=13.5)
peak02kbr <- readChrom("KBRF02A.txt", do.plot = TRUE, t1=6, t2=6.7)
peak03kbr <- readChrom("KBRF03A.txt", do.plot = TRUE, t1=3.95, t2=4.4)
peak04kbr <- readChrom("KBRF04A.txt", do.plot = TRUE, t1= 2.95, t2=3.35)
peak06kbr <- readChrom("KBRF06A.txt", do.plot = TRUE, t1=1.95, t2=2.2)
peak08kbr <- readChrom("KBRF08A.txt", do.plot = TRUE, t1=1.455, t2=1.65)
peak10kbr <- readChrom("KBRF10A.txt", do.plot = TRUE, t1=1.16, t2=1.32)
peak12kbr <- readChrom("KBRF12A.txt", do.plot = TRUE, t1= 0.98, t2=1.10)
peak14kbr <- readChrom("KBRF14A.txt", do.plot = TRUE, t1=0.85, t2=0.95)
peak16kbr <- readChrom("KBRF16A.txt", do.plot = TRUE, t1=0.72, t2=0.84)
peak20kbr <- readChrom("KBRF20A.txt", do.plot = TRUE, t1=0.59, t2=0.66)
peak25kbr <- readChrom("KBRF25A.txt", do.plot = TRUE, t1=0.46, t2=0.54)


#peak30kbr <- readChrom("KBRF30A.txt", do.plot = TRUE, t1=0.389, t2=0.455)

peak30kbr <- readChrom("KBRF30A.txt", do.plot = TRUE, t1=0.399, t2=0.435)



peaks_kbr <- list(peak01kbr,peak02kbr, peak03kbr, peak04kbr,  peak06kbr, 
                  peak08kbr, peak10kbr, peak12kbr, peak14kbr, peak16kbr, 
                  peak20kbr, peak25kbr, peak30kbr)

flows <- c(0.1,0.2, 0.3, 0.4,  0.6, 0.8, 1, 1.2, 1.4, 1.6, 2, 2.5, 3)
compounds <- rep("Kbr", length(flows))

parameters_dead <- data.frame()
for (i in 1:length(flows)){
  tryCatch({
    p <- processPeak(peaks_kbr[[i]], baseline=FALSE, flow=flows[i],
                     method="pvmg", compound=compounds[i], area=FALSE)
    parameters_dead <- rbind(parameters_dead, p)
  }, error=function(e){print(paste("File", i, "failed. Modify retention times with readChrom function"))})
}

parameters_dead


## Extracolumn


setwd("H:/Projects/RpeakChrom/DATA/EXTRACOLUMN")

a<-0.2
b<-1

peak01kbr <- readChrom("KBREF01.txt", do.plot = TRUE,t1=a,t2=b)
peak02kbr <- readChrom("KBREF02.txt", do.plot = TRUE, t1=a/2, t2=b/2)
peak03kbr <- readChrom("KBREF03.txt", do.plot = TRUE, t1=a/3, t2=b/3)
peak04kbr <- readChrom("KBREF04.txt", do.plot = TRUE, t1= a/4, t2=b/4)
peak06kbr <- readChrom("KBREF06.txt", do.plot = TRUE, t1=a/6, t2=b/6)
peak08kbr <- readChrom("KBREF08.txt", do.plot = TRUE, t1=a/8, t2=b/8)
peak10kbr <- readChrom("KBREF10.txt", do.plot = TRUE, t1=a/10, t2=b/10)
peak12kbr <- readChrom("KBREF12.txt", do.plot = TRUE, t1= a/12, t2=b/12)
peak14kbr <- readChrom("KBREF14.txt", do.plot = TRUE, t1=a/14+0.0045 ,t2=b/14-0.01)
peak16kbr <- readChrom("KBREF16.txt", do.plot = TRUE, t1=a/16, t2=b/16)
peak20kbr <- readChrom("KBREF20.txt", do.plot = TRUE, t1=a/20+0.0065, t2=b/20-0.015)
peak25kbr <- readChrom("KBREF25.txt", do.plot = TRUE, t1=a/25+0.0065, t2=b/25-0.001)
peak30kbr <- readChrom("KBREF30.txt", do.plot = TRUE, t1=0.0025+0.011, t2=0.04-0.007)


peaks_kbr_e <- list(peak01kbr,peak02kbr, peak03kbr, peak04kbr,  peak06kbr, 
                  peak08kbr, peak10kbr, peak12kbr, peak14kbr, peak16kbr,  
                  peak20kbr, peak25kbr, peak30kbr)

flows <- c(0.1,0.2, 0.3, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 2, 2.5, 3)
compounds <- rep("Kbre", length(flows))

parameters_ext<- data.frame()
for (i in 1:length(flows)){
  tryCatch({
    p <- processPeak(peaks_kbr_e[[i]], baseline=FALSE, flow=flows[i],
                     method="pvmg", compound=compounds[i], area=TRUE)
    parameters_ext <- rbind(parameters_ext, p)
  }, error=function(e){print(paste("File", i, "failed. Modify retention times with readChrom function"))})
}


setwd("H:/Projects/RpeakChrom/DATA/ZORBAX/MERACOL")

peak01mera <- readChrom("MERF01A.txt", do.plot = TRUE,  t1=26.5, t2=29)
peak02mera <- readChrom("MERF02A.txt", do.plot = TRUE, t1=13.2, t2=14.5)
peak03mera <- readChrom("MERF03A.txt", do.plot = TRUE, t1=8.8, t2=9.5)
peak04mera <- readChrom("MERF04A.txt", do.plot = TRUE, t1=6.5, t2=7.2)
peak06mera <- readChrom("MERF06A.txt", do.plot = TRUE, t1=4.3, t2=4.7)
peak08mera <- readChrom("MERF08A.txt", do.plot = TRUE, t1=3.25, t2=3.55)
peak10mera <- readChrom("MERF10A.txt", do.plot = TRUE, t1=2.57, t2=2.85)
peak12mera <- readChrom("MERF12A.txt", do.plot = TRUE, t1=2.15, t2=2.35)
peak14mera <- readChrom("MERF14A.txt", do.plot = TRUE, t1=1.85, t2=2.02)
peak16mera <- readChrom("MERF16A.txt", do.plot = TRUE, t1=1.61, t2=1.78)
peak20mera <- readChrom("MERF20A.txt", do.plot = TRUE, t1=1.29, t2=1.43)
peak25mera <- readChrom("MERF25A.txt", do.plot = TRUE, t1=1.03, t2=1.15)
peak30mera <- readChrom("MERF30A.txt", do.plot = TRUE, t1=0.87, t2=0.98)

peaks_mera <- list(peak01mera, peak02mera, peak03mera, peak04mera,  peak06mera, 
                   peak08mera, peak10mera, peak12mera, peak14mera, peak16mera, 
                   peak20mera, peak25mera, peak30mera)


compounds <- rep("Sulphamerazine", length(flows))

parameters_col_mera <- data.frame()
for (i in 1:length(flows)){
  tryCatch({
    p <- processPeak(peaks_mera[[i]], baseline=FALSE, flow=flows[i],
                     method="pvmg", compound=compounds[i], area=TRUE)
    parameters_col_mera <- rbind(parameters_col_mera, p)
  }, error=function(e){print(paste("File", i, "failed. Modify retention times with readChrom function"))})
}





setwd("H:/Projects/RpeakChrom/DATA/ZORBAX/MIXCOL")


peak01dimeti <- readChrom("MIXF01A.txt", do.plot = TRUE, t1=28, t2=31)
peak01cloro <- readChrom("MIXF01A.txt", do.plot = TRUE, t1=41.9, t2=47)
peak01sulfi <- readChrom("MIXF01A.txt", do.plot = TRUE, t1=57.5, t2=62.5)
peak01metoxi <- readChrom("MIXF01B.txt", do.plot = TRUE, t1=84, t2=93)

peak02dimeti <- readChrom("MIXF02B.txt", do.plot = TRUE, t1=28/2, t2=31/2)
peak02cloro <- readChrom("MIXF02B.txt", do.plot = TRUE, t1=41.9/2+0.05, t2=23)
peak02sulfi <- readChrom("MIXF02B.txt", do.plot = TRUE, t1=57.5/2, t2=31)
peak02metoxi <- readChrom("MIXF02B.txt", do.plot = TRUE, t1=84/2, t2=93/2)


peak03dimeti <- readChrom("MIXF03B.txt", do.plot = TRUE, t1=28/3, t2=10.2)
peak03cloro <- readChrom("MIXF03B.txt", do.plot = TRUE, t1=41.9/3, t2=15.2)
peak03sulfi <- readChrom("MIXF03B.txt", do.plot = TRUE, t1=19, t2=20.4)
peak03metoxi <- readChrom("MIXF03A.txt", do.plot = TRUE, t1=28.8, t2=30)


peak04dimeti <- readChrom("MIXF04A.txt", do.plot = TRUE, t1=28/4, t2=7.6)
peak04cloro <- readChrom("MIXF04A.txt", do.plot = TRUE, t1=41.9/4-0.03, t2=11.4)
peak04sulfi <- readChrom("MIXF04A.txt", do.plot = TRUE, t1=14.2, t2=15.3)
peak04metoxi <- readChrom("MIXF04B.txt", do.plot = TRUE, t1=20.8, t2=22.6)


peak06dimeti <- readChrom("MIXF06A.txt", do.plot = TRUE, t1=28/6, t2=5)
peak06cloro <- readChrom("MIXF06A.txt", do.plot = TRUE, t1=41.9/6-0.08, t2=7.4+0.09)
peak06sulfi <- readChrom("MIXF06A.txt", do.plot = TRUE, t1=9.4, t2=10)

peak06metoxi <- readChrom("MIXF06A.txt", do.plot = TRUE, t1=13.7, t2=15.)
#peak06metoxi <- readChrom("MIXF06A.txt", do.plot = TRUE, t1=14, t2=14.9)

peak08dimeti <- readChrom("MIXF08A.txt", do.plot = TRUE, t1=3.49, t2=3.70)
peak08cloro <- readChrom("MIXF08A.txt", do.plot = TRUE, t1=41.9/8-0.15, t2=5.60)
peak08sulfi <- readChrom("MIXF08A.txt", do.plot = TRUE, t1=6.95, t2=7.4)
peak08metoxi <- readChrom("MIXF08B.txt", do.plot = TRUE, t1=10.25, t2=11.1)


peak10dimeti <- readChrom("MIXF10A.txt", do.plot = TRUE, t1=28/10, t2=2.98)
peak10cloro <- readChrom("MIXF10B.txt", do.plot = TRUE, t1=41.9/10-0.15, t2=4.5-0.05)
peak10sulfi <- readChrom("MIXF10A.txt", do.plot = TRUE, t1=5.56, t2=5.95)
#peak10metoxi <- readChrom("MIXF10A.txt", do.plot = TRUE, t1=8.3, t2=8.9)
peak10metoxi <- readChrom("MIXF10B.txt", do.plot = TRUE, t1=8.10, t2=8.85)


peak12dimeti <- readChrom("MIXF12A.txt", do.plot = TRUE, t1=2.3, t2=2.5)
peak12cloro <- readChrom("MIXF12A.txt", do.plot = TRUE, t1=41.9/12-0.12, t2=3.73)
peak12sulfi <- readChrom("MIXF12A.txt", do.plot = TRUE, t1=4.5, t2=5)
peak12metoxi <- readChrom("MIXF12B.txt", do.plot = TRUE, t1=6.65, t2=7.4)

peak14dimeti <- readChrom("MIXF14A.txt", do.plot = TRUE, t1=1.98, t2=2.11)
peak14cloro <- readChrom("MIXF14A.txt", do.plot = TRUE, t1=2.85, t2=3.2)
peak14sulfi <- readChrom("MIXF14B.txt", do.plot = TRUE, t1=3.82, t2=4.25)
peak14metoxi <- readChrom("MIXF14B.txt", do.plot = TRUE, t1=5.7, t2=6.3)


peak16dimeti <- readChrom("MIXF16A.txt", do.plot = TRUE, t1=1.73, t2=1.86)
peak16cloro <- readChrom("MIXF16A.txt", do.plot = TRUE, t1=2.50, t2=2.8)
peak16sulfi <- readChrom("MIXF16B.txt", do.plot = TRUE, t1=3.35, t2=3.7)
peak16metoxi <- readChrom("MIXF16B.txt", do.plot = TRUE, t1=4.95, t2=5.5)


peak20dimeti <- readChrom("MIXF20A.txt", do.plot = TRUE, t1=1.38, t2=1.5)
peak20cloro <- readChrom("MIXF20A.txt", do.plot = TRUE, t1=2.03, t2=2.22)
peak20sulfi <- readChrom("MIXF20A.txt", do.plot = TRUE, t1=2.7, t2=3)
peak20metoxi <- readChrom("MIXF20A.txt", do.plot = TRUE, t1=4.05, t2=4.5)


peak25dimeti <- readChrom("MIXF25A.txt", do.plot = TRUE, t1=1.1, t2=1.2)
peak25cloro <- readChrom("MIXF25A.txt", do.plot = TRUE, t1=1.62, t2=1.8)
peak25sulfi <- readChrom("MIXF25A.txt", do.plot = TRUE, t1=2.15, t2=2.45)
peak25metoxi <- readChrom("MIXF25B.txt", do.plot = TRUE, t1=3.20, t2=3.65)


peak30dimeti <- readChrom("MIXF30A.txt", do.plot = TRUE, t1=0.9, t2=31/30)
peak30cloro <- readChrom("MIXF30A.txt", do.plot = TRUE, t1=1.33, t2=1.52)
peak30sulfi <- readChrom("MIXF30A.txt", do.plot = TRUE, t1=1.82, t2=2.04)
#peak30metoxi <- readChrom("MIXF30A.txt", do.plot = TRUE, t1=2.8, t2=2.98)
peak30metoxi <- readChrom("MIXF30A.txt", do.plot = TRUE, t1=2.70, t2=3.05)

peaks_dimeti <- list(peak01dimeti, peak02dimeti, peak03dimeti, peak04dimeti,  peak06dimeti, 
                   peak08dimeti, peak10dimeti, peak12dimeti, peak14dimeti, peak16dimeti, 
                   peak20dimeti, peak25dimeti, peak30dimeti)

compounds <- rep("dimeti", length(flows))
parameters_col_dimeti <- data.frame()
for (i in 1:length(flows)){
  tryCatch({
    p <- processPeak(peaks_dimeti[[i]], baseline=FALSE, flow=flows[i],
                     method="pvmg", compound=compounds[i], area=TRUE)
    parameters_col_dimeti <- rbind(parameters_col_dimeti, p)
  }, error=function(e){print(paste("File", i, "failed. Modify retention times with readChrom function"))})
}




peaks_cloro <- list(peak01cloro, peak02cloro, peak03cloro, peak04cloro, peak06cloro, 
                     peak08cloro, peak10cloro, peak12cloro, peak14cloro, peak16cloro, 
                     peak20cloro, peak25cloro, peak30cloro)

compounds <- rep("cloro", length(flows))
parameters_col_cloro <- data.frame()
for (i in 1:length(flows)){
  tryCatch({
    p <- processPeak(peaks_cloro[[i]], baseline=FALSE, flow=flows[i],
                     method="pvmg", compound=compounds[i], area=TRUE)
    parameters_col_cloro <- rbind(parameters_col_cloro, p)
  }, error=function(e){print(paste("File", i, "failed. Modify retention times with readChrom function"))})
}


peaks_sulfi <- list(peak01sulfi, peak02sulfi, peak03sulfi, peak04sulfi,  peak06sulfi, 
                    peak08sulfi, peak10sulfi, peak12sulfi, peak14sulfi, peak16sulfi, 
                    peak20sulfi, peak25sulfi, peak30sulfi)


compounds <- rep("sulfi", length(flows))
parameters_col_sulfi <- data.frame()
for (i in 1:length(flows)){
  tryCatch({
    p <- processPeak(peaks_sulfi[[i]], baseline=FALSE, flow=flows[i],
                     method="pvmg", compound=compounds[i], area=TRUE)
    parameters_col_sulfi <- rbind(parameters_col_sulfi, p)
  }, error=function(e){print(paste("File", i, "failed. Modify retention times with readChrom function"))})
}

peaks_metoxi <- list(peak01metoxi, peak02metoxi, peak03metoxi, peak04metoxi,  peak06metoxi, 
                    peak08metoxi,peak10metoxi, peak12metoxi, peak14metoxi, peak16metoxi, 
                    peak20metoxi, peak25metoxi, peak30metoxi)


compounds <- rep("metoxi", length(flows))
parameters_col_metoxi <- data.frame()
for (i in 1:length(flows)){
  tryCatch({
    p <- processPeak(peaks_metoxi[[i]], baseline=FALSE, flow=flows[i],
                     method="pvmg", compound=compounds[i], area=TRUE)
    parameters_col_metoxi <- rbind(parameters_col_metoxi, p)
  }, error=function(e){print(paste("File", i, "failed. Modify retention times with readChrom function"))})
}




col<-rbind(parameters_col_metoxi,
parameters_col_sulfi,
parameters_col_cloro,
parameters_col_dimeti,
parameters_col_mera)


### Peak Assymetry
parameters_col_dimeti
parameters_col_mera
assy_dimeti_10<-parameters_col_dimeti$B10/parameters_col_dimeti$A10
assy_dimeti_60<-parameters_col_dimeti$B60/parameters_col_dimeti$A60

assy_dimeti_10<-mean(assy_dimeti_10)
assy_dimeti_60<-mean(assy_dimeti_60)

assy_mera_10<-parameters_col_mera$B10/parameters_col_mera$A10
assy_mera_60<-parameters_col_mera$B60/parameters_col_mera$A60

assy_mera_10<-mean(assy_mera_10)
assy_mera_60<-mean(assy_mera_60)

assy_sulfi_10<-parameters_col_sulfi$B10/parameters_col_sulfi$A10
assy_sulfi_60<-parameters_col_sulfi$B60/parameters_col_sulfi$A60

assy_sulfi_10<-mean(assy_sulfi_10)
assy_sulfi_60<-mean(assy_sulfi_60)

assy_cloro_10<-parameters_col_cloro$B10/parameters_col_cloro$A10
assy_cloro_60<-parameters_col_cloro$B60/parameters_col_cloro$A60

assy_cloro_10<-mean(assy_cloro_10)
assy_cloro_60<-mean(assy_cloro_60)

assy_metoxi_10<-parameters_col_cloro$B10/parameters_col_cloro$A10
assy_metoxi_60<-parameters_col_cloro$B60/parameters_col_cloro$A60

assy_metoxi_10<-mean(assy_cloro_10)
assy_metoxi_60<-mean(assy_cloro_60)


#### vanDeemter with FOLEY=TRUE



png(file="fmera.png",res = 800, width = 6, height = 4, units = 'in')

fmera<- vanDeemter(col=parameters_col_mera, ext=parameters_ext, dead=parameters_dead ,length=150,
                   A=3,B=150,C=0.02, GG=FALSE, Foley=TRUE, do.plot=TRUE)
fmera$MeanError


dev.off()
png(file="fmetoxi.png",res = 800, width = 6, height = 4, units = 'in')

fmetoxi<- vanDeemter(col=parameters_col_metoxi, ext=parameters_ext, dead=parameters_dead ,length=150,
                     A=3,B=150,C=0.02, GG=FALSE, Foley=TRUE, do.plot=TRUE)
fmetoxi$MeanError


dev.off()
png(file="fsulfi.png",res = 800, width = 6, height = 4, units = 'in')

fsulfi<- vanDeemter(col=parameters_col_sulfi, ext=parameters_ext, dead=parameters_dead ,length=150,
                    A=3,B=150,C=0.02, GG=FALSE, Foley=TRUE, do.plot=TRUE)
fsulfi$MeanError


dev.off()
png(file="fcloro.png",res = 800, width = 6, height = 4, units = 'in')

fcloro<- vanDeemter(col=parameters_col_cloro, ext=parameters_ext, dead=parameters_dead ,length=150,
                    A=3,B=150,C=0.02, GG=FALSE, Foley=TRUE, do.plot=TRUE)

fcloro$MeanError


dev.off()
png(file="fdimeti.png",res = 800, width = 6, height = 4, units = 'in')

fdimeti<- vanDeemter(col=parameters_col_dimeti, ext=parameters_ext, dead=parameters_dead ,length=150,
                     A=3,B=150,C=0.02, GG=FALSE, Foley=TRUE, do.plot=TRUE)
fdimeti$MeanError
dev.off()


fall<-rbind(fdimeti[[1]],fcloro[[1]],fsulfi[[1]],fmetoxi[[1]],fmera[[1]])
A<-mean(fall[,1])
B<-mean(fall[,2])
C<-mean(fall[,3])

## The coefficients parameters corresponding to Van Deemter Classic (Foley)
A
B
C

MEvandeemter<-rbind(fdimeti$MeanError,fcloro$MeanError,fsulfi$MeanError,fmera$MeanError,fmetoxi$MeanError)
MEvandeemter<-mean(MEvandeemter)
MEvandeemter


####  vanDeemter function with GG=TRUE




png(file="ggmera.png",res = 800, width = 6, height = 4, units = 'in')


ggmera<- vanDeemter(col=parameters_col_mera, ext=parameters_ext, dead=parameters_dead ,length=150,
                    A=6,B=200,C=0.04, Foley=FALSE, GG=TRUE, do.plot=TRUE)



dev.off()
png(file="ggmetoxi.png",res = 800, width = 6, height = 4, units = 'in')



ggmetoxi<- vanDeemter(col=parameters_col_metoxi, ext=parameters_ext, dead=parameters_dead ,length=150,
                      A=6,B=200,C=0.04, Foley=FALSE, GG=TRUE, do.plot=TRUE)



dev.off()
png(file="ggsulfi.png",res = 800, width = 6, height = 4, units = 'in')

ggsulfi<- vanDeemter(col=parameters_col_sulfi, ext=parameters_ext, dead=parameters_dead ,length=150,
                     A=6,B=200,C=0.04, Foley=FALSE, GG=TRUE, do.plot=TRUE)






dev.off()
png(file="ggcloro.png",res = 800, width = 6, height = 4, units = 'in')

ggcloro<- vanDeemter(col=parameters_col_cloro, ext=parameters_ext, dead=parameters_dead ,length=150,
                     A=6,B=200,C=0.04, Foley=FALSE, GG=TRUE, do.plot=TRUE)






dev.off()
png(file="ggdimeti.png",res = 800, width = 6, height = 4, units = 'in')

ggdimeti<- vanDeemter(col=parameters_col_dimeti, ext=parameters_ext, dead=parameters_dead ,length=150,
                      A=6,B=200,C=0.04, Foley=FALSE, GG=TRUE, do.plot=TRUE)


dev.off()

ggall<-rbind(ggdimeti[[1]],ggcloro[[1]],ggsulfi[[1]],ggmetoxi[[1]],ggmera[[1]])


Aa<-mean(ggall[,1])
Bb<-mean(ggall[,2])
Cc<-mean(ggall[,3])

## The coefficients parameters  corresponding to Van Deemter Classic (GG)
Aa
Bb
Cc


MEvandeemter<-rbind(ggdimeti$MeanError,ggcloro$MeanError,ggsulfi$MeanError,ggmera$MeanError,ggmetoxi$MeanError)
MEvandeemter<-mean(MEvandeemter)
MEvandeemter

#### vanDeemterAlternative function 


png(file="ApproachI_zor.png",res = 800, width = 6, height = 4, units = 'in')

ApproachI <-  vanDeemterAlternative(col=col, ext=parameters_ext, dead=parameters_dead, 
                                    length=150, approachI=TRUE, A=6, B=200, C=0.04, approachII=FALSE)

dev.copy(png,"myplot.png",res=800,width=6,height=4,units='in')

dev.off()



##
ApproachI$correlation
ApproachI$coefficients
ApproachI$step1
ApproachI$MeanError
ApproachI$RSE







png(file="ApproachIIb_zor.png",res = 800, width = 6, height = 4, units = 'in')

ApproachII <- vanDeemterAlternative(col=col, ext=parameters_ext, dead=parameters_dead, 
                                    length=150, approachI=FALSE, A=6, B=200, C=0.04, approachII=TRUE)


dev.off()
##  The coefficients parameters  corresponding to Approach II

ApproachII$coefficients
MREA<-ApproachII$MREA
MREB<-ApproachII$MREB
MREC<-ApproachII$MREC
(MREA+MREB+MREC)/3
MREA
MREB
MREC


#################### END THE SCRIPT

