
library(RpeakChrom)



## I Reading and Processing Peaks

##column
setwd("H:/Projects/RpeakChrom/DATA/CHROMOLITH/50/KBR")


t1=6.5
t2=8

peak01kbr <- readChrom("KBRF01T50A.txt", do.plot = TRUE,t1=t1,t2=t2)
peak02kbr <- readChrom("KBRF02T50A.txt", do.plot = TRUE, t1=t1/2, t2=t2/2)
peak03kbr <- readChrom("KBRF03T50A.txt", do.plot = TRUE, t1=t1/3, t2=t2/3)
#peak04kbr <- readChrom("KBRF04T50A.txt", do.plot = TRUE, t1= t1/4, t2=t2/4)
peak04kbr <- readChrom("KBRF04T50A.txt", do.plot = TRUE, t1= t1/4, t2=t2/4)
peak06kbr <- readChrom("KBRF06T50A.txt", do.plot = TRUE, t1=t1/6, t2=t2/6)
peak08kbr <- readChrom("KBRF08T50A.txt", do.plot = TRUE, t1=t1/8, t2=t2/8)
peak10kbr <- readChrom("KBRF10T50A.txt", do.plot = TRUE, t1=t1/10, t2=t2/10)
peak12kbr <- readChrom("KBRF12T50A.txt", do.plot = TRUE, t1= t1/12, t2=t2/12)
peak14kbr <- readChrom("KBRF14T50A.txt", do.plot = TRUE, t1=t1/14, t2=t2/14)
peak16kbr <- readChrom("KBRF16T50A.txt", do.plot = TRUE, t1=t1/16, t2=t2/16)
#peak20kbr <- readChrom("KBRF20T50B.txt", do.plot = TRUE, t1=t1/20+0.008, t2=t2/20-0.01)
peak20kbr <- readChrom("KBRF20T50B.txt", do.plot = TRUE, t1=t1/20, t2=t2/20-0.01)
peak25kbr <- readChrom("KBRF25T50A.txt", do.plot = TRUE, t1=t1/25, t2=t2/25)
peak30kbr <- readChrom("KBRF30T50B.txt", do.plot = TRUE, t1=0.215, t2=0.27)
#peak30kbr <- readChrom("KBRF30T50A.txt", do.plot = TRUE, t1=0.225, t2=0.26)


peaks_kbr <- list(peak01kbr,peak02kbr, peak03kbr, peak04kbr,  peak06kbr, 
                  peak08kbr, peak10kbr, peak12kbr, peak14kbr, peak16kbr, 
                  peak20kbr, peak25kbr, peak30kbr)

flows <- c(0.1,0.2, 0.3, 0.4,  0.6, 0.8, 1, 1.2, 1.4, 1.6, 2, 2.5, 3)
compounds <- rep("Kbr", length(flows))

parameters_dead <- data.frame()
for (i in 1:length(flows)){
  tryCatch({
    p <- processPeak(peaks_kbr[[i]], baseline=FALSE, flow=flows[i],
                     method="pvmg", compound=compounds[i], area=TRUE)
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



setwd("H:/Projects/RpeakChrom/DATA/CHROMOLITH/50/MER")

a=9.5
b=11.4

peak01mera <- readChrom("MERF01T50B.txt", do.plot = TRUE,t1=a,t2=b)
peak02mera <- readChrom("MERF02T50A.txt", do.plot = TRUE, t1=a/2, t2=b/2)
peak03mera <- readChrom("MERF03T50A.txt", do.plot = TRUE, t1=a/3+0.08, t2=b/3-0.2)
peak04mera <- readChrom("MERF04T50A.txt", do.plot = TRUE, t1=a/4, t2=b/4-0.08)
peak06mera <- readChrom("MERF06T50A.txt", do.plot = TRUE, t1=a/6, t2=b/6-0.08)
peak08mera <- readChrom("MERF08T50A.txt", do.plot = TRUE, t1=a/8-0.01, t2=b/8-0.04)
peak10mera <- readChrom("MERF10T50A.txt", do.plot = TRUE, t1=a/10, t2=b/10-0.04)
peak12mera <- readChrom("MERF12T50A.txt", do.plot = TRUE, t1=a/12, t2=b/12-0.04)
peak14mera <- readChrom("MERF14T50A.txt", do.plot = TRUE, t1=a/14, t2=b/14-0.04)
peak16mera <- readChrom("MERF16T50A.txt", do.plot = TRUE, t1=a/16, t2=b/16-0.03)
#peak20mera <- readChrom("MERF20T50A.txt", do.plot = TRUE, t1=a/20+0.008, t2=b/20-0.03)
peak20mera <- readChrom("MERF20T50B.txt", do.plot = TRUE, t1=a/20-0.008, t2=b/20)
peak25mera <- readChrom("MERF25T50A.txt", do.plot = TRUE, t1=a/25, t2=b/25-0.01)
peak30mera <- readChrom("MERF30T50A.txt", do.plot = TRUE, t1=a/30, t2=b/30)

peaks_mera <- list(peak01mera, peak02mera, peak03mera, peak04mera,  peak06mera, 
                   peak08mera, peak10mera, peak12mera, peak14mera, peak16mera,peak20mera,
                    peak25mera, peak30mera)


compounds <- rep("Mera", length(flows))

parameters_col_mera <- data.frame()
for (i in 1:length(flows)){
  tryCatch({
    p <- processPeak(peaks_mera[[i]], baseline=FALSE, flow=flows[i],
                     method="pvmg", compound=compounds[i], area=TRUE)
    parameters_col_mera <- rbind(parameters_col_mera, p)
  }, error=function(e){print(paste("File", i, "failed. Modify retention times with readChrom function"))})
}


parameters_col_mera


setwd("H:/Projects/RpeakChrom/DATA/CHROMOLITH/50/MIX")

a=9.8
b=11.5
c=12.3
d=14.5
e=15
f=17.5
g=19.5
h=23
peak01dimeti <- readChrom("MIXF01T50A.txt", do.plot = TRUE,t1=a,t2=b)
peak01cloro <- readChrom("MIXF01T50A.txt", do.plot = TRUE, t1=c, t2=d)
peak01sulfi <- readChrom("MIXF01T50A.txt", do.plot = TRUE, t1=e, t2=f)
peak01metoxi <- readChrom("MIXF01T50A.txt", do.plot = TRUE, t1=g, t2=h)

peak02dimeti <- readChrom("MIXF02T50A.txt", do.plot = TRUE, t1=a/2, t2=b/2)
peak02cloro <- readChrom("MIXF02T50A.txt", do.plot = TRUE, t1=c/2, t2=d/2)
peak02sulfi <- readChrom("MIXF02T50A.txt", do.plot = TRUE, t1=e/2, t2=f/2)
peak02metoxi <- readChrom("MIXF02T50A.txt", do.plot = TRUE, t1=g/2, t2=h/2)


peak03dimeti <- readChrom("MIXF03T50A.txt", do.plot = TRUE, t1=a/3, t2=b/3)
peak03cloro <- readChrom("MIXF03T50B.txt", do.plot = TRUE, t1=c/3, t2=d/3)
peak03sulfi <- readChrom("MIXF03T50A.txt", do.plot = TRUE, t1=e/3, t2=f/3)
peak03metoxi <- readChrom("MIXF03T50B.txt", do.plot = TRUE, t1=g/3, t2=h/3)


peak04dimeti <- readChrom("MIXF04T50A.txt", do.plot = TRUE, t1=a/4, t2=b/4)
peak04cloro <- readChrom("MIXF04T50B.txt", do.plot = TRUE, t1=c/4, t2=d/4)
peak04sulfi <- readChrom("MIXF04T50A.txt", do.plot = TRUE, t1=e/4, t2=f/4)
peak04metoxi <- readChrom("MIXF04T50A.txt", do.plot = TRUE, t1=g/4, t2=h/4)


peak06dimeti <- readChrom("MIXF06T50A.txt", do.plot = TRUE, t1=a/6, t2=b/6-0.01)
peak06cloro <- readChrom("MIXF06T50B.txt", do.plot = TRUE, t1=c/6, t2=d/6-0.01)
peak06sulfi <- readChrom("MIXF06T50A.txt", do.plot = TRUE, t1=e/6, t2=f/6-0.01)
peak06metoxi <- readChrom("MIXF06T50B.txt", do.plot = TRUE, t1=g/6+0.08, t2=h/6)


peak08dimeti <- readChrom("MIXF08T50B.txt", do.plot = TRUE, t1=a/8, t2=b/8-0.01)
peak08cloro <- readChrom("MIXF08T50A.txt", do.plot = TRUE, t1=c/8, t2=d/8-0.02)
peak08sulfi <- readChrom("MIXF08T50A.txt", do.plot = TRUE, t1=e/8, t2=f/8-0.02)
peak08metoxi <- readChrom("MIXF08T50B.txt", do.plot = TRUE, t1=g/8, t2=h/8)


peak10dimeti <- readChrom("MIXF10T50B.txt", do.plot = TRUE, t1=a/10, t2=b/10)
peak10cloro <- readChrom("MIXF10T50A.txt", do.plot = TRUE, t1=c/10, t2=d/10-0.03)
peak10sulfi <- readChrom("MIXF10T50A.txt", do.plot = TRUE, t1=e/10, t2=f/10-0.03)
peak10metoxi <- readChrom("MIXF10T50B.txt", do.plot = TRUE, t1=g/10-0.03, t2=h/10)

peak12dimeti <- readChrom("MIXF12T50A.txt", do.plot = TRUE, t1=a/12, t2=b/12-0.03)
peak12cloro <- readChrom("MIXF12T50A.txt", do.plot = TRUE, t1=c/12, t2=d/12)
peak12sulfi <- readChrom("MIXF12T50B.txt", do.plot = TRUE, t1=e/12, t2=f/12-0.03)
peak12metoxi <- readChrom("MIXF12T50B.txt", do.plot = TRUE, t1=g/12-0.03, t2=h/12)

peak14dimeti <- readChrom("MIXF14T50A.txt", do.plot = TRUE, t1=a/14, t2=b/14-0.03)
peak14cloro <- readChrom("MIXF14T50A.txt", do.plot = TRUE, t1=c/14, t2=d/14-0.03)
peak14sulfi <- readChrom("MIXF14T50A.txt", do.plot = TRUE, t1=e/14, t2=f/14-0.03)
peak14metoxi <- readChrom("MIXF14T50B.txt", do.plot = TRUE, t1=g/14-0.03, t2=h/14)


peak16dimeti <- readChrom("MIXF16T50B.txt", do.plot = TRUE, t1=a/16-0.02, t2=b/16+0.02)
peak16cloro <- readChrom("MIXF16T50A.txt", do.plot = TRUE, t1=c/16, t2=d/16-0.04)
peak16sulfi <- readChrom("MIXF16T50B.txt", do.plot = TRUE, t1=e/16-0.02, t2=f/16)
peak16metoxi <- readChrom("MIXF16T50B.txt", do.plot = TRUE, t1=g/16-0.03, t2=h/16)


peak20dimeti <- readChrom("MIXF20T50B.txt", do.plot = TRUE, t1=a/20+0.01, t2=b/20-0.028)
peak20cloro <- readChrom("MIXF20T50B.txt", do.plot = TRUE, t1=c/20-0.025, t2=d/20+0.01)
peak20sulfi <- readChrom("MIXF20T50B.txt", do.plot = TRUE, t1=e/20-0.02, t2=f/20)
peak20metoxi <- readChrom("MIXF20T50B.txt", do.plot = TRUE, t1=g/20-0.03, t2=h/20)


peak25dimeti <- readChrom("MIXF25T50A.txt", do.plot = TRUE, t1=a/25, t2=b/25-0.015)
peak25cloro <- readChrom("MIXF25T50A.txt", do.plot = TRUE, t1=c/25, t2=d/25-0.01)
peak25sulfi <- readChrom("MIXF25T50A.txt", do.plot = TRUE, t1=e/25, t2=f/25-0.02)
peak25metoxi <- readChrom("MIXF25T50A.txt", do.plot = TRUE, t1=g/25-0.03, t2=h/25)


peak30dimeti <- readChrom("MIXF30T50A.txt", do.plot = TRUE, t1=a/30-0.02, t2=b/30+0.03)
peak30cloro <- readChrom("MIXF30T50A.txt", do.plot = TRUE, t1=c/30, t2=d/30-0.015)
peak30sulfi <- readChrom("MIXF30T50B.txt", do.plot = TRUE, t1=e/30, t2=f/30-0.015)
peak30metoxi <- readChrom("MIXF30T50B.txt", do.plot = TRUE, t1=g/30-0.02, t2=h/30)


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



####  vanDeemter with FOLEY=TRUE



png(file="fmera_chr.png",res = 800, width = 6, height = 4, units = 'in')

fmera<- vanDeemter(col=parameters_col_mera, ext=parameters_ext, dead=parameters_dead ,length=50,
                          A=3,B=150,C=0.02, GG=FALSE, Foley=TRUE, do.plot=TRUE)
fmera$MeanError


dev.off()
png(file="fmetoxi_chr.png",res = 800, width = 6, height = 4, units = 'in')

fmetoxi<- vanDeemter(col=parameters_col_metoxi, ext=parameters_ext, dead=parameters_dead ,length=50,
                             A=3,B=150,C=0.02, GG=FALSE, Foley=TRUE, do.plot=TRUE)
fmetoxi$MeanError


dev.off()
png(file="fsulfi_chr.png",res = 800, width = 6, height = 4, units = 'in')

fsulfi<- vanDeemter(col=parameters_col_sulfi, ext=parameters_ext, dead=parameters_dead ,length=50,
                           A=3,B=150,C=0.02, GG=FALSE, Foley=TRUE, do.plot=TRUE)
fsulfi$MeanError


dev.off()
png(file="fcloro_chr.png",res = 800, width = 6, height = 4, units = 'in')

fcloro<- vanDeemter(col=parameters_col_cloro, ext=parameters_ext, dead=parameters_dead ,length=50,
                           A=3,B=150,C=0.02, GG=FALSE, Foley=TRUE, do.plot=TRUE)

fcloro$MeanError


dev.off()
png(file="fdimeti_chr.png",res = 800, width = 6, height = 4, units = 'in')

fdimeti<- vanDeemter(col=parameters_col_dimeti, ext=parameters_ext, dead=parameters_dead ,length=50,
                            A=3,B=150,C=0.02, GG=FALSE, Foley=TRUE, do.plot=TRUE)

dev.off()





fdimeti$MeanError

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


#### example vanDeemter function with GG=TRUE




dev.off()
png(file="ggmera_chr.png",res = 800, width = 6, height = 4, units = 'in')


ggmera<- vanDeemter(col=parameters_col_mera, ext=parameters_ext, dead=parameters_dead ,length=50,
                          A=6,B=200,C=0.04, Foley=FALSE, GG=TRUE, do.plot=TRUE)



dev.off()
png(file="ggmetoxi_chr.png",res = 800, width = 6, height = 4, units = 'in')



ggmetoxi<- vanDeemter(col=parameters_col_metoxi, ext=parameters_ext, dead=parameters_dead ,length=50,
                            A=6,B=200,C=0.04, Foley=FALSE, GG=TRUE, do.plot=TRUE)



dev.off()
png(file="ggsulfi_chr.png",res = 800, width = 6, height = 4, units = 'in')

ggsulfi<- vanDeemter(col=parameters_col_sulfi, ext=parameters_ext, dead=parameters_dead ,length=50,
                           A=6,B=200,C=0.04, Foley=FALSE, GG=TRUE, do.plot=TRUE)






dev.off()
png(file="ggcloro_chr.png",res = 800, width = 6, height = 4, units = 'in')

ggcloro<- vanDeemter(col=parameters_col_cloro, ext=parameters_ext, dead=parameters_dead ,length=50,
                           A=6,B=200,C=0.04, Foley=FALSE, GG=TRUE, do.plot=TRUE)






dev.off()
png(file="ggdimeti_chr.png",res = 800, width = 6, height = 4, units = 'in')

ggdimeti<- vanDeemter(col=parameters_col_dimeti, ext=parameters_ext, dead=parameters_dead ,length=50,
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

png(file="ApproachIb.png",res = 800, width = 6, height = 4, units = 'in')

ApproachI <-  vanDeemterAlternative(col=col, ext=parameters_ext, dead=parameters_dead, 
                                    length=50, approachI=TRUE, A=6, B=200, C=0.04, approachII=FALSE)

dev.off()


##  The coefficients parameters  corresponding to Approach I
ApproachI$correlation
ApproachI$coefficients
ApproachI$step1
ApproachI$MeanError
ApproachI$RSE

png(file="ApproachIIb.png",res = 800, width = 6, height = 4, units = 'in')

ApproachII <- vanDeemterAlternative(col=col, ext=parameters_ext, dead=parameters_dead, 
                                    length=50, approachI=FALSE, A=6, B=200, C=0.04, approachII=TRUE)

dev.off()
##  The coefficients parameters corresponding to Approach II
ApproachII$r2C
ApproachII$coefficients
MREA<-ApproachII$MREA
MREB<-ApproachII$MREB
MREC<-ApproachII$MREC
(MREA+MREB+MREC)/3
MREA
MREB
MREC
#################### end

