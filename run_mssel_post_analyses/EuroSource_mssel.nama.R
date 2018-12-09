setwd("//share/hennlab/data/SA_sim/sim/nama/Euro_mig/")
library(rehh.data, lib.loc="/share/hennlab/R/x86_64-pc-linux-gnu-library/3.5")
library(gplots, lib.loc="/share/hennlab/R/x86_64-pc-linux-gnu-library/3.5")
library(rehh,  lib.loc="/share/hennlab/R/x86_64-pc-linux-gnu-library/3.5")


################
#   Function   #
################
var.test<-function(locus){ # test if a locus is segregating in samples
  if (length(levels(as.factor(locus)))==1) {
    seg<-0
  }else(seg<-1)
  return(seg)
}

SegsiteOut <- function(haps, sel.site){ # output haps that only contain variants
  nsam<-dim(haps)[1]
  index.var<-apply(haps, FUN=var.test, MARGIN = 2)
  var.num <- sum(index.var)
  sel.new.pos <- sum(index.var[c(1:sel.site)])
  haps.sub <- haps[,as.logical(index.var)]
  return(list(var.num, sel.new.pos, haps.sub)) # number of true seg sites, selected allele's new position(a fake number), and variants
}

SeqErr <- function(haps, sel.site){ 
  nsam <- dim(haps)[1]
  nsite<-dim(haps)[2]
  store.selsite <- haps[,sel.site] # temporary storage of selected site
  index.var <- apply(haps[which(haps[,sel.site]=="1"),], FUN=var.test, MARGIN=2)
  #nvar <- sum(index.var) # number of variant sites at derived haplotypes 
  #lambda.FN <- 0.04*nvar
  #lambda.FP <- 0.02*nvar
  for (i in 1:nsam){
    index.ref <- which(haps[i,]=="0")
    index.alt <- which(haps[i,]=="1")
    n.alt <- length(index.alt)
    if(n.alt==0){
    next
    }
    lambda.FN <- 0.04*n.alt
    lambda.FP <- 0.02*n.alt
    n.FP <- rpois(1, lambda.FP)
    n.FN <-rpois(1, lambda.FN)
    
    if(n.FP>length(index.ref)){
    index.flip.FP <- index.ref
    }else{
    index.flip.FP <- index.ref[sample(length(index.ref), size=n.FP, replace=F)]
    }
    
    if(n.FN > length(index.alt)){
    index.flip.FN <- index.alt
    }else{
    index.flip.FN <- index.alt[sample(length(index.alt), size=n.FN, replace=F)]
    }
    index.flip <- c(index.flip.FP, index.flip.FN)
    
    if(length(index.flip)!=0){
    haps[i,index.flip] = as.character(abs(as.numeric(haps[i,index.flip])-1))
    }
  }
  haps[,sel.site] <- store.selsite # making sure selsite is not flipped
  return(haps) 
}


CalcPi <- function(haps){ # calculate Tajima's theta: mean pairwise differences
  nsam <- dim(haps)[1]
  #freq.vec <- unlist(lapply(apply(haps, FUN=table, MARGIN = 2),`[[`,1))/nsam# allele freq at all loci
  try(count.vec <- unlist(lapply(strsplit(as.vector(summary(haps)[1,]), split="0:"),`[[`,2)))
  if(exists("count.vec")){
  freq.vec <- as.numeric(count.vec)/nsam
  pi <- sum(nsam / (nsam-1) * 2*freq.vec*(1-freq.vec))
  }else{
       nsam <-length(haps)
          count <- length(which(haps=="0"))
          freq <- count/nsam
          pi <- nsam/(nsam-1)*2*freq*(1-freq)
  }
  return(pi)
}

RecRate<-function(k, mean){
theta <- mean/(k-1)
rate <- rgamma(1, shape=k, scale=theta)
return(rate)
}

RecodeRefAlt <- function(haps,selsite){ # recode 0 and 1 ans ref and alt, as compared to a Euro hap, the last one is assumed to be euro
 n.hap <- dim(haps)[1]
 index.recode <- which(haps[n.hap,]=="1") # which alleles in the Euro hap are coded as 1, while we want them to be present as just 0 (ref)
 for (i in c(1:(n.hap-1))) {
    haps[i,index.recode] <- as.character(abs(as.numeric(haps[i,index.recode])-1)) # flip all alleles at these loci, for all San haplotypes
 }
 #attention: the derived rs142 is recoded as 0 now, flip this site back
 haps[,selsite] <- as.character(abs(as.numeric(haps[,selsite])-1))
 
 return(haps[c(1:(n.hap-1)),]) #only return San haplotypes
}



SumStat <- function(TrajFile){
  # before running mssel, get demographic parameters
  name.seg <- strsplit(TrajFile, split="_")
  sel <- as.numeric(name.seg[[1]][2])
  rm(name.seg)
  
  Ne = 10000 # current Ne of KhoeSan #!
  #r.rate <- RecRate(2.5, 6e-9) # recombination rate per site per generation
  r.rate <- RecRate(1.1, 6e-10)
  mu <- 1.9e-8 # mutation rate per generation
  nanc.K <- 161*2 - 172 #number of ancestral copies in KhoeSan #!
  nsel.K <- 172 #!
  nanc.E <- 0
  nsel.E <- 198
  nanc.B <- 198
  nsel.B <- 0
  nanc.total <-nanc.K +nanc.E +nanc.B
  nsel.total <- nsel.K + nsel.E + nsel.B
  nsam.total <- nanc.total + nsel.total
  
  #other fixed demographic paramters
  tEK <- 10  #!
  mEK <- 0.17 #!
  tBK <- 14
  mBK <- 0.02 #!
  NE <- 10000
  NB <- 17000
  
  
  
  
  #bottleneck timing: 20/4/10000
  #run mssel
  system(paste("/share/hennlab/progs/Meng_mssel_programs/mssel ", nsam.total, " 1 ", nanc.total, " ", nsel.total, " /share/hennlab/data/SA_sim/trajectories/nama/Nama_EuroMigOnly_traj/", TrajFile, 
  " 23316 -r ", 4*Ne*31700*r.rate, " 31700 -t ", 4*Ne*31700*mu, " -I 3 ", nanc.E, " ", nsel.E, " ", nanc.K, " ", nsel.K, " ", nanc.B, " ", nsel.B, 
  " -n 1 ", NE/Ne, " -n 2 1.0 -n 3 ", NB/Ne, 
  " -es ", tEK/4/Ne, " 2 ", 1-mEK, " -ej ", tEK/4/Ne, " 4 1 -en ", tEK/4/Ne, " 1 ", NE/Ne, 
  " -es ", tBK/4/Ne, " 2 ", 1-mBK, " -ej ", tBK/4/Ne, " 5 3 -en ", tBK/4/Ne, " 3 ", NB/Ne, 
  " -en 0.0005 2 2.0", #nama bottleneck #!
  " -en ", 1800/4/Ne, " 1 ", 1000/Ne, #OOA bottleneck
  " -ej ", 75000/30/4/Ne, " 1 3", # divergence between Bantu and Europeans,
  " -ej ", 110000/30/4/Ne, " 3 2 > ", TrajFile, ".mssel", sep="")) # divergence between San and others
  
  #split mssel to derived and ancestral for calc on s and pi, EHH not affected
  ms.all <- read.table(paste(TrajFile,".mssel",sep=""), sep = "", header = F, colClasses = "character", skip = 8)
  ms.output<- as.matrix(ms.all[c(c((nanc.E+nsel.E+1):(nanc.E+nsel.E+nanc.K+nsel.K)), 1),]) # Keep San haplotypes, and one European haplotype as reference
  ms.intermediate<-apply(data.frame(as.character(ms.output[,1])),MARGIN=1,FUN = strsplit, split="")
  haps<-matrix(NA, nrow=(nanc.K+nsel.K+1),ncol=nchar(ms.output[1,1]))
  ##reformating haplotypes to a matrix
  for(i in 1:(nanc.K+nsel.K+1)) {
    haps[i,]<- ms.intermediate[[i]][[1]]
  }
  rm(ms.intermediate)
  ##find the selected site
  system(paste("head -n 6 ", TrajFile, ".mssel | tail -n 1 > ", TrajFile, ".selsite", sep="")) # store selected site
  selsite.file<-read.table(paste(TrajFile,".selsite", sep=""))
  selsite <- selsite.file$V2[1] +1
  
  #now seperate San haplotypes and the one Euro reference hap, and recode it as ref and alt (except for the selected site)
  haps.san <- RecodeRefAlt(haps, selsite)
  haps <- haps.san 
  rm(haps.san)
  
  ##adding sequencing errors 
  err.haps <- SeqErr(haps, selsite)
  haps <- err.haps
  rm (err.haps)
  
  ##keep only derived haplotypes
  derived_haplotypes <- haps[which(haps[,selsite]=="1"),]
  ##shorten derived haplotypes to only variants 
  temp <- SegsiteOut(derived_haplotypes, selsite)
  haps.der.var <- temp[[3]]
  s <- temp[[1]] # number of segregating sites on derived haplotypes
  rm(temp)
  pi <- CalcPi(haps.der.var) # pi on derived haplotyeps
  
  #*********converting to rehh input************
  
  system(paste("head -n 8 ", TrajFile, ".mssel | tail -n 1 > ", TrajFile, ".position", sep="")) # store positions
  ## creating rehh hap file
  write.table(haps,paste(TrajFile,".rehh.hap",sep=""),col.names=F,row.names = T,quote=F,sep=' ')
  ## creating rehh map file
  position<-read.table(paste(TrajFile,".position", sep=""))
  #fake.map<- data.frame(SNP=c(1:nchar(temp.hap[1,1])),CHR=15,POS=round(31599*as.numeric(position[1,-1])) + 48403201, anc=0,der=1)
  fake.map<- data.frame(SNP=c(1:length(haps[1,])),CHR=15,POS=round(31701*as.numeric(position[1,-1])) + 48403169, anc=0,der=1)
  fake.map$POS<-as.numeric(fake.map$POS)
  fake.map$POS[which(duplicated(fake.map$POS,fromLast = T))]<-fake.map$POS[which(duplicated(fake.map$POS,fromLast = T))]-1
  write.table(fake.map, paste(TrajFile, ".rehh.map", sep=""), col.names = F,row.names = F,quote=F,sep=' ')
  ## do ihh calc
  hap <- data2haplohh(hap_file=paste(TrajFile, ".rehh.hap", sep=""), map_file = paste(TrajFile, ".rehh.map", sep=""),recode.allele = T)
  hap.ehh <-calc_ehh(hap, mrk = selsite, plotehh = F)
  ihh <- as.numeric(hap.ehh$ihh[2])
  ehh.up<- as.numeric(hap.ehh$ehh[2,1])
  ehh.down <- as.numeric(hap.ehh$ehh[2,dim(hap.ehh$ehh)[2]])
  ## extract s2 from traj name
  temp <- strsplit(as.character(TrajFile), split="_")[[1]]
  sel <- temp[2]
  
  #clean up
  name.random <- round(500000*runif(1))
  system(paste("mv ",TrajFile,".mssel ./mssel_output/",TrajFile,"_id",name.random, ".mssel", sep="")) 
  system(paste("rm ",TrajFile,".*", sep=""))
  return(c(sel, pi,s,ihh, ehh.up, ehh.down, r.rate))
}



###################
#  main function  #
###################
args <- commandArgs(T)
traj.list.init<-read.table("traj.list")
#seg.size <- nrow(traj.list.init) %/% Njob + 1 
seg.size <- 50
traj.list <- data.frame(V1=traj.list.init[c((seg.size*(as.numeric(args[1])-1)+1):min(seg.size*(as.numeric(args[1])), nrow(traj.list.init))),])
sum.stats <- matrix(nrow=10*nrow(traj.list),ncol=7)
colnames(sum.stats) <- c("s2", "pi", "ss","ihh","ehh.up","ehh.down","RecRate")

for(i in 1:(10*nrow(traj.list))){
  TrajFile <- as.character(traj.list$V1[(i-1)%/%10+1])
  sum.stats[i,] <- SumStat(TrajFile)
}

write.table(sum.stats, paste("job", as.character(args[1]), "_sumstats.out", sep=""),col.names = T,row.names = F, quote=F,sep='\t')

