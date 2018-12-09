setwd("/home/mlin/data/Coalescent/4pop_complex/nama/trajectories/")

args=commandArgs(TRUE)
job = as.numeric(args[1]) # potentially have 1,000 jobs, each job corresponding to 10 s2

s2_range = round(runif(10, min=0.005, max=0.2), 2)
mEP = 0.3 # Eurasian to eastern pastoralists
mBK = 0.02 # Bantu to KhoeSan
nsam = 322
nder = 172

for (s2 in s2_range){
for (i in 1:100 ){ # each s2 has 100 trajectories
  #sampling demographic parameters 
  mPK = round(runif(1, min=0.02, max=0.1),2)
  tEP = round(runif(1, min=167, max=333),0)
  tPK = round(runif(1, min=30, max=100),0)

  random = runif(1, min =1, max=1e5)%/%1
  #traj nreps h nsamK nderK s t mEP tEP mEK tEK mPK tPK nsamP nderP mBK tBK seed 
  system(paste("timeout 300 /home/mlin/progs/Meng_mssel_programs/traj4pop_nama 1 0.5 ", nsam, " ", nder, " ", 2*s2, " 446 ", mEP, " ", tEP, " 0.17 10 ", mPK, " " ,tPK, " 100 38 0.02 14 ", random, " | /home/mlin/progs/Meng_mssel_programs/stepftn4pop > s2_", s2, "_i",i,"_R",random,"_task",job,".traj", sep=""))
  system(paste("head -n 5 s2_", s2, "_i",i,"_R",random,"_task",job,".traj | tail -n 1 | cut -f 2 > temp_s2_", s2, "_i", i, "_R",random, "_task", job,sep=""))
  try(freq <- read.table(paste("temp_s2_", s2, "_i",i,"_R",random, "_task", job, sep=""), colClasses = "numeric"), silent = T)
  if(exists("freq")) {
    if(freq!= 0.525){# Situation - freq is not within the range, attention: check bin file
      system(paste("rm temp_s2_", s2, "_i",i,"_R",random, "_task", job, sep=""))
      system(paste("rm s2_", s2, "_i",i,"_R",random, "_task", job, ".traj", sep=""))
    }else{ #Situation - freq is what we want, so we keep .traj
      system(paste("rm temp_s2_", s2, "_i",i,"_R",random, "_task", job, sep=""))
    }
    rm(freq)
  }else{# situation - timeout, so .traj is empty
      system(paste("rm temp_s2_", s2, "_i",i,"_R",random, "_task", job, sep=""))
      system(paste("rm s2_", s2, "_i",i,"_R",random, "_task", job, ".traj", sep=""))
    }
}
}
