setwd("/home/mlin/data/Coalescent/EuroMig_only/Khm/trajectories/")

args=commandArgs(TRUE)
job = as.numeric(args[1]) # potentially have 1,000 jobs, each job corresponding to 1 s2

s2 = round(runif(1, min=0.005, max=1), 2)

for (i in 1:100 ){ # each s2 has 100 trajectories
  random = runif(1, min =1, max=1e5)%/%1
  system(paste("timeout 300 /home/mlin/progs/Meng_mssel_programs/traj2pop_wBantu 1 0.16 ", 2*s2, " 0.5 10000 20000 446 7 0.12 538 175 ", random, " | /home/mlin/progs/Meng_mssel_programs/stepftn3pop > s2_", s2, "_i",i,"_R",random,"_task",job,".traj", sep=""))
  system(paste("head -n 5 s2_", s2, "_i",i,"_R",random,"_task",job,".traj | tail -n 1 | cut -f 3 > temp_s2_", s2, "_i", i, "_R",random, "_task", job,sep=""))
  try(freq <- read.table(paste("temp_s2_", s2, "_i",i,"_R",random, "_task", job, sep=""), colClasses = "numeric"), silent = T)
  if(exists("freq")) {
    if(freq!= 0.325){# Situation - freq is not within the range, attention: check bin file
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
