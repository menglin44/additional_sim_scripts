# additional_sim_scripts

### 1) replacing_simulator_scripts/

These are scripts that replace the corresponding original files in Mssel software directory, or need to be put directly under Mssel software directory. 

popsize*.c: Ne files required for trajectory generator. Europeans (E), Khomani (K), Nama (Nama), Eastern Pastoralists (P).

2pop_wBantu_wTimer.c: Trajectory generator under European Source Model, with a timeout feature based on number of iterations it attempts to generate a traj matching the desired allele freq (maximum allowed 1M). It generates a trajectory for European, San, and Bantu, under the model that only a recent migration from Europeans to San happened. Compiling instrucstions are in the header.

Trajectory_Euro_EP_Bantu_wTimer.c: Trajectory generator under East African Source Model. Timer function same as above. generates a trajectory for San, Eastern Pastoralists, Europeans, and Bantu, under the East African Source scenario. Compiling instructions are in the header.

stepftn_3pop.c: replace stepftn.c for European Source Model

stepftn_4pop.c: replace stepftn.c for East African Source Model

freqints.h: uses different bins than original freqints.h in mssel.

### 2) generating_trajectories

Pipeline scripts to call trajectory generators and condition them on accepted final allele frequency (that matches the observed Khomani / Nama allele frequency). 

### 3) mssel_ABC

Pipeline scripts to run coalescent simulations after trajectories are generated. Followed by post simulation analyses (adding sequencing errors, calculating summary statistics etc.) 

