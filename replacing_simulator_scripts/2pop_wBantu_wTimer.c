// Modification: this simulates a scenario where migration happens between pop1 and pop2,
// where beneficial allele has been selected to fixation in pop1 when migration occurs.
// Before migration, the locus has never been segregating in pop2.
// nreps: number of trajectories to generate
// s1: selection coefficient in the first population (source)
// s2: selection coefficient in the second population (target)
// h: dominance coefficient
// N1: diploid population size 1 (Euro)
// N2: diploid population size 2 (KS)
// t1: an onset of selection in pop1 (generations ago)
// t2: an onset of selection in / migration to pop2 (generations ago)
// m: migration rate (one pulse)
// nsam: total number of chromosomes in your samples (pop2)
// nder: total number of chromosomes carrying a focal allele (pop2)
// seed: for the random number generator

//Compilation:  [NO: gcc -o trajc trajconst.c binomial.c rand1.c -lm]
// gcc -o trajc trajpopsel_wBantu.c binomial_mt.c mt.c
// usage:
// trajc nrep s1 s2 h N1 N2 t1 t2 m nsam nder seed > my.out


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

main(int argc, char **argv )
{
  double s1, s2, h, m,  pfinal, wbar1, wbar2, pprime1, pprime2, prob1, prob2, *freq1, *freq2, p1, p2;
  int N1, N2, nreps, gen, t1, t2, maxgen, rep, i, nsam, nder, seed, count;
  float bnldev(float , int);
  int  success1, success2  ;
 
 //record timing of running
 double time_spent;
 //float begin, end;
 
  //double ran1() ;
  double genrand_real3(void)  ;
  float bnldev(float,int) ;
  unsigned short seedv[3], *seed48()  ;


clock_t begin=clock();



  if( argc < 11 ){
    fprintf(stderr,"trajc nrep s1 s2 h N1 N2 t1 t2 m nsam nder seed \n");
    exit(1);
  }
  
  printf("// ");
  for( i=0; i<argc ; i++) printf(" %s",argv[i]);
  printf("\n");

  
  maxgen = 5000 ;
  freq1 = (double *)malloc( (unsigned)maxgen*sizeof(double) ) ; // freq1 and freq2 are pointers to address storing frequency
  freq2 = (double *)malloc( (unsigned)maxgen*sizeof(double) ) ;
  nreps = strtol(argv[1],NULL,10); // str to long, base10
  s1 = strtod( argv[2],NULL ) ; // str to double
  s2 = strtod(argv[3],NULL);
  h = strtod( argv[4],NULL ) ;
  N1 = strtol( argv[5],NULL,10);
  N2 = strtol(argv[6], NULL, 10); 
  t1 = strtol(argv[7], NULL, 10);
  t2 = strtol(argv[8], NULL, 10);
  m = strtod(argv[9], NULL);
  //if( m >= 1 ) m = m/(2.0*N2) ;
  nsam = strtol(argv[10],NULL,10);
  nder = strtol(argv[11],NULL,10);
  pfinal = nder/(double)nsam;


  fprintf(stderr, "Parameters input:\n") ;
  fprintf(stderr,"%d %lf %lf %lf %d %d %d %d %lf %d %d\n",nreps, s1, s2, h, N1, N2, t1, t2, m, nsam, nder);
  //#if( argc > 12 ) seedv[0] = strtol( argv[12] ,NULL,10) ;
  //#else seedv[0] = 3579 ;
  //#seedv[1] = 27011;
  //#seedv[2] = 59243 ;
  //#seed48( seedv);
   sscanf(argv[12], "%x", &seed) ;
 // Initialize Mersenne Twistor.
        void init_by_array(unsigned long init_key[], int key_length)    ;
        unsigned long init[4] = {seed, 0x265, 0x367, 0x569}, length=4 ;
        init_by_array(init, length) ;



  //#printf("%d\n",nreps);
  printf("%d N0: %d\n", nreps, N2);
  for( rep =0; rep < nreps; rep++){ //within one trajectory
      // Population 1
      success1 = 0 ;
      while(success1 == 0 ){
         p1 = 1/(2.0*N1) ; // Hard sweep model for population 1
         p2 = 0.0 ; // beneficial allele absent in population 2
         gen = 0 ;
         freq1[0] = 0.0 ;
         freq2[0] = 0.0 ;
	     gen = 1 ;
         freq1[gen] = p1 ;
         freq2[gen] = 0.0 ;
 
         // before migration happens 
         while( (gen < (t1-t2)) && (p1 > 0. ) && (p1 < 1.0 ) ){
            gen++;
		    if( gen > maxgen-1 ){
			    maxgen += 10000 ;
			    freq1 = (double *)realloc(freq1,(unsigned)(maxgen*sizeof(double)) ); // needs bigger space 
			    freq2 = (double *)realloc(freq2,(unsigned)(maxgen*sizeof(double)) );
		    }
            wbar1 = p1*p1*(1.+s1) + 2.*p1*(1.-p1)*(1.+s1*h)  + (1.-p1)*(1.-p1) ; // mean fitness
            pprime1 = ( p1*p1*(1.+s1) + p1*(1.-p1)*(1.+s1*h) ) / wbar1 ; // expected p at t+1: p(t)*w1 / wbar
			p1 = bnldev(pprime1,2*N1)/(2*N1) ; // binomial sampling in 2N1, with expected p(t+1)
	        freq1[gen] = p1 ;
	        freq2[gen] = 0.0 ;
	     }
	     // prob1 = pow(p1/0.99, (double)(nsam-1))*pow((1.0-p1)/(1.0-0.99), 1.0); 
	     //if(ran1() < prob1){ /* this situation does not include fixation */
	     //success1 =1 ; //other wise have to re-initiate this trajectory for pop1 from scratch
	     //} 
	     if(p1==1.0){/* another situation able to end while loop for pop1: fixation reached before migration happens */
	     success1 =1 ;
	     for (i = gen; i <= t1-t2; i++){ /* fill in freq traj between fixation time and migration time*/
            freq1[i+1] = p1;
            freq2[i+1] = 0.0;	     
	        }
	     } 
	  }
	
	//add!!!
	count = 0;  
	  // migration happens at t1-t2+1 generation forward in time
	  success2 = 0 ;
	  while(success2 == 0){
	     count++;
	     gen = t1-t2+1;
	     freq2[gen] = m*freq1[gen] ;
	     p2 = freq2[gen];
	     if( freq1[gen] < 1.0){
	     fprintf(stderr, "Population 1 didn't go to fixation. Please reconsider your t1 or s1.\n");
	     fprintf(stderr, "Your current allele frequency in population 1 is: %lf.\n", freq1[gen]);
	     exit(1);
	     } 
	     while( (gen < t1) && (p2 >0. ) && (p2 < 1.0) ){
	        gen++;
	        if( gen > maxgen-1 ){
			    maxgen += 10000 ;
			    freq1 = (double *)realloc(freq1,(unsigned)(maxgen*sizeof(double)) ); // needs bigger space 
			    freq2 = (double *)realloc(freq2,(unsigned)(maxgen*sizeof(double)) );
		    }
		    freq1[gen] = 1.0 ; // pop1 has gone to fixation before migration
	        wbar2 = p2*p2*(1.+s2) + 2.*p2*(1.-p2)*(1.+s2*h)  + (1.-p2)*(1.-p2) ;
	        pprime2 = ( p2*p2*(1.+s2) + p2*(1.-p2)*(1.+s2*h) ) / wbar2 ;
	        //p2 = bnldev(pprime2, 2*N2)/(2*N2) + m ; // each generation there's migration into pop2 at rate m
	        p2 = bnldev(pprime2, 2*N2)/(2*N2);
		freq2[gen] = p2 ;
	     }
         prob2 = pow(p2/pfinal, (double)nder)*pow((1.0-p2)/(1.0-pfinal), (double)(nsam-nder)); // in case p(t+1) is too small to cause while loop quits
         //#if(ran1() < prob2){ // ran1() is a random number generated between 0 and 1, compiled from rand1.c
         
         // if taking too much time... Exit
         if(count == 1e4){
         clock_t end = clock();
         time_spent = (double) (end-begin) / CLOCKS_PER_SEC;
         fprintf(stderr, "Number of attempts to generate a trajectory exceeds 10,000. Terminating program.\n");
         fprintf(stderr, "Total time(s): %lf.\n", time_spent);
         exit(1);
         }
         
         if(genrand_real3() < prob2){
            success2 = 1 ;
         }
      }

	  //#printf("#\n");
	  printf("# s:%lf age: %lf freq: %lf\n", s2, (double)t2 / (4.0 * N2), m);
      for( i = 0 ; i <=  gen; i++){
         //printf("%lf\t%lf\t%lf\n",i/(double)(4*(N1+N2)),freq1[i],freq2[i] );
	 printf("%lf\t%lf\t%lf\t%lf\n",i/(double)(4*N2),freq1[i],freq2[i], 0.0 ); // Ne 1.0 is set the same as the target pop
      }
   }

fprintf(stderr," !%.10f\n", nreps / (double)count) ;

clock_t end = clock();
time_spent = (double) (end-begin) / CLOCKS_PER_SEC;

fprintf(stderr, "Total time(s): %lf.\n", time_spent);


}
  
