#include <stdio.h>
#include <math.h>

int naK = 20000;
int naBottleneck = 10000;
	
int popsizeK(int t) {
    int pop ;
    if (t > 20){
	pop = naK ;
	}else{
	pop = naBottleneck;
	}
    return(pop) ;
}

