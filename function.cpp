#include <math.h>
#include "function.h"

int getRelfromK(double ibd1, double ibd2, double bkg, double tot_genome, int maxDeg){
    double K = std::max((ibd1/4.0 + ibd2/2.0 - bkg/4.0)/tot_genome, 0.0); 
    if (K == 0){return -1;}
    int deg = (int)(-log2(K) + 0.5) - 1; // 四舍五入
    if (deg <= maxDeg){return deg;}
    else{return -1;}
}

