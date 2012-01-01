#include "ChangePointModelKS.h"
#include <cmath>
#include <cstdlib>
#include <R.h>

ChangePointModelKS::ChangePointModelKS() {
	m_startup=20;
}


ChangePointModelKS::ChangePointModelKS(const std::vector<double> &thresholds, int startup) {
	m_thresholds = thresholds;
	m_startup = startup;
}


void ChangePointModelKS::cpmMLEaux(std::vector <double> &Us) {
	int i,j;
	double n0,n1,a,b,statistic, temp,z,correction;
	double N = (double) m_statistics[0].size();

	double *cumsums;
	cumsums = (double*) malloc(N * sizeof(double));
	
    Us.push_back(0);
	
	for (i = 1; i < N-2 ;i++) {
		n0 = i+1;
		n1 = N-i;
		
		a = 1/n0;
		b = -1/n1;

		
		for (j = 0 ; j < N ; j++) {
			if (m_statistics[1][j] <= n0) {
				cumsums[j] = a;
			}
			else {
				cumsums[j] = b;
			}			
		}

		for (j = 1 ; j < N; j++) {
			cumsums[j] = cumsums[j-1] + cumsums[j]; //i dont know how, but cumsums[j] must be difference between the two edfs for the j th point
		}
		statistic = 0;
		for (j = 0; j < N; j++) {
			temp = fabs(cumsums[j]);
			if(temp > statistic) {
				statistic=temp;
			}
		}
		
		//now we compute p value
		//this is a continuity correction for the KS statistic so I can more accurately use the asymptotic distribution. See Kim 1969
		correction = 0;
        
			//rearrange so n0 >= n1
			if (n1 > n0) {
				a = n1;
				n1 = n0;
				n0 = a;
			} 
		
			if (n0 > 2*n1) {
				correction=1/(2*sqrt(n0));
			}
			else {
				if ((int) n0 % (int) n1 == 0) {
					correction=2/(3*sqrt(n0));
				}
				else {
					correction=2/(5*sqrt(n0));
				}
			}
		//}
		z = statistic*sqrt((n0*n1)/(n0+n1)) + correction;

		z = z*z;
		Us.push_back(1-2*(exp(-2*z) - exp(-8 * z))); 
	}
    Us.push_back(0);
    Us.push_back(0);
    
	free(cumsums);
}
