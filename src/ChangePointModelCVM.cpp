#include "ChangePointModelCVM.h"
#include <cmath>
#include <cstdlib>
#include <R.h>

ChangePointModelCVM::ChangePointModelCVM() {
	m_startup=20;
}


ChangePointModelCVM::ChangePointModelCVM(const std::vector<double> &thresholds, int startup) {
	m_thresholds = thresholds;
	m_startup = startup;
}


void ChangePointModelCVM::cpmMLEaux(std::vector <double> &Us) {
	int i,j;
	double a,b,statistic,mu,sigma,N,prod,n0,n1;
	
	N = m_statistics[1].size();
	
	double *cumsums;
	cumsums = (double*) malloc(N * sizeof(double));


	mu = (double)1/6 + 1/(6*N);
    Us.push_back(0);
	for (i = 1; i < N-2 ;i++) {
		n0 = i+1;
		n1 = N-i-1;
		
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
			cumsums[j] = cumsums[j-1] + cumsums[j]; 
		}
		statistic = 0;
		for (j = 0; j < N; j++) {
			statistic+= cumsums[j]*cumsums[j];
		}
		prod = n0*n1;
		sigma = sqrt(  (double) 1/45 * (N+1)/(N*N)  * (4*prod*N - 3*(n1*n1 + n0*n0) - 2*prod)/(4*prod));

		Us.push_back((statistic * (prod)/(N*N) - mu)/sigma);

	}
	free(cumsums);
    Us.push_back(0); Us.push_back(0);
}

