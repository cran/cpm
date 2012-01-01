#include "ChangePointModelMW.h"
#include <cmath>
#include <cstdlib>

ChangePointModelMW::ChangePointModelMW() {
	m_startup=20;
}


ChangePointModelMW::ChangePointModelMW(const std::vector<double> &thresholds, int startup) {
	m_thresholds = thresholds;
	m_startup = startup;
}


void ChangePointModelMW::cpmMLEaux(std::vector <double> &Us) {
	int i;
	double n0,n1,R1,U,mu,sd;
	
	double sz = m_statistics[0].size();
	double *cumsums;
	cumsums = (double *) malloc(sz * sizeof(double));
	cumsums[0] =  m_statistics[1][0];
	
	for (i = 1 ; i < sz; i++) {
		cumsums[i] = cumsums[i-1] + m_statistics[1][i];
	}
	
    Us.push_back(0);
	for (i = 1; i < sz-2 ;i++) {
		n0 = i+1;
		n1 = sz-n0;
		R1 = cumsums[i];
		
		U = R1 - ( n0*(n0+1)/2) ;		
				
		mu = n0*n1/2;			
		sd = sqrt(n0*n1*(n0+n1+1)/12);
		Us.push_back(fabs((U-mu)/sd));		
	}
    Us.push_back(0);
    Us.push_back(0);
	free(cumsums);
}

