#include "ChangePointModelMood.h"
#include <cmath>
#include <cstdlib>

ChangePointModelMood::ChangePointModelMood() {
	m_startup=20;
}


ChangePointModelMood::ChangePointModelMood(const std::vector<double> &thresholds, int startup) {
	m_thresholds = thresholds;
	m_startup = startup;
}


void ChangePointModelMood::cpmMLEaux(std::vector <double> &Us) {
	int i;
    double n0,n1,M,mu,sd,med;
	
	double sz = m_statistics[0].size();
	med = (sz+1)/2;
        
	double *cumsums;
	cumsums = (double *) malloc(sz * sizeof(double));
	cumsums[0] =  (m_statistics[1][0]-med)*(m_statistics[1][0]-med);
	
    Us.push_back(0);
	for (i = 1 ; i < sz; i++) {
		cumsums[i] = cumsums[i-1] + (m_statistics[1][i]-med)*(m_statistics[1][i]-med);
	}
	
	for (i = 1; i < sz-2 ;i++) {
		n0 = i+1;
		n1 = sz-n0;
        M = cumsums[i];
        mu= n0*(sz*sz - 1)/12;
        sd = sqrt(n0*n1*(sz+1)*(sz*sz-4)/180);
		Us.push_back(fabs((M-mu)/sd));	
	}
    Us.push_back(0);
    Us.push_back(0);
	free(cumsums);
}

