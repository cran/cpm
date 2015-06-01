#include "ChangePointModelExponential.h"
#include <cmath> 

ChangePointModelExponential::ChangePointModelExponential() {
	m_startup=20;
}


ChangePointModelExponential::ChangePointModelExponential(const std::vector<double> &thresholds, int startup) {
	m_thresholds = thresholds;
	m_startup=startup;
}

//poisson GLR is based only on sum of observations in each sample
void ChangePointModelExponential::updateStatistics(const double &obs) {
	double S = obs; 

	if (m_statistics[0].size() > 0) {
		S = m_statistics[0].back() + obs;
	}
	m_statistics[0].push_back(S);
}
	

void ChangePointModelExponential::cpmMLEaux(std::vector <double> &Us) {
	int i,sz;
	double n0,n1,s0,s1,K,temp, dsz;
	Us.reserve(m_statistics[0].size());	
	sz = m_statistics[0].size();
  dsz = (double) sz;
    //Us.push_back(0);
    
	//for (i = 1 ; i < sz-2 ; i++) {
	for (i = 0 ; i < sz-1 ; i++) {
		n0 = (double) i+1;
		n1 = (double) sz-n0;
		s0 = (double) m_statistics[0][i];
		s1 = (double) m_statistics[0].back() - s0;
		K = sz*log(dsz) - n0*log(n0)-n1*log(n1);
		temp = -sz*log(s0+s1) + n0*log(s0) + n1*log(s1)+K;
		Us.push_back(-2*temp);    
	}
    Us.push_back(0);
    //Us.push_back(0);
}

void ChangePointModelExponential::reset() {
	m_statistics[0].clear();
	n=0;
}



	
