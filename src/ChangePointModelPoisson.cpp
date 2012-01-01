#include "ChangePointModelPoisson.h"
#include <cmath> 

ChangePointModelPoisson::ChangePointModelPoisson() {
	m_startup=20;
}


ChangePointModelPoisson::ChangePointModelPoisson(const std::vector<double> &thresholds, int startup) {
	m_thresholds = thresholds;
	m_startup=startup;
}

//poisson GLR is based only on sum of observations in each sample
void ChangePointModelPoisson::updateStatistics(const double &obs) {
	double S = obs; 

	if (m_statistics[0].size() > 0) {
		S = m_statistics[0].back() + obs;
	}
	m_statistics[0].push_back(S);
}
	

void ChangePointModelPoisson::cpmMLEaux(std::vector <double> &Us) {
	int i,sz;
	double n0,n1,s0,s1,d, temp;
	Us.reserve(m_statistics[0].size());	
	sz = m_statistics[0].size();
    Us.push_back(0);
    
	for (i = 1 ; i < sz-2 ; i++) {
		n0 = (double) i+1;
		n1 = (double) sz-n0;
		s0 = (double) m_statistics[0][i];
		s1 = (double) m_statistics[0].back() - s0;
		d = n1/n0;
		
		if (s0==0) {s0=0.5;}
		if (s1==0) {s1=0.5;}
		
		temp =  fabs((log(s1/s0)  - log(d))/sqrt(1/s0 + 1/s1));
		Us.push_back(temp);    
	}
    Us.push_back(0);
    Us.push_back(0);
}

void ChangePointModelPoisson::reset() {
	m_statistics[0].clear();
	n=0;
}



	
