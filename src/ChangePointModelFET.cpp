#include "ChangePointModelFET.h"
#include <R.h>
#include <Rmath.h>

ChangePointModelFET::ChangePointModelFET() {
	m_startup=20;
	m_lambda=1.0;
}


ChangePointModelFET::ChangePointModelFET(const std::vector<double> &thresholds, int startup, double lambda) {
	m_thresholds = thresholds;
	m_startup=startup;
	m_lambda=lambda;
}


void ChangePointModelFET::updateStatistics(const double &obs) {
	double S = obs; 

	if (m_statistics[0].size() > 0) {
		S = m_statistics[0].back() + obs;
	}
	m_statistics[0].push_back(S);
}
	
void ChangePointModelFET::cpmMLEaux(std::vector <double> &Us) {
	int i,sz,n0,n1,s0,s1;
	
	Us.reserve(m_statistics[0].size());	
	sz = m_statistics[0].size();
    Us.push_back(0);
	for (i = 1 ; i < sz-2 ; i++) {
		n0 = i+1;
		n1 = sz-n0;
		s0 = m_statistics[0][i];
		s1 = m_statistics[0].back() - s0;
		//Us.push_back(1.0 - my_gsl_cdf_hypergeometric_P(s0,s0+s1,n0+n1-s0-s1,n0));	
		Us.push_back(1.0 - phyper(s0,s0+s1,n0+n1-s0-s1,n0,1,0));	        
	}
    Us.push_back(0);
    Us.push_back(0);
    
	if (sz > 3 && m_lambda < 1.0) {
		for (i = 2 ; i < sz-2; i++) { 
			Us[i] = (1 - m_lambda)*Us[i-1] + m_lambda * Us[i];
		}
	}
}
		
		
		
void ChangePointModelFET::reset() {
	m_statistics[0].clear();
	m_statistics[1].clear();
	n=0;
}



	
