#include "ChangePointModelT.h"
#include <cmath> 
#include <R.h>
ChangePointModelT::ChangePointModelT() {
	m_startup=20;
}


ChangePointModelT::ChangePointModelT(const std::vector<double> &thresholds, int startup) {
	m_thresholds = thresholds;
	m_startup=startup;
}


void ChangePointModelT::updateStatistics(const double &obs) {
	double S = obs;
	double W = 0;
	if (m_statistics[0].size() > 0) {
		S = obs + m_statistics[0].back(); //last element of the S list
		double temp =  (n-1)*obs - m_statistics[0].back();
		W = m_statistics[1].back() + temp*temp / (n*(n-1));	
	}
	m_statistics[0].push_back(S);
	m_statistics[1].push_back(W);
}
	

void ChangePointModelT::cpmMLEaux(std::vector <double> &Us) {
	int i;
    double j,nn = (double) n;
	double sigma = n-2; sigma=sqrt(sigma/(sigma-2));
	double temp,E;	
	
	Us.reserve(m_statistics[0].size());	

	int sz = m_statistics[0].size();	
    Us.push_back(0);
    
	for (i = 1 ; i < sz-2 ; i++) {
		j=(double)i+1;
		temp = n*m_statistics[0][i] - j*m_statistics[0].back();
		E = temp*temp / (nn*j*(nn-j));
		Us.push_back(sqrt( (nn-2)*E / (m_statistics[1].back()-E))/sigma);
	}
    Us.push_back(0);
    Us.push_back(0);
}

void ChangePointModelT::reset() {
	m_statistics[0].clear();
	m_statistics[1].clear();
	n=0;
}



	
