#include "ChangePointModelJointNormal.h"
#include <cmath> 

ChangePointModelJointNormal::ChangePointModelJointNormal() {
	m_startup=20;
}


ChangePointModelJointNormal::ChangePointModelJointNormal(const std::vector<double> &thresholds, int startup) {
	m_thresholds = thresholds;
	m_startup=startup;
}
	
//nb is 0, only higher if we are using window
void ChangePointModelJointNormal::cpmMLEaux(std::vector<double> &Us) {
	int i,n1;
	double mu1,mu2, Sok, Son, Skn,C,G,n2;

	Us.reserve(m_statistics[0].size());	
	
	int sz = m_statistics[0].size();
    Us.push_back(0);	
	for (i = 1 ; i < sz-2 ; i++) {
		n1=i+1; n2 = n-n1;
		mu1 = m_statistics[0][n1-1]/(double)n1;  
        mu2 = (m_statistics[0][n-1]-m_statistics[0][n1-1])/(n2);
		
		Sok = m_statistics[1][n1-1]/(double)n1;
		Son = m_statistics[1][n-1]/(double)n;
		Skn = m_statistics[1][n-1] - m_statistics[1][n1-1] - (n1*(n-n1)*(mu1-mu2)*(mu1-mu2))/n;
		Skn = Skn/(double)n2;

		C = 1 + 11/(double)12 * (1/(double)n1 + 1/(double)n2 - 1/(double)n) + (1/(double)(n1*n1) + 1/(double)(n2*n2) - 1/(double)(n*n));
		
		G = (n1*log(Son/Sok) + n2*log(Son/Skn))/C;
		Us.push_back(G);	
	}
	Us.push_back(0); Us.push_back(0);
}

void ChangePointModelJointNormal::reset() {
	m_statistics[0].clear();
	m_statistics[1].clear();
	n=0;
}

	
