#include "ChangePointModelF.h"
#include <cmath> 

ChangePointModelF::ChangePointModelF() {
	m_startup=20;
}


ChangePointModelF::ChangePointModelF(const std::vector<double> &thresholds, int startup) {
	m_thresholds = thresholds;
	m_startup=startup;
}
	
//nb is 0, only higher if we are using window
void ChangePointModelF::cpmMLEaux(std::vector<double> &Us) {
	int i,n1;
	double mu1,mu2,V1,V2,sigma1,sigma2,sigma,C,G,n2;

	Us.reserve(m_statistics[0].size());	
	
	int sz = m_statistics[0].size();
    Us.push_back(0);	
	for (i = 1 ; i < sz-2 ; i++) {
		n1=i+1; n2 = n-n1;
		mu1 = m_statistics[0][n1-1]/(double)n1;  
        mu2 = (m_statistics[0][n-1]-m_statistics[0][n1-1])/(n2);
		
		V1 = m_statistics[1][n1-1];
		V2 = m_statistics[1][n-1] - m_statistics[1][n1-1] - (n1*(n-n1)*(mu1-mu2)*(mu1-mu2))/n;
		
		sigma1 = V1/(n1-1);//cast to double
		sigma2 = V2/(n-n1-1);
		sigma = (V1+V2)/(n-2);	
		C = 1 + ( 1/(double)(n1-1) + 1/(double)(n-n1-1) - 1/(double)(n-2))/(double)3;
		G = ((n1-1)*log(sigma/sigma1) + (n-n1-1)*log(sigma/sigma2))/C;
		Us.push_back(G);	
	}
	Us.push_back(0); Us.push_back(0);
}

void ChangePointModelF::reset() {
	m_statistics[0].clear();
	m_statistics[1].clear();
	n=0;
}



	
