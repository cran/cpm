#include "ChangePointModelJointNormalHawkins.h"
#include <cmath> 

ChangePointModelJointNormalHawkins::ChangePointModelJointNormalHawkins() {
	m_startup=20;  
	meanAdjustments[0] = 2.2989; meanAdjustments[1] = 2.0814; meanAdjustments[2] = 2.0335;
    sdAdjustments[0] = 2.3151; sdAdjustments[1] = 2.0871; sdAdjustments[2] = 2.0368;
}


ChangePointModelJointNormalHawkins::ChangePointModelJointNormalHawkins(const std::vector<double> &thresholds, int startup) {
	m_thresholds = thresholds;
	m_startup=startup;

    meanAdjustments[0] = 2.2989; meanAdjustments[1] = 2.0814; meanAdjustments[2] = 2.0335;
    sdAdjustments[0] = 2.3151; sdAdjustments[1] = 2.0871; sdAdjustments[2] = 2.0368;
}
	
//nb is 0, only higher if we are using window
void ChangePointModelJointNormalHawkins::cpmMLEaux(std::vector<double> &Us) {
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

    //now do the adjustment
    int len = Us.size(); 
    if (len < 10) {return;}

    //needs to be n-3, n-4, n-5 since C indexes from 0, so if length=20 last element is only index 19
    Us[1] = (Us[1] - meanAdjustments[0])/sdAdjustments[0];
    Us[len-3] = (Us[len-3] - meanAdjustments[0])/sdAdjustments[0];
    Us[2] = (Us[2] - meanAdjustments[1])/sdAdjustments[1];
    Us[len-4] = (Us[len-4] - meanAdjustments[1])/sdAdjustments[1];
    Us[3] = (Us[3] - meanAdjustments[2])/sdAdjustments[2];
    Us[len-5] = (Us[len-5] - meanAdjustments[2])/sdAdjustments[2];

    //now convert back to chi-square
    double truemean = 2;
    double truesd = 2;

    Us[1] = (Us[1] * truesd) + truemean;
    Us[2] = (Us[2] * truesd) + truemean;
    Us[3] = (Us[3] * truesd) + truemean;
    Us[n-3] = (Us[n-3] * truesd) + truemean;
    Us[n-4] = (Us[n-4] * truesd) + truemean;
    Us[n-5] = (Us[n-5] * truesd) + truemean;


}




void ChangePointModelJointNormalHawkins::reset() {
	m_statistics[0].clear();
	m_statistics[1].clear();
	n=0;
}

	
