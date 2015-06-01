#include "ChangePointModel.h"
#include <algorithm>
#include <vector>
#include <cmath>
#include <R.h>

ChangePointModel::ChangePointModel() {
	std::vector<double> vec;
	vec.reserve(5000);		
	m_statistics.push_back(vec);
	m_statistics.push_back(vec);
	n=0;
}

ChangePointModel::~ChangePointModel() {}

void ChangePointModel::processPoint(const double &obs) {
	n++;		
	updateStatistics(obs);	//add to window
}	 


void ChangePointModel::processStream(const std::vector<double> &x, std::vector<double> &Us, std::vector<int> &cps, std::vector<int> &dts) {
	double U;
  int maxindex, dt, cp, i;
	int lastChange = 0;
	int sz = x.size();
	
	double threshold = 0;
	int thresholdLength = m_thresholds.size();
    i = -1;
    while (i < sz-1) {
    	++i;
        //Rprintf("i:%d len:%d\n",i,(int)x.size());
        this->processPoint(x[i]);

		  if (n >= m_startup) {
  			this->cpmMLE(U,maxindex);
  			//Us.push_back(U);
  
  			if (thresholdLength==0) {threshold=99999.0;}
  			else if (n >= thresholdLength) {threshold = m_thresholds[thresholdLength-1];}
  			else {threshold = m_thresholds[n-1];}
  
  			//change found, we must reset
  			if (U > threshold) {
  				dt = i + 1;
  				cp = lastChange + maxindex + 1;
  				lastChange = cp;
  				i = cp-1; //start from observation after change
  				dts.push_back(dt);
          cps.push_back(cp);
  				this->reset();			
        }
		} else {
			//Us.push_back(0);
		}
	}
}	

void ChangePointModel::detectChange(const std::vector<double> &x, std::vector<double> &Us, int &cp, int &dt) {
	double U;
  int maxind;

	int sz = x.size();

	double threshold = 0;
	int thresholdLength = m_thresholds.size();
	for (int i = 0; i < sz; i++) {
		this->processPoint(x[i]);
		if (n >= m_startup) {
			this->cpmMLE(U, maxind);
			Us.push_back(U);

			if (thresholdLength==0) {threshold=9999999.0;}
			else if (n >= thresholdLength) {threshold = m_thresholds[thresholdLength-1];}
			else {threshold = m_thresholds[n-1];}

        
			if (U > threshold) {
				dt = i+1;
                cp = maxind+1;
				return;			
			}
		} else {
			Us.push_back(0);
		}
	}
	cp=0;
}	

//returns the maximum value, and its index
void ChangePointModel::cpmMLE(double &maxvalue, int &maxindex) {
	if (n < m_startup) {return;}
	
	std::vector<double> Us;
	Us.reserve(n);
	cpmMLEaux(Us);

	//Us[0]=0; Us[Us.size()-1]=0; Us[Us.size()-2]=0;
    
    maxvalue = 0;
    maxindex = 0;

	int sz = Us.size();
    for (int i = 1; i < sz-2; i++) {
        if (Us[i] > maxvalue) {
            maxvalue = Us[i];
            maxindex = i;
        }
    }
}

