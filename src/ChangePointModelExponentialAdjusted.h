
#ifndef CHANGE_POINT_MODEL_EXPONENTIAL_ADJUSTED
#define CHANGE_POINT_MODEL_EXPONENTIAL_ADJUSTED

#include "ChangePointModel.h"
#include "ChangePointModelExponential.h"
#include <vector>
using namespace std;


class ChangePointModelExponentialAdjusted : public ChangePointModelExponential {
	public:
		ChangePointModelExponentialAdjusted();	
		ChangePointModelExponentialAdjusted(const std::vector<double> &thresholds, int startup=20);	
		void cpmMLEaux(std::vector <double> &Us);
		void initarray();
		
        double meanAdjustments[5];
        double sdAdjustments[5];
		double digamma[10001];
		
	private:
};

#endif

