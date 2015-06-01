

#ifndef CHANGE_POINT_MODEL_EXPONENTIAL
#define CHANGE_POINT_MODEL_EXPONENTIAL

#include "ChangePointModel.h"
#include <vector>
using namespace std;

class ChangePointModelExponential : public ChangePointModel {
	public:
		ChangePointModelExponential();	
		ChangePointModelExponential(const std::vector<double> &thresholds, int startup=20);	
		void cpmMLEaux(std::vector <double> &Us);
		void updateStatistics(const double &obs);
		void reset();
	private:
};

#endif

