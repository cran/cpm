
#ifndef CHANGE_POINT_MODEL_POISSON
#define CHANGE_POINT_MODEL_POISSON

#include "ChangePointModel.h"
#include <vector>
using namespace std;

class ChangePointModelPoisson : public ChangePointModel {
	public:
		ChangePointModelPoisson();	
		ChangePointModelPoisson(const std::vector<double> &thresholds, int startup=20);	
		void cpmMLEaux(std::vector <double> &Us);
		void updateStatistics(const double &obs);
		void reset();
	private:
};

#endif

