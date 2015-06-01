
#ifndef CHANGE_POINT_MODEL_T
#define CHANGE_POINT_MODEL_T

#include "ChangePointModel.h"
#include <vector>
using namespace std;

class ChangePointModelT : public ChangePointModel {
	public:
		ChangePointModelT();	
		ChangePointModelT(const std::vector<double> &thresholds, int startup=20);	
		void cpmMLEaux(std::vector <double> &Us);
		void updateStatistics(const double &obs);
		void reset();
	private:
};

#endif

