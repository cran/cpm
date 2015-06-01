
#ifndef CHANGE_POINT_MODEL_FET
#define CHANGE_POINT_MODEL_FET

#include "ChangePointModel.h"
#include <vector>
using namespace std;


class ChangePointModelFET : public ChangePointModel {
	public:
		ChangePointModelFET();	
		ChangePointModelFET(const std::vector<double> &thresholds, int startup=20, double lambda=1.0);	
		void cpmMLEaux(std::vector <double> &Us);
		void updateStatistics(const double &obs);
		void reset();
		double m_lambda;
	private:
};

#endif

