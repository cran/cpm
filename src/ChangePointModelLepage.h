
#ifndef CHANGE_POINT_MODEL_LEPAGE
#define CHANGE_POINT_MODEL_LEPAGE

#include "ChangePointModelMW.h"
#include "ChangePointModelMood.h"
#include <vector>
using namespace std;

class ChangePointModelLepage : public ChangePointModel {
	public:
		ChangePointModelMW m_cpmMW;
		ChangePointModelMood m_cpmMood;
		
		ChangePointModelLepage();	
		ChangePointModelLepage(const std::vector<double> &thresholds, int startup=20);	

		void updateStatistics(const double &obs);
		void cpmMLEaux(std::vector <double> &Us);
		void reset();
	private:
};

#endif

