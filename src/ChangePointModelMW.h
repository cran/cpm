
#ifndef CHANGE_POINT_MODEL_MW
#define CHANGE_POINT_MODEL_MW

#include "ChangePointModelWithRanks.h"
#include <vector>
using namespace std;

class ChangePointModelMW : public ChangePointModelWithRanks {
	public:
		ChangePointModelMW();	
		ChangePointModelMW(const std::vector<double> &thresholds, int startup=20);	
		void cpmMLEaux(std::vector <double> &Us);
	private:
};

#endif

