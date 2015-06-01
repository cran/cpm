
#ifndef CHANGE_POINT_MODEL_MOOD
#define CHANGE_POINT_MODEL_MOOD

#include "ChangePointModelWithRanks.h"
#include <vector>
using namespace std;

class ChangePointModelMood : public ChangePointModelWithRanks {
	public:
		ChangePointModelMood();	
		ChangePointModelMood(const std::vector<double> &thresholds, int startup=20);	
		void cpmMLEaux(std::vector <double> &Us);
	private:
};

#endif

