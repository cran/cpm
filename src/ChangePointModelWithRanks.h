
#ifndef CHANGE_POINT_MODEL_WITH_RANKS
#define CHANGE_POINT_MODEL_WITH_RANKS

#include "ChangePointModel.h"
#include <vector>
using namespace std;

class ChangePointModelWithRanks : public ChangePointModel {
	public:
		ChangePointModelWithRanks() {}
        	void updateStatistics(const double &obs);
        	void reset(); 
	private:
};

#endif

