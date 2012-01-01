using namespace std;

#ifndef CHANGE_POINT_MODEL_WITH_RANKS
#define CHANGE_POINT_MODEL_WITH_RANKS

#include "ChangePointModel.h"
#include <vector>

class ChangePointModelWithRanks : public ChangePointModel {
	public:
		ChangePointModelWithRanks() {}
        	void updateStatistics(const double &obs);
        	void reset(); 
	private:
};

#endif

