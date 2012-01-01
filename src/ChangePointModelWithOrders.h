using namespace std;

#ifndef CHANGE_POINT_MODEL_WITH_ORDERS
#define CHANGE_POINT_MODEL_WITH_ORDERS

#include "ChangePointModel.h"
#include <vector>

class ChangePointModelWithOrders : public ChangePointModel {
	public:
		ChangePointModelWithOrders() {}
        	void updateStatistics(const double &obs);
        	void reset(); 
	private:
};

#endif

