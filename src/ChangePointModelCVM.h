#include <vector>
using namespace std;

#ifndef CHANGE_POINT_MODEL_CVM
#define CHANGE_POINT_MODEL_CVM

#include "ChangePointModelWithOrders.h"


class ChangePointModelCVM : public ChangePointModelWithOrders {
	public:
		ChangePointModelCVM();	
		ChangePointModelCVM(const std::vector<double> &thresholds, int startup=20);	
		void cpmMLEaux(std::vector <double> &Us);
	private:
};

#endif


