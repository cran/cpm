
#ifndef CHANGE_POINT_MODEL_KS
#define CHANGE_POINT_MODEL_KS

#include "ChangePointModelWithOrders.h"
#include <vector>
using namespace std;

class ChangePointModelKS : public ChangePointModelWithOrders {
	public:
		ChangePointModelKS();	
		ChangePointModelKS(const std::vector<double> &thresholds, int startup=20);	
		void cpmMLEaux(std::vector <double> &Us);
	private:
};

#endif

