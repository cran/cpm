
#ifndef CHANGE_POINT_MODEL_F
#define CHANGE_POINT_MODEL_F

#include "ChangePointModel.h"
#include "ChangePointModelT.h"
#include <vector>
using namespace std;

class ChangePointModelF : public ChangePointModelT {
	public:
		ChangePointModelF();	
		ChangePointModelF(const std::vector<double> &thresholds, int startup=20);	
		void cpmMLEaux(std::vector <double> &Us);
		void reset();
	private:
};

#endif

