using namespace std;

#ifndef CHANGE_POINT_MODEL_JOINT_NORMAL
#define CHANGE_POINT_MODEL_JOINT_NORMAL

#include "ChangePointModel.h"
#include "ChangePointModelT.h"
#include <vector>

class ChangePointModelJointNormal : public ChangePointModelT {
	public:
		ChangePointModelJointNormal();	
		ChangePointModelJointNormal(const std::vector<double> &thresholds, int startup=20);	
		void cpmMLEaux(std::vector <double> &Us);
		void reset();
	private:
};

#endif

