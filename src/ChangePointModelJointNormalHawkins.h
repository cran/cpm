
#ifndef CHANGE_POINT_MODEL_JOINT_NORMAL_HAWKINS
#define CHANGE_POINT_MODEL_JOINT_NORMAL_HAWKINS

#include "ChangePointModel.h"
#include "ChangePointModelT.h"
#include <vector>
using namespace std;

class ChangePointModelJointNormalHawkins : public ChangePointModelT {
	
    public:
		ChangePointModelJointNormalHawkins();	
		ChangePointModelJointNormalHawkins(const std::vector<double> &thresholds, int startup=2);	
		void cpmMLEaux(std::vector <double> &Us);
		void reset();
        double meanAdjustments[3];
        double sdAdjustments[3];

    
	private:
};

#endif

