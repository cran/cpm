
#ifndef CHANGE_POINT_MODEL_JOINT_NORMAL_ADJUSTED
#define CHANGE_POINT_MODEL_JOINT_NORMAL_ADJUSTED

#include "ChangePointModel.h"
#include "ChangePointModelT.h"
#include <vector>
using namespace std;

class ChangePointModelJointNormalAdjusted : public ChangePointModelT {
	
    public:
		ChangePointModelJointNormalAdjusted();	
		ChangePointModelJointNormalAdjusted(const std::vector<double> &thresholds, int startup=2);	
		void cpmMLEaux(std::vector <double> &Us);
		void reset();
		void initarray();
		double digamma[10001]; //these are half integer values

        double meanAdjustments[3];
        double sdAdjustments[3];

    
	private:
};

#endif

