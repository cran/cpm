
#include <R.h>
#include <vector>
#include <algorithm>
#include <cstring>
#include "ChangePointModel.h"
#include "ChangePointModelWithRanks.h"
#include "ChangePointModelWithOrders.h"
#include "ChangePointModelT.h"
#include "ChangePointModelF.h"
#include "ChangePointModelMW.h"
#include "ChangePointModelMood.h"
#include "ChangePointModelFET.h"
#include "ChangePointModelLepage.h"
#include "ChangePointModelJointNormal.h"
#include "ChangePointModelJointNormalAdjusted.h"
#include "ChangePointModelJointNormalHawkins.h"
#include "ChangePointModelCVM.h"
#include "ChangePointModelKS.h"
#include "ChangePointModelPoisson.h"
#include "ChangePointModelExponential.h"
#include "ChangePointModelExponentialAdjusted.h"
using namespace std;


extern "C" {
	void cpmDetectChange(char **cpmType, double *x, int* nx, double *thresholds, int *nthresholds, int *startup, double *Ds, int *cp, int *dt, double *lambda) {
		std::vector<double> vthresholds(thresholds,thresholds+*nthresholds);
		
		std::vector<double> xv(x,x + *nx);
        std::vector<double> Dsv;
        Dsv.reserve(*nx);
        
        ChangePointModel* cpm;
		
        if (strcmp(*cpmType, "Student") ==0) {
			cpm = new ChangePointModelT(vthresholds,*startup);
		} else if (strcmp(*cpmType, "Bartlett") ==0) {
			cpm = new ChangePointModelF(vthresholds,*startup);
		} else if (strcmp(*cpmType, "MW") ==0) {
			cpm = new ChangePointModelMW(vthresholds,*startup);
		} else if (strcmp(*cpmType, "Mood") ==0) {
			cpm = new ChangePointModelMood(vthresholds,*startup);
		} else if (strcmp(*cpmType, "FET") ==0) {
			cpm = new ChangePointModelFET(vthresholds,*startup,*lambda);
		} else if (strcmp(*cpmType, "LP") ==0) {
			cpm = new ChangePointModelLepage(vthresholds,*startup);
		} else if (strcmp(*cpmType, "Joint") ==0) {
			cpm = new ChangePointModelJointNormal(vthresholds,*startup);
		} else if (strcmp(*cpmType, "JointAdjusted") ==0) {
    		cpm = new ChangePointModelJointNormalAdjusted(vthresholds,*startup);
		} else if (strcmp(*cpmType, "CVM") ==0) {
			cpm = new ChangePointModelCVM(vthresholds,*startup);
		} else if (strcmp(*cpmType, "KS") ==0) {
			cpm = new ChangePointModelKS(vthresholds,*startup);
		}  else if (strcmp(*cpmType, "Poisson") ==0) {
			cpm = new ChangePointModelPoisson(vthresholds,*startup);
		}  else if (strcmp(*cpmType, "Exponential") ==0) {
			cpm = new ChangePointModelExponential(vthresholds,*startup);
		}  else if (strcmp(*cpmType, "ExponentialAdjusted") ==0) {
			cpm = new ChangePointModelExponentialAdjusted(vthresholds,*startup);
		} else if (strcmp(*cpmType, "JointHawkins") ==0) {
    		cpm = new ChangePointModelJointNormalHawkins(vthresholds,*startup);
		} else {
			Rprintf("Change point model type not supported\n");
			return;
		}
		
        int cp2=0, dt2=0;
		cpm->detectChange(xv, Dsv, cp2, dt2);	
        std::copy(Dsv.begin(),Dsv.end(),Ds);
        *cp = cp2;
        *dt = dt2;
		delete cpm;
		return;
	}
	
	void cpmProcessStream(char **cpmType, double *x, int* nx, double *thresholds, int *nthresholds, int *startup, double *Ds, int *cps, int *dts, int *numChanges, double *lambda) {
		std::vector<double> vthresholds(thresholds,thresholds+*nthresholds);
		
		std::vector<double> xv(x,x + *nx);
        std::vector<double> Dsv;
        Dsv.reserve(*nx);
        
        ChangePointModel* cpm;
		
        if (strcmp(*cpmType, "Student") ==0) {
			cpm = new ChangePointModelT(vthresholds,*startup);
		} else if (strcmp(*cpmType, "Bartlett") ==0) {
			cpm = new ChangePointModelF(vthresholds,*startup);
		} else if (strcmp(*cpmType, "MW") ==0) {
			cpm = new ChangePointModelMW(vthresholds,*startup);
		} else if (strcmp(*cpmType, "Mood") ==0) {
			cpm = new ChangePointModelMood(vthresholds,*startup);
		} else if (strcmp(*cpmType, "LP") ==0) {
			cpm = new ChangePointModelLepage(vthresholds,*startup);
		} else if (strcmp(*cpmType, "FET") ==0) {
			cpm = new ChangePointModelFET(vthresholds,*startup,*lambda);
		} else if (strcmp(*cpmType, "Joint") ==0) {
			cpm = new ChangePointModelJointNormal(vthresholds,*startup);
		}  else if (strcmp(*cpmType, "JointAdjusted") ==0) {
        	cpm = new ChangePointModelJointNormalAdjusted(vthresholds,*startup);
		} else if (strcmp(*cpmType, "CVM") ==0) {
			cpm = new ChangePointModelCVM(vthresholds,*startup);
		} else if (strcmp(*cpmType, "KS") ==0) {
			cpm = new ChangePointModelKS(vthresholds,*startup);
		}  else if (strcmp(*cpmType, "Poisson") ==0) {
			cpm = new ChangePointModelPoisson(vthresholds,*startup);
		} else if (strcmp(*cpmType, "Exponential") ==0) {
			cpm = new ChangePointModelExponential(vthresholds,*startup);
		} else if (strcmp(*cpmType, "ExponentialAdjusted") ==0) {
			cpm = new ChangePointModelExponentialAdjusted(vthresholds,*startup);
		}
			
		else {
			Rprintf("Error: Change point model type not supported\n");
			return;
		}
        std::vector<int> cpsv, dtsv;
		cpm->processStream(xv, Dsv, cpsv,dtsv);	
        *numChanges = dtsv.size();
        //std::copy(Dsv.begin(),Dsv.end(),Ds);
        std::copy(cpsv.begin(),cpsv.end(),cps);
        std::copy(dtsv.begin(),dtsv.end(),dts);
		delete cpm;
		return;
	}
    
    
    void cpmDetectChangeBatch(char **cpmType, double *x, int* nx, double *Ds, double *lambda) {
        std::vector<double> Dsv, vthresholds;
        int startup=20;
        ChangePointModel* cpm;
		
        if (strcmp(*cpmType, "Student") ==0) {
			cpm = new ChangePointModelT(vthresholds,startup);
		} else if (strcmp(*cpmType, "Bartlett") ==0) {
			cpm = new ChangePointModelF(vthresholds,startup);
		} else if (strcmp(*cpmType, "MW") ==0) {
			cpm = new ChangePointModelMW(vthresholds,startup);
		} else if (strcmp(*cpmType, "Mood") ==0) {
			cpm = new ChangePointModelMood(vthresholds,startup);
		} else if (strcmp(*cpmType, "FET") ==0) {
			cpm = new ChangePointModelFET(vthresholds,startup,*lambda);
		} else if (strcmp(*cpmType, "LP") ==0) {
			cpm = new ChangePointModelLepage(vthresholds,startup);
		} else if (strcmp(*cpmType, "Joint") ==0) {
			cpm = new ChangePointModelJointNormal(vthresholds,startup);
		} else if (strcmp(*cpmType, "JointAdjusted") ==0) {
    		cpm = new ChangePointModelJointNormalAdjusted(vthresholds,startup);
		} else if (strcmp(*cpmType, "CVM") ==0) {
			cpm = new ChangePointModelCVM(vthresholds,startup);
		} else if (strcmp(*cpmType, "KS") ==0) {
			cpm = new ChangePointModelKS(vthresholds,startup);
		} else if (strcmp(*cpmType, "Poisson") ==0) {
			cpm = new ChangePointModelPoisson(vthresholds,startup);
		} else if (strcmp(*cpmType, "Exponential") ==0) {
			cpm = new ChangePointModelExponential(vthresholds,startup);
		} else if (strcmp(*cpmType, "ExponentialAdjusted") ==0) {
			cpm = new ChangePointModelExponentialAdjusted(vthresholds,startup);
		} else if (strcmp(*cpmType, "JointHawkins") ==0) {
			cpm = new ChangePointModelJointNormalHawkins(vthresholds,startup);
		} else {
			Rprintf("Change point model type not supported\n");
			return;
		}
		
        for (int i = 0; i < *nx; i++) {
            cpm->processPoint(x[i]);
        }
        cpm->cpmMLEaux(Dsv);
        std::copy(Dsv.begin(),Dsv.end(),Ds);
		delete cpm;
		return;
	}
}		
		
