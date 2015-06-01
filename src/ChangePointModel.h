#include <vector>
using namespace std;
 
#ifndef CHANGE_POINT_MODEL
#define CHANGE_POINT_MODEL


class ChangePointModel {
	public:

		std::vector<std::vector<double> > m_statistics; 
		std::vector<double> m_thresholds;
		int m_startup;
		long int n;		

    	ChangePointModel();
    	virtual ~ChangePointModel() = 0;
		
		void processPoint(const double &obs);
		void processStream(const std::vector<double> &x, std::vector<double> &Us, std::vector<int> &cps, std::vector<int> &dts);
		void detectChange(const std::vector<double> &x, std::vector<double> &Us, int &cp, int &dt);

		virtual void reset() {};
		void cpmMLE(double &maxvalue, int &maxindex);
		virtual void cpmMLEaux(std::vector <double> &Us) {};
		virtual void updateStatistics(const double &obs) {};	
		
        //fix to return reference
        std::vector<double> getThresholds() {return(m_thresholds);}
	private:
	
};

#endif


