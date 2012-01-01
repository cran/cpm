#include "ChangePointModelLepage.h"

ChangePointModelLepage::ChangePointModelLepage() {
	m_startup=20;
}


ChangePointModelLepage::ChangePointModelLepage(const std::vector<double> &thresholds, int startup) {
	m_thresholds = thresholds;
	m_startup = startup;
}


void ChangePointModelLepage::updateStatistics(const double &obs) {
	m_cpmMW.updateStatistics(obs);
	m_cpmMood.updateStatistics(obs);
}

void ChangePointModelLepage::cpmMLEaux(std::vector <double> &Us) {
	std::vector <double> UsMood;
	
	m_cpmMW.cpmMLEaux(Us);
	m_cpmMood.cpmMLEaux(UsMood);
	
	int sz = Us.size();
	for (int i = 1; i < sz-2; i++) {
		Us[i] = Us[i]*Us[i] + UsMood[i]*UsMood[i];
	}
}


void ChangePointModelLepage::reset() {
	m_cpmMW.reset();
	m_cpmMood.reset();
	n=0;
}

