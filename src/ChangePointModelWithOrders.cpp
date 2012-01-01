#include "ChangePointModelWithOrders.h"
#include <R.h>

void ChangePointModelWithOrders::updateStatistics(const double &obs) {
	int sz = m_statistics[0].size();
	if (sz==0) {
		m_statistics[0].push_back(obs);
		m_statistics[1].push_back(1);
		return;
	}

	int rank = 0;
		
	for (int i = 0; i < sz; i++) {
		if (obs > m_statistics[0][i]) {
			++rank;
		}	
	}
	
	if (rank==sz) {
		m_statistics[1].push_back(sz+1);
	} else {
	 	std::vector<double>::iterator it = m_statistics[1].begin() + rank;
		m_statistics[1].insert(it,sz+1);
	}
	m_statistics[0].push_back(obs);	
}
    
void ChangePointModelWithOrders::reset() {
    m_statistics[0].clear();
    m_statistics[1].clear();
    n=0;
}
    
	

