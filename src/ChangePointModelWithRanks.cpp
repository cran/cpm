#include "ChangePointModelWithRanks.h"
#include <vector>
#include <R.h>
void ChangePointModelWithRanks::updateStatistics(const double &obs) {
	m_statistics[0].push_back(obs);
	double rank = 1;
	int sz = m_statistics[1].size();

    std::vector<int> ties;

	for (int i = 0; i < sz; i++) {
		if (m_statistics[0][i] > obs) {
			m_statistics[1][i]++;
		} else if (obs > m_statistics[0][i]) {
			rank++;
		} else {
			ties.push_back(i); //there is a tied rank
		}
    }

    //if there were any ties, then assign each observation to the average rank
    int tiessz = (int) ties.size();
    if (tiessz > 0) {
        rank = (rank + rank + tiessz)/2.0;
        for (int i = 0; i < tiessz; i++) {
            m_statistics[1][ties[i]] = rank;
        }
	}

	m_statistics[1].push_back(rank);
}
    
void ChangePointModelWithRanks::reset() {
    m_statistics[0].clear();
    m_statistics[1].clear();
    n=0;
}
    
	

