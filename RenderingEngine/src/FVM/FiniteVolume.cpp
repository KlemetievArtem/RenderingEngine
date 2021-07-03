#include "FiniteVolume.h"

bool operator==(const Line& l1, const Line& l2) {
	if (((l1.m_p1 == l2.m_p1) and (l1.m_p1 == l2.m_p2)) or
		((l1.m_p2 == l2.m_p1) and (l1.m_p2 == l2.m_p2))) {
		return true;
	}
	else return false;
}

bool operator==(const Surface& s1, const Surface& s2) {
	if (s1.m_lineArray.size() != s2.m_lineArray.size()) return false;
	int a_size = s1.m_lineArray.size();
	int nrOfSameLines = 0;
	for (int i = 0;i < a_size;i++) {
		for (int ii = 0;ii < a_size;ii++) {
			if (s1.m_lineArray[i] == s2.m_lineArray[ii]) {
				nrOfSameLines++;
				break;
			}
		}
	}
	if (nrOfSameLines == a_size) return true;
	return false;
}
