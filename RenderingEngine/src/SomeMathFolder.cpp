
#include "SomeMathFolder.h"

double SurfaceArea(glm::vec3 p1, glm::vec3 p2, glm::vec3 p3) {
	double retVal;
	Matrix m1(3, 3);
	m1(0, 0) = 1.0; m1(1, 0) = p1.x; m1(2, 0) = p1.y;
	m1(0, 1) = 1.0; m1(1, 1) = p2.x; m1(2, 1) = p2.y;
	m1(0, 2) = 1.0; m1(1, 2) = p3.x; m1(2, 2) = p3.y;
	Matrix m2(3, 3);
	m2(0, 0) = 1.0; m2(1, 0) = p1.x; m2(2, 0) = p1.z;
	m2(0, 1) = 1.0; m2(1, 1) = p2.x; m2(2, 1) = p2.z;
	m2(0, 2) = 1.0; m2(1, 2) = p3.x; m2(2, 2) = p3.z;
	Matrix m3(3, 3);
	m3(0, 0) = 1.0; m3(1, 0) = p1.y; m3(2, 0) = p1.z;
	m3(0, 1) = 1.0; m3(1, 1) = p2.y; m3(2, 1) = p2.z;
	m3(0, 2) = 1.0; m3(1, 2) = p3.y; m3(2, 2) = p3.z;
	return retVal = 0.5*sqrt(pow(determinant(m1), 2) + pow(determinant(m2), 2) + pow(determinant(m3), 2));
}