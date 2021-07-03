#include "ComputationalDomain.h"


float getRandomNumber(float min, float max) {
	static const double fraction = 1.0 / (static_cast<double>(RAND_MAX) + 1.0);
	return (rand()*fraction*(max - min) + min);
}

glm::mat3x3 Xrot(float radians) {
	glm::mat3x3 Retval;
	Retval[0][0] = 1; Retval[0][1] = 0;            Retval[0][2] = 0;
	Retval[1][0] = 0; Retval[1][1] = cos(radians); Retval[1][2] = -sin(radians);
	Retval[2][0] = 0; Retval[2][1] = sin(radians); Retval[2][2] = cos(radians);
	return Retval;
};
glm::mat3x3 Yrot(float radians) {
	glm::mat3x3 Retval;
	Retval[0][0] = cos(radians);     Retval[0][1] = 0; Retval[0][2] = sin(radians);
	Retval[1][0] = 0;                Retval[1][1] = 1; Retval[1][2] = 0;
	Retval[2][0] = 0 - sin(radians); Retval[2][1] = 0; Retval[2][2] = cos(radians);
	return Retval;
};
glm::mat3x3 Zrot(float radians) {
	glm::mat3x3 Retval;
	Retval[0][0] = cos(radians); Retval[0][1] = -sin(radians); Retval[0][2] = 0;
	Retval[1][0] = sin(radians); Retval[1][1] = cos(radians);  Retval[1][2] = 0;
	Retval[2][0] = 0;            Retval[2][1] = 0;             Retval[2][2] = 1;
	return Retval;
};



void FVM_CD::Initilization(glm::vec3 velocity, std::vector<BoundaryBase*>* activeBoundaries) {
	if (computationalDomain_mode == MODE_CD::CD_DEBUG) {
		std::cout << "FVM_CD::Initilization\n";
	}

	//setStatMinScale(InitialSmR*m_options.smoothingKernelLengthCoefficient);
	//  d(Vx)/dt d(Density)/dt

	// Перенести в функцию EquationsInitialization
	//EquationsInitialization();



	cd_prec minX = (CD_Boundaries)[0]->getMinX();
	cd_prec minY = (CD_Boundaries)[0]->getMinY();
	cd_prec minZ = (CD_Boundaries)[0]->getMinZ();
	for (size_t i = 1; i < CD_Boundaries.size(); i++) {
		if ((CD_Boundaries)[i]->getMinX() < minX) minX = (CD_Boundaries)[i]->getMinX();
		if ((CD_Boundaries)[i]->getMinY() < minY) minY = (CD_Boundaries)[i]->getMinY();
		if ((CD_Boundaries)[i]->getMinZ() < minZ) minZ = (CD_Boundaries)[i]->getMinZ();
	}
	cd_prec maxX = (CD_Boundaries)[0]->getMaxX();
	cd_prec maxY = (CD_Boundaries)[0]->getMaxY();
	cd_prec maxZ = (CD_Boundaries)[0]->getMaxZ();
	for (size_t i = 1; i < activeBoundaries->size(); i++) {
		if ((CD_Boundaries)[i]->getMaxX() > maxX) maxX = (CD_Boundaries)[i]->getMaxX();
		if ((CD_Boundaries)[i]->getMaxY() > maxY) maxY = (CD_Boundaries)[i]->getMaxY();
		if ((CD_Boundaries)[i]->getMaxZ() > maxZ) maxZ = (CD_Boundaries)[i]->getMaxZ();
	}

	glm::vec3 positionMin(minX, minY, minZ);
	glm::vec3 positionMax(maxX, maxY, maxZ);



	MeshInitilization(positionMin, positionMax);
	BoundaryMeshInitialization(activeBoundaries);



}


void FVM_CD::MeshInitilization(glm::vec3 positionMin, glm::vec3 positionMax) {
	if (computationalDomain_mode == MODE_CD::CD_DEBUG) {
		std::cout << "FVM_CD::MeshInitilization\n";
	}
	int count = 0;
	fv_prec pressure = 10E05;
	switch (m_options.mesh_s) {
	case(STRUCTED):
	{
		switch (nrOfDim) {
		case(D1):
			m_options.average_dim_steps.x = (positionMax.x - positionMin.x) / m_options.nrOfFVinDir.x;
			for (int x = 0;x < m_options.nrOfFVinDir.x;x++) {
				fv_prec_3 p1({ positionMin.x + m_options.average_dim_steps.x*x, positionMin.y, positionMin.z });
				fv_prec_3 p2({ positionMin.x + m_options.average_dim_steps.x*(x + 1), positionMin.y, positionMin.z });

				Line tempLine(p1, p2);

				FVM.FiniteVolumes.push_back(new FiniteVolume(count, FINITE_VOLUME_TYPE::FVT_DEFAULT, tempLine, glm::vec3(0.0f, 0.f, 0.f), pressure));
				FVM.FiniteVolumeMesh.push_back(new FVMeshNodeEmpty(FVM.FiniteVolumes[count]));
				for (int fvc = 0;fvc < count;fvc++) {
					if (FVM.FiniteVolumes[count] == FVM.FiniteVolumes[fvc]) {
						if (FVM.FiniteVolumes[count]->m_l.ConnectedWith(FVM.FiniteVolumes[fvc]->m_l)) {
							std::cout << "adding FiniteVolumeMesh\n";
							FVM.FiniteVolumeMesh[count]->addNeighbour(FVM.FiniteVolumes[fvc]);
							FVM.FiniteVolumeMesh[fvc]->addNeighbour(FVM.FiniteVolumes[count]);
						}
					}
				}
				count++;
			}
			m_options.nrOfFV = count;
			break;
		case(D2):
			m_options.average_dim_steps.x = (positionMax.x - positionMin.x) / m_options.nrOfFVinDir.x;
			m_options.average_dim_steps.y = (positionMax.y - positionMin.y) / m_options.nrOfFVinDir.y;
			for (int y = 0;y < m_options.nrOfFVinDir.y;y++) {
				//std::cout << "x = " << x << "\n";
				for (int x = 0;x < m_options.nrOfFVinDir.x;x++) {
					//std::cout << "y = " << y << "\n";
					fv_prec_3 p1({ positionMin.x + m_options.average_dim_steps.x*x, positionMin.y + m_options.average_dim_steps.y*y, (positionMin.z+ positionMax.z)/2.0});
					fv_prec_3 p2({ positionMin.x + m_options.average_dim_steps.x*(x + 1), positionMin.y + m_options.average_dim_steps.y*y, (positionMin.z + positionMax.z) / 2.0 });
					fv_prec_3 p3({ positionMin.x + m_options.average_dim_steps.x*x, positionMin.y + m_options.average_dim_steps.y*(y + 1), (positionMin.z + positionMax.z) / 2.0 });
					fv_prec_3 p4({ positionMin.x + m_options.average_dim_steps.x*(x + 1), positionMin.y + m_options.average_dim_steps.y*(y + 1), (positionMin.z + positionMax.z) / 2.0 });
					//std::cout << "{" << p1.x << "," << p1.y << "}, {" << p2.x << "," << p2.y << "}, {" << p3.x << "," << p3.y << "}, {" << p4.x << "," << p4.y << "}\n";
					Surface tempSurface(Line(p1, p2), Line(p2, p4), Line(p4, p3)); tempSurface.addLine(Line(p3, p1));
					tempSurface.calcCentre();
					FVM.FiniteVolumes.push_back(new FiniteVolume(count, FINITE_VOLUME_TYPE::FVT_DEFAULT, tempSurface, glm::vec3(0.0f, 0.f, 0.f), pressure));
					FVM.FiniteVolumeMesh.push_back(new FVMeshNodeEmpty(FVM.FiniteVolumes[count]));
					for (int fvc = 0;fvc < count;fvc++) {
						//std::cout << "TUT";
						//std::cout << FVM.FiniteVolumes[count]->m_s.getCentre().x << ", " << FVM.FiniteVolumes[count]->m_s.getCentre().y << "\n";
						if (FVM.FiniteVolumes[count]->m_s.ConnectedWith(FVM.FiniteVolumes[fvc]->m_s)) {
							std::cout << "adding FiniteVolumeMesh\n";
							FVM.FiniteVolumeMesh[count]->addNeighbour(FVM.FiniteVolumes[fvc]);
							FVM.FiniteVolumeMesh[fvc]->addNeighbour(FVM.FiniteVolumes[count]);
						}
					}
					count++;
				}
			}
			m_options.nrOfFV = count;
			break;
		case(D3):
			m_options.average_dim_steps.x = (positionMax.x - positionMin.x) / m_options.nrOfFVinDir.x;
			m_options.average_dim_steps.y = (positionMax.y - positionMin.y) / m_options.nrOfFVinDir.y;
			m_options.average_dim_steps.z = (positionMax.z - positionMin.z) / m_options.nrOfFVinDir.z;
			for (int z = 0;z < m_options.nrOfFVinDir.z;z++) {
				for (int y = 0;y < m_options.nrOfFVinDir.y;y++) {
					for (int x = 0;x < m_options.nrOfFVinDir.x;x++) {
						fv_prec_3 p1({ positionMin.x + m_options.average_dim_steps.x*x, positionMin.y + m_options.average_dim_steps.y*y, positionMin.z + m_options.average_dim_steps.z*z });
						fv_prec_3 p2({ positionMin.x + m_options.average_dim_steps.x*(x + 1), positionMin.y + m_options.average_dim_steps.y*y, positionMin.z + m_options.average_dim_steps.z*z });
						fv_prec_3 p3({ positionMin.x + m_options.average_dim_steps.x*x, positionMin.y + m_options.average_dim_steps.y*(y + 1), positionMin.z + m_options.average_dim_steps.z*z });
						fv_prec_3 p4({ positionMin.x + m_options.average_dim_steps.x*(x + 1), positionMin.y + m_options.average_dim_steps.y*(y + 1), positionMin.z + m_options.average_dim_steps.z*z });
						fv_prec_3 p5({ positionMin.x + m_options.average_dim_steps.x*x, positionMin.y + m_options.average_dim_steps.y*y, positionMin.z + m_options.average_dim_steps.z*(z + 1) });
						fv_prec_3 p6({ positionMin.x + m_options.average_dim_steps.x*(x + 1), positionMin.y + m_options.average_dim_steps.y*y, positionMin.z + m_options.average_dim_steps.z*(z + 1) });
						fv_prec_3 p7({ positionMin.x + m_options.average_dim_steps.x*x, positionMin.y + m_options.average_dim_steps.y*(y + 1), positionMin.z + m_options.average_dim_steps.z*(z + 1) });
						fv_prec_3 p8({ positionMin.x + m_options.average_dim_steps.x*(x + 1), positionMin.y + m_options.average_dim_steps.y*(y + 1), positionMin.z + m_options.average_dim_steps.z*(z + 1) });

						Surface s1(Line(p1, p2), Line(p2, p4), Line(p4, p3)); s1.addLine(Line(p3, p1));
						Surface s2(Line(p1, p5), Line(p5, p6), Line(p6, p2)); s1.addLine(Line(p2, p1));
						Surface s3(Line(p2, p6), Line(p6, p8), Line(p8, p4)); s1.addLine(Line(p4, p2));
						Surface s4(Line(p4, p8), Line(p8, p7), Line(p7, p3)); s1.addLine(Line(p3, p4));
						Surface s5(Line(p3, p7), Line(p7, p5), Line(p5, p1)); s1.addLine(Line(p1, p3));
						Surface s6(Line(p5, p6), Line(p6, p8), Line(p8, p7)); s1.addLine(Line(p7, p5));

						Volume tempVolume(s1, s2, s3, s4); tempVolume.addSurface(s5); tempVolume.addSurface(s6); tempVolume.calcCentre();

						FVM.FiniteVolumes.push_back(new FiniteVolume(count, FINITE_VOLUME_TYPE::FVT_DEFAULT, tempVolume, glm::vec3(0.0f, 0.f, 0.f), pressure));
						FVM.FiniteVolumeMesh.push_back(new FVMeshNodeEmpty(FVM.FiniteVolumes[count]));
						for (int fvc = 0;fvc < count;fvc++) {
							if (FVM.FiniteVolumes[count] == FVM.FiniteVolumes[fvc]) {
								if (FVM.FiniteVolumes[count]->m_v.ConnectedWith(FVM.FiniteVolumes[fvc]->m_v)) {
									std::cout << "adding FiniteVolumeMesh\n";
									FVM.FiniteVolumeMesh[count]->addNeighbour(FVM.FiniteVolumes[fvc]);
								}
							}
						}
						count++;
					}
				}
			}
			m_options.nrOfFV = count;
			break;
		default:
			break;
		}
	}
		break;
	default:
		break;
	}

	if (computationalDomain_mode == MODE_CD::CD_DEBUG) {
		std::cout << "number of Finite Volume = " << m_options.nrOfFV << "\n";
	}
}
/*
void FVM_CD::BoundaryMeshInitialization(std::vector<BoundaryBase*>* activeBoundaries) {


	float Xstart;
	float Ystart;

	glm::vec3 AverageNormal;
	int BoundaryVoluemCountOld = FVM.FiniteVolumes.size();
	int BoundaryVoluemCount = FVM.FiniteVolumes.size();
	fv_prec pressure = 10E05;
	for (auto*& i : (*activeBoundaries)) {
		int nrOfPolygons = static_cast<int>(i->ReturnMesh()->getNumberOfIndice()) / 3;

		glm::vec3 OXYnormal = { 0.f,0.f,1.f };
		float PlaneCoefA;
		float PlaneCoefB;
		float PlaneCoefC;
		float PlaneCoefD;

		glm::vec3 Plane1;
		glm::vec3 Plane2;
		glm::vec3 Plane3;

		for (int p = 0;p < nrOfPolygons;p++) {
			//WE WANT TO ROTATE EVERY POLYGON IN A WAY, THAT TWO OF HIS SIDES COULD BE ANALISED AS LINEAR FUNCTIONS 
			GLuint* indexarray = i->ReturnMesh()->getIndexArray();
			int p1 = indexarray[p * 3 + 0];
			int p2 = indexarray[p * 3 + 1];
			int p3 = indexarray[p * 3 + 2];
			Vertex* vertexArray = i->ReturnMesh()->getVertexArray();
			glm::vec3 pos1 = vertexArray[p1].position;
			glm::vec3 pos2 = vertexArray[p2].position;
			glm::vec3 pos3 = vertexArray[p3].position;

			AverageNormal = (vertexArray[p1].normal + vertexArray[p2].normal + vertexArray[p3].normal);
			AverageNormal[0] /= 3;
			AverageNormal[1] /= 3;
			AverageNormal[2] /= 3;

			Plane1 = pos1;
			Plane2 = pos2;
			Plane3 = pos3;
			//Зная три точки, можно определить функцию плоскости
			PlaneCoefA = (Plane2.y - Plane1.y)*(Plane3.z - Plane1.z) - (Plane3.y - Plane1.y)*(Plane2.z - Plane1.z);
			PlaneCoefB = -(Plane2.x - Plane1.x)*(Plane3.z - Plane1.z) - (Plane3.x - Plane1.x)*(Plane2.z - Plane1.z);
			PlaneCoefC = (Plane2.x - Plane1.x)*(Plane3.y - Plane1.y) - (Plane3.x - Plane1.x)*(Plane2.y - Plane1.y);
			PlaneCoefD = -Plane1.x*PlaneCoefA - Plane1.y*PlaneCoefB - Plane1.z*PlaneCoefC;
			//std::cout << PlaneCoefA << "*x + " << PlaneCoefB << "*y + " << PlaneCoefC << "*z + " << PlaneCoefD << " =  0\n";
			//Далее определим, линию пересечение полученной плоскости с плоскостью OXY  (z=0):
			// x*PlaneCoefA + y*PlaneCoefB + PlaneCoefD = 0
			//Определим пересечение этой линии с каждой линией, из которых состоит полигон
			glm::vec3 interP;
			std::vector<glm::vec3> PolyLine0 = { pos1 ,pos2 ,pos3 };
			std::vector<glm::vec3> PolyLine1 = { pos2 ,pos3 ,pos1 };
			std::vector<glm::vec3> interPoints;

			float t;
			float a;
			for (int side = 0;side < 3;side++) {
				if (PlaneCoefA != 0) {
					a = -PolyLine0[side].z / (PolyLine1[side].z - PolyLine0[side].z);
					t = a * (PolyLine1[side].y - PolyLine0[side].y) + PolyLine0[side].y;
					if (-PlaneCoefD / PlaneCoefA - t * PlaneCoefB / PlaneCoefA == a * (PolyLine1[side].x - PolyLine0[side].x) + PolyLine0[side].x) {

						interP.x = -PlaneCoefD / PlaneCoefA - t * PlaneCoefB / PlaneCoefA;
						interP.y = t;
						interP.z = 0.f;
						interPoints.push_back(interP);
					}
				}
				else {
					a = -PolyLine0[side].z / (PolyLine1[side].z - PolyLine0[side].z);
					t = a * (PolyLine1[side].x - PolyLine0[side].x) + PolyLine0[side].x;
					if (-PlaneCoefD / PlaneCoefB - t * PlaneCoefA / PlaneCoefB == a * (PolyLine1[side].y - PolyLine0[side].y) + PolyLine0[side].y) {

						interP.x = t;
						interP.y = -PlaneCoefD / PlaneCoefB - t * PlaneCoefA / PlaneCoefB;
						interP.z = 0.f;
						interPoints.push_back(interP);
					}
				}
			}
			//std::cout << "GLOBAL COORDINATS\n";
			glm::vec3 newCoordOrigin = pos1;
			pos1 += -newCoordOrigin;
			pos2 += -newCoordOrigin;
			pos3 += -newCoordOrigin;

			cd_prec fi_x, fi_y, fi_z = 0.f;
			bool rotationFlags[] = { false,false, false };
			// Z-ROTATION
			if ((pos2.y != 0) or (pos2.x < 0)) {
				fi_z = (asin(pos2.y / sqrt(pow(pos2.y, 2) + pow(pos2.x, 2))));
				if ((pos2.x < 0))
					fi_z = M_PI - fi_z;
				pos2 = Zrot(fi_z) * pos2;
				pos3 = Zrot(fi_z) * pos3;

				rotationFlags[2] = true;
			}
			// Y-ROTATION
			if ((pos2.z != 0) or (pos2.x < 0)) {
				fi_y = (asin(pos2.z / sqrt(pow(pos2.x, 2) + pow(pos2.z, 2))));
				if ((pos2.x < 0))
					fi_y += M_PI;
				pos2 = Yrot(-fi_y) * pos2;
				pos3 = Yrot(-fi_y) * pos3;

				rotationFlags[1] = true;
			}
			// X-ROTATION
			if ((pos3.z != 0) or (pos3.y < 0)) {
				fi_x = (asin(pos3.z / sqrt(pow(pos3.z, 2) + pow(pos3.y, 2))));
				if (pos3.y*pos3.z < 0) {
					fi_x = -fi_x;
				}
				if ((pos3.y < 0))
					fi_x += M_PI;
				pos2 = Xrot(fi_x) * pos2;
				pos3 = Xrot(fi_x) * pos3;
				rotationFlags[0] = true;
			}
			for (auto &ip : interPoints) {
				ip += -newCoordOrigin;
				if (rotationFlags[2] == true)
					ip = Zrot(fi_z) * ip;
				if (rotationFlags[1] == true)
					ip = Yrot(-fi_y) * ip;
				if (rotationFlags[0] == true)
					ip = Xrot(fi_x) * ip;
			}
			LinearFunc LF1_3 = LinearFuncCoefficients(glm::vec2(pos1.x, pos1.y), glm::vec2(pos3.x, pos3.y));
			LinearFunc LF2_3 = LinearFuncCoefficients(glm::vec2(pos2.x, pos2.y), glm::vec2(pos3.x, pos3.y));
			std::cout << "{" << pos1.x << "," << pos1.y << "," << pos1.z << "},{" << pos2.x << "," << pos2.y << "," << pos2.z << "},{" << pos3.x << "," << pos3.y << "," << pos3.z << "}\n";
			LinearFunc LFIntersection = LinearFuncCoefficients(glm::vec2(interPoints[0].x, interPoints[0].y), glm::vec2(interPoints[1].x, interPoints[1].y));

			glm::vec3 mincoord = pos1;
			glm::vec3 maxcoord = pos1;
			// MAX AND MIN VALUES PER POLYGON
			for (int dim = 0;dim < 3;dim++) {
				if (mincoord[dim] > pos2[dim])
					mincoord[dim] = pos2[dim];
				if (mincoord[dim] > pos3[dim])
					mincoord[dim] = pos3[dim];
				if (maxcoord[dim] < pos2[dim])
					maxcoord[dim] = pos2[dim];
				if (maxcoord[dim] < pos3[dim])
					maxcoord[dim] = pos3[dim];
			}
			fv_prec_3 rand_pos_min;
			fv_prec_3 rand_pos;
			fv_prec_3 rand_pos_max;
			cd_prec partB_x = 0.f;
			cd_prec partB_y = 0.f;
			cd_prec partB_z = 0.f;

			int BNx = (m_options.nrOfFVinDir.x / nrOfPolygons);
			int BNy = (m_options.nrOfFVinDir.y / nrOfPolygons);
			fv_prec_3 Bd;
			fv_prec_3 Bd_rot;
			Bd.x = (maxcoord[0] - mincoord[0]) / BNx;
			Bd.y = (maxcoord[1] - mincoord[1]) / BNy;
			Bd.z = 0.0;

			//FUTURE WARNING
			if (Bd.x > Bd.y * 2)
				Bd.x /= 2;
			else if (Bd.y > Bd.x * 2)
				Bd.y /= 2;
			//FUTURE WARNING
			fv_prec_3 min2D = interPoints[0];
			fv_prec_3 max2D = interPoints[0];
			if (nrOfDim == D2) {
				if (min2D.x > interPoints[1].x)  min2D.x = interPoints[1].x;
				if (max2D.x < interPoints[1].x)  max2D.x = interPoints[1].x;
				if (min2D.y > interPoints[1].y)  min2D.y = interPoints[1].y;
				if (max2D.y < interPoints[1].y)  max2D.y = interPoints[1].y;
			}
			//std::cout << "max2D_x:" << max2D.x << ", min2D_x:" << min2D.x << "\n";
			//std::cout << "max2D_y:" << max2D.y << ", min2D_y:" << min2D.y << "\n";
			if (nrOfDim == D2) {
				BNx = (m_options.nrOfFVinDir.x / nrOfPolygons);
				BNy = BNx;
				Bd.x = (max2D.x - min2D.x) / ((BNx));
				Bd.y = (max2D.y - min2D.y) / ((BNy));
			}
			//std::cout << "Bdx:" << Bd.x << ", Bdy:" << Bd.y << "\n";
			//std::cout << "Bnx:" << BNx << ", Bny:" << BNy << "\n";
			int Bnx = 0;
			int Bny = 0;
			for (int j = 0;j < BNx;) {
				if (nrOfDim == D3) {
					Xstart = (partB_y - LF1_3.b) / LF1_3.k;
					Ystart = 0.f;
					partB_x = Xstart + Bnx * Bd.x;
					partB_y = Ystart + Bny * Bd.y;
					partB_z = 0;

					Bnx++;
					if (partB_x > (partB_y - LF2_3.b) / LF2_3.k) {
						Bnx = 0;
						Bny++;
						if (partB_y > maxcoord[1]) {
							partB_y = 0.f;
							j = static_cast<int>(m_options.nrOfFVinDir.y / nrOfPolygons);
							break;
						}
					}
				}
				if (nrOfDim == D2) {
					partB_z = 0.f;
					Xstart = min2D.x + Bd.x * 0.5;
					Ystart = min2D.y + Bd.y * 0.5;
					partB_x = Xstart + Bnx * Bd.x;

					if (abs(min2D.x - max2D.x) < 10E-08) {
						partB_y = Ystart + Bnx * Bd.y;
					}
					else {
						partB_y = LFIntersection.k*partB_x + LFIntersection.b;
					}
					Bnx++;
				}
				rand_pos_min = { partB_x - Bd.x / 2.f,partB_y - Bd.y/ 2.f, partB_z };
				rand_pos = { partB_x, partB_y, partB_z };
				rand_pos_max = { partB_x + Bd.x / 2.f,partB_y + Bd.y / 2.f, partB_z };

				//std::cout << "{" << rand_pos_min.x << "," << rand_pos_min.y << "," << rand_pos_min.z << "} ";
				//std::cout << "{" << rand_pos.x << "," << rand_pos.y << "," << rand_pos.z << "} ";
				//std::cout << "{" << rand_pos_max.x << "," << rand_pos_max.y << "," << rand_pos_max.z << "}\n";
				

				bool flagIn = false;
				Bd_rot = Bd;
				if ((rand_pos.x < pos3.x) and (rand_pos.x > pos1.x)) {
					if (rand_pos.y < LF1_3.k*rand_pos.x + LF1_3.b) {
						flagIn = true;
						//std::cout << "under pos3\n";
					}
				}
				else if ((rand_pos.x > pos3.x) and (rand_pos.x < pos2.x)) {
					if (rand_pos.y < LF2_3.k*rand_pos.x + LF2_3.b) {
						flagIn = true;
						//std::cout << "over pos3\n";
					}
				}
				if (flagIn) {
					if (rotationFlags[0]) {
						//rand_pos_min = Xrot(-fi_x)* rand_pos_min;
						rand_pos = Xrot(-fi_x)* rand_pos;
						//rand_pos_max = Xrot(-fi_x)* rand_pos_max;
						Bd_rot = Xrot(-fi_x)*Bd_rot;
					}
					if (rotationFlags[1]) {
						//rand_pos_min = Yrot(fi_y)* rand_pos_min;
						rand_pos = Yrot(fi_y)* rand_pos;
						//rand_pos_max = Yrot(fi_y)* rand_pos_max;
						Bd_rot = Yrot(fi_y)* Bd_rot;
					}
					if (rotationFlags[2]) {
						//rand_pos_min = Zrot(-fi_z)* rand_pos_min;
						rand_pos = Zrot(-fi_z)* rand_pos;
						//rand_pos_max = Zrot(-fi_z)* rand_pos_max;
						Bd_rot = Zrot(-fi_z)* Bd_rot;
					}
					Bd_rot = sqrt(Bd_rot * Bd_rot);
					//std::cout << "{" << Bd_rot.x << "," << Bd_rot.y << "," << Bd_rot.z << "}\n";
					rand_pos += newCoordOrigin;
					rand_pos_min = rand_pos - Bd_rot / 2.f;
					rand_pos_max = rand_pos + Bd_rot / 2.f;
					if (abs(rand_pos_min.x) < 10E-08) rand_pos_min.x = 0.0;
					if (abs(rand_pos_min.y) < 10E-08) rand_pos_min.y = 0.0;
					if (abs(rand_pos_min.z) < 10E-08) rand_pos_min.z = 0.0;
					if (abs(rand_pos_max.x) < 10E-08) rand_pos_max.x = 0.0;
					if (abs(rand_pos_max.y) < 10E-08) rand_pos_max.y = 0.0;
					if (abs(rand_pos_max.z) < 10E-08) rand_pos_max.z = 0.0;

					std::cout << "{" << rand_pos_min.x << "," << rand_pos_min.y << "," << rand_pos_min.z << "} ";
					//std::cout << "{" << rand_pos.x << "," << rand_pos.y << "," << rand_pos.z << "} ";
					std::cout << "{" << rand_pos_max.x << "," << rand_pos_max.y << "," << rand_pos_max.z << "}\n";

					glm::vec3 bc_velocity = i->getVelocity();
					//std::cout << "vel:" << bc_velocity.x << "," << bc_velocity.y << "," << bc_velocity.z << "\n";
					if (i->isSource()) i->setPositions(rand_pos);

					if (nrOfDim == D2) {
						Line tempLine(rand_pos_min, rand_pos_max);
						FVM.FiniteVolumes.push_back(new FiniteVolume(BoundaryVoluemCount, FINITE_VOLUME_TYPE::FVT_BOUNDARY, tempLine, bc_velocity, pressure));
						FVM.FiniteVolumeMesh.push_back(new FVMeshNodeEmpty(FVM.FiniteVolumes[BoundaryVoluemCount]));
						for (int fvc = 0;fvc < BoundaryVoluemCount;fvc++) {
							for (auto l : FVM.FiniteVolumes[fvc]->m_s.m_lineArray) {
								if (FVM.FiniteVolumes[BoundaryVoluemCount]->m_l == l) {
									std::cout << "adding FiniteVolumeMesh\n";
									FVM.FiniteVolumeMesh[BoundaryVoluemCount]->addNeighbour(FVM.FiniteVolumes[fvc]);
									FVM.FiniteVolumeMesh[fvc]->addNeighbour(FVM.FiniteVolumes[BoundaryVoluemCount]);
								}
							}
						}
					}
					j++; BoundaryVoluemCount++;

				}
			}
			vertexArray = nullptr;
			indexarray = nullptr;
		}
	}
	std::cout << "FVM_CD::Initilization::BoundaryMeshInitialization::" << BoundaryVoluemCount - BoundaryVoluemCountOld << "\n";
}
*/
void FVM_CD::BoundaryMeshInitialization(std::vector<BoundaryBase*>* activeBoundaries) {
	int BoundaryVoluemCountOld = FVM.FiniteVolumes.size();
	int BoundaryVoluemCount = FVM.FiniteVolumes.size();
	fv_prec pressure = 10E05;
	bool has_a_nighbour;
	for (auto fvm1 : FVM.FiniteVolumeMesh) {
		for (auto fvm1_l : fvm1->FVPtr()->m_s.m_lineArray) {
			//std::cout << "fvm1_l: {" << fvm1_l.m_p1.x << "," << fvm1_l.m_p1.y << "},{" << fvm1_l.m_p2.x << "," << fvm1_l.m_p2.y << "}\n";
			has_a_nighbour = false;
			for (auto fvm2 : FVM.FiniteVolumeMesh) {
				if (fvm1 == fvm2) { continue; }
				for (auto fvm2_l : fvm2->FVPtr()->m_s.m_lineArray) {
					//std::cout << "fvm2_l: {" << fvm2_l.m_p1.x << "," << fvm2_l.m_p1.y << "},{" << fvm2_l.m_p2.x << "," << fvm2_l.m_p2.y << "}\n";
					if (fvm1_l == fvm2_l) {
						has_a_nighbour = true;
						break;
					}
				}
				if (has_a_nighbour) {
					break;
				}
			}
			if (!has_a_nighbour) {
				glm::vec3 AverageNormal;
				fv_prec_3 bc_velocity;
				bool middle_is_on_a_polygon;
				for (auto i : (*activeBoundaries)) {
					int nrOfPolygons = static_cast<int>(i->ReturnMesh()->getNumberOfIndice()) / 3;

					glm::vec3 OXYnormal = { 0.f,0.f,1.f };
					float PlaneCoefA;
					float PlaneCoefB;
					float PlaneCoefC;
					float PlaneCoefD;

					glm::vec3 Plane1;
					glm::vec3 Plane2;
					glm::vec3 Plane3;

					for (int p = 0;p < nrOfPolygons;p++) {
						middle_is_on_a_polygon = false;
						//WE WANT TO ROTATE EVERY POLYGON IN A WAY, THAT TWO OF HIS SIDES COULD BE ANALISED AS LINEAR FUNCTIONS 
						GLuint* indexarray = i->ReturnMesh()->getIndexArray();
						int p1 = indexarray[p * 3 + 0];
						int p2 = indexarray[p * 3 + 1];
						int p3 = indexarray[p * 3 + 2];
						Vertex* vertexArray = i->ReturnMesh()->getVertexArray();
						glm::vec3 pos1 = vertexArray[p1].position;
						glm::vec3 pos2 = vertexArray[p2].position;
						glm::vec3 pos3 = vertexArray[p3].position;

						AverageNormal = (vertexArray[p1].normal + vertexArray[p2].normal + vertexArray[p3].normal);
						AverageNormal[0] /= 3;
						AverageNormal[1] /= 3;
						AverageNormal[2] /= 3;

						Plane1 = pos1;
						Plane2 = pos2;
						Plane3 = pos3;
						//Зная три точки, можно определить функцию плоскости
						PlaneCoefA = (Plane2.y - Plane1.y)*(Plane3.z - Plane1.z) - (Plane3.y - Plane1.y)*(Plane2.z - Plane1.z);
						PlaneCoefB = -(Plane2.x - Plane1.x)*(Plane3.z - Plane1.z) - (Plane3.x - Plane1.x)*(Plane2.z - Plane1.z);
						PlaneCoefC = (Plane2.x - Plane1.x)*(Plane3.y - Plane1.y) - (Plane3.x - Plane1.x)*(Plane2.y - Plane1.y);
						PlaneCoefD = -Plane1.x*PlaneCoefA - Plane1.y*PlaneCoefB - Plane1.z*PlaneCoefC;
						if (abs(PlaneCoefA*fvm1_l.getCentre().x + PlaneCoefB * fvm1_l.getCentre().y + PlaneCoefC * fvm1_l.getCentre().z + PlaneCoefD) < 10E-05) {
							middle_is_on_a_polygon = true;
							break;
						}
					}
					if (middle_is_on_a_polygon) {
						bc_velocity = i->getVelocity();
						Line tempLine(fvm1_l);
						FVM.FiniteVolumes.push_back(new FiniteVolume(BoundaryVoluemCount, FINITE_VOLUME_TYPE::FVT_BOUNDARY, tempLine, bc_velocity, pressure));
						FVM.FiniteVolumes[BoundaryVoluemCount]->assignToBoundary(i);
						std::cout << FVM.FiniteVolumes[BoundaryVoluemCount]->fv_boundary << "\n";
						FVM.FiniteVolumeMesh.push_back(new FVMeshNodeEmpty(FVM.FiniteVolumes[BoundaryVoluemCount]));
						for (int fvc = 0;fvc < BoundaryVoluemCount;fvc++) {
							for (auto l : FVM.FiniteVolumes[fvc]->m_s.m_lineArray) {
								if (FVM.FiniteVolumes[BoundaryVoluemCount]->m_l == l) {
									std::cout << "adding FiniteVolumeMesh\n";
									FVM.FiniteVolumeMesh[BoundaryVoluemCount]->addNeighbour(FVM.FiniteVolumes[fvc]);
									FVM.FiniteVolumeMesh[fvc]->addNeighbour(FVM.FiniteVolumes[BoundaryVoluemCount]);
								}
							}
						}
						BoundaryVoluemCount++;
						break;
					}
				}
			}
		}
	}
}

void FVM_CD::PRB_refresh() {
	if (computationalDomain_mode == MODE_CD::CD_DEBUG) {
		std::cout << "FVM_CD::PRB_refresh\n";
	}
	PRB.clear();
	fv_prec_3 offset;
	fv_prec_3 size;
	fv_prec_3 rotation = { 0.0, 0.0, 0.0 };
	for (auto*& i : FVM.FiniteVolumes) {
		switch (i->FVdim) {
		case(FVA_S):
			offset = 0.5f*m_options.average_dim_steps;
			size = m_options.average_dim_steps;
			break;
		case(FVA_L):
			//rotation =  СЛОЖНО ЕСЛИ ПО-ПРАВИЛЬНОМУ
			if (abs((i->m_l.m_p2 - i->m_l.m_p1).x) >= abs((i->m_l.m_p2 - i->m_l.m_p1).y)) {
				offset.x = 0.5f*m_options.average_dim_steps.x;
				size.x = m_options.average_dim_steps.x;
				offset.y = 0.125f*m_options.average_dim_steps.y;
				size.y = 0.25f*m_options.average_dim_steps.y;
			}
			else {
				offset.y = 0.5f*m_options.average_dim_steps.y;
				size.y = m_options.average_dim_steps.y;
				offset.x = 0.125f*m_options.average_dim_steps.x;
				size.x = 0.25f*m_options.average_dim_steps.x;
			}
			break;
		default:
			break;
		}
		PRB.push_back({ i->m_Centre - offset, i->m_color, size });
	}
}

void FVM_CD::punctualColorChange(int number, glm::vec3 color)
{
}

float FVM_CD::getGlobalStats(int number)
{
	return 0.0f;
}

std::string FVM_CD::getLocalStats(int number)
{
	return std::string();
}


void FVM_CD::UpdateRendering(std::vector<Model*>* models, Texture* tex, Texture* tex_specualar, std::vector<Material*>* materials) {
	if (computationalDomain_mode == MODE_CD::CD_DEBUG) {
		std::cout << "FVM_CD::UpdateRendering\n";
	}
	Mesh* mesh;
	Material* material;

	if (m_options.firstCycle == true) {
		m_options.firstCycle = false;

		switch (m_options.mesh_s) {
		case(STRUCTED):
		{
			for (int i = 0;i < PRB.size();i++) {
				//std::cout << i << "\n";
				//(*models)[ModelId]->meshes.push_back(new Mesh(&Quad(glm::vec3(-0.5f), Quad::QuadNormal(Quad::Z, Quad::PLUS), 1.f, 1.f, PRB[i].color), PRB[i].position, glm::vec3(0.f), glm::vec3(0.f), glm::vec3(PRB[i].size)));
				
				switch (nrOfDim) {
				case(D1):
					mesh = (new Mesh(&Quad(glm::vec3(-0.5f), Quad::QuadNormal(Quad::Z, Quad::PLUS), 1.f, 1.f, glm::vec3(0.3f)), PRB[i].position, glm::vec3(0.f), PRB[i].rotation, glm::vec3(PRB[i].size)));
					break;
				case(D2):
					mesh = (new Mesh(&Quad(glm::vec3(0.f), Quad::QuadNormal(Quad::Z, Quad::PLUS), 1.0, 1.0, glm::vec3(0.3f)), PRB[i].position, glm::vec3(0.f), PRB[i].rotation, PRB[i].size));
					break;
				case(D3):
					mesh = (new Mesh(&Qube(PRB[i].color), PRB[i].position, glm::vec3(0.f), PRB[i].rotation, glm::vec3(PRB[i].size)));
					break;
				default:
					mesh = (new Mesh(&Quad(glm::vec3(-0.5f), Quad::QuadNormal(Quad::Z, Quad::PLUS), 1.f, 1.f, glm::vec3(0.3f)), PRB[i].position, glm::vec3(0.f), PRB[i].rotation, glm::vec3(PRB[i].size)));
					break;
				}
				material = (new Material(PRB[i].color, glm::vec3(0.9f), glm::vec3(1.f), 0, 1));
				(*models).push_back(new Model(PRB[i].position, material, tex, tex_specualar, mesh));
				//int last_model = models->size()-1;
				(*models)[ModelId + i]->getMaterial()->ChangeLighting(PRB[i].color, glm::vec3(0.9f), glm::vec3(1.f));
				(*models)[ModelId + i]->moveTo(PRB[i].position);
				(*models)[ModelId + i]->scaleUpTo(glm::vec3(PRB[i].size));
				//std::cout << ModelId + i << "   " << (*models).size()-1 << "\n";
				delete material;
				delete mesh;
			}
			//std::cout << "after: " << models->size() << "\n";

			assignNewModels(ModelId + m_options.nrOfFV);


		}
		break;
		default:
			break;
		}


		//std::cout << "before: " << models->size() << "\n";

		//for (auto*& i : (*models)[ModelId]->meshes) {
		//	delete i;
		//}
		//(*models)[ModelId]->meshes.resize(0);
	}
	else {
		for (int i = 0;i < PRB.size();i++) {
			//(*models)[ModelId]->meshes[i]->changeColorTo(PRB[i].color);
			//(*models)[ModelId]->meshes[i]->moveTo(PRB[i].position);
			//(*models)[ModelId]->meshes[i]->changeScaleTo(glm::vec3(PRB[i].size));
			//(*models)[ModelId]->meshes[i]->update();
			//std::cout << ModelId + i << "\n";
			(*models)[ModelId + i]->getMaterial()->ChangeLighting(PRB[i].color, glm::vec3(0.9f), glm::vec3(1.f));
			(*models)[ModelId + i]->moveTo(PRB[i].position);
			(*models)[ModelId + i]->scaleUpTo(glm::vec3(PRB[i].size));
			//std::cout << ModelId + i << "   " << (*models).size() - 1 << "\n";
		}
		//std::cout << "after: " << models->size() << "\n";
	}
}



void FVM_CD::AfterRendering(std::vector<Model*>* models) {
}

void FVM_CD::timeStep(cd_prec dt) {
	if (computationalDomain_mode == MODE_CD::CD_DEBUG) {
		std::cout << "FVM_CD::timeStep\n";
	}
	{
		//RecalcPrep();
		//ColoringBtType();
		Coloring();
		//BoundariesUpdate(InitialSmR, Dens0);
		PRB_refresh();
	}
	timeUpdate(dt);
}

void FVM_CD::timeStep_thread(cd_prec dt, std::atomic<bool>& dataReadyForRender, std::atomic<bool>& dataIsRendering) {
	if (computationalDomain_mode == MODE_CD::CD_DEBUG) {
		std::cout << "FVM_CD::timeStep_thread\n";
	}
	{
		//RecalcPrep();
		//ColoringBtType();
		Coloring();
		dataReadyForRender.store(false);
		while (dataIsRendering) {}
		//std::cout << "SPH_CD::timeStepEnd::dataReadyForRender:" << dataReadyForRender << "\n";
		PRB_refresh();
		dataReadyForRender.store(true);
		//std::cout << "SPH_CD::timeStepEnd::dataReadyForRender:" << dataReadyForRender << "\n";
	}
	timeUpdate(dt);
}



void FVM_CD::Solve(std::vector<FVMeshNode*>* Mesh, std::vector<fv_prec>* phi_param) {

}

void FVM_CD::Solve(std::vector<FVMeshNode*>* Mesh, std::vector<fv_prec_3>* phi_param) {
	(*Mesh)[0]->getPhiName();
	
	std::vector<cd_prec> Sources; Sources.resize((*Mesh).size(), 0.0);
	Matrix mat((*Mesh).size(), (*Mesh).size());
	for (int i = 0; i < (*Mesh).size(); i++) {
		mat(i, i) = (*Mesh)[i]->getCoef((*Mesh)[i]->getVFPointers());
		Sources[i] = (*Mesh)[i]->getSource();
		for (int j = 0; j < (*Mesh).size(); j++) {
			for (auto ptr : (*Mesh)[i]->getNbPointers()) {
				if (j == (*ptr).m_id) {
					mat(i, j) = (*Mesh)[i]->getCoef(ptr);
				}
			}
		}
	}
	//std::cout << mat;
	//std::cout << "\n";

	std::vector<cd_prec> tempVec; tempVec.resize((*Mesh).size(), 0.0);
	int lvl = 0;
	for (int i = 1; i < (*Mesh).size(); i++) {
		for (int j = 0; j < i; j++) {
			if (mat(i, j) != 0.0) {
				//if(mat(i - 1-j, j) != 0.0){
					for (int k = 0;k < (*Mesh).size(); k++) {
						tempVec[k] = mat(i - 1 - lvl+j, k)*mat(i, j) / mat(i - 1 - lvl+j, j);
					}
					for (int k = 0;k < (*Mesh).size(); k++) {
						mat(i, k) -= tempVec[k];
					}
				//}
			}
			//std::cout << "\n";
			//std::cout << mat;
			//std::cout << "\n";
		}
		lvl++;
	}
	std::cout << mat;
	std::cout << "\n";
	for (int i = 0; i < (*Mesh).size(); i++) {
		std::cout << Sources[i] << " ";
	}
	std::cout << "\n";

	std::vector<cd_prec> Results; Results.resize((*Mesh).size(), 0.0);
	for (int i = (*Mesh).size() - 1; i > 0; i--) {
		for (int j = (*Mesh).size() - 1; j > i; j--) {
			if (i == (*Mesh).size() - 1) {
				Results[i] = Sources[j] / mat(i, j);
			}
			else {
				Results[i] = (Sources[j] - Results[i + 1] * mat(i + 1, j)) / mat(i, j);
			}
		}
	}
	for (int i = 0; i < (*Mesh).size(); i++) {
		std::cout << Results[i] << " ";
	}
	std::cout << "\n\n";


	switch((*Mesh)[0]->getPhiName()){
	case(PHI_PARAMS::UX):
		break;
	case(PHI_PARAMS::UY):
		break;
	default:
		break;
	}

}


void FVM_CD::SIMPLE_ALGO() {
	std::vector<fv_prec> PApr;
	PApr.resize(FVM.FiniteVolumes.size(), 0.0);
	std::vector<fv_prec_3> VelApr;
	VelApr.resize(FVM.FiniteVolumes.size(), { 0.0,0.0,0.0 });
	std::vector<fv_prec> Pcor;
	Pcor.resize(FVM.FiniteVolumes.size(), 0.0);

	std::vector<Ch_<fv_prec>> Eps;

	std::vector<FVMeshNode*> FVMesh_VxA;
	std::vector<FVMeshNode*> FVMesh_VyA;
	std::vector<FVMeshNode*> FVMesh_PC;
	Ch_<fv_prec> sum = 0;
	for (auto e : Eps) { sum = sum + e; }

	for (auto i : FVM.FiniteVolumes) {
		PApr[i->m_id] = i->m_pressure.val;
	}
	while (sum.dval - sum.val > 0.001) {


		// ae ue* = SUM(anb unb*) + bu + dy(PP*-PE*)
		for (int i = 0;i < FVM.FiniteVolumeMesh.size();i++) {
			if (FVM.FiniteVolumeMesh[i]->FVPtr()->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
				FVMesh_VxA.push_back(new FVMeshNode(*(FVM.FiniteVolumeMesh[i]), PHI_PARAMS::UX));
			}
		}



		std::cout << "PHI_PARAMS::UX\n";
		Solve(&FVMesh_VxA, &VelApr);

		// an vn* = SUM(anb vnb*) + bv + dx(PP*-PN*)
		for (int i = 0;i < FVM.FiniteVolumeMesh.size();i++) {
			if (FVM.FiniteVolumeMesh[i]->FVPtr()->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
				FVMesh_VyA.push_back(new FVMeshNode(*(FVM.FiniteVolumeMesh[i]), PHI_PARAMS::UY));
			}
		}



		std::cout << "PHI_PARAMS::UY\n";
		Solve(&FVMesh_VyA, &VelApr);

		// ap PP' = SUM(anb PNB') + bp* , bp* = -{(Rhop-Rhop0)/dt dVp + SUM(Fnb*)} 
		for (int i = 0;i < FVM.FiniteVolumeMesh.size();i++) {
			if (FVM.FiniteVolumeMesh[i]->FVPtr()->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
				FVMesh_PC.push_back(new FVMeshNode(*(FVM.FiniteVolumeMesh[i]), PHI_PARAMS::P));
			}
		}



		std::cout << "PHI_PARAMS::P\n";
		Solve(&FVMesh_PC, &Pcor);
		/*
		// ue = ue* +de(PP'- PE')
		for (auto i : FVM.FiniteVolumes) {
			i->m_velocity.dval.x = VelApr[i->m_id].x + de * (Pcor[i->m_id] - Pcor[FVM.FiniteVolumeMesh[i->m_id].E]);
		}

		// vn = vn* +dn(PP'- PN')
		for (auto i : FVM.FiniteVolumes) {
			i->m_velocity.dval.y = VelApr[i->m_id].y + dn * (Pcor[i->m_id] - Pcor[FVM.FiniteVolumeMesh[i->m_id].N]);
		}
		*/
		// P = P* + P'
		for (auto i : FVM.FiniteVolumes) {
			if (i->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
				i->m_pressure.dval = Pcor[i->m_id] + PApr[i->m_id];
			}
		}

		// P* = P
		for (auto i : FVM.FiniteVolumes) {
			if (i->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
				PApr[i->m_id] = i->m_pressure.dval;
			}
		}


	}
	for (auto i : FVMesh_VxA) { delete i; }
	FVMesh_VxA.resize(0);
	for (auto i : FVMesh_VyA) { delete i; }
	FVMesh_VyA.resize(0);
	for (auto i : FVMesh_PC) { delete i; }
	FVMesh_PC.resize(0);



	for (auto i : FVM.FiniteVolumes) {
		if (i->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
			i->m_pressure.val = i->m_pressure.dval;
			i->m_velocity.val = i->m_velocity.dval;
		}
	}




		sum = 0;
		for (auto e : Eps) { sum = sum + e; }
	}

}


void FVM_CD::timeUpdate(cd_prec dt){
	if (m_options.firstCycle == false) {
		if (computationalDomain_mode == MODE_CD::CD_DEBUG) {
			std::cout << "SPH_CD::timeIntegration\n";
		}
		// We know P* field
		std::vector<fv_prec> PApr;
		PApr.resize(FVM.FiniteVolumes.size(), 0.0);
		std::vector<fv_prec_3> VelApr;
		VelApr.resize(FVM.FiniteVolumes.size(), { 0.0,0.0,0.0 });
		std::vector<fv_prec> Pcor;
		Pcor.resize(FVM.FiniteVolumes.size(), 0.0);

		std::vector<Ch_<fv_prec>> Eps;
		if (m_options.Steady) {
			
			if (m_options.mesh_s == STRUCTED) {
				std::vector<FVMeshNode*> FVMesh_VxA;
				std::vector<FVMeshNode*> FVMesh_VyA;
				std::vector<FVMeshNode*> FVMesh_PC;
				Ch_<fv_prec> sum = 0;
				for (auto e : Eps) { sum = sum + e; }

				for (auto i : FVM.FiniteVolumes) {
					PApr[i->m_id] = i->m_pressure.val;
				}

				//while(sum.dval - sum.val < 0.001)
				{
					// ae ue* = SUM(anb unb*) + bu + dy(PP*-PE*)
					for (int i = 0;i < FVM.FiniteVolumeMesh.size();i++) {
						if (FVM.FiniteVolumeMesh[i]->FVPtr()->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
							FVMesh_VxA.push_back(new FVMeshNode(*(FVM.FiniteVolumeMesh[i]), PHI_PARAMS::UX));
						}
					}
					for (auto i : FVMesh_VxA) {
						i->calcCoefsAndSources();
					}
					std::cout << "PHI_PARAMS::UX\n";
					Solve(&FVMesh_VxA, &VelApr);

					// an vn* = SUM(anb vnb*) + bv + dx(PP*-PN*)
					for (int i = 0;i < FVM.FiniteVolumeMesh.size();i++) {
						if (FVM.FiniteVolumeMesh[i]->FVPtr()->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
							FVMesh_VyA.push_back(new FVMeshNode(*(FVM.FiniteVolumeMesh[i]), PHI_PARAMS::UY));
						}
					}
					for (auto i : FVMesh_VyA) {
						i->calcCoefsAndSources();
					}
					std::cout << "PHI_PARAMS::UY\n";
					Solve(&FVMesh_VyA, &VelApr);

					// ap PP' = SUM(anb PNB') + bp* , bp* = -{(Rhop-Rhop0)/dt dVp + SUM(Fnb*)} 
					for (int i = 0;i < FVM.FiniteVolumeMesh.size();i++) {
						if (FVM.FiniteVolumeMesh[i]->FVPtr()->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
							FVMesh_PC.push_back(new FVMeshNode(*(FVM.FiniteVolumeMesh[i]), PHI_PARAMS::P));
						}
					}
					for (auto i : FVMesh_PC) {
						i->calcCoefsAndSources();
					}
					std::cout << "PHI_PARAMS::P\n";
					Solve(&FVMesh_PC, &Pcor);
					/*
					// ue = ue* +de(PP'- PE')
					for (auto i : FVM.FiniteVolumes) {
						i->m_velocity.dval.x = VelApr[i->m_id].x + de * (Pcor[i->m_id] - Pcor[FVM.FiniteVolumeMesh[i->m_id].E]);
					}

					// vn = vn* +dn(PP'- PN')
					for (auto i : FVM.FiniteVolumes) {
						i->m_velocity.dval.y = VelApr[i->m_id].y + dn * (Pcor[i->m_id] - Pcor[FVM.FiniteVolumeMesh[i->m_id].N]);
					}
					*/
					// P = P* + P'
					for (auto i : FVM.FiniteVolumes) {
						if (i->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
							i->m_pressure.dval = Pcor[i->m_id] + PApr[i->m_id];
						}
					}

					// P* = P
					for (auto i : FVM.FiniteVolumes) {
						if (i->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
							PApr[i->m_id] = i->m_pressure.dval;
						}
					}


				}
				for (auto i : FVMesh_VxA) { delete i; }
				FVMesh_VxA.resize(0);
				for (auto i : FVMesh_VyA) { delete i; }
				FVMesh_VyA.resize(0);
				for (auto i : FVMesh_PC) { delete i; }
				FVMesh_PC.resize(0);



				for (auto i : FVM.FiniteVolumes) {
					if (i->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
						i->m_pressure.val = i->m_pressure.dval;
						i->m_velocity.val = i->m_velocity.dval;
					}
				}
			}
		}
		else {

		}




		/*
		switch (m_options.timeIntegrationScheme) {
		case(EXPLICIT):
			////Полностью явная схема
		break;
		case(IMPLICIT):
			////Полностью неявная схема
		break;
		case(SEMI_IMPICIT):
			////Полунеявная схема
		break;
		case(NONE):
			break;
		default:
			break;
		}
		*/

	}
	else {
	}
}



inline glm::vec3 FVM_CD::ColoringGradient(float maxVal, float minVal, float currentVal) {

	float Red = 0.f;
	float Green = 0.f;
	float Blue = 0.f;
	float aveVal = (minVal + maxVal) / 2.f;
	if ((currentVal >= minVal) and (currentVal < aveVal)) {
		Red = 0.f;
		Green = (1.f - 0.f) / (aveVal - minVal)*(currentVal - minVal) + 0.f;
		Blue = 1.f - (1.f - 0.f) / (aveVal - minVal)*(currentVal - minVal);
	}
	else if ((currentVal >= aveVal) and (currentVal < maxVal)) {
		Red = (1.f - 0.f) / (maxVal - aveVal)*(currentVal - aveVal) + 0.f;
		Green = 1.f - (1.f - 0.f) / (maxVal - aveVal)*(currentVal - aveVal);
		Blue = 0.f;
	}
	if (currentVal == maxVal) {
		Red = 1.f;
		Green = 0.f;
		Blue = 0.f;
	}
	return glm::vec3(Red, Green, Blue);
}
void FVM_CD::Coloring() {
	if (computationalDomain_mode == MODE_CD::CD_DEBUG) {
		std::cout << "FVM_CD::Coloring\n";
	}
	float currentVal;
	std::string ColoringParameter = getColorParam();
	ColoringUnknow = false;
	if (ColoringParameter == "Vx") {
		changeColorParamTo("Vx");
		float maxVal = FVM.FiniteVolumes[0]->m_velocity.val.x;
		float minVal = FVM.FiniteVolumes[0]->m_velocity.val.x;
		for (auto*& i : FVM.FiniteVolumes) {
			//if (ColoringCondition(*i)) {
				if (i->m_velocity.val.x < minVal) { minVal = i->m_velocity.val.x; }
				if (i->m_velocity.val.x > maxVal) { maxVal = i->m_velocity.val.x; }
			//}
		}
		m_maxVal = maxVal;
		m_minVal = minVal;

		for (auto*& i : FVM.FiniteVolumes) {

			//if (ColoringCondition(*i)) {
				currentVal = i->m_velocity.val.x;
				i->m_color = ColoringGradient(m_maxVal, m_minVal, currentVal);
				continue;
			//}
			//i->m_color = glm::vec3(0.9f, 0.9f, 0.9f);

		}
	}
	else if (ColoringParameter == "Vy") {
		ColoringUnknow = false;
		changeColorParamTo("Vy");
		float maxVal = FVM.FiniteVolumes[0]->m_velocity.val.y;
		float minVal = FVM.FiniteVolumes[0]->m_velocity.val.y;
		for (auto*& i : FVM.FiniteVolumes) {
			//if (ColoringCondition(*i)) {
				if (i->m_velocity.val.y < minVal) { minVal = i->m_velocity.val.y; }
				if (i->m_velocity.val.y > maxVal) { maxVal = i->m_velocity.val.y; }
			//}
		}
		m_maxVal = maxVal;
		m_minVal = minVal;

		for (auto*& i : FVM.FiniteVolumes) {

			//if (ColoringCondition(*i)) {
				currentVal = i->m_velocity.val.y;
				i->m_color = ColoringGradient(m_maxVal, m_minVal, currentVal);
				continue;
			//}
			//i->m_color = glm::vec3(0.9f, 0.9f, 0.9f);

		}
	}
	else if (ColoringParameter == "abs(V)") {
		ColoringUnknow = false;
		changeColorParamTo("abs(V)");
		float maxVal = sqrt(pow(FVM.FiniteVolumes[0]->m_velocity.val.x, 2) + pow(FVM.FiniteVolumes[0]->m_velocity.val.y, 2) + pow(FVM.FiniteVolumes[0]->m_velocity.val.z, 2));
		float minVal = sqrt(pow(FVM.FiniteVolumes[0]->m_velocity.val.x, 2) + pow(FVM.FiniteVolumes[0]->m_velocity.val.y, 2) + pow(FVM.FiniteVolumes[0]->m_velocity.val.z, 2));
		for (auto*& i : FVM.FiniteVolumes) {
			//if (ColoringCondition(*i)) {
				if (sqrt(pow(i->m_velocity.val.x, 2) + pow(i->m_velocity.val.y, 2) + pow(i->m_velocity.val.z, 2)) < minVal) { minVal = sqrt(pow(i->m_velocity.val.x, 2) + pow(i->m_velocity.val.y, 2) + pow(i->m_velocity.val.z, 2)); }
				if (sqrt(pow(i->m_velocity.val.x, 2) + pow(i->m_velocity.val.y, 2) + pow(i->m_velocity.val.z, 2)) > maxVal) { maxVal = sqrt(pow(i->m_velocity.val.x, 2) + pow(i->m_velocity.val.y, 2) + pow(i->m_velocity.val.z, 2)); }
			//}
		}
		m_maxVal = maxVal;
		m_minVal = minVal;

		for (auto*& i : FVM.FiniteVolumes) {
			//if (ColoringCondition(*i)) {
				currentVal = sqrt(pow(i->m_velocity.val.x, 2) + pow(i->m_velocity.val.y, 2) + pow(i->m_velocity.val.z, 2));
				i->m_color = ColoringGradient(maxVal, minVal, currentVal);
				continue;
			//}
			//i->m_color = glm::vec3(0.9f, 0.9f, 0.9f);
		}
	}
	else if (ColoringParameter == "d(Vx)/dt") {
		ColoringUnknow = false;
		changeColorParamTo("d(Vx)/dt");
		float maxVal = FVM.FiniteVolumes[0]->m_velocity.dval.x;
		float minVal = FVM.FiniteVolumes[0]->m_velocity.dval.x;
		for (auto*& i : FVM.FiniteVolumes) {
			//if (ColoringCondition(*i)) {
				if (i->m_velocity.dval.x < minVal) { minVal = i->m_velocity.dval.x; }
				if (i->m_velocity.dval.x > maxVal) { maxVal = i->m_velocity.dval.x; }
			//}
		}
		m_maxVal = maxVal;
		m_minVal = minVal;

		for (auto*& i : FVM.FiniteVolumes) {
			//if (ColoringCondition(*i)) {
				currentVal = i->m_velocity.dval.x;
				i->m_color = ColoringGradient(maxVal, minVal, currentVal);
				continue;
			//}
			//i->m_color = glm::vec3(0.9f, 0.9f, 0.9f);
		}
	}
	else if (ColoringParameter == "d(Vy)/dt") {
		ColoringUnknow = false;
		changeColorParamTo("d(Vy)/dt");
		float maxVal = FVM.FiniteVolumes[0]->m_velocity.dval.y;
		float minVal = FVM.FiniteVolumes[0]->m_velocity.dval.y;
		for (auto*& i : FVM.FiniteVolumes) {
			//if (ColoringCondition(*i)) {
				if (i->m_velocity.dval.y < minVal) { minVal = i->m_velocity.dval.y; }
				if (i->m_velocity.dval.y > maxVal) { maxVal = i->m_velocity.dval.y; }
			//}
		}
		m_maxVal = maxVal;
		m_minVal = minVal;

		for (auto*& i : FVM.FiniteVolumes) {
			//if (ColoringCondition(*i)) {
				currentVal = i->m_velocity.dval.y;
				i->m_color = ColoringGradient(maxVal, minVal, currentVal);
				continue;
			//}
			//i->m_color = glm::vec3(0.9f, 0.9f, 0.9f);
		}
	}
	else if (ColoringParameter == "dV/dt") {
		ColoringUnknow = false;
		changeColorParamTo("dV/dt");
		float maxVal = sqrt(pow(FVM.FiniteVolumes[0]->m_velocity.dval.x, 2) + pow(FVM.FiniteVolumes[0]->m_velocity.dval.y, 2) + pow(FVM.FiniteVolumes[0]->m_velocity.dval.z, 2));
		float minVal = sqrt(pow(FVM.FiniteVolumes[0]->m_velocity.dval.x, 2) + pow(FVM.FiniteVolumes[0]->m_velocity.dval.y, 2) + pow(FVM.FiniteVolumes[0]->m_velocity.dval.z, 2));
		for (auto*& i : FVM.FiniteVolumes) {
			//if (ColoringCondition(*i)) {
				if (sqrt(pow(i->m_velocity.dval.x, 2) + pow(i->m_velocity.dval.y, 2) + pow(i->m_velocity.dval.z, 2)) < minVal) { minVal = sqrt(pow(i->m_velocity.dval.x, 2) + pow(i->m_velocity.dval.y, 2) + pow(i->m_velocity.dval.z, 2)); }
				if (sqrt(pow(i->m_velocity.dval.x, 2) + pow(i->m_velocity.dval.y, 2) + pow(i->m_velocity.dval.z, 2)) > maxVal) { maxVal = sqrt(pow(i->m_velocity.dval.x, 2) + pow(i->m_velocity.dval.y, 2) + pow(i->m_velocity.dval.z, 2)); }
			//}
		}
		m_maxVal = maxVal;
		m_minVal = minVal;

		for (auto*& i : FVM.FiniteVolumes) {
			//if (ColoringCondition(*i)) {
				currentVal = sqrt(pow(i->m_velocity.dval.x, 2) + pow(i->m_velocity.dval.y, 2) + pow(i->m_velocity.dval.z, 2));
				i->m_color = ColoringGradient(maxVal, minVal, currentVal);
				continue;
			//}
			//i->m_color = glm::vec3(0.9f, 0.9f, 0.9f);
		}
	}
	else if (ColoringParameter == "Density") {
		ColoringUnknow = false;
		changeColorParamTo("Density");
		float maxVal = FVM.FiniteVolumes[0]->m_density.val;
		float minVal = FVM.FiniteVolumes[0]->m_density.val;
		for (auto*& i : FVM.FiniteVolumes) {
			//if (ColoringCondition(*i)) {
				if (i->m_density.val < minVal) { minVal = i->m_density.val; }
				if (i->m_density.val > maxVal) { maxVal = i->m_density.val; }
			//}
		}
		m_maxVal = maxVal;
		m_minVal = minVal;

		for (auto*& i : FVM.FiniteVolumes) {
			//if (ColoringCondition(*i)) {
				currentVal = i->m_density.val;
				i->m_color = ColoringGradient(m_maxVal, m_minVal, currentVal);
				continue;
			//}
			//i->m_color = glm::vec3(0.9f, 0.9f, 0.9f);
		}
	}
	else if (ColoringParameter == "d(Density)/dt") {
		ColoringUnknow = false;
		changeColorParamTo("d(Density)/dt");
		float maxVal = FVM.FiniteVolumes[0]->m_density.dval;
		float minVal = FVM.FiniteVolumes[0]->m_density.dval;
		for (auto*& i : FVM.FiniteVolumes) {
			//if (ColoringCondition(*i)) {
				if (i->m_density.dval < minVal) { minVal = i->m_density.dval; }
				if (i->m_density.dval > maxVal) { maxVal = i->m_density.dval; }
			//}
		}
		m_maxVal = maxVal;
		m_minVal = minVal;

		for (auto*& i : FVM.FiniteVolumes) {
			//if (ColoringCondition(*i)) {
				currentVal = i->m_density.dval;
				i->m_color = ColoringGradient(maxVal, minVal, currentVal);
				continue;
			//}
			//i->m_color = glm::vec3(0.9f, 0.9f, 0.9f);
		}
	}
	else if (ColoringParameter == "P") {
		ColoringUnknow = false;
		changeColorParamTo("P");
		float maxVal = FVM.FiniteVolumes[0]->m_pressure.val;
		float minVal = FVM.FiniteVolumes[0]->m_pressure.val;
		for (auto*& i : FVM.FiniteVolumes) {
			//if (ColoringCondition(*i)) {
				if (i->m_pressure.val < minVal) { minVal = i->m_pressure.val; }
				if (i->m_pressure.val > maxVal) { maxVal = i->m_pressure.val; }
			//}
		}
		m_maxVal = maxVal;
		m_minVal = minVal;

		for (auto*& i : FVM.FiniteVolumes) {
			//if (ColoringCondition(*i)) {
				currentVal = i->m_pressure.val;
				i->m_color = ColoringGradient(maxVal, minVal, currentVal);
				continue;
			//}
			//i->m_color = glm::vec3(0.9f, 0.9f, 0.9f);
		}
	}
	else if (ColoringParameter == "dP/dt") {
		ColoringUnknow = false;
		changeColorParamTo("dP/dt");
		float maxVal = FVM.FiniteVolumes[0]->m_pressure.dval;
		float minVal = FVM.FiniteVolumes[0]->m_pressure.dval;
		for (auto*& i : FVM.FiniteVolumes) {
			//if (ColoringCondition(*i)) {
			if (i->m_pressure.dval < minVal) { minVal = i->m_pressure.dval; }
			if (i->m_pressure.dval > maxVal) { maxVal = i->m_pressure.dval; }
			//}
		}
		m_maxVal = maxVal;
		m_minVal = minVal;

		for (auto*& i : FVM.FiniteVolumes) {
			//if (ColoringCondition(*i)) {
			currentVal = i->m_pressure.dval;
			i->m_color = ColoringGradient(maxVal, minVal, currentVal);
			continue;
			//}
			//i->m_color = glm::vec3(0.9f, 0.9f, 0.9f);
		}
	}
	else if (ColoringParameter == "SurfaceNormal_x") {
		changeColorParamTo("mass");
		float maxVal = FVM.FiniteVolumes[0]->m_normalToSurface.x;
		float minVal = FVM.FiniteVolumes[0]->m_normalToSurface.x;
		for (auto*& i : FVM.FiniteVolumes) {
			//if (ColoringCondition(*i)) {
				if (i->m_normalToSurface.x < minVal) { minVal = i->m_normalToSurface.x; }
				if (i->m_normalToSurface.x > maxVal) { maxVal = i->m_normalToSurface.x; }
			//}
		}
		float aveVal = (minVal + maxVal) / 2;

		m_maxVal = maxVal;
		m_minVal = minVal;

		for (auto*& i : FVM.FiniteVolumes) {
			//if (ColoringCondition(*i)) {
				currentVal = i->m_normalToSurface.x;
				i->m_color = ColoringGradient(maxVal, minVal, currentVal);
				continue;
			//}
			//i->m_color = glm::vec3(0.9f, 0.9f, 0.9f);
		}
	}
	else if (ColoringParameter == "SurfaceNormal_y") {
		changeColorParamTo("SurfaceNormal_y");
		float maxVal = FVM.FiniteVolumes[0]->m_normalToSurface.y;
		float minVal = FVM.FiniteVolumes[0]->m_normalToSurface.y;
		for (auto*& i : FVM.FiniteVolumes) {
			//if (ColoringCondition(*i)) {
				if (i->m_normalToSurface.y < minVal) { minVal = i->m_normalToSurface.y; }
				if (i->m_normalToSurface.y > maxVal) { maxVal = i->m_normalToSurface.y; }
			//}
		}
		float aveVal = (minVal + maxVal) / 2;

		m_maxVal = maxVal;
		m_minVal = minVal;

		for (auto*& i : FVM.FiniteVolumes) {
			//if (ColoringCondition(*i)) {
				currentVal = i->m_normalToSurface.y;
				i->m_color = ColoringGradient(maxVal, minVal, currentVal);
				continue;
			//}
			//i->m_color = glm::vec3(0.9f, 0.9f, 0.9f);
		}
	}
	else {
		ColoringUnknow = true;
		changeColorParamTo("id");
		float maxVal = FVM.FiniteVolumes[0]->m_id;
		float minVal = FVM.FiniteVolumes[0]->m_id;
		for (auto*& i : FVM.FiniteVolumes) {
			//if (ColoringCondition(*i)) {
				if (i->m_id < minVal) { minVal = i->m_id; }
				if (i->m_id > maxVal) { maxVal = i->m_id; }
			//}
		}
		m_maxVal = maxVal;
		m_minVal = minVal;

		for (auto*& i : FVM.FiniteVolumes) {
			//if (ColoringCondition(*i)) {
				currentVal = i->m_id;
				i->m_color = ColoringGradient(maxVal, minVal, currentVal);
				continue;
			//}
			//i->m_color = glm::vec3(0.9f, 0.9f, 0.9f);
		}
	}
}



FVM_CD::~FVM_CD() {}