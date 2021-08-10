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

	// œÂÂÌÂÒÚË ‚ ÙÛÌÍˆË˛ EquationsInitialization
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
	int meshCount = 0;
	//fv_prec pressure = 10E05;
	fv_prec pressure = 0;
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
							//std::cout << "adding FiniteVolumeMesh\n";
							FVM.FiniteVolumeMesh[count]->addNeighbour(FVM.FiniteVolumes[fvc]);
							meshCount++;
							FVM.FiniteVolumeMesh[fvc]->addNeighbour(FVM.FiniteVolumes[count]);
							meshCount++;
						}
					}
				}
				count++;
			}
			m_options.nrOfFV[FVT_DEFAULT] = count;
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
							//std::cout << "adding FiniteVolumeMesh\n";
							FVM.FiniteVolumeMesh[count]->addNeighbour(FVM.FiniteVolumes[fvc]);
							meshCount++;
							FVM.FiniteVolumeMesh[fvc]->addNeighbour(FVM.FiniteVolumes[count]);
							meshCount++;
						}
					}
					count++;
				}
			}
			m_options.nrOfFV[FVT_DEFAULT] = count;
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
									//std::cout << "adding FiniteVolumeMesh\n";
									FVM.FiniteVolumeMesh[count]->addNeighbour(FVM.FiniteVolumes[fvc]);
									meshCount++;
									FVM.FiniteVolumeMesh[fvc]->addNeighbour(FVM.FiniteVolumes[count]);
									meshCount++;
								}
							}
						}
						count++;
					}
				}
			}
			m_options.nrOfFV[FVT_DEFAULT] = count;
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
	std::cout << count << " finite volumes added\n";
	std::cout << meshCount << " neighbours in mesh\n";
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
			//«Ì‡ˇ ÚË ÚÓ˜ÍË, ÏÓÊÌÓ ÓÔÂ‰ÂÎËÚ¸ ÙÛÌÍˆË˛ ÔÎÓÒÍÓÒÚË
			PlaneCoefA = (Plane2.y - Plane1.y)*(Plane3.z - Plane1.z) - (Plane3.y - Plane1.y)*(Plane2.z - Plane1.z);
			PlaneCoefB = -(Plane2.x - Plane1.x)*(Plane3.z - Plane1.z) - (Plane3.x - Plane1.x)*(Plane2.z - Plane1.z);
			PlaneCoefC = (Plane2.x - Plane1.x)*(Plane3.y - Plane1.y) - (Plane3.x - Plane1.x)*(Plane2.y - Plane1.y);
			PlaneCoefD = -Plane1.x*PlaneCoefA - Plane1.y*PlaneCoefB - Plane1.z*PlaneCoefC;
			//std::cout << PlaneCoefA << "*x + " << PlaneCoefB << "*y + " << PlaneCoefC << "*z + " << PlaneCoefD << " =  0\n";
			//ƒ‡ÎÂÂ ÓÔÂ‰ÂÎËÏ, ÎËÌË˛ ÔÂÂÒÂ˜ÂÌËÂ ÔÓÎÛ˜ÂÌÌÓÈ ÔÎÓÒÍÓÒÚË Ò ÔÎÓÒÍÓÒÚ¸˛ OXY  (z=0):
			// x*PlaneCoefA + y*PlaneCoefB + PlaneCoefD = 0
			//ŒÔÂ‰ÂÎËÏ ÔÂÂÒÂ˜ÂÌËÂ ˝ÚÓÈ ÎËÌËË Ò Í‡Ê‰ÓÈ ÎËÌËÂÈ, ËÁ ÍÓÚÓ˚ı ÒÓÒÚÓËÚ ÔÓÎË„ÓÌ
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
	int meshCount = 0;
	int BoundaryVoluemCountOld = FVM.FiniteVolumes.size();
	int BoundaryVoluemCount = FVM.FiniteVolumes.size();
	//fv_prec pressure = 10E05;
	fv_prec pressure = 0;
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
						//«Ì‡ˇ ÚË ÚÓ˜ÍË, ÏÓÊÌÓ ÓÔÂ‰ÂÎËÚ¸ ÙÛÌÍˆË˛ ÔÎÓÒÍÓÒÚË
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
						FVM.FiniteVolumeMesh.push_back(new FVMeshNodeEmpty(FVM.FiniteVolumes[BoundaryVoluemCount]));
						for (int fvc = 0;fvc < BoundaryVoluemCount;fvc++) {
							for (auto l : FVM.FiniteVolumes[fvc]->m_s.m_lineArray) {
								if (FVM.FiniteVolumes[BoundaryVoluemCount]->m_l == l) {
									//std::cout << "adding FiniteVolumeMesh\n";
									FVM.FiniteVolumeMesh[BoundaryVoluemCount]->addNeighbour(FVM.FiniteVolumes[fvc]);
									meshCount++;
									FVM.FiniteVolumeMesh[fvc]->addNeighbour(FVM.FiniteVolumes[BoundaryVoluemCount]);
									meshCount++;
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

	m_options.nrOfFV[FVT_BOUNDARY] = BoundaryVoluemCount - BoundaryVoluemCountOld;
	std::cout << BoundaryVoluemCount - BoundaryVoluemCountOld << " finite volumes added\n";
	std::cout << meshCount << " neighbours in mesh\n";
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
			//rotation =  —ÀŒ∆ÕŒ ≈—À» œŒ-œ–¿¬»À‹ÕŒÃ”
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

			assignNewModels(ModelId + FVM.FiniteVolumes.size());


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



//void FVM_CD::Solve(std::vector<FVMeshNode*>* Mesh, std::vector<fv_prec>* phi_param) {
//
//}

void FVM_CD::SolveWithFirstZero(std::vector<FVMeshNode*>* Mesh, std::vector<fv_prec>* phi_param) {
	(*Mesh)[0]->getPhiName();
	phi_param->resize(0);
	bool print_mat = false;
	bool print_source = false;
	bool print_res = false;
	int m_size = 0;
	for (auto fvm : *Mesh) {
		if (fvm->fv.fv->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
			m_size++;
		}
	}

	std::vector<cd_prec> Sources; Sources.resize(m_size, 0.0);
	Matrix mat(m_size, m_size);
	for (int i = 0; i < m_size; i++) {
		for (int j = 0; j < m_size; j++) {
			//for (auto ptr : (*Mesh)[i]->nb_fvs) {
			//	if (j == ptr.innerid) {
			//		mat(i, j) = -(*Mesh)[i]->getCoef(ptr.fv);
			//	}
			//}
			for (auto ptr : (*Mesh)[i]->getNbPointers()) {
				if (j == (*ptr).m_id) {
					mat(i, j) = -(*Mesh)[i]->getCoef(ptr);
				}
			}
		}
		mat(i, i) = (*Mesh)[i]->El().fv_coef;
		Sources[i] = (*Mesh)[i]->getSource();
	}
	
	fv_prec def_press = 0.0;
	for (int i = 0;i < m_size;i++) {
		Sources[i] -= mat(i, 0) * def_press;
		mat(i, 0) = 0.0;
	}


	if (print_mat) {
		std::cout << mat;
		std::cout << "\n";
	}
	if (print_source) {
		std::cout << "Source: ";
		for (int i = 0; i < m_size; i++) {
			std::cout << Sources[i] << " ";
		}
		std::cout << "\n";
	}
	Matrix tempMat(mat);
	std::vector<cd_prec> tempVec1; tempVec1.resize(m_size, 0.0);
	cd_prec tempSourceChange;
	for (int j = 0; j < m_size; j++) { //  ÓÎÎÓÌ‡
		for (int i = j + 1; i < m_size; i++) { // –ˇ‰
			//fv_prec Koef = tempMat(i, j) / tempMat(j, j);
			//std::cout << tempMat(i, j) << "/" << tempMat(j, j) << "\n";
			if (tempMat(i, j) != 0) {
				for (int k = 0; k < m_size; k++) { //  ÓÎÎÓÌ‡
					//std::cout << tempMat(i, k) - tempMat(j, k)*tempMat(i, j) / tempMat(j, j) << " ";
					tempVec1[k] = tempMat(i, k) - tempMat(j, k)*tempMat(i, j) / tempMat(j, j);
				}
				tempSourceChange = Sources[i] - Sources[j] * tempMat(i, j) / tempMat(j, j);
				for (int k = 0; k < m_size; k++) { //  ÓÎÎÓÌ‡
					tempMat(i, k) = tempVec1[k];
				}
				Sources[i] = tempSourceChange;
				//std::cout << "\n";
				//std::cout << tempMat;
				//std::cout << "\n";
			}
		}
	}
	mat = tempMat;
	if (print_mat) {
		std::cout << mat;
		std::cout << "\n";
	}
	if (print_source) {
		std::cout << "Source: ";
		for (int i = 0; i < m_size; i++) {
			std::cout << Sources[i] << " ";
		}
		std::cout << "\n";
	}

	std::vector<cd_prec> Results; Results.resize(m_size, 0.0);
	for (int i = m_size - 1; i > 0; i--) {
		Results[i] = Sources[i];
		if (i == m_size - 1) { }
		for (int j = m_size - 1; j > i; j--) {
			if (i == m_size - 1) { }
			else { Results[i] -= Results[j] * mat(i, j); }
		}
		Results[i] /= mat(i, i);

	}
	for (int i = 0; i < m_size; i++) {
		phi_param->push_back(Results[i]);
	}
	if (print_res) {
		std::cout << "Reaults: ";
		for (int i = 0; i < m_size; i++) {
			std::cout << Results[i] << " ";
		}
		std::cout << "\n\n";
	}


	switch ((*Mesh)[0]->getPhiName()) {
	case(PHI_PARAMS::UX):
		break;
	case(PHI_PARAMS::UY):
		break;
	default:
		break;
	}

}

Matrix FVM_CD::Solve(std::vector<FVMeshNode*>* Mesh, std::vector<fv_prec>* phi_param) {
	(*Mesh)[0]->getPhiName();
	phi_param->resize(0);
	bool print_mat = false;
	bool print_source = false;
	bool print_res = false;
	int m_size = 0;
	for (auto fvm : *Mesh) {
		if (fvm->fv.fv->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
			m_size++;
		}
	}


	std::vector<cd_prec> Sources; Sources.resize(m_size, 0.0);
	Matrix mat(m_size, m_size);
	for (int i = 0; i < m_size; i++) {
		for (int j = 0; j < m_size; j++) {
			//for (auto ptr : (*Mesh)[i]->nb_fvs) {
			//	if (j == ptr.innerid) {
			//		mat(i, j) = -(*Mesh)[i]->getCoef(ptr.fv);
			//	}
			//}
			for (auto ptr : (*Mesh)[i]->getNbPointers()) {
				if (j == (*ptr).m_id) {
					mat(i, j) = -(*Mesh)[i]->getCoef(ptr);
				}
			}
		}
		mat(i, i) = (*Mesh)[i]->El().fv_coef;
		Sources[i] = (*Mesh)[i]->getSource();
	}
	if(print_mat){
		std::cout << mat;
		std::cout << "\n";
	}
	if (print_source) {
		std::cout << "Source: ";
		for (int i = 0; i < m_size; i++) {
			std::cout << Sources[i] << " ";
		}
		std::cout << "\n";
	}
	Matrix RetValmat(mat);


	//std::vector<cd_prec> tempVec; tempVec.resize((*Mesh).size(), 0.0);
	//int lvl = 1;
	//for (int j = 0; j < lvl; j++) {
	//	for (int i = lvl; i < (*Mesh).size(); i++) {
	//		if (mat(i, j) != 0.0) {
	//			//if(mat(i - lvl+j, j) != 0.0){
	//				for (int k = 0;k < (*Mesh).size(); k++) {
	//					tempVec[k] = mat(j, k)*mat(i, j) / mat(j, j);
	//				}
	//				for (int k = 0;k < (*Mesh).size(); k++) {
	//					mat(i, k) -= tempVec[k];
	//				}
	//			//}
	//  		//std::cout << "\n";
	//  		//std::cout << mat;
	//  		//std::cout << "\n";
	//		}
	//	}
	//	lvl++;
	//}



	Matrix tempMat(mat);
	std::vector<cd_prec> tempVec1; tempVec1.resize(m_size, 0.0);
	cd_prec tempSourceChange;
	for (int j = 0; j < m_size; j++) { //  ÓÎÎÓÌ‡
		for (int i = j+1; i < m_size; i++) { // –ˇ‰
			//fv_prec Koef = tempMat(i, j) / tempMat(j, j);
			//std::cout << tempMat(i, j) << "/" << tempMat(j, j) << "\n";
			if(tempMat(i, j) != 0){
				for (int k = 0; k < m_size; k++) { //  ÓÎÎÓÌ‡
					//std::cout << tempMat(i, k) - tempMat(j, k)*tempMat(i, j) / tempMat(j, j) << " ";
					tempVec1[k] = tempMat(i, k) - tempMat(j, k)*tempMat(i, j) / tempMat(j, j);
				}
				tempSourceChange = Sources[i] - Sources[j] * tempMat(i, j) / tempMat(j, j);
				for (int k = 0; k < m_size; k++) { //  ÓÎÎÓÌ‡
					tempMat(i, k) = tempVec1[k];
				}
				Sources[i] = tempSourceChange;
				//std::cout << "\n";
				//std::cout << tempMat;
				//std::cout << "\n";
			}
		}
	}
	mat = tempMat;

	//std::cout << "\n";
	//std::cout << tempMat;
	//std::cout << "\n";

	fv_prec Koef;

	//std::vector<cd_prec> tempVec; tempVec.resize((*Mesh).size(), 0.0);
	//int lvl = 0;
	//for (int i = 1; i < (*Mesh).size(); i++) {
	//	for (int j = 0; j < i; j++) {
	//		if (mat(i, j) != 0.0) {
	//			if(mat(i - 1 - lvl + j, j) != 0.0){
	//				Koef = mat(i, j) / mat(i - 1 - lvl + j, j);
	//				for (int k = 0;k < (*Mesh).size(); k++) {
	//					tempVec[k] = mat(i - 1 - lvl+j, k)*mat(i, j) / mat(i - 1 - lvl + j, j);
	//					Sources[i] -= Sources[i - 1 - lvl + j] * mat(i, j) / mat(i - 1 - lvl + j, j);
	//				}
	//				for (int k = 0;k < (*Mesh).size(); k++) {
	//					mat(i, k) -= tempVec[k];
	//				}
	//			}
	//		}
	//		//std::cout << "\n";
			//std::cout << mat;
			//std::cout << "\n";
			//for (int i = 0; i < (*Mesh).size(); i++) {
			//	std::cout << Sources[i] << " ";
			//}
			//std::cout << "\n";
	//	}
	//	lvl++;
	//}

	if (print_mat) {
		std::cout << mat;
		std::cout << "\n";
	}
	if (print_source) {
		std::cout << "Source: ";
		for (int i = 0; i < m_size; i++) {
			std::cout << Sources[i] << " ";
		}
		std::cout << "\n";
	}

	std::vector<cd_prec> Results; Results.resize(m_size, 0.0);
	for (int i = m_size - 1; i >= 0; i--) {
		Results[i] = Sources[i];
		if (i == m_size - 1) { }
		for (int j = m_size - 1; j > i; j--) {
			if (i == m_size - 1) { }
			else { Results[i] -= Results[j] * mat(i, j); }
		}
		Results[i] /= mat(i, i);
	}
	for (int i = 0; i < m_size; i++) {
		phi_param->push_back(Results[i]);
	}
	if (print_res) {
		std::cout << "Results: ";
		for (int i = 0; i < m_size; i++) {
			std::cout << Results[i] << " ";
		}
		std::cout << "\n\n";
	}


	return RetValmat;

}


inline fv_prec vec_mod(fv_prec_3 vec) {
	return (sqrt(pow(vec.x, 2) + pow(vec.y, 2) + pow(vec.z, 2)));
}

inline fv_prec FVM_CD::PeFunc(fv_prec Pe) {
	switch (m_options.convectionDiffusionScheme) {
	case(CONV_DIFF_SCHEME::CENTRAL_DIFFERENSING): return 1.0 - 0.5*abs(Pe);
	case(CONV_DIFF_SCHEME::UPWIND): return 1.0;
	case(CONV_DIFF_SCHEME::HYBRID): return glm::max(glm::max(-1.0*Pe, (1-0.5*Pe)), 0.0);
	case(CONV_DIFF_SCHEME::POWER_LAW): return glm::max(0.0, pow(1.0 - 0.1*abs(Pe),5));
	case(CONV_DIFF_SCHEME::EXPONENTIAL): return Pe/ (exp(Pe) - 1);
	default: assert("FVM_CD::PeFunc" && 0);
	}
}


void FVM_CD::SIMPLE_ALGO(cd_prec dt) {
	std::cout << "FVM_CD::SIMPLE_ALGO(" << dt << ")\n";
	int AlgoIteration = 1;

	fv_prec AlphaP = 0.8;
	fv_prec AlphaU = 0.5;
	fv_prec AlphaV = 0.5;
	fv_prec AlphaW = 0.5;
	fv_prec residualScale = 1.0;
	fv_prec Eps = 0.1;
	fv_prec maxIteration = 1;


	std::vector<FiniteVolume*> UVolumes;
	std::vector<FVMeshNode*> FVMesh_U;
	std::cout << "Creating U mesh \n";
	//std::vector<FiniteVolumeMesh*> UVolumeMesh;

	int meshCount = 0;
	int UVc = 0;
	for (auto fvm : FVM.FiniteVolumeMesh) {
		if ((fvm->getPtrToSide(NB_E)->fv_boundary == nullptr) and(fvm->FVPtr()->fv_boundary == nullptr)) {
			fv_prec_3 p1(fvm->FVPtr()->m_s.getCentre() + glm::vec3(0.0, -m_options.average_dim_steps.y / 2.0, 0.0));
			fv_prec_3 p2(fvm->getPtrToSide(NB_E)->m_s.getCentre() + glm::vec3(0.0, -m_options.average_dim_steps.y / 2.0, 0.0));
			fv_prec_3 p3(fvm->FVPtr()->m_s.getCentre() + glm::vec3(0.0, m_options.average_dim_steps.y / 2.0, 0.0));
			fv_prec_3 p4(fvm->getPtrToSide(NB_E)->m_s.getCentre() + glm::vec3(0.0, m_options.average_dim_steps.y / 2.0, 0.0));
			if (fvm->getPtrToSide(NB_W)->fv_boundary != nullptr) {
				p1 =(fvm->FVPtr()->m_s.getCentre() + glm::vec3(-m_options.average_dim_steps.x / 2.0, -m_options.average_dim_steps.y / 2.0, 0.0));
				p3 =(fvm->FVPtr()->m_s.getCentre() + glm::vec3(-m_options.average_dim_steps.x / 2.0, m_options.average_dim_steps.y / 2.0, 0.0));
			}
			if (FVM.FiniteVolumeMesh[fvm->getPtrToSide(NB_E)->m_id]->getPtrToSide(NB_E)->fv_boundary != nullptr) {
				p2 =(fvm->getPtrToSide(NB_E)->m_s.getCentre() + glm::vec3(m_options.average_dim_steps.x / 2.0, -m_options.average_dim_steps.y / 2.0, 0.0));
				p4 =(fvm->getPtrToSide(NB_E)->m_s.getCentre() + glm::vec3(m_options.average_dim_steps.x / 2.0, m_options.average_dim_steps.y / 2.0, 0.0));
			}
			//std::cout << UVc << ": {" << p1.x << "," << p1.y << "}, {" << p2.x << "," << p2.y << "}, {" << p3.x << "," << p3.y << "}, {" << p4.x << "," << p4.y << "}\n";
			Surface tempSurface(Line(p1, p2), Line(p2, p4), Line(p4, p3)); tempSurface.addLine(Line(p3, p1));
			tempSurface.calcCentre();
			UVolumes.push_back(new FiniteVolume(UVc, FINITE_VOLUME_TYPE::FVT_DEFAULT, tempSurface, glm::vec3(0.0f, 0.f, 0.f), 0.0));
			FVMesh_U.push_back(new FVMeshNode(UVolumes[UVc], FPARAM::FU));

			for (int fvc = 0;fvc < UVc;fvc++) {
				//if (UVc == 4) {
				//	std::cout << "{" << UVolumes[0]->m_s.m_lineArray[2].m_p1.x << ", " << UVolumes[0]->m_s.m_lineArray[2].m_p1.y << "},{" << UVolumes[0]->m_s.m_lineArray[2].m_p2.x << ", " << UVolumes[0]->m_s.m_lineArray[2].m_p2.y << "}-{" << UVolumes[UVc]->m_s.m_lineArray[0].m_p1.x << ", " << UVolumes[UVc]->m_s.m_lineArray[0].m_p1.y << "},{" << UVolumes[UVc]->m_s.m_lineArray[0].m_p2.x << ", " << UVolumes[UVc]->m_s.m_lineArray[0].m_p2.y << "}" << (UVolumes[0]->m_s.m_lineArray[2] == UVolumes[UVc]->m_s.m_lineArray[0]) << "\n";
				//	std::cout << "" << UVolumes[0]->m_s.m_lineArray[2].m_p1.x << "-" << UVolumes[UVc]->m_s.m_lineArray[0].m_p2.x<< "= "  << UVolumes[0]->m_s.m_lineArray[2].m_p1.x - UVolumes[UVc]->m_s.m_lineArray[0].m_p2.x << " " << (UVolumes[0]->m_s.m_lineArray[2].m_p1.x == UVolumes[UVc]->m_s.m_lineArray[0].m_p2.x) << "\n";
				//	std::cout << "" << UVolumes[0]->m_s.m_lineArray[2].m_p1.y << "-" << UVolumes[UVc]->m_s.m_lineArray[0].m_p2.y << "= " << UVolumes[0]->m_s.m_lineArray[2].m_p1.y - UVolumes[UVc]->m_s.m_lineArray[0].m_p2.y << " " << (UVolumes[0]->m_s.m_lineArray[2].m_p1.y == UVolumes[UVc]->m_s.m_lineArray[0].m_p2.y) << "\n";
				//	std::cout << "" << UVolumes[0]->m_s.m_lineArray[2].m_p2.x << "-" << UVolumes[UVc]->m_s.m_lineArray[0].m_p1.x << "= " << UVolumes[0]->m_s.m_lineArray[2].m_p2.x - UVolumes[UVc]->m_s.m_lineArray[0].m_p1.x << " " << (UVolumes[0]->m_s.m_lineArray[2].m_p2.x == UVolumes[UVc]->m_s.m_lineArray[0].m_p1.x) << "\n";
				//	std::cout << "" << UVolumes[0]->m_s.m_lineArray[2].m_p2.y << "-" << UVolumes[UVc]->m_s.m_lineArray[0].m_p1.y << "= " << UVolumes[0]->m_s.m_lineArray[2].m_p2.y - UVolumes[UVc]->m_s.m_lineArray[0].m_p1.y << " "  << (UVolumes[0]->m_s.m_lineArray[2].m_p2.y == UVolumes[UVc]->m_s.m_lineArray[0].m_p1.y) << "\n";
				//
				//}
				if (UVolumes[UVc]->m_s.ConnectedWith(UVolumes[fvc]->m_s)) {
					//std::cout << "(" << UVolumes[UVc]->m_Centre.x << "," << UVolumes[UVc]->m_Centre.y << " - " << UVolumes[fvc]->m_Centre.x << "," << UVolumes[fvc]->m_Centre.y << ") ";
					//std::cout << "adding FiniteVolumeMesh \n";
					FVMesh_U[UVc]->addNeighbour(UVolumes[fvc]);
					meshCount++;
					FVMesh_U[fvc]->addNeighbour(UVolumes[UVc]);
					meshCount++;
				}
			}
			//std::cout << FVMesh_U[UVc]->nb_fvs.size() << "\n";
			UVc++;
		}
	}
	int Usize = UVc;
	for (auto fvm : FVM.FiniteVolumeMesh) {
		fv_prec_3 p1;
		fv_prec_3 p2;
		Line l;
		fv_prec_3 vel;
		fv_prec pressure;
		if ((fvm->FVPtr()->fv_boundary == nullptr)) {
			if ((fvm->getPtrToSide(NB_E)->fv_boundary == nullptr)) {
				if (fvm->getPtrToSide(NB_W)->fv_boundary != nullptr) {
					l = fvm->getPtrToSide(NB_W)->m_l;
					vel = fvm->getPtrToSide(NB_W)->m_velocity.val;
					pressure = fvm->getPtrToSide(NB_W)->m_pressure.val;
					//std::cout << UVc << ": {" << l.m_p1.x << "," << l.m_p1.y << "}, {" << l.m_p2.x << "," << l.m_p2.y << "}" << " NB_W ";
					UVolumes.push_back(new FiniteVolume(UVc, FINITE_VOLUME_TYPE::FVT_BOUNDARY, l, vel, pressure));
					UVolumes[UVc]->assignToBoundary(fvm->getPtrToSide(NB_W)->fv_boundary);
					FVMesh_U.push_back(new FVMeshNode(UVolumes[UVc], FPARAM::FU));
					for (int fvc = 0;fvc < UVc;fvc++) {
						for (auto l : UVolumes[fvc]->m_s.m_lineArray) {
							if (UVolumes[UVc]->m_l == l) {
								//std::cout << "adding FiniteVolumeMesh\n";
								FVMesh_U[UVc]->addNeighbour(UVolumes[fvc]);
								meshCount++;
								FVMesh_U[fvc]->addNeighbour(UVolumes[UVc]);
								meshCount++;
							}
						}
					}
					//std::cout << UVc << ": " << UVolumes[UVc]->fv_boundary << " NB_W ";
					UVc++;

				}
				if (fvm->getPtrToSide(NB_N)->fv_boundary != nullptr) {					
					p1 = fvm->getPtrToSide(NB_N)->m_Centre;
					p2 = FVM.FiniteVolumeMesh[fvm->getPtrToSide(NB_E)->m_id]->getPtrToSide(NB_N)->m_Centre;
					
					vel = 0.5f*(fvm->getPtrToSide(NB_N)->m_velocity.val + FVM.FiniteVolumeMesh[fvm->getPtrToSide(NB_E)->m_id]->getPtrToSide(NB_N)->m_velocity.val);
					pressure = 0.5f*(fvm->getPtrToSide(NB_N)->m_pressure.val + FVM.FiniteVolumeMesh[fvm->getPtrToSide(NB_E)->m_id]->getPtrToSide(NB_N)->m_pressure.val);
					
					if (fvm->getPtrToSide(NB_W)->fv_boundary != nullptr)
						p1 = (fvm->FVPtr()->m_s.getCentre() + glm::vec3(-m_options.average_dim_steps.x / 2.0, m_options.average_dim_steps.y / 2.0, 0.0));
					if (FVM.FiniteVolumeMesh[fvm->getPtrToSide(NB_E)->m_id]->getPtrToSide(NB_E)->fv_boundary != nullptr) 
						p2 = (fvm->getPtrToSide(NB_E)->m_s.getCentre() + glm::vec3(m_options.average_dim_steps.x / 2.0, m_options.average_dim_steps.y / 2.0, 0.0));
					l = Line(p1, p2);
					//std::cout << UVc << ": {" << l.m_p1.x << "," << l.m_p1.y << "}, {" << l.m_p2.x << "," << l.m_p2.y << "}" << " NB_N ";
					UVolumes.push_back(new FiniteVolume(UVc, FINITE_VOLUME_TYPE::FVT_BOUNDARY, l, vel, pressure));
					UVolumes[UVc]->assignToBoundary(fvm->getPtrToSide(NB_N)->fv_boundary);
					FVMesh_U.push_back(new FVMeshNode(UVolumes[UVc], FPARAM::FU));
					for (int fvc = 0;fvc < UVc;fvc++) {
						for (auto l : UVolumes[fvc]->m_s.m_lineArray) {
							if (UVolumes[UVc]->m_l == l) {
								//std::cout << "adding FiniteVolumeMesh\n";
								FVMesh_U[UVc]->addNeighbour(UVolumes[fvc]);
								meshCount++;
								FVMesh_U[fvc]->addNeighbour(UVolumes[UVc]);
								meshCount++;
							}
						}
					}
					//std::cout << UVc << ": " << UVolumes[UVc]->fv_boundary << " NB_N ";
					UVc++;

				}
				if (fvm->getPtrToSide(NB_S)->fv_boundary != nullptr) {
					p1 = (fvm->FVPtr()->m_s.getCentre() + glm::vec3(0.0, -m_options.average_dim_steps.y / 2.0, 0.0));
					p2 = (fvm->getPtrToSide(NB_E)->m_s.getCentre() + glm::vec3(0.0, -m_options.average_dim_steps.y / 2.0, 0.0));

					vel = 0.5f*(fvm->getPtrToSide(NB_S)->m_velocity.val + FVM.FiniteVolumeMesh[fvm->getPtrToSide(NB_E)->m_id]->getPtrToSide(NB_S)->m_velocity.val);
					pressure = 0.5f*(fvm->getPtrToSide(NB_S)->m_pressure.val + FVM.FiniteVolumeMesh[fvm->getPtrToSide(NB_E)->m_id]->getPtrToSide(NB_S)->m_pressure.val);

					if (fvm->getPtrToSide(NB_W)->fv_boundary != nullptr)
						p1 = (fvm->FVPtr()->m_s.getCentre() + glm::vec3(-m_options.average_dim_steps.x / 2.0, -m_options.average_dim_steps.y / 2.0, 0.0));
					if (FVM.FiniteVolumeMesh[fvm->getPtrToSide(NB_E)->m_id]->getPtrToSide(NB_E)->fv_boundary != nullptr)
						p2 = (fvm->getPtrToSide(NB_E)->m_s.getCentre() + glm::vec3(m_options.average_dim_steps.x / 2.0, -m_options.average_dim_steps.y / 2.0, 0.0));
					l = Line(p1, p2);

					//std::cout << UVc << ": {" << l.m_p1.x << "," << l.m_p1.y << "}, {" << l.m_p2.x << "," << l.m_p2.y << "}" << " NB_S ";
					UVolumes.push_back(new FiniteVolume(UVc, FINITE_VOLUME_TYPE::FVT_BOUNDARY, l, vel, pressure));
					UVolumes[UVc]->assignToBoundary(fvm->getPtrToSide(NB_S)->fv_boundary);
					FVMesh_U.push_back(new FVMeshNode(UVolumes[UVc], FPARAM::FU));
					for (int fvc = 0;fvc < UVc;fvc++) {
						for (auto l : UVolumes[fvc]->m_s.m_lineArray) {
							if (UVolumes[UVc]->m_l == l) {
								//std::cout << "adding FiniteVolumeMesh\n";
								FVMesh_U[UVc]->addNeighbour(UVolumes[fvc]);
								meshCount++;
								FVMesh_U[fvc]->addNeighbour(UVolumes[UVc]);
								meshCount++;
							}
						}
					}
					//std::cout << UVc << ": " << UVolumes[UVc]->fv_boundary << " NB_S ";
					UVc++;


				}
			}
			else {
				l = fvm->getPtrToSide(NB_E)->m_l;
				vel = fvm->getPtrToSide(NB_E)->m_velocity.val;
				pressure = fvm->getPtrToSide(NB_E)->m_pressure.val;
				//std::cout << UVc << ": {" << l.m_p1.x << "," << l.m_p1.y << "}, {" << l.m_p2.x << "," << l.m_p2.y << "}" << " NB_E ";
				UVolumes.push_back(new FiniteVolume(UVc, FINITE_VOLUME_TYPE::FVT_BOUNDARY, l, vel, pressure));
				UVolumes[UVc]->assignToBoundary(fvm->getPtrToSide(NB_E)->fv_boundary);
				FVMesh_U.push_back(new FVMeshNode(UVolumes[UVc], FPARAM::FU));
				for (int fvc = 0;fvc < UVc;fvc++) {
					 for (auto l : UVolumes[fvc]->m_s.m_lineArray) {
						 if (UVolumes[UVc]->m_l == l) {
							 //std::cout << "adding FiniteVolumeMesh\n";
							 FVMesh_U[UVc]->addNeighbour(UVolumes[fvc]);
							 meshCount++;
							 FVMesh_U[fvc]->addNeighbour(UVolumes[UVc]);
							 meshCount++;
						 }
					 }
				}
				//std::cout << UVc << ": " << UVolumes[UVc]->fv_boundary << " NB_E ";
				UVc++;
			}
			//std::cout << "\n";
		}
	}
	std::cout << UVc << " finite volumes added\n";
	std::cout << meshCount << " neighbours in mesh\n";
	//int UVbc = UVc;
	//for (auto fvm : FVMesh_U) {
	//	fv_prec_3 p1;
	//	fv_prec_3 p2;
	//	if (fvm->getPtrToSide(NB_E) == nullptr) {
	//		p1 = fv_prec_3(fvm->fv.fv->m_Centre.x + m_options.average_dim_steps.x / 2.0, fvm->fv.fv->m_Centre.y - m_options.average_dim_steps.y / 2.0, 0.0);
	//		p2 = fv_prec_3(fvm->fv.fv->m_Centre.x + m_options.average_dim_steps.x / 2.0, fvm->fv.fv->m_Centre.y + m_options.average_dim_steps.y / 2.0, 0.0);
	//	}
	//	if (fvm->getPtrToSide(NB_W) == nullptr) {
	//		p1 = fv_prec_3(fvm->fv.fv->m_Centre.x - m_options.average_dim_steps.x / 2.0, fvm->fv.fv->m_Centre.y - m_options.average_dim_steps.y / 2.0, 0.0);
	//		p2 = fv_prec_3(fvm->fv.fv->m_Centre.x - m_options.average_dim_steps.x / 2.0, fvm->fv.fv->m_Centre.y + m_options.average_dim_steps.y / 2.0, 0.0);
	//	}
	//	if (fvm->getPtrToSide(NB_N) == nullptr) {
	//		p1 = fv_prec_3(fvm->fv.fv->m_Centre.x - m_options.average_dim_steps.x / 2.0, fvm->fv.fv->m_Centre.y + m_options.average_dim_steps.y / 2.0, 0.0);
	//		p2 = fv_prec_3(fvm->fv.fv->m_Centre.x + m_options.average_dim_steps.x / 2.0, fvm->fv.fv->m_Centre.y + m_options.average_dim_steps.y / 2.0, 0.0);
	//	}
	//	if (fvm->getPtrToSide(NB_S) == nullptr) {
	//		p1 = fv_prec_3(fvm->fv.fv->m_Centre.x - m_options.average_dim_steps.x / 2.0, fvm->fv.fv->m_Centre.y - m_options.average_dim_steps.y / 2.0, 0.0);
	//		p2 = fv_prec_3(fvm->fv.fv->m_Centre.x + m_options.average_dim_steps.x / 2.0, fvm->fv.fv->m_Centre.y - m_options.average_dim_steps.y / 2.0, 0.0);
	//	}
	//	Line tempLine(p1,p2);
	//	UVolumes.push_back(new FiniteVolume(UVbc, FINITE_VOLUME_TYPE::FVT_BOUNDARY, tempLine, {0.0,0.0,0.0}, 0.0));
	//	FVMesh_U.push_back(new FVMeshNode(UVolumes[UVbc]));
	//	for (int fvc = 0;fvc < UVbc;fvc++) {
	//		for (auto l : FVM.FiniteVolumes[fvc]->m_s.m_lineArray) {
	//			if (UVolumes[UVbc]->m_l == l) {
	//				std::cout << "adding FiniteVolumeMesh\n";
	//				FVMesh_U[UVbc]->addNeighbour(UVolumes[fvc]);
	//				FVMesh_U[fvc]->addNeighbour(UVolumes[UVbc]);
	//			}
	//		}
	//	}
	//	UVbc++;
	//}
	fv_prec Uresidual = 0.0;
	fv_prec residualScaleU = 0.0;
	fv_prec MaxUresidual = 0.0;

	std::vector<FiniteVolume*> VVolumes;
	std::vector<FVMeshNode*> FVMesh_V;
	std::cout << "Creating V mesh \n";
	//std::vector<FiniteVolumeMesh*> VVolumeMesh;
	meshCount = 0;
	int VVc = 0;
	for (auto fvm : FVM.FiniteVolumeMesh) {
		if ((fvm->getPtrToSide(NB_N)->fv_boundary == nullptr) and (fvm->FVPtr()->fv_boundary == nullptr)) {
			fv_prec_3 p1(fvm->FVPtr()->m_s.getCentre() + glm::vec3(-m_options.average_dim_steps.x / 2.0, 0.0, 0.0));
			fv_prec_3 p2(fvm->FVPtr()->m_s.getCentre() + glm::vec3(m_options.average_dim_steps.x / 2.0, 0.0, 0.0));
			fv_prec_3 p3(fvm->getPtrToSide(NB_N)->m_s.getCentre() + glm::vec3(-m_options.average_dim_steps.x / 2.0, 0.0, 0.0));
			fv_prec_3 p4(fvm->getPtrToSide(NB_N)->m_s.getCentre() + glm::vec3(m_options.average_dim_steps.x / 2.0, 0.0, 0.0));
			if (fvm->getPtrToSide(NB_S)->fv_boundary != nullptr) {
				p1 = (fvm->FVPtr()->m_s.getCentre() + glm::vec3(-m_options.average_dim_steps.x / 2.0, -m_options.average_dim_steps.y / 2.0, 0.0));
				p2 = (fvm->FVPtr()->m_s.getCentre() + glm::vec3(m_options.average_dim_steps.x / 2.0, -m_options.average_dim_steps.y / 2.0, 0.0));
			}
			if (FVM.FiniteVolumeMesh[fvm->getPtrToSide(NB_N)->m_id]->getPtrToSide(NB_N)->fv_boundary != nullptr) {
				p3 = (fvm->getPtrToSide(NB_N)->m_s.getCentre() + glm::vec3(-m_options.average_dim_steps.x / 2.0, m_options.average_dim_steps.y / 2.0, 0.0));
				p4 = (fvm->getPtrToSide(NB_N)->m_s.getCentre() + glm::vec3(m_options.average_dim_steps.x / 2.0, m_options.average_dim_steps.y / 2.0, 0.0));
			}
			//std::cout << VVc << ": {" << p1.x << "," << p1.y << "}, {" << p2.x << "," << p2.y << "}, {" << p3.x << "," << p3.y << "}, {" << p4.x << "," << p4.y << "}\n";
			Surface tempSurface(Line(p1, p2), Line(p2, p4), Line(p4, p3)); tempSurface.addLine(Line(p3, p1));
			tempSurface.calcCentre();
			VVolumes.push_back(new FiniteVolume(VVc, FINITE_VOLUME_TYPE::FVT_DEFAULT, tempSurface, glm::vec3(0.0f, 0.f, 0.f), 0.0));
			FVMesh_V.push_back(new FVMeshNode(VVolumes[VVc], FPARAM::FV));

			for (int fvc = 0;fvc < VVc;fvc++) {
				if (VVolumes[VVc]->m_s.ConnectedWith(VVolumes[fvc]->m_s)) {
					//std::cout << "adding FiniteVolumeMesh\n";
					FVMesh_V[VVc]->addNeighbour(VVolumes[fvc]);
					meshCount++;
					FVMesh_V[fvc]->addNeighbour(VVolumes[VVc]);
					meshCount++;
				}
			}
			VVc++;
		}
	}
	int Vsize = VVc;
	for (auto fvm : FVM.FiniteVolumeMesh) {
		fv_prec_3 p1;
		fv_prec_3 p2;
		Line l;
		fv_prec_3 vel;
		fv_prec pressure;
		if ((fvm->FVPtr()->fv_boundary == nullptr)) {
			//std::wcout << fvm->FVPtr()->m_id << " ";
			if ((fvm->getPtrToSide(NB_N)->fv_boundary == nullptr)) {
				if (fvm->getPtrToSide(NB_S)->fv_boundary != nullptr) {
					l = fvm->getPtrToSide(NB_S)->m_l;
					vel = fvm->getPtrToSide(NB_S)->m_velocity.val;
					pressure = fvm->getPtrToSide(NB_S)->m_pressure.val;
					//std::cout << VVc << ": {" << l.m_p1.x << "," << l.m_p1.y << "}, {" << l.m_p2.x << "," << l.m_p2.y << "}" << " NB_S ";
					VVolumes.push_back(new FiniteVolume(VVc, FINITE_VOLUME_TYPE::FVT_BOUNDARY, l, vel, pressure));
					VVolumes[VVc]->assignToBoundary(fvm->getPtrToSide(NB_S)->fv_boundary);
					FVMesh_V.push_back(new FVMeshNode(VVolumes[VVc], FPARAM::FU));
					for (int fvc = 0;fvc < VVc;fvc++) {
						for (auto l : VVolumes[fvc]->m_s.m_lineArray) {
							if (VVolumes[VVc]->m_l == l) {
								//std::cout << "adding FiniteVolumeMesh\n";
								FVMesh_V[VVc]->addNeighbour(VVolumes[fvc]);
								meshCount++;
								FVMesh_V[fvc]->addNeighbour(VVolumes[VVc]);
								meshCount++;
							}
						}
					}
					VVc++;

				}
				if (fvm->getPtrToSide(NB_E)->fv_boundary != nullptr) {
					p1 = fvm->getPtrToSide(NB_E)->m_Centre;
					p2 = FVM.FiniteVolumeMesh[fvm->getPtrToSide(NB_N)->m_id]->getPtrToSide(NB_E)->m_Centre;

					vel = 0.5f*(fvm->getPtrToSide(NB_E)->m_velocity.val + FVM.FiniteVolumeMesh[fvm->getPtrToSide(NB_N)->m_id]->getPtrToSide(NB_E)->m_velocity.val);
					pressure = 0.5f*(fvm->getPtrToSide(NB_E)->m_pressure.val + FVM.FiniteVolumeMesh[fvm->getPtrToSide(NB_N)->m_id]->getPtrToSide(NB_E)->m_pressure.val);

					if (fvm->getPtrToSide(NB_S)->fv_boundary != nullptr)
						p1 = (fvm->FVPtr()->m_s.getCentre() + glm::vec3(m_options.average_dim_steps.x / 2.0, -m_options.average_dim_steps.y / 2.0, 0.0));
					if (FVM.FiniteVolumeMesh[fvm->getPtrToSide(NB_N)->m_id]->getPtrToSide(NB_N)->fv_boundary != nullptr)
						p2 = (fvm->getPtrToSide(NB_N)->m_s.getCentre() + glm::vec3(m_options.average_dim_steps.x / 2.0, m_options.average_dim_steps.y / 2.0, 0.0));
					l = Line(p1, p2);
					//std::cout << VVc << ": {" << l.m_p1.x << "," << l.m_p1.y << "}, {" << l.m_p2.x << "," << l.m_p2.y << "}" << " NB_E ";
					VVolumes.push_back(new FiniteVolume(VVc, FINITE_VOLUME_TYPE::FVT_BOUNDARY, l, vel, pressure));
					VVolumes[VVc]->assignToBoundary(fvm->getPtrToSide(NB_E)->fv_boundary);
					FVMesh_V.push_back(new FVMeshNode(VVolumes[VVc], FPARAM::FU));
					for (int fvc = 0;fvc < VVc;fvc++) {
						for (auto l : VVolumes[fvc]->m_s.m_lineArray) {
							if (VVolumes[VVc]->m_l == l) {
								//std::cout << "adding FiniteVolumeMesh\n";
								FVMesh_V[VVc]->addNeighbour(VVolumes[fvc]);
								meshCount++;
								FVMesh_V[fvc]->addNeighbour(VVolumes[VVc]);
								meshCount++;
							}
						}
					}
					VVc++;

				}
				if (fvm->getPtrToSide(NB_W)->fv_boundary != nullptr) {
					p1 = fvm->getPtrToSide(NB_W)->m_Centre;
					p2 = FVM.FiniteVolumeMesh[fvm->getPtrToSide(NB_N)->m_id]->getPtrToSide(NB_W)->m_Centre;

					vel = 0.5f*(fvm->getPtrToSide(NB_W)->m_velocity.val + FVM.FiniteVolumeMesh[fvm->getPtrToSide(NB_N)->m_id]->getPtrToSide(NB_W)->m_velocity.val);
					pressure = 0.5f*(fvm->getPtrToSide(NB_W)->m_pressure.val + FVM.FiniteVolumeMesh[fvm->getPtrToSide(NB_N)->m_id]->getPtrToSide(NB_W)->m_pressure.val);

					if (fvm->getPtrToSide(NB_S)->fv_boundary != nullptr)
						p1 = (fvm->FVPtr()->m_s.getCentre() + glm::vec3(-m_options.average_dim_steps.x / 2.0, -m_options.average_dim_steps.y / 2.0, 0.0));
					if (FVM.FiniteVolumeMesh[fvm->getPtrToSide(NB_N)->m_id]->getPtrToSide(NB_N)->fv_boundary != nullptr)
						p2 = (fvm->getPtrToSide(NB_N)->m_s.getCentre() + glm::vec3(-m_options.average_dim_steps.x / 2.0, m_options.average_dim_steps.y / 2.0, 0.0));
					l = Line(p1, p2);
					//std::cout << VVc << ": {" << l.m_p1.x << "," << l.m_p1.y << "}, {" << l.m_p2.x << "," << l.m_p2.y << "}" << " NB_W ";
					VVolumes.push_back(new FiniteVolume(VVc, FINITE_VOLUME_TYPE::FVT_BOUNDARY, l, vel, pressure));
					VVolumes[VVc]->assignToBoundary(fvm->getPtrToSide(NB_W)->fv_boundary);
					FVMesh_V.push_back(new FVMeshNode(VVolumes[VVc], FPARAM::FU));
					for (int fvc = 0;fvc < VVc;fvc++) {
						for (auto l : VVolumes[fvc]->m_s.m_lineArray) {
							if (VVolumes[VVc]->m_l == l) {
								//std::cout << "adding FiniteVolumeMesh\n";
								FVMesh_V[VVc]->addNeighbour(VVolumes[fvc]);
								meshCount++;
								FVMesh_V[fvc]->addNeighbour(VVolumes[VVc]);
								meshCount++;
							}
						}
					}
					VVc++;


				}
			}
			else {
				l = fvm->getPtrToSide(NB_N)->m_l;
				vel = fvm->getPtrToSide(NB_N)->m_velocity.val;
				pressure = fvm->getPtrToSide(NB_N)->m_pressure.val;
				//std::cout << VVc << ": {" << l.m_p1.x << "," << l.m_p1.y << "}, {" << l.m_p2.x << "," << l.m_p2.y << "}" << " NB_N ";
				VVolumes.push_back(new FiniteVolume(VVc, FINITE_VOLUME_TYPE::FVT_BOUNDARY, l, vel, pressure));
				VVolumes[VVc]->assignToBoundary(fvm->getPtrToSide(NB_N)->fv_boundary);
				FVMesh_V.push_back(new FVMeshNode(VVolumes[VVc], FPARAM::FU));
				for (int fvc = 0;fvc < VVc;fvc++) {
					for (auto l : VVolumes[fvc]->m_s.m_lineArray) {
						if (VVolumes[VVc]->m_l == l) {
							//std::cout << "adding FiniteVolumeMesh\n";
							FVMesh_V[VVc]->addNeighbour(VVolumes[fvc]);
							meshCount++;
							FVMesh_V[fvc]->addNeighbour(VVolumes[VVc]);
							meshCount++;
						}
					}
				}
				VVc++;
			}
			//std::cout << "\n";
		}
	}
	std::cout << VVc << " finite volumes added\n";
	std::cout << meshCount << " neighbours in mesh\n";
	fv_prec Vresidual = 0.0;
	fv_prec residualScaleV = 0.0;
	fv_prec MaxVresidual = 0.0;

	std::vector<int> last_x_fvs;
	std::vector<int> last_y_fvs;
	fv_prec max_x = 0.0; fv_prec max_y = 0.0; fv_prec max_z = 0.0;
	for (int i = 0;i < FVM.FiniteVolumeMesh.size();i++) {
		if (max_x < FVM.FiniteVolumeMesh[i]->FVPtr()->m_Centre.x) max_x = FVM.FiniteVolumeMesh[i]->FVPtr()->m_Centre.x;
		if (max_y < FVM.FiniteVolumeMesh[i]->FVPtr()->m_Centre.y) max_y = FVM.FiniteVolumeMesh[i]->FVPtr()->m_Centre.y;
		if (max_z < FVM.FiniteVolumeMesh[i]->FVPtr()->m_Centre.z) max_z = FVM.FiniteVolumeMesh[i]->FVPtr()->m_Centre.z;
	}
	max_x = 1.00 - m_options.average_dim_steps.x;
	max_y = 1.00 - m_options.average_dim_steps.y;
	for (int i = 0;i < FVM.FiniteVolumeMesh.size();i++) {
		if (FVM.FiniteVolumeMesh[i]->FVPtr()->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
			if (FVM.FiniteVolumeMesh[i]->FVPtr()->m_Centre.x >= max_x) { last_x_fvs.push_back(i); }
			if (FVM.FiniteVolumeMesh[i]->FVPtr()->m_Centre.y >= max_y) { last_y_fvs.push_back(i); }
		}
	}

	fv_prec sigma = 1.0;

	Ch_<fv_prec> sum;
	sum.dval = 10.0;
	sum.val = 0.0;

	fv_prec Presidual = 0.0;
	fv_prec MaxPresidual = 0.0;


	std::vector<fv_prec> PApr;
	PApr.resize(FVM.FiniteVolumes.size(), 0.0);
	std::vector<fv_prec> Pcor;
	Pcor.resize(FVM.FiniteVolumes.size(), 0.0);


	Ch_<std::vector<fv_prec>> UApr;
	UApr.val.resize(Usize, 0.0);
	UApr.dval.resize(Usize, 0.0);
	std::vector<fv_prec> UCor_ed;
	UCor_ed.resize(Usize, 0.0);
	std::vector<fv_prec> UCor0;
	UCor0.resize(Usize, 0.0);

	Ch_<std::vector<fv_prec>> VApr;
	VApr.val.resize(Vsize, 0.0);
	VApr.dval.resize(Vsize, 0.0);
	std::vector<fv_prec> VCor_ed;
	VCor_ed.resize(Vsize, 0.0);
	std::vector<fv_prec> VCor0;
	VCor0.resize(Vsize, 0.0);



	std::vector<FVMeshNode*> FVMesh_PC;
	for (int i = 0;i < FVM.FiniteVolumeMesh.size();i++) {
		if (FVM.FiniteVolumeMesh[i]->FVPtr()->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
			FVMesh_PC.push_back(new FVMeshNode(*(FVM.FiniteVolumeMesh[i]), PHI_PARAMS::P));
		}
	}

	for (int i = 0;i < FVMesh_PC.size();i++) {
		FVMesh_PC[i]->fv.fv->m_density.dval = FVMesh_PC[i]->fv.fv->m_density.val;
		FVMesh_PC[i]->fv.fv->m_velocity.dval = FVMesh_PC[i]->fv.fv->m_velocity.val;
		FVMesh_PC[i]->fv.fv->m_pressure.dval = FVMesh_PC[i]->fv.fv->m_pressure.val;
	}


	bool whileCondition = true;
	while (whileCondition) {
		if (AlgoIteration > maxIteration) break;

		std::wcout << "Iteration number " << AlgoIteration << "\n";
		Presidual = 0.0;
		MaxPresidual = 0.0;

		Uresidual = 0.0;
		MaxUresidual = 0.0;
		Vresidual = 0.0;
		MaxVresidual = 0.0;

		residualScale = 0.0;
		residualScaleU = 0.0;
		residualScaleV = 0.0;
		for (auto fvm : FVMesh_U) {
			//  |_____|_____|_____|_____|
			//  |     |     |     |     |
			//  |  0  |  1  |  2  |  3  |
			//  |_____|_____|_____|_____|
			//
			//  |________|_____|________|
			//  |        |     |        |
			//  |     0' |  1' |  2'    |
			//  |________|_____|________|
			//  D0' = (D1(x0'-x0)+D0(x1-x0'))/(x1-x0)

			if (fvm->fv.fv->m_type == FVT_DEFAULT) {
				//std::cout << "(" << FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + 1 + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->FVPtr()->m_id << "*(" << fvm->fv.fv->m_id << "," << FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->FVPtr()->m_id << ")+" << FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->FVPtr()->m_id << "*(" << FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + 1 + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->FVPtr()->m_id << "," << fvm->fv.fv->m_id << "))/(" << FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + 1 + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->FVPtr()->m_id << "-" << FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->FVPtr()->m_id << ")\n";
				fvm->fv.fv->m_DVisc = (FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + 1 + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->FVPtr()->m_DVisc* fvm->fv.fv->distance(FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->FVPtr()) + FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->FVPtr()->m_DVisc* FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + 1 + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->FVPtr()->distance(fvm->fv.fv)) / (FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + 1 + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->FVPtr()->distance(FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->FVPtr()));
				fvm->fv.fv->m_density.dval = (FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + 1 + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->FVPtr()->m_density.dval* fvm->fv.fv->distance(FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->FVPtr()) + FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->FVPtr()->m_density.dval* FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + 1 + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->FVPtr()->distance(fvm->fv.fv)) / (FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + 1 + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->FVPtr()->distance(FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->FVPtr()));
			}
		}
		for (auto fvm : FVMesh_V) {
			//   _____ _     _____ _
			//  |     |     |     |
			//  |  12 |     |  8' |
			//  |_____|_    |     |
			//  |     |     |_____|_
			//  |  8  |     |     |
			//  |_____|_    |  4' |
			//  |     |     |     |
			//  |  4  |     |_____|_
			//  |_____|_    |     |
			//  |     |     |  0' |
			//  |  0  |     |     |
			//  |_____|_    |_____|_
			//  D0' = (D4(y0'-y0)+D0(y4-y0'))/(y4-y0)

			if (fvm->fv.fv->m_type == FVT_DEFAULT) {
				//std::cout << "(" << FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + last_x_fvs.size()]->FVPtr()->m_id << "*(" << fvm->fv.fv->m_id << "," << FVM.FiniteVolumeMesh[fvm->fv.fv->m_id]->FVPtr()->m_id << ")+" << FVM.FiniteVolumeMesh[fvm->fv.fv->m_id]->FVPtr()->m_id << "*(" << FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + last_x_fvs.size()]->FVPtr()->m_id << "," << fvm->fv.fv->m_id << "))/(" << FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + last_x_fvs.size()]->FVPtr()->m_id << "-" << FVM.FiniteVolumeMesh[fvm->fv.fv->m_id]->FVPtr()->m_id << ")=";
				fvm->fv.fv->m_DVisc = (FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + last_x_fvs.size()]->FVPtr()->m_DVisc* fvm->fv.fv->distance(FVM.FiniteVolumeMesh[fvm->fv.fv->m_id]->FVPtr()) + FVM.FiniteVolumeMesh[fvm->fv.fv->m_id]->FVPtr()->m_DVisc* FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + last_x_fvs.size()]->FVPtr()->distance(fvm->fv.fv)) / (FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + last_x_fvs.size()]->FVPtr()->distance(FVM.FiniteVolumeMesh[fvm->fv.fv->m_id]->FVPtr()));
				fvm->fv.fv->m_density.dval = (FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + last_x_fvs.size()]->FVPtr()->m_density.dval* fvm->fv.fv->distance(FVM.FiniteVolumeMesh[fvm->fv.fv->m_id]->FVPtr()) + FVM.FiniteVolumeMesh[fvm->fv.fv->m_id]->FVPtr()->m_density.dval* FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + last_x_fvs.size()]->FVPtr()->distance(fvm->fv.fv)) / (FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + last_x_fvs.size()]->FVPtr()->distance(FVM.FiniteVolumeMesh[fvm->fv.fv->m_id]->FVPtr()));
			}
		}


		// ‘P*(aP0 + sigma*SUM(aNB) - Spm*Vol*sigma) = SUM(‘NB * aNB) + ‘P0*(aP0 - (1-sigma)*SUM(aNB) - (1-sigma)*Spm0*Vol) + (1-sigma)*SUM(‘NB0 * aNB) + (1-sigma)*Spp0*Vol + sigma*Spp*Vol + dy*(P*P - P*E)
		//      _          _             _               _              _              _                    _                      _                          _                  _              -            
		for (auto fvm : FVMesh_U) {
			fv_prec_3 Area_centre;
			fv_prec Area;
			fv_prec Volume;
			fv_prec_3 InNormal;
			switch (fvm->fv.fv->FVdim) {
			case(FVA_S):
				Volume = fvm->fv.fv->m_s.size();
				break;
			default:
				break;
			}

			fv_prec effD;
			fv_prec_3 vec_fv_to_nb;
			fv_prec_3 vec_fv_to_b;
			fv_prec_3 vec_b_to_nb;
			fvm->fv.fv_coef = 0.0;
			fvm->fv_source = 0.0;

			fv_prec boundary_coefs = 0.0;
			fv_prec b_coef = 0.0;
			fv_prec b_vel = 0.0;

			if (fvm->fv.fv->m_type == FVT_DEFAULT) {

				for (int nb = 0;nb < fvm->nb_fvs.size();nb++) {
					FV_NB relativ_pos;
					fv_prec Dnb_dir = 0.0;
					fv_prec Fnb_dir = 0.0;
					switch (fvm->fv.fv->FVdim) {
					case(FVA_S):
						for (auto fv_p : fvm->fv.fv->m_s.m_lineArray) {
							bool NBIsFound = false;

							if (fvm->nb_fvs[nb].fv->FVdim == FVA_L) {
								//std::cout << "{" << fv_p.m_p1.x << "," << fv_p.m_p1.y << "},{" << fv_p.m_p2.x << "," << fv_p.m_p2.y << "}-{" << nb.fv->m_l.m_p1.x << "," << nb.fv->m_l.m_p1.y << "},{" << nb.fv->m_l.m_p2.x << "," << nb.fv->m_l.m_p2.y << "}" << nb.fv->m_type << "\n";
								if (fv_p == fvm->nb_fvs[nb].fv->m_l) {
									Area = fv_p.size();
									Area_centre = fv_p.getCentre();
									if (Area_centre.x > fvm->El().fv->m_Centre.x) relativ_pos = NB_E;
									if (Area_centre.x < fvm->El().fv->m_Centre.x) relativ_pos = NB_W;
									if (Area_centre.y > fvm->El().fv->m_Centre.y) relativ_pos = NB_N;
									if (Area_centre.y < fvm->El().fv->m_Centre.y) relativ_pos = NB_S;
									NBIsFound = true;
								}
							}
							else if (fvm->nb_fvs[nb].fv->FVdim == FVA_S) {
								for (auto fv_nb : fvm->nb_fvs[nb].fv->m_s.m_lineArray) {
									if (fv_p == fv_nb) {
										Area = fv_p.size();
										Area_centre = fv_p.getCentre();
										if (Area_centre.x > fvm->fv.fv->m_Centre.x) relativ_pos = NB_E;
										if (Area_centre.x < fvm->fv.fv->m_Centre.x) relativ_pos = NB_W;
										if (Area_centre.y > fvm->fv.fv->m_Centre.y) relativ_pos = NB_N;
										if (Area_centre.y < fvm->fv.fv->m_Centre.y) relativ_pos = NB_S;
										NBIsFound = true;
										break;
									}
								}
							}
							if (!NBIsFound) continue;

							vec_fv_to_nb = fvm->fv.fv->m_Centre - fvm->nb_fvs[nb].fv->m_Centre;
							vec_fv_to_b = fvm->fv.fv->m_Centre - Area_centre;
							vec_b_to_nb = Area_centre - fvm->nb_fvs[nb].fv->m_Centre;
							InNormal = glm::normalize(-vec_fv_to_b);

							//std::cout << fvm->nb_fvs[nb].fv->m_type << "  ";
							//std::cout << fvm->fv.fv->m_id << " (";
							//switch (relativ_pos) {
							//case(NB_E):
							//	std::cout << "NB_E";
							//	break;
							//case(NB_W):
							//	std::cout << "NB_W";
							//	break;
							//case(NB_N):
							//	std::cout << "NB_N";
							//	break;
							//case(NB_S):
							//	std::cout << "NB_S";
							//	break;
							//default:
							//	break;
							//}
							//std::cout << "): " << fvm->nb_fvs[nb].fv->m_id << " ";
							if (fvm->nb_fvs[nb].fv->m_type == FVT_DEFAULT) {
								//std::cout << "\n";

								fv_prec p_coef = fvm->nb_fvs[nb].fv->m_DVisc*vec_mod(vec_fv_to_b) / vec_mod(vec_fv_to_nb);;
								fv_prec nb_coef = fvm->fv.fv->m_DVisc*vec_mod(vec_b_to_nb)/ vec_mod(vec_fv_to_nb);
								effD = fvm->nb_fvs[nb].fv->m_DVisc*fvm->fv.fv->m_DVisc / (p_coef + nb_coef);
								Dnb_dir = effD/ vec_mod(vec_fv_to_nb);
								Fnb_dir = 0.5*(fvm->fv.fv->m_density.dval*fvm->fv.fv->m_velocity.dval.x + fvm->nb_fvs[nb].fv->m_density.dval*fvm->nb_fvs[nb].fv->m_velocity.dval.x);
								switch (relativ_pos) {
								case NB_E:
								case NB_W:
									Fnb_dir = 0.5*(fvm->fv.fv->m_density.dval*fvm->fv.fv->m_velocity.dval.x + fvm->nb_fvs[nb].fv->m_density.dval*fvm->nb_fvs[nb].fv->m_velocity.dval.x);
									break;
								case NB_N:
									Fnb_dir = 0.5*(VVolumes[FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->FVPtr()->m_id]->m_density.dval*VVolumes[FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->FVPtr()->m_id]->m_velocity.dval.y + VVolumes[FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + 1 + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->FVPtr()->m_id]->m_density.dval *VVolumes[FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + 1 + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->FVPtr()->m_id]->m_velocity.dval.y);
									break;
								case NB_S:
									Fnb_dir = 0.5*(VVolumes[FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->getPtrToSide(NB_S)->m_id]->m_density.dval*VVolumes[FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->getPtrToSide(NB_S)->m_id]->m_velocity.dval.y + VVolumes[FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + 1 + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->getPtrToSide(NB_S)->m_id]->m_density.dval*VVolumes[FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + 1 + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->getPtrToSide(NB_S)->m_id]->m_velocity.dval.y);
									break;
								default:
									break;
								}
								Dnb_dir *= Area;
								Fnb_dir *= Area;
								//std::cout << "Fnb_dir: " << Fnb_dir << " Dnb_dir:" << Dnb_dir << " - ";
								// aNB							
								{
									// Diffusion
									fvm->nb_fvs[nb].fv_coef = PeFunc(Fnb_dir / Dnb_dir) * Dnb_dir;
									//std::cout << fvm->nb_fvs[nb].fv_coef << " - >";
									// Convention
									//fvm->nb_fvs[nb].fv_coef += glm::max(0.0f, Fnb_dir);
									if (relativ_pos == NB_E) fvm->nb_fvs[nb].fv_coef += glm::max(-Fnb_dir, 0.0f);
									if (relativ_pos == NB_W) fvm->nb_fvs[nb].fv_coef += glm::max(Fnb_dir, 0.0f);
									if (relativ_pos == NB_N) fvm->nb_fvs[nb].fv_coef += glm::max(-Fnb_dir, 0.0f);
									if (relativ_pos == NB_S) fvm->nb_fvs[nb].fv_coef += glm::max(Fnb_dir, 0.0f);
									//std::cout << fvm->nb_fvs[nb].fv_coef << " ";
								}
								//fvm->nb_fvs[nb].fv_coef *= Area;
								//boundary_coefs = 0.0;
							}
							else if (fvm->nb_fvs[nb].fv->m_type == FVT_BOUNDARY) {

								b_coef = fvm->nb_fvs[nb].fv->fv_boundary->getVelocityCoefs().x;
								//std::cout << b_coef << " ";
								b_vel = fvm->nb_fvs[nb].fv->fv_boundary->getVelocity().x;
								//std::cout << b_vel << " ";

								fv_prec p_coef = fvm->El().fv->m_DVisc / vec_mod(vec_fv_to_b);
								boundary_coefs = 1.0 / (1.0 / p_coef + 1.0 / b_coef)*Area;
								//boundary_coefs = 1.0 / (1.0 / (effD / (vec_mod(vec_fv_to_b))) + 1.0 / b_coef)*Area;

								//if ((relativ_pos == NB_E) or (relativ_pos == NB_W)) boundary_coefs -= 0.001;
								//std::cout << boundary_coefs;

								fvm->fv.fv_coef += boundary_coefs * sigma;
								fvm->fv_source -= fvm->fv.fv->m_velocity.dval.x * boundary_coefs * (1.0 - sigma);
								fvm->fv_source += boundary_coefs * b_vel * sigma; // *(1 - sigma)
								//std::cout << fvm->fv_source << "->";
								fvm->fv_source += boundary_coefs * b_vel * (1.0 - sigma); //old boundary velocity
								//std::cout << fvm->fv_source << "->";

							}
							/*
							else if (fvm->nb_fvs[nb].fv->m_type == FVT_BOUNDARY) {
								if (Area_centre.x > fvm->fv.fv->m_Centre.x) relativ_pos = NB_E;
								if (Area_centre.x < fvm->fv.fv->m_Centre.x) relativ_pos = NB_W;
								if (Area_centre.y > fvm->fv.fv->m_Centre.y) relativ_pos = NB_N;
								if (Area_centre.y < fvm->fv.fv->m_Centre.y) relativ_pos = NB_S;
								vec_fv_to_b = fvm->fv.fv->m_Centre - Area_centre;

								std::cout << fvm->fv.fv->m_id << " (";
								switch (relativ_pos) {
								case(NB_E):
									std::cout << "NB_E";
									break;
								case(NB_W):
									std::cout << "NB_W";
									break;
								case(NB_N):
									std::cout << "NB_N";
									break;
								case(NB_S):
									std::cout << "NB_S";
									break;
								default:
									break;
								}
								std::cout << "): ";

								switch (relativ_pos) {
								case(NB_E):
									std::cout << FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + 1 + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->getPtrToSide(NB_E)->m_id << "\n";
									effD = fvm->El().fv->m_DVisc;
									b_coef = FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + 1 + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->getPtrToSide(NB_E)->fv_boundary->getVelocityCoefs().x;
									b_vel = FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + 1 + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->getPtrToSide(NB_E)->fv_boundary->getVelocity().x;
									boundary_coefs += 1.0 / (1.0 / (effD / (2.0*vec_mod(vec_fv_to_b))) + 1.0 / b_coef)*Area;
									break;
								case(NB_W):
									std::cout << FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->getPtrToSide(NB_W)->m_id << "\n";
									effD = fvm->El().fv->m_DVisc;
									b_coef = FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->getPtrToSide(NB_W)->fv_boundary->getVelocityCoefs().x;
									b_vel = FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->getPtrToSide(NB_W)->fv_boundary->getVelocity().x;
									boundary_coefs += 1.0 / (1.0 / (effD / (2.0*vec_mod(vec_fv_to_b))) + 1.0 / b_coef)*Area;
									break;
								case(NB_N):
									std::cout << FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->getPtrToSide(NB_N)->m_id << "  " << FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + 1 + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->getPtrToSide(NB_N)->m_id << "\n";
									effD = fvm->El().fv->m_DVisc;
									b_coef = 0.5*(FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->getPtrToSide(NB_N)->fv_boundary->getVelocityCoefs().x + FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + 1 + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->getPtrToSide(NB_N)->fv_boundary->getVelocityCoefs().x);
									b_vel = 0.5*(FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->getPtrToSide(NB_N)->fv_boundary->getVelocity().x + FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + 1 + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->getPtrToSide(NB_N)->fv_boundary->getVelocity().x);
									boundary_coefs += 1.0 / (1.0 / (effD / vec_mod(vec_fv_to_b)) + 1.0 / b_coef)*Area;
									break;
								case(NB_S):
									std::cout << FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->getPtrToSide(NB_S)->m_id << "  " << FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + 1 + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->getPtrToSide(NB_S)->m_id << "\n";
									effD = fvm->El().fv->m_DVisc;
									b_coef = 0.5*(FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->getPtrToSide(NB_S)->fv_boundary->getVelocityCoefs().x + FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + 1 + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->getPtrToSide(NB_N)->fv_boundary->getVelocityCoefs().x);
									b_vel = 0.5*(FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->getPtrToSide(NB_S)->fv_boundary->getVelocity().x + FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + 1 + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->getPtrToSide(NB_N)->fv_boundary->getVelocity().x);
									boundary_coefs += 1.0 / (1.0 / (effD / vec_mod(vec_fv_to_b)) + 1.0 / b_coef)*Area;
									break;
								default:
									break;
								}
							}
							*/
							//std::cout << "\n";
							//fvm->fv.fv_coef -= Fnb_dir * Area;
							if ((relativ_pos == NB_E) or (relativ_pos == NB_W)) fvm->fv.fv_coef += Fnb_dir * InNormal.x * Area;
							if ((relativ_pos == NB_N) or (relativ_pos == NB_S)) fvm->fv.fv_coef += Fnb_dir * InNormal.y * Area;

							//if (relativ_pos == NB_E) fvm->fv.fv_coef += Fnb_dir * Area;
							//if (relativ_pos == NB_W) fvm->fv.fv_coef += -Fnb_dir * Area;
							//if (relativ_pos == NB_N) fvm->fv.fv_coef += Fnb_dir * Area;
							//if (relativ_pos == NB_S) fvm->fv.fv_coef += -Fnb_dir * Area;
						}
						break;
					default:
						break;
					}
					// sigma*SUM(aNB)
					fvm->fv.fv_coef += fvm->nb_fvs[nb].fv_coef * sigma;
					//
					// ‘P0*(1-sigma)*SUM(aNB)
					fvm->fv_source -= fvm->fv.fv->m_velocity.dval.x * fvm->nb_fvs[nb].fv_coef * (1.0 - sigma); // *(1 - sigma)
					// (1-sigma)*SUM(‘NB0 * aNB)
					fvm->fv_source += fvm->nb_fvs[nb].fv->m_velocity.dval.x * fvm->nb_fvs[nb].fv_coef * (1.0 - sigma); // *(1 - sigma)







				}
				// Pressure
				for (auto Pcor_fv_ls : FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->FVPtr()->m_s.m_lineArray) {
					for (auto Pcor_nb_ls : FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->getPtrToSide(NB_E)->m_s.m_lineArray) {
						if (Pcor_fv_ls == Pcor_nb_ls) {
							Area = Pcor_fv_ls.size();
						}
					}
				}
				// dy*(P*P - P*E)
				fvm->fv_source += Area * (PApr[fvm->fv.fv->m_id + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)] - PApr[fvm->fv.fv->m_id + 1 + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]);
				//std::cout << fvm->fv.fv->m_id << " -{" << fvm->fv.fv->m_id +(fvm->fv.fv->m_id)/(last_y_fvs.size()-1) << ", " << fvm->fv.fv->m_id + 1 + (fvm->fv.fv->m_id ) / (last_y_fvs.size() - 1) << "}\n";

				// Spm*Vol*sigma
				fvm->fv.fv_coef -= fvm->sf.Spm.val * Volume*sigma;
				// sigma*Spp*Vol
				fvm->fv_source += fvm->sf.Spp.val * Volume*sigma;
				// ‘P0*(1-sigma)*Spm0*Vol
				fvm->fv_source -= fvm->fv.fv->m_velocity.val.x * fvm->sf.Spm.dval * Volume*(1.0 - sigma);
				// (1-sigma)*Spp0*Vol
				fvm->fv_source += fvm->sf.Spp.dval * Volume*(1.0 - sigma);


				if (dt == 0) {}
				else {
					// aP0
					fvm->fv.fv_coef += fvm->fv.fv->m_density.dval * Volume / dt;
					// ‘P0*aP0
					fvm->fv_source += fvm->fv.fv->m_velocity.val.x* fvm->fv.fv->m_density.val * Volume / dt;
				}


				//std::cout << fvm->getPtrToSide(NB_E) << " " << fvm->getPtrToSide(NB_W) << " " << fvm->getPtrToSide(NB_N) << " " << fvm->getPtrToSide(NB_S) <<"\n";
				/*
				fv_prec boundary_coefs = 0.0;
				fv_prec b_coef = 0.0;
				fv_prec b_vel = 0.0;
				if (fvm->getPtrToSide(NB_E) == nullptr) {
					bool ex_flag = false;
					for (auto fv_p : fvm->fv.fv->m_s.m_lineArray) {
						if (fv_p.getCentre().x > fvm->fv.fv->m_Centre.x) {
							Area = fv_p.size();
							Area_centre = fv_p.getCentre();
							ex_flag = true;
							break;
						}
					}
					if (!ex_flag) continue;
					vec_fv_to_b = fvm->fv.fv->m_Centre - Area_centre;
					//std::cout << fvm->fv.fv->m_id << " (" << "NB_E-" << Area_centre.x << "," << Area_centre.y << "," << Area_centre.z <<"): " << FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + 1 + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->getPtrToSide(NB_E)->m_id << "\n";

					effD = fvm->El().fv->m_DVisc;
					b_coef = fvm->getPtrToSide(NB_E)->fv_boundary->getVelocityCoefs().x;
					b_vel = fvm->getPtrToSide(NB_E)->fv_boundary->getVelocity().x;
					//b_coef = FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + 1 + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->getPtrToSide(NB_E)->fv_boundary->getVelocityCoefs().x;
					//b_vel = FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + 1 + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->getPtrToSide(NB_E)->fv_boundary->getVelocity().x;
					//boundary_coefs += 1.0 / (1.0 / (effD / (2.0*vec_mod(vec_fv_to_b))) + 1.0 / b_coef)*Area;
					boundary_coefs = 1.0 / (1.0 / (effD / (vec_mod(vec_fv_to_b))) + 1.0 / b_coef)*Area;

					fvm->fv.fv_coef += boundary_coefs * sigma;
					fvm->fv_source -= fvm->fv.fv->m_velocity.dval.x * boundary_coefs * (1.0 - sigma);
					fvm->fv_source += boundary_coefs * b_vel * sigma; // *(1 - sigma)
					//std::cout << fvm->fv_source << "->";
					fvm->fv_source += boundary_coefs * b_vel * (1.0 - sigma); //old boundary velocity
					//std::cout << fvm->fv_source << "->";
				}
				if (fvm->getPtrToSide(NB_W) == nullptr) {
					bool ex_flag = false;
					for (auto fv_p : fvm->fv.fv->m_s.m_lineArray) {
						if (fv_p.getCentre().x < fvm->fv.fv->m_Centre.x) {
							Area = fv_p.size();
							Area_centre = fv_p.getCentre();
							ex_flag = true;
							break;
						}
					}
					if (!ex_flag) continue;
					vec_fv_to_b = fvm->fv.fv->m_Centre - Area_centre;
					//std::cout << fvm->fv.fv->m_id << " (" << "NB_W-" << Area_centre.x << "," << Area_centre.y << "," << Area_centre.z << "): " << FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->getPtrToSide(NB_W)->m_id << "\n";
					effD = fvm->El().fv->m_DVisc;
					b_coef = fvm->getPtrToSide(NB_W)->fv_boundary->getVelocityCoefs().x;
					b_vel = fvm->getPtrToSide(NB_W)->fv_boundary->getVelocity().x;
					//b_coef = FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->getPtrToSide(NB_W)->fv_boundary->getVelocityCoefs().x;
					//b_vel = FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->getPtrToSide(NB_W)->fv_boundary->getVelocity().x;
					//boundary_coefs += 1.0 / (1.0 / (effD / (2.0*vec_mod(vec_fv_to_b))) + 1.0 / b_coef)*Area;
					boundary_coefs = 1.0 / (1.0 / (effD / (vec_mod(vec_fv_to_b))) + 1.0 / b_coef)*Area;

					fvm->fv.fv_coef += boundary_coefs * sigma;
					fvm->fv_source -= fvm->fv.fv->m_velocity.dval.x * boundary_coefs * (1.0 - sigma);
					fvm->fv_source += boundary_coefs * b_vel * sigma; // *(1 - sigma)
					//std::cout << fvm->fv_source << "->";
					fvm->fv_source += boundary_coefs * b_vel * (1.0 - sigma); //old boundary velocity
					//std::cout << fvm->fv_source << "->";
				}
				if (fvm->getPtrToSide(NB_N) == nullptr) {
					bool ex_flag = false;
					for (auto fv_p : fvm->fv.fv->m_s.m_lineArray) {
						if (fv_p.getCentre().y > fvm->fv.fv->m_Centre.y) {
							Area = fv_p.size();
							Area_centre = fv_p.getCentre();
							ex_flag = true;
							break;
						}
					}
					if (!ex_flag) continue;
					vec_fv_to_b = fvm->fv.fv->m_Centre - Area_centre;
					//std::cout << fvm->fv.fv->m_id << " (" << "NB_N-" << Area_centre.x << "," << Area_centre.y << "," << Area_centre.z << "): " << FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->getPtrToSide(NB_N)->m_id << "  " << FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + 1 + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->getPtrToSide(NB_N)->m_id << "\n";
					effD = fvm->El().fv->m_DVisc;
					b_coef = fvm->getPtrToSide(NB_N)->fv_boundary->getVelocityCoefs().x;
					b_vel = fvm->getPtrToSide(NB_N)->fv_boundary->getVelocity().x;
					//b_coef = 0.5*(FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->getPtrToSide(NB_N)->fv_boundary->getVelocityCoefs().x + FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + 1 + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->getPtrToSide(NB_N)->fv_boundary->getVelocityCoefs().x);
					//b_vel = 0.5*(FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->getPtrToSide(NB_N)->fv_boundary->getVelocity().x + FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + 1 + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->getPtrToSide(NB_N)->fv_boundary->getVelocity().x);
					boundary_coefs = 1.0 / (1.0 / (effD / vec_mod(vec_fv_to_b)) + 1.0 / b_coef)*Area;
					//std::cout << boundary_coefs << "\n";
					fvm->fv.fv_coef += boundary_coefs * sigma;
					fvm->fv_source -= fvm->fv.fv->m_velocity.dval.x * boundary_coefs * (1.0 - sigma);
					fvm->fv_source += boundary_coefs * b_vel * sigma; // *(1 - sigma)
					//std::cout << fvm->fv_source << "->";
					fvm->fv_source += boundary_coefs * b_vel * (1.0 - sigma); //old boundary velocity
					//std::cout << fvm->fv_source << "->";
				}
				if (fvm->getPtrToSide(NB_S) == nullptr) {
					bool ex_flag = false;
					for (auto fv_p : fvm->fv.fv->m_s.m_lineArray) {
						if (fv_p.getCentre().y < fvm->fv.fv->m_Centre.y) {
							Area = fv_p.size();
							Area_centre = fv_p.getCentre();
							ex_flag = true;
							break;
						}
					}
					if (!ex_flag) continue;
					vec_fv_to_b = fvm->fv.fv->m_Centre - Area_centre;
					//std::cout << fvm->fv.fv->m_id << " (" << "NB_S-" << Area_centre.x << "," << Area_centre.y << "," << Area_centre.z << "): " << FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->getPtrToSide(NB_S)->m_id << "  " << FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + 1 + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->getPtrToSide(NB_S)->m_id << "\n";
					effD = fvm->El().fv->m_DVisc;
					b_coef = fvm->getPtrToSide(NB_S)->fv_boundary->getVelocityCoefs().x;
					b_vel = fvm->getPtrToSide(NB_S)->fv_boundary->getVelocity().x;
					//b_coef = 0.5*(FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->getPtrToSide(NB_S)->fv_boundary->getVelocityCoefs().x + FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + 1 + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->getPtrToSide(NB_S)->fv_boundary->getVelocityCoefs().x);
					//b_vel = 0.5*(FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->getPtrToSide(NB_S)->fv_boundary->getVelocity().x + FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + 1 + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->getPtrToSide(NB_S)->fv_boundary->getVelocity().x);
					boundary_coefs = 1.0 / (1.0 / (effD / vec_mod(vec_fv_to_b)) + 1.0 / b_coef)*Area;
					//std::cout << boundary_coefs << "\n";
					fvm->fv.fv_coef += boundary_coefs * sigma;
					fvm->fv_source -= fvm->fv.fv->m_velocity.dval.x * boundary_coefs * (1.0 - sigma);
					fvm->fv_source += boundary_coefs * b_vel * sigma; // *(1 - sigma)
					//std::cout << fvm->fv_source << "->";
					fvm->fv_source += boundary_coefs * b_vel * (1.0 - sigma); //old boundary velocity
					//std::cout << fvm->fv_source << "->";
				}
				*/


				fvm->fv.fv_coef /= (AlphaU);
				fvm->fv_source += (1.0 - AlphaU)*fvm->fv.fv_coef*fvm->fv.fv->m_velocity.dval.x;
			}
		}
		std::cout << "PHI_PARAMS::UX\n";
		Matrix Umat(Solve(&FVMesh_U, &UApr.dval));



		// ‘P*(aP0 + sigma*SUM(aNB) - Spm*Vol) = SUM(‘NB * aNB) + ‘P0*(aP0 - (1-sigma)*SUM(aNB) - (1-sigma)*Spm0*Vol) + (1-sigma)*SUM(‘NB0 * aNB) + (1-sigma)*Spp0*Vol + sigma*Spp*Vol + dx*(P*P - P*N)
		//      _          _             _               _              _              _                    _                      _                          _                  _              -            
		for (auto fvm : FVMesh_V) {
			fv_prec_3 Area_centre;
			fv_prec Area;
			fv_prec Volume;
			fv_prec_3 InNormal;
			switch (fvm->fv.fv->FVdim) {
			case(FVA_S):
				Volume = fvm->fv.fv->m_s.size();
				break;
			default:
				break;
			}

			fv_prec effD;
			fv_prec_3 vec_fv_to_nb;
			fv_prec_3 vec_fv_to_b;
			fv_prec_3 vec_b_to_nb;
			fvm->fv.fv_coef = 0.0;
			fvm->fv_source = 0.0;

			fv_prec boundary_coefs = 0.0;
			fv_prec b_coef = 0.0;
			fv_prec b_vel = 0.0;

			if (fvm->fv.fv->m_type == FVT_DEFAULT) {
				for (int nb = 0;nb < fvm->nb_fvs.size();nb++) {
					FV_NB relativ_pos;
					fv_prec Dnb_dir = 0.0;
					fv_prec Fnb_dir = 0.0;
					switch (fvm->fv.fv->FVdim) {
					case(FVA_S):
						for (auto fv_p : fvm->fv.fv->m_s.m_lineArray) {
							bool NBIsFound = false;

							if (fvm->nb_fvs[nb].fv->FVdim == FVA_L) {
								//std::cout << "{" << fv_p.m_p1.x << "," << fv_p.m_p1.y << "},{" << fv_p.m_p2.x << "," << fv_p.m_p2.y << "}-{" << nb.fv->m_l.m_p1.x << "," << nb.fv->m_l.m_p1.y << "},{" << nb.fv->m_l.m_p2.x << "," << nb.fv->m_l.m_p2.y << "}" << nb.fv->m_type << "\n";
								if (fv_p == fvm->nb_fvs[nb].fv->m_l) {
									Area = fv_p.size();
									Area_centre = fv_p.getCentre();
									if (Area_centre.x > fvm->El().fv->m_Centre.x) relativ_pos = NB_E;
									if (Area_centre.x < fvm->El().fv->m_Centre.x) relativ_pos = NB_W;
									if (Area_centre.y > fvm->El().fv->m_Centre.y) relativ_pos = NB_N;
									if (Area_centre.y < fvm->El().fv->m_Centre.y) relativ_pos = NB_S;
									NBIsFound = true;
								}
							}
							else if (fvm->nb_fvs[nb].fv->FVdim == FVA_S) {
								for (auto fv_nb : fvm->nb_fvs[nb].fv->m_s.m_lineArray) {
									if (fv_p == fv_nb) {
										Area = fv_p.size();
										Area_centre = fv_p.getCentre();
										if (Area_centre.x > fvm->fv.fv->m_Centre.x) relativ_pos = NB_E;
										if (Area_centre.x < fvm->fv.fv->m_Centre.x) relativ_pos = NB_W;
										if (Area_centre.y > fvm->fv.fv->m_Centre.y) relativ_pos = NB_N;
										if (Area_centre.y < fvm->fv.fv->m_Centre.y) relativ_pos = NB_S;
										NBIsFound = true;
										break;
									}
								}
							}
							if (!NBIsFound) continue;
							vec_fv_to_nb = fvm->fv.fv->m_Centre - fvm->nb_fvs[nb].fv->m_Centre;
							vec_fv_to_b = fvm->fv.fv->m_Centre - Area_centre;
							vec_b_to_nb = Area_centre - fvm->nb_fvs[nb].fv->m_Centre;
							InNormal = glm::normalize(-vec_fv_to_b);


							//std::cout << fvm->fv.fv->m_id << " (";
							//switch (relativ_pos) {
							//case(NB_E):
							//	std::cout << "NB_E";
							//	break;
							//case(NB_W):
							//	std::cout << "NB_W";
							//	break;
							//case(NB_N):
							//	std::cout << "NB_N";
							//	break;
							//case(NB_S):
							//	std::cout << "NB_S";
							//	break;
							//default:
							//	break;
							//}
							//std::cout << "): " << fvm->nb_fvs[nb].fv->m_id << " ";
							if (fvm->nb_fvs[nb].fv->m_type == FVT_DEFAULT) {
								fv_prec p_coef = fvm->nb_fvs[nb].fv->m_DVisc*vec_mod(vec_fv_to_b) / vec_mod(vec_fv_to_nb);;
								fv_prec nb_coef = fvm->fv.fv->m_DVisc*vec_mod(vec_b_to_nb) / vec_mod(vec_fv_to_nb);
								effD = fvm->nb_fvs[nb].fv->m_DVisc*fvm->fv.fv->m_DVisc / (p_coef + nb_coef);
								Dnb_dir = effD / vec_mod(vec_fv_to_nb);

								switch (relativ_pos)
								{
								case NB_E:
									Fnb_dir = 0.5*(UVolumes[fvm->fv.fv->m_id - (fvm->fv.fv->m_id) / (last_y_fvs.size()) + (last_x_fvs.size() - 1)]->m_density.dval*UVolumes[fvm->fv.fv->m_id - (fvm->fv.fv->m_id) / (last_y_fvs.size()) + (last_x_fvs.size() - 1)]->m_velocity.dval.x + UVolumes[fvm->fv.fv->m_id - (fvm->fv.fv->m_id) / (last_y_fvs.size())]->m_density.dval*UVolumes[fvm->fv.fv->m_id - (fvm->fv.fv->m_id) / (last_y_fvs.size())]->m_velocity.dval.x);
									break;
								case NB_W:
									Fnb_dir = 0.5*(UVolumes[fvm->fv.fv->m_id - (fvm->fv.fv->m_id) / (last_y_fvs.size()) - 1 + (last_x_fvs.size() - 1)]->m_density.dval*UVolumes[fvm->fv.fv->m_id - (fvm->fv.fv->m_id) / (last_y_fvs.size()) - 1 + (last_x_fvs.size() - 1)]->m_velocity.dval.x + UVolumes[fvm->fv.fv->m_id - (fvm->fv.fv->m_id) / (last_y_fvs.size()) - 1]->m_density.dval*UVolumes[fvm->fv.fv->m_id - (fvm->fv.fv->m_id) / (last_y_fvs.size()) - 1]->m_velocity.dval.x);
									break;
								case NB_N:
								case NB_S:
									Fnb_dir = 0.5*(fvm->fv.fv->m_density.dval*fvm->fv.fv->m_velocity.dval.y + fvm->nb_fvs[nb].fv->m_density.dval*fvm->nb_fvs[nb].fv->m_velocity.dval.y);
									break;
								default:
									break;
								}
								Dnb_dir *= Area;
								Fnb_dir *= Area;
								// aNB							
								{
									//std::cout << "Fnb_dir: " << Fnb_dir << " Dnb_dir:" << Dnb_dir << " - ";
									// Diffusion
									fvm->nb_fvs[nb].fv_coef = PeFunc(Fnb_dir / Dnb_dir) * Dnb_dir;
									//std::cout << fvm->nb_fvs[nb].fv_coef << " - >";
									// Convention
									//fvm->nb_fvs[nb].fv_coef += glm::max(0.0f, Fnb_dir);
									if (relativ_pos == NB_E) fvm->nb_fvs[nb].fv_coef += glm::max(-Fnb_dir, 0.0f);
									if (relativ_pos == NB_W) fvm->nb_fvs[nb].fv_coef += glm::max(Fnb_dir, 0.0f);
									if (relativ_pos == NB_N) fvm->nb_fvs[nb].fv_coef += glm::max(-Fnb_dir, 0.0f);
									if (relativ_pos == NB_S) fvm->nb_fvs[nb].fv_coef += glm::max(Fnb_dir, 0.0f);
									//std::cout << fvm->nb_fvs[nb].fv_coef << " ";
								}
								//fvm->nb_fvs[nb].fv_coef *= Area;
								//boundary_coefs = 0.0;

							}
							else if (fvm->nb_fvs[nb].fv->m_type == FVT_BOUNDARY) {

								b_coef = fvm->nb_fvs[nb].fv->fv_boundary->getVelocityCoefs().y;
								//std::cout << b_coef << " ";
								b_vel = fvm->nb_fvs[nb].fv->fv_boundary->getVelocity().y;
								//std::cout << b_vel << " ";
								effD = fvm->El().fv->m_DVisc;

								fv_prec p_coef = effD / vec_mod(vec_fv_to_b);

								//std::cout << effD << " " << vec_mod(vec_fv_to_b) << " ";
								boundary_coefs = 1.0 / (1.0 / p_coef + 1.0 / b_coef)*Area;
								//boundary_coefs = 1.0 / (1.0 / (effD / vec_mod(vec_fv_to_b)) + 1.0 / b_coef)*Area;

								//if ((relativ_pos == NB_N) or (relativ_pos == NB_S)) boundary_coefs -= 0.001;

								//std::cout << boundary_coefs;


								fvm->fv.fv_coef += boundary_coefs * sigma;
								fvm->fv_source -= fvm->fv.fv->m_velocity.dval.x * boundary_coefs * (1.0 - sigma);
								fvm->fv_source += boundary_coefs * b_vel * sigma; // *(1 - sigma)
								//std::cout << fvm->fv_source << "->";
								fvm->fv_source += boundary_coefs * b_vel * (1.0 - sigma); //old boundary velocity
								//std::cout << fvm->fv_source << "->";




							}
							//std::cout << "\n";
							//fvm->fv.fv_coef -= Fnb_dir * Area;
							if ((relativ_pos == NB_E) or (relativ_pos == NB_W)) fvm->fv.fv_coef += Fnb_dir * InNormal.x * Area;
							if ((relativ_pos == NB_N) or (relativ_pos == NB_S)) fvm->fv.fv_coef += Fnb_dir * InNormal.y * Area;

							//if (relativ_pos == NB_E) fvm->fv.fv_coef += Fnb_dir * Area;
							//if (relativ_pos == NB_W) fvm->fv.fv_coef += -Fnb_dir * Area;
							//if (relativ_pos == NB_N) fvm->fv.fv_coef += Fnb_dir * Area;
							//if (relativ_pos == NB_S) fvm->fv.fv_coef += -Fnb_dir * Area;
						}
						break;
					default:
						break;
					}
					// sigma*SUM(aNB)
					fvm->fv.fv_coef += fvm->nb_fvs[nb].fv_coef * sigma;
					//
					// ‘P0*(1-sigma)*SUM(aNB)
					fvm->fv_source -= fvm->fv.fv->m_velocity.dval.y * fvm->nb_fvs[nb].fv_coef * (1.0 - sigma); // *(1 - sigma)
					// (1-sigma)*SUM(‘NB0 * aNB)
					fvm->fv_source += fvm->nb_fvs[nb].fv->m_velocity.dval.y * fvm->nb_fvs[nb].fv_coef * (1.0 - sigma); // *(1 - sigma)




				}
				// Pressure
				for (auto Pcor_fv_ls : FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->FVPtr()->m_s.m_lineArray) {
					for (auto Pcor_nb_ls : FVM.FiniteVolumeMesh[fvm->fv.fv->m_id + (fvm->fv.fv->m_id) / (last_y_fvs.size() - 1)]->getPtrToSide(NB_E)->m_s.m_lineArray) {
						if (Pcor_fv_ls == Pcor_nb_ls) {
							Area = Pcor_fv_ls.size();
						}
					}
				}
				// dx*(P*P - P*N)
				fvm->fv_source += Area * (PApr[fvm->fv.fv->m_id] - PApr[fvm->fv.fv->m_id + last_x_fvs.size()]);
				//std::cout << fvm->fv.fv->m_id << " -{" << fvm->fv.fv->m_id  << ", " << fvm->fv.fv->m_id + last_x_fvs.size() << "}\n";

				// Spm*Vol*sigma
				fvm->fv.fv_coef -= fvm->sf.Spm.val * Volume*sigma;
				// sigma*Spp*Vol
				fvm->fv_source += fvm->sf.Spp.val * Volume*sigma;
				// ‘P0*(1-sigma)*Spm0*Vol
				fvm->fv_source -= fvm->fv.fv->m_velocity.val.y * fvm->sf.Spm.dval * Volume*(1.0 - sigma);
				// (1-sigma)*Spp0*Vol
				fvm->fv_source += fvm->sf.Spp.dval * Volume*(1.0 - sigma);


				if (dt == 0) {}
				else {
					// aP0
					fvm->fv.fv_coef += fvm->fv.fv->m_density.dval * Volume / dt;
					// ‘P0*aP0
					fvm->fv_source += fvm->fv.fv->m_velocity.val.y* fvm->fv.fv->m_density.val * Volume / dt;
				}

				fvm->fv.fv_coef /= (AlphaV);
				fvm->fv_source += (1.0 - AlphaV)*fvm->fv.fv_coef*fvm->fv.fv->m_velocity.dval.y;
			}
		}
		std::cout << "PHI_PARAMS::UY\n";
		Matrix Vmat(Solve(&FVMesh_V, &VApr.dval));


		for (int fv = 0; fv < UApr.dval.size(); fv++) {
			FVMesh_U[fv]->fv.fv->m_velocity.dval.x = UApr.dval[fv];
		}
		for (auto fvm : FVMesh_U) {
			fvm->RecalcResidual();
			Uresidual += fvm->residual;
			if (MaxUresidual < abs(fvm->residual)) MaxUresidual = abs(fvm->residual);
			for (auto fvnb : fvm->nb_fvs) {
				residualScaleU += abs(fvnb.fv_coef*(fvm->fv.fv->m_velocity.dval.x - fvnb.fv->m_velocity.dval.x));
			}

		}
		std::cout << "Uresidual= " << Uresidual << "\n";
		std::cout << "MaxUresidual= " << MaxUresidual << "\n";
		std::cout << "residualScaleU= " << residualScaleU << "\n";

		for (int fv = 0; fv < VApr.dval.size(); fv++) {
			FVMesh_V[fv]->fv.fv->m_velocity.dval.y = VApr.dval[fv];
		}

		for (auto fvm : FVMesh_V) {
			fvm->RecalcResidual();
			Vresidual += fvm->residual;
			if (MaxVresidual < abs(fvm->residual)) MaxVresidual = abs(fvm->residual);
			for (auto fvnb : fvm->nb_fvs) {
				residualScaleV += abs(fvnb.fv_coef*(fvm->fv.fv->m_velocity.dval.y - fvnb.fv->m_velocity.dval.y));
			}
		}
		std::cout << "Vresidual= " << Vresidual << "\n";
		std::cout << "MaxVresidual= " << MaxVresidual << "\n";
		std::cout << "residualScaleV= " << residualScaleV << "\n";


		for (auto fvm : FVMesh_PC) {
			fv_prec_3 Area_centre;
			fv_prec Area;
			fv_prec Volume;
			fvm->fv.fv_coef = 0.0;
			fv_prec de = 0.0; fv_prec dw = 0.0; fv_prec dn = 0.0; fv_prec ds = 0.0;
			fv_prec ue = 0.0; fv_prec uw = 0.0; fv_prec vn = 0.0; fv_prec vs = 0.0;
			for (int nb = 0;nb < fvm->nb_fvs.size();nb++) {
				FV_NB relativ_pos;
				//boundary_coefs = 0.0;
				fv_prec dens = fvm->El().fv->m_density.val;
				//std::cout << dens << "\n";
				switch (fvm->El().fv->FVdim) {
				case(FVA_L):
					break;
				case(FVA_S):
					Volume = fvm->El().fv->m_s.size();
					for (auto fv_p : fvm->El().fv->m_s.m_lineArray) {
						bool NBIsFound = false;
						if (fvm->nb_fvs[nb].fv->FVdim == FVA_L) {
							//std::cout << "{" << fv_p.m_p1.x << "," << fv_p.m_p1.y << "},{" << fv_p.m_p2.x << "," << fv_p.m_p2.y << "}-{" << nb.fv->m_l.m_p1.x << "," << nb.fv->m_l.m_p1.y << "},{" << nb.fv->m_l.m_p2.x << "," << nb.fv->m_l.m_p2.y << "}" << nb.fv->m_type << "\n";
							if (fv_p == fvm->nb_fvs[nb].fv->m_l) {
								Area = fv_p.size();
								Area_centre = fv_p.getCentre();
								if (Area_centre.x > fvm->El().fv->m_Centre.x) relativ_pos = NB_E;
								if (Area_centre.x < fvm->El().fv->m_Centre.x) relativ_pos = NB_W;
								if (Area_centre.y > fvm->El().fv->m_Centre.y) relativ_pos = NB_N;
								if (Area_centre.y < fvm->El().fv->m_Centre.y) relativ_pos = NB_S;
								NBIsFound = true;
							}
						}
						else if (fvm->nb_fvs[nb].fv->FVdim == FVA_S) {
							for (auto fv_nb : fvm->nb_fvs[nb].fv->m_s.m_lineArray) {
								if (fv_p == fv_nb) {
									Area = fv_p.size();
									Area_centre = fv_p.getCentre();
									if (Area_centre.x > fvm->El().fv->m_Centre.x) relativ_pos = NB_E;
									if (Area_centre.x < fvm->El().fv->m_Centre.x) relativ_pos = NB_W;
									if (Area_centre.y > fvm->El().fv->m_Centre.y) relativ_pos = NB_N;
									if (Area_centre.y < fvm->El().fv->m_Centre.y) relativ_pos = NB_S;
									NBIsFound = true;
									//std::cout << "{" << fv_p.m_p1.x << "," << fv_p.m_p1.y << "},{" << fv_p.m_p2.x << "," << fv_p.m_p2.y << "}-{" << fv_nb.m_p1.x << "," << fv_nb.m_p1.y << "},{" << fv_nb.m_p2.x << "," << fv_nb.m_p2.y << "}" << "\n";
									//std::cout << fvm->fv.fv->m_id << " " <<  fvm->nb_fvs[nb].fv->m_id << " TUT\n";
									break;
								}
							}
						}
						if (!NBIsFound) continue;

						//std::cout << fvm->fv.fv->m_id << " (";
						//switch (relativ_pos) {
						//case(NB_E):
						//	std::cout << "NB_E";
						//	break;
						//case(NB_W):
						//	std::cout << "NB_W";
						//	break;
						//case(NB_N):
						//	std::cout << "NB_N";
						//	break;
						//case(NB_S):
						//	std::cout << "NB_S";
						//	break;
						//default:
						//	break;
						//}
						//std::cout << "): ";
						if (relativ_pos == NB_E) {
							if (fvm->nb_fvs[nb].fv->m_type == FINITE_VOLUME_TYPE::FVT_BOUNDARY) {
								//std::cout << fvm->nb_fvs[nb].fv->m_id << "\n";
								ue = fvm->nb_fvs[nb].fv->fv_boundary->getVelocity().x;
								//de = Area / FVMesh_U[]->fv.fv_coef;
								de = 0.0;
							}
							else {
								//std::cout << fvm->fv.fv->m_id - (fvm->fv.fv->m_id) / (last_y_fvs.size()) << "\n";
								ue = UVolumes[fvm->fv.fv->m_id - (fvm->fv.fv->m_id) / (last_y_fvs.size())]->m_velocity.dval.x;
								de = Area / FVMesh_U[fvm->fv.fv->m_id - (fvm->fv.fv->m_id) / (last_y_fvs.size())]->fv.fv_coef;
								dens = 0.5*(fvm->fv.fv->m_density.dval + fvm->nb_fvs[nb].fv->m_density.dval);
								fvm->nb_fvs[nb].fv_coef = dens * Area*de;
								//std::cout << FVMesh_U[fvm->fv.fv->m_id - (fvm->fv.fv->m_id) / (last_y_fvs.size())]->fv.fv_coef << " ";
							}
							//std::cout << dens * Area*de << " ";
						}
						else if (relativ_pos == NB_W) {
							if (fvm->nb_fvs[nb].fv->m_type == FINITE_VOLUME_TYPE::FVT_BOUNDARY) {
								//std::cout << fvm->nb_fvs[nb].fv->m_id << "\n";
								uw = fvm->nb_fvs[nb].fv->fv_boundary->getVelocity().x;
								dw = 0.0;}
							else {
								//std::cout << fvm->fv.fv->m_id - (fvm->fv.fv->m_id) / (last_y_fvs.size()) - 1 << "\n";
								uw = UVolumes[fvm->fv.fv->m_id - (fvm->fv.fv->m_id) / (last_y_fvs.size()) - 1]->m_velocity.dval.x;
								dw = Area / FVMesh_U[fvm->fv.fv->m_id - (fvm->fv.fv->m_id) / (last_y_fvs.size()) - 1]->fv.fv_coef;
								dens = 0.5*(fvm->fv.fv->m_density.dval + fvm->nb_fvs[nb].fv->m_density.dval);
								fvm->nb_fvs[nb].fv_coef = dens * Area*dw;
								//std::cout << FVMesh_U[fvm->fv.fv->m_id - (fvm->fv.fv->m_id) / (last_y_fvs.size()) - 1]->fv.fv_coef << " ";
							}
							//std::cout << dens * Area*dw << " ";
						}
						else if (relativ_pos == NB_N) {
							if (fvm->nb_fvs[nb].fv->m_type == FINITE_VOLUME_TYPE::FVT_BOUNDARY) {
								//std::cout << fvm->nb_fvs[nb].fv->m_id << "\n";
								vn = fvm->nb_fvs[nb].fv->fv_boundary->getVelocity().y;
								dn = 0.0;
							}
							else {
								//std::cout << fvm->fv.fv->m_id << "\n";
								vn = VVolumes[fvm->fv.fv->m_id]->m_velocity.dval.y;
								dn = Area / FVMesh_V[fvm->fv.fv->m_id]->fv.fv_coef;
								dens = 0.5*(fvm->fv.fv->m_density.dval + fvm->nb_fvs[nb].fv->m_density.dval);
								fvm->nb_fvs[nb].fv_coef = dens * Area*dn;
								//std::cout << FVMesh_V[fvm->fv.fv->m_id]->fv.fv_coef << " ";
							}
							//std::cout << dens * Area*dn << " ";
						}
						else if (relativ_pos == NB_S) {
							if (fvm->nb_fvs[nb].fv->m_type == FINITE_VOLUME_TYPE::FVT_BOUNDARY) {
								//std::cout << fvm->nb_fvs[nb].fv->m_id << "\n";
								vs = fvm->nb_fvs[nb].fv->fv_boundary->getVelocity().y;
								ds = 0.0;
							}
							else {
								//std::cout << fvm->fv.fv->m_id - last_x_fvs.size() << "\n";
								vs = VVolumes[fvm->fv.fv->m_id - last_x_fvs.size()]->m_velocity.dval.y;
								ds = Area / FVMesh_V[fvm->fv.fv->m_id - last_x_fvs.size()]->fv.fv_coef;
								dens = 0.5*(fvm->fv.fv->m_density.dval + fvm->nb_fvs[nb].fv->m_density.dval);
								fvm->nb_fvs[nb].fv_coef = dens * Area*ds;
								//std::cout << FVMesh_V[fvm->fv.fv->m_id - last_x_fvs.size()]->fv.fv_coef << " ";
							}
							//std::cout << dens * Area*ds << " ";
						}
						//std::cout << "\n";
					}
					break;
				default:
					break;
				}
				//if (relativ_pos == NB_E) { std::cout << "de=" << de << " - " << "ue=" << ue << "\n"; }
				//else if (relativ_pos == NB_W) { std::cout << "dw=" << dw << " - " << "uw=" << uw << "\n"; }
				//else if (relativ_pos == NB_N) { std::cout << "dn=" << dn << " - " << "vn=" << vn << "\n"; }
				//else if (relativ_pos == NB_S) { std::cout << "ds=" << ds << " - " << "vs=" << vs << "\n"; }

				fvm->fv.fv_coef += fvm->nb_fvs[nb].fv_coef;
				fvm->fv_source += dens * Area * (uw - ue + vs - vn);

				de = 0.0;  dw = 0.0;  dn = 0.0;  ds = 0.0;
				ue = 0.0;  uw = 0.0;  vn = 0.0;  vs = 0.0;
			}

			if (dt == 0) { }
			else { fvm->fv_source -= (fvm->fv.fv->m_density.dval - fvm->fv.fv->m_density.val) * Volume / dt; }
		}
		std::cout << "PHI_PARAMS::P\n";

		for (auto err : FVMesh_PC) {
			Presidual += abs(err->fv_source);
			if (MaxPresidual < abs(err->fv_source)) MaxPresidual = abs(err->fv_source);
		}

		SolveWithFirstZero(&FVMesh_PC, &Pcor);


		fv_prec Area;
		// ue = ue* +de(PP'- PE')
		//std::cout << "UCor_ed: ";
		for (auto fvmU : FVMesh_U) {
			if (fvmU->fv.fv->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
				for (auto fv_m : FVM.FiniteVolumeMesh[fvmU->fv.fv->m_id + (fvmU->fv.fv->m_id) / (last_y_fvs.size() - 1)]->FVPtr()->m_s.m_lineArray) {
					for (auto fv_nb : FVM.FiniteVolumeMesh[fvmU->fv.fv->m_id + (fvmU->fv.fv->m_id) / (last_y_fvs.size() - 1)]->getPtrToSide(NB_E)->m_s.m_lineArray) {
						if (fv_m == fv_nb) {
							Area = fv_m.size();
						}
					}
				}
				//std::cout << fvmU->fv.fv->m_id << " - " << fvmU->fv.fv->m_id + (fvmU->fv.fv->m_id) / (last_y_fvs.size() - 1) << " " << fvmU->fv.fv->m_id + 1 + (fvmU->fv.fv->m_id) / (last_y_fvs.size() - 1) << "\n";
				fvmU->fv.fv->m_velocity.dval.x = fvmU->fv.fv->m_velocity.dval.x + Area/ fvmU->fv.fv_coef*(Pcor[fvmU->fv.fv->m_id + (fvmU->fv.fv->m_id) / (last_y_fvs.size() - 1)] - Pcor[fvmU->fv.fv->m_id + 1 + (fvmU->fv.fv->m_id) / (last_y_fvs.size() - 1)]);
				//std::cout << fvmU->fv.fv->m_velocity.dval.x << " ";
			}
		}
		//std::cout << "\n";

		// vn = vn* +dn(PP'- PN')
		//std::cout << "VCor_ed: ";
		for (auto fvmV : FVMesh_V) {
			if (fvmV->fv.fv->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
				for (auto fv_m : FVM.FiniteVolumeMesh[fvmV->fv.fv->m_id]->FVPtr()->m_s.m_lineArray) {
					for (auto fv_nb : FVM.FiniteVolumeMesh[fvmV->fv.fv->m_id]->getPtrToSide(NB_N)->m_s.m_lineArray) {
						if (fv_m == fv_nb) {
							Area = fv_m.size();
						}
					}
				}
				//std::cout << fvmV->fv.fv->m_id << " - " << fvmV->fv.fv->m_id << " " << fvmV->fv.fv->m_id + last_y_fvs.size() << "\n";
				fvmV->fv.fv->m_velocity.dval.y = fvmV->fv.fv->m_velocity.dval.y + Area / fvmV->fv.fv_coef*(Pcor[fvmV->fv.fv->m_id] - Pcor[fvmV->fv.fv->m_id + last_y_fvs.size()]);
				//std::cout << fvmV->fv.fv->m_velocity.dval.y << " ";
			}
		}
		//std::cout << "\n";

		/*
		//ÒÍÓÓÒÚË
		std::cout << "\n";
		std::cout << "Velocity.x: ";
		for (auto fvm : FVM.FiniteVolumeMesh) {
			if (fvm->FVPtr()->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
				fv_prec velE = 0.0;
				fv_prec velW = 0.0;
				//std::cout << (fvm->FVPtr()->m_id) / (last_x_fvs.size()) << " ";
				//std::cout << fvm->FVPtr()->m_id << ": ";
				if (fvm->getPtrToSide(NB_W)->fv_boundary != nullptr) {
					velW = fvm->getPtrToSide(NB_W)->fv_boundary->getVelocity().x;
					//std::cout << fvm->getPtrToSide(NB_W)->m_id << " ";
				}
				else {
					velW = UCor_ed[fvm->getPtrToSide(NB_W)->m_id - (fvm->FVPtr()->m_id) / (last_x_fvs.size())];
					//std::cout << fvm->getPtrToSide(NB_W)->m_id- (fvm->FVPtr()->m_id) / (last_x_fvs.size()) << " ";
				}
				if (fvm->getPtrToSide(NB_E)->fv_boundary != nullptr) {
					velE = fvm->getPtrToSide(NB_E)->fv_boundary->getVelocity().x;
					//std::cout << fvm->getPtrToSide(NB_E)->m_id << " ";
				}
				else {
					velE = UCor_ed[fvm->FVPtr()->m_id - (fvm->FVPtr()->m_id) / (last_x_fvs.size())];
					//std::cout << fvm->FVPtr()->m_id - (fvm->FVPtr()->m_id) / (last_x_fvs.size()) << " ";
				}
				fvm->FVPtr()->m_velocity.dval.x = 0.5*(velW + velE);
				//std::cout << "\n";
				std::cout << fvm->FVPtr()->m_velocity.dval.x << " ";
			}
		}
		std::cout << "\n";
		std::cout << "Velocity.y: ";
		for (auto fvm : FVM.FiniteVolumeMesh) {
			if (fvm->FVPtr()->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
				fv_prec velN = 0.0;
				fv_prec velS = 0.0;
				//std::cout << (fvm->FVPtr()->m_id) / (last_x_fvs.size()) << " ";
				//std::cout << fvm->FVPtr()->m_id << ": ";
				if (fvm->getPtrToSide(NB_S)->fv_boundary != nullptr) {
					velS = fvm->getPtrToSide(NB_S)->fv_boundary->getVelocity().y;
					//std::cout << fvm->getPtrToSide(NB_S)->m_id << " ";
				}
				else {
					velS = VCor_ed[fvm->getPtrToSide(NB_S)->m_id];
					//std::cout << fvm->getPtrToSide(NB_S)->m_id << " ";
				}
				if (fvm->getPtrToSide(NB_N)->fv_boundary != nullptr) {
					velN = fvm->getPtrToSide(NB_N)->fv_boundary->getVelocity().y;
					//std::cout << fvm->getPtrToSide(NB_N)->m_id << " ";
				}
				else {
					velN = VCor_ed[fvm->FVPtr()->m_id];
					//std::cout << fvm->FVPtr()->m_id << " ";
				}
				fvm->FVPtr()->m_velocity.dval.y = 0.5*(velS + velN);
				std::cout << fvm->FVPtr()->m_velocity.dval.y << " ";
				//std::cout << "\n";
				//std::cout << velS << " " << velN << "\n";
			}
		}
		std::cout << "\n";
		*/



		// P = P* + P'
		//std::cout << "Pressure: ";
		for (auto i : FVM.FiniteVolumes) {
			if (i->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
				//i->m_pressure.val = i->m_pressure.dval;
				i->m_pressure.dval = PApr[i->m_id] + AlphaP * Pcor[i->m_id];
				//std::cout << i->m_pressure.dval << " ";
			}
		}
		for(int i =0; i < FVM.FiniteVolumeMesh.size(); i++) {
			if (FVM.FiniteVolumeMesh[i]->FVPtr()->m_type == FINITE_VOLUME_TYPE::FVT_BOUNDARY) {
				int m_id = FVM.FiniteVolumeMesh[i]->FVPtr()->m_id;
				int nb_id = FVM.FiniteVolumeMesh[i]->NB_FVSPrt()[0]->m_id;
				FV_NB nb_relative_pos;
				if (FVM.FiniteVolumeMesh[m_id]->FVPtr()->m_Centre.x < FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_Centre.x) nb_relative_pos = NB_E;
				if (FVM.FiniteVolumeMesh[m_id]->FVPtr()->m_Centre.x > FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_Centre.x) nb_relative_pos = NB_W;
				if (FVM.FiniteVolumeMesh[m_id]->FVPtr()->m_Centre.y < FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_Centre.y) nb_relative_pos = NB_N;
				if (FVM.FiniteVolumeMesh[m_id]->FVPtr()->m_Centre.y > FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_Centre.y) nb_relative_pos = NB_S;

				int nb_nb_id = FVM.FiniteVolumeMesh[nb_id]->getPtrToSide(nb_relative_pos)->m_id;
				fv_prec diff31 = FVM.FiniteVolumeMesh[nb_nb_id]->FVPtr()->distance(FVM.FiniteVolumeMesh[i]->FVPtr());
				fv_prec diff21 = FVM.FiniteVolumeMesh[nb_id]->FVPtr()->distance(FVM.FiniteVolumeMesh[i]->FVPtr());
				fv_prec diff32= FVM.FiniteVolumeMesh[nb_nb_id]->FVPtr()->distance(FVM.FiniteVolumeMesh[nb_id]->FVPtr());

				FVM.FiniteVolumeMesh[i]->FVPtr()->m_pressure.dval = (FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_pressure.dval*(diff31) - FVM.FiniteVolumeMesh[nb_nb_id]->FVPtr()->m_pressure.dval*(diff21))/(diff32);
			}
		}


		for (auto fvm : FVMesh_PC) {
			for (auto fvnb : fvm->nb_fvs) {
				residualScale += abs(fvnb.fv_coef*(fvm->fv.fv->m_pressure.dval - fvnb.fv->m_pressure.dval));
			}
		}
		std::cout << "\n";
		std::cout << "residualScale= " << residualScale << "\n";

		// P* = P
		for (auto i : FVM.FiniteVolumes) {
			if (i->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
				PApr[i->m_id] = i->m_pressure.dval;
			}
		}
		
		//std::cout << "U velocity\n";
		//std::cout << "Max residual: " << MaxUresidual / residualScaleU * UVolumes.size() << "\n";
		//std::cout << "Average residual: " << Uresidual / residualScaleU << "\n";
		//
		//std::cout << "V velocity\n";
		//std::cout << "Max residual: " << MaxVresidual / residualScaleV * VVolumes.size() << "\n";
		//std::cout << "Average residual: " << Vresidual / residualScaleV << "\n";
		//
		//std::cout << "Pressure\n";
		//std::cout << "Max residual: " << MaxPresidual / residualScale * m_options.nrOfFV[FINITE_VOLUME_TYPE::FVT_DEFAULT] << "\n";
		//std::cout << "Average residual: " << Presidual / residualScale << "\n";
		//
		//std::cout << "_________________________________________________ \n";


		fv_prec MaxResidual = 0.0;
		fv_prec AverageResidual = 0.0;

		if (residualScaleU != 0.0) {
			if (AverageResidual < Uresidual / residualScaleU) {
				AverageResidual = Uresidual / residualScaleU;
				MaxResidual = MaxUresidual / residualScaleU * UVolumes.size();
			}
		}
		if (residualScaleV != 0.0) {
			if (AverageResidual < Vresidual / residualScaleV) {
				AverageResidual = Vresidual / residualScaleV;
				MaxResidual = MaxVresidual / residualScaleV * VVolumes.size();
			}
		}
		if (residualScale != 0.0) {
			if (AverageResidual < Presidual / residualScale) {
				AverageResidual = Presidual / residualScale;
				MaxResidual = MaxPresidual / residualScale * m_options.nrOfFV[FINITE_VOLUME_TYPE::FVT_DEFAULT];
			}
		}








		std::cout << "Average residual: " << AverageResidual << ", MaxResidual: " << MaxResidual << "\n";
		std::cout << "_________________________________________________ \n";
		// (FVMesh_PC[i]->fv.fv->m_density.dval)_Recalculation

		if (AverageResidual < Eps) {
			whileCondition = false;
		}

		AlgoIteration++;
	}

	//std::cout << "velocity.x ";
	//std::cout << "velocity.y ";
	for (auto fvm : FVM.FiniteVolumeMesh) {
		if(fvm->FVPtr()->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT){
			if (fvm->getPtrToSide(NB_E)->fv_boundary != nullptr) {
				//std::cout << "NB_E ";
				//std::cout << fvm->FVPtr()->m_id << " " << fvm->FVPtr()->m_id - (fvm->FVPtr()->m_id) / (last_y_fvs.size()) - 1 << " " << fvm->getPtrToSide(NB_E)->m_id << "\n";
				fvm->FVPtr()->m_velocity.val.x = 0.5*(FVMesh_U[fvm->FVPtr()->m_id - (fvm->FVPtr()->m_id) / (last_y_fvs.size()) - 1]->fv.fv->m_velocity.dval.x + fvm->getPtrToSide(NB_E)->fv_boundary->getVelocity().x);
			}
			else if (fvm->getPtrToSide(NB_W)->fv_boundary != nullptr) {
				//std::cout << "NB_W ";
				//std::cout << fvm->FVPtr()->m_id << " " << fvm->getPtrToSide(NB_W)->m_id << " " << fvm->FVPtr()->m_id - (fvm->FVPtr()->m_id) / (last_y_fvs.size()) << "\n";
				fvm->FVPtr()->m_velocity.val.x = 0.5*(fvm->getPtrToSide(NB_W)->fv_boundary->getVelocity().x + FVMesh_U[fvm->FVPtr()->m_id - (fvm->FVPtr()->m_id) / (last_y_fvs.size())]->fv.fv->m_velocity.dval.x);

			}
			else {
				//std::cout << fvm->FVPtr()->m_id << " " << fvm->FVPtr()->m_id - (fvm->FVPtr()->m_id) / (last_y_fvs.size()) - 1 << " " << fvm->FVPtr()->m_id - (fvm->FVPtr()->m_id) / (last_y_fvs.size()) << "\n";
				fvm->FVPtr()->m_velocity.val.x = 0.5*(FVMesh_U[fvm->FVPtr()->m_id - (fvm->FVPtr()->m_id) / (last_y_fvs.size()) - 1]->fv.fv->m_velocity.dval.x + FVMesh_U[fvm->FVPtr()->m_id - (fvm->FVPtr()->m_id) / (last_y_fvs.size())]->fv.fv->m_velocity.dval.x);

			}
			//std::cout << fvm->FVPtr()->m_velocity.val.x << " ";


			if (fvm->getPtrToSide(NB_N)->fv_boundary != nullptr) {
				//std::cout << "NB_N ";
				//std::cout << fvm->FVPtr()->m_id << " " << fvm->FVPtr()->m_id - last_y_fvs.size() << " " << fvm->getPtrToSide(NB_N)->m_id << "\n";
				fvm->FVPtr()->m_velocity.val.y = 0.5*(VVolumes[fvm->FVPtr()->m_id - last_y_fvs.size()]->m_velocity.dval.y + fvm->getPtrToSide(NB_N)->fv_boundary->getVelocity().y);

			}
			else if (fvm->getPtrToSide(NB_S)->fv_boundary != nullptr) {
				//std::cout << "NB_S ";
				//std::cout << fvm->FVPtr()->m_id << " " << fvm->getPtrToSide(NB_S)->m_id << " " << fvm->FVPtr()->m_id << "\n";
				fvm->FVPtr()->m_velocity.val.y = 0.5*(fvm->getPtrToSide(NB_S)->fv_boundary->getVelocity().y + UVolumes[fvm->FVPtr()->m_id]->m_velocity.dval.y);

			}
			else {
				//std::cout << fvm->FVPtr()->m_id << " " << fvm->FVPtr()->m_id - last_y_fvs.size() << " " << fvm->FVPtr()->m_id << "\n";
				fvm->FVPtr()->m_velocity.val.y = 0.5*(VVolumes[fvm->FVPtr()->m_id - last_y_fvs.size()]->m_velocity.dval.y + VVolumes[fvm->FVPtr()->m_id]->m_velocity.dval.y);

			}
			//std::cout << fvm->FVPtr()->m_velocity.val.y << " ";
		}

	}
	//std::cout << "\n";

	for (auto i : FVM.FiniteVolumes) {
		if (i->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
			i->m_pressure.val = PApr[i->m_id];
		}
	}


	for (auto i : FVMesh_U) { delete i; }
	FVMesh_U.resize(0);
	for (auto i : FVMesh_V) { delete i; }
	FVMesh_V.resize(0);
	for (auto i : FVMesh_PC) { delete i; }
	FVMesh_PC.resize(0);


	/*
	std::vector<fv_prec> PApr;
	PApr.resize(FVM.FiniteVolumes.size(), 0.0);
	std::vector<fv_prec> Pcor;
	Pcor.resize(FVM.FiniteVolumes.size(), 0.0);


	fv_prec sigma = 1.0;



	fv_prec AlphaP = 0.8;
	fv_prec AlphaU = 0.5;
	fv_prec AlphaV = 0.5;
	fv_prec AlphaW = 0.5;

	std::vector<FVMeshNode*> FVMesh_VxA;
	std::vector<FVMeshNode*> FVMesh_VyA;
	std::vector<FVMeshNode*> FVMesh_PC;

	for (auto i : FVM.FiniteVolumes) {
		//i->m_pressure.dval = i->m_pressure.val;
		PApr[i->m_id] = i->m_pressure.val;
		i->m_pressure.dval = i->m_pressure.val;
		i->m_velocity.dval = i->m_velocity.val;
		//  ÕŒ Õ¿ƒŒ œ≈–≈—“–Œ»“‹ «Õ¿◊≈Õ»ﬂ Œ“ ÕŒ–Ã¿À‹Õ€’  Œ   —Ã≈Ÿ≈ÕÕ€Ã ƒÀﬂ — Œ–Œ—“», » “Œ√ƒ¿ ”Ã≈Õ‹ÿ»“‹ –¿«Ã≈–€ UCor_ed, VCor_ed, WCor_ed
		//  Õ” “Œ√ƒ¿ œŒÃ≈Õﬂ“‹ » –¿«Ã≈–€ UApr, VApr, WApr
	}


	std::vector<int> last_x_fvs;
	std::vector<int> last_y_fvs;
	fv_prec max_x = 0.0; fv_prec max_y = 0.0; fv_prec max_z = 0.0;
	for (int i = 0;i < FVM.FiniteVolumeMesh.size();i++) {
		if (max_x < FVM.FiniteVolumeMesh[i]->FVPtr()->m_Centre.x) max_x = FVM.FiniteVolumeMesh[i]->FVPtr()->m_Centre.x;
		if (max_y < FVM.FiniteVolumeMesh[i]->FVPtr()->m_Centre.y) max_y = FVM.FiniteVolumeMesh[i]->FVPtr()->m_Centre.y;
		if (max_z < FVM.FiniteVolumeMesh[i]->FVPtr()->m_Centre.z) max_z = FVM.FiniteVolumeMesh[i]->FVPtr()->m_Centre.z;
	}
	max_x = 1.00 - m_options.average_dim_steps.x;
	max_y = 1.00 - m_options.average_dim_steps.y;
	for (int i = 0;i < FVM.FiniteVolumeMesh.size();i++) {
		if (FVM.FiniteVolumeMesh[i]->FVPtr()->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
			if (FVM.FiniteVolumeMesh[i]->FVPtr()->m_Centre.x >= max_x) { last_x_fvs.push_back(i); }
			if (FVM.FiniteVolumeMesh[i]->FVPtr()->m_Centre.y >= max_y) { last_y_fvs.push_back(i); }
		}
	}

	int num_of_def = 0;
	for (int i = 0;i < FVM.FiniteVolumeMesh.size();i++) {
		if (FVM.FiniteVolumeMesh[i]->FVPtr()->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
			num_of_def++;
		}
	}


	Ch_<std::vector<fv_prec>> UApr;
	UApr.val.resize(num_of_def - last_x_fvs.size(), 0.0);
	UApr.dval.resize(num_of_def - last_x_fvs.size(), 0.0);
	std::vector<fv_prec> UCor_ed;
	UCor_ed.resize(num_of_def - last_x_fvs.size(), 0.0);
	std::vector<fv_prec> UCor0;
	UCor0.resize(num_of_def - last_x_fvs.size(), 0.0);

	Ch_<std::vector<fv_prec>> VApr;
	VApr.val.resize(num_of_def - last_y_fvs.size(), 0.0);
	VApr.dval.resize(num_of_def - last_y_fvs.size(), 0.0);
	std::vector<fv_prec> VCor_ed;
	VCor_ed.resize(num_of_def - last_y_fvs.size(), 0.0);
	std::vector<fv_prec> VCor0;
	VCor0.resize(num_of_def - last_y_fvs.size(), 0.0);


	std::vector<fv_prec> WApr;
	WApr.resize(FVM.FiniteVolumes.size(), 0.0);
	std::vector<fv_prec> WCor_ed;
	WCor_ed.resize(FVM.FiniteVolumes.size(), 0.0);


	std::cout << "UCor0[n]: ";
	int n = 0;
	for (auto fvm : FVM.FiniteVolumeMesh) {
		if (((fvm->FVPtr()->m_id + 1) % last_y_fvs.size() != 0) and (n < num_of_def - last_x_fvs.size())) {
			UCor0[n] = 0.5*(fvm->FVPtr()->m_velocity.val.x + fvm->getPtrToSide(NB_E)->m_velocity.val.x);
			UCor_ed[n] = 0.5*(fvm->FVPtr()->m_velocity.dval.x + fvm->getPtrToSide(NB_E)->m_velocity.dval.x);
			std::cout << UCor0[n] << " ";
			n++;
		}
	}
	std::cout << "\n";


	std::cout << "VCor0[n]: ";
	n = 0;
	for (auto fvm : FVM.FiniteVolumeMesh) {
		if (n < num_of_def - last_y_fvs.size()) {
			VCor0[n] = 0.5*(fvm->FVPtr()->m_velocity.val.y + fvm->getPtrToSide(NB_N)->m_velocity.val.y);
			VCor_ed[n] = 0.5*(fvm->FVPtr()->m_velocity.dval.y + fvm->getPtrToSide(NB_N)->m_velocity.dval.y);
			std::cout << VCor0[n] << " ";
			n++;
		}
	}
	std::cout << "\n";

	//for (int i = 0; i < VCor0.size();i++) {	std::cout << "{" << UCor_ed[i] << ", " << VCor_ed[i] << "}\n"; }

	//n = 0;
	//for (auto fvm : FVM.FiniteVolumeMesh) {
	//	if (n < num_of_def - last_x_fvs.size()) {
	//		VCor_ed[n] = fvm->FVPtr()->m_velocity.val.x + fvm->getPtrToSide(NB_N)->m_velocity.val.x;
	//		n++;
	//	}
	//}

	std::cout << "Err: " << abs(sum.dval - sum.val) << "\n";
	while (abs(sum.dval - sum.val) > 0.00001) {


		// ae ue* = SUM(anb unb*) + bu + dy(PP*-PE*)
		for (int i = 0;i < FVM.FiniteVolumeMesh.size();i++) {
			if (FVM.FiniteVolumeMesh[i]->FVPtr()->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
				float skip_flag = false;
				for (auto x : last_x_fvs) { if (x == i) { skip_flag = true; } }
				if (skip_flag) { continue; }
				FVMesh_VxA.push_back(new FVMeshNode(*(FVM.FiniteVolumeMesh[i]), PHI_PARAMS::UX));
			}
		}
		fv_prec effD;
		fv_prec_3 vec_fv_to_nb;
		fv_prec_3 vec_fv_to_b;
		fv_prec_3 vec_b_to_nb;
		fv_prec Area;
		fv_prec Volume;
		fv_prec_3 Area_centre;
		fv_prec boundary_coefs;
		fv_prec Dnb_dir = 0.0;
		fv_prec Fnb_dir = 0.0;
		fv_prec FP0;
		fv_prec FNB0;
		fv_prec Dens = 1000.0;
		fv_prec Dens0;
		fv_prec DensNB;

		// ‘P*(aP0 + sigma*SUM(aNB) - Spm*Vol) = SUM(‘NB * aNB) + ‘P0*(aP0 - (1-sigma)*SUM(aNB) - (1-sigma)*Spm0*Vol) + (1-sigma)*SUM(‘NB0 * aNB) + (1-sigma)*Spp0*Vol + sigma*Spp*Vol + dy*(P*P - P*E)
		//      +          +             +               +              +              +                    +                      +                                                            +            
		glm::max(0.0, 0.1);

		int main_innerid = 0;
		for (auto fvm : FVMesh_VxA) {
			fv_prec ap0;
			fvm->fv.fv_coef = 0.0;
			Dens0 = fvm->El().fv->m_density.val;
			fvm->fv.innerid = main_innerid;
			//std::cout << fvm->El().fv->m_id << " ";
			for (int nb = 0;nb < fvm->nb_fvs.size();nb++) {
				FV_NB relativ_pos;
				boundary_coefs = 0.0;
				fvm->nb_fvs[nb].fv_coef = 0.0;

				switch (fvm->El().fv->FVdim) {
				case(FVA_L):
					break;
				case(FVA_S):
					Volume = fvm->El().fv->m_s.size();
					for (auto fv_p : fvm->El().fv->m_s.m_lineArray) {
						bool NBIsFound = false;
						if (fvm->nb_fvs[nb].fv->FVdim == FVA_L) {
							//std::cout << "{" << fv_p.m_p1.x << "," << fv_p.m_p1.y << "},{" << fv_p.m_p2.x << "," << fv_p.m_p2.y << "}-{" << nb.fv->m_l.m_p1.x << "," << nb.fv->m_l.m_p1.y << "},{" << nb.fv->m_l.m_p2.x << "," << nb.fv->m_l.m_p2.y << "}" << nb.fv->m_type << "\n";
							if(fv_p == fvm->nb_fvs[nb].fv->m_l) {
								Area = fv_p.size();
								Area_centre = fv_p.getCentre();
								if (Area_centre.x > fvm->El().fv->m_Centre.x) relativ_pos = NB_E;
								if (Area_centre.x < fvm->El().fv->m_Centre.x) relativ_pos = NB_W;
								if (Area_centre.y > fvm->El().fv->m_Centre.y) relativ_pos = NB_N;
								if (Area_centre.y < fvm->El().fv->m_Centre.y) relativ_pos = NB_S;
								NBIsFound = true;
							}
						}
						else if(fvm->nb_fvs[nb].fv->FVdim == FVA_S){
							for (auto fv_nb : fvm->nb_fvs[nb].fv->m_s.m_lineArray) {
								if (fv_p==fv_nb) {
									Area = fv_p.size();
									Area_centre = fv_p.getCentre();
									if (Area_centre.x > fvm->El().fv->m_Centre.x) relativ_pos = NB_E;
									if (Area_centre.x < fvm->El().fv->m_Centre.x) relativ_pos = NB_W;
									if (Area_centre.y > fvm->El().fv->m_Centre.y) relativ_pos = NB_N;
									if (Area_centre.y < fvm->El().fv->m_Centre.y) relativ_pos = NB_S;
									NBIsFound = true;
									break;
								}
							}
						}
						if (!NBIsFound) continue;
						//std::cout << relativ_pos << "\n";
						int nb_id = fvm->nb_fvs[nb].fv->m_id;
						if ((relativ_pos == NB_E) and ((fvm->nb_fvs[nb].fv->m_id + 1) % (last_x_fvs.size()) == 0)) {
							nb_id = FVM.FiniteVolumeMesh[fvm->nb_fvs[nb].fv->m_id]->getPtrToSide(NB_E)->m_id;
						}


						int innerid;				
						if (fvm->nb_fvs[nb].fv->m_id >= m_options.nrOfFVinDir.x* m_options.nrOfFVinDir.y) {
							innerid = fvm->nb_fvs[nb].fv->m_id;
						}
						else {
							if ((fvm->nb_fvs[nb].fv->m_Centre.x < max_x)) innerid = fvm->nb_fvs[nb].fv->m_id - fvm->nb_fvs[nb].fv->m_id / last_x_fvs.size();
							else innerid = FVM.FiniteVolumeMesh[fvm->nb_fvs[nb].fv->m_id - fvm->nb_fvs[nb].fv->m_id / last_x_fvs.size()]->getPtrToSide(NB_E)->m_id;
						}
						fvm->nb_fvs[nb].innerid = innerid;
						//std::cout << innerid << "\n";

						// Pressure
						// dy*(P*P - P*E)
						if (relativ_pos == NB_E) {
							fvm->fv_source += Area * (PApr[fvm->El().fv->m_id] - PApr[fvm->nb_fvs[nb].fv->m_id]);
							//std::cout << fvm->fv_source << "->";
						}

						if (FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
							vec_fv_to_nb = fvm->El().fv->m_Centre - FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_Centre;;
							vec_fv_to_b = fvm->El().fv->m_Centre - Area_centre;
							vec_b_to_nb = Area_centre - FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_Centre;


							if (relativ_pos == NB_E) {}
							if (relativ_pos == NB_W) {}
							if (relativ_pos == NB_N) {}
							if (relativ_pos == NB_S) {}

							//;

						}

						if (FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
							vec_fv_to_nb = fvm->El().fv->m_Centre - FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_Centre;;
							vec_fv_to_b = fvm->El().fv->m_Centre - Area_centre;
							vec_b_to_nb = Area_centre - FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_Centre;

							effD = (FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_DVisc * fvm->El().fv->m_DVisc) / ((FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_DVisc* vec_mod(vec_fv_to_b) / vec_mod(vec_fv_to_nb)) + (fvm->El().fv->m_DVisc *vec_mod(vec_b_to_nb) / vec_mod(vec_fv_to_nb)));
							Dnb_dir = (effD * Area) / vec_mod(vec_fv_to_nb);

							fv_prec_3 vel({0.0,0.0,0.0}); // {UCor_ed[i->m_id],VCor_ed[i->m_id],WCor_ed[i->m_id]}
							vel.x = 0.5*(UCor_ed[fvm->fv.innerid] + UCor_ed[fvm->nb_fvs[nb].innerid]);

							if (relativ_pos == NB_E) vel.y = FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_velocity.dval.y;
							if (relativ_pos == NB_W) vel.y = fvm->El().fv->m_velocity.dval.y;
							//›“Œ —ÀŒ∆ÕŒ
							if (relativ_pos == NB_N) vel.y = 0.25*(fvm->El().fv->m_velocity.dval.y + FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_velocity.dval.y + FVM.FiniteVolumeMesh[fvm->El().fv->m_id]->getPtrToSide(NB_N)->m_velocity.dval.y + FVM.FiniteVolumeMesh[FVM.FiniteVolumeMesh[nb_id]->getPtrToSide(NB_N)->m_id]->FVPtr()->m_velocity.dval.y);
							if (relativ_pos == NB_S) vel.y = 0.25*(fvm->El().fv->m_velocity.dval.y + FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_velocity.dval.y + FVM.FiniteVolumeMesh[fvm->El().fv->m_id]->getPtrToSide(NB_S)->m_velocity.dval.y + FVM.FiniteVolumeMesh[FVM.FiniteVolumeMesh[nb_id]->getPtrToSide(NB_S)->m_id]->FVPtr()->m_velocity.dval.y);

							


							//Fnb_dir = Dens0 * (sqrt(pow(vel.x, 2) + pow(vel.y, 2) + pow(vel.z, 2))) * Area;
							//Fnb_dir = Dens0 * vel.x * Area;
							//// aNB							
							//{
							//	// Diffusion
							//	fvm->nb_fvs[nb].fv_coef += PeFunc(abs(Fnb_dir / Dnb_dir)) * Dnb_dir;
							//	// Convention
							//	if ((relativ_pos == NB_E) and Fnb_dir < 0) fvm->nb_fvs[nb].fv_coef -= Fnb_dir;
							//	if ((relativ_pos == NB_W) and Fnb_dir > 0) fvm->nb_fvs[nb].fv_coef += Fnb_dir;
							//}
							fvm->nb_fvs[nb].fv_coef += Dnb_dir;
						}
						else if (FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_type == FINITE_VOLUME_TYPE::FVT_BOUNDARY) {
							switch (relativ_pos)
							{
							case NB_E:
							case NB_W:
								// œÓÍ‡ ÌÂ ÔÓÌˇÎ Ì‡‰Ó, Ë ÂÒÎË Ì‡‰Ó ÚÓ Í‡Í ÚÓ˜ÌÓ ÔÓÒ˜ËÚ‡¸? ≈ÒÎË ÌÂ„‰Â ÚÓ˜ÌÓ„Ó ‡Ò˜ÂÚ‡ ÌÂÚ, ÚÓ Ì‡‚ÂÌÓÂ ÌÂ Ì‡‰Ó
								break;
							case NB_N:
							case NB_S:
								// ›ÚÓ ‚ÓÁÏÓÊÌÓ ÚÓÊÂ ÌÂ Ì‡‰Ó
								vec_fv_to_nb = fvm->El().fv->m_Centre - FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_Centre;
								vec_fv_to_b = fvm->El().fv->m_Centre - Area_centre;
								effD = fvm->El().fv->m_DVisc;
								boundary_coefs += 1.0 / (1.0 / (effD* vec_mod(vec_fv_to_b) / vec_mod(vec_fv_to_nb)) + 1.0 / fvm->nb_fvs[nb].fv->fv_boundary->getVelocityCoefs().x)*Area;
								//std::cout << "nb.fv->fv_boundary->getVelocityCoefs().x:" << nb.fv->fv_boundary->getVelocityCoefs().x << "\n";
								break;
							default:
								break;
							}
						
						}
						//Ã˚ Ì‡ „‡ÌËˆÂ Ë ÁÌ‡ÂÏ Â∏ Ô‡‡ÏÂÚ˚
					}
					break;
				default:
					break;
				}
				FP0 = UCor0[fvm->fv.innerid]; // UCor_ed[fvm->El().fv->m_id]
				FNB0 = UCor0[fvm->nb_fvs[nb].innerid]; // UCor_ed[fvm->nb_fvs[nb].fv->m_id]
				// sigma*SUM(aNB)
				fvm->fv.fv_coef += fvm->nb_fvs[nb].fv_coef * sigma;  // 
				if (fvm->nb_fvs[nb].fv->m_type == FINITE_VOLUME_TYPE::FVT_BOUNDARY) {
					fvm->fv.fv_coef += boundary_coefs * sigma;
					fvm->fv.fv_coef -= FP0 * boundary_coefs * (1.0 - sigma);
					fvm->fv_source += boundary_coefs * fvm->nb_fvs[nb].fv->fv_boundary->getVelocity().x * sigma; // *(1 - sigma)
					//std::cout << fvm->fv_source << "->";
					fvm->fv_source += boundary_coefs * fvm->nb_fvs[nb].fv->fv_boundary->getVelocity().x * (1.0 - sigma);
					//std::cout << fvm->fv_source << "->";
				}
				// ‘P0*(1-sigma)*SUM(aNB)
				fvm->fv_source -= FP0 * fvm->nb_fvs[nb].fv_coef * (1.0 - sigma); // *(1 - sigma)
				// (1-sigma)*SUM(‘NB0 * aNB)
				fvm->fv_source += FNB0 * fvm->nb_fvs[nb].fv_coef * (1.0 - sigma); // *(1 - sigma)
			}

			// Spm*Vol
			fvm->fv.fv_coef -= fvm->sf.Spm.val * Volume;
			// ‘P0*(1-sigma)*Spm0*Vol
			fvm->fv_source -= FP0 * fvm->sf.Spm.dval * Volume*(1.0 - sigma);
			// (1-sigma)*Spp0*Vol
			fvm->fv_source += fvm->sf.Spp.dval * Volume*(1.0 - sigma);// *(1 - sigma)
			// sigma*Spp*Vol
			fvm->fv_source += fvm->sf.Spp.val * Volume*sigma;// *sigma


			// aP0
			if (dt == 0) { ap0 = 0.0; }
			else { ap0 = Dens0 * Volume / dt; }
			fvm->fv.fv_coef += ap0;
			// ‘P0*aP0
			FP0 = UCor0[fvm->fv.innerid]; // UCor_ed[fvm->El().fv->m_id]
			fvm->fv_source += FP0 * ap0;

			//std::cout <<  "\n";
			main_innerid++;
		}
		std::cout << "PHI_PARAMS::UX\n";
		Solve(&FVMesh_VxA, &UApr.dval);

		// an vn* = SUM(anb vnb*) + bv + dx(PP*-PN*)
		for (int i = 0;i < FVM.FiniteVolumeMesh.size();i++) {
			if (FVM.FiniteVolumeMesh[i]->FVPtr()->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
				float skip_flag = false;
				for (auto y : last_y_fvs) { if (y == i) { skip_flag = true; } }
				if (skip_flag) { continue; }
				FVMesh_VyA.push_back(new FVMeshNode(*(FVM.FiniteVolumeMesh[i]), PHI_PARAMS::UY));
			}
		}
		Dnb_dir = 0.0;
		Fnb_dir = 0.0;

		// ‘P*(aP0 + sigma*SUM(aNB) - Spm*Vol) = SUM(‘NB * aNB) + ‘P0*(aP0 - (1-sigma)*SUM(aNB) - (1-sigma)*Spm0*Vol) + (1-sigma)*SUM(‘NB0 * aNB) + (1-sigma)*Spp0*Vol + sigma*Spp*Vol + dx*(P*P - P*N)
		//      +          +             +               +              +              +                    +                      +                                                            +            
		for (auto fvm : FVMesh_VyA) {
			fv_prec ap0;
			fvm->fv.fv_coef = 0.0;
			Dens0 = fvm->El().fv->m_density.val;
			for (int nb = 0;nb < fvm->nb_fvs.size();nb++) {
				FV_NB relativ_pos;
				boundary_coefs = 0.0;
				fvm->nb_fvs[nb].fv_coef = 0.0;
				switch (fvm->El().fv->FVdim) {
				case(FVA_L):
					break;
				case(FVA_S):
					Volume = fvm->El().fv->m_s.size();
					for (auto fv_p : fvm->El().fv->m_s.m_lineArray) {
						bool NBIsFound = false;
						if (fvm->nb_fvs[nb].fv->FVdim == FVA_L) {
							//std::cout << "{" << fv_p.m_p1.x << "," << fv_p.m_p1.y << "},{" << fv_p.m_p2.x << "," << fv_p.m_p2.y << "}-{" << nb.fv->m_l.m_p1.x << "," << nb.fv->m_l.m_p1.y << "},{" << nb.fv->m_l.m_p2.x << "," << nb.fv->m_l.m_p2.y << "}" << nb.fv->m_type << "\n";
							if (fv_p == fvm->nb_fvs[nb].fv->m_l) {
								Area = fv_p.size();
								Area_centre = fv_p.getCentre();
								if (Area_centre.x > fvm->El().fv->m_Centre.x) relativ_pos = NB_E;
								if (Area_centre.x < fvm->El().fv->m_Centre.x) relativ_pos = NB_W;
								if (Area_centre.y > fvm->El().fv->m_Centre.y) relativ_pos = NB_N;
								if (Area_centre.y < fvm->El().fv->m_Centre.y) relativ_pos = NB_S;
								NBIsFound = true;
							}
						}
						else if (fvm->nb_fvs[nb].fv->FVdim == FVA_S) {
							for (auto fv_nb : fvm->nb_fvs[nb].fv->m_s.m_lineArray) {
								if (fv_p==fv_nb) {
									Area = fv_p.size();
									Area_centre = fv_p.getCentre();
									if (Area_centre.x > fvm->El().fv->m_Centre.x) relativ_pos = NB_E;
									if (Area_centre.x < fvm->El().fv->m_Centre.x) relativ_pos = NB_W;
									if (Area_centre.y > fvm->El().fv->m_Centre.y) relativ_pos = NB_N;
									if (Area_centre.y < fvm->El().fv->m_Centre.y) relativ_pos = NB_S;
									NBIsFound = true;
									break;
								}
							}
						}
						if (!NBIsFound) continue;
						//std::cout << relativ_pos << "\n";
						//std::cout << VCor_ed.size() << "  " << VCor_ed.size() + last_y_fvs.size() << "\n";
						int nb_id = fvm->nb_fvs[nb].fv->m_id;
						if ((relativ_pos == NB_N) and (fvm->El().fv->m_id >= FVMesh_VyA.size())) {
							nb_id = FVM.FiniteVolumeMesh[fvm->nb_fvs[nb].fv->m_id]->getPtrToSide(NB_N)->m_id;
						}
						if ((nb_id >= VCor_ed.size()) and (nb_id < VCor_ed.size() + last_y_fvs.size())) {
							nb_id = FVM.FiniteVolumeMesh[nb_id]->getPtrToSide(NB_N)->m_id;
						}
						//std::cout << " " << fvm->El().fv->m_id << " " << nb_id << " " << relativ_pos  << " "<< "\n"; // —“–¿ÕÕŒ


						int innerid = fvm->nb_fvs[nb].fv->m_id;
						fvm->nb_fvs[nb].innerid = innerid;


						// Pressure
						// dx*(P*P - P*N)
						if (relativ_pos == NB_N) {
							//std::cout << fvm->El().fv->m_id <<" - " << fvm->nb_fvs[nb].fv->m_id << "\n";
							fvm->fv_source += Area * (PApr[fvm->El().fv->m_id] - PApr[fvm->nb_fvs[nb].fv->m_id]);
						}

						if (FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
							vec_fv_to_nb = fvm->El().fv->m_Centre - FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_Centre;
							vec_fv_to_b = fvm->El().fv->m_Centre - Area_centre;
							vec_b_to_nb = Area_centre - FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_Centre;

							effD = (FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_DVisc * fvm->El().fv->m_DVisc) / ((FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_DVisc* vec_mod(vec_fv_to_b) / vec_mod(vec_fv_to_nb)) + (fvm->El().fv->m_DVisc *vec_mod(vec_b_to_nb) / vec_mod(vec_fv_to_nb)));
							Dnb_dir = (effD * Area) / vec_mod(vec_fv_to_nb);
							fv_prec_3 vel(fvm->El().fv->m_velocity.val); // {UCor_ed[i->m_id],VCor_ed[i->m_id],WCor_ed[i->m_id]}
							//std::cout << fvm->fv.innerid << " " << fvm->nb_fvs[nb].innerid << "\n";
							//std::cout << VCor_ed[fvm->fv.innerid] <<", " << VCor_ed[fvm->nb_fvs[nb].innerid] << " ";
							vel.y = 0.5*(VCor_ed[fvm->fv.innerid] + VCor_ed[fvm->nb_fvs[nb].innerid]);
							//std::cout << vel.y << "\n";

							if (relativ_pos == NB_N) vel.x = FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_velocity.dval.x;
							if (relativ_pos == NB_S) vel.x = fvm->El().fv->m_velocity.dval.x;
							//›“Œ —ÀŒ∆ÕŒ
							if (relativ_pos == NB_E) vel.x = 0.25*(fvm->El().fv->m_velocity.dval.x + FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_velocity.dval.x + FVM.FiniteVolumeMesh[fvm->El().fv->m_id]->getPtrToSide(NB_E)->m_velocity.dval.x + FVM.FiniteVolumeMesh[FVM.FiniteVolumeMesh[nb_id]->getPtrToSide(NB_E)->m_id]->FVPtr()->m_velocity.dval.x);
							if (relativ_pos == NB_W) vel.x = 0.25*(fvm->El().fv->m_velocity.dval.x + FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_velocity.dval.x + FVM.FiniteVolumeMesh[fvm->El().fv->m_id]->getPtrToSide(NB_W)->m_velocity.dval.x + FVM.FiniteVolumeMesh[FVM.FiniteVolumeMesh[nb_id]->getPtrToSide(NB_W)->m_id]->FVPtr()->m_velocity.dval.x);

							//Fnb_dir = Dens0 * (sqrt(pow(vel.x, 2) + pow(vel.y, 2) + pow(vel.z, 2))) * Area;
							Fnb_dir = Dens0 * vel.y * Area;
							// aNB							
							//{
							//	// Diffusion
							//	fvm->nb_fvs[nb].fv_coef += PeFunc(abs(Fnb_dir / Dnb_dir)) * Dnb_dir;
							//	// Convention
							//	if ((relativ_pos == NB_N) and Fnb_dir < 0) fvm->nb_fvs[nb].fv_coef -= Fnb_dir;
							//	if ((relativ_pos == NB_S) and Fnb_dir > 0) fvm->nb_fvs[nb].fv_coef += Fnb_dir;
							//}
							fvm->nb_fvs[nb].fv_coef += Dnb_dir;
							//std::cout << "fvm->nb_fvs[nb].fv_coef" << " " << fvm->nb_fvs[nb].fv_coef << "\n"; // —“–¿ÕÕŒ
						}
						else if (FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_type == FINITE_VOLUME_TYPE::FVT_BOUNDARY) {
							switch (relativ_pos)
							{
							case NB_E:
							case NB_W:
								// ›ÚÓ ‚ÓÁÏÓÊÌÓ ÚÓÊÂ ÌÂ Ì‡‰Ó
								vec_fv_to_nb = fvm->El().fv->m_Centre - FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_Centre;
								vec_fv_to_b = fvm->El().fv->m_Centre - Area_centre;
								effD = fvm->El().fv->m_DVisc;
								boundary_coefs += 1.0 / (1.0 / (effD* vec_mod(vec_fv_to_b) / vec_mod(vec_fv_to_nb)) + 1.0 / fvm->nb_fvs[nb].fv->fv_boundary->getVelocityCoefs().y)*Area;
								break;
							case NB_N:
							case NB_S:
								// œÓÍ‡ ÌÂ ÔÓÌˇÎ Ì‡‰Ó, Ë ÂÒÎË Ì‡‰Ó ÚÓ Í‡Í ÚÓ˜ÌÓ ÔÓÒ˜ËÚ‡¸? ≈ÒÎË ÌÂ„‰Â ÚÓ˜ÌÓ„Ó ‡Ò˜ÂÚ‡ ÌÂÚ, ÚÓ Ì‡‚ÂÌÓÂ ÌÂ Ì‡‰Ó
								break;
							default:
								break;
							}
						}
					}
					break;
				default:
					break;
				}
				FP0 = VCor0[fvm->fv.innerid]; // VCor_ed[fvm->El().fv->m_id]
				FNB0 = VCor0[fvm->nb_fvs[nb].innerid]; // VCor_ed[fvm->nb_fvs[nb].fv->m_id]
				// sigma*SUM(aNB)
				fvm->fv.fv_coef += fvm->nb_fvs[nb].fv_coef * sigma;  // 
				if (fvm->nb_fvs[nb].fv->m_type == FINITE_VOLUME_TYPE::FVT_BOUNDARY) {
					fvm->fv.fv_coef += boundary_coefs * sigma;
					fvm->fv.fv_coef -= FP0 * boundary_coefs * (1.0 - sigma);
					fvm->fv_source += boundary_coefs * fvm->nb_fvs[nb].fv->fv_boundary->getVelocity().y * sigma; // *(1 - sigma)
					fvm->fv_source += boundary_coefs * fvm->nb_fvs[nb].fv->fv_boundary->getVelocity().y * (1 - sigma);
				}
				// ‘P0*(1-sigma)*SUM(aNB)
				fvm->fv_source -= FP0 * fvm->nb_fvs[nb].fv_coef * (1 - sigma); // *(1 - sigma)
				// (1-sigma)*SUM(‘NB0 * aNB)
				fvm->fv_source += FNB0 * fvm->nb_fvs[nb].fv_coef * (1 - sigma); // *(1 - sigma)
			}
			// Spm*Vol
			fvm->fv.fv_coef -= fvm->sf.Spm.val * Volume;
			// ‘P0*(1-sigma)*Spm0*Vol
			fvm->fv_source -= FP0 * fvm->sf.Spm.dval * Volume*(1 - sigma);
			// (1-sigma)*Spp0*Vol
			fvm->fv_source += fvm->sf.Spp.dval * Volume*(1 - sigma);// *(1 - sigma)
			// sigma*Spp*Vol
			fvm->fv_source += fvm->sf.Spp.val * Volume*sigma;// *sigma

			// aP0
			if (dt == 0) { ap0 = 0.0; }
			else { ap0 = Dens0 * Volume / dt; }
			fvm->fv.fv_coef += ap0;
			// ‘P0*aP0
			FP0 = VCor0[fvm->fv.innerid]; // VCor_ed[fvm->El().fv->m_id]
			fvm->fv_source += FP0 * ap0;
		}
		std::cout << "PHI_PARAMS::UY\n";
		Solve(&FVMesh_VyA, &VApr.dval);


		std::cout << "UApr.dval: ";
		// Fn = AlfaF*F + (1-AlfaF)*Fn-1
		for (int i = 0;i < UApr.val.size();i++) {
			UApr.dval[i] = AlphaU * UApr.dval[i] + (1 - AlphaU)*UApr.val[i];
			std::cout << UApr.dval[i] << " ";
		}
		std::cout << "\n";
		std::cout << "VApr.dval: ";
		for (int i = 0;i < VApr.val.size();i++) { 
			VApr.dval[i] = AlphaV * VApr.dval[i] + (1 - AlphaV)*VApr.val[i];
			std::cout << VApr.dval[i] << " ";
		}
		std::cout << "\n";


		// ap PP' = SUM(anb PNB') + bp* , bp* = -{(Rhop-Rhop0)/dt dVp + SUM(Fnb*)} 
		for (int i = 0;i < FVM.FiniteVolumeMesh.size();i++) {
			if (FVM.FiniteVolumeMesh[i]->FVPtr()->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
				FVMesh_PC.push_back(new FVMeshNode(*(FVM.FiniteVolumeMesh[i]), PHI_PARAMS::P));
			}
		}
		Dnb_dir = 0.0;
		Fnb_dir = 0.0;

		// ‘P*SUM(dens*dAnb*dnb) = SUM(‘NB*dens*dAnb*dnb) + (-(aP0 + dens * U*e *dAe - dens * U*w *dAw + dens * U*n *dAn - dens * U*s *dAs)
		//      +                           +                   +      +                +                  +                 +                                                        
		for (auto fvm : FVMesh_PC) {
			fv_prec ap0;
			fvm->fv.fv_coef = 0.0;
			Dens0 = fvm->El().fv->m_density.val;
			fv_prec de = 0.0; fv_prec dw = 0.0; fv_prec dn = 0.0; fv_prec ds = 0.0;
			fv_prec ue = 0.0; fv_prec uw = 0.0; fv_prec vn = 0.0; fv_prec vs = 0.0;
			for (int nb = 0;nb < fvm->nb_fvs.size();nb++) {
				FV_NB relativ_pos;
				boundary_coefs = 0.0;
				fv_prec dens = Dens0;
				switch (fvm->El().fv->FVdim) {
				case(FVA_L):
					break;
				case(FVA_S):
					Volume = fvm->El().fv->m_s.size();
					for (auto fv_p : fvm->El().fv->m_s.m_lineArray) {
						bool NBIsFound = false;
						if (fvm->nb_fvs[nb].fv->FVdim == FVA_L) {
							//std::cout << "{" << fv_p.m_p1.x << "," << fv_p.m_p1.y << "},{" << fv_p.m_p2.x << "," << fv_p.m_p2.y << "}-{" << nb.fv->m_l.m_p1.x << "," << nb.fv->m_l.m_p1.y << "},{" << nb.fv->m_l.m_p2.x << "," << nb.fv->m_l.m_p2.y << "}" << nb.fv->m_type << "\n";
							if (fv_p == fvm->nb_fvs[nb].fv->m_l) {
								Area = fv_p.size();
								Area_centre = fv_p.getCentre();
								if (Area_centre.x > fvm->El().fv->m_Centre.x) relativ_pos = NB_E;
								if (Area_centre.x < fvm->El().fv->m_Centre.x) relativ_pos = NB_W;
								if (Area_centre.y > fvm->El().fv->m_Centre.y) relativ_pos = NB_N;
								if (Area_centre.y < fvm->El().fv->m_Centre.y) relativ_pos = NB_S;
								NBIsFound = true;
							}
						}
						else if (fvm->nb_fvs[nb].fv->FVdim == FVA_S) {
							for (auto fv_nb : fvm->nb_fvs[nb].fv->m_s.m_lineArray) {
								if (fv_p == fv_nb) {
									Area = fv_p.size();
									Area_centre = fv_p.getCentre();
									if (Area_centre.x > fvm->El().fv->m_Centre.x) relativ_pos = NB_E;
									if (Area_centre.x < fvm->El().fv->m_Centre.x) relativ_pos = NB_W;
									if (Area_centre.y > fvm->El().fv->m_Centre.y) relativ_pos = NB_N;
									if (Area_centre.y < fvm->El().fv->m_Centre.y) relativ_pos = NB_S;
									NBIsFound = true;
									//std::cout << "{" << fv_p.m_p1.x << "," << fv_p.m_p1.y << "},{" << fv_p.m_p2.x << "," << fv_p.m_p2.y << "}-{" << fv_nb.m_p1.x << "," << fv_nb.m_p1.y << "},{" << fv_nb.m_p2.x << "," << fv_nb.m_p2.y << "}" << "\n";
									//std::cout << fvm->fv.fv->m_id << " " <<  fvm->nb_fvs[nb].fv->m_id << " TUT\n";
									break;
								}
							}
						}
						if (!NBIsFound) continue;
						fvm->nb_fvs[nb].innerid = fvm->nb_fvs[nb].fv->m_id;
						if (relativ_pos == NB_E) {
							int id = fvm->fv.fv->m_id;
							if (fvm->nb_fvs[nb].fv->m_type == FINITE_VOLUME_TYPE::FVT_BOUNDARY) {
								//std::cout << fvm->fv.fv->m_id << "e: " << fvm->nb_fvs[nb].fv->m_id << "\n";
								ue = fvm->nb_fvs[nb].fv->fv_boundary->getVelocity().x;
								de = 0.0;
							}
							else {
								//std::cout << fvm->fv.fv->m_id << "e: " << fvm->getPtrToSide(NB_E)->m_id - 1 - id / last_x_fvs.size() << "\n";
								ue = UApr.dval[fvm->getPtrToSide(NB_E)->m_id - 1 - id / last_x_fvs.size()];
								de = Area / FVMesh_VxA[fvm->getPtrToSide(NB_E)->m_id - 1 - id / last_x_fvs.size()]->fv.fv_coef;
								DensNB = fvm->getPtrToSide(NB_E)->m_density.val;
								//std::cout << FVMesh_VxA[fvm->getPtrToSide(NB_E)->m_id - 1 - id / last_x_fvs.size()]->fv.fv->m_id << " " << FVMesh_VxA[fvm->getPtrToSide(NB_E)->m_id - 1 - id / last_x_fvs.size()]->fv.fv_coef << "\n";
							}
							dens = 0.5*(Dens0 + DensNB);
							fvm->nb_fvs[nb].fv_coef = dens *Area*de;
						}
						else if (relativ_pos == NB_W) {
							if (fvm->fv.fv->m_id % last_x_fvs.size() == 0) {
								//std::cout << fvm->fv.fv->m_id << "w: " << FVMesh_PC[fvm->fv.fv->m_id]->getPtrToSide(NB_W)->m_id << "\n";
								uw = FVMesh_PC[fvm->fv.fv->m_id]->getPtrToSide(NB_W)->fv_boundary->getVelocity().x;
								dw = 0.0;
							}
							else {
								int id = fvm->getPtrToSide(NB_W)->m_id;
								//std::cout << fvm->fv.fv->m_id << "w: " << id - static_cast<int>(id / (last_x_fvs.size())) << "\n";
								uw = UApr.dval[id- static_cast<int>(id / (last_x_fvs.size()))];
								dw = Area / FVMesh_VxA[id- static_cast<int>(id / (last_x_fvs.size()))]->fv.fv_coef;
								DensNB = fvm->getPtrToSide(NB_W)->m_density.val;
								//std::cout << FVMesh_VxA[id - static_cast<int>(id / (last_x_fvs.size()))]->fv.fv->m_id << " " << FVMesh_VxA[id - static_cast<int>(id / (last_x_fvs.size()))]->fv.fv_coef << "\n";
							}
							dens = 0.5*(Dens0 + DensNB);
							fvm->nb_fvs[nb].fv_coef = dens *Area*dw;
						}
						else if (relativ_pos == NB_N) {
							if (fvm->fv.fv->m_id >= VApr.dval.size()) {
								//std::cout << fvm->fv.fv->m_id << "n: " << FVMesh_PC[fvm->fv.fv->m_id]->getPtrToSide(NB_N)->m_id << "\n";
								vn = FVMesh_PC[fvm->fv.fv->m_id]->getPtrToSide(NB_N)->fv_boundary->getVelocity().y;
								dn = 0.0;
							}
							else {
								//std::cout << fvm->fv.fv->m_id << "n: " << fvm->fv.fv->m_id << "\n";
								vn = VApr.dval[fvm->fv.fv->m_id];
								dn = Area / FVMesh_VxA[fvm->fv.fv->m_id]->fv.fv_coef;
								DensNB = fvm->getPtrToSide(NB_N)->m_density.val;
								//std::cout << FVMesh_VxA[fvm->fv.fv->m_id]->fv.fv->m_id << " " << FVMesh_VxA[fvm->fv.fv->m_id]->fv.fv_coef << "\n";
							}
							dens = 0.5*(Dens0 + DensNB);
							fvm->nb_fvs[nb].fv_coef = dens *Area*dn;
						}
						else if (relativ_pos == NB_S) {
							int id = fvm->fv.fv->m_id;
							if (id < last_y_fvs.size()) {
								//std::cout << fvm->fv.fv->m_id << "s: " << FVMesh_PC[id]->getPtrToSide(NB_S)->m_id << "\n";
								vs = FVMesh_PC[id]->getPtrToSide(NB_S)->fv_boundary->getVelocity().y;
								ds = 0.0;
							}
							else {
								//std::cout << fvm->fv.fv->m_id << "s: " << id - last_y_fvs.size() << "\n";
								vs = VApr.dval[id- last_y_fvs.size()];
								ds = Area / FVMesh_VxA[id- last_y_fvs.size()]->fv.fv_coef;
								DensNB = fvm->getPtrToSide(NB_S)->m_density.val;
								//std::cout << FVMesh_VxA[id - last_y_fvs.size()]->fv.fv->m_id << " " << FVMesh_VxA[id - last_y_fvs.size()]->fv.fv_coef << "\n";
							}
							dens = 0.5*(Dens0 + DensNB);
							fvm->nb_fvs[nb].fv_coef = dens * Area*ds;
						}
					}
					break;
				default:
					break;
				}
				//if (relativ_pos == NB_E) { std::cout << "de=" << de << " - " << "ue=" << ue << "\n"; }
				//else if (relativ_pos == NB_W) { std::cout << "dw=" << dw << " - " << "uw=" << uw << "\n"; }
				//else if (relativ_pos == NB_N) { std::cout << "dn=" << dn << " - " << "vn=" << vn << "\n"; }
				//else if (relativ_pos == NB_S) { std::cout << "ds=" << ds << " - " << "vs=" << vs << "\n"; }

				fvm->fv.fv_coef += fvm->nb_fvs[nb].fv_coef;
				fvm->fv_source += dens *Area * (uw - ue + vs - vn);

				de = 0.0;  dw = 0.0;  dn = 0.0;  ds = 0.0;
				ue = 0.0;  uw = 0.0;  vn = 0.0;  vs = 0.0;
			}

			if (dt == 0) { ap0 = 0.0; }
			else { ap0 = (Dens - Dens0) * Volume / dt; }
			fvm->fv_source -= ap0;
		}
		std::cout << "PHI_PARAMS::P\n";


		sum.val = sum.dval;
		sum.dval = 0.0;
		for (auto err : FVMesh_PC) {
			sum.dval += err->fv_source;
		}
		std::cout << "sum.val: " << sum.val << ", sum.dval: " << sum.dval << "\n";

		SolveWithFirstZero(&FVMesh_PC, &Pcor);

		// ue = ue* +de(PP'- PE')
		std::cout << "UCor_ed: ";
		for (int i = 0;i < UApr.dval.size();i++) {
			//int num = i + i / (last_x_fvs.size() - 1);
			//FVM.FiniteVolumes[num]->m_velocity.dval.x = UApr[i] + FVMesh_VxA[i]->getCoef(FVMesh_VxA[i]->fv.fv) * (Pcor[i] - Pcor[FVMesh_VyA[i]->getPtrToSide(NB_E)->m_id]);
			UCor_ed[i] = UApr.dval[i] + FVMesh_VxA[i]->getCoef(FVMesh_VxA[i]->fv.fv) * (Pcor[i] - Pcor[FVMesh_VyA[i]->getPtrToSide(NB_E)->m_id]);
			std::cout << UCor_ed[i] << " ";
			UApr.val[i] = UApr.dval[i];
		}
		std::cout << "\n";
		// vn = vn* +dn(PP'- PN')
		std::cout << "VCor_ed: ";
		for (int i = 0;i < VApr.dval.size();i++) {
			//FVM.FiniteVolumes[i]->m_velocity.dval.y = VApr[i] + FVMesh_VyA[i]->getCoef(FVMesh_VyA[i]->fv.fv) * (Pcor[i] - Pcor[FVMesh_VyA[i]->getPtrToSide(NB_N)->m_id]);
			VCor_ed[i] = VApr.dval[i] + FVMesh_VyA[i]->getCoef(FVMesh_VyA[i]->fv.fv) * (Pcor[i] - Pcor[FVMesh_VyA[i]->getPtrToSide(NB_N)->m_id]);
			std::cout << VCor_ed[i] << " ";
			VApr.val[i] = VApr.dval[i];
		}
		std::cout << "\n";


		//ÒÍÓÓÒÚË
		std::cout << "\n";
		std::cout << "Velocity.x: ";
		for (auto fvm : FVM.FiniteVolumeMesh) {
			if(fvm->FVPtr()->m_type ==  FINITE_VOLUME_TYPE::FVT_DEFAULT){
				fv_prec velE = 0.0;
				fv_prec velW = 0.0;
				//std::cout << (fvm->FVPtr()->m_id) / (last_x_fvs.size()) << " ";
				//std::cout << fvm->FVPtr()->m_id << ": ";
				if (fvm->getPtrToSide(NB_W)->fv_boundary != nullptr) {
					velW = fvm->getPtrToSide(NB_W)->fv_boundary->getVelocity().x;
					//std::cout << fvm->getPtrToSide(NB_W)->m_id << " ";
				}
				else {
					velW = UCor_ed[fvm->getPtrToSide(NB_W)->m_id - (fvm->FVPtr()->m_id) / (last_x_fvs.size())];
					//std::cout << fvm->getPtrToSide(NB_W)->m_id- (fvm->FVPtr()->m_id) / (last_x_fvs.size()) << " ";
				}
				if (fvm->getPtrToSide(NB_E)->fv_boundary != nullptr) { 
					velE = fvm->getPtrToSide(NB_E)->fv_boundary->getVelocity().x;
					//std::cout << fvm->getPtrToSide(NB_E)->m_id << " ";
				}
				else { 
					velE = UCor_ed[fvm->FVPtr()->m_id - (fvm->FVPtr()->m_id) / (last_x_fvs.size())];
					//std::cout << fvm->FVPtr()->m_id - (fvm->FVPtr()->m_id) / (last_x_fvs.size()) << " ";
				}
				fvm->FVPtr()->m_velocity.dval.x = 0.5*(velW + velE);
				//std::cout << "\n";
				std::cout << fvm->FVPtr()->m_velocity.dval.x << " ";
			}
		}
		std::cout << "\n";
		std::cout << "Velocity.y: ";
		for (auto fvm : FVM.FiniteVolumeMesh) {
			if (fvm->FVPtr()->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
				fv_prec velN = 0.0;
				fv_prec velS = 0.0;
				//std::cout << (fvm->FVPtr()->m_id) / (last_x_fvs.size()) << " ";
				//std::cout << fvm->FVPtr()->m_id << ": ";
				if (fvm->getPtrToSide(NB_S)->fv_boundary != nullptr) {
					velS = fvm->getPtrToSide(NB_S)->fv_boundary->getVelocity().y;
					//std::cout << fvm->getPtrToSide(NB_S)->m_id << " ";
				}
				else {
					velS = VCor_ed[fvm->getPtrToSide(NB_S)->m_id];
					//std::cout << fvm->getPtrToSide(NB_S)->m_id << " ";
				}
				if (fvm->getPtrToSide(NB_N)->fv_boundary != nullptr) { velN = fvm->getPtrToSide(NB_N)->fv_boundary->getVelocity().y;
					//std::cout << fvm->getPtrToSide(NB_N)->m_id << " ";
				}
				else { velN = VCor_ed[fvm->FVPtr()->m_id];
					//std::cout << fvm->FVPtr()->m_id << " ";
				}
				fvm->FVPtr()->m_velocity.dval.y = 0.5*(velS + velN);
				std::cout << fvm->FVPtr()->m_velocity.dval.y << " ";
				//std::cout << "\n";
				//std::cout << velS << " " << velN << "\n";
			}
		}
		std::cout << "\n";


		// P = P* + P'

		std::cout << "Pressure: ";
		for (auto i : FVM.FiniteVolumes) {
			if (i->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
				//i->m_pressure.val = i->m_pressure.dval;
				i->m_pressure.dval = PApr[i->m_id] + AlphaP * Pcor[i->m_id];
				std::cout << i->m_pressure.dval << " ";
			}
		}
		std::cout << "\n";
		std::cout << "_________________________________________________ \n";

		// P* = P
		for (auto i : FVM.FiniteVolumes) {
			if (i->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
				PApr[i->m_id] = i->m_pressure.dval;
			}
		}
		for (auto i : FVMesh_VxA) { delete i; }
		FVMesh_VxA.resize(0);
		for (auto i : FVMesh_VyA) { delete i; }
		FVMesh_VyA.resize(0);
		for (auto i : FVMesh_PC) { delete i; }
		FVMesh_PC.resize(0);

		std::cout << "Err: " << abs(sum.dval - sum.val) << "\n";

	}




	// U0 = U; ÕŒ Õ¿ƒŒ œ≈–≈—“–Œ»“‹ «Õ¿◊≈Õ»ﬂ Œ“ —Ã≈Ÿ≈ÕÕ€’  Œ   ÕŒ–Ã¿À‹Õ€Ã ƒÀﬂ — Œ–Œ—“»
	// P0 = P
	for (auto i : FVM.FiniteVolumes) {
		if (i->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
			i->m_pressure.val = i->m_pressure.dval;
			i->m_velocity.val = i->m_velocity.dval;
		}
	}
	*/

	std::cout << "\n";
}

/*
void FVM_CD::SIMPLE_ALGO(cd_prec dt) {
	   

	std::vector<fv_prec> PApr;
	PApr.resize(FVM.FiniteVolumes.size(), 0.0);
	std::vector<fv_prec> Pcor;
	Pcor.resize(FVM.FiniteVolumes.size(), 0.0);


	std::vector<Ch_<fv_prec>> Eps;
	fv_prec sigma = 1.0;



	fv_prec AlphaP = 0.8;
	fv_prec AlphaU = 0.5;
	fv_prec AlphaV = 0.5;
	fv_prec AlphaW = 0.5;

	std::vector<FVMeshNode*> FVMesh_VxA;
	std::vector<FVMeshNode*> FVMesh_VyA;
	std::vector<FVMeshNode*> FVMesh_PC;

	Ch_<fv_prec> sum;
	sum.dval = 10;
	sum.val = 0;
	for (auto e : Eps) { sum = sum + e; }

	for (auto i : FVM.FiniteVolumes) {
		//i->m_pressure.dval = i->m_pressure.val;
		PApr[i->m_id] = i->m_pressure.val;
		i->m_pressure.dval = i->m_pressure.val;
		i->m_velocity.dval = i->m_velocity.val;
		//  ÕŒ Õ¿ƒŒ œ≈–≈—“–Œ»“‹ «Õ¿◊≈Õ»ﬂ Œ“ ÕŒ–Ã¿À‹Õ€’  Œ   —Ã≈Ÿ≈ÕÕ€Ã ƒÀﬂ — Œ–Œ—“», » “Œ√ƒ¿ ”Ã≈Õ‹ÿ»“‹ –¿«Ã≈–€ UCor_ed, VCor_ed, WCor_ed
		//  Õ” “Œ√ƒ¿ œŒÃ≈Õﬂ“‹ » –¿«Ã≈–€ UApr, VApr, WApr
	}


	std::vector<int> last_x_fvs;
	std::vector<int> last_y_fvs;
	fv_prec max_x = 0.0; fv_prec max_y = 0.0; fv_prec max_z = 0.0;
	for (int i = 0;i < FVM.FiniteVolumeMesh.size();i++) {
		if (max_x < FVM.FiniteVolumeMesh[i]->FVPtr()->m_Centre.x) max_x = FVM.FiniteVolumeMesh[i]->FVPtr()->m_Centre.x;
		if (max_y < FVM.FiniteVolumeMesh[i]->FVPtr()->m_Centre.y) max_y = FVM.FiniteVolumeMesh[i]->FVPtr()->m_Centre.y;
		if (max_z < FVM.FiniteVolumeMesh[i]->FVPtr()->m_Centre.z) max_z = FVM.FiniteVolumeMesh[i]->FVPtr()->m_Centre.z;
	}
	max_x = 1.00 - m_options.average_dim_steps.x;
	max_y = 1.00 - m_options.average_dim_steps.y;
	for (int i = 0;i < FVM.FiniteVolumeMesh.size();i++) {
		if (FVM.FiniteVolumeMesh[i]->FVPtr()->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
			if (FVM.FiniteVolumeMesh[i]->FVPtr()->m_Centre.x >= max_x) { last_x_fvs.push_back(i); }
			if (FVM.FiniteVolumeMesh[i]->FVPtr()->m_Centre.y >= max_y) { last_y_fvs.push_back(i); }
		}
	}

	int num_of_def = 0;
	for (int i = 0;i < FVM.FiniteVolumeMesh.size();i++) {
		if (FVM.FiniteVolumeMesh[i]->FVPtr()->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
			num_of_def++;
		}
	}


	Ch_<std::vector<fv_prec>> UApr;
	UApr.val.resize(num_of_def - last_x_fvs.size(), 0.0);
	UApr.dval.resize(num_of_def - last_x_fvs.size(), 0.0);
	std::vector<fv_prec> UCor_ed;
	UCor_ed.resize(num_of_def - last_x_fvs.size(), 0.0);
	std::vector<fv_prec> UCor0;
	UCor0.resize(num_of_def - last_x_fvs.size(), 0.0);

	Ch_<std::vector<fv_prec>> VApr;
	VApr.val.resize(num_of_def - last_y_fvs.size(), 0.0);
	VApr.dval.resize(num_of_def - last_y_fvs.size(), 0.0);
	std::vector<fv_prec> VCor_ed;
	VCor_ed.resize(num_of_def - last_y_fvs.size(), 0.0);
	std::vector<fv_prec> VCor0;
	VCor0.resize(num_of_def - last_y_fvs.size(), 0.0);


	std::vector<fv_prec> WApr;
	WApr.resize(FVM.FiniteVolumes.size(), 0.0);
	std::vector<fv_prec> WCor_ed;
	WCor_ed.resize(FVM.FiniteVolumes.size(), 0.0);


	std::cout << "UCor0[n]: ";
	int n = 0;
	for (auto fvm : FVM.FiniteVolumeMesh) {
		if (((fvm->FVPtr()->m_id + 1) % last_y_fvs.size() != 0) and (n < num_of_def - last_x_fvs.size())) {
			UCor0[n] = 0.5*(fvm->FVPtr()->m_velocity.val.x + fvm->getPtrToSide(NB_E)->m_velocity.val.x);
			UCor_ed[n] = 0.5*(fvm->FVPtr()->m_velocity.dval.x + fvm->getPtrToSide(NB_E)->m_velocity.dval.x);
			std::cout << UCor0[n] << " ";
			n++;
		}
	}
	std::cout << "\n";


	std::cout << "VCor0[n]: ";
	n = 0;
	for (auto fvm : FVM.FiniteVolumeMesh) {
		if (n < num_of_def - last_y_fvs.size()) {
			VCor0[n] = 0.5*(fvm->FVPtr()->m_velocity.val.y + fvm->getPtrToSide(NB_N)->m_velocity.val.y);
			VCor_ed[n] = 0.5*(fvm->FVPtr()->m_velocity.dval.y + fvm->getPtrToSide(NB_N)->m_velocity.dval.y);
			std::cout << VCor0[n] << " ";
			n++;
		}
	}
	std::cout << "\n";

	//for (int i = 0; i < VCor0.size();i++) {	std::cout << "{" << UCor_ed[i] << ", " << VCor_ed[i] << "}\n"; }

	//n = 0;
	//for (auto fvm : FVM.FiniteVolumeMesh) {
	//	if (n < num_of_def - last_x_fvs.size()) {
	//		VCor_ed[n] = fvm->FVPtr()->m_velocity.val.x + fvm->getPtrToSide(NB_N)->m_velocity.val.x;
	//		n++;
	//	}
	//}

	std::cout << "Err: " << abs(sum.dval - sum.val) << "\n";
	while (abs(sum.dval - sum.val) > 0.00001) {


		// ae ue* = SUM(anb unb*) + bu + dy(PP*-PE*)
		for (int i = 0;i < FVM.FiniteVolumeMesh.size();i++) {
			if (FVM.FiniteVolumeMesh[i]->FVPtr()->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
				float skip_flag = false;
				for (auto x : last_x_fvs) { if (x == i) { skip_flag = true; } }
				if (skip_flag) { continue; }
				FVMesh_VxA.push_back(new FVMeshNode(*(FVM.FiniteVolumeMesh[i]), PHI_PARAMS::UX));
			}
		}
		fv_prec effD;
		fv_prec_3 vec_fv_to_nb;
		fv_prec_3 vec_fv_to_b;
		fv_prec_3 vec_b_to_nb;
		fv_prec Area;
		fv_prec Volume;
		fv_prec_3 Area_centre;
		fv_prec boundary_coefs;
		fv_prec Dnb_dir = 0.0;
		fv_prec Fnb_dir = 0.0;
		fv_prec FP0;
		fv_prec FNB0;
		fv_prec Dens = 1000.0;
		fv_prec Dens0;
		fv_prec DensNB;

		// ‘P*(aP0 + sigma*SUM(aNB) - Spm*Vol) = SUM(‘NB * aNB) + ‘P0*(aP0 - (1-sigma)*SUM(aNB) - (1-sigma)*Spm0*Vol) + (1-sigma)*SUM(‘NB0 * aNB) + (1-sigma)*Spp0*Vol + sigma*Spp*Vol + dy*(P*P - P*E)
		//      +          +             +               +              +              +                    +                      +                                                            +            
		glm::max(0.0, 0.1);

		int main_innerid = 0;
		for (auto fvm : FVMesh_VxA) {
			fv_prec ap0;
			fvm->fv.fv_coef = 0.0;
			Dens0 = fvm->El().fv->m_density.val;
			fvm->fv.innerid = main_innerid;
			//std::cout << fvm->El().fv->m_id << " ";
			for (int nb = 0;nb < fvm->nb_fvs.size();nb++) {
				FV_NB relativ_pos;
				boundary_coefs = 0.0;
				fvm->nb_fvs[nb].fv_coef = 0.0;

				switch (fvm->El().fv->FVdim) {
				case(FVA_L):
					break;
				case(FVA_S):
					Volume = fvm->El().fv->m_s.size();
					for (auto fv_p : fvm->El().fv->m_s.m_lineArray) {
						bool NBIsFound = false;
						if (fvm->nb_fvs[nb].fv->FVdim == FVA_L) {
							//std::cout << "{" << fv_p.m_p1.x << "," << fv_p.m_p1.y << "},{" << fv_p.m_p2.x << "," << fv_p.m_p2.y << "}-{" << nb.fv->m_l.m_p1.x << "," << nb.fv->m_l.m_p1.y << "},{" << nb.fv->m_l.m_p2.x << "," << nb.fv->m_l.m_p2.y << "}" << nb.fv->m_type << "\n";
							if (fv_p == fvm->nb_fvs[nb].fv->m_l) {
								Area = fv_p.size();
								Area_centre = fv_p.getCentre();
								if (Area_centre.x > fvm->El().fv->m_Centre.x) relativ_pos = NB_E;
								if (Area_centre.x < fvm->El().fv->m_Centre.x) relativ_pos = NB_W;
								if (Area_centre.y > fvm->El().fv->m_Centre.y) relativ_pos = NB_N;
								if (Area_centre.y < fvm->El().fv->m_Centre.y) relativ_pos = NB_S;
								NBIsFound = true;
							}
						}
						else if (fvm->nb_fvs[nb].fv->FVdim == FVA_S) {
							for (auto fv_nb : fvm->nb_fvs[nb].fv->m_s.m_lineArray) {
								if (fv_p == fv_nb) {
									Area = fv_p.size();
									Area_centre = fv_p.getCentre();
									if (Area_centre.x > fvm->El().fv->m_Centre.x) relativ_pos = NB_E;
									if (Area_centre.x < fvm->El().fv->m_Centre.x) relativ_pos = NB_W;
									if (Area_centre.y > fvm->El().fv->m_Centre.y) relativ_pos = NB_N;
									if (Area_centre.y < fvm->El().fv->m_Centre.y) relativ_pos = NB_S;
									NBIsFound = true;
									break;
								}
							}
						}
						if (!NBIsFound) continue;
						//std::cout << relativ_pos << "\n";
						int nb_id = fvm->nb_fvs[nb].fv->m_id;
						if ((relativ_pos == NB_E) and ((fvm->nb_fvs[nb].fv->m_id + 1) % (last_x_fvs.size()) == 0)) {
							nb_id = FVM.FiniteVolumeMesh[fvm->nb_fvs[nb].fv->m_id]->getPtrToSide(NB_E)->m_id;
						}


						int innerid;
						if (fvm->nb_fvs[nb].fv->m_id >= m_options.nrOfFVinDir.x* m_options.nrOfFVinDir.y) {
							innerid = fvm->nb_fvs[nb].fv->m_id;
						}
						else {
							if ((fvm->nb_fvs[nb].fv->m_Centre.x < max_x)) innerid = fvm->nb_fvs[nb].fv->m_id - fvm->nb_fvs[nb].fv->m_id / last_x_fvs.size();
							else innerid = FVM.FiniteVolumeMesh[fvm->nb_fvs[nb].fv->m_id - fvm->nb_fvs[nb].fv->m_id / last_x_fvs.size()]->getPtrToSide(NB_E)->m_id;
						}
						fvm->nb_fvs[nb].innerid = innerid;
						//std::cout << innerid << "\n";

						// Pressure
						// dy*(P*P - P*E)
						if (relativ_pos == NB_E) {
							fvm->fv_source += Area * (PApr[fvm->El().fv->m_id] - PApr[fvm->nb_fvs[nb].fv->m_id]);
							//std::cout << fvm->fv_source << "->";
						}

						if (FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
							vec_fv_to_nb = fvm->El().fv->m_Centre - FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_Centre;;
							vec_fv_to_b = fvm->El().fv->m_Centre - Area_centre;
							vec_b_to_nb = Area_centre - FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_Centre;


							if (relativ_pos == NB_E) {}
							if (relativ_pos == NB_W) {}
							if (relativ_pos == NB_N) {}
							if (relativ_pos == NB_S) {}

							//;

						}

						if (FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
							vec_fv_to_nb = fvm->El().fv->m_Centre - FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_Centre;;
							vec_fv_to_b = fvm->El().fv->m_Centre - Area_centre;
							vec_b_to_nb = Area_centre - FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_Centre;

							effD = (FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_DVisc * fvm->El().fv->m_DVisc) / ((FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_DVisc* vec_mod(vec_fv_to_b) / vec_mod(vec_fv_to_nb)) + (fvm->El().fv->m_DVisc *vec_mod(vec_b_to_nb) / vec_mod(vec_fv_to_nb)));
							Dnb_dir = (effD * Area) / vec_mod(vec_fv_to_nb);

							fv_prec_3 vel({ 0.0,0.0,0.0 }); // {UCor_ed[i->m_id],VCor_ed[i->m_id],WCor_ed[i->m_id]}
							vel.x = 0.5*(UCor_ed[fvm->fv.innerid] + UCor_ed[fvm->nb_fvs[nb].innerid]);

							if (relativ_pos == NB_E) vel.y = FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_velocity.dval.y;
							if (relativ_pos == NB_W) vel.y = fvm->El().fv->m_velocity.dval.y;
							//›“Œ —ÀŒ∆ÕŒ
							if (relativ_pos == NB_N) vel.y = 0.25*(fvm->El().fv->m_velocity.dval.y + FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_velocity.dval.y + FVM.FiniteVolumeMesh[fvm->El().fv->m_id]->getPtrToSide(NB_N)->m_velocity.dval.y + FVM.FiniteVolumeMesh[FVM.FiniteVolumeMesh[nb_id]->getPtrToSide(NB_N)->m_id]->FVPtr()->m_velocity.dval.y);
							if (relativ_pos == NB_S) vel.y = 0.25*(fvm->El().fv->m_velocity.dval.y + FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_velocity.dval.y + FVM.FiniteVolumeMesh[fvm->El().fv->m_id]->getPtrToSide(NB_S)->m_velocity.dval.y + FVM.FiniteVolumeMesh[FVM.FiniteVolumeMesh[nb_id]->getPtrToSide(NB_S)->m_id]->FVPtr()->m_velocity.dval.y);




							//Fnb_dir = Dens0 * (sqrt(pow(vel.x, 2) + pow(vel.y, 2) + pow(vel.z, 2))) * Area;
							//Fnb_dir = Dens0 * vel.x * Area;
							//// aNB							
							//{
							//	// Diffusion
							//	fvm->nb_fvs[nb].fv_coef += PeFunc(abs(Fnb_dir / Dnb_dir)) * Dnb_dir;
							//	// Convention
							//	if ((relativ_pos == NB_E) and Fnb_dir < 0) fvm->nb_fvs[nb].fv_coef -= Fnb_dir;
							//	if ((relativ_pos == NB_W) and Fnb_dir > 0) fvm->nb_fvs[nb].fv_coef += Fnb_dir;
							//}
							fvm->nb_fvs[nb].fv_coef += Dnb_dir;
						}
						else if (FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_type == FINITE_VOLUME_TYPE::FVT_BOUNDARY) {
							switch (relativ_pos)
							{
							case NB_E:
							case NB_W:
								// œÓÍ‡ ÌÂ ÔÓÌˇÎ Ì‡‰Ó, Ë ÂÒÎË Ì‡‰Ó ÚÓ Í‡Í ÚÓ˜ÌÓ ÔÓÒ˜ËÚ‡¸? ≈ÒÎË ÌÂ„‰Â ÚÓ˜ÌÓ„Ó ‡Ò˜ÂÚ‡ ÌÂÚ, ÚÓ Ì‡‚ÂÌÓÂ ÌÂ Ì‡‰Ó
								break;
							case NB_N:
							case NB_S:
								// ›ÚÓ ‚ÓÁÏÓÊÌÓ ÚÓÊÂ ÌÂ Ì‡‰Ó
								vec_fv_to_nb = fvm->El().fv->m_Centre - FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_Centre;
								vec_fv_to_b = fvm->El().fv->m_Centre - Area_centre;
								effD = fvm->El().fv->m_DVisc;
								boundary_coefs += 1.0 / (1.0 / (effD* vec_mod(vec_fv_to_b) / vec_mod(vec_fv_to_nb)) + 1.0 / fvm->nb_fvs[nb].fv->fv_boundary->getVelocityCoefs().x)*Area;
								//std::cout << "nb.fv->fv_boundary->getVelocityCoefs().x:" << nb.fv->fv_boundary->getVelocityCoefs().x << "\n";
								break;
							default:
								break;
							}

						}
					}
					break;
				default:
					break;
				}
				FP0 = UCor0[fvm->fv.innerid]; // UCor_ed[fvm->El().fv->m_id]
				FNB0 = UCor0[fvm->nb_fvs[nb].innerid]; // UCor_ed[fvm->nb_fvs[nb].fv->m_id]
				// sigma*SUM(aNB)
				fvm->fv.fv_coef += fvm->nb_fvs[nb].fv_coef * sigma;  // 
				if (fvm->nb_fvs[nb].fv->m_type == FINITE_VOLUME_TYPE::FVT_BOUNDARY) {
					fvm->fv.fv_coef += boundary_coefs * sigma;
					fvm->fv.fv_coef -= FP0 * boundary_coefs * (1.0 - sigma);
					fvm->fv_source += boundary_coefs * fvm->nb_fvs[nb].fv->fv_boundary->getVelocity().x * sigma; // *(1 - sigma)
					//std::cout << fvm->fv_source << "->";
					fvm->fv_source += boundary_coefs * fvm->nb_fvs[nb].fv->fv_boundary->getVelocity().x * (1.0 - sigma);
					//std::cout << fvm->fv_source << "->";
				}
				// ‘P0*(1-sigma)*SUM(aNB)
				fvm->fv_source -= FP0 * fvm->nb_fvs[nb].fv_coef * (1.0 - sigma); // *(1 - sigma)
				// (1-sigma)*SUM(‘NB0 * aNB)
				fvm->fv_source += FNB0 * fvm->nb_fvs[nb].fv_coef * (1.0 - sigma); // *(1 - sigma)
			}

			// Spm*Vol
			fvm->fv.fv_coef -= fvm->sf.Spm.val * Volume;
			// ‘P0*(1-sigma)*Spm0*Vol
			fvm->fv_source -= FP0 * fvm->sf.Spm.dval * Volume*(1.0 - sigma);
			// (1-sigma)*Spp0*Vol
			fvm->fv_source += fvm->sf.Spp.dval * Volume*(1.0 - sigma);// *(1 - sigma)
			// sigma*Spp*Vol
			fvm->fv_source += fvm->sf.Spp.val * Volume*sigma;// *sigma


			// aP0
			if (dt == 0) { ap0 = 0.0; }
			else { ap0 = Dens0 * Volume / dt; }
			fvm->fv.fv_coef += ap0;
			// ‘P0*aP0
			FP0 = UCor0[fvm->fv.innerid]; // UCor_ed[fvm->El().fv->m_id]
			fvm->fv_source += FP0 * ap0;

			//std::cout <<  "\n";
			main_innerid++;
		}
		std::cout << "PHI_PARAMS::UX\n";
		Solve(&FVMesh_VxA, &UApr.dval);

		// an vn* = SUM(anb vnb*) + bv + dx(PP*-PN*)
		for (int i = 0;i < FVM.FiniteVolumeMesh.size();i++) {
			if (FVM.FiniteVolumeMesh[i]->FVPtr()->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
				float skip_flag = false;
				for (auto y : last_y_fvs) { if (y == i) { skip_flag = true; } }
				if (skip_flag) { continue; }
				FVMesh_VyA.push_back(new FVMeshNode(*(FVM.FiniteVolumeMesh[i]), PHI_PARAMS::UY));
			}
		}
		Dnb_dir = 0.0;
		Fnb_dir = 0.0;

		// ‘P*(aP0 + sigma*SUM(aNB) - Spm*Vol) = SUM(‘NB * aNB) + ‘P0*(aP0 - (1-sigma)*SUM(aNB) - (1-sigma)*Spm0*Vol) + (1-sigma)*SUM(‘NB0 * aNB) + (1-sigma)*Spp0*Vol + sigma*Spp*Vol + dx*(P*P - P*N)
		//      +          +             +               +              +              +                    +                      +                                                            +            
		for (auto fvm : FVMesh_VyA) {
			fv_prec ap0;
			fvm->fv.fv_coef = 0.0;
			Dens0 = fvm->El().fv->m_density.val;
			for (int nb = 0;nb < fvm->nb_fvs.size();nb++) {
				FV_NB relativ_pos;
				boundary_coefs = 0.0;
				fvm->nb_fvs[nb].fv_coef = 0.0;
				switch (fvm->El().fv->FVdim) {
				case(FVA_L):
					break;
				case(FVA_S):
					Volume = fvm->El().fv->m_s.size();
					for (auto fv_p : fvm->El().fv->m_s.m_lineArray) {
						bool NBIsFound = false;
						if (fvm->nb_fvs[nb].fv->FVdim == FVA_L) {
							//std::cout << "{" << fv_p.m_p1.x << "," << fv_p.m_p1.y << "},{" << fv_p.m_p2.x << "," << fv_p.m_p2.y << "}-{" << nb.fv->m_l.m_p1.x << "," << nb.fv->m_l.m_p1.y << "},{" << nb.fv->m_l.m_p2.x << "," << nb.fv->m_l.m_p2.y << "}" << nb.fv->m_type << "\n";
							if (fv_p == fvm->nb_fvs[nb].fv->m_l) {
								Area = fv_p.size();
								Area_centre = fv_p.getCentre();
								if (Area_centre.x > fvm->El().fv->m_Centre.x) relativ_pos = NB_E;
								if (Area_centre.x < fvm->El().fv->m_Centre.x) relativ_pos = NB_W;
								if (Area_centre.y > fvm->El().fv->m_Centre.y) relativ_pos = NB_N;
								if (Area_centre.y < fvm->El().fv->m_Centre.y) relativ_pos = NB_S;
								NBIsFound = true;
							}
						}
						else if (fvm->nb_fvs[nb].fv->FVdim == FVA_S) {
							for (auto fv_nb : fvm->nb_fvs[nb].fv->m_s.m_lineArray) {
								if (fv_p == fv_nb) {
									Area = fv_p.size();
									Area_centre = fv_p.getCentre();
									if (Area_centre.x > fvm->El().fv->m_Centre.x) relativ_pos = NB_E;
									if (Area_centre.x < fvm->El().fv->m_Centre.x) relativ_pos = NB_W;
									if (Area_centre.y > fvm->El().fv->m_Centre.y) relativ_pos = NB_N;
									if (Area_centre.y < fvm->El().fv->m_Centre.y) relativ_pos = NB_S;
									NBIsFound = true;
									break;
								}
							}
						}
						if (!NBIsFound) continue;
						//std::cout << relativ_pos << "\n";
						//std::cout << VCor_ed.size() << "  " << VCor_ed.size() + last_y_fvs.size() << "\n";
						int nb_id = fvm->nb_fvs[nb].fv->m_id;
						if ((relativ_pos == NB_N) and (fvm->El().fv->m_id >= FVMesh_VyA.size())) {
							nb_id = FVM.FiniteVolumeMesh[fvm->nb_fvs[nb].fv->m_id]->getPtrToSide(NB_N)->m_id;
						}
						if ((nb_id >= VCor_ed.size()) and (nb_id < VCor_ed.size() + last_y_fvs.size())) {
							nb_id = FVM.FiniteVolumeMesh[nb_id]->getPtrToSide(NB_N)->m_id;
						}
						//std::cout << " " << fvm->El().fv->m_id << " " << nb_id << " " << relativ_pos  << " "<< "\n"; // —“–¿ÕÕŒ


						int innerid = fvm->nb_fvs[nb].fv->m_id;
						fvm->nb_fvs[nb].innerid = innerid;


						// Pressure
						// dx*(P*P - P*N)
						if (relativ_pos == NB_N) {
							//std::cout << fvm->El().fv->m_id <<" - " << fvm->nb_fvs[nb].fv->m_id << "\n";
							fvm->fv_source += Area * (PApr[fvm->El().fv->m_id] - PApr[fvm->nb_fvs[nb].fv->m_id]);
						}

						if (FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
							vec_fv_to_nb = fvm->El().fv->m_Centre - FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_Centre;
							vec_fv_to_b = fvm->El().fv->m_Centre - Area_centre;
							vec_b_to_nb = Area_centre - FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_Centre;

							effD = (FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_DVisc * fvm->El().fv->m_DVisc) / ((FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_DVisc* vec_mod(vec_fv_to_b) / vec_mod(vec_fv_to_nb)) + (fvm->El().fv->m_DVisc *vec_mod(vec_b_to_nb) / vec_mod(vec_fv_to_nb)));
							Dnb_dir = (effD * Area) / vec_mod(vec_fv_to_nb);
							fv_prec_3 vel(fvm->El().fv->m_velocity.val); // {UCor_ed[i->m_id],VCor_ed[i->m_id],WCor_ed[i->m_id]}
							//std::cout << fvm->fv.innerid << " " << fvm->nb_fvs[nb].innerid << "\n";
							//std::cout << VCor_ed[fvm->fv.innerid] <<", " << VCor_ed[fvm->nb_fvs[nb].innerid] << " ";
							vel.y = 0.5*(VCor_ed[fvm->fv.innerid] + VCor_ed[fvm->nb_fvs[nb].innerid]);
							//std::cout << vel.y << "\n";

							if (relativ_pos == NB_N) vel.x = FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_velocity.dval.x;
							if (relativ_pos == NB_S) vel.x = fvm->El().fv->m_velocity.dval.x;
							//›“Œ —ÀŒ∆ÕŒ
							if (relativ_pos == NB_E) vel.x = 0.25*(fvm->El().fv->m_velocity.dval.x + FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_velocity.dval.x + FVM.FiniteVolumeMesh[fvm->El().fv->m_id]->getPtrToSide(NB_E)->m_velocity.dval.x + FVM.FiniteVolumeMesh[FVM.FiniteVolumeMesh[nb_id]->getPtrToSide(NB_E)->m_id]->FVPtr()->m_velocity.dval.x);
							if (relativ_pos == NB_W) vel.x = 0.25*(fvm->El().fv->m_velocity.dval.x + FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_velocity.dval.x + FVM.FiniteVolumeMesh[fvm->El().fv->m_id]->getPtrToSide(NB_W)->m_velocity.dval.x + FVM.FiniteVolumeMesh[FVM.FiniteVolumeMesh[nb_id]->getPtrToSide(NB_W)->m_id]->FVPtr()->m_velocity.dval.x);

							//Fnb_dir = Dens0 * (sqrt(pow(vel.x, 2) + pow(vel.y, 2) + pow(vel.z, 2))) * Area;
							Fnb_dir = Dens0 * vel.y * Area;
							// aNB							
							//{
							//	// Diffusion
							//	fvm->nb_fvs[nb].fv_coef += PeFunc(abs(Fnb_dir / Dnb_dir)) * Dnb_dir;
							//	// Convention
							//	if ((relativ_pos == NB_N) and Fnb_dir < 0) fvm->nb_fvs[nb].fv_coef -= Fnb_dir;
							//	if ((relativ_pos == NB_S) and Fnb_dir > 0) fvm->nb_fvs[nb].fv_coef += Fnb_dir;
							//}
							fvm->nb_fvs[nb].fv_coef += Dnb_dir;
							//std::cout << "fvm->nb_fvs[nb].fv_coef" << " " << fvm->nb_fvs[nb].fv_coef << "\n"; // —“–¿ÕÕŒ
						}
						else if (FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_type == FINITE_VOLUME_TYPE::FVT_BOUNDARY) {
							switch (relativ_pos)
							{
							case NB_E:
							case NB_W:
								// ›ÚÓ ‚ÓÁÏÓÊÌÓ ÚÓÊÂ ÌÂ Ì‡‰Ó
								vec_fv_to_nb = fvm->El().fv->m_Centre - FVM.FiniteVolumeMesh[nb_id]->FVPtr()->m_Centre;
								vec_fv_to_b = fvm->El().fv->m_Centre - Area_centre;
								effD = fvm->El().fv->m_DVisc;
								boundary_coefs += 1.0 / (1.0 / (effD* vec_mod(vec_fv_to_b) / vec_mod(vec_fv_to_nb)) + 1.0 / fvm->nb_fvs[nb].fv->fv_boundary->getVelocityCoefs().y)*Area;
								break;
							case NB_N:
							case NB_S:
								// œÓÍ‡ ÌÂ ÔÓÌˇÎ Ì‡‰Ó, Ë ÂÒÎË Ì‡‰Ó ÚÓ Í‡Í ÚÓ˜ÌÓ ÔÓÒ˜ËÚ‡¸? ≈ÒÎË ÌÂ„‰Â ÚÓ˜ÌÓ„Ó ‡Ò˜ÂÚ‡ ÌÂÚ, ÚÓ Ì‡‚ÂÌÓÂ ÌÂ Ì‡‰Ó
								break;
							default:
								break;
							}
						}
					}
					break;
				default:
					break;
				}
				FP0 = VCor0[fvm->fv.innerid]; // VCor_ed[fvm->El().fv->m_id]
				FNB0 = VCor0[fvm->nb_fvs[nb].innerid]; // VCor_ed[fvm->nb_fvs[nb].fv->m_id]
				// sigma*SUM(aNB)
				fvm->fv.fv_coef += fvm->nb_fvs[nb].fv_coef * sigma;  // 
				if (fvm->nb_fvs[nb].fv->m_type == FINITE_VOLUME_TYPE::FVT_BOUNDARY) {
					fvm->fv.fv_coef += boundary_coefs * sigma;
					fvm->fv.fv_coef -= FP0 * boundary_coefs * (1.0 - sigma);
					fvm->fv_source += boundary_coefs * fvm->nb_fvs[nb].fv->fv_boundary->getVelocity().y * sigma; // *(1 - sigma)
					fvm->fv_source += boundary_coefs * fvm->nb_fvs[nb].fv->fv_boundary->getVelocity().y * (1 - sigma);
				}
				// ‘P0*(1-sigma)*SUM(aNB)
				fvm->fv_source -= FP0 * fvm->nb_fvs[nb].fv_coef * (1 - sigma); // *(1 - sigma)
				// (1-sigma)*SUM(‘NB0 * aNB)
				fvm->fv_source += FNB0 * fvm->nb_fvs[nb].fv_coef * (1 - sigma); // *(1 - sigma)
			}
			// Spm*Vol
			fvm->fv.fv_coef -= fvm->sf.Spm.val * Volume;
			// ‘P0*(1-sigma)*Spm0*Vol
			fvm->fv_source -= FP0 * fvm->sf.Spm.dval * Volume*(1 - sigma);
			// (1-sigma)*Spp0*Vol
			fvm->fv_source += fvm->sf.Spp.dval * Volume*(1 - sigma);// *(1 - sigma)
			// sigma*Spp*Vol
			fvm->fv_source += fvm->sf.Spp.val * Volume*sigma;// *sigma

			// aP0
			if (dt == 0) { ap0 = 0.0; }
			else { ap0 = Dens0 * Volume / dt; }
			fvm->fv.fv_coef += ap0;
			// ‘P0*aP0
			FP0 = VCor0[fvm->fv.innerid]; // VCor_ed[fvm->El().fv->m_id]
			fvm->fv_source += FP0 * ap0;
		}
		std::cout << "PHI_PARAMS::UY\n";
		Solve(&FVMesh_VyA, &VApr.dval);


		std::cout << "UApr.dval: ";
		// Fn = AlfaF*F + (1-AlfaF)*Fn-1
		for (int i = 0;i < UApr.val.size();i++) {
			UApr.dval[i] = AlphaU * UApr.dval[i] + (1 - AlphaU)*UApr.val[i];
			std::cout << UApr.dval[i] << " ";
		}
		std::cout << "\n";
		std::cout << "VApr.dval: ";
		for (int i = 0;i < VApr.val.size();i++) {
			VApr.dval[i] = AlphaV * VApr.dval[i] + (1 - AlphaV)*VApr.val[i];
			std::cout << VApr.dval[i] << " ";
		}
		std::cout << "\n";


		// ap PP' = SUM(anb PNB') + bp* , bp* = -{(Rhop-Rhop0)/dt dVp + SUM(Fnb*)} 
		for (int i = 0;i < FVM.FiniteVolumeMesh.size();i++) {
			if (FVM.FiniteVolumeMesh[i]->FVPtr()->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
				FVMesh_PC.push_back(new FVMeshNode(*(FVM.FiniteVolumeMesh[i]), PHI_PARAMS::P));
			}
		}
		Dnb_dir = 0.0;
		Fnb_dir = 0.0;

		// ‘P*SUM(dens*dAnb*dnb) = SUM(‘NB*dens*dAnb*dnb) + (-(aP0 + dens * U*e *dAe - dens * U*w *dAw + dens * U*n *dAn - dens * U*s *dAs)
		//      +                           +                   +      +                +                  +                 +                                                        
		for (auto fvm : FVMesh_PC) {
			fv_prec ap0;
			fvm->fv.fv_coef = 0.0;
			Dens0 = fvm->El().fv->m_density.val;
			fv_prec de = 0.0; fv_prec dw = 0.0; fv_prec dn = 0.0; fv_prec ds = 0.0;
			fv_prec ue = 0.0; fv_prec uw = 0.0; fv_prec vn = 0.0; fv_prec vs = 0.0;
			for (int nb = 0;nb < fvm->nb_fvs.size();nb++) {
				FV_NB relativ_pos;
				boundary_coefs = 0.0;
				fv_prec dens = Dens0;
				switch (fvm->El().fv->FVdim) {
				case(FVA_L):
					break;
				case(FVA_S):
					Volume = fvm->El().fv->m_s.size();
					for (auto fv_p : fvm->El().fv->m_s.m_lineArray) {
						bool NBIsFound = false;
						if (fvm->nb_fvs[nb].fv->FVdim == FVA_L) {
							//std::cout << "{" << fv_p.m_p1.x << "," << fv_p.m_p1.y << "},{" << fv_p.m_p2.x << "," << fv_p.m_p2.y << "}-{" << nb.fv->m_l.m_p1.x << "," << nb.fv->m_l.m_p1.y << "},{" << nb.fv->m_l.m_p2.x << "," << nb.fv->m_l.m_p2.y << "}" << nb.fv->m_type << "\n";
							if (fv_p == fvm->nb_fvs[nb].fv->m_l) {
								Area = fv_p.size();
								Area_centre = fv_p.getCentre();
								if (Area_centre.x > fvm->El().fv->m_Centre.x) relativ_pos = NB_E;
								if (Area_centre.x < fvm->El().fv->m_Centre.x) relativ_pos = NB_W;
								if (Area_centre.y > fvm->El().fv->m_Centre.y) relativ_pos = NB_N;
								if (Area_centre.y < fvm->El().fv->m_Centre.y) relativ_pos = NB_S;
								NBIsFound = true;
							}
						}
						else if (fvm->nb_fvs[nb].fv->FVdim == FVA_S) {
							for (auto fv_nb : fvm->nb_fvs[nb].fv->m_s.m_lineArray) {
								if (fv_p == fv_nb) {
									Area = fv_p.size();
									Area_centre = fv_p.getCentre();
									if (Area_centre.x > fvm->El().fv->m_Centre.x) relativ_pos = NB_E;
									if (Area_centre.x < fvm->El().fv->m_Centre.x) relativ_pos = NB_W;
									if (Area_centre.y > fvm->El().fv->m_Centre.y) relativ_pos = NB_N;
									if (Area_centre.y < fvm->El().fv->m_Centre.y) relativ_pos = NB_S;
									NBIsFound = true;
									//std::cout << "{" << fv_p.m_p1.x << "," << fv_p.m_p1.y << "},{" << fv_p.m_p2.x << "," << fv_p.m_p2.y << "}-{" << fv_nb.m_p1.x << "," << fv_nb.m_p1.y << "},{" << fv_nb.m_p2.x << "," << fv_nb.m_p2.y << "}" << "\n";
									//std::cout << fvm->fv.fv->m_id << " " <<  fvm->nb_fvs[nb].fv->m_id << " TUT\n";
									break;
								}
							}
						}
						if (!NBIsFound) continue;
						fvm->nb_fvs[nb].innerid = fvm->nb_fvs[nb].fv->m_id;
						if (relativ_pos == NB_E) {
							int id = fvm->fv.fv->m_id;
							if (fvm->nb_fvs[nb].fv->m_type == FINITE_VOLUME_TYPE::FVT_BOUNDARY) {
								//std::cout << fvm->fv.fv->m_id << "e: " << fvm->nb_fvs[nb].fv->m_id << "\n";
								ue = fvm->nb_fvs[nb].fv->fv_boundary->getVelocity().x;
								de = 0.0;
							}
							else {
								//std::cout << fvm->fv.fv->m_id << "e: " << fvm->getPtrToSide(NB_E)->m_id - 1 - id / last_x_fvs.size() << "\n";
								ue = UApr.dval[fvm->getPtrToSide(NB_E)->m_id - 1 - id / last_x_fvs.size()];
								de = Area / FVMesh_VxA[fvm->getPtrToSide(NB_E)->m_id - 1 - id / last_x_fvs.size()]->fv.fv_coef;
								DensNB = fvm->getPtrToSide(NB_E)->m_density.val;
								//std::cout << FVMesh_VxA[fvm->getPtrToSide(NB_E)->m_id - 1 - id / last_x_fvs.size()]->fv.fv->m_id << " " << FVMesh_VxA[fvm->getPtrToSide(NB_E)->m_id - 1 - id / last_x_fvs.size()]->fv.fv_coef << "\n";
							}
							dens = 0.5*(Dens0 + DensNB);
							fvm->nb_fvs[nb].fv_coef = dens * Area*de;
						}
						else if (relativ_pos == NB_W) {
							if (fvm->fv.fv->m_id % last_x_fvs.size() == 0) {
								//std::cout << fvm->fv.fv->m_id << "w: " << FVMesh_PC[fvm->fv.fv->m_id]->getPtrToSide(NB_W)->m_id << "\n";
								uw = FVMesh_PC[fvm->fv.fv->m_id]->getPtrToSide(NB_W)->fv_boundary->getVelocity().x;
								dw = 0.0;
							}
							else {
								int id = fvm->getPtrToSide(NB_W)->m_id;
								//std::cout << fvm->fv.fv->m_id << "w: " << id - static_cast<int>(id / (last_x_fvs.size())) << "\n";
								uw = UApr.dval[id - static_cast<int>(id / (last_x_fvs.size()))];
								dw = Area / FVMesh_VxA[id - static_cast<int>(id / (last_x_fvs.size()))]->fv.fv_coef;
								DensNB = fvm->getPtrToSide(NB_W)->m_density.val;
								//std::cout << FVMesh_VxA[id - static_cast<int>(id / (last_x_fvs.size()))]->fv.fv->m_id << " " << FVMesh_VxA[id - static_cast<int>(id / (last_x_fvs.size()))]->fv.fv_coef << "\n";
							}
							dens = 0.5*(Dens0 + DensNB);
							fvm->nb_fvs[nb].fv_coef = dens * Area*dw;
						}
						else if (relativ_pos == NB_N) {
							if (fvm->fv.fv->m_id >= VApr.dval.size()) {
								//std::cout << fvm->fv.fv->m_id << "n: " << FVMesh_PC[fvm->fv.fv->m_id]->getPtrToSide(NB_N)->m_id << "\n";
								vn = FVMesh_PC[fvm->fv.fv->m_id]->getPtrToSide(NB_N)->fv_boundary->getVelocity().y;
								dn = 0.0;
							}
							else {
								//std::cout << fvm->fv.fv->m_id << "n: " << fvm->fv.fv->m_id << "\n";
								vn = VApr.dval[fvm->fv.fv->m_id];
								dn = Area / FVMesh_VxA[fvm->fv.fv->m_id]->fv.fv_coef;
								DensNB = fvm->getPtrToSide(NB_N)->m_density.val;
								//std::cout << FVMesh_VxA[fvm->fv.fv->m_id]->fv.fv->m_id << " " << FVMesh_VxA[fvm->fv.fv->m_id]->fv.fv_coef << "\n";
							}
							dens = 0.5*(Dens0 + DensNB);
							fvm->nb_fvs[nb].fv_coef = dens * Area*dn;
						}
						else if (relativ_pos == NB_S) {
							int id = fvm->fv.fv->m_id;
							if (id < last_y_fvs.size()) {
								//std::cout << fvm->fv.fv->m_id << "s: " << FVMesh_PC[id]->getPtrToSide(NB_S)->m_id << "\n";
								vs = FVMesh_PC[id]->getPtrToSide(NB_S)->fv_boundary->getVelocity().y;
								ds = 0.0;
							}
							else {
								//std::cout << fvm->fv.fv->m_id << "s: " << id - last_y_fvs.size() << "\n";
								vs = VApr.dval[id - last_y_fvs.size()];
								ds = Area / FVMesh_VxA[id - last_y_fvs.size()]->fv.fv_coef;
								DensNB = fvm->getPtrToSide(NB_S)->m_density.val;
								//std::cout << FVMesh_VxA[id - last_y_fvs.size()]->fv.fv->m_id << " " << FVMesh_VxA[id - last_y_fvs.size()]->fv.fv_coef << "\n";
							}
							dens = 0.5*(Dens0 + DensNB);
							fvm->nb_fvs[nb].fv_coef = dens * Area*ds;
						}
					}
					break;
				default:
					break;
				}
				//if (relativ_pos == NB_E) { std::cout << "de=" << de << " - " << "ue=" << ue << "\n"; }
				//else if (relativ_pos == NB_W) { std::cout << "dw=" << dw << " - " << "uw=" << uw << "\n"; }
				//else if (relativ_pos == NB_N) { std::cout << "dn=" << dn << " - " << "vn=" << vn << "\n"; }
				//else if (relativ_pos == NB_S) { std::cout << "ds=" << ds << " - " << "vs=" << vs << "\n"; }

				fvm->fv.fv_coef += fvm->nb_fvs[nb].fv_coef;
				fvm->fv_source += dens * Area * (uw - ue + vs - vn);

				de = 0.0;  dw = 0.0;  dn = 0.0;  ds = 0.0;
				ue = 0.0;  uw = 0.0;  vn = 0.0;  vs = 0.0;
			}

			if (dt == 0) { ap0 = 0.0; }
			else { ap0 = (Dens - Dens0) * Volume / dt; }
			fvm->fv_source -= ap0;
		}
		std::cout << "PHI_PARAMS::P\n";


		sum.val = sum.dval;
		sum.dval = 0.0;
		for (auto err : FVMesh_PC) {
			sum.dval += err->fv_source;
		}
		std::cout << "sum.val: " << sum.val << ", sum.dval: " << sum.dval << "\n";

		SolveWithFirstZero(&FVMesh_PC, &Pcor);

		// ue = ue* +de(PP'- PE')
		std::cout << "UCor_ed: ";
		for (int i = 0;i < UApr.dval.size();i++) {
			//int num = i + i / (last_x_fvs.size() - 1);
			//FVM.FiniteVolumes[num]->m_velocity.dval.x = UApr[i] + FVMesh_VxA[i]->getCoef(FVMesh_VxA[i]->fv.fv) * (Pcor[i] - Pcor[FVMesh_VyA[i]->getPtrToSide(NB_E)->m_id]);
			UCor_ed[i] = UApr.dval[i] + FVMesh_VxA[i]->getCoef(FVMesh_VxA[i]->fv.fv) * (Pcor[i] - Pcor[FVMesh_VyA[i]->getPtrToSide(NB_E)->m_id]);
			std::cout << UCor_ed[i] << " ";
			UApr.val[i] = UApr.dval[i];
		}
		std::cout << "\n";
		// vn = vn* +dn(PP'- PN')
		std::cout << "VCor_ed: ";
		for (int i = 0;i < VApr.dval.size();i++) {
			//FVM.FiniteVolumes[i]->m_velocity.dval.y = VApr[i] + FVMesh_VyA[i]->getCoef(FVMesh_VyA[i]->fv.fv) * (Pcor[i] - Pcor[FVMesh_VyA[i]->getPtrToSide(NB_N)->m_id]);
			VCor_ed[i] = VApr.dval[i] + FVMesh_VyA[i]->getCoef(FVMesh_VyA[i]->fv.fv) * (Pcor[i] - Pcor[FVMesh_VyA[i]->getPtrToSide(NB_N)->m_id]);
			std::cout << VCor_ed[i] << " ";
			VApr.val[i] = VApr.dval[i];
		}
		std::cout << "\n";


		//ÒÍÓÓÒÚË
		std::cout << "\n";
		std::cout << "Velocity.x: ";
		for (auto fvm : FVM.FiniteVolumeMesh) {
			if (fvm->FVPtr()->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
				fv_prec velE = 0.0;
				fv_prec velW = 0.0;
				//std::cout << (fvm->FVPtr()->m_id) / (last_x_fvs.size()) << " ";
				//std::cout << fvm->FVPtr()->m_id << ": ";
				if (fvm->getPtrToSide(NB_W)->fv_boundary != nullptr) {
					velW = fvm->getPtrToSide(NB_W)->fv_boundary->getVelocity().x;
					//std::cout << fvm->getPtrToSide(NB_W)->m_id << " ";
				}
				else {
					velW = UCor_ed[fvm->getPtrToSide(NB_W)->m_id - (fvm->FVPtr()->m_id) / (last_x_fvs.size())];
					//std::cout << fvm->getPtrToSide(NB_W)->m_id- (fvm->FVPtr()->m_id) / (last_x_fvs.size()) << " ";
				}
				if (fvm->getPtrToSide(NB_E)->fv_boundary != nullptr) {
					velE = fvm->getPtrToSide(NB_E)->fv_boundary->getVelocity().x;
					//std::cout << fvm->getPtrToSide(NB_E)->m_id << " ";
				}
				else {
					velE = UCor_ed[fvm->FVPtr()->m_id - (fvm->FVPtr()->m_id) / (last_x_fvs.size())];
					//std::cout << fvm->FVPtr()->m_id - (fvm->FVPtr()->m_id) / (last_x_fvs.size()) << " ";
				}
				fvm->FVPtr()->m_velocity.dval.x = 0.5*(velW + velE);
				//std::cout << "\n";
				std::cout << fvm->FVPtr()->m_velocity.dval.x << " ";
			}
		}
		std::cout << "\n";
		std::cout << "Velocity.y: ";
		for (auto fvm : FVM.FiniteVolumeMesh) {
			if (fvm->FVPtr()->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
				fv_prec velN = 0.0;
				fv_prec velS = 0.0;
				//std::cout << (fvm->FVPtr()->m_id) / (last_x_fvs.size()) << " ";
				//std::cout << fvm->FVPtr()->m_id << ": ";
				if (fvm->getPtrToSide(NB_S)->fv_boundary != nullptr) {
					velS = fvm->getPtrToSide(NB_S)->fv_boundary->getVelocity().y;
					//std::cout << fvm->getPtrToSide(NB_S)->m_id << " ";
				}
				else {
					velS = VCor_ed[fvm->getPtrToSide(NB_S)->m_id];
					//std::cout << fvm->getPtrToSide(NB_S)->m_id << " ";
				}
				if (fvm->getPtrToSide(NB_N)->fv_boundary != nullptr) {
					velN = fvm->getPtrToSide(NB_N)->fv_boundary->getVelocity().y;
					//std::cout << fvm->getPtrToSide(NB_N)->m_id << " ";
				}
				else {
					velN = VCor_ed[fvm->FVPtr()->m_id];
					//std::cout << fvm->FVPtr()->m_id << " ";
				}
				fvm->FVPtr()->m_velocity.dval.y = 0.5*(velS + velN);
				std::cout << fvm->FVPtr()->m_velocity.dval.y << " ";
				//std::cout << "\n";
				//std::cout << velS << " " << velN << "\n";
			}
		}
		std::cout << "\n";


		// P = P* + P'

		std::cout << "Pressure: ";
		for (auto i : FVM.FiniteVolumes) {
			if (i->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
				//i->m_pressure.val = i->m_pressure.dval;
				i->m_pressure.dval = PApr[i->m_id] + AlphaP * Pcor[i->m_id];
				std::cout << i->m_pressure.dval << " ";
			}
		}
		std::cout << "\n";
		std::cout << "_________________________________________________ \n";

		// P* = P
		for (auto i : FVM.FiniteVolumes) {
			if (i->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
				PApr[i->m_id] = i->m_pressure.dval;
			}
		}
		for (auto i : FVMesh_VxA) { delete i; }
		FVMesh_VxA.resize(0);
		for (auto i : FVMesh_VyA) { delete i; }
		FVMesh_VyA.resize(0);
		for (auto i : FVMesh_PC) { delete i; }
		FVMesh_PC.resize(0);

		std::cout << "Err: " << abs(sum.dval - sum.val) << "\n";

	}




	// U0 = U; ÕŒ Õ¿ƒŒ œ≈–≈—“–Œ»“‹ «Õ¿◊≈Õ»ﬂ Œ“ —Ã≈Ÿ≈ÕÕ€’  Œ   ÕŒ–Ã¿À‹Õ€Ã ƒÀﬂ — Œ–Œ—“»
	// P0 = P
	for (auto i : FVM.FiniteVolumes) {
		if (i->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
			i->m_pressure.val = i->m_pressure.dval;
			i->m_velocity.val = i->m_velocity.dval;
		}
	}
}
*/


void FVM_CD::timeUpdate(cd_prec dt){
	if (m_options.firstCycle == false) {
		if (computationalDomain_mode == MODE_CD::CD_DEBUG) {
			std::cout << "SPH_CD::timeIntegration\n";
		}


		SIMPLE_ALGO(0.0);

		/*
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
					
					//// ue = ue* +de(PP'- PE')
					//for (auto i : FVM.FiniteVolumes) {
					//	i->m_velocity.dval.x = VelApr[i->m_id].x + de * (Pcor[i->m_id] - Pcor[FVM.FiniteVolumeMesh[i->m_id].E]);
					//}
					//
					//// vn = vn* +dn(PP'- PN')
					//for (auto i : FVM.FiniteVolumes) {
					//	i->m_velocity.dval.y = VelApr[i->m_id].y + dn * (Pcor[i->m_id] - Pcor[FVM.FiniteVolumeMesh[i->m_id].N]);
					//}
					
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