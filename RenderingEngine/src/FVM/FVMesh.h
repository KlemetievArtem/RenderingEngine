#pragma once
enum MESH_POSITON {
	ORIGNAL,
	SHIFTED_RIGHT,
	SHIFTED_UP
};
class FVMeshNode;


enum FV_NB {
	NB_E,
	NB_W,
	NB_N,
	NB_S,
};

class FVMeshNodeEmpty {
	FiniteVolume* fv;
	std::vector<FiniteVolume*> nb_fvs;
public:
	FVMeshNodeEmpty(FiniteVolume* fvp) : fv(fvp) { }
	void addNeighbour(FiniteVolume* fvp) {
		nb_fvs.push_back(fvp);
	}
	friend class FVMeshNode;
	FiniteVolume* FVPtr() {
		return fv;
	}
	std::vector<FiniteVolume*> NB_FVSPrt() {
		return nb_fvs;
	}
	FiniteVolume* getPtrToSide(FV_NB s) {
		for (auto nb : nb_fvs) {
			if ((nb->m_Centre.x > fv->m_Centre.x) and s == NB_E) return nb;
			if ((nb->m_Centre.x < fv->m_Centre.x) and s == NB_W) return nb;
			if ((nb->m_Centre.y > fv->m_Centre.y) and s == NB_N) return nb;
			if ((nb->m_Centre.y < fv->m_Centre.y) and s == NB_S) return nb;
		}
	}
};




struct SourceFunction {
	// S(Ô) = Spp - Spm*Ô 
	Ch_<fv_prec> Spm;
	Ch_<fv_prec> Spp;
	SourceFunction(fv_prec spm = 0.0 , fv_prec spp = 0.0) {
		Spm.val = spm;
		Spp.val = spp;
	}
};

struct MeshElement {
	FiniteVolume* fv;
	fv_prec fv_coef ;
	int innerid;
};


enum FPARAM {
	FP,
	FU,
	FV,
	FW
};

class FVMeshNode {
public:
	FPARAM type;
	MeshElement fv;
	std::vector<MeshElement> nb_fvs;
	PHI_PARAMS phi;

	fv_prec fv_source = 0.0;

	SourceFunction sf;
	FVMeshNode(FiniteVolume* fvp, FPARAM t) { fv.fv = fvp; }
	void addNeighbour(FiniteVolume* fvp) {
		nb_fvs.push_back(MeshElement({ fvp,0.0,0 }));
	}


	void clear() {
		fv.fv_coef = 0.0;
		for (auto nb : nb_fvs) {
			nb.fv_coef = 0.0;
		}
		fv_source = 0.0;
	}


	FVMeshNode(FVMeshNodeEmpty FVMesh, PHI_PARAMS add_phi) : phi(add_phi) {
		fv.fv = FVMesh.fv;
		for (int i = 0; i < FVMesh.nb_fvs.size(); i++) {
			nb_fvs.push_back({FVMesh.nb_fvs[i],0.0});
		}
	}
	PHI_PARAMS getPhiName() {
		return phi;
	}
	/*
	void calcCoefsAndSources(cd_prec dt = 0.0, fv_prec sigma = 1.0) {
		fv_prec effD;
		fv_prec_3 vec_distance_to_nb;
		fv_prec_3 vec_distance_to_border;
		fv_prec Area;
		fv_prec Volume;
		fv_prec_3 Area_centre;
		fv_prec boundary_coefs;

		for (int i = 0;i < nb_fvs.size();i++) {
			fv_prec ap0;
			if (dt == 0) { ap0 = 0.0; }
			else { ap0 = fv->m_density.val *Volume / dt; }
			fv_coef = ap0;
			switch (phi) {
			case(PHI_PARAMS::UX):
				fv_source = fv->m_velocity.val.x * ap0;
				break;
			case(PHI_PARAMS::UY):
				fv_source = fv->m_velocity.val.y * ap0;
				break;
			case(PHI_PARAMS::P):
				fv_source = fv->m_pressure.val * ap0;
				break;
			default:
				break;
			}

			fv_prec Dnb_dir = 0.0;
			fv_prec Fnb_dir = 0.0;
			boundary_coefs = 0.0;
			switch (fv->FVdim) {
			case(FVA_S):
				Volume = fv->m_s.size();
				for (int dim = 1; dim <= 2;dim++) {
					for (auto fv_p : fv->m_s.m_lineArray) {
						for (auto fv_nb : nb_fvs[i]->m_s.m_lineArray) {
							if (fv_p.ConnectedWith(fv_nb)) {
								Area = fv_p.size();
								Area_centre = fv_p.getCentre();
							}
						}
					}
					switch (phi) {
					case(PHI_PARAMS::UX):
					case(PHI_PARAMS::UY):
						if (nb_fvs[i]->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
							effD = nb_fvs[i]->m_DVisc;
							vec_distance_to_nb = fv->m_Centre - nb_fvs[i]->m_Centre;
							Dnb_dir = (effD * Area) / (sqrt(pow(vec_distance_to_nb.x, 2) + pow(vec_distance_to_nb.y, 2) + pow(vec_distance_to_nb.z, 2)));
							Fnb_dir = fv->m_density.val * fv->m_velocity.val[dim - 1] * Area;
							//std::cout << "Dnb_dir:" << Dnb_dir << ",Fnb_dir:"  << Fnb_dir << " \n";
							//Diffusion
							nb_fv_coefs[i] += PeFunc(abs(Fnb_dir / Dnb_dir)) * Dnb_dir;
							//Convection
							if ((fv->m_Centre[dim - 1] > Area_centre[dim - 1]) and Fnb_dir > 0) {
								nb_fv_coefs[i] += Fnb_dir;
							}
							else if ((fv->m_Centre[dim - 1] < Area_centre[dim - 1]) and Fnb_dir < 0) {
								nb_fv_coefs[i] += Fnb_dir;
							}
						}
						if (nb_fvs[i]->m_type == FINITE_VOLUME_TYPE::FVT_BOUNDARY) {
							vec_distance_to_nb = fv->m_Centre - nb_fvs[i]->m_Centre;
							vec_distance_to_border = fv->m_Centre - Area_centre;
							boundary_coefs += 1.0 / (1.0 / (effD*(sqrt(pow(vec_distance_to_border.x, 2) + pow(vec_distance_to_border.y, 2) + pow(vec_distance_to_border.z, 2))) / (sqrt(pow(vec_distance_to_nb.x, 2) + pow(vec_distance_to_nb.y, 2) + pow(vec_distance_to_nb.z, 2)))) + 1.0 / nb_fvs[i]->fv_boundary->getVelocityCoefs()[dim - 1])*Area;
						}
						break;
					case(PHI_PARAMS::P):
						break;
					default:
						break;
					}
				}
				switch (phi) {
				case(PHI_PARAMS::UX):
				case(PHI_PARAMS::UY):
					nb_fv_coefs[i] += Area * (fv->m_pressure.val - nb_fvs[i]->m_pressure.val);
					break;
				case(PHI_PARAMS::P):
					break;
				default:
					break;
				}



				break;
			default:
				break;
			}
			fv_coef += nb_fv_coefs[i] * sigma;  // 
			if (nb_fvs[i]->m_type == FINITE_VOLUME_TYPE::FVT_BOUNDARY) {
				fv_coef += boundary_coefs * sigma;
			}
			fv_prec FP0;
			fv_prec FNB0;
			switch (phi) {
			case(PHI_PARAMS::UX):
				if (nb_fvs[i]->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
					FP0 = fv->m_velocity.val.x;
					FNB0 = nb_fvs[i]->m_velocity.val.x;
					fv_source -= FP0 * nb_fv_coefs[i] * (1 - sigma); // *(1 - sigma)
					fv_source += FNB0 * nb_fv_coefs[i] * (1 - sigma); // *(1 - sigma)
				}
				if (nb_fvs[i]->m_type == FINITE_VOLUME_TYPE::FVT_BOUNDARY) {
					FP0 = nb_fvs[i]->fv_boundary->getVelocity().x;
					fv_source += FP0 * boundary_coefs;
				}
				break;
			case(PHI_PARAMS::UY):
				if (nb_fvs[i]->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
					FP0 = fv->m_velocity.val.y;
					FNB0 = nb_fvs[i]->m_velocity.val.y;
					fv_source -= FP0 * nb_fv_coefs[i] * (1 - sigma); // *(1 - sigma)
					fv_source += FNB0 * nb_fv_coefs[i] * (1 - sigma); // *(1 - sigma)
				}
				if (nb_fvs[i]->m_type == FINITE_VOLUME_TYPE::FVT_BOUNDARY) {
					FP0 = nb_fvs[i]->fv_boundary->getVelocity().y;
					fv_source += FP0 * boundary_coefs;
				}
				break;
			case(PHI_PARAMS::P):
				if (nb_fvs[i]->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
					FP0 = fv->m_pressure.val;
					FNB0 = nb_fvs[i]->m_pressure.val;
					fv_source -= FP0 * nb_fv_coefs[i] * (1 - sigma); // *(1 - sigma)
					fv_source += FNB0 * nb_fv_coefs[i] * (1 - sigma); // *(1 - sigma)
				}
				if (nb_fvs[i]->m_type == FINITE_VOLUME_TYPE::FVT_BOUNDARY) {
					FP0 = nb_fvs[i]->m_pressure.val;
					fv_source += FP0 * boundary_coefs;
				}
				break;
			default:
				break;
			}


			fv_coef += sf.Spm.val * Volume;
			switch (phi) {
			case(PHI_PARAMS::UX):
				fv_source -= fv->m_velocity.val.x * sf.Spm.dval * Volume;
				break;
			case(PHI_PARAMS::UY):
				fv_source -= fv->m_velocity.val.x * sf.Spm.dval * Volume;
				break;
			case(PHI_PARAMS::P):
				fv_source -= fv->m_pressure.val * sf.Spm.dval * Volume;
				break;
			default:
				break;
			}
			fv_source += sf.Spp.dval * Volume*(1 - sigma);// *(1 - sigma)
			fv_source += sf.Spp.val * Volume*sigma;// *sigma

		}



		
		//switch (phi) {
		//case(PHI_PARAMS::UX):
		//case(PHI_PARAMS::UY):
		//	for (int i = 0;i < nb_fvs.size();i++) {
		//		switch (fv->FVATYPE) {
		//		case(FVA_P):
		//			break;
		//		case(FVA_L):
		//			for (auto fv_p : fv->m_s.m_lineArray) {
		//				for (auto fv_nb : nb_fvs[i]->m_s.m_lineArray) {
		//					if (fv_p.ConnectedWith(fv_nb)) {
		//						Area = fv_p.size();
		//					}
		//				}
		//			}
		//			break;
		//		case(FVA_S):
		//			break;
		//		default:
		//			break;
		//		}
		//		effD = nb_fvs[i]->m_DVisc;
		//		vec_distance = fv->m_Centre - nb_fvs[i]->m_Centre;
		//		//std::cout << (sqrt(pow(vec_distance.x, 2) + pow(vec_distance.y, 2) + pow(vec_distance.z, 2))) << "\n";
		//		nb_fv_coefs[i] = (effD * Area) / (sqrt(pow(vec_distance.x, 2) + pow(vec_distance.y, 2) + pow(vec_distance.z, 2)));
		//		//std::cout << nb_fv_coefs[i] << "\n";
		//		fv_coef += nb_fv_coefs[i];
		//	}
		//	break;
		//default:
		//	break;
		//}
		
	}

	void calcCoefsAndSources(cd_prec dt = 0.0, fv_prec sigma = 1.0) {
		fv_prec effD;
		fv_prec_3 vec_distance_to_nb;
		fv_prec_3 vec_distance_to_border;
		fv_prec Area;
		fv_prec Volume;
		fv_prec_3 Area_centre;
		fv_prec boundary_coefs;

		fv_prec Dnb_dir = 0.0;
		fv_prec Fnb_dir = 0.0;

		fv_prec FP0;
		fv_prec FNB0;

		for (int i = 0;i < nb_fvs.size();i++) {
			fv_prec ap0;
			if (dt == 0) { ap0 = 0.0; }
			else { ap0 = fv->m_density.val *Volume / dt; }
			fv_coef = ap0;
			switch (phi) {
			case(PHI_PARAMS::UX):
				fv_source = fv->m_velocity.val.x * ap0;
				FP0 = fv->m_velocity.val.x;
				FNB0 = nb_fvs[i]->m_velocity.val.x;
			case(PHI_PARAMS::UY):
				fv_source = fv->m_velocity.val.x * ap0;

				boundary_coefs = 0.0;
				switch (fv->FVdim) {
				case(FVA_S):
					Volume = fv->m_s.size();
					for (int dim = 1; dim <= 2;dim++) {
						for (auto fv_p : fv->m_s.m_lineArray) {
							for (auto fv_nb : nb_fvs[i]->m_s.m_lineArray) {
								if (fv_p.ConnectedWith(fv_nb)) {
									Area = fv_p.size();
									Area_centre = fv_p.getCentre();
								}
							}
						}
						if (nb_fvs[i]->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
							effD = nb_fvs[i]->m_DVisc;
							vec_distance_to_nb = fv->m_Centre - nb_fvs[i]->m_Centre;
							Dnb_dir = (effD * Area) / (sqrt(pow(vec_distance_to_nb.x, 2) + pow(vec_distance_to_nb.y, 2) + pow(vec_distance_to_nb.z, 2)));
							Fnb_dir = fv->m_density.val * fv->m_velocity.val[dim - 1] * Area;
							//std::cout << "Dnb_dir:" << Dnb_dir << ",Fnb_dir:"  << Fnb_dir << " \n";
							//Diffusion
							nb_fv_coefs[i] += PeFunc(abs(Fnb_dir / Dnb_dir)) * Dnb_dir;
							//Convection
							if ((fv->m_Centre[dim - 1] > Area_centre[dim - 1]) and Fnb_dir > 0) {
								nb_fv_coefs[i] += Fnb_dir;
							}
							else if ((fv->m_Centre[dim - 1] < Area_centre[dim - 1]) and Fnb_dir < 0) {
								nb_fv_coefs[i] += Fnb_dir;
							}
						}
						if (nb_fvs[i]->m_type == FINITE_VOLUME_TYPE::FVT_BOUNDARY) {
							vec_distance_to_nb = fv->m_Centre - nb_fvs[i]->m_Centre;
							vec_distance_to_border = fv->m_Centre - Area_centre;
							boundary_coefs += 1.0 / (1.0 / (effD*(sqrt(pow(vec_distance_to_border.x, 2) + pow(vec_distance_to_border.y, 2) + pow(vec_distance_to_border.z, 2))) / (sqrt(pow(vec_distance_to_nb.x, 2) + pow(vec_distance_to_nb.y, 2) + pow(vec_distance_to_nb.z, 2)))) + 1.0 / nb_fvs[i]->fv_boundary->getVelocityCoefs()[dim - 1])*Area;
						}
					}
					if (phi == PHI_PARAMS::UX) {
						if(fv->m_Centre.x < nb_fvs[i]->m_Centre.x) fv_source += Area * (fv->m_pressure.val - nb_fvs[i]->m_pressure.val);
						else if (fv->m_Centre.x > nb_fvs[i]->m_Centre.x) fv_source += Area * (nb_fvs[i]->m_pressure.val - fv->m_pressure.val);
					}
					else if(phi == PHI_PARAMS::UY){
						if (fv->m_Centre.y < nb_fvs[i]->m_Centre.y) fv_source += Area * (fv->m_pressure.val - nb_fvs[i]->m_pressure.val);
						else if (fv->m_Centre.y > nb_fvs[i]->m_Centre.y) fv_source += Area * (nb_fvs[i]->m_pressure.val - fv->m_pressure.val);
					}
					break;
				default:
					break;
				}
				break;
				fv_coef += nb_fv_coefs[i] * sigma;  // 
				if (nb_fvs[i]->m_type == FINITE_VOLUME_TYPE::FVT_BOUNDARY) {
					fv_coef += boundary_coefs * sigma;
				}
				FP0 = fv->m_velocity.val.y;
				FNB0 = nb_fvs[i]->m_velocity.val.y;
				fv_source -= FP0 * nb_fv_coefs[i] * (1 - sigma); // *(1 - sigma)
				fv_source += FNB0 * nb_fv_coefs[i] * (1 - sigma); // *(1 - sigma)

				fv_coef += sf.Spm.val * Volume;
				fv_source -= FP0 * sf.Spm.dval * Volume;

				fv_source += sf.Spp.dval * Volume*(1 - sigma);// *(1 - sigma)
				fv_source += sf.Spp.val * Volume*sigma;// *sigma
			case(PHI_PARAMS::P):
				fv_source = fv->m_pressure.val * ap0;

				boundary_coefs = 0.0;
				switch (fv->FVdim) {
				case(FVA_S):
					Volume = fv->m_s.size();
					for (int dim = 1; dim <= 2;dim++) {
						for (auto fv_p : fv->m_s.m_lineArray) {
							for (auto fv_nb : nb_fvs[i]->m_s.m_lineArray) {
								if (fv_p.ConnectedWith(fv_nb)) {
									Area = fv_p.size();
									Area_centre = fv_p.getCentre();
								}
							}
						}
						if (nb_fvs[i]->m_type == FINITE_VOLUME_TYPE::FVT_DEFAULT) {
							effD = nb_fvs[i]->m_DVisc;
							vec_distance_to_nb = fv->m_Centre - nb_fvs[i]->m_Centre;
							Dnb_dir = (effD * Area) / (sqrt(pow(vec_distance_to_nb.x, 2) + pow(vec_distance_to_nb.y, 2) + pow(vec_distance_to_nb.z, 2)));
							Fnb_dir = fv->m_density.val * fv->m_velocity.val[dim - 1] * Area;
							//std::cout << "Dnb_dir:" << Dnb_dir << ",Fnb_dir:"  << Fnb_dir << " \n";
							//Diffusion
							nb_fv_coefs[i] += PeFunc(abs(Fnb_dir / Dnb_dir)) * Dnb_dir;
							//Convection
							if ((fv->m_Centre[dim - 1] > Area_centre[dim - 1]) and Fnb_dir > 0) {
								nb_fv_coefs[i] += Fnb_dir;
							}
							else if ((fv->m_Centre[dim - 1] < Area_centre[dim - 1]) and Fnb_dir < 0) {
								nb_fv_coefs[i] += Fnb_dir;
							}
						}
						if (nb_fvs[i]->m_type == FINITE_VOLUME_TYPE::FVT_BOUNDARY) {
							vec_distance_to_nb = fv->m_Centre - nb_fvs[i]->m_Centre;
							vec_distance_to_border = fv->m_Centre - Area_centre;
							boundary_coefs += 1.0 / (1.0 / (effD*(sqrt(pow(vec_distance_to_border.x, 2) + pow(vec_distance_to_border.y, 2) + pow(vec_distance_to_border.z, 2))) / (sqrt(pow(vec_distance_to_nb.x, 2) + pow(vec_distance_to_nb.y, 2) + pow(vec_distance_to_nb.z, 2)))) + 1.0 / nb_fvs[i]->fv_boundary->getVelocityCoefs()[dim - 1])*Area;
						}
					}
					if (phi == PHI_PARAMS::UX) {
						if (fv->m_Centre.x < nb_fvs[i]->m_Centre.x) fv_source += Area * (fv->m_pressure.val - nb_fvs[i]->m_pressure.val);
						else if (fv->m_Centre.x > nb_fvs[i]->m_Centre.x) fv_source += Area * (nb_fvs[i]->m_pressure.val - fv->m_pressure.val);
					}
					else if (phi == PHI_PARAMS::UY) {
						if (fv->m_Centre.y < nb_fvs[i]->m_Centre.y) fv_source += Area * (fv->m_pressure.val - nb_fvs[i]->m_pressure.val);
						else if (fv->m_Centre.y > nb_fvs[i]->m_Centre.y) fv_source += Area * (nb_fvs[i]->m_pressure.val - fv->m_pressure.val);
					}
					break;
				default:
					break;
				}


				break;
			default:
				break;
			}
		}
	}

	void SpecialSource(){
		fv_prec Area;
		fv_prec PP = fv->m_pressure.val;
		fv_prec PNB;
		fv_prec_3 distance;
		fv_prec_3 OX(1.0, 0.0, 0.0);
		fv_prec_3 OY(0.0, 1.0, 0.0);
		fv_prec_3 OZ(0.0, 0.0, 1.0);
		for (int i = 0;i < nb_fvs.size();i++) {
			switch (fv->FVdim) {
			case(FVA_L):
				break;
			case(FVA_S):
				for (auto fv_p : fv->m_s.m_lineArray) {
					for (auto fv_nb : nb_fvs[i]->m_s.m_lineArray) {
						if (fv_p.ConnectedWith(fv_nb)) {
							distance = nb_fvs[i]->m_Centre;
							Area = fv_p.size();
							PNB = nb_fvs[i]->m_pressure.val; 
						}
					}
				}
				break;
			case(FVA_V):
				break;
			default:
				break;
			}
			distance = -fv->m_Centre;
			fv_source = Area * (PP - PNB);
		}
		cd_prec cos_OX_dist_ = ((OX.x*(distance.x)) + (OX.y*(distance.y)) + (OX.z*(distance.z))) / (sqrt(pow(OX.x, 2) + pow(OX.y, 2) + pow(OX.z, 2))*sqrt(pow(distance.x, 2) + pow(distance.y, 2) + pow(distance.z, 2)));

		switch (phi) {
		case(PHI_PARAMS::UX):
			fv_source *= cos_OX_dist_;
			break;
		case(PHI_PARAMS::UY):
			fv_source *= sqrt(1-pow(cos_OX_dist_,2));
			break;
		default:
			break;
		}
		//fv_source = 1.0;
	}

	*/

	FiniteVolume* getPtrToSide(FV_NB s) {
		for (auto nb : nb_fvs) {
			if ((nb.fv->m_Centre.x > fv.fv->m_Centre.x) and s == NB_E) return nb.fv;
			if ((nb.fv->m_Centre.x < fv.fv->m_Centre.x) and s == NB_W) return nb.fv;
			if ((nb.fv->m_Centre.y > fv.fv->m_Centre.y) and s == NB_N) return nb.fv;
			if ((nb.fv->m_Centre.y < fv.fv->m_Centre.y) and s == NB_S) return nb.fv;
		}
		return nullptr;
	}

	fv_prec getCoefToSide(FV_NB s) {
		for (auto nb : nb_fvs) {
			if ((nb.fv->m_Centre.x > fv.fv->m_Centre.x) and s == NB_E) return nb.fv_coef;
			if ((nb.fv->m_Centre.x < fv.fv->m_Centre.x) and s == NB_W) return nb.fv_coef;
			if ((nb.fv->m_Centre.y > fv.fv->m_Centre.y) and s == NB_N) return nb.fv_coef;
			if ((nb.fv->m_Centre.y < fv.fv->m_Centre.y) and s == NB_S) return nb.fv_coef;
		}
	}
	

	fv_prec getSource(){
		return fv_source;
	}

	std::vector<MeshElement>* NbEls() {
		return &nb_fvs;
	}
	MeshElement El() {
		return fv;
	}

	std::vector<FiniteVolume*> getNbPointers() {
		std::vector<FiniteVolume*> retVal;
		for (auto nb : nb_fvs) {
			retVal.push_back(nb.fv);
		}
		return retVal;
	}
	FiniteVolume* getVFPointers() {
		return fv.fv;
	}
	void setInnerId(FiniteVolume* fv_ptr, int id) {
		if (fv_ptr == fv.fv) return;
		for (int i = 0;i < nb_fvs.size(); i++) {
			if (fv_ptr == nb_fvs[i].fv) {
				nb_fvs[i].innerid = id;
			}
		}
	}

	fv_prec getCoef(FiniteVolume* fv_ptr) {
		if (fv_ptr == fv.fv) return fv.fv_coef;
		for (int i = 0;i < nb_fvs.size(); i++) {
			if (fv_ptr == nb_fvs[i].fv) {
				return nb_fvs[i].fv_coef;
			}
		}
		return 0.0;
	}
	fv_prec residual = 0.0;
	//fv_prec residualScale = 0.0;
	void RecalcResidual() {
		residual = fv_source;
		switch (type){
		case(FPARAM::FU):
			for (auto fvnb : nb_fvs) residual += fvnb.fv_coef * fvnb.fv->m_velocity.dval.x;
			residual += fv.fv_coef * fv.fv->m_velocity.dval.x;
			break;
		case(FPARAM::FV):
			for (auto fvnb : nb_fvs) residual += fvnb.fv_coef * fvnb.fv->m_velocity.dval.y;
			residual += fv.fv_coef * fv.fv->m_velocity.dval.y;
			break;
		default:
			break;
		}
	}

};


class tempFVMesh {

};