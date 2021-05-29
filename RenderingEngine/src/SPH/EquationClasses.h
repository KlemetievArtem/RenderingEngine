#pragma once

class ConcreteEquation {
public:
	virtual void Func(std::vector<part_prec>* val, ParticlePair * pp, PARTICLEPAIR_NUMERATIN ppn) = 0;
	void Calc(std::vector<part_prec>* val, ParticlePair * pp, PAIR_INTERACTION_TYPE pit) {
		switch (pit) {
		case PIT_MAIN:
			this->Func(val, pp, PPN_MAIN);
			break;
		case PIT_NEIGHBOUR:
			this->Func(val, pp, PPN_NEIGHBOUR);
			break;
		case PIT_BOTH:
			this->Func(val, pp, PPN_MAIN);
			this->Func(val, pp, PPN_NEIGHBOUR);
			break;
		case PIT_NONE:
			break;
		default:
			break;
		}
	}
	virtual void Func(std::vector<part_prec_3>* vector_val, ParticlePair * pp, DIMENSIONS dim, PARTICLEPAIR_NUMERATIN ppn) = 0;
	void Calc(std::vector<part_prec_3>* vector_val, ParticlePair * pp, DIMENSIONS dim, PAIR_INTERACTION_TYPE pit) {
		switch (pit) {
		case PIT_MAIN:
			this->Func(vector_val, pp, dim, PPN_MAIN);
			break;
		case PIT_NEIGHBOUR:
			this->Func(vector_val, pp, dim, PPN_NEIGHBOUR);
			break;
		case PIT_BOTH:
			this->Func(vector_val, pp, dim, PPN_MAIN);
			this->Func(vector_val, pp, dim, PPN_NEIGHBOUR);
			break;
		case PIT_NONE:
			break;
		default:
			break;
		}
	}

};
class ConcreteTimeUpdateEquation {
public:
	virtual void FuncT(Particle* p, part_prec t) = 0;
	void Calc(Particle* p, part_prec t) {
		this->FuncT(p, t);
	}

	virtual void FuncT(Particle* p, part_prec_3 val) = 0;
	void Calc3(Particle* p, part_prec_3 val) {
		this->FuncT(p, val);
	}

};
class TimeEquation : public EquationsBuilder {
public:
	TimeEquation(ConcreteTimeUpdateEquation* t) : _concreteEquation(t) { }
	void changeConcreteEquation(ConcreteTimeUpdateEquation* t) {
		delete _concreteEquation;
		_concreteEquation = t;
	}
	void Calc(Particle* p, cd_prec t) {
		if (ParticleTypesInUse == 1) {
			if (usedParticlesActivity[0] == ACTIVE) {
				if (p->m_type == usedParticles[0]) {
					if (_concreteEquation != nullptr) _concreteEquation->Calc(p, t);
				}
			}
		}
		if (ParticleTypesInUse == 2) {
			bool condition = false;
			for (int k = 0; k < ParticleTypesInUse; k++) {
				if (usedParticlesActivity[k] == ACTIVE) {
					condition = condition or (p->m_type == usedParticlesActivity[k]);
				}
				else if (usedParticlesActivity[k] == PASSIVE) {
					condition = condition and (p->m_type == usedParticlesActivity[k]);
				}
			}
			if (condition) {
				if (_concreteEquation != nullptr) _concreteEquation->Calc(p, t);
			}
		}
	}
	void Calc(Particle* p, part_prec_3 val) {
		if (ParticleTypesInUse == 1) {
			if (usedParticlesActivity[0] == ACTIVE) {
				if (p->m_type == usedParticles[0]) {
					if (_concreteEquation != nullptr) _concreteEquation->Calc3(p, val);
				}
			}
		}
		if (ParticleTypesInUse == 2) {
			bool condition = false;
			for (int k = 0; k < ParticleTypesInUse; k++) {
				if (usedParticlesActivity[k] == ACTIVE) {
					condition = condition or (p->m_type == usedParticlesActivity[k]);
				}
				else if (usedParticlesActivity[k] == PASSIVE) {
					condition = condition and (p->m_type == usedParticlesActivity[k]);
				}
			}
			if (condition) {
				if (_concreteEquation != nullptr) _concreteEquation->Calc3(p, val);
			}
		}
	}

	ConcreteTimeUpdateEquation* _concreteEquation;
	~TimeEquation() {
		delete _concreteEquation;
	}
};
//  LOGIC CUALITY OF LIFE (SPH.ParticlePairs-loop)
// MAIN PARTICLE IS 
#define PPL_MAIN_PARTICLE_IS(type) (i->Mpart_ptr->m_type == type)
#define PPL_MAIN_PARTICLE_IS_NOT(type) (i->Mpart_ptr->m_type != type)
	// NEIGBOUR PARTICLE IS 
#define PPL_NEIGHBOUR_PARTICLE_IS(type) (i->NBpart_ptr->m_type == type)
#define PPL_NEIGHBOUR_PARTICLE_IS_NOT(type) (i->NBpart_ptr->m_type != type)


#define PPL_BOTH_PARTICLES_ARE(type) ((i->NBpart_ptr->m_type == type) and (i->Mpart_ptr->m_type == type))
#define PPL_BOTH_PARTICLES_ARE_NOT(type) ((i->NBpart_ptr->m_type != type) and (i->Mpart_ptr->m_type != type))

#define PPL_ONE_OF_PARTICLE_IS_1_OTHER_IS_2(type1,type2) (((i->NBpart_ptr->m_type == type1) and (i->Mpart_ptr->m_type == type2)) or ((i->NBpart_ptr->m_type == type2) and (i->Mpart_ptr->m_type == type1)))


class Equation : public EquationsBuilder {
public:
	Equation(ConcreteEquation* c) : _concreteEquation(c) { }
	void changeConcreteEquation(ConcreteEquation* c) {
		delete _concreteEquation;
		_concreteEquation = c;
	}

	void Calc(std::vector<part_prec>* val, ParticlePair* i) {
		if (_concreteEquation == nullptr) assert("Equation::_concreteEquation==" && 0);
		PAIR_INTERACTION_TYPE interactionType;
		if (ParticleTypesInUse == 1) {
			switch (usedParticlesActivity[0]) {
			case(ACTIVE):
				//std::cout << "A    ";
				interactionType = PIT_BOTH;
				break;
			case(REACTIVE):
				//std::cout << "R    ";
				interactionType = PIT_NONE;
				break;
			case(PASSIVE):
				//std::cout << "P    ";
				interactionType = PIT_NONE;
				break;
			default:
				break;
			}
			_concreteEquation->Calc(val, i, interactionType);
			//std::cout << i->Mpart_ptr->m_type << "  " << i->NBpart_ptr->m_type << "  " << interactionType << "\n";
		}
		else if (ParticleTypesInUse == 2) {
			bool exit_condition = false;
			for (int k = 0; k < ParticleTypesInUse; k++) {
				if (exit_condition) break;
				switch (usedParticlesActivity[k]) {
				case(ACTIVE):
					//std::cout << "ACTIVE ";
					for (int kk = k; kk < ParticleTypesInUse; kk++) {
						switch (usedParticlesActivity[kk]) {
						case(ACTIVE):
							if ((i->Mpart_ptr->m_type == usedParticles[kk]) and (i->NBpart_ptr->m_type == usedParticles[k])) {
								//std::cout << "A A    ";
								interactionType = PIT_BOTH;
								exit_condition = true;
							}
							break;
						case(REACTIVE):
							if ((i->Mpart_ptr->m_type == usedParticles[kk]) and (i->NBpart_ptr->m_type == usedParticles[k])) {
								//std::cout << "A R    ";
								interactionType = PIT_MAIN;
								exit_condition = true;
							}
							if ((i->NBpart_ptr->m_type == usedParticles[kk]) and (i->Mpart_ptr->m_type == usedParticles[k])) {
								//std::cout << "A R    ";
								interactionType = PIT_NEIGHBOUR;
								exit_condition = true;
							}
							break;
						case(PASSIVE):
							if ((i->Mpart_ptr->m_type == usedParticles[kk]) and (i->NBpart_ptr->m_type == usedParticles[k])) {
								//std::cout << "A P    ";
								interactionType = PIT_NEIGHBOUR;
								exit_condition = true;
							}
							if ((i->NBpart_ptr->m_type == usedParticles[kk]) and (i->Mpart_ptr->m_type == usedParticles[k])) {
								//std::cout << "A P    ";
								interactionType = PIT_MAIN;
								exit_condition = true;
							}
							break;
						default:
							break;
						}
					}
					break;
				case(REACTIVE):
					//std::cout << "REACTIVE ";
					for (int kk = k; kk < ParticleTypesInUse; kk++) {
						switch (usedParticlesActivity[kk]) {
						case(ACTIVE):
							//std::cout << "R A    ";
							if ((i->Mpart_ptr->m_type == usedParticles[kk]) and (i->NBpart_ptr->m_type == usedParticles[k])) {
								interactionType = PIT_NEIGHBOUR;
								exit_condition = true;
							}
							if ((i->NBpart_ptr->m_type == usedParticles[kk]) and (i->Mpart_ptr->m_type == usedParticles[k])) {
								interactionType = PIT_MAIN;
								exit_condition = true;
							}
							break;
						case(REACTIVE):
							//std::cout << "R R    ";
							if ((i->Mpart_ptr->m_type == usedParticles[kk]) and (i->NBpart_ptr->m_type == usedParticles[k])) {
								interactionType = PIT_NONE;
								exit_condition = true;
							}
							break;
						case(PASSIVE):
							//std::cout << "R P    ";
							if ((i->Mpart_ptr->m_type == usedParticles[kk]) and (i->NBpart_ptr->m_type == usedParticles[k])) {
								//interactionType = PIT_NONE;
								interactionType = PIT_NEIGHBOUR;
								exit_condition = true;
							}
							if ((i->NBpart_ptr->m_type == usedParticles[kk]) and (i->Mpart_ptr->m_type == usedParticles[k])) {
								//interactionType = PIT_NONE;
								interactionType = PIT_MAIN;
								exit_condition = true;
							}
							break;
						default:
							break;
						}
					}
					break;
				case(PASSIVE):
					//std::cout << "PASSIVE ";
					for (int kk = k; kk < ParticleTypesInUse; kk++) {
						switch (usedParticlesActivity[kk]) {
						case(ACTIVE):
							//std::cout << "P A    ";
							if ((i->Mpart_ptr->m_type == usedParticles[kk]) and (i->NBpart_ptr->m_type == usedParticles[k])) {
								interactionType = PIT_MAIN;
								exit_condition = true;
							}
							if ((i->NBpart_ptr->m_type == usedParticles[kk]) and (i->Mpart_ptr->m_type == usedParticles[k])) {
								interactionType = PIT_NEIGHBOUR;
								exit_condition = true;
							}
							break;
						case(REACTIVE):
							//std::cout << "P R    ";
							if ((i->Mpart_ptr->m_type == usedParticles[kk]) and (i->NBpart_ptr->m_type == usedParticles[k])) {
								//interactionType = PIT_NONE;
								interactionType = PIT_MAIN;
								exit_condition = true;
							}
							if ((i->NBpart_ptr->m_type == usedParticles[kk]) and (i->Mpart_ptr->m_type == usedParticles[k])) {
								//interactionType = PIT_NONE;
								interactionType = PIT_NEIGHBOUR;
								exit_condition = true;
							}
							break;
						case(PASSIVE):
							//std::cout << "P P    ";
							if ((i->Mpart_ptr->m_type == usedParticles[kk]) and (i->NBpart_ptr->m_type == usedParticles[k])) {
								interactionType = PIT_NONE;
								exit_condition = true;
							}
							break;
						default:
							break;
						}
					}
					break;
				default:
					break;
				}
			}
			_concreteEquation->Calc(val, i, interactionType);
		}
		else if (ParticleTypesInUse == 3) {
			bool exit_condition = false;
			for (int k = 0; k < ParticleTypesInUse; k++) {
				if (exit_condition) break;
				switch (usedParticlesActivity[k]) {
				case(ACTIVE):
					//std::cout << "ACTIVE ";
					for (int kk = k; kk < ParticleTypesInUse; kk++) {
						switch (usedParticlesActivity[kk]) {
						case(ACTIVE):
							if ((i->Mpart_ptr->m_type == usedParticles[kk]) and (i->NBpart_ptr->m_type == usedParticles[k])) {
								//std::cout << "A A    ";
								interactionType = PIT_BOTH;
								exit_condition = true;
							}
							break;
						case(REACTIVE):
							if ((i->Mpart_ptr->m_type == usedParticles[kk]) and (i->NBpart_ptr->m_type == usedParticles[k])) {
								//std::cout << "A R    ";
								interactionType = PIT_MAIN;
								exit_condition = true;
							}
							if ((i->NBpart_ptr->m_type == usedParticles[kk]) and (i->Mpart_ptr->m_type == usedParticles[k])) {
								//std::cout << "A R    ";
								interactionType = PIT_NEIGHBOUR;
								exit_condition = true;
							}
							break;
						case(PASSIVE):
							if ((i->Mpart_ptr->m_type == usedParticles[kk]) and (i->NBpart_ptr->m_type == usedParticles[k])) {
								//std::cout << "A P    ";
								interactionType = PIT_NEIGHBOUR;
								exit_condition = true;
							}
							if ((i->NBpart_ptr->m_type == usedParticles[kk]) and (i->Mpart_ptr->m_type == usedParticles[k])) {
								//std::cout << "A P    ";
								interactionType = PIT_MAIN;
								exit_condition = true;
							}
							break;
						default:
							break;
						}
					}
					break;
				case(REACTIVE):
					//std::cout << "REACTIVE ";
					for (int kk = k; kk < ParticleTypesInUse; kk++) {
						switch (usedParticlesActivity[kk]) {
						case(ACTIVE):
							//std::cout << "R A    ";
							if ((i->Mpart_ptr->m_type == usedParticles[kk]) and (i->NBpart_ptr->m_type == usedParticles[k])) {
								interactionType = PIT_NEIGHBOUR;
								exit_condition = true;
							}
							if ((i->NBpart_ptr->m_type == usedParticles[kk]) and (i->Mpart_ptr->m_type == usedParticles[k])) {
								interactionType = PIT_MAIN;
								exit_condition = true;
							}
							break;
						case(REACTIVE):
							//std::cout << "R R    ";
							if ((i->Mpart_ptr->m_type == usedParticles[kk]) and (i->NBpart_ptr->m_type == usedParticles[k])) {
								interactionType = PIT_NONE;
								exit_condition = true;
							}
							break;
						case(PASSIVE):
							//std::cout << "R P    ";
							if ((i->Mpart_ptr->m_type == usedParticles[kk]) and (i->NBpart_ptr->m_type == usedParticles[k])) {
								//interactionType = PIT_NONE;
								interactionType = PIT_NEIGHBOUR;
								exit_condition = true;
							}
							if ((i->NBpart_ptr->m_type == usedParticles[kk]) and (i->Mpart_ptr->m_type == usedParticles[k])) {
								//interactionType = PIT_NONE;
								interactionType = PIT_MAIN;
								exit_condition = true;
							}
							break;
						default:
							break;
						}
					}
					break;
				case(PASSIVE):
					//std::cout << "PASSIVE ";
					for (int kk = k; kk < ParticleTypesInUse; kk++) {
						switch (usedParticlesActivity[kk]) {
						case(ACTIVE):
							//std::cout << "P A    ";
							if ((i->Mpart_ptr->m_type == usedParticles[kk]) and (i->NBpart_ptr->m_type == usedParticles[k])) {
								interactionType = PIT_MAIN;
								exit_condition = true;
							}
							if ((i->NBpart_ptr->m_type == usedParticles[kk]) and (i->Mpart_ptr->m_type == usedParticles[k])) {
								interactionType = PIT_NEIGHBOUR;
								exit_condition = true;
							}
							break;
						case(REACTIVE):
							//std::cout << "P R    ";
							if ((i->Mpart_ptr->m_type == usedParticles[kk]) and (i->NBpart_ptr->m_type == usedParticles[k])) {
								//interactionType = PIT_NONE;
								interactionType = PIT_MAIN;
								exit_condition = true;
							}
							if ((i->NBpart_ptr->m_type == usedParticles[kk]) and (i->Mpart_ptr->m_type == usedParticles[k])) {
								//interactionType = PIT_NONE;
								interactionType = PIT_NEIGHBOUR;
								exit_condition = true;
							}
							break;
						case(PASSIVE):
							//std::cout << "P P    ";
							if ((i->Mpart_ptr->m_type == usedParticles[kk]) and (i->NBpart_ptr->m_type == usedParticles[k])) {
								interactionType = PIT_NONE;
								exit_condition = true;
							}
							break;
						default:
							break;
						}
					}
					break;
				default:
					break;
				}
			}
			_concreteEquation->Calc(val, i, interactionType);
		}
	}

	void Calc(std::vector<part_prec_3>* vector_val, ParticlePair* i, DIMENSIONS dim) {
		if (_concreteEquation == nullptr) assert("Equation::_concreteEquation==" && 0);
		PAIR_INTERACTION_TYPE interactionType;
		if (ParticleTypesInUse == 1) {
			switch (usedParticlesActivity[0]) {
			case(ACTIVE):
				//std::cout << "A    ";
				interactionType = PIT_BOTH;
				break;
			case(REACTIVE):
				//std::cout << "R    ";
				interactionType = PIT_NONE;
				break;
			case(PASSIVE):
				//std::cout << "P    ";
				interactionType = PIT_NONE;
				break;
			default:
				break;
			}
			_concreteEquation->Calc(vector_val, i, dim, interactionType);
		}
		else if (ParticleTypesInUse == 2) {
			bool exit_condition = false;
			for (int k = 0; k < ParticleTypesInUse; k++) {
				if (exit_condition) break;
				switch (usedParticlesActivity[k]) {
				case(ACTIVE):
					//std::cout << "ACTIVE ";
					for (int kk = k; kk < ParticleTypesInUse; kk++) {
						switch (usedParticlesActivity[kk]) {
						case(ACTIVE):
							if ((i->Mpart_ptr->m_type == usedParticles[kk]) and (i->NBpart_ptr->m_type == usedParticles[k])) {
								//std::cout << "A A    ";
								interactionType = PIT_BOTH;
								exit_condition = true;
							}
							break;
						case(REACTIVE):
							if ((i->Mpart_ptr->m_type == usedParticles[kk]) and (i->NBpart_ptr->m_type == usedParticles[k])) {
								//std::cout << "A R    ";
								interactionType = PIT_MAIN;
								exit_condition = true;
							}
							if ((i->NBpart_ptr->m_type == usedParticles[kk]) and (i->Mpart_ptr->m_type == usedParticles[k])) {
								//std::cout << "A R    ";
								interactionType = PIT_NEIGHBOUR;
								exit_condition = true;
							}
							break;
						case(PASSIVE):
							if ((i->Mpart_ptr->m_type == usedParticles[kk]) and (i->NBpart_ptr->m_type == usedParticles[k])) {
								//std::cout << "A P    ";
								interactionType = PIT_NEIGHBOUR;
								exit_condition = true;
							}
							if ((i->NBpart_ptr->m_type == usedParticles[kk]) and (i->Mpart_ptr->m_type == usedParticles[k])) {
								//std::cout << "A P    ";
								interactionType = PIT_MAIN;
								exit_condition = true;
							}
							break;
						default:
							break;
						}
					}
					break;
				case(REACTIVE):
					//std::cout << "REACTIVE ";
					for (int kk = k; kk < ParticleTypesInUse; kk++) {
						switch (usedParticlesActivity[kk]) {
						case(ACTIVE):
							//std::cout << "R A    ";
							if ((i->Mpart_ptr->m_type == usedParticles[kk]) and (i->NBpart_ptr->m_type == usedParticles[k])) {
								interactionType = PIT_NEIGHBOUR;
								exit_condition = true;
							}
							if ((i->NBpart_ptr->m_type == usedParticles[kk]) and (i->Mpart_ptr->m_type == usedParticles[k])) {
								interactionType = PIT_MAIN;
								exit_condition = true;
							}
							break;
						case(REACTIVE):
							//std::cout << "R R    ";
							if ((i->Mpart_ptr->m_type == usedParticles[kk]) and (i->NBpart_ptr->m_type == usedParticles[k])) {
								interactionType = PIT_NONE;
								exit_condition = true;
							}
							break;
						case(PASSIVE):
							//std::cout << "R P    ";
							if ((i->Mpart_ptr->m_type == usedParticles[kk]) and (i->NBpart_ptr->m_type == usedParticles[k])) {
								interactionType = PIT_NONE;
								exit_condition = true;
							}
							//if (i->Mpart_ptr->m_type == usedParticles[kk]) interactionType = PIT_NEIGHBOUR;
							//if (i->NBpart_ptr->m_type == usedParticles[kk]) interactionType = PIT_MAIN;
							break;
						default:
							break;
						}
					}
					break;
				case(PASSIVE):
					//std::cout << "PASSIVE ";
					for (int kk = k; kk < ParticleTypesInUse; kk++) {
						switch (usedParticlesActivity[kk]) {
						case(ACTIVE):
							//std::cout << "P A    ";
							if ((i->Mpart_ptr->m_type == usedParticles[kk]) and (i->NBpart_ptr->m_type == usedParticles[k])) {
								interactionType = PIT_MAIN;
								exit_condition = true;
							}
							if ((i->NBpart_ptr->m_type == usedParticles[kk]) and (i->Mpart_ptr->m_type == usedParticles[k])) {
								interactionType = PIT_NEIGHBOUR;
								exit_condition = true;
							}
							break;
						case(REACTIVE):
							//std::cout << "P R    ";
							if ((i->Mpart_ptr->m_type == usedParticles[kk]) and (i->NBpart_ptr->m_type == usedParticles[k])) {
								interactionType = PIT_NONE;
								exit_condition = true;
							}
							//if (i->Mpart_ptr->m_type == usedParticles[kk]) interactionType = PIT_MAIN;
							//if (i->NBpart_ptr->m_type == usedParticles[kk]) interactionType = PIT_NEIGHBOUR;
							break;
						case(PASSIVE):
							//std::cout << "P P    ";
							if ((i->Mpart_ptr->m_type == usedParticles[kk]) and (i->NBpart_ptr->m_type == usedParticles[k])) {
								interactionType = PIT_NONE;
								exit_condition = true;
							}
							break;
						default:
							break;
						}
					}
					break;
				default:
					break;
				}
			}
			_concreteEquation->Calc(vector_val, i, dim, interactionType);
		}
		else if (ParticleTypesInUse == 3) {
			bool exit_condition = false;
			for (int k = 0; k < ParticleTypesInUse; k++) {
				if (exit_condition) break;
				switch (usedParticlesActivity[k]) {
				case(ACTIVE):
					//std::cout << "ACTIVE ";
					for (int kk = k; kk < ParticleTypesInUse; kk++) {
						switch (usedParticlesActivity[kk]) {
						case(ACTIVE):
							if ((i->Mpart_ptr->m_type == usedParticles[kk]) and (i->NBpart_ptr->m_type == usedParticles[k])) {
								//std::cout << "A A    ";
								interactionType = PIT_BOTH;
								exit_condition = true;
							}
							break;
						case(REACTIVE):
							if ((i->Mpart_ptr->m_type == usedParticles[kk]) and (i->NBpart_ptr->m_type == usedParticles[k])) {
								//std::cout << "A R    ";
								interactionType = PIT_MAIN;
								exit_condition = true;
							}
							if ((i->NBpart_ptr->m_type == usedParticles[kk]) and (i->Mpart_ptr->m_type == usedParticles[k])) {
								//std::cout << "A R    ";
								interactionType = PIT_NEIGHBOUR;
								exit_condition = true;
							}
							break;
						case(PASSIVE):
							if ((i->Mpart_ptr->m_type == usedParticles[kk]) and (i->NBpart_ptr->m_type == usedParticles[k])) {
								//std::cout << "A P    ";
								interactionType = PIT_NEIGHBOUR;
								exit_condition = true;
							}
							if ((i->NBpart_ptr->m_type == usedParticles[kk]) and (i->Mpart_ptr->m_type == usedParticles[k])) {
								//std::cout << "A P    ";
								interactionType = PIT_MAIN;
								exit_condition = true;
							}
							break;
						default:
							break;
						}
					}
					break;
				case(REACTIVE):
					//std::cout << "REACTIVE ";
					for (int kk = k; kk < ParticleTypesInUse; kk++) {
						switch (usedParticlesActivity[kk]) {
						case(ACTIVE):
							//std::cout << "R A    ";
							if ((i->Mpart_ptr->m_type == usedParticles[kk]) and (i->NBpart_ptr->m_type == usedParticles[k])) {
								interactionType = PIT_NEIGHBOUR;
								exit_condition = true;
							}
							if ((i->NBpart_ptr->m_type == usedParticles[kk]) and (i->Mpart_ptr->m_type == usedParticles[k])) {
								interactionType = PIT_MAIN;
								exit_condition = true;
							}
							break;
						case(REACTIVE):
							//std::cout << "R R    ";
							if ((i->Mpart_ptr->m_type == usedParticles[kk]) and (i->NBpart_ptr->m_type == usedParticles[k])) {
								interactionType = PIT_NONE;
								exit_condition = true;
							}
							break;
						case(PASSIVE):
							//std::cout << "R P    ";
							if ((i->Mpart_ptr->m_type == usedParticles[kk]) and (i->NBpart_ptr->m_type == usedParticles[k])) {
								interactionType = PIT_NONE;
								exit_condition = true;
							}
							//if (i->Mpart_ptr->m_type == usedParticles[kk]) interactionType = PIT_NEIGHBOUR;
							//if (i->NBpart_ptr->m_type == usedParticles[kk]) interactionType = PIT_MAIN;
							break;
						default:
							break;
						}
					}
					break;
				case(PASSIVE):
					//std::cout << "PASSIVE ";
					for (int kk = k; kk < ParticleTypesInUse; kk++) {
						switch (usedParticlesActivity[kk]) {
						case(ACTIVE):
							//std::cout << "P A    ";
							if ((i->Mpart_ptr->m_type == usedParticles[kk]) and (i->NBpart_ptr->m_type == usedParticles[k])) {
								interactionType = PIT_MAIN;
								exit_condition = true;
							}
							if ((i->NBpart_ptr->m_type == usedParticles[kk]) and (i->Mpart_ptr->m_type == usedParticles[k])) {
								interactionType = PIT_NEIGHBOUR;
								exit_condition = true;
							}
							break;
						case(REACTIVE):
							//std::cout << "P R    ";
							if ((i->Mpart_ptr->m_type == usedParticles[kk]) and (i->NBpart_ptr->m_type == usedParticles[k])) {
								interactionType = PIT_NONE;
								exit_condition = true;
							}
							//if (i->Mpart_ptr->m_type == usedParticles[kk]) interactionType = PIT_MAIN;
							//if (i->NBpart_ptr->m_type == usedParticles[kk]) interactionType = PIT_NEIGHBOUR;
							break;
						case(PASSIVE):
							//std::cout << "P P    ";
							if ((i->Mpart_ptr->m_type == usedParticles[kk]) and (i->NBpart_ptr->m_type == usedParticles[k])) {
								interactionType = PIT_NONE;
								exit_condition = true;
							}
							break;
						default:
							break;
						}
					}
					break;
				default:
					break;
				}
			}
			_concreteEquation->Calc(vector_val, i, dim, interactionType);
		}
		//std::cout << i->Mpart_ptr->m_type << "  " << i->NBpart_ptr->m_type << "  " << interactionType << "\n";
	}


	ConcreteEquation* _concreteEquation;
	~Equation() {
		delete _concreteEquation;
	}
};



