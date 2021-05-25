#pragma once

#include "virtualComputationalDomain.h"
#include "Particle.h"
//#include "FiniteParticleMethodFULL.h"

enum NEIGBOURS_SEARCH_ALGORITHM {
	DIRECT,
	UNIFORM_GRID
};
enum DISTRIBUTION {
	RANDOM,
	UNIFORM
};
enum TIME_INTEGRATION_SCHEME {
	EXPLICIT,
	IMPLICIT,
	SEMI_IMPICIT,
	SECOND_ORDER_SCEME,
	NONE
};
enum VAL_CHANGE_USING{
	VAL,
	DVAL_DT
};

enum GLOBAL_METHOD {
	STANDART_SPH,
	CORRECTIVE_SPH,
	FINITE_PARTICLE_METHOD,
	THE_ONLY_RIGHT
};

enum BOUNDARIES_TYPES {
	FIRSTORDER_FIRSTORDER,
	PERIODIC_PERIODIC,
	PERIODIC_FIRSTORDER,
	FIRSTORDER_PERIODIC
};

enum BOUNDARY_HANDLING {
	DUMMY_PARTICLES,
	MIRROR_PARTICLES,
	REPULSIVE_FORCE,
	RENORMALIZATION
};


struct ParticleRendererBuffer {
	part_prec_3 position;
	glm::vec3 color;
	part_prec size;
};

enum PARTICLE_TYPE_ACTIVITY {
	ACTIVE,
	PASSIVE,
	REACTIVE
};

enum PAIR_INTERACTION_TYPE {
	PIT_NONE,
	PIT_MAIN,
	PIT_NEIGHBOUR,
	PIT_BOTH
};

class EquationsBuilder {
public:
	PARTICLETYPE usedParticles[ALLTYPES];
	PARTICLE_TYPE_ACTIVITY usedParticlesActivity[ALLTYPES];
	unsigned int ParticleTypesInUse=0;

	void addParticleType(PARTICLETYPE pt, PARTICLE_TYPE_ACTIVITY pta) {
		usedParticles[ParticleTypesInUse] = pt;
		usedParticlesActivity[ParticleTypesInUse] = pta;
		ParticleTypesInUse++;
	}
};
class Equation;
class TimeEquation;

// Передавать одни и теже опции в SPH и в пары
struct SPH_OPTIONS {
	unsigned int nrOfParticles[ALLTYPES];

	NEIGBOURS_SEARCH_ALGORITHM NBSAlg = UNIFORM_GRID; // DIRECT	UNIFORM_GRID
	//DIMENSIONS KernelDimension = D1;

	part_prec smoothingKernelLengthCoefficient = 2.0;

	GLOBAL_METHOD SPH_algorithm = THE_ONLY_RIGHT;
	GLOBAL_METHOD SPH_iteration_algorithm = FINITE_PARTICLE_METHOD;


	VAL_CHANGE_USING densityChangeUsing = DVAL_DT; // VAL DVAL_DT
	int densityChangeAlgorithm = 2;
	int velocityChangeAlgorithm_pressurePart = 2;
	int velocityChangeAlgorithm_viscosityPart = 2;

	DISTRIBUTION distributionREAL = UNIFORM; // RANDOM UNIFORM
	DISTRIBUTION distributionBOUNDARY = UNIFORM; // RANDOM UNIFORM
	TIME_INTEGRATION_SCHEME timeIntegrationScheme = SECOND_ORDER_SCEME; // EXPLICIT IMPLICIT SEMI_IMPICIT SECOND_ORDER_SCEME NONE


	BOUNDARY_HANDLING boundary_handling = MIRROR_PARTICLES; // MIRROR_PARTICLES RENORMALIZATION

	// true  false
	bool firstCycle = true;
	bool cornerVP = true; 
	bool IsStab = false;
	bool Density_Diffusion_term = false;
	bool BoundaryCalculation = false;


	part_prec_3 average_dim_steps;

	Equation* ContinuityEquation;
	TimeEquation* DensityUpdate;
	Equation* VelocityEquation_PressurePart;
	Equation* VelocityEquation_ViscosityPart;
	TimeEquation* VelocityUpdate;
	TimeEquation* PositionUpdate;

	//Equation* BoundaryDensity;
	//Equation* CorrectionEquation;
	TimeEquation* BoundaryDensityUpdate;
	

};

#include "ParticlePair.h"
#include "BoundaryMentor.h"
#include "UniformTreeGrid.h"

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
			this->Func(vector_val, pp, dim,PPN_MAIN);
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
		//std::cout << i->Mpart_ptr->m_type << "  " << i->NBpart_ptr->m_type << "  " << interactionType << "\n";
	}


	ConcreteEquation* _concreteEquation;
	~Equation() {
		delete _concreteEquation;
	}
};

#include "EquationClasses.h"

struct SPH_ESENTIALS {
	std::vector<Particle*> Particles;
	std::vector<ParticlePair*> ParticlePairs;
	~SPH_ESENTIALS() {
		for (auto*&i : Particles)
			delete i;
		for (auto*&i : ParticlePairs)
			delete i;
	}
};

class SPH_CD : public CompDomain {
private:
	SPH_ESENTIALS SPH;
	SPH_OPTIONS m_options;

	bool FirstSycle = true;

	part_prec BH_D;
	part_prec BH_r;
	BoundaryMentor BM;
	//std::unique_ptr<BoundaryMentor> BM;
	std::vector<ParticleRendererBuffer> PRB;

	uniformTreeGrid_2D* UG;

public:
	SPH_CD(SPH_OPTIONS options, DIMENSIONS dim) {
		setTypeTo(MICROCOSME::MC_SPH);
		nrOfDim = dim;
		m_options = options;
	}
	void Initilization();
	//void InitialRendering(std::vector<Mesh*>* meshes);

	//void UpdateRendering(std::vector<Model*>* models);
	void UpdateRendering(std::vector<Model*>* models, Texture* tex, Texture* tex_specualar, std::vector<Material*>* materials);
	void AfterRendering(std::vector<Model*>* models);

	void timeUpdate(cd_prec dt);

	void timeStep(cd_prec dt);
	void timeStep_thread(cd_prec dt, std::atomic<bool>& dataReadyForRender, std::atomic<bool>& dataIsRendering);

	void neighbourSearch();

	void uniformGridPartitionInitialization();
	void creatingVirtualParticles();
	void addingVirtualTOPairs();
	void deletingParticlePairs();

	void DeletingVirtualParticles();


	void DensityVariation();

	void ExternalForces();

	void SaveMaxVelocity();

	void Boundary_Handling();


	bool chekParticlePairsFor(Particle* p1, Particle* p2);

	void PPC_Density(ParticlePair* pp, std::vector<part_prec>* drhodt);
	void PPC_InternalForces(ParticlePair* pp, std::vector<part_prec_3>* dvdt);
	
	// В идеале
	//  function approximation                1-st derivative approximation           2-nd derivative approximation
	//
	//  Sum(Wij*Vj) = 1                       Sum(Wij,x*Vj) = 0                       Sum(Wij,xx*Vj) = 0
	//  Sum((xj-xi)*Wij*Vj) = 0				  Sum((xj-xi)*Wij,x*Vj) = 1				  Sum((xj-xi)*Wij,xx*Vj) = 0
	//  Sum((xj-xi)^2*Wij*Vj) = 0			  Sum((xj-xi)^2*Wij,x*Vj) = 0			  Sum((xj-xi)^2*Wij,xx*Vj) = 2
	//         ...							        ...								        ...
	//         ...							        ...								        ...
	//  Sum((xj-xi)^n*Wij*Vj) = 0			  Sum((xj-xi)^n*Wij,x*Vj) = 0             Sum((xj-xi)^n*Wij,xx*Vj) = 0
	//
	//FINITE PARTICLE METHOD
	//STANDART                                                                                                          в идеале матрица должна быть следующего вида:
	//  | Sum(fj*Wij*Vj)   |     | fi   |   | Sum(Wij*Vj)      Sum((xj-xi)*Wij*Vj)      Sum((yj-yi)*Wij*Vj)    |                      | 1    0   0  |
	//  | Sum(fj*Wij,x*Vj) |  =  | fi,x | * | Sum(Wij,x*Vj)    Sum((xj-xi)*Wij,x*Vj)    Sum((yj-yi)*Wij,x*Vj)  |					  | 0    1   ?  |
	//  | Sum(fj*Wij,y*Vj) |     | fi,y |   | Sum(Wij,y*Vj)    Sum((xj-xi)*Wij,y*Vj)    Sum((yj-yi)*Wij,y*Vj)  |					  | 0    ?   1  |
	void PPC_FPMStandart_General(ParticlePair* pp, std::vector<Matrix>* FPM_matrix);
	void PPC_FPMStandart_vec(ParticlePair* pp, std::vector<std::vector<part_prec>>* FPM_matrix, std::string param);

	//Difference
	//  | Sum((fj-fi)*Wij*Vj)   |     |   0  |   | 1    Sum((xj-xi)*Wij*Vj)      Sum((yj-yi)*Wij*Vj)    |
	//  | Sum((fj-fi)*Wij,x*Vj) |  =  | fi,x | * | 1    Sum((xj-xi)*Wij,x*Vj)    Sum((yj-yi)*Wij,x*Vj)  |
	//  | Sum((fj-fi)*Wij,y*Vj) |     | fi,y |   | 1    Sum((xj-xi)*Wij,y*Vj)    Sum((yj-yi)*Wij,y*Vj)  |
	void PPC_FPMDifference_General(ParticlePair* pp, std::vector<Matrix>* FPM_matrix);
	void PPC_FPMDifference_vec(ParticlePair* pp, std::vector<std::vector<part_prec>>* FPM_matrix, std::string param);

	//Difference + Trim
	//  | Sum((fj-fi)*Wij,x*Vj) |  =  | fi,x | * | 1    Sum((xj-xi)*Wij,x*Vj)    Sum((yj-yi)*Wij,x*Vj)  |
	//  | Sum((fj-fi)*Wij,y*Vj) |     | fi,y |   | 1    Sum((xj-xi)*Wij,y*Vj)    Sum((yj-yi)*Wij,y*Vj)  |
	//void PPC_FPMDifTrim_General(ParticlePair* pp, std::vector<Matrix>* FPM_matrix);
	//void PPC_FPMDifTrim_vec(ParticlePair* pp, std::vector<std::vector<part_prec>>* FPM_matrix, std::string param);


	//FOR CORRECTIVE SPH
	//  | Sum((fj-fi)*Wij,x*Vj) |  =  | fi,x | * | Sum((xj-xi)*Wij,x*Vj)  Sum((yj-yi)*Wij,x*Vj)  Sum((zj-zi)*Wij,x*Vj)  |
	//  | Sum((fj-fi)*Wij,y*Vj) |  =  | fi,y | * | Sum((xj-xi)*Wij,y*Vj)  Sum((yj-yi)*Wij,y*Vj)  Sum((zj-zi)*Wij,y*Vj)  |
	//  | Sum((fj-fi)*Wij,z*Vj) |  =  | fi,z | * | Sum((xj-xi)*Wij,z*Vj)  Sum((yj-yi)*Wij,z*Vj)  Sum((zj-zi)*Wij,z*Vj)  |
	void PPC_CorSPH_General(ParticlePair* pp, std::vector<Matrix>* FPM_matrix);
	void PPC_CorSPH_vec(ParticlePair* pp, std::vector<std::vector<part_prec>>* FPM_matrix, std::string param);



	void DensityRecalculation();
	void VelocityRecalculation();
	void DensityAndVelocityRecalculation();
	void RecalcPrep();


	void DensityDiffusiveTerm(std::vector<part_prec>* drhodt, ParticlePair* pp);


	void RenormFactorCalc(std::vector<part_prec>* gamma, ParticlePair* pp);
	void RenormFactorDerivCalc(std::vector<part_prec_3>* gamma_deriv, ParticlePair* pp);
	void RenormPressurePartCalculation(std::vector<part_prec_3>* dvdt, ParticlePair* pp, std::vector<part_prec>* gamma);
	
	void BoundaryCalculation();
	



	~SPH_CD();


	//COLORING
	void Coloring();
	void ColoringBtType();
	bool ColoringCondition(Particle& p);
	glm::vec3 ColoringGradient(float maxVal, float minVal, float currentVal);

	void punctualColorChange(int number, glm::vec3 color);

	void PRB_refresh();


	int getNrOfParticles(PARTICLETYPE type);
	float getGlobalStats(int number);
	std::string getLocalStats(int number);
};
