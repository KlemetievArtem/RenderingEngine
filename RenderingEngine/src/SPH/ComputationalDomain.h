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


#include "ParticlePair.h"
#include "BoundaryMentor.h"
#include "UniformTreeGrid.h"
#include "EquationClasses.h"
#include "ParticleSource.h"


// Передавать одни и теже опции в SPH и в пары
struct SPH_OPTIONS {
	unsigned int nrOfParticles[ALLTYPES];
	part_prec smoothingKernelLengthCoefficient = 2.0;

	NEIGBOURS_SEARCH_ALGORITHM NBSAlg = UNIFORM_GRID; // DIRECT	UNIFORM_GRID
	//DIMENSIONS KernelDimension = D1;

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

	std::vector<Equation*> ContinuityEquations;
	std::vector<TimeEquation*> DensityUpdates;

	std::vector<Equation*> MomentumConservation_PressureParts;
	std::vector<Equation*> MomentumConservation_ViscosityParts;
	std::vector<TimeEquation*> VelocityUpdates;

	std::vector<TimeEquation*> PositionUpdates;

		Equation* RenormalizationEquation;
	Equation* RenormalizationPressuerPartEquation;

	//Equation* BoundaryDensity;
	//Equation* CorrectionEquation;
	

};
#include "ConcreteEquationClasses.h"

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
private:
	void EquationsInitialization();
	void RealParticlesInitialization(int nrOfParticlesForVolume, glm::vec3 positionMin, glm::vec3 domainSize, glm::vec3 velocity, part_prec smR, part_prec density);
	void BoundaryParticlesInitialization(std::vector<BoundaryBase*>* activeBoundaries, part_prec smR, part_prec density);//(glm::vec3 positionMin, glm::vec3 domainSize, glm::vec3 velocity, part_prec smR, part_prec density);

	void Initilization(glm::vec3 velocity, std::vector<BoundaryBase*>* activeBoundaries);
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
	void RenormPressurePartCalculation(std::vector<part_prec_3>* dvdt, ParticlePair* pp);
	





	std::vector<ParticleSourcesContainer> particleSourceContainers;

	void BoundariesUpdate(part_prec smR, part_prec density);


	void FirstMiddleLastCMoutput();

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
