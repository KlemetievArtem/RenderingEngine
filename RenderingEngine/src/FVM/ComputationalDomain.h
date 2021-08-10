#pragma once

#include "virtualComputationalDomain.h"
#include "FiniteVolume.h"

enum PHI_PARAMS {
	UX,
	UY,
	UZ,
	P,
	T //ЧЕ
};

#include "FVMesh.h"

enum MESH_STRUCTURE {
	STRUCTED,
	UNSTRUCTED
};
enum TIME_INTEGRATION_SCHEME {
	EXPLICIT,
	IMPLICIT,
	SEMI_IMPICIT,
	NONE
};
enum CONV_DIFF_SCHEME {
	CENTRAL_DIFFERENSING,
	UPWIND,
	HYBRID,
	POWER_LAW,
	EXPONENTIAL
};


struct RendererBuffer {
	fv_prec_3 position;
	glm::vec3 color;
	fv_prec_3 size;
	fv_prec_3 rotation = {0.0, 0.0, 0.0};
};

// Передавать одни и теже опции в SPH и в пары
struct FVM_OPTIONS {

	MESH_STRUCTURE mesh_s = STRUCTED;
	int nrOfFV[FVT_ALL];
	fv_prec_3 nrOfFVinDir = {1,1,1};
	//int densityChangeAlgorithm = 2;
	//int velocityChangeAlgorithm_pressurePart = 2;
	//int velocityChangeAlgorithm_viscosityPart = 2;
	CONV_DIFF_SCHEME convectionDiffusionScheme = POWER_LAW;
	TIME_INTEGRATION_SCHEME timeIntegrationScheme = EXPLICIT; // EXPLICIT IMPLICIT SEMI_IMPICIT SECOND_ORDER_SCEME NONE
	// true  false
	bool firstCycle = true;
	bool cornerVP = true;
	bool IsStab = false;
	bool Density_Diffusion_term = false;
	bool BoundaryCalculation = false;
	bool CalcGamma = true;

	bool PressurePartCalc = true;
	bool ViscosityPartCalc = true;
	bool SurfaceTensionPartCalc = false;

	bool TimeStepSizing = true;

	bool Steady = true;

	fv_prec_3 average_dim_steps;
};

struct FVM_ESENTIALS {
	std::vector<FiniteVolume*> FiniteVolumes;
	std::vector<FVMeshNodeEmpty*> FiniteVolumeMesh;
	~FVM_ESENTIALS() {
		for (auto*&i : FiniteVolumes)
			delete i;
		for (auto*&i : FiniteVolumeMesh)
			delete i;
	}
};


class FVM_CD : public CompDomain {
private:
	FVM_ESENTIALS FVM;
	FVM_OPTIONS m_options;

	bool FirstSycle = true;

	//part_prec BH_D;
	//part_prec BH_r;
	//BoundaryMentor BM;
	//std::unique_ptr<BoundaryMentor> BM;
	std::vector<RendererBuffer> PRB;

	
public:
	FVM_CD(FVM_OPTIONS options, DIMENSIONS dim) {
		setTypeTo(MICROCOSME::MC_SPH);
		nrOfDim = dim;
		m_options = options;
	}
private:
	void Initilization(glm::vec3 velocity, std::vector<BoundaryBase*>* activeBoundaries);
	//void InitialRendering(std::vector<Mesh*>* meshes);
	void MeshInitilization(glm::vec3 positionMin, glm::vec3 positionMax);
	void BoundaryMeshInitialization(std::vector<BoundaryBase*>* activeBoundaries);
	void PHIEquationInitializatin();

	//void UpdateRendering(std::vector<Model*>* models);
	void UpdateRendering(std::vector<Model*>* models, Texture* tex, Texture* tex_specualar, std::vector<Material*>* materials);
	void AfterRendering(std::vector<Model*>* models);

	void timeStep(cd_prec dt);
	void timeStep_thread(cd_prec dt, std::atomic<bool>& dataReadyForRender, std::atomic<bool>& dataIsRendering);
	void timeUpdate(cd_prec dt);
	void SIMPLE_ALGO(cd_prec dt);

	inline fv_prec PeFunc(fv_prec Pe);



	Matrix Solve(std::vector<FVMeshNode*>* Mesh, std::vector<fv_prec>* phi_param);
	void SolveWithFirstZero(std::vector<FVMeshNode*>* Mesh, std::vector<fv_prec>* phi_param);

	void SaveMaxVelocity();


	void FirstMiddleLastCMoutput();

	~FVM_CD();


	//COLORING
	void Coloring();

	glm::vec3 ColoringGradient(float maxVal, float minVal, float currentVal);

	void PRB_refresh();

	void punctualColorChange(int number, glm::vec3 color);

	float getGlobalStats(int number);
	std::string getLocalStats(int number);
};
