#pragma once

enum BOUNDARYCONDITION {
	FIRSTORDER_BC,
	SECONDORDER_BC
};
enum BCPARAMETER {
	VELOCITY_X,
	VELOCITY_Y,
	VELOCITY_Z,
	ALL_PARAMETRS
};

class CD_Boundary {
private:
	Mesh* B_mesh;
	std::vector<BOUNDARYCONDITION> BC_types;
	std::vector<BCPARAMETER> BC_parameters;
	std::vector<float> BC_value;
	bool periodic = false;
	Mesh* periodic_mesh = nullptr;

	bool source = false;

public:
	   
	CD_Boundary(Mesh* mesh, BOUNDARYCONDITION BC_type, BCPARAMETER BC_param, float BC_val)
	: B_mesh(new Mesh(mesh)) {
		periodic = false;
		BC_parameters.resize(ALL_PARAMETRS);
		BC_types.resize(ALL_PARAMETRS);
		BC_value.resize(ALL_PARAMETRS);

		for (int i = 0;i < ALL_PARAMETRS;i++) {
			BC_parameters[i] = static_cast<BCPARAMETER>(i);
			BC_types[i] = FIRSTORDER_BC;
			BC_value[i] = 0.f;
		}

		BC_parameters[BC_param] = BC_param;
		BC_types[BC_param] = BC_type;
		BC_value[BC_param] = BC_val;
	}

	void addBC(BOUNDARYCONDITION BC_type, BCPARAMETER BC_param, float BC_val) {
		if (periodic = true) {
			assert("ÑD_Boundary::addBC BOUNDARY IS PERIODIC" && 0);
		}
		else {
			BC_parameters[BC_param] = BC_param;
			BC_types[BC_param] = BC_type;
			BC_value[BC_param] = BC_val;
		}
	}

	CD_Boundary(Mesh* mesh, Mesh* mesh_to_loop) : B_mesh(new Mesh(mesh)), periodic_mesh(new Mesh(mesh_to_loop)) {
		periodic = true;
	}

	void transformToSource() {
		this->source = true;
		for (int i = 0;i < ALL_PARAMETRS;i++) {
			if (BC_types[i] != FIRSTORDER_BC) {
				BC_types[i] = FIRSTORDER_BC;
				std::cout << "BC_types were changed? you better chek \n";
			}
		}
		if (sqrt(pow(BC_value[BCPARAMETER::VELOCITY_X], 2) + pow(BC_value[BCPARAMETER::VELOCITY_Y], 2) + pow(BC_value[BCPARAMETER::VELOCITY_Z], 2)) == 0) {
			assert("Source boundary velocity is equal to" && 0);
		}
	}


	CD_Boundary(CD_Boundary* newBoundary)
		: B_mesh(newBoundary->B_mesh) {

		for (auto i : newBoundary->BC_types) {
			this->BC_types.push_back(i);
		}
		for (auto i : newBoundary->BC_parameters) {
			this->BC_parameters.push_back(i);
		}
		for (auto i : newBoundary->BC_value) {
			this->BC_value.push_back(i);
		}
		std::cout << this->periodic << "\n";
		newBoundary->periodic = this->periodic;
		if (this->periodic == true) {
			this->periodic_mesh = nullptr;
		}
		else {
			periodic_mesh = newBoundary->periodic_mesh;
		}
	}

	Mesh* ReturnMesh() {
		return B_mesh;
	}
	Mesh* ReturnMeshToLoop() {
		return periodic_mesh;
	}

	std::vector<BOUNDARYCONDITION> getBCtypes() { return this->BC_types; }
	std::vector<BCPARAMETER> getBCparams() { return this->BC_parameters; }
	std::vector<float> getBCvals() { return this->BC_value; }
	bool isPeriodic() { return this->periodic; }
	bool isSource() { return this->source; }



	friend bool operator==(const CD_Boundary& CDB1, const CD_Boundary& CDB2);
	friend bool operator!=(const CD_Boundary& CDB1, const CD_Boundary& CDB2);






	~CD_Boundary() {
		safeDelete(&B_mesh);
		if (periodic_mesh != nullptr) {
			safeDelete(&periodic_mesh);
		}

	}
};
