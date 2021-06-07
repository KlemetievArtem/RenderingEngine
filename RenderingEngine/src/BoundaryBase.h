#pragma once
#include "libs.h"

enum BCPARAMETER {
	VELOCITY_X,
	VELOCITY_Y,
	VELOCITY_Z,
	ALL_PARAMETRS
};
class BoundaryBase;

class BoundaryCondition {
	cd_prec BC_val;
public:
	friend class BoundaryBase;
	BoundaryCondition(cd_prec val): BC_val(val) { }
	BCPARAMETER BC_type;
	BCPARAMETER getBCtype() { return this->BC_type; }
	cd_prec getBCval() { return this->BC_val; }

	//virtual void BoundariesUpdate() = 0;
	virtual Mesh* ReturnMesh() = 0;

//Periodic
	bool b_periodic = false;
//Source
	std::vector<glm::vec3> part_pos;
	bool b_source = false;
	bool alreadyCreatedForTimeInterval = false;


protected:
	virtual void addPosition(glm::vec3 pp) = 0;
	virtual glm::vec3 CalcVelocity() = 0;
	virtual cd_prec getSourceTimeStep() = 0;
	cd_prec timer = 0.0;
	cd_prec timeStep = 0.0; //[ñ]
public:
	virtual ~BoundaryCondition(){}
};



class FirstOrderBC: public BoundaryCondition {
public:
	FirstOrderBC(cd_prec val) :BoundaryCondition(val) { }
	void changeType(BCPARAMETER t) {
		BC_type = t;
	}
	void addPosition(glm::vec3 pp) override {}
	Mesh* ReturnMesh() override { return nullptr; }
	glm::vec3 CalcVelocity() override { assert("FirstOrderBC::CalcVelocity()" && 0); return glm::vec3(0.0); }
	cd_prec getSourceTimeStep() override { return 0.0; }
};

class Velcity_X_BC : public FirstOrderBC {
	//cd_prec BC_val; //[ì/ñ]
public:
	Velcity_X_BC(cd_prec val) : FirstOrderBC(val) {
		changeType(BCPARAMETER::VELOCITY_X);
	}
};
class Velcity_Y_BC : public FirstOrderBC {
	//cd_prec BC_val; //[ì/ñ]
public:
	Velcity_Y_BC(cd_prec val) : FirstOrderBC(val) {
		changeType(BCPARAMETER::VELOCITY_Y);
	}

};
class Velcity_Z_BC : public FirstOrderBC {
	//cd_prec BC_val; //[ì/ñ]
public:
	Velcity_Z_BC(cd_prec val) : FirstOrderBC(val) {
		changeType(BCPARAMETER::VELOCITY_Z);
	}

};


class PeriodicBC: public BoundaryCondition {
public:
	Mesh* optional_mesh = nullptr;
	void addPosition(glm::vec3 pp) override {}
	Mesh* ReturnMesh() override {
		return optional_mesh;
	}
	glm::vec3 CalcVelocity() override { assert("PeriodicBC::CalcVelocity()"&&0); return glm::vec3(0.0); }
	cd_prec getSourceTimeStep() override { return 0.0; }
	PeriodicBC(Mesh* mesh_to_loop) :BoundaryCondition(0.0) {
		b_periodic = true;
		optional_mesh = new Mesh(mesh_to_loop);
	}
};

class SourceBC : public BoundaryCondition {
	Mesh* optional_mesh = nullptr;

	//cd_prec BC_val; //[êã/ñ]
	cd_prec particlesMass;//[êã]
	glm::vec3 normal;
	cd_prec density; //[êã/ñ]
	cd_prec ave_velocity;
public:
	void addPosition(glm::vec3 pp) override { part_pos.push_back(pp); }
	Mesh* ReturnMesh() override { return optional_mesh; }

	SourceBC(cd_prec val, cd_prec mass, cd_prec rho, glm::vec3 source_normal, Mesh* mesh) : BoundaryCondition(val), particlesMass(mass) {
		b_source = true;
		optional_mesh = new Mesh(mesh);
		density = rho;
		timeStep = mass / val;
		ave_velocity = val / (density * optional_mesh->Area());
		normal = source_normal;
	}
	glm::vec3 CalcVelocity() { return static_cast<float>(ave_velocity)*normal; }
	cd_prec getSourceTimeStep() override { return timeStep; }

};




class BoundaryBase {
	Mesh* B_mesh;
	std::vector<BoundaryCondition*> BCarray;
public:
	BoundaryBase(Mesh* mesh) :B_mesh(new Mesh(mesh)) { }
	void addBC(BoundaryCondition* BC) {
		BCarray.push_back(BC);
	}
	void addWallBC(glm::vec3 vel) {
		BCarray.push_back(new Velcity_X_BC(vel.x));
		BCarray.push_back(new Velcity_Y_BC(vel.y));
		BCarray.push_back(new Velcity_Z_BC(vel.z));
	}
	glm::vec3 getVelocity() {
		glm::vec3 retVal(0.0, 0.0, 0.0);
		if (BCarray.size() < 3) {
			if ((!this->isPeriodic()) and (!this->isSource())) assert("BoundaryBase::getVelocity() not enough assigned conditions for velocity" && 0);
			if (this->isSource()) return BCarray[0]->CalcVelocity();
		}
		else {
			for (auto bc : BCarray) {
				switch (bc->getBCtype()) {
				case(BCPARAMETER::VELOCITY_X):
					retVal.x = bc->getBCval();
					break;
				case(BCPARAMETER::VELOCITY_Y):
					retVal.y = bc->getBCval();
					break;
				case(BCPARAMETER::VELOCITY_Z):
					retVal.z = bc->getBCval();
					break;
				default:
					break;
				}
			}
		}
		return retVal;
	}
	cd_prec getTimeStep() {
		cd_prec retVal(0.0);
		if (this->isSource()) retVal = BCarray[0]->getSourceTimeStep();
		return retVal;
	}
	void setPositions(glm::vec3 pp) {
		BCarray[0]->addPosition(pp);
	}
	std::vector<glm::vec3> getPositions() {
		return BCarray[0]->part_pos;
	}
	Mesh* ReturnMesh() {
		return B_mesh;
	}

	void UpdateBoundaryTimers(cd_prec dt) {
		for (auto bc : BCarray) {
			bc->timer += dt;
		}
		//if (BCarray[0]->b_source) std::cout << "BC_Time:" << BCarray[0]->timer << "\n";

	}

	bool isPeriodic() { return BCarray[0]->b_periodic; }
	Mesh* ReturnMeshToLoop() { return BCarray[0]->ReturnMesh(); }

	bool Active() {
		if (BCarray[0]->timer > 0) {
			Deactivate();
			return true;
		}
		else {
			return false;
		}
	}
	void Deactivate() {
		BCarray[0]->timer -= BCarray[0]->timeStep;
	}

	bool isSource() {
		return BCarray[0]->b_source;
	}


	cd_prec getMinX() { return B_mesh->getMinX(); }
	cd_prec getMaxX() { return B_mesh->getMaxX(); }
	cd_prec getMinY() { return B_mesh->getMinY(); }
	cd_prec getMaxY() { return B_mesh->getMaxY(); }
	cd_prec getMinZ() { return B_mesh->getMinZ(); }
	cd_prec getMaxZ() { return B_mesh->getMaxZ(); }

};