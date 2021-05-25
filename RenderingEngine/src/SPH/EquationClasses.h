#pragma once

class K0_ContinuityEquation_dval : public ConcreteEquation {
public:
	K0_ContinuityEquation_dval() {};
	void Func(std::vector<part_prec>* val, ParticlePair * pp, PARTICLEPAIR_NUMERATIN ppn) override {
		switch (ppn)
		{
		case PPN_MAIN:
			(*val)[pp->Mpart_ptr->m_id] += pp->Mpart_ptr->m_density.val*pp->NBpart_ptr->m_mass / pp->NBpart_ptr->m_density.val*glm::dot(pp->Mpart_ptr->dV(pp->NBpart_ptr), pp->dWd(R, pp->Mpart_ptr)*pp->vec_e(pp->Mpart_ptr));
			break;
		case PPN_NEIGHBOUR:
			(*val)[pp->NBpart_ptr->m_id] += pp->NBpart_ptr->m_density.val*pp->Mpart_ptr->m_mass / pp->Mpart_ptr->m_density.val*glm::dot(pp->NBpart_ptr->dV(pp->Mpart_ptr), pp->dWd(R, pp->NBpart_ptr)*pp->vec_e(pp->NBpart_ptr));
			break;
		default:
			break;
		}
	}

	void Func(std::vector<part_prec_3>* vector_val, ParticlePair * pp, DIMENSIONS dim, PARTICLEPAIR_NUMERATIN ppn) override { 
		assert("K0_ContinuityEquation_dval::wrong_func" && 0);
	}
};

class K1_ContinuityEquation_dval : public ConcreteEquation {
public:
	K1_ContinuityEquation_dval() {};
	void Func(std::vector<part_prec>* val, ParticlePair * pp, PARTICLEPAIR_NUMERATIN ppn) override {
		switch (ppn)
		{
		case PPN_MAIN:
			(*val)[pp->Mpart_ptr->m_id] += pp->NBpart_ptr->m_mass * glm::dot(pp->Mpart_ptr->dV(pp->NBpart_ptr), pp->dWd(R, pp->Mpart_ptr)* pp->vec_e(pp->Mpart_ptr));
			break;
		case PPN_NEIGHBOUR:
			(*val)[pp->NBpart_ptr->m_id] += pp->Mpart_ptr->m_mass * glm::dot(pp->NBpart_ptr->dV(pp->Mpart_ptr), pp->dWd(R, pp->NBpart_ptr)*pp->vec_e(pp->NBpart_ptr));
			break;
		default:
			break;
		}
	}

	void Func(std::vector<part_prec_3>* vector_val, ParticlePair * pp, DIMENSIONS dim, PARTICLEPAIR_NUMERATIN ppn) override {
		assert("K1_ContinuityEquation_dval::wrong_func" && 0);
	}
};

class K2_ContinuityEquation_dval : public ConcreteEquation {
public:
	K2_ContinuityEquation_dval() {};
	void Func(std::vector<part_prec>* val, ParticlePair * pp, PARTICLEPAIR_NUMERATIN ppn) override {
		switch (ppn)
		{
		case PPN_MAIN:
			(*val)[pp->Mpart_ptr->m_id] += pp->Mpart_ptr->m_density.val*pp->NBpart_ptr->m_mass / pp->NBpart_ptr->m_density.val*glm::dot(pp->Mpart_ptr->dV(pp->NBpart_ptr), pp->dWd(R, pp->Mpart_ptr)*pp->vec_e(pp->Mpart_ptr));
			break;
		case PPN_NEIGHBOUR:
			(*val)[pp->NBpart_ptr->m_id] += pp->NBpart_ptr->m_density.val*pp->Mpart_ptr->m_mass / pp->Mpart_ptr->m_density.val*glm::dot(pp->NBpart_ptr->dV(pp->Mpart_ptr), pp->dWd(R, pp->NBpart_ptr)*pp->vec_e(pp->NBpart_ptr));
			break;
		default:
			break;
		}
	}

	void Func(std::vector<part_prec_3>* vector_val, ParticlePair * pp, DIMENSIONS dim, PARTICLEPAIR_NUMERATIN ppn) override {
		assert("K2_ContinuityEquation_dval::wrong_func" && 0);
	}
};


class K0_ContinuityEquation_val : public ConcreteEquation {
public:
	K0_ContinuityEquation_val() {};
	void Func(std::vector<part_prec>* val, ParticlePair * pp, PARTICLEPAIR_NUMERATIN ppn) override {
		switch (ppn)
		{
		case PPN_MAIN:
			(*val)[pp->Mpart_ptr->m_id] += pp->NBpart_ptr->m_mass * pp->Mpart_ptr->m_density.val / pp->NBpart_ptr->m_density.val * pp->W(pp->Mpart_ptr);
			break;
		case PPN_NEIGHBOUR:
			(*val)[pp->NBpart_ptr->m_id] += pp->Mpart_ptr->m_mass * pp->NBpart_ptr->m_density.val / pp->Mpart_ptr->m_density.val *pp->W(pp->NBpart_ptr);
			break;
		default:
			break;
		}
	}

	void Func(std::vector<part_prec_3>* vector_val, ParticlePair * pp, DIMENSIONS dim, PARTICLEPAIR_NUMERATIN ppn) override {
		assert("K0_ContinuityEquation_val::wrong_func" && 0);
	}
};

class K1_ContinuityEquation_val : public ConcreteEquation {
public:
	K1_ContinuityEquation_val() {};
	void Func(std::vector<part_prec>* val, ParticlePair * pp, PARTICLEPAIR_NUMERATIN ppn) override {
		switch (ppn)
		{
		case PPN_MAIN:
			(*val)[pp->Mpart_ptr->m_id] += pp->NBpart_ptr->m_mass * pp->W(pp->Mpart_ptr);
			break;
		case PPN_NEIGHBOUR:
			(*val)[pp->NBpart_ptr->m_id] += pp->Mpart_ptr->m_mass * pp->W(pp->NBpart_ptr);
			break;
		default:
			break;
		}
	}

	void Func(std::vector<part_prec_3>* vector_val, ParticlePair * pp, DIMENSIONS dim, PARTICLEPAIR_NUMERATIN ppn) override {
		assert("K1_ContinuityEquation_val::wrong_func" && 0);
	}
};

class K2_ContinuityEquation_val : public ConcreteEquation {
public:
	K2_ContinuityEquation_val() {};
	void Func(std::vector<part_prec>* val, ParticlePair * pp, PARTICLEPAIR_NUMERATIN ppn) override {
		switch (ppn)
		{
		case PPN_MAIN:
			(*val)[pp->Mpart_ptr->m_id] += pp->NBpart_ptr->m_mass * pp->NBpart_ptr->m_density.val / pp->Mpart_ptr->m_density.val * pp->W(pp->Mpart_ptr);
			break;
		case PPN_NEIGHBOUR:
			(*val)[pp->NBpart_ptr->m_id] += pp->Mpart_ptr->m_mass * pp->Mpart_ptr->m_density.val / pp->NBpart_ptr->m_density.val *pp->W(pp->NBpart_ptr);
			break;
		default:
			break;
		}
	}

	void Func(std::vector<part_prec_3>* vector_val, ParticlePair * pp, DIMENSIONS dim, PARTICLEPAIR_NUMERATIN ppn) override {
		assert("K2_ContinuityEquation_val::wrong_func" && 0);
	}

};





class K0_VelPressurePart_dval : public ConcreteEquation {
public:
	K0_VelPressurePart_dval() {};

	void Func(std::vector<part_prec>* val, ParticlePair * pp, PARTICLEPAIR_NUMERATIN ppn) override {
		assert("K0_VelPressurePart_dval::wrong_func" && 0);
	}

	void Func(std::vector<part_prec_3>* vector_val, ParticlePair * pp, DIMENSIONS dim, PARTICLEPAIR_NUMERATIN ppn) override {
		switch (ppn)
		{
		case PPN_MAIN:
			(*vector_val)[pp->Mpart_ptr->m_id] += (-pp->NBpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val*pp->NBpart_ptr->m_density.val)*(pp->Mpart_ptr->m_pressure.val + pp->NBpart_ptr->m_pressure.val)*pp->dWd(R, pp->Mpart_ptr) * pp->vec_e(pp->Mpart_ptr);
			break;
		case PPN_NEIGHBOUR:
			(*vector_val)[pp->NBpart_ptr->m_id] += (-pp->Mpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val*pp->Mpart_ptr->m_density.val)*(pp->NBpart_ptr->m_pressure.val + pp->Mpart_ptr->m_pressure.val)*pp->dWd(R, pp->NBpart_ptr) * pp->vec_e(pp->NBpart_ptr);
			break;
		default:
			break;
		}
		
		
	}
};

class K1_VelPressurePart_dval : public ConcreteEquation {
public:
	K1_VelPressurePart_dval() {};

	void Func(std::vector<part_prec>* val, ParticlePair * pp, PARTICLEPAIR_NUMERATIN ppn) override {
		assert("K1_VelPressurePart_dval::wrong_func" && 0);
	}

	void Func(std::vector<part_prec_3>* vector_val, ParticlePair * pp, DIMENSIONS dim, PARTICLEPAIR_NUMERATIN ppn) override {
		switch (ppn)
		{
		case PPN_MAIN:
			(*vector_val)[pp->Mpart_ptr->m_id] += (-pp->NBpart_ptr->m_mass)*(pp->Mpart_ptr->m_pressure.val / pow(pp->Mpart_ptr->m_density.val, 2) + pp->NBpart_ptr->m_pressure.val / pow(pp->NBpart_ptr->m_density.val, 2))*pp->dWd(R, pp->Mpart_ptr) * pp->vec_e(pp->Mpart_ptr);
			break;
		case PPN_NEIGHBOUR:
			(*vector_val)[pp->NBpart_ptr->m_id] += (-pp->Mpart_ptr->m_mass)*(pp->NBpart_ptr->m_pressure.val / pow(pp->NBpart_ptr->m_density.val, 2) + pp->Mpart_ptr->m_pressure.val / pow(pp->Mpart_ptr->m_density.val, 2))*pp->dWd(R, pp->NBpart_ptr) * pp->vec_e(pp->NBpart_ptr);
			break;
		default:
			break;
		}		
	}
};

class K2_VelPressurePart_dval : public ConcreteEquation {
public:
	K2_VelPressurePart_dval() {};

	void Func(std::vector<part_prec>* val, ParticlePair * pp, PARTICLEPAIR_NUMERATIN ppn) override {
		assert("K2_VelPressurePart_dval::wrong_func" && 0);
	}

	void Func(std::vector<part_prec_3>* vector_val, ParticlePair * pp, DIMENSIONS dim, PARTICLEPAIR_NUMERATIN ppn) override {
		switch (ppn)
		{
		case PPN_MAIN:
			(*vector_val)[pp->Mpart_ptr->m_id] += (-pp->NBpart_ptr->m_mass)*(pp->Mpart_ptr->m_pressure.val*pp->NBpart_ptr->m_density.val / pow(pp->Mpart_ptr->m_density.val, 3) + pp->NBpart_ptr->m_pressure.val*pp->Mpart_ptr->m_density.val / pow(pp->NBpart_ptr->m_density.val, 3))*pp->dWd(R, pp->Mpart_ptr) * pp->vec_e(pp->Mpart_ptr);
			break;
		case PPN_NEIGHBOUR:
			(*vector_val)[pp->NBpart_ptr->m_id] += (-pp->Mpart_ptr->m_mass)*(pp->NBpart_ptr->m_pressure.val*pp->Mpart_ptr->m_density.val / pow(pp->NBpart_ptr->m_density.val, 3) + pp->Mpart_ptr->m_pressure.val*pp->NBpart_ptr->m_density.val / pow(pp->Mpart_ptr->m_density.val, 3))*pp->dWd(R, pp->NBpart_ptr) * pp->vec_e(pp->NBpart_ptr);
			break;
		default:
			break;
		}
	}
};






class K0_VelViscosityPart_dval : public ConcreteEquation {
public:
	K0_VelViscosityPart_dval() {};

	void Func(std::vector<part_prec>* val, ParticlePair * pp, PARTICLEPAIR_NUMERATIN ppn) override {
		assert("K0_VelViscosityPart_dval::wrong_func" && 0);
	}

	void Func(std::vector<part_prec_3>* vector_val, ParticlePair * pp, DIMENSIONS dim, PARTICLEPAIR_NUMERATIN ppn) override {
		switch (ppn)
		{
		case PPN_MAIN:
			(*vector_val)[pp->Mpart_ptr->m_id] += (static_cast<part_prec>(2.0)*pp->Mpart_ptr->m_DVisc)*pp->Mpart_ptr->dV(pp->NBpart_ptr) / (pp->getDist())*pp->dWd(R, pp->Mpart_ptr);
			break;
		case PPN_NEIGHBOUR:
			(*vector_val)[pp->NBpart_ptr->m_id] += (static_cast<part_prec>(2.0)*pp->NBpart_ptr->m_DVisc)*pp->Mpart_ptr->dV(pp->Mpart_ptr) / (pp->getDist())*pp->dWd(R, pp->NBpart_ptr);
			break;
		default:
			break;
		}
	}
};

class K1_VelViscosityPart_dval : public ConcreteEquation {
public:
	K1_VelViscosityPart_dval() {};

	void Func(std::vector<part_prec>* val, ParticlePair * pp, PARTICLEPAIR_NUMERATIN ppn) override {
		assert("K1_VelViscosityPart_dval::wrong_func" && 0);
	}

	void Func(std::vector<part_prec_3>* vector_val, ParticlePair * pp, DIMENSIONS dim, PARTICLEPAIR_NUMERATIN ppn) override {
		switch (ppn)
		{
		case PPN_MAIN:
			(*vector_val)[pp->Mpart_ptr->m_id] += (static_cast<part_prec>(2.0)*pp->Mpart_ptr->m_DVisc)*static_cast<part_prec>(dim + 2)*glm::dot(pp->Mpart_ptr->dV(pp->NBpart_ptr), pp->vec_e(pp->Mpart_ptr)) / (pp->Mpart_ptr->m_density.val*pp->NBpart_ptr->m_density.val*pp->getDist())*pp->dWd(R, pp->Mpart_ptr) * pp->vec_e(pp->Mpart_ptr);
			break;
		case PPN_NEIGHBOUR:
			(*vector_val)[pp->NBpart_ptr->m_id] += (static_cast<part_prec>(2.0)*pp->NBpart_ptr->m_DVisc)*static_cast<part_prec>(dim + 2)*glm::dot(pp->Mpart_ptr->dV(pp->Mpart_ptr), pp->vec_e(pp->NBpart_ptr)) / (pp->NBpart_ptr->m_density.val*pp->Mpart_ptr->m_density.val*pp->getDist())*pp->dWd(R, pp->NBpart_ptr) * pp->vec_e(pp->NBpart_ptr);
			break;
		default:
			break;
		}
	}
};

class K2_VelViscosityPart_dval : public ConcreteEquation {
public:
	K2_VelViscosityPart_dval() {};

	void Func(std::vector<part_prec>* val, ParticlePair * pp, PARTICLEPAIR_NUMERATIN ppn) override {
		assert("K2_VelViscosityPart_dval::wrong_func" && 0);
	}

	void Func(std::vector<part_prec_3>* vector_val, ParticlePair * pp, DIMENSIONS dim, PARTICLEPAIR_NUMERATIN ppn) override {
		switch (ppn)
		{
		case PPN_MAIN:
			(*vector_val)[pp->Mpart_ptr->m_id] += ((pp->Mpart_ptr->m_DVisc + pp->NBpart_ptr->m_DVisc)) / (static_cast<part_prec>(2.0)*pp->Mpart_ptr->m_density.val*pp->NBpart_ptr->m_density.val*pp->getDist())*(static_cast<part_prec>(dim + 2)*glm::dot(pp->Mpart_ptr->dV(pp->NBpart_ptr), pp->vec_e(pp->Mpart_ptr))*pp->vec_e(pp->Mpart_ptr) + pp->Mpart_ptr->dV(pp->NBpart_ptr))*pp->dWd(R, pp->Mpart_ptr);
			break;
		case PPN_NEIGHBOUR:
			(*vector_val)[pp->NBpart_ptr->m_id] += ((pp->NBpart_ptr->m_DVisc + pp->Mpart_ptr->m_DVisc)) / (static_cast<part_prec>(2.0)*pp->NBpart_ptr->m_density.val*pp->Mpart_ptr->m_density.val*pp->getDist())*(static_cast<part_prec>(dim + 2)*glm::dot(pp->NBpart_ptr->dV(pp->Mpart_ptr), pp->vec_e(pp->NBpart_ptr))*pp->vec_e(pp->NBpart_ptr) + pp->NBpart_ptr->dV(pp->Mpart_ptr))*pp->dWd(R, pp->NBpart_ptr);
			break;
		default:
			break;
		}
	}
};



class CorrectionFactor : public ConcreteEquation {
public:
	CorrectionFactor() {};
	void Func(std::vector<part_prec>* val, ParticlePair * pp, PARTICLEPAIR_NUMERATIN ppn) override {
		switch (ppn)
		{
		case PPN_MAIN:
			(*val)[pp->Mpart_ptr->m_id] += pp->NBpart_ptr->m_mass / pp->NBpart_ptr->m_density.val * pp->W(pp->Mpart_ptr);
			break;
		case PPN_NEIGHBOUR:
			(*val)[pp->NBpart_ptr->m_id] += pp->Mpart_ptr->m_mass / pp->Mpart_ptr->m_density.val * pp->W(pp->NBpart_ptr);
			break;
		default:
			break;
		}
	}

	void Func(std::vector<part_prec_3>* vector_val, ParticlePair * pp, DIMENSIONS dim, PARTICLEPAIR_NUMERATIN ppn) override {
		assert("CorrectionFactor::wrong_func" && 0);
	}
};







class DensityUpdate_VAL : public ConcreteTimeUpdateEquation {
public:
	DensityUpdate_VAL() {};

	void FuncT(Particle* p, part_prec t) override {
		p->m_density.val = p->m_density.dval;
	}
	void FuncT(Particle* p, part_prec_3 val) override { }
};

class DensityUpdate_DVAL : public ConcreteTimeUpdateEquation {
public:
	DensityUpdate_DVAL() {};

	void FuncT(Particle* p, part_prec t) override {
		p->m_density.val += p->m_density.dval * t;
	}
	void FuncT(Particle* p, part_prec_3 val) override { }
};



class VelocityUpdate_VAL : public ConcreteTimeUpdateEquation {
public:
	VelocityUpdate_VAL() {};

	void FuncT(Particle* p, part_prec t) override {
		p->m_velocity.val = p->m_velocity.dval;
	}
	void FuncT(Particle* p, part_prec_3 val) override { }

};

class VelocityUpdate_DVAL : public ConcreteTimeUpdateEquation {
public:
	VelocityUpdate_DVAL() {};

	void FuncT(Particle* p, part_prec t) override {
		p->m_velocity.val += p->m_velocity.dval * t;
	}
	void FuncT(Particle* p, part_prec_3 val) override { }
};



class PositionUpdate_add : public ConcreteTimeUpdateEquation {
public:
	PositionUpdate_add() {};

	void FuncT(Particle* p, part_prec t) override { }

	void FuncT(Particle* p, part_prec_3 val) override {
		p->m_position.val += val;
	}
};







class BoundaryDensityUpdate_VAL : public ConcreteEquation {
public:
	BoundaryDensityUpdate_VAL() {};
	void Func(std::vector<part_prec>* val, ParticlePair * pp, PARTICLEPAIR_NUMERATIN ppn) override {
		switch (ppn)
		{
		case PPN_MAIN:
			(*val)[pp->Mpart_ptr->m_id] += pp->NBpart_ptr->m_mass * pp->W(pp->Mpart_ptr);
			break;
		case PPN_NEIGHBOUR:
			(*val)[pp->NBpart_ptr->m_id] += pp->Mpart_ptr->m_mass * pp->W(pp->NBpart_ptr);
			break;
		default:
			break;
		}
	}

	void Func(std::vector<part_prec_3>* vector_val, ParticlePair * pp, DIMENSIONS dim, PARTICLEPAIR_NUMERATIN ppn) override {
		assert("BoundaryDensityUpdate_VAL::wrong_func" && 0);
	}
};



class BoundaryVelocityUpdate_VAL : public ConcreteEquation {
public:
	BoundaryVelocityUpdate_VAL() {};

	void Func(std::vector<part_prec>* val, ParticlePair * pp, PARTICLEPAIR_NUMERATIN ppn) override {
		assert("BoundaryVelocityUpdate_VAL::wrong_func" && 0);
	}

	void Func(std::vector<part_prec_3>* vector_val, ParticlePair * pp, DIMENSIONS dim, PARTICLEPAIR_NUMERATIN ppn) override {
		switch (ppn){
		case PPN_MAIN:
			(*vector_val)[pp->Mpart_ptr->m_id] += pp->NBpart_ptr->m_velocity.val * pp->NBpart_ptr->m_mass / pp->NBpart_ptr->m_density.val * pp->W(pp->Mpart_ptr);
			break;
		case PPN_NEIGHBOUR:
			(*vector_val)[pp->NBpart_ptr->m_id] += pp->Mpart_ptr->m_velocity.val * pp->Mpart_ptr->m_mass / pp->Mpart_ptr->m_density.val * pp->W(pp->NBpart_ptr);
			break;
		default:
			break;
		}
	}
};



class BoundaryCorrection : public ConcreteEquation {
public:
	BoundaryCorrection() {};
	void Func(std::vector<part_prec>* val, ParticlePair * pp, PARTICLEPAIR_NUMERATIN ppn) override {
		switch (ppn)
		{
		case PPN_MAIN:
			(*val)[pp->Mpart_ptr->m_id] += pp->NBpart_ptr->m_mass / pp->NBpart_ptr->m_density.val * pp->W(pp->Mpart_ptr);
			break;
		case PPN_NEIGHBOUR:
			(*val)[pp->NBpart_ptr->m_id] += pp->Mpart_ptr->m_mass / pp->Mpart_ptr->m_density.val * pp->W(pp->NBpart_ptr);
			break;
		default:
			break;
		}
	}

	void Func(std::vector<part_prec_3>* vector_val, ParticlePair * pp, DIMENSIONS dim, PARTICLEPAIR_NUMERATIN ppn) override {
		assert("BoundaryCorrection::wrong_func" && 0);
	}
};


