#pragma once


//typedef double part_prec;
//typedef glm::highp_dvec3 part_prec_3;
typedef float part_prec;
typedef glm::vec3 part_prec_3;



enum PARTICLETYPE
{
	REAL,
	BOUNDARY,
	VIRTUAL,
	ALLTYPES
};

template <class T>
struct Ch_ {
public:
	Ch_(T v=static_cast<T>(0)) : val(v),dval(static_cast<T>(0.0)) {}
	T val;
	T dval;
};



struct Particle {
public:
	int m_id;
	PARTICLETYPE m_type;
	Ch_ <part_prec_3> m_position;  //part_prec_3
	Ch_ <part_prec_3> m_velocity;  //part_prec_3

	part_prec absVelocity;           //part_prec
	part_prec m_SmR;                 //part_prec
	int nrOfNeighbours = 0;

	glm::vec3 m_color{ 0.9f, 0.9f, 0.9f };

	//float correctiveKernel;
	part_prec m_mass;                //part_prec
	Ch_<part_prec> m_density;
	part_prec m_pressure;       //part_prec

	part_prec m_SoundVelocity;
	part_prec m_Temperature = 357.0;
	part_prec m_DVisc;


	part_prec gamma;
	part_prec_3 grad_gamma;
	

	part_prec dummyParameter;

	Particle(int id, PARTICLETYPE type, glm::vec3 pos, glm::vec3 vel, part_prec SmR, part_prec density, float mass)
		: m_id(id), m_type(type), m_SmR(SmR), m_mass(mass) {
		m_position.val = pos; m_position.dval = part_prec_3(0.0);
		m_velocity.val = vel; m_velocity.dval = part_prec_3(0.0);
		m_density.val = density; m_density.dval = part_prec(0.0);
		//if (m_type == PARTICLETYPE::REAL)
		p_art_water();
		//else
			//p_gas();
		setDVisc();
	}

	Particle(Particle* part): Particle(part->m_id, part->m_type, part->m_position.val, part->m_velocity.val, part->m_SmR, part->m_density.val, part->m_mass){
		this->m_position.dval = part->m_position.dval;
		this->m_velocity.dval = part->m_velocity.dval;
		this->m_density.dval = part->m_density.dval;
		this->InitPolygonNormal = part->InitPolygonNormal;
		this->m_pressure = part->m_pressure;
		this->m_DVisc = part->m_DVisc;
	}

	~Particle() {
		//VirtualCounterpartNormals.resize(0);
		//VirtualCounterpartFlags.resize(0);
	}

	CD_Boundary* particle_boundary;
	void assignToBoundary(CD_Boundary* cd_boundaty) { particle_boundary = cd_boundaty; }

	part_prec dx(Particle* other) { return this->m_position.val.x - other->m_position.val.x; }
	part_prec dy(Particle* other) { return this->m_position.val.y - other->m_position.val.y; }
	part_prec dz(Particle* other) { return this->m_position.val.z - other->m_position.val.z; }
	part_prec_3 dr(Particle* other) { return this->m_position.val - other->m_position.val; }
	part_prec distance(Particle* other) { return std::sqrt(dot(dr(other), dr(other))); }



	part_prec dVx(Particle* other) { return this->m_velocity.val.x - other->m_velocity.val.x; }
	part_prec dVy(Particle* other) { return this->m_velocity.val.y - other->m_velocity.val.y; }
	part_prec dVz(Particle* other) { return this->m_velocity.val.z - other->m_velocity.val.z; }

	part_prec_3 dV(Particle* other) { return this->m_velocity.val - other->m_velocity.val; }


	part_prec dx(glm::vec3 pos) { return this->m_position.val.x - pos.x; }
	part_prec dy(glm::vec3 pos) { return this->m_position.val.y - pos.y; }
	part_prec dz(glm::vec3 pos) { return this->m_position.val.z - pos.z; }
	part_prec_3 dr(glm::vec3 pos) { return { dx(pos), dy(pos), dz(pos) };}
	part_prec distance(glm::vec3 pos) { return std::sqrt(dot(dr(pos), dr(pos))); }
	

	void p_art_water() {
		int gamma = 7;
		setSoundVel(1480.0);
		part_prec Pressure0 = 1.0E05;
		part_prec b = 10000;

		part_prec Dens0 = 1.0;

		//p = b * ((m_dens / Dens0)**gamma - 1)
		//m_pressure.val = (b*pow((m_density.val / Dens0), gamma) - 1);
		//m_pressure.val = Pressure0 + Dens0 * pow(m_SoundVelocity, 2) / static_cast<part_prec>(gamma) * (pow(m_density.val / Dens0, gamma) - static_cast<part_prec>(1.0));
		//m_pressure.val = pow(m_SoundVelocity, 2)  * (m_density.val - Dens0);

		m_pressure =  b* (pow(m_density.val / Dens0, gamma) - static_cast<part_prec>(1.0));

		//m_pressure.val = b * ((m_density.val- Dens0) / Dens0);

		//if(m_id == 961 - 31){
		//	std::cout << "density = " <<m_density.val << "\n";
		//	std::cout << "pressure = " << Dens0 << "*" << pow(m_SoundVelocity, 2) << "/" << gamma << "*(" << "(" << m_density.val <<"/"<< Dens0 << ")^" << gamma << "-" << 1.f << ")" << "=" << m_pressure.val << "\n";
		//}
	}
	void p_gas() {
		float gamma = 1.4f;
		m_pressure = ((gamma - 1)*m_density.val*m_Temperature);
		setSoundVel(sqrt((gamma - 1)*m_Temperature));

	}

	void setPressureTo(part_prec pressure) {
		this->m_pressure = pressure;
	}

	void addNeighbour() {
		nrOfNeighbours++;
	}
	void refreshNeighbours() {
		nrOfNeighbours = 0;
	}

private:
	void setSoundVel(float val) {
		m_SoundVelocity = val;
	}


	void setDVisc() {
		if (m_type == PARTICLETYPE::REAL)
			m_DVisc = 0;
		if (m_type == PARTICLETYPE::VIRTUAL)
			m_DVisc = 0;
		if (m_type == PARTICLETYPE::BOUNDARY)
			m_DVisc = 1.0E-03;


		m_DVisc = 1000E-06;
	}
public: 
	glm::vec3 InitPolygonNormal;// = glm::vec3(0.f);
	glm::vec3 GetNormal() { return this->InitPolygonNormal; }
};



