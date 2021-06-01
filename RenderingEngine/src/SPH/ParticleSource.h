#pragma once


class ParticleSource {
private:
	part_prec_3 m_position;
	part_prec_3 m_velocity;  //part_prec_3

	part_prec m_SmR;
	//glm::vec3 m_color{ 0.9f, 0.9f, 0.9f };
	part_prec m_mass;                //part_prec
	part_prec m_density;

	part_prec m_Temperature = 357.0;
	PARTICLETYPE m_type;
public:
	ParticleSource(part_prec_3 pos, part_prec_3 vel, Particle* p)
		: m_position(pos), m_velocity(vel) {
		m_type = p->m_type;
		m_SmR = p->m_SmR;
		m_mass = p->m_mass;
		m_density = p->m_density.val;
	}
};


class ParticleSourcesContainer {
	std::vector<ParticleSource> particleSources;
	Particle copyParticle;
public:
	ParticleSourcesContainer(Particle p)
		:copyParticle(p) {
		
	}
};