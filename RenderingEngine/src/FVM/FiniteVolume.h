#pragma once

typedef float fv_prec;
typedef glm::vec3 fv_prec_3;


template <class T>
struct Ch_ {
public:
	Ch_(T v = static_cast<T>(0)) : val(v), dval(static_cast<T>(0.0)) {}
	T val;
	T dval;

	friend Ch_<T> operator+(const Ch_<T>& v1, Ch_<T>& v2) {
		Ch_<T> retVal(v1.val); retVal.dval = v1.dval;
		retVal.val = retVal.val + v2.val;
		retVal.dval = retVal.dval + v2.dval;
		return retVal;
	}
};

struct Line {
	fv_prec_3 m_p1;
	fv_prec_3 m_p2;
	Line() {}
	Line(fv_prec_3 p1, fv_prec_3 p2): m_p1(p1), m_p2(p2){}
	fv_prec_3* getFirstPointPtr() { return &m_p1; }
	fv_prec_3* getSecondPointPtr() { return &m_p2; }
	fv_prec_3 getCentre() {
		return { ((m_p1.x + m_p2.x) / 2.0), ((m_p1.y + m_p2.y) / 2.0), ((m_p1.z + m_p2.z) / 2.0) };
	}
	bool ConnectedWith(Line& l) {
		if (((m_p1 == l.m_p1) or (m_p1 == l.m_p2)) or
			((m_p2 == l.m_p1) or (m_p2 == l.m_p2))) {
			return true;
		}
		else return false;
	}


	inline fv_prec loc_mod(fv_prec_3 vec) {
		return (sqrt(pow(vec.x, 2) + pow(vec.y, 2) + pow(vec.z, 2)));
	}

	friend bool operator==(const Line& l1, const Line& l2) {
		fv_prec Err = 1E-05;
		fv_prec diff1 = sqrt(pow((l1.m_p1 - l2.m_p1).x, 2) + pow((l1.m_p1 - l2.m_p1).y, 2) + pow((l1.m_p1 - l2.m_p1).z, 2));
		fv_prec diff2 = sqrt(pow((l1.m_p2 - l2.m_p2).x, 2) + pow((l1.m_p2 - l2.m_p2).y, 2) + pow((l1.m_p2 - l2.m_p2).z, 2));
		fv_prec diff3 = sqrt(pow((l1.m_p1 - l2.m_p2).x, 2) + pow((l1.m_p1 - l2.m_p2).y, 2) + pow((l1.m_p1 - l2.m_p2).z, 2));
		fv_prec diff4 = sqrt(pow((l1.m_p2 - l2.m_p1).x, 2) + pow((l1.m_p2 - l2.m_p1).y, 2) + pow((l1.m_p1 - l2.m_p1).z, 2));
		if (((diff1 < Err) and (diff2 < Err)) or
			((diff3 < Err) and (diff4 < Err))) {
			return true;
		}
		else return false;
	}

	fv_prec size() {
		fv_prec_3 diff = m_p1 - m_p2;
		return (sqrt(pow(diff.x, 2) + pow(diff.y, 2) + pow(diff.z, 2)));
	}

};

struct Surface {
	std::vector<Line> m_lineArray;
	fv_prec_3 m_Centre;
	Surface() {}
	Surface(Line l1, Line l2, Line l3) {
		m_lineArray.push_back(l1);
		m_lineArray.push_back(l2);
		m_lineArray.push_back(l3);
	}
	fv_prec_3 getCentre() {
		return m_Centre;
	}
	void calcCentre() {
		fv_prec_3 SUM = {0.0, 0.0, 0.0};
		for (auto l : m_lineArray) {
			SUM = SUM + l.getCentre();
		}
		m_Centre = 1.f/ m_lineArray.size() * SUM;
	}

	void addLine(Line l) {
		m_lineArray.push_back(l);
	}

	bool chekIfClosed() {

		return true;
	}

	bool ConnectedWith(Surface& s) {
		for (int i = 0;i < this->m_lineArray.size();i++) {
			for (int ii = 0;ii < s.m_lineArray.size();ii++) {
				//std::cout <<"{"<< this->m_lineArray[i].m_p1.x << "," <<this->m_lineArray[i].m_p1.y  << "}, {" << this->m_lineArray[i].m_p2.x <<","<< this->m_lineArray[i].m_p2.y << "} and {" << s.m_lineArray[ii].m_p1.x <<"," << s.m_lineArray[ii].m_p1.y << "}, {" << s.m_lineArray[ii].m_p2.x<<"," << s.m_lineArray[ii].m_p2.y << "} " << (this->m_lineArray[i] == s.m_lineArray[ii]) <<"\n";
				if (this->m_lineArray[i] == s.m_lineArray[ii]) return true;
			}
		} 
		return false;
	}

	friend bool operator==(const Surface& s1, const Surface& s2) {
		if (s1.m_lineArray.size() != s2.m_lineArray.size()) return false;
		int a_size = s1.m_lineArray.size();
		int nrOfSameLines = 0;
		for (int i = 0;i < a_size;i++) {
			for (int ii = 0;ii < a_size;ii++) {
				if (s1.m_lineArray[i] == s2.m_lineArray[ii]) {
					nrOfSameLines++;
					break;
				}
			}
		}
		if (nrOfSameLines == a_size) return true;
		return false;
	}


	fv_prec size() {
		fv_prec Sum = 0.0;
		std::vector<fv_prec_3*> uniq_points;
		bool no_match_1;
		bool no_match_2;
		// Creating array with all points
		for (auto l : m_lineArray) {
			no_match_1 = true;
			no_match_2 = true;
			if (uniq_points.size() == 0) {
				uniq_points.push_back(&l.m_p1);
				uniq_points.push_back(&l.m_p2);
			}
			else {
				for (auto up : uniq_points) {
					if (l.m_p1 == *up) {
						no_match_1 = false;
					}
					if (l.m_p2 == *up) {
						no_match_2 = false;
					}
					if ((no_match_1 == false) and (no_match_1 == false)) {
						break;
					}
				}
				if (no_match_1) { uniq_points.push_back(&l.m_p1); }
				if (no_match_2) { uniq_points.push_back(&l.m_p2); }
			}
		}
		for (int p = 1;p < uniq_points.size()-1;p++) {
			Sum += SurfaceArea(*(uniq_points[0]), *(uniq_points[p]), *(uniq_points[p+1]));
		}
		return Sum;
	}

};

struct Volume {
	std::vector<Surface> m_surfaceArray;
	fv_prec_3 m_Centre;
	Volume(){}
	Volume(Surface s1, Surface s2, Surface s3, Surface s4) {
		m_surfaceArray.push_back(s1);
		m_surfaceArray.push_back(s2);
		m_surfaceArray.push_back(s3);
		m_surfaceArray.push_back(s4);
	}
	fv_prec_3 getCentre() {
		return m_Centre;
	}
	void calcCentre() {
		fv_prec_3 SUM = { 0.0, 0.0, 0.0 };
		for (auto s : m_surfaceArray) {
			SUM = SUM + s.m_Centre;
		}
		m_Centre = 1.f / m_surfaceArray.size() * SUM;
	}

	void addSurface(Surface s) {
		m_surfaceArray.push_back(s);
	}

	bool ConnectedWith(Volume& v) {
		int a_size = this->m_surfaceArray.size();
		for (int i = 0;i < a_size;i++) {
			for (int ii = 0;ii < v.m_surfaceArray.size();ii++) {
				if (this->m_surfaceArray[i] == v.m_surfaceArray[ii]) {
					return true;
				}
			}
		}
		return false;
	}

};

enum FINITE_VOLUME_DIM {
	FVA_L,
	FVA_S,
	FVA_V
};

enum FINITE_VOLUME_TYPE {
	FVT_DEFAULT,
	FVT_BOUNDARY,
	FVT_ALL
};

class FiniteVolume {
public:
	int m_id;
	FINITE_VOLUME_TYPE m_type;
	fv_prec_3 m_Centre;

	FINITE_VOLUME_DIM FVdim;
	Volume m_v;
	Surface m_s;
	Line m_l;
	
	Ch_<fv_prec_3> m_velocity;  //part_prec_3

	fv_prec absVelocity;           //part_prec
	//int nrOfNeighbours = 0;

	glm::vec3 m_color{ 0.9f, 0.9f, 0.9f };

	//float correctiveKernel;
	Ch_<fv_prec> m_density;
	Ch_<fv_prec> m_pressure;       //part_prec

	fv_prec m_SoundVelocity;
	fv_prec m_Temperature = 357.0;
	fv_prec m_DVisc;
	void setDVisc() { m_DVisc = 1.0E-03; }

	fv_prec m_gamma;
	fv_prec_3 m_grad_gamma;


	fv_prec m_SurfaceTension;
	void setSurfTen() { m_SurfaceTension = 60E-03; }
	fv_prec_3 m_normalToSurface;

	fv_prec dummyParameter;
public:
	FiniteVolume(int id, FINITE_VOLUME_TYPE type, Volume v, glm::vec3 vel, fv_prec pressure)
		: m_id(id), m_type(type), m_v(v) {
		FVdim = FVA_V;
		m_v.calcCentre();
		m_Centre = m_v.getCentre();
		m_velocity.val = vel; m_velocity.dval = fv_prec_3(0.0);
		m_pressure.val = pressure; m_pressure.dval = fv_prec(0.0);
		m_density.val = 1000.0;
		//if (m_type == PARTICLETYPE::REAL){
		//} else { p_gas(); }
		setDVisc();
		setSurfTen();
	}
	FiniteVolume(int id, FINITE_VOLUME_TYPE type, Surface s, glm::vec3 vel, fv_prec pressure)
		: m_id(id), m_type(type), m_s(s) {
		FVdim = FVA_S;
		m_s.calcCentre();
		m_Centre = m_s.getCentre();
		m_velocity.val = vel; m_velocity.dval = fv_prec_3(0.0);
		m_pressure.val = pressure; m_pressure.dval = fv_prec(0.0);
		m_density.val = 1000.0;
		//if (m_type == PARTICLETYPE::REAL){
		//} else { p_gas(); }
		setDVisc();
		setSurfTen();
	}
	FiniteVolume(int id, FINITE_VOLUME_TYPE type, Line l, glm::vec3 vel, fv_prec pressure)
		: m_id(id), m_type(type), m_l(l) {
		FVdim = FVA_L;
		m_Centre = m_l.getCentre();
		m_velocity.val = vel; m_velocity.dval = fv_prec_3(0.0);
		m_pressure.val = pressure; m_pressure.dval = fv_prec(0.0);
		m_density.val = 1000.0;
		//if (m_type == PARTICLETYPE::REAL){
		//} else { p_gas(); }
		setDVisc();
		setSurfTen();
	}


	BoundaryBase* fv_boundary;
	void assignToBoundary(BoundaryBase* cd_boundary) { fv_boundary = cd_boundary; }

	fv_prec distance(FiniteVolume* fv_nb) {
		fv_prec_3 diff = m_Centre - fv_nb->m_Centre;
		return sqrt(pow(diff.x, 2) + pow(diff.y, 2) + pow(diff.z, 2));
	}

	~FiniteVolume() {
		//VirtualCounterpartNormals.resize(0);
		//VirtualCounterpartFlags.resize(0);
	}

};