#pragma once
#include <math.h>
#include <unordered_set>
#include <functional>
#include <vector>
#include <unordered_map>
#include "geos/geos.h"
#include <geos/index/kdtree/KdTree.h>

#define PI acos(-1)
#define SectEps 1.0e-4
#define SectHashEps 1.0e4
#define DllEXPORT __declspec(dllexport)

using namespace std;
using namespace geos::geom;
using namespace geos::index::kdtree;

enum class facie
{
	_default = 0,
	channel_sand = 1,
	levee = 2,
	crevasse = 3,
	migrate_mud = 4,
};

class DllEXPORT rect3d {
private:
	double xmin;
	double xmax;
	double ymin;
	double ymax;
	double zmin;
	double zmax;
public:
	rect3d(const double& xmin, const double& xmax, const double& ymin, const double& ymax, const double& zmin, const double& zmax);
	rect3d(const rect3d& other);
	double get_xmin() const;
	double get_xmax() const;
	double get_ymin() const;
	double get_ymax() const;
	double get_zmin() const;
	double get_zmax() const;
};

class DllEXPORT well3d {
private:
	facie f;
	int wellid;
	Coordinate pos;
	double top;
	double bot;
public:
	well3d(const int& wellid, const  Coordinate& pos, const facie& facie, const  double& top, const  double& bot);

	well3d(const well3d& other);

	const Coordinate& get_pos() const;
	double get_top() const;
	double get_bot() const;
	double get_thick() const;
	int get_id() const;
	facie get_facie() const;
};

class DllEXPORT well3d_sequence {
private:
	unique_ptr<vector<well3d>> inwells;
	int mid;
	unique_ptr<KdTree> tree;

	well3d_sequence(const well3d_sequence& other) = delete;

	well3d_sequence& operator=(const well3d_sequence&) = delete;
public:
	well3d_sequence(const vector<well3d>* wells);

	size_t get_size() const;

	const well3d& get_at(size_t index) const;

	const well3d& get_with_id(const int& wellid) const;

	bool is_empty() const;

	unique_ptr<vector<well3d>> query(const Envelope& range) const;

	unique_ptr<vector<well3d>> query(const geos::geom::Polygon* range) const;
};

class DllEXPORT well2d {
private:
	facie f;
	int wellid;
	Coordinate pos;
public:
	well2d(const int& wellid, const  Coordinate& pos, const facie& facie);

	well2d(const well2d& other);

	const Coordinate& get_pos() const;

	int get_id() const;

	facie get_facie() const;
};

class DllEXPORT well2d_sequence {
private:
	unique_ptr<vector<well2d>> wells;
	int maxid;
	int minid;
	unique_ptr<KdTree> tree;

	well2d_sequence(const well2d_sequence& other) = delete;

	well2d_sequence& operator=(const well2d_sequence&) = delete;
public:
	well2d_sequence(vector<well3d>* wells);

	size_t get_size() const;

	const well2d& get_at(size_t index) const;

	const well2d& get_with_id(const int& wellid) const;

	bool is_empty() const;

	unique_ptr<vector<well2d>> query(const Envelope& range) const;

	unique_ptr<vector<well2d>> query(const geos::geom::Polygon* range) const;
};

class DllEXPORT modelpar {
private:
	double xy_tol;
	double z_tol;
	double ctrl_dis;
	int undefinevalue;
public:
	modelpar(const double& xy_tol, const double& z_tol, const double& ctrl_dis, const int& undefine);

	modelpar(const modelpar& other);

	double get_xytol() const;

	double get_ztol() const;

	double get_ctrl_dis() const;

	int get_undefinevalue() const { return undefinevalue; }

	void set_ctrl_dis(const double& dis);
};

class DllEXPORT gridpoint {
private:
	Coordinate p;
	void* data;
public:
	gridpoint(double p_x, double p_y, double p_z, void* p_data);
	gridpoint(const Coordinate& p_p, void* p_data);

	double getX() const { return p.x; }
	double getY() const { return p.y; }
	double getZ() const { return p.z; }
	const Coordinate& getCoordinate() { return p; }
	void* getData() { return data; }
};

class DllEXPORT grid3d {
private:
	unique_ptr<modelpar> pt_par;
	//unique_ptr<vector<gridpoint>> pt_points;
	unique_ptr<int[]> facies;
	rect3d inrange;
	Coordinate insize;
	int xcount;
	int ycount;
	int zcount;
	grid3d(const grid3d& other) = delete;

	grid3d& operator=(const grid3d&) = delete;
public:
	grid3d(const modelpar& m_par, const rect3d& range, const Coordinate& size);

	const rect3d& get_range() const { return inrange; }

	const Coordinate& get_size() const { return insize; }

	double get_ratio(int f) const;

	int get_xcount() const {
		return xcount;
	}

	int get_ycount() const {
		return ycount;
	}

	int get_zcount() const {
		return zcount;
	}

	size_t get_gridunits_count() const;

	int* const get_facies() const { return facies.get(); }

	void add_unit(const Coordinate& pos3d, int data);
};

class DllEXPORT grid2d {
private:
	vector<int> facies;
	Envelope inrange;
	Coordinate insize;
	int xcount;
	int ycount;
	grid2d(const grid2d& other) = delete;

	grid2d& operator=(const grid2d&) = delete;
public:
	grid2d(const Envelope& range, const Coordinate& size);

	const Envelope& get_range() const { return inrange; }

	const Coordinate& get_size() const { return insize; }

	int get_xcount() const {
		return xcount;
	}

	int get_ycount() const {
		return ycount;
	}

	size_t get_gridunits_count() const;

	const vector<int>& get_facies() const { return facies; }

	void add_unit(const Coordinate& pos, int data);

	void add_unit(const size_t& index, int data);
};

class DllEXPORT channelzone {
private:
	double angel;
	Coordinate start_dot;
	double mainlen;
public:
	channelzone(const double& angel, const Coordinate& startdot, const double& mainlen);

	channelzone(const channelzone& zone);

	double get_angel() const;

	const Coordinate& get_start_dot() const;

	double get_mainlen() const;
};

class DllEXPORT channelpar {
private:
	channelzone zone;
	modelpar m_par;
	double width;
	double thick;
	double wave;
	double amp;
	double levee_w_ratio;
	double levee_b_ratio;
	double levee_c_ratio;
public:
	channelpar(const channelzone& chzone, const modelpar& modelpar, const double& width, const double& thick, const double& wave, const double& amp
	);

	channelpar(const channelpar& chpar);

	double get_width() const;
	double get_thick() const;
	double get_wave() const;
	double get_amp() const;
	double get_levee_w_ratio() const;
	double get_levee_b_ratio() const;
	double get_levee_c_ratio() const;
	const channelzone& get_zone() const;
	const modelpar& get_modelpar() const;
};