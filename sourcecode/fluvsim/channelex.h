#pragma once
#include "section.h"

class __declspec(dllexport) curvature {
private:
	
	double p;
	
	Coordinate dir;
	
	Coordinate targent;
public:
	curvature();

	curvature(const double& p, const Coordinate& dir, const Coordinate& targent);

	const Coordinate& get_dir() const;

	const Coordinate& get_targent() const;

	double get_p() const;

	void set_p(const double& p);

	bool isleft() const;

	double get_ay() const;
};


class __declspec(dllexport) nodeinfo {
private:
	
	double levee_b;
	
	double levee_c;
	
	double levee_w;
	
	double ch_z;
	
	curvature cv;
	
	double ch_width;
	
	double ch_thick;
	nodeinfo(const double& ch_width, const double& ch_thick, const double& ch_z, const double& levee_b, const double& levee_c, const double& levee_w);
public:

	nodeinfo(const double& ch_width, const double& ch_thick, const double& ch_z);

	nodeinfo(const double& ch_width, const double& ch_thick, const double& ch_z, const channelpar* m_par);

	nodeinfo static interpolation(const nodeinfo& info1, const nodeinfo& info2, double k);

	double get_levee_width() const;

	double get_levee_b_thick() const;

	double get_levee_c_thick() const;

	const curvature& get_curv() const;

	curvature& get_change_curv();

	void set_curv(const curvature& value);

	double get_z() const;

	double get_width() const;

	double get_thick() const;

	void set_thick(double thick);

	void set_width(double width);

	void set_z(double z);
};

struct crevasse_info {
private:

	Coordinate pos;

	Coordinate dir;

	double radius;
public:

	crevasse_info(const Coordinate& pos, const Coordinate& dir, const double& radius);

	const Coordinate& get_pos() const;

	const Coordinate& get_dir() const;

	double get_radius() const;
};

class trans_well {
private:
	Coordinate tran_pos;
	int well_id;
	double top;
	double bot;
	facie _facie;
public:
	trans_well(const well3d& well, const Coordinate& tran_pos);

	const Coordinate& get_pos()const;

	int get_wellid() const;

	double get_thick() const { return top - bot; }

	double get_top() const {
		return top;
	}

	double get_bot() const {
		return bot;
	}

	const facie& get_facie() const;
};

class ch_node {
private:
	Coordinate pos;
	double thick;
	double width;
public:
	ch_node(const Coordinate& pos, const double& width, const double& thick);

	const Coordinate& get_pos() const;

	double get_thick() const;

	double get_width() const;
};

enum class move_near_state {

	move_closer = 0,
	
	too_near = 1,

	too_far = 2
};

class move_near_info {
private:
	unique_ptr<vector<Coordinate>> new_center;
public:

	double move_dis;

	move_near_state state;

	move_near_info(vector<Coordinate>* new_center, const double& move_dis, const move_near_state& state)
		:new_center(new_center), move_dis(move_dis), state(state) {
	}

	move_near_info(const move_near_info&) = delete;

	const vector<Coordinate>* get_center() const {
		return new_center.get();
	}
};

class node_change_info {
public:
	size_t line_index;
	double newvalue;
	node_change_info(const size_t& line_index, const double& newvalue) {
		this->line_index = line_index;
		this->newvalue = newvalue;
	}
};

class __declspec(dllexport) wellsfilter {
private:

	unique_ptr<vector<int>> ignoreids;

	unique_ptr<vector<facie>> facies;
public:

	wellsfilter(const vector<int>& ignoreids, const vector<facie>& fs);

	wellsfilter(const wellsfilter& other) = delete;

	wellsfilter& operator=(const wellsfilter&) = delete;

	virtual const vector<int>* get_ignoreids()const;

	virtual const vector<facie>* get_facie_filter()const;
};

class channelex {
public:

	static double interp_width(const vector<nodeinfo>* infos, const LineString* line, const size_t& index, const Coordinate& cross);

	static double interp_thick(const vector<nodeinfo>* infos, const LineString* line, const size_t& index, const Coordinate& cross);

	static nodeinfo interp_nodeinfo(const vector<nodeinfo>* infos, const LineString* line, const size_t& index, const Coordinate& cross);

	static shared_ptr<section> get_channel_section_inzone(LineString* line, vector<nodeinfo>* nodeinfos, vector<crevasse_info>* crevasse_infos, const channelzone& zone, const Coordinate& trans_pos);

	static shared_ptr<section> get_channel_section(LineString* line, vector<nodeinfo>* nodeinfos, vector<crevasse_info>* crevasse_infos, const channelzone& zone, const Coordinate& pos);

	static double mindis_to_nosand_wells(const Coordinate& pos, const vector<trans_well>* nosand_wells);

	static unique_ptr<move_near_info> move_to_point(const Coordinate& point, LineString* line, vector<nodeinfo>* nodeinfos, const double& radius, const double& xy_tol);

	static CoordinateArraySequence* move_away_point(const Coordinate& point, LineString* line, vector<nodeinfo>* nodeinfos, const double& xy_tol, const double& maxradius);

	static bool channel_thick_near(const double& well_thick, const ch_section* section, const double& ratio);

	static std::unique_ptr<vector<Coordinate>> create_ctrl_points(const vector<trans_well>& wells, const channelpar* m_par);
};
