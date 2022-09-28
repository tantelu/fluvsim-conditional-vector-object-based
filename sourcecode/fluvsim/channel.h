#pragma once
#include "channelex.h"
class __declspec(dllexport) channel
{
protected:
	channelpar* m_par = NULL;

	LineString* center = NULL;

	vector<nodeinfo>* nodeinfos = NULL;

	vector<crevasse_info>* crevasse_infos;
public:
	channel(const channelpar& _par);
	channel(const string& file, const modelpar& modelpar);
	channel(const channel& other);
	~channel();
	channel();

	int save_to_file(string file);

	vector<channel> cut_center_with_others(const vector<channel>& others, const size_t& ignoreindex);

	void build_center_nodes(const double& decrease_rate);

	void build_center_nodes(const well3d_sequence* wells);

	void cut_center_with_others(const vector<channel>& others, const well3d& sandwell_pos);

	void move_center_nodes(const well3d_sequence* need_obey_wells, const vector<int>* match_center_sand_ids, const double* maxDis, const double* radius);

	unique_ptr<CoordinateArraySequence> get_cross_chsectionline(const Coordinate* pos, const double& spacing);

	void move_center_nodes(const well3d_sequence* need_obey_wells, const double* wellSpace, const double* stepDis, const int* sandStepNums);

	void center_migration(const double& lateral_dis, const double& forward_dis, const double& up_dis, const double& restrain);

	void adjust_thick_match_sand(const well3d_sequence* need_obey_wells);

	unique_ptr<vector<int>> get_obey_wells(const well3d_sequence* all_wells) const;

	bool is_well_match(const well3d& well) const;

	double get_wells_dis(const well3d_sequence* all_wells, const double& maxdis) const;

	unique_ptr<CoordinateArraySequence> get_centerline() const;

	unique_ptr<CoordinateArraySequence> buffer_centerline(const double& distance) const;

	unique_ptr<vector<curvature>> get_curvatures() const;

	const LineString* get_centerline_internal() const;

	const channelpar* get_par() const;

	void set_par(const channelpar& m_par);

	virtual void get_in_channel_pos(grid3d* inunits) const;

	virtual void get_in_channel_pos(grid2d* inunits) const;

	void get_migrate_in_channel_pos(grid3d* inindexs) const;

	virtual facie get_channel_facie(Coordinate pos) const;

	void add_crevasse_pos(const Coordinate& pos, double radius);

	unique_ptr<vector<well3d>> in_ch_segment(const well3d_sequence& wells, const wellsfilter& filter) const;

	double match_ratio(const well3d_sequence& wells, const wellsfilter& filter) const;

	void well_cross_info(const well3d_sequence& wells, const double& maxdis, double& welllen, double& cross, double& sand, double& mud) const;

	const size_t get_nodes_count() const;

	Coordinate get_node_at(size_t index) const;

	nodeinfo get_info_at(size_t index) const;

private:
	void clear();

	void set_center(LineString* line);

	void set_nodeinfos(vector<nodeinfo>* nodeinfos);

	bool is_conflict_with_nosand_wells(const vector<trans_well>* nosand_obey_wells);

	bool is_sandwell_incenter(const vector<trans_well>* sand_wells);

	void adjust_width_match_nosand(const vector<trans_well>* nosand_wells);

	void adjust_thick_match_sand(const vector<trans_well>* need_obey_wells, const bool& bigchange);

	void adjust_ch_z(const vector<trans_well>* need_obey_wells);

	void re_computer_curvature();

	void re_interpo_center_nodes();
};

