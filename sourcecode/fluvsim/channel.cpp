#include "channel.h"
#include "pch.h"
#include "bezier.h"
#include <geos\algorithm\Angle.h>
#include <geos/index/kdtree/KdTree.h>
#include "geos\operation\union\UnaryUnionOp.h"
#include "mathhelp.h"
#include "f2c.h"
#include "blaswrap.h"
#include <stdlib.h>
#include <time.h>
#include <map>
#include <iostream>
#include <fstream>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>

extern "C"
{
#include <clapack.h>
}

using namespace geos::index::kdtree;
using namespace geos::operation::geounion;

channel::~channel()
{
	clear();
	if (m_par != NULL) {
		delete m_par;
		m_par = NULL;
	}
}

channel::channel()
{
	this->nodeinfos = new vector<nodeinfo>();
	crevasse_infos = new vector<crevasse_info>();
}

int channel::save_to_file(string file)
{
	ofstream out_file(file);
	if (out_file.is_open()) {
		//zone
		auto& zone = this->get_par()->get_zone();
		out_file << zone.get_start_dot().x << "," << zone.get_start_dot().y << "," << zone.get_start_dot().z << ",";
		out_file << zone.get_mainlen() << "," << zone.get_angel() << "\n";
		//par
		out_file << this->get_par()->get_width() << "," << this->get_par()->get_thick() << "," << this->get_par()->get_wave() << "," << this->get_par()->get_amp() << "\n";
		//numscount
		int nums = this->center->getNumPoints();
		out_file << nums << endl;
		//dots
		auto dots = this->center->getCoordinatesRO();
		for (size_t i = 0; i < nums; i++)
		{
			auto& dot = dots->getAt(i);
			out_file << dot.x << "," << dot.y << "," << dot.z << "\n";
		}
		out_file.flush();
		//infos
		for (size_t i = 0; i < nums; i++)
		{
			out_file << this->nodeinfos->at(i).get_z() << "," << this->nodeinfos->at(i).get_width() << "," << this->nodeinfos->at(i).get_thick() << ",";
			//cv
			auto& cv = this->nodeinfos->at(i).get_curv();
			out_file << cv.get_ay() << "," << cv.get_p() << "," << cv.get_dir().x << "," << cv.get_dir().y << "," << cv.get_targent().x << "," << cv.get_targent().y << "\n";
		}
		out_file.close();
		return 1;
	}
	else {
		return 0;
	}
}

channel::channel(const channelpar& _par)
{
	this->nodeinfos = new vector<nodeinfo>();
	m_par = new channelpar(_par);
	crevasse_infos = new vector<crevasse_info>();
}

channel::channel(const string& file, const modelpar& modelpar)
{
	ifstream infile;
	infile.open(file, ios::in);
	if (infile.is_open())
	{
		string buf;
		//zone
		channelzone* zone = NULL;
		{
			getline(infile, buf);
			vector<string> res;
			boost::split(res, buf, boost::is_any_of(","), boost::token_compress_on);
			if (res.size() == 5) {
				auto x = std::stod(res[0]);
				auto y = std::stod(res[1]);
				auto z = std::stod(res[2]);
				auto len = std::stod(res[3]);
				auto angel = std::stod(res[4]);
				zone = new channelzone(angel, Coordinate(x, y, z), len);
			}
		}
		//河道宽度 最大厚度 波长 振幅
		double width = 0, thick = 0, wave = 0, amp = 0;
		{
			getline(infile, buf);
			vector<string> res;
			boost::split(res, buf, boost::is_any_of(","), boost::token_compress_on);
			if (res.size() == 4) {
				width = std::stod(res[0]);
				thick = std::stod(res[1]);
				wave = std::stod(res[2]);
				amp = std::stod(res[3]);
			}
		}
		//dotnums
		getline(infile, buf);
		int dotnums = std::stoi(buf);
		auto factory = GeometryFactory::getDefaultInstance();
		CoordinateArraySequence* cas = new CoordinateArraySequence();
		for (int i = 0; i < dotnums; i++)
		{
			getline(infile, buf);
			vector<string> res;
			boost::split(res, buf, boost::is_any_of(","), boost::token_compress_on);
			if (res.size() == 3) {
				if (res[2] != "nan") {
					cas->add(Coordinate(std::stod(res[0]), std::stod(res[1]), std::stod(res[2])));
				}
				else {
					cas->add(Coordinate(std::stod(res[0]), std::stod(res[1])));
				}
			}
		}
		set_center(factory->createLineString(cas));

		auto newinfos = new vector<nodeinfo>();
		newinfos->reserve(dotnums);
		for (int i = 0; i < dotnums; i++)
		{
			getline(infile, buf);
			vector<string> res;
			boost::split(res, buf, boost::is_any_of(","), boost::token_compress_on);
			if (res.size() == 9) {
				auto chz = std::stod(res[0]);
				auto width = std::stod(res[1]);
				auto thick = std::stod(res[2]);
				auto ay = std::stod(res[3]);
				auto p = std::stod(res[4]);
				double dirx = NAN;
				double diry = NAN;
				double targentx = NAN;
				double targenty = NAN;
				if (res[5] != "nan") {
					dirx = std::stod(res[5]);
				}
				if (res[6] != "nan") {
					diry = std::stod(res[6]);
				}if (res[7] != "nan") {
					targentx = std::stod(res[7]);
				}if (res[8] != "nan") {
					targenty = std::stod(res[8]);
				}
				curvature cur(p, Coordinate(dirx, diry), Coordinate(targentx, targenty));
				newinfos->emplace_back(width, thick, chz);
				newinfos->at(newinfos->size() - 1).set_curv(cur);
			}
		}
		set_nodeinfos(newinfos);
		this->m_par = new channelpar(*zone, modelpar, width, thick, wave, amp);
		delete zone;
		zone = NULL;
	}
	crevasse_infos = new vector<crevasse_info>();
}

channel::channel(const channel& other)
{
	this->nodeinfos = new vector<nodeinfo>(*other.nodeinfos);
	m_par = new channelpar(*other.get_par());
	crevasse_infos = new vector<crevasse_info>(*other.crevasse_infos);
	auto copycenter = GeometryFactory::getDefaultInstance()->createLineString(*other.center->getCoordinatesRO());
	set_center(copycenter);
}

void channel::move_center_nodes(const well3d_sequence* need_obey_wells, const vector<int>* match_center_sand_ids, const double* maxDis, const double* radius)
{
	if (need_obey_wells != nullptr && need_obey_wells->get_size() > 0 && *maxDis > 0 && *radius > 0) {
		vector<trans_well> nosandwell;//需要移动来匹配的非砂岩井
		vector<trans_well> sandwell;//需要移动来匹配的砂岩井
		vector<trans_well> match_center_sandwell;//需要锚定在中线上的砂岩井
		auto factory = GeometryFactory::getDefaultInstance();
		auto distance = *radius + this->m_par->get_width();
		auto pt = this->center->buffer(distance, 8, 2);
		auto polygon = (geos::geom::Polygon*)pt.get();
		for (size_t i = 0; i < need_obey_wells->get_size(); i++) {

			Coordinate dot = geo_math::cor_translate(need_obey_wells->get_at(i).get_pos(), this->m_par->get_zone().get_start_dot(), this->m_par->get_zone().get_angel());
			auto dot2 = factory->createPoint(dot);
			if (polygon->contains(dot2)) {
				if (need_obey_wells->get_at(i).get_facie() == facie::channel_sand) {
					if (match_center_sand_ids != nullptr && match_center_sand_ids->size() > 0) {
						if (count(match_center_sand_ids->begin(), match_center_sand_ids->end(), need_obey_wells->get_at(i).get_id()) > 0) {
							match_center_sandwell.emplace_back(need_obey_wells->get_at(i), dot);
						}
						else {
							sandwell.emplace_back(need_obey_wells->get_at(i), dot);
						}
					}
					else {
						sandwell.emplace_back(need_obey_wells->get_at(i), dot);
					}
				}
				else if (need_obey_wells->get_at(i).get_facie() != facie::channel_sand) {
					nosandwell.emplace_back(need_obey_wells->get_at(i), dot);
				}
			}
			factory->destroyGeometry(dot2);
		}
		double acc_dis = 0;
		while (acc_dis < *maxDis) {//先处理非砂岩
			acc_dis += 0.5 * this->m_par->get_modelpar().get_xytol();
			for (size_t i = 0; i < nosandwell.size(); i++) {
				double movedis = 0;
				auto cas = channelex::move_away_point(nosandwell[i].get_pos(), center, this->nodeinfos, this->m_par->get_modelpar().get_xytol(), *maxDis);
				if (cas != NULL) {
					set_center(factory->createLineString(cas));
				}
			}
		}
		if (is_conflict_with_nosand_wells(&nosandwell)) {//如果冲突修改宽度
			adjust_width_match_nosand(&nosandwell);
		}
		if (is_conflict_with_nosand_wells(&nosandwell)) {
			//throw exception("非河道没有条件化");
		}
		if (match_center_sandwell.size() == 0) {
			for (size_t i = 0; i < sandwell.size(); i++) {
				//判断如何移动 砂岩需要靠近和远离
				auto pt = channelex::get_channel_section_inzone(this->center, this->nodeinfos, this->crevasse_infos, this->m_par->get_zone(), sandwell[i].get_pos());
				auto sect = pt.get();
				auto copy = factory->createLineString(*this->center->getCoordinatesRO());
				if (sect->get_face() != facie::channel_sand || ((ch_section*)sect)->get_thick() < sandwell[i].get_thick()) {//靠近砂岩井点
					double acc_dis = 0;
					while (acc_dis < *maxDis) {
						double movedis = 0;
						auto pt = channelex::move_to_point(sandwell[i].get_pos(), center, this->nodeinfos, *radius, this->m_par->get_modelpar().get_xytol());
						auto move_res = pt.get();
						if (move_res->state == move_near_state::move_closer) {
							acc_dis += move_res->move_dis;
							auto moved = move_res->get_center();
							CoordinateArraySequence* cas = new CoordinateArraySequence();
							for (size_t index = 0; index < moved->size(); index++)
							{
								cas->add(Coordinate(moved->at(index).x, moved->at(index).y));
							}
							set_center(factory->createLineString(cas));
							auto temppt = channelex::get_channel_section_inzone(this->center, this->nodeinfos, this->crevasse_infos, this->m_par->get_zone(), sandwell[i].get_pos());
							auto tempsect = temppt.get();
							if (tempsect->get_face() == facie::channel_sand && ((ch_section*)tempsect)->get_thick() >= sandwell[i].get_thick()) {
								break;//井点处河道厚度已经大于井点 结束靠近
							}
						}
						else if (move_res->state == move_near_state::too_near) {//在中线上
							break;
						}
						else {//(move_res->state == move_near_state::too_far)
							break;//不能移动跳出
						}
					}
				}
				else {//远离
					double acc_dis = 0;
					while (acc_dis < *maxDis) {
						acc_dis += 0.5 * this->m_par->get_modelpar().get_xytol();
						auto cas = channelex::move_away_point(sandwell[i].get_pos(), center, this->nodeinfos, this->m_par->get_modelpar().get_xytol(), *maxDis);
						set_center(factory->createLineString(cas));
						auto temppt = channelex::get_channel_section_inzone(this->center, this->nodeinfos, this->crevasse_infos, this->m_par->get_zone(), sandwell[i].get_pos());
						auto tempsect = temppt.get();
						if (tempsect->get_face() == facie::channel_sand && ((ch_section*)tempsect)->get_thick() <= sandwell[i].get_thick()) {
							break;//井点处河道厚度已经小于井点 结束远离
						}
					}
				}
				if (is_conflict_with_nosand_wells(&nosandwell)) {
					set_center(copy);
				}
				else {
					auto pt2 = channelex::get_channel_section_inzone(this->center, this->nodeinfos, this->crevasse_infos, this->m_par->get_zone(), sandwell[i].get_pos());
					auto sect2 = pt2.get();
					if (sect2->get_face() == facie::channel_sand) {
						auto subs = ((ch_section*)sect2)->get_thick() - sandwell[i].get_thick();
						if (abs(subs) > 0.5 * this->m_par->get_modelpar().get_ztol()) {
							set_center(copy);
							continue;
						}
					}
					factory->destroyGeometry(copy);
				}
			}
		}
		re_interpo_center_nodes();
		adjust_ch_z(&sandwell);
		adjust_thick_match_sand(&sandwell, false);
	}
}

unique_ptr<CoordinateArraySequence> channel::get_cross_chsectionline(const Coordinate* pos, const double& spacing)
{
	//转换坐标
	Coordinate trans_pos = geo_math::cor_translate(*pos, this->m_par->get_zone().get_start_dot(), this->m_par->get_zone().get_angel());
	size_t index;
	Coordinate cross;
	auto dis = geo_math::point_min_pos_to_line(trans_pos, this->center, index, cross);
	auto nodeinfo = nodeinfos->at(index);
	if (index < this->center->getCoordinatesRO()->size() - 1) {
		nodeinfo = channelex::interp_nodeinfo(nodeinfos, this->center, index, cross);
	}
	auto ch_z = nodeinfos->at(index).get_z();
	auto cv = nodeinfos->at(index).get_curv();
	bool transisleft = geo_math::cross(cv.get_targent(), Coordinate(trans_pos.x - cross.x, trans_pos.y - cross.y)) > 0;
	double wid_pos = nodeinfo.get_width() * 0.5 + (transisleft ? -1 : 1) * dis;
	ch_section sect(nodeinfo.get_width(), nodeinfo.get_thick(), cv.get_ay(), ch_z, wid_pos);
	return sect.get_bezier_lines(spacing);
}

void channel::move_center_nodes(const well3d_sequence* need_obey_wells, const double* wellSpace, const double* stepDis, const int* sandStepNums)
{
	if (need_obey_wells != nullptr && need_obey_wells->get_size() > 0 && *stepDis > 0 && *sandStepNums >= 0) {
		vector<trans_well> wells;
		auto factory = GeometryFactory::getDefaultInstance();
		auto Influencedis = *wellSpace;
		auto checkDis = *wellSpace / 2.0;
		//砂岩移动步长通过参数传递
		for (size_t i = 0; i < need_obey_wells->get_size(); i++) {

			Coordinate dot = geo_math::cor_translate(need_obey_wells->get_at(i).get_pos(), this->m_par->get_zone().get_start_dot(), this->m_par->get_zone().get_angel());
			auto dis = geo_math::point_min_pos_to_line(dot, this->center);
			if (dis <= checkDis) {//检查范围
				wells.emplace_back(need_obey_wells->get_at(i), dot);
			}
		}
		for (size_t stepi = 0; stepi < *sandStepNums; stepi++)
		{
			for (size_t welli = 0; welli < wells.size(); welli++) {
				if (wells[welli].get_facie() == facie::channel_sand) {
					auto temppt = channelex::get_channel_section_inzone(this->center, this->nodeinfos, this->crevasse_infos, this->m_par->get_zone(), wells[welli].get_pos());
					auto tempsect = temppt.get();
					if (tempsect->get_face() == facie::channel_sand && (((ch_section*)tempsect)->is_in_maxthick(*stepDis) || ((ch_section*)tempsect)->get_thick() >= wells[welli].get_thick())) {
						continue;//井在最厚处或者已经大于河道厚度了
					}
					auto move_res = channelex::move_to_point(wells[welli].get_pos(), center, this->nodeinfos, Influencedis, *stepDis);
					if (move_res->state == move_near_state::move_closer) {
						auto moved = move_res->get_center();
						CoordinateArraySequence* cas = new CoordinateArraySequence();
						for (size_t index = 0; index < moved->size(); index++)
						{
							cas->add(Coordinate(moved->at(index).x, moved->at(index).y));
						}
						set_center(factory->createLineString(cas));
						auto temppt = channelex::get_channel_section_inzone(this->center, this->nodeinfos, this->crevasse_infos, this->m_par->get_zone(), wells[welli].get_pos());
						auto tempsect = temppt.get();
					}
				}
			}
		}
		for (size_t welli = 0; welli < wells.size(); welli++) {
			if (wells[welli].get_facie() == facie::_default) {
				size_t index;
				Coordinate cross;
				auto mindis = geo_math::point_min_pos_to_line(wells[welli].get_pos(), this->center, index, cross);
				auto width = nodeinfos->at(index).get_width();
				if (index < this->center->getNumPoints() - 1) {
					width = channelex::interp_width(nodeinfos, this->center, index, cross);
				}
				if (mindis < width / 2.0) {
					//河道宽度应该小于井距的一半 否则后面的井移动结果会影响前面已经移动好的
					int needmovesteps = (int)(geo_math::random(width / 1.8 - mindis, checkDis - mindis) / (*stepDis)); //移动距离在半宽和检测距离之间,
					for (size_t stepi = 0; stepi < needmovesteps; stepi++)
					{
						auto cas = channelex::move_away_point(wells[welli].get_pos(), center, this->nodeinfos, *stepDis, Influencedis);
						if (cas != NULL) {
							set_center(factory->createLineString(cas));
						}
					}
					auto mindis2 = geo_math::point_min_pos_to_line(wells[welli].get_pos(), this->center, index, cross);
					auto width2 = nodeinfos->at(index).get_width();
					if (index < this->center->getNumPoints() - 1) {
						width2 = channelex::interp_width(nodeinfos, this->center, index, cross);
					}
				}
			}
		}
		re_interpo_center_nodes();
	}
}

void channel::center_migration(const double& lateral_dis, const double& forward_dis, const double& up_dis, const double& restrain)
{
	CoordinateArraySequence* cas = new CoordinateArraySequence();
	auto cors = this->center->getCoordinatesRO();
	vector<tuple<size_t, double>> ps;
	ps.reserve(cas->size() / 3);
	ps.emplace_back(0, 0);
	int curinte = 0;
	for (size_t i = 1; i < cors->size() - 1; i++)
	{
		curinte++;
		auto p1 = this->nodeinfos->at(i - 1).get_curv().get_p();
		auto p2 = this->nodeinfos->at(i).get_curv().get_p();
		auto p3 = this->nodeinfos->at(i + 1).get_curv().get_p();
		if ((p1 - p2 < 0 && p2 - p3 > 0 && p2 > 0.5)) {
			auto dir = this->nodeinfos->at(i).get_curv().get_dir();
			auto tagent = geo_math::sub(cors->getAt(i), cors->getAt(i - 1));
			ps.emplace_back(i, lateral_dis * /*p2 **/ (geo_math::cross(tagent, dir) > 0 ? 1 : -1));
			curinte = 0;
		}
		else if (curinte >= 8) {
			auto dir = this->nodeinfos->at(i).get_curv().get_dir();
			auto tagent = geo_math::sub(cors->getAt(i), cors->getAt(i - 1));
			ps.emplace_back(i, lateral_dis * p2 * (geo_math::cross(tagent, dir) > 0 ? 1 : -1));
			curinte = 0;
		}
	}
	ps.emplace_back(cors->size() - 1, this->nodeinfos->at(cors->size() - 1).get_curv().get_p());
	cas->add(cors->getAt(0));
	bool crossout = false;
	for (size_t i = 0; i < ps.size() - 1; i++)
	{
		size_t start = get<0>(ps[i]);
		size_t end = get<0>(ps[i + 1]);
		double startp = get<1>(ps[i]);
		double endp = get<1>(ps[i + 1]);
		size_t index = start;
		while (index < end) {
			if (index > 0) {
				double move_dis = (startp + (index - start) * (endp - startp) / (end - start));
				auto tagent = geo_math::normal(cors->getAt(index - 1), cors->getAt(index));
				auto dir = this->nodeinfos->at(index).get_curv().get_dir();
				auto p = this->nodeinfos->at(index).get_curv().get_p();
				auto k = geo_math::cross(tagent, dir) * move_dis;
				auto x = cors->getAt(index).x + dir.x * (k > 0 ? 1 : -1) * abs(move_dis) + forward_dis;
				if (x > this->m_par->get_zone().get_mainlen()) {
					i = ps.size() - 1;
					crossout = true;
					break;
				}
				double yk = 1.0;
				if (restrain > 0) {
					double cory = abs(cors->getAt(index).y);
					if (cory > (restrain - this->m_par->get_width())) {
						yk = 0.0;
					}
					else {
						yk = (restrain - cory) / restrain;
					}
				}
				auto y = cors->getAt(index).y + dir.y * yk * (k > 0 ? 1 : -1) * abs(move_dis);
				cas->add(Coordinate(x, y));
			}
			index++;
		}
	}
	if (!crossout)
	{
		auto index = get<0>(ps[ps.size() - 1]);
		auto move_dis = get<1>(ps[ps.size() - 1]);
		auto tagent = geo_math::normal(cors->getAt(index - 1), cors->getAt(index));
		auto dir = this->nodeinfos->at(index).get_curv().get_dir();
		auto p = this->nodeinfos->at(index).get_curv().get_p();
		auto k = geo_math::cross(tagent, dir) * move_dis;
		auto x = cors->getAt(index).x + dir.x * (k > 0 ? 1 : -1) * abs(move_dis);
		auto y = cors->getAt(index).y + dir.y * (k > 0 ? 1 : -1) * abs(move_dis);
		cas->add(Coordinate(x, y));
	}

	double cutdis = this->m_par->get_width() * 1.0;
	KdTree tree;
	for (size_t i = 0; i < cas->size(); i++)
	{
		tree.insert(cas->getAt(i), (void*)i);
	}
	vector<KdNode*> result;
	vector<size_t> cutoff;
	int min_index_delt = max(this->m_par->get_width() * 2.4 / this->m_par->get_modelpar().get_ctrl_dis(), 5);
	CoordinateArraySequence* cutcas = new CoordinateArraySequence();
	vector<nodeinfo>* newinfos = new vector<nodeinfo>();
	for (size_t i = 0; i < cas->size(); i++)
	{
		auto curpt = cas->getAt(i);
		newinfos->push_back(this->nodeinfos->at(i));
		Envelope range(curpt.x - cutdis, curpt.x + cutdis, curpt.y - cutdis, curpt.y + cutdis);
		tree.query(range, result);
		for (size_t index = 0; index < result.size(); index++)
		{
			auto findindex = (int)result[index]->getData();
			if (findindex > i && (findindex - i) > min_index_delt) {
				i = findindex + 1;
				break;
			}
		}
		cutcas->add(curpt);
	}
	if (cutcas->getAt(cutcas->getSize() - 1).x < (this->m_par->get_zone().get_mainlen() - this->m_par->get_modelpar().get_ctrl_dis())) {
		cutcas->add(Coordinate(this->m_par->get_zone().get_mainlen(), cutcas->getAt(cutcas->getSize() - 1).y));
		newinfos->push_back(newinfos->at(newinfos->size() - 1));
	}
	auto factory = GeometryFactory::getDefaultInstance();
	set_center(factory->createLineString(cutcas));
	set_nodeinfos(newinfos);
	re_interpo_center_nodes();
	for (size_t i = 0; i < this->nodeinfos->size(); i++)
	{
		this->nodeinfos->at(i).set_z(this->nodeinfos->at(i).get_z() + up_dis);
	}
}

//生成过程中避开非河道井
void channel::build_center_nodes(const well3d_sequence* wells)
{
	this->clear();
	double inte = (this->m_par->get_wave() / 2.0);
	vector<trans_well> trans;
	auto factory = GeometryFactory::getDefaultInstance();
	if (wells != nullptr) {
		trans.reserve(wells->get_size());
		for (size_t i = 0; i < wells->get_size(); i++) {
			if (this->m_par->get_zone().get_start_dot().z <= wells->get_at(i).get_top() &&
				this->m_par->get_zone().get_start_dot().z > wells->get_at(i).get_bot())
			{
				Coordinate dot = geo_math::cor_translate(wells->get_at(i).get_pos(), this->m_par->get_zone().get_start_dot(), this->m_par->get_zone().get_angel());
				trans.emplace_back(wells->get_at(i), dot);
			}
		}
	}
	auto linept = channelex::create_ctrl_points(trans, this->m_par);
	auto line = linept.get();
	if (line->size() > 1) {
		CoordinateArraySequence* cas = new CoordinateArraySequence();
		auto newinfos = new vector<nodeinfo>();
		for (size_t i = 0; i < line->size(); i++)
		{
			cas->add(line->at(i));
			auto width = this->m_par->get_width() * (1 - 0.2 * (rand() / double(RAND_MAX) * 2 - 1));
			auto thick = this->m_par->get_thick() * (1 - 0.2 * (rand() / double(RAND_MAX) * 2 - 1));
			newinfos->emplace_back(width, thick, this->m_par->get_zone().get_start_dot().z, this->m_par);
		}
		set_nodeinfos(newinfos);
		set_center(factory->createLineString(cas));
	}
	else {
		CoordinateArraySequence* cas = new CoordinateArraySequence();
		auto newinfos = new vector<nodeinfo>();
		cas->add(Coordinate(0, 0));
		cas->add(Coordinate(this->m_par->get_zone().get_mainlen(), 0));
		newinfos->emplace_back(this->m_par->get_width(), this->m_par->get_thick(), this->m_par->get_zone().get_start_dot().z, this->m_par);
		newinfos->emplace_back(this->m_par->get_width(), this->m_par->get_thick(), this->m_par->get_zone().get_start_dot().z, this->m_par);
		set_nodeinfos(newinfos);
		set_center(factory->createLineString(cas));
	}
	re_interpo_center_nodes();
}

void channel::cut_center_with_others(const vector<channel>& others, const well3d& sandwell_pos)
{
	if (others.size() > 0) {

		Coordinate tran_pos = geo_math::cor_translate(sandwell_pos.get_pos(), this->m_par->get_zone().get_start_dot(), this->m_par->get_zone().get_angel());
		Coordinate start_dot(this->m_par->get_zone().get_start_dot().x, this->m_par->get_zone().get_start_dot().y);

		auto factory = GeometryFactory::getDefaultInstance();
		size_t cur = 0;
		bool indexchange = true;
		for (size_t i = 0; i < others.size(); i++)
		{
			LocationIndexedLine loc(this->center);
			if (indexchange) {//剪切之后重新计算井点的索引 
				//loc = LocationIndexedLine(this->center);
				Coordinate cross;
				auto dis = geo_math::point_min_pos_to_line(tran_pos, this->center, cur, cross);
			}
			if (cur > 0 && cur < this->center->getNumPoints() - 1)
			{
				auto otherpt = others[i].get_centerline();
				auto otherdots = otherpt.get();
				CoordinateArraySequence* cas = new CoordinateArraySequence();
				for (size_t otheri = 0; otheri < otherdots->size(); otheri++)
				{
					auto otherdot = otherdots->getAt(otheri);
					auto tran_otherdot = geo_math::cor_translate(otherdot, start_dot, this->m_par->get_zone().get_angel());
					cas->add(tran_otherdot);
				}
				auto otherline = factory->createLineString(cas);
				auto crosspt = this->center->intersection(otherline);
				factory->destroyGeometry(otherline);
				auto cross = crosspt.get();
				if (cross != NULL && !cross->isEmpty()) {
					vector<Coordinate> cross_dots;
					if (cross->getGeometryTypeId() == GeometryTypeId::GEOS_MULTIPOINT) {
						MultiPoint* mul_point = (MultiPoint*)cross;
						for (size_t j = 0; j < mul_point->getNumGeometries(); j++)
						{
							auto point = mul_point->getGeometryN(j)->getCoordinate();
							cross_dots.push_back(*point);
						}
					}
					else if (cross->getGeometryTypeId() == GeometryTypeId::GEOS_POINT)
					{
						Point* point = (Point*)cross;
						cross_dots.push_back(*point->getCoordinate());
					}
					else {
						throw exception("交点类型错误");
					}
					size_t minindex = 0;
					size_t maxindex = this->center->getNumPoints() - 1;
					Coordinate* firstcross = NULL;
					Coordinate* secondcross = NULL;
					for (size_t j = 0; j < cross_dots.size(); j++)
					{
						auto pos = loc.indexOf(cross_dots[j]);
						auto temp = pos.getSegmentIndex();
						if (temp < cur) {
							if (minindex < temp) {
								minindex = temp;
								firstcross = &cross_dots[j];
							}
						}
						else if (temp > cur) {
							if (maxindex > temp) {
								maxindex = temp;
								secondcross = &cross_dots[j];
							}
						}
					}
					if (minindex > 0 || maxindex < this->center->getNumPoints() - 1) {//切断中线
						if (cur<minindex || cur>maxindex)
						{
							throw exception("");
						}
						auto corseq = this->center->getCoordinatesRO();
						CoordinateArraySequence* cas = new CoordinateArraySequence();
						auto newinfos = new vector<nodeinfo>();
						newinfos->reserve(maxindex - minindex + 2);
						//添加第一个交点
						if (firstcross != NULL) {
							cas->add(*firstcross);
							newinfos->emplace_back(this->nodeinfos->at(minindex));
						}
						else {
							cas->add(corseq->getAt(minindex));
							newinfos->emplace_back(this->nodeinfos->at(minindex));
						}
						for (size_t index = minindex + 1; index <= maxindex; index++)
						{
							cas->add(corseq->getAt(index));
							newinfos->emplace_back(this->nodeinfos->at(index));
						}
						//添加第二个交点
						if (secondcross != NULL) {
							cas->add(*secondcross);
							newinfos->emplace_back(this->nodeinfos->at(maxindex));
						}
						set_center(factory->createLineString(cas));
						set_nodeinfos(newinfos);
						indexchange = true;
					}
					else {
						indexchange = false;
					}
				}
			}
		}
	}
}

vector<channel> channel::cut_center_with_others(const vector<channel>& others, const size_t& ignoreindex)
{
	vector<channel> chs;
	if (others.size() > 0) {
		auto factory = GeometryFactory::getDefaultInstance();
		Coordinate start_dot(this->m_par->get_zone().get_start_dot().x, this->m_par->get_zone().get_start_dot().y);

		set<size_t> cutindexs;
		for (size_t i = 0; i < others.size(); i++)
		{
			if (i == ignoreindex) {//自己不切自己
				continue;
			}
			LocationIndexedLine loc(this->center);

			auto otherpt = others.at(i).get_centerline();//不同中线坐标系不同 需要转到当前坐标系下
			auto otherdots = otherpt.get();
			CoordinateArraySequence* cas = new CoordinateArraySequence();
			for (size_t otheri = 0; otheri < otherdots->size(); otheri++)
			{
				auto otherdot = otherdots->getAt(otheri);
				auto tran_otherdot = geo_math::cor_translate(otherdot, start_dot, this->m_par->get_zone().get_angel());
				cas->add(tran_otherdot);
			}
			auto otherline = factory->createLineString(cas);
			auto crosspt = this->center->intersection(otherline);
			factory->destroyGeometry(otherline);
			auto cross = crosspt.get();
			if (cross != NULL && !cross->isEmpty()) {
				vector<Coordinate> cross_dots;
				if (cross->getGeometryTypeId() == GeometryTypeId::GEOS_MULTIPOINT) {
					MultiPoint* mul_point = (MultiPoint*)cross;
					for (size_t j = 0; j < mul_point->getNumGeometries(); j++)
					{
						auto point = mul_point->getGeometryN(j)->getCoordinate();
						cross_dots.push_back(*point);
					}
				}
				else if (cross->getGeometryTypeId() == GeometryTypeId::GEOS_POINT)
				{
					Point* point = (Point*)cross;
					cross_dots.push_back(*point->getCoordinate());
				}
				/*else {
					throw exception("交点类型错误");
				}*/
				for (size_t j = 0; j < cross_dots.size(); j++)//添加所有交点
				{
					auto pos = loc.indexOf(cross_dots[j]);
					auto temp = pos.getSegmentIndex();
					cutindexs.insert(temp);
				}
			}
		}
		int lastcutindex = 0;
		auto corseq = this->center->getCoordinatesRO();
		for (auto i = cutindexs.begin(); i != cutindexs.end(); i++)
		{
			if (*i - lastcutindex < 5) {//节点数太少不切
				continue;
			}
			else {//构建
				channel cutch(*this->get_par());
				CoordinateArraySequence* cas = new CoordinateArraySequence();
				auto newinfos = new vector<nodeinfo>();
				newinfos->reserve(*i - lastcutindex + 1);
				for (int curi = lastcutindex; curi < *i; curi++)
				{
					cas->add(corseq->getAt(curi));
					newinfos->emplace_back(this->nodeinfos->at(curi));
				}
				cutch.set_center(factory->createLineString(cas));
				cutch.set_nodeinfos(newinfos);
				chs.push_back(cutch);
				lastcutindex = *i;
			}
		}
		if (this->center->getNumPoints() - lastcutindex > 4) {
			channel cutch(*this->get_par());
			CoordinateArraySequence* cas = new CoordinateArraySequence();
			auto newinfos = new vector<nodeinfo>();
			newinfos->reserve(this->center->getNumPoints() - lastcutindex + 1);
			for (int curi = lastcutindex; curi < this->center->getNumPoints(); curi++)
			{
				cas->add(corseq->getAt(curi));
				newinfos->emplace_back(this->nodeinfos->at(curi));
			}
			cutch.set_center(factory->createLineString(cas));
			cutch.set_nodeinfos(newinfos);
			chs.push_back(cutch);
		}
	}
	if (chs.size() == 0) {
		chs.push_back(*this);
	}
	return chs;
}

void channel::build_center_nodes(const double& decrease_rate)
{
	this->clear();
	auto factory = GeometryFactory::getDefaultInstance();
	double curx = 0;
	double max_x = this->get_par()->get_zone().get_mainlen();
	double inte = this->get_par()->get_wave();
	CoordinateArraySequence* cas = new CoordinateArraySequence();
	auto newinfos = new vector<nodeinfo>();
	while (curx < max_x) {
		auto randvalue = rand() / double(RAND_MAX) * 2 - 1;
		double y = randvalue * this->get_par()->get_amp();
		cas->add(Coordinate(curx, y));
		double width = this->m_par->get_width() * max(0.1, (1 - curx * decrease_rate / 100.0));
		auto thick = this->m_par->get_thick() * max(0.1, (1 - curx * decrease_rate / 100.0));
		newinfos->emplace_back(width, thick, this->m_par->get_zone().get_start_dot().z, this->m_par);
		curx += inte;
	}
	cas->setAt(Coordinate(0, 0), 0);
	set_nodeinfos(newinfos);
	set_center(factory->createLineString(cas));
	re_interpo_center_nodes();
}

void channel::re_interpo_center_nodes()
{
	auto nodes = this->center->getCoordinatesRO();
	if (nodes->size() != this->nodeinfos->size()) {
		throw exception("长度不一致");
	}
	auto len = this->center->getLength();
	double inte = min(this->m_par->get_modelpar().get_ctrl_dis(), len / 5);
	CoordinateArraySequence* cas = new CoordinateArraySequence();
	vector<Coordinate> points;
	points.reserve(nodes->size());
	for (size_t i = 0; i < nodes->size(); i++)
	{
		points.emplace_back(nodes->getAt(i));
	}
	vector<Coordinate> first;
	vector<Coordinate> second;
	auto res = geo_math::calculate_control_point(points, first, second);
	auto newinfos = new vector<nodeinfo>();
	if (res) {
		newinfos->reserve(nodes->size() * 2);
		double cur_t = 0;//线段的当前位置
		for (size_t i = 0; i < nodes->size() - 1; i++)
		{
			auto len = points[i].distance(points[i + 1]);//线段总长度
			auto restdis = len;//线段剩余
			while (len > cur_t)//当前位置还在线段内
			{
				vector<Bezier::Point> four = { {Bezier::Point(points[i].x,points[i].y)}, {Bezier::Point(first[i].x,first[i].y)}
				, {Bezier::Point(second[i].x,second[i].y)} , {Bezier::Point(points[i + 1].x,points[i + 1].y)} };
				Bezier::Bezier<3> bezier(four);
				auto x = bezier.valueAt(cur_t / len, 0);
				auto y = bezier.valueAt(cur_t / len, 1);
				cas->add(Coordinate(x, y));
				auto k = cur_t / len;
				auto width = this->nodeinfos->at(i).get_width() + (this->nodeinfos->at(i + 1).get_width() - this->nodeinfos->at(i).get_width()) * k;
				auto thick = this->nodeinfos->at(i).get_thick() + (this->nodeinfos->at(i + 1).get_thick() - this->nodeinfos->at(i).get_thick()) * k;
				auto ch_z = this->nodeinfos->at(i).get_z() + (this->nodeinfos->at(i + 1).get_z() - this->nodeinfos->at(i).get_z()) * k;
				newinfos->emplace_back(width, thick, ch_z, this->m_par);
				restdis = len - cur_t;
				cur_t += this->m_par->get_modelpar().get_ctrl_dis();
			}
			cur_t = this->m_par->get_modelpar().get_ctrl_dis() - restdis;//减去剩下的距离
		}
		cas->add(nodes->getAt(nodes->size() - 1));
		newinfos->emplace_back(this->nodeinfos->at(nodes->size() - 1));
	}
	set_center(GeometryFactory::getDefaultInstance()->createLineString(cas));
	set_nodeinfos(newinfos);
	re_computer_curvature();
}

double channel::match_ratio(const well3d_sequence& wells, const wellsfilter& filter) const
{
	double alllenth = 0;
	double macth = 0;
	auto ignore = filter.get_ignoreids();
	auto facies = filter.get_facie_filter();
	for (size_t i = 0; i < wells.get_size(); i++)
	{
		if (std::find(facies->begin(), facies->end(), wells.get_at(i).get_facie()) != facies->end()) {
			if (std::find(ignore->begin(), ignore->end(), wells.get_at(i).get_id()) == ignore->end()) {
				auto sectpt = channelex::get_channel_section(this->center, this->nodeinfos, this->crevasse_infos, this->m_par->get_zone(), wells.get_at(i).get_pos());
				auto sect = sectpt.get();
				if (sect->get_face() == facie::channel_sand) {
					auto ch_sect = (ch_section*)sect;
					auto cross = ch_sect->cross_lenth(wells.get_at(i).get_top(), wells.get_at(i).get_bot());
					if (cross > 0.0001) {
						alllenth += cross;
						if (sect->get_face() == wells.get_at(i).get_facie()) {
							macth += cross;
						}
					}
				}
			}
		}
	}
	return macth < 0.000001 ? 0.0 : macth / alllenth;
}

void channel::well_cross_info(const well3d_sequence& wells, const double& maxdis, double& welllen, double& cross, double& sand, double& mud) const
{
	cross = 0;
	sand = 0;
	mud = 0;
	welllen = 0;
	for (size_t i = 0; i < wells.get_size(); i++)
	{
		Coordinate trans = geo_math::cor_translate(wells.get_at(i).get_pos(), this->m_par->get_zone().get_start_dot(), this->m_par->get_zone().get_angel());
		size_t index;
		Coordinate crossdot;
		auto dis = geo_math::point_min_pos_to_line(trans, this->center, index, crossdot);
		if (dis > maxdis) {
			continue;
		}
		welllen += wells.get_at(i).get_thick();
		auto nodeinfo = nodeinfos->at(index);
		if (index < this->center->getCoordinatesRO()->size() - 1) {
			nodeinfo = channelex::interp_nodeinfo(nodeinfos, this->center, index, crossdot);
		}
		if (dis <= nodeinfo.get_width() / 2) {
			auto ch_z = nodeinfos->at(index).get_z();
			auto cv = nodeinfos->at(index).get_curv();
			bool transisleft = geo_math::cross(cv.get_targent(), Coordinate(trans.x - crossdot.x, trans.y - crossdot.y)) > 0;
			double wid_pos = nodeinfo.get_width() * 0.5 + (transisleft ? -1 : 1) * dis;
			ch_section sect(nodeinfo.get_width(), nodeinfo.get_thick(), cv.get_ay(), ch_z, wid_pos);
			auto inte = sect.cross_lenth(wells.get_at(i).get_top(), wells.get_at(i).get_bot());
			if (wells.get_at(i).get_facie() == facie::channel_sand) {
				sand += inte;
			}
			else {
				mud += inte;
			}
			cross += inte;
		}
	}
}

const size_t channel::get_nodes_count() const
{
	return this->center->getNumPoints();
}

Coordinate channel::get_node_at(size_t index) const
{
	auto orignalpos = geo_math::get_original_cor(Coordinate(this->m_par->get_zone().get_start_dot().x, this->m_par->get_zone().get_start_dot().y), this->m_par->get_zone().get_angel());
	auto trans = geo_math::cor_translate(this->center->getCoordinateN(index), orignalpos, -this->m_par->get_zone().get_angel());
	return trans;
}

nodeinfo channel::get_info_at(size_t index) const
{
	return this->nodeinfos->at(index);
}

void channel::clear()
{
	set_center(NULL);
	if (nodeinfos != NULL) {
		nodeinfos->clear();
	}
	if (crevasse_infos != NULL) {
		crevasse_infos->clear();
	}
}

void channel::set_center(LineString* line)
{
	if (center != NULL) {
		GeometryFactory::getDefaultInstance()->destroyGeometry(center);
		center = NULL;
	}
	if (line != NULL) {
		center = line;
	}
}

void channel::set_nodeinfos(vector<nodeinfo>* nodeinfos)
{
	if (this->nodeinfos != NULL) {
		this->nodeinfos->clear();
		delete this->nodeinfos;
	}
	this->nodeinfos = nodeinfos;
}

double channel::get_wells_dis(const well3d_sequence* all_wells, const double& maxdis) const
{
	double acc_dis = 0;
	//double maxwidth = 0;
	/*vector<int>* obeys = new vector<int>();
	auto cors = center->getCoordinatesRO();
	for (size_t i = 0; i < cors->size(); i++)
	{
		maxwidth = max(maxwidth, this->nodeinfos->at(i).get_width());
	}
	auto pt = this->center->buffer(maxdis + maxwidth / 2, 8, 2);
	auto polygon = (geos::geom::Polygon*)pt.get();*/
	if (all_wells->get_size() > 0) {
		auto factory = GeometryFactory::getDefaultInstance();
		vector<trans_well> trans_wells;
		for (size_t i = 0; i < all_wells->get_size(); i++) {
			auto getfacie = get_channel_facie(all_wells->get_at(i).get_pos());
			if (all_wells->get_at(i).get_facie() == facie::channel_sand) {
				if (getfacie == facie::channel_sand) {//砂岩匹配
					acc_dis -= 1;
				}
				else {//砂岩在外

				}
			}
			else {
				if (getfacie == facie::channel_sand) {//泥岩在内
					acc_dis += 2;
				}
				else {//泥岩匹配

				}
			}
		}
	}
	return acc_dis;
}

unique_ptr<CoordinateArraySequence> channel::get_centerline() const
{
	CoordinateArraySequence* res = new CoordinateArraySequence();
	auto cors = center->getCoordinatesRO();
	auto orignalpos = geo_math::get_original_cor(Coordinate(this->m_par->get_zone().get_start_dot().x, this->m_par->get_zone().get_start_dot().y), this->m_par->get_zone().get_angel());
	for (size_t i = 0; i < cors->size(); i++)
	{
		auto& pos = cors->getAt(i);
		auto trans = geo_math::cor_translate(cors->getAt(i), orignalpos, -this->m_par->get_zone().get_angel());
		res->add(trans);
	}
	return unique_ptr<CoordinateArraySequence>(res);
}

unique_ptr<CoordinateArraySequence> channel::buffer_centerline(const double& distance) const
{
	CoordinateArraySequence* res = new CoordinateArraySequence();
	auto pt = this->center->buffer(distance, 8, 2);
	auto polygon = (geos::geom::Polygon*)pt.get();
	auto dotspt = polygon->getCoordinates();
	auto dots = dotspt.get();
	auto orignalpos = geo_math::get_original_cor(Coordinate(this->m_par->get_zone().get_start_dot().x, this->m_par->get_zone().get_start_dot().y), this->m_par->get_zone().get_angel());
	for (size_t i = 0; i < dots->size(); i++)
	{
		auto pos = dots->getAt(i);
		auto trans = geo_math::cor_translate(dots->getAt(i), orignalpos, -this->m_par->get_zone().get_angel());
		res->add(trans);
	}
	return unique_ptr<CoordinateArraySequence>(res);
}

unique_ptr<vector<curvature>> channel::get_curvatures() const
{
	vector<curvature>* curs = new vector<curvature>();
	curs->reserve(this->nodeinfos->size());
	for (size_t i = 0; i < this->nodeinfos->size(); i++)
	{
		curs->emplace_back(this->nodeinfos->at(i).get_curv());
	}
	return unique_ptr<vector<curvature>>(curs);
}

const LineString* channel::get_centerline_internal() const
{
	return this->center;
}

const channelpar* channel::get_par() const
{
	return this->m_par;
}

void channel::set_par(const channelpar& m_par)
{
	if (this->m_par != NULL) {
		delete this->m_par;
	}
	this->m_par = new channelpar(m_par);
}

void channel::get_in_channel_pos(grid3d* inunits) const
{
	if (inunits != NULL && this->center != NULL) {
		const rect3d gridRange = inunits->get_range();
		const Coordinate gridSize = inunits->get_size();
		auto factory = GeometryFactory::getDefaultInstance();
		double maxwidth = 0;
		auto cor_pos = Coordinate(this->m_par->get_zone().get_start_dot().x, this->m_par->get_zone().get_start_dot().y);
		auto orignalpos = geo_math::get_original_cor(cor_pos, this->m_par->get_zone().get_angel());
		auto cors = center->getCoordinatesRO();
		CoordinateArraySequence* cas = new CoordinateArraySequence();
		for (size_t i = 0; i < cors->size(); i++)
		{
			auto trans = geo_math::cor_translate(cors->getAt(i), orignalpos, -this->m_par->get_zone().get_angel());
			cas->add(trans);
			maxwidth = max(maxwidth, this->nodeinfos->at(i).get_width());
		}
		auto originalcenter = factory->createLineString(cas);
		auto pt = originalcenter->buffer(maxwidth, 8, 2);
		factory->destroyGeometry(originalcenter);
		auto polygon = (geos::geom::Polygon*)pt.get();
		auto rect = polygon->getEnvelopeInternal();
		auto xmin = gridRange.get_xmin() + max(0, gridSize.x * floor((rect->getMinX() - gridRange.get_xmin()) / gridSize.x));
		auto xmax = gridRange.get_xmax() - max(0, gridSize.x * floor((gridRange.get_xmax() - rect->getMaxX()) / gridSize.x));
		auto ymin = gridRange.get_ymin() + max(0, gridSize.y * floor((rect->getMinY() - gridRange.get_ymin()) / gridSize.y));
		auto ymax = gridRange.get_ymax() - max(0, gridSize.y * floor((gridRange.get_ymax() - rect->getMaxY()) / gridSize.y));
		vector<Coordinate> indots;
		for (double y = ymin; y <= ymax; y += gridSize.y)
		{
			CoordinateArraySequence* segdots = new CoordinateArraySequence();
			segdots->add(Coordinate(xmin, y));
			segdots->add(Coordinate(xmax, y));
			auto segline = factory->createLineString(segdots);
			auto temppt = polygon->intersection(segline);
			factory->destroyGeometry(segline);
			auto cross = temppt.get();
			if (cross != NULL && !cross->isEmpty()) {
				set<double> x;
				if (cross->getGeometryTypeId() == GeometryTypeId::GEOS_LINESTRING) {
					auto crosssequence = ((LineString*)cross)->getCoordinatesRO();
					for (size_t crossindex = 0; crossindex < crosssequence->size(); crossindex++)
					{
						x.insert(crosssequence->getAt(crossindex).x);
					}
				}
				else if (cross->getGeometryTypeId() == GeometryTypeId::GEOS_MULTILINESTRING) {
					MultiLineString* mul_line = (MultiLineString*)cross;
					for (size_t k = 0; k < mul_line->getNumGeometries(); k++)
					{
						auto crosssequence = ((LineString*)mul_line->getGeometryN(k))->getCoordinatesRO();
						for (size_t crossindex = 0; crossindex < crosssequence->size(); crossindex++)
						{
							x.insert(crosssequence->getAt(crossindex).x);
						}
					}
				}
				if (x.size() >= 2) {
					vector<double> sortedx;
					sortedx.insert(sortedx.begin(), x.begin(), x.end());
					geo_math::inpolygon_dots(polygon, sortedx, gridRange.get_xmin(), y, gridSize.x, indots);
				}
			}
		}
		for (size_t i = 0; i < indots.size(); i++)
		{
			auto trans = geo_math::cor_translate(indots[i], cor_pos, this->m_par->get_zone().get_angel());
			auto sectpt = channelex::get_channel_section_inzone(this->center, this->nodeinfos, this->crevasse_infos, this->m_par->get_zone(), Coordinate(trans.x, trans.y));
			auto sect = sectpt.get();
			if (sect->get_face() != facie::_default) {
				for (double z = gridRange.get_zmin(); z <= gridRange.get_zmax(); z += gridSize.z)
				{
					if (sect->z_is_in(z)) {
						inunits->add_unit(Coordinate(indots[i].x, indots[i].y, z), (int)sect->get_face());
					}
				}
			}
		}
	}
}

void channel::get_in_channel_pos(grid2d* inunits) const
{
	if (inunits != NULL && this->center != NULL) {
		const Envelope gridRange = inunits->get_range();
		const Coordinate gridSize = inunits->get_size();
		auto factory = GeometryFactory::getDefaultInstance();
		double maxwidth = 0;
		auto cor_pos = Coordinate(this->m_par->get_zone().get_start_dot().x, this->m_par->get_zone().get_start_dot().y);
		auto orignalpos = geo_math::get_original_cor(cor_pos, this->m_par->get_zone().get_angel());
		auto cors = center->getCoordinatesRO();
		CoordinateArraySequence* cas = new CoordinateArraySequence();
		for (size_t i = 0; i < cors->size(); i++)
		{
			auto trans = geo_math::cor_translate(cors->getAt(i), orignalpos, -this->m_par->get_zone().get_angel());
			cas->add(trans);
			maxwidth = max(maxwidth, this->nodeinfos->at(i).get_width());
		}
		auto originalcenter = factory->createLineString(cas);
		auto pt = originalcenter->buffer(maxwidth, 8, 1);
		factory->destroyGeometry(originalcenter);
		auto polygon = (geos::geom::Polygon*)pt.get();
		auto rect = polygon->getEnvelopeInternal();

		auto xmin = gridRange.getMinX() + max(0, gridSize.x * floor((rect->getMinX() - gridRange.getMinX()) / gridSize.x));
		auto xmax = gridRange.getMaxX() - max(0, gridSize.x * floor((gridRange.getMaxX() - rect->getMaxX()) / gridSize.x));
		auto ymin = gridRange.getMinY() + max(0, gridSize.y * floor((rect->getMinY() - gridRange.getMinY()) / gridSize.y));
		auto ymax = gridRange.getMaxY() - max(0, gridSize.y * floor((gridRange.getMaxY() - rect->getMaxY()) / gridSize.y));
		vector<Coordinate> indots;
		for (double y = ymin; y <= ymax; y += gridSize.y)
		{
			CoordinateArraySequence* segdots = new CoordinateArraySequence();
			segdots->add(Coordinate(xmin, y));
			segdots->add(Coordinate(xmax, y));
			auto segline = factory->createLineString(segdots);
			auto temppt = polygon->intersection(segline);
			factory->destroyGeometry(segline);
			auto cross = temppt.get();
			if (cross != NULL && !cross->isEmpty()) {
				set<double> x;
				if (cross->getGeometryTypeId() == GeometryTypeId::GEOS_LINESTRING) {
					auto crosssequence = ((LineString*)cross)->getCoordinatesRO();
					for (size_t crossindex = 0; crossindex < crosssequence->size(); crossindex++)
					{
						x.insert(crosssequence->getAt(crossindex).x);
					}
				}
				else if (cross->getGeometryTypeId() == GeometryTypeId::GEOS_MULTILINESTRING) {
					MultiLineString* mul_line = (MultiLineString*)cross;
					for (size_t k = 0; k < mul_line->getNumGeometries(); k++)
					{
						auto crosssequence = ((LineString*)mul_line->getGeometryN(k))->getCoordinatesRO();
						for (size_t crossindex = 0; crossindex < crosssequence->size(); crossindex++)
						{
							x.insert(crosssequence->getAt(crossindex).x);
						}
					}
				}
				if (x.size() >= 2) {
					vector<double> sortedx;
					sortedx.insert(sortedx.begin(), x.begin(), x.end());
					geo_math::inpolygon_dots(polygon, sortedx, gridRange.getMinX(), y, gridSize.x, indots);
				}
			}
		}
		for (size_t i = 0; i < indots.size(); i++)
		{
			auto trans = geo_math::cor_translate(indots[i], cor_pos, this->m_par->get_zone().get_angel());
			auto sectpt = channelex::get_channel_section_inzone(this->center, this->nodeinfos, this->crevasse_infos, this->m_par->get_zone(), Coordinate(trans.x, trans.y));
			auto sect = sectpt.get();
			if (sect->get_face() != facie::_default) {
				inunits->add_unit(Coordinate(indots[i].x, indots[i].y), (int)sect->get_face());
			}
		}
	}
}

void channel::get_migrate_in_channel_pos(grid3d* inunits) const
{
	if (inunits != NULL && this->center != NULL) {
		const rect3d gridRange = inunits->get_range();
		const Coordinate gridSize = inunits->get_size();
		auto factory = GeometryFactory::getDefaultInstance();
		double maxwidth = 0;
		auto cor_pos = Coordinate(this->m_par->get_zone().get_start_dot().x, this->m_par->get_zone().get_start_dot().y);
		auto orignalpos = geo_math::get_original_cor(cor_pos, this->m_par->get_zone().get_angel());
		auto cors = center->getCoordinatesRO();
		CoordinateArraySequence* cas = new CoordinateArraySequence();
		for (size_t i = 0; i < cors->size(); i++)
		{
			auto trans = geo_math::cor_translate(cors->getAt(i), orignalpos, -this->m_par->get_zone().get_angel());
			cas->add(trans);
			maxwidth = max(maxwidth, this->nodeinfos->at(i).get_width());
		}
		auto originalcenter = factory->createLineString(cas);
		auto pt = originalcenter->buffer(maxwidth, 8, 2);
		factory->destroyGeometry(originalcenter);
		auto polygon = (geos::geom::Polygon*)pt.get();
		auto rect = polygon->getEnvelopeInternal();
		auto xmin = gridRange.get_xmin() + max(0, gridSize.x * floor((rect->getMinX() - gridRange.get_xmin()) / gridSize.x));
		auto xmax = gridRange.get_xmax() - max(0, gridSize.x * floor((gridRange.get_xmax() - rect->getMaxX()) / gridSize.x));
		auto ymin = gridRange.get_ymin() + max(0, gridSize.y * floor((rect->getMinY() - gridRange.get_ymin()) / gridSize.y));
		auto ymax = gridRange.get_ymax() - max(0, gridSize.y * floor((gridRange.get_ymax() - rect->getMaxY()) / gridSize.y));
		vector<Coordinate> indots;
		for (double y = ymin; y <= ymax; y += gridSize.y)
		{
			CoordinateArraySequence* segdots = new CoordinateArraySequence();
			segdots->add(Coordinate(xmin, y));
			segdots->add(Coordinate(xmax, y));
			auto segline = factory->createLineString(segdots);
			auto temppt = polygon->intersection(segline);
			factory->destroyGeometry(segline);
			auto cross = temppt.get();
			if (cross != NULL && !cross->isEmpty()) {
				set<double> x;
				if (cross->getGeometryTypeId() == GeometryTypeId::GEOS_LINESTRING) {
					auto crosssequence = ((LineString*)cross)->getCoordinatesRO();
					for (size_t crossindex = 0; crossindex < crosssequence->size(); crossindex++)
					{
						x.insert(crosssequence->getAt(crossindex).x);
					}
				}
				else if (cross->getGeometryTypeId() == GeometryTypeId::GEOS_MULTILINESTRING) {
					MultiLineString* mul_line = (MultiLineString*)cross;
					for (size_t k = 0; k < mul_line->getNumGeometries(); k++)
					{
						auto crosssequence = ((LineString*)mul_line->getGeometryN(k))->getCoordinatesRO();
						for (size_t crossindex = 0; crossindex < crosssequence->size(); crossindex++)
						{
							x.insert(crosssequence->getAt(crossindex).x);
						}
					}
				}
				if (x.size() >= 2) {
					vector<double> sortedx;
					sortedx.insert(sortedx.begin(), x.begin(), x.end());
					geo_math::inpolygon_dots(polygon, sortedx, gridRange.get_xmin(), y, gridSize.x, indots);
				}
			}
		}
		for (size_t i = 0; i < indots.size(); i++)
		{
			auto trans = geo_math::cor_translate(indots[i], cor_pos, this->m_par->get_zone().get_angel());

			size_t index;
			Coordinate cross;
			auto dis = geo_math::point_min_pos_to_line(trans, this->center, index, cross);
			auto nodeinfo = nodeinfos->at(index);
			if (index < this->center->getCoordinatesRO()->size() - 1) {
				nodeinfo = channelex::interp_nodeinfo(nodeinfos, this->center, index, cross);
			}
			auto width = nodeinfo.get_width();
			auto thick = nodeinfo.get_thick();
			auto ch_z = nodeinfo.get_z();
			auto cv = nodeinfo.get_curv();
			if (dis <= width / 2) {

				bool transisleft = geo_math::cross(cv.get_targent(), Coordinate(trans.x - cross.x, trans.y - cross.y)) > 0;
				double wid_pos = width * 0.5 + (transisleft ? -1 : 1) * dis;
				ch_section sect(width, thick, cv.get_ay(), ch_z, wid_pos);
				double sect_thick = sect.get_thick();
				double z_start = min(ch_z - sect_thick, ch_z);
				{
					inunits->add_unit(Coordinate(indots[i].x, indots[i].y, z_start - gridSize.z), (int)facie::migrate_mud);
				}
				for (double z = z_start; z <= ch_z; z += gridSize.z)
				{
					inunits->add_unit(Coordinate(indots[i].x, indots[i].y, z), (int)facie::channel_sand);
				}
			}
		}
	}
}

facie channel::get_channel_facie(Coordinate pos)const
{
	auto sectpt = channelex::get_channel_section(this->center, this->nodeinfos, this->crevasse_infos, this->m_par->get_zone(), pos);
	auto sect = sectpt.get();
	return sect->get_face();
}

void channel::add_crevasse_pos(const Coordinate& pos, double radius)
{
	auto trans = geo_math::cor_translate(pos, this->m_par->get_zone().get_start_dot(), this->m_par->get_zone().get_angel());
	size_t index;
	Coordinate cross;
	geo_math::point_min_pos_to_line(trans, this->center, index, cross);
	auto nodeinfo = channelex::interp_nodeinfo(this->nodeinfos, this->center, index, cross);
	crevasse_info info(cross, nodeinfo.get_curv().get_dir(), radius);
	this->crevasse_infos->push_back(info);
}

unique_ptr<vector<well3d>> channel::in_ch_segment(const well3d_sequence& wells, const wellsfilter& filter) const
{
	vector<well3d>* inwell = new vector<well3d>();
	auto ignore = filter.get_ignoreids();
	auto facies = filter.get_facie_filter();
	for (size_t i = 0; i < wells.get_size(); i++)
	{
		if (std::find(facies->begin(), facies->end(), wells.get_at(i).get_facie()) != facies->end()) {
			if (std::find(ignore->begin(), ignore->end(), wells.get_at(i).get_id()) == ignore->end()) {
				auto sectpt = channelex::get_channel_section(this->center, this->nodeinfos, this->crevasse_infos, this->m_par->get_zone(), wells.get_at(i).get_pos());
				auto sect = sectpt.get();
				if (sect->get_face() == facie::channel_sand) {
					auto ch_sect = (ch_section*)sect;
					auto top = min(wells.get_at(i).get_top(), ch_sect->get_z());
					auto bot = max(wells.get_at(i).get_bot(), ch_sect->get_z() - ch_sect->get_thick());
					if (top > bot) {
						inwell->push_back(well3d(wells.get_at(i).get_id(), wells.get_at(i).get_pos(),
							wells.get_at(i).get_facie(), top, bot));
					}
				}
			}
		}
	}
	return unique_ptr<vector<well3d>>(inwell);
}

unique_ptr<vector<int>> channel::get_obey_wells(const well3d_sequence* all_wells) const {
	vector<int>* obeys = new vector<int>();
	if (all_wells->get_size() > 0) {
		auto factory = GeometryFactory::getDefaultInstance();
		vector<trans_well> trans_wells;
		for (size_t i = 0; i < all_wells->get_size(); i++) {
			Coordinate dot = geo_math::cor_translate(all_wells->get_at(i).get_pos(), this->m_par->get_zone().get_start_dot(), this->m_par->get_zone().get_angel());
			trans_wells.emplace_back(all_wells->get_at(i), dot);
		}
		for (size_t i = 0; i < trans_wells.size(); i++) {
			auto pt = channelex::get_channel_section_inzone(this->center, this->nodeinfos, this->crevasse_infos, this->m_par->get_zone(), trans_wells[i].get_pos());
			auto sect = pt.get();
			if (sect->get_face() == trans_wells[i].get_facie()) {
				if (trans_wells[i].get_facie() == facie::channel_sand) {
					auto subs = ((ch_section*)sect)->get_thick() - trans_wells[i].get_thick();
					if (abs(subs) <= 0.5 * this->m_par->get_modelpar().get_ztol()) {
						obeys->push_back(trans_wells[i].get_wellid());
					}
					else {
						subs = subs;
					}
				}
				else {
					obeys->push_back(trans_wells[i].get_wellid());
				}
			}
		}
	}
	return unique_ptr<vector<int>>(obeys);
}

bool channel::is_well_match(const well3d& well) const
{
	Coordinate trans = geo_math::cor_translate(well.get_pos(), this->m_par->get_zone().get_start_dot(), this->m_par->get_zone().get_angel());
	size_t index;
	Coordinate cross;
	auto dis = geo_math::point_min_pos_to_line(trans, this->center, index, cross);

	if (index == 0 || index >= this->center->getCoordinatesRO()->size() - 1) {//首尾不算
		return false;
	}

	auto nodeinfo = nodeinfos->at(index);
	if (index < this->center->getCoordinatesRO()->size() - 1) {
		nodeinfo = channelex::interp_nodeinfo(nodeinfos, this->center, index, cross);
	}
	if (dis <= nodeinfo.get_width() / 2) {//井在河道内
		if (well.get_facie() == facie::channel_sand) {//砂岩井
			auto ch_z = nodeinfos->at(index).get_z();
			auto& cv = nodeinfos->at(index).get_curv();
			bool transisleft = geo_math::cross(cv.get_targent(), Coordinate(trans.x - cross.x, trans.y - cross.y)) > 0;
			double wid_pos = nodeinfo.get_width() * 0.5 + (transisleft ? -1 : 1) * dis;
			ch_section sect(nodeinfo.get_width(), nodeinfo.get_thick(), cv.get_ay(), ch_z, wid_pos);
			auto gt = sect.get_thick() - well.get_thick();
			return gt >= 0;
		}
		else {
			return false;
		}
	}
	else {
		if (well.get_facie() == facie::channel_sand)
		{
			return false;
		}
		else {
			return true;
		}
	}
}

bool channel::is_conflict_with_nosand_wells(const vector<trans_well>* nosand_obey_wells)
{
	bool conflict = false;
	for (size_t i = 0; i < nosand_obey_wells->size(); i++)
	{
		auto pt = channelex::get_channel_section_inzone(this->center, this->nodeinfos, this->crevasse_infos, this->m_par->get_zone(), nosand_obey_wells->at(i).get_pos());
		auto facie = pt.get();
		if (facie->get_face() != nosand_obey_wells->at(i).get_facie()) {
			conflict = true;
			break;
		}
	}
	return conflict;
}

void channel::adjust_thick_match_sand(const vector<trans_well>* sandwells, const bool& bigchange)
{
	vector<node_change_info> thick_info;
	vector<int> delids;
	thick_info.emplace_back(0, this->nodeinfos->at(0).get_thick());
	for (size_t i = 0; i < sandwells->size(); i++)
	{
		size_t index;
		Coordinate cross;
		auto trans = sandwells->at(i).get_pos();
		auto dis = geo_math::point_min_pos_to_line(trans, this->center, index, cross);
		auto width = nodeinfos->at(index).get_width();
		if (index < this->center->getCoordinatesRO()->size() - 1) {
			width = channelex::interp_width(nodeinfos, this->center, index, cross);
		}
		auto thick = nodeinfos->at(index).get_thick();
		if (dis <= width / 2) {
			auto cv = nodeinfos->at(index).get_curv();
			auto ch_z = nodeinfos->at(index).get_z();
			auto ay = cv.get_ay();
			bool transisleft = geo_math::cross(cv.get_targent(), Coordinate(trans.x - cross.x, trans.y - cross.y)) > 0;
			double wid_pos = width * 0.5 + (transisleft ? -1 : 1) * dis;
			ch_section sect(width, thick, ay, ch_z, wid_pos);
			auto pos_thick = sect.get_thick();
			if (pos_thick > 0) {
				auto thick_k = sandwells->at(i).get_thick() / pos_thick;//厚度变换比例
				if (bigchange && dis < (0.5 * this->m_par->get_modelpar().get_xytol())) {
					thick_info.emplace_back(index, thick_k * thick);
				}
				else {
					if (thick_k < 2 && thick > 0.5) {
						thick_info.emplace_back(index, thick_k * thick);
					}
				}
			}
		}
	}
	thick_info.emplace_back(this->nodeinfos->size() - 1, this->nodeinfos->at(this->nodeinfos->size() - 1).get_thick());
	sort(thick_info.begin(), thick_info.end(), [=](const node_change_info& info1, const node_change_info& info2) -> bool {return info1.line_index < info2.line_index; });
	for (size_t i = 1; i < thick_info.size(); i++)
	{
		size_t cur = thick_info[i].line_index;
		size_t last = thick_info[i - 1].line_index;
		if (cur == last) {
			continue;
		}
		auto k = (thick_info[i].newvalue - thick_info[i - 1].newvalue) / (cur - last);
		for (size_t index = last; index < cur; index++)
		{
			this->nodeinfos->at(index).set_thick(thick_info[i - 1].newvalue + (index - last) * k);
		}
	}
}

void channel::adjust_ch_z(const vector<trans_well>* sandwells)
{
	vector<node_change_info> z_info;
	vector<int> delids;
	z_info.emplace_back(0, this->nodeinfos->at(0).get_thick());
	for (size_t i = 0; i < sandwells->size(); i++)
	{
		size_t index;
		Coordinate cross;
		auto trans = sandwells->at(i).get_pos();
		auto dis = geo_math::point_min_pos_to_line(trans, this->center, index, cross);
		auto width = nodeinfos->at(index).get_width();
		if (index < this->center->getCoordinatesRO()->size() - 1) {
			width = channelex::interp_width(nodeinfos, this->center, index, cross);
		}
		auto thick = nodeinfos->at(index).get_thick();
		if (dis <= width / 2) {
			auto cv = nodeinfos->at(index).get_curv();
			auto ch_z = nodeinfos->at(index).get_z();
			z_info.emplace_back(index, sandwells->at(i).get_top());
		}
	}
	z_info.emplace_back(this->nodeinfos->size() - 1, this->nodeinfos->at(this->nodeinfos->size() - 1).get_thick());
	sort(z_info.begin(), z_info.end(), [=](const node_change_info& info1, const node_change_info& info2) -> bool {return info1.line_index < info2.line_index; });
	for (size_t i = 1; i < z_info.size(); i++)
	{
		size_t cur = z_info[i].line_index;
		size_t last = z_info[i - 1].line_index;
		if (cur == last) {
			continue;
		}
		auto k = (z_info[i].newvalue - z_info[i - 1].newvalue) / (cur - last);
		for (size_t index = last; index < cur; index++)
		{
			this->nodeinfos->at(index).set_z(z_info[i - 1].newvalue + (index - last) * k);
		}
	}
}

void channel::re_computer_curvature()
{
	auto pts = center->getCoordinatesRO();
	//计算曲率
	double leftmaxcv = 0;
	double rightmaxcv = 0;
	auto count = center->getNumPoints();
	if (count > 2) {
		for (size_t i = 1; i < count - 1; i++)
		{
			auto dot1 = pts->getAt(i - 1);
			auto dot2 = pts->getAt(i);
			auto dot3 = pts->getAt(i + 1);
			Coordinate pre(dot2.x - dot1.x, dot2.y - dot1.y);
			Coordinate last(dot3.x - dot2.x, dot3.y - dot2.y);
			double prelen = pre.distance(Coordinate(0, 0));
			double lastlen = last.distance(Coordinate(0, 0));
			auto theta = asin((pre.x * last.y - pre.y * last.x) / (prelen * lastlen));
			int mul = 1;//曲率方向向左
			if (theta < 0)
			{
				mul = -1;//曲率方向向右
			}
			auto dir = geo_math::rotate(pre, mul * 0.5 * (PI - 0.5 * theta));
			geo_math::normal(dir);
			theta = theta / prelen;
			geo_math::normal(pre);
			this->nodeinfos->at(i).set_curv(curvature(abs(theta), dir, pre));
			if (theta < 0) {
				rightmaxcv = max(abs(theta), rightmaxcv);
			}
			else {
				leftmaxcv = max(abs(theta), leftmaxcv);
			}
		}
		if (leftmaxcv > 0.0001 && rightmaxcv > 0.0001) {
			this->nodeinfos->at(0).set_curv(this->nodeinfos->at(1).get_curv());
			this->nodeinfos->at(center->getNumPoints() - 1).set_curv(this->nodeinfos->at(center->getNumPoints() - 2).get_curv());
			for (size_t i = 0; i < center->getNumPoints(); i++)
			{
				if (this->nodeinfos->at(i).get_curv().isleft()) {
					auto newvalue = log((this->nodeinfos->at(i).get_curv().get_p() / leftmaxcv) * (exp(1) - 1) + 1);
					if (isnan(newvalue)) {
						newvalue = 0;
					}
					this->nodeinfos->at(i).get_change_curv().set_p(newvalue);
				}
				else {
					auto newvalue = log((this->nodeinfos->at(i).get_change_curv().get_p() / rightmaxcv) * (exp(1) - 1) + 1);
					if (isnan(newvalue)) {
						newvalue = 0;
					}
					this->nodeinfos->at(i).get_change_curv().set_p(newvalue);
				}
			}
		}
	}
	else {
		throw exception("节点数不够");
	}
}

void channel::adjust_thick_match_sand(const well3d_sequence* need_obey_wells)
{
	auto factory = GeometryFactory::getDefaultInstance();
	auto max_ele = max_element(this->nodeinfos->begin(), this->nodeinfos->end(), [](const nodeinfo& a, const nodeinfo& b) -> bool { return a.get_width() < b.get_width(); });
	auto maxwidth = (*max_ele).get_width();
	auto pt = this->center->buffer(maxwidth);
	auto polygon = (geos::geom::Polygon*)(pt.get());

	vector<trans_well> sandwells;
	for (size_t i = 0; i < need_obey_wells->get_size(); i++) {
		Coordinate dot = geo_math::cor_translate(need_obey_wells->get_at(i).get_pos(), this->m_par->get_zone().get_start_dot(), this->m_par->get_zone().get_angel());
		auto temppoint = factory->createPoint(Coordinate(dot.x, dot.y));
		if (polygon->contains(temppoint)) {
			if (need_obey_wells->get_at(i).get_facie() == facie::channel_sand) {
				sandwells.emplace_back(need_obey_wells->get_at(i), dot);
			}
		}
		factory->destroyGeometry(temppoint);
	}
	adjust_thick_match_sand(&sandwells, true);
}

bool channel::is_sandwell_incenter(const vector<trans_well>* sand_wells)
{
	for (size_t i = 0; i < sand_wells->size(); i++)
	{
		size_t index;
		Coordinate cross;
		auto dis = geo_math::point_min_pos_to_line(sand_wells->at(i).get_pos(), this->center, index, cross);
		if (dis > 0.5 * this->m_par->get_modelpar().get_xytol()) {
			return false;
		}
	}
	return true;
}

void channel::adjust_width_match_nosand(const vector<trans_well>* nosand_wells)
{
	size_t shift = 10;
	for (size_t i = 0; i < nosand_wells->size(); i++)
	{
		if (nosand_wells->at(i).get_facie() != facie::channel_sand)
		{
			size_t index;
			Coordinate cross;
			auto dis = geo_math::point_min_pos_to_line(nosand_wells->at(i).get_pos(), this->center, index, cross);
			auto width = nodeinfos->at(index).get_width();
			if (dis <= (width / 2)) {
				size_t minindex = max(0, index - shift);
				size_t maxindex = min(nodeinfos->size() - 1, index + shift);
				auto maxchange = (width / 2 - dis) * 2.5;
				for (size_t j = minindex; j <= maxindex; j++)
				{
					auto k = ((double)shift - abs((double)j - (double)index)) / ((double)shift);
					this->nodeinfos->at(j).set_width(this->nodeinfos->at(j).get_width() - maxchange * k);
				}
			}
		}
	}
}