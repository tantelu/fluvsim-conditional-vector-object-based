#include "channelex.h"
#include "pch.h"
#include "bezier.h"
#include <geos\algorithm\Angle.h>
#include "geos\operation\union\UnaryUnionOp.h"
#include "mathhelp.h"
#include "f2c.h"
#include "blaswrap.h"
#include <stdlib.h>
#include <time.h>
#include <map>

extern "C"
{
#include <clapack.h>
}

using namespace std;
using namespace geos::operation::geounion;

curvature::curvature() {
	dir = Coordinate::getNull();
	targent = Coordinate::getNull();
	p = 0;
}

curvature::curvature(const double& p, const Coordinate& dir, const Coordinate& targent) {
	this->p = p;
	this->dir = dir;
	this->targent = targent;
}

const Coordinate& curvature::get_dir() const {
	return this->dir;
}

const Coordinate& curvature::get_targent() const {
	return this->targent;
}

double curvature::get_p() const {
	return p;
}

void curvature::set_p(const double& p) {
	this->p = p;
}

bool curvature::isleft() const {
	return geo_math::cross(targent, dir) > 0;
}

double curvature::get_ay() const {
	auto ay = 0.5 * (1 + (isleft() ? -1 : 1) * p);
	ay = min(0.8, ay);
	ay = max(0.2, ay);
	return ay;
}


nodeinfo::nodeinfo(const double& ch_width, const double& ch_thick, const double& ch_z, const double& levee_b, const double& levee_c, const double& levee_w) {
	this->ch_thick = ch_thick;
	this->ch_width = ch_width;
	this->ch_z = ch_z;
	this->levee_b = levee_b;
	this->levee_c = levee_c;
	this->levee_w = levee_w;
}

nodeinfo::nodeinfo(const double& ch_width, const double& ch_thick, const double& ch_z) :nodeinfo(ch_width, ch_thick, ch_z, ch_thick * 0.2, ch_thick * 0.2, ch_width * 0.0) {

}

nodeinfo::nodeinfo(const double& ch_width, const double& ch_thick, const double& ch_z, const channelpar* m_par) : nodeinfo(ch_width, ch_thick, ch_z)
{
	if (m_par->get_levee_b_ratio() >= 0 && m_par->get_levee_b_ratio() <= 1) {
		this->levee_b = ch_width * m_par->get_levee_b_ratio();
	}
	if (m_par->get_levee_c_ratio() >= 0 && m_par->get_levee_c_ratio() <= 1) {
		this->levee_c = ch_width * m_par->get_levee_c_ratio();
	}
	if (m_par->get_levee_w_ratio() >= 0 && m_par->get_levee_w_ratio() <= 1) {
		this->levee_w = ch_width * m_par->get_levee_w_ratio();
	}
}

nodeinfo nodeinfo::interpolation(const nodeinfo& info1, const nodeinfo& info2, double k) {
	double ch_z = info1.ch_z + (info2.ch_z - info1.ch_z) * k;
	double ch_width = info1.ch_width + (info2.ch_width - info1.ch_width) * k;
	double ch_thick = info1.ch_thick + (info2.ch_thick - info1.ch_thick) * k;
	double levee_b = info1.levee_b + (info2.levee_b - info1.levee_b) * k;;
	double levee_c = info1.levee_c + (info2.levee_c - info1.levee_c) * k;;
	double levee_w = info1.levee_w + (info2.levee_w - info1.levee_w) * k;;
	double p = info1.cv.get_p() + ((info1.cv.isleft() == info2.cv.isleft() ? 1 : -1) * info2.cv.get_p() - info1.cv.get_p()) * k;
	nodeinfo info = nodeinfo(ch_width, ch_thick, ch_z, levee_b, levee_c, levee_w);
	Coordinate tar1 = info1.cv.get_targent();
	Coordinate tar2 = info2.cv.get_targent();
	Coordinate dir1 = info1.cv.get_dir();
	Coordinate dir2 = info2.cv.get_dir();
	double tarx = tar1.x + (tar2.x - tar1.x) * k;
	double tary = tar1.y + (tar2.y - tar1.y) * k;
	double dirx = dir1.x + (dir2.x - dir1.x) * k;
	double diry = dir1.y + (dir2.y - dir1.y) * k;
	auto tardis = sqrt(tarx * tarx + (tary * tary));
	auto dirdis = sqrt(dirx * dirx + (diry * diry));
	if (tardis >= 0.0001 && dirdis > 0.0001)
	{
		tarx = tarx / tardis;
		tary = tary / tardis;
		dirx = dirx / dirdis;
		diry = diry / dirdis;
		info.set_curv(curvature(p, Coordinate(dirx, diry), Coordinate(tarx, tary)));
	}
	else {
		info.set_curv(curvature(0, dir1, tar1));
	}
	return info;
}

double nodeinfo::get_levee_width() const {
	return this->levee_w;
}

double nodeinfo::get_levee_b_thick() const {
	return  this->levee_b;
}

double nodeinfo::get_levee_c_thick() const {
	return  this->levee_c;
}

const curvature& nodeinfo::get_curv() const {
	return cv;
}

curvature& nodeinfo::get_change_curv() {
	return cv;
}

void nodeinfo::set_curv(const curvature& value) {
	this->cv = value;
}

double nodeinfo::get_z() const {
	return this->ch_z;
}

double nodeinfo::get_width() const {
	return this->ch_width;
}

double nodeinfo::get_thick() const {
	return this->ch_thick;
}

void nodeinfo::set_thick(double thick) {
	this->ch_thick = thick;
}

void nodeinfo::set_width(double width) {
	this->ch_width = width;
}

void nodeinfo::set_z(double z) {
	this->ch_z = z;
}

crevasse_info::crevasse_info(const Coordinate& pos, const Coordinate& dir, const double& radius) {
	this->pos = pos;
	this->dir = dir;
	this->radius = radius;
}

const Coordinate& crevasse_info::get_pos() const {
	return pos;
}

const Coordinate& crevasse_info::get_dir() const {
	return dir;
}

double crevasse_info::get_radius() const {
	return radius;
}

trans_well::trans_well(const well3d& well, const Coordinate& tran_pos) {
	this->tran_pos = tran_pos;
	this->well_id = well.get_id();
	this->top = well.get_top();
	this->bot = well.get_bot();
	this->_facie = well.get_facie();
}

const Coordinate& trans_well::get_pos()const {
	return tran_pos;
}

int trans_well::get_wellid() const {
	return well_id;
}


const facie& trans_well::get_facie() const {
	return _facie;
}

ch_node::ch_node(const Coordinate& pos, const double& width, const double& thick) {
	if (width <= 0 || thick <= 0) {
		assert("宽厚必须大于0");
	}
	this->pos = pos;
	this->thick = thick;
	this->width = width;
}

const Coordinate& ch_node::get_pos() const {
	return pos;
}

double ch_node::get_thick() const {
	return thick;
}

double ch_node::get_width() const {
	return width;
}

double channelex::interp_width(const vector<nodeinfo>* infos, const LineString* line, const size_t& index, const Coordinate& cross)
{
	auto factory = GeometryFactory::getDefaultInstance();
	Point* startp = factory->createPoint(line->getCoordinatesRO()->getAt(index));
	Point* endp = factory->createPoint(line->getCoordinatesRO()->getAt(index + 1));
	Point* crossp = factory->createPoint(cross);
	auto k = startp->distance(crossp) / startp->distance(endp);
	factory->destroyGeometry(startp);
	factory->destroyGeometry(endp);
	factory->destroyGeometry(crossp);
	auto width = infos->at(index).get_width() + (infos->at(index + 1).get_width() - infos->at(index).get_width()) * k;
	return width;
}

double channelex::interp_thick(const vector<nodeinfo>* infos, const LineString* line, const size_t& index, const Coordinate& cross)
{
	auto factory = GeometryFactory::getDefaultInstance();
	Point* startp = factory->createPoint(line->getCoordinatesRO()->getAt(index));
	Point* endp = factory->createPoint(line->getCoordinatesRO()->getAt(index + 1));
	Point* crossp = factory->createPoint(cross);
	auto k = startp->distance(crossp) / startp->distance(endp);
	factory->destroyGeometry(startp);
	factory->destroyGeometry(endp);
	factory->destroyGeometry(crossp);
	auto thick = infos->at(index).get_thick() + (infos->at(index + 1).get_thick() - infos->at(index).get_thick()) * k;
	return thick;
}

nodeinfo channelex::interp_nodeinfo(const vector<nodeinfo>* infos, const LineString* line, const size_t& index, const Coordinate& cross)
{
	auto factory = GeometryFactory::getDefaultInstance();
	Point* startp = factory->createPoint(line->getCoordinatesRO()->getAt(index));
	Point* endp = factory->createPoint(line->getCoordinatesRO()->getAt(index + 1));
	Point* crossp = factory->createPoint(cross);
	auto k = startp->distance(crossp) / startp->distance(endp);
	factory->destroyGeometry(startp);
	factory->destroyGeometry(endp);
	factory->destroyGeometry(crossp);
	return nodeinfo::interpolation(infos->at(index), infos->at(index + 1), k);
}

shared_ptr<section> channelex::get_channel_section_inzone(LineString* line, vector<nodeinfo>* nodeinfos, vector<crevasse_info>* crevasse_infos, const channelzone& zone, const Coordinate& trans_pos)
{
	size_t index;
	Coordinate cross;
	auto dis = geo_math::point_min_pos_to_line(trans_pos, line, index, cross);
	auto nodeinfo = nodeinfos->at(index);
	if (index < line->getCoordinatesRO()->size() - 1) {
		nodeinfo = interp_nodeinfo(nodeinfos, line, index, cross);
	}
	auto levee_wid = nodeinfo.get_levee_width();
	if (dis <= nodeinfo.get_width() / 2) {
		auto ch_z = nodeinfos->at(index).get_z();
		auto cv = nodeinfos->at(index).get_curv();
		bool transisleft = geo_math::cross(cv.get_targent(), Coordinate(trans_pos.x - cross.x, trans_pos.y - cross.y)) > 0;
		double wid_pos = nodeinfo.get_width() * 0.5 + (transisleft ? -1 : 1) * dis;
		ch_section* sect = new ch_section(nodeinfo.get_width(), nodeinfo.get_thick(), cv.get_ay(), ch_z, wid_pos);
		return shared_ptr<section>(sect);
	}
	return shared_ptr<section>(new default_section());
}

shared_ptr<section> channelex::get_channel_section(LineString* line, vector<nodeinfo>* nodeinfos, vector<crevasse_info>* crevasse_infos, const channelzone& zone, const Coordinate& pos)
{
	if (line != NULL) {
		Coordinate trans = geo_math::cor_translate(pos, zone.get_start_dot(), zone.get_angel());
		auto pt = get_channel_section_inzone(line, nodeinfos, crevasse_infos, zone, trans);
		auto res = pt;
		return res;
	}
	return shared_ptr<section>(NULL);
}

double channelex::mindis_to_nosand_wells(const Coordinate& pos, const vector<trans_well>* nosand_wells)
{
	double dis = DBL_MAX;
	for (size_t i = 0; i < nosand_wells->size(); i++)
	{
		auto temp = pos.distanceSquared(nosand_wells->at(i).get_pos());
		if (temp < dis) {
			dis = temp;
		}
	}
	return sqrt(dis);
}

unique_ptr<move_near_info> channelex::move_to_point(const Coordinate& point, LineString* line, vector<nodeinfo>* nodeinfos, const double& radius, const double& xy_tol)
{
	if (radius > 0 && line != NULL) {
		vector<Coordinate>* moveDots = new vector<Coordinate>();
		moveDots->reserve(line->getNumPoints() + 1);
		double cur_move = 0;
		size_t index;
		Coordinate cross;
		auto mindis = geo_math::point_min_pos_to_line(point, line, index, cross);
		auto width = nodeinfos->at(index).get_width();
		if (index < line->getNumPoints() - 1) {
			width = interp_width(nodeinfos, line, index, cross);
		}
		if (mindis < xy_tol) {//在线上
			return unique_ptr<move_near_info>(new move_near_info(moveDots, 0, move_near_state::too_near));
		}
		else if (mindis - radius <= -0.01) {
			auto move_radius = 1.6 * radius;
			Coordinate dir1(point.x - cross.x, point.y - cross.y);
			geo_math::normal(dir1);
			Coordinate targetdot(point.x + dir1.x * radius, point.y + dir1.y * radius);
			auto move_dis = xy_tol;
			double k = move_dis / (mindis - move_radius);
			auto cors = line->getCoordinatesRO();
			for (size_t i = 0; i < cors->size(); i++)
			{
				auto _pos = cors->getAt(i);
				Coordinate pos(_pos.x, _pos.y);
				Coordinate dir2(targetdot.x - pos.x, targetdot.y - pos.y);
				geo_math::normal(dir2);
				auto dis = point.distance(pos);
				if (dis < move_radius) {//搜索范围内的才需要移动
					auto needmove = (dis - move_radius) * k;
					pos.x = pos.x + dir2.x * needmove;
					pos.y = pos.y + dir2.y * needmove;
					moveDots->push_back(pos);
				}
				else {
					moveDots->push_back(pos);
				}
			}
			return unique_ptr<move_near_info>(new move_near_info(moveDots, move_dis, move_near_state::move_closer));
		}//在范围内
		else {
			return unique_ptr<move_near_info>(new move_near_info(moveDots, 0, move_near_state::too_far));
		}
	}
	throw exception("参数为空，必须指定半径和河道中线");
}

CoordinateArraySequence* channelex::move_away_point(const Coordinate& point, LineString* line, vector<nodeinfo>* nodeinfos, const double& xy_tol, const double& maxradius)
{
	if (line != NULL) {
		size_t index;
		Coordinate cross;
		auto mindis = geo_math::point_min_pos_to_line(point, line, index, cross);
		auto width = nodeinfos->at(index).get_width();
		if (index < line->getNumPoints() - 1) {
			width = interp_width(nodeinfos, line, index, cross);
		}
		if (mindis > (maxradius / 2.0)) {
			return NULL;
		}
		else {//否则搜索范围内必有点
			double maxmovedis = xy_tol;
			Coordinate dir(cross.x - point.x, cross.y - point.y);
			if (cross.equals2D(point, 0.0001)) {
				dir.x = geo_math::random(0.5, 1);
				dir.y = geo_math::random(0.5, 1);
			}
			geo_math::normal(dir);
			auto cors = line->getCoordinatesRO();
			CoordinateArraySequence* cas = new CoordinateArraySequence();
			for (size_t i = 0; i < cors->size(); i++)
			{
				auto& _pos = cors->getAt(i);
				auto disSquar = _pos.distanceSquared(point);
				if (disSquar < (maxradius * maxradius)) {
					auto dis = sqrt(disSquar);
					auto movedis = (1 - ((dis - mindis) / (maxradius - mindis))) * maxmovedis;
					Coordinate newpos(_pos.x + dir.x * movedis, _pos.y + dir.y * movedis);
					cas->add(newpos);
				}
				else {
					cas->add(_pos);
				}
			}
			return cas;
		}
	}
	return NULL;
}

bool channelex::channel_thick_near(const double& well_thick, const ch_section* section, const double& ratio) {
	auto pos_thick = section->get_thick();
	if (abs(pos_thick - well_thick) <= pos_thick * ratio) {
		return true;
	}
	return false;
}

std::unique_ptr<vector<Coordinate>> channelex::create_ctrl_points(const vector<trans_well>& wells, const channelpar* m_par) {
	double inte = (m_par->get_wave() / 2.0);
	double amp = m_par->get_amp();
	double halfwidth = m_par->get_width() / 2.0;
	double max_x = m_par->get_zone().get_mainlen();
	auto factory = GeometryFactory::getDefaultInstance();
	unique_ptr<Geometry> unionpt;
	if (wells.size() > 0) {
		vector<geos::geom::Polygon*> nosand_polygons;
		nosand_polygons.reserve(wells.size());
		for (size_t i = 0; i < wells.size(); i++)
		{
			if (wells[i].get_facie() != facie::channel_sand) {
				auto temppt = factory->createPoint(wells[i].get_pos());
				auto polygonpt = temppt->buffer(halfwidth);
				auto polygon = polygonpt.get();
				if (polygon->getGeometryTypeId() == GeometryTypeId::GEOS_POLYGON) {
					nosand_polygons.emplace_back(((geos::geom::Polygon*)factory->createGeometry(polygon)));
				}
				factory->destroyGeometry(temppt);
			}
		}
		UnaryUnionOp op(nosand_polygons);
		unionpt = op.Union();
		for (size_t i = 0; i < nosand_polygons.size(); i++)
		{
			factory->destroyGeometry(nosand_polygons[i]);
		}
		nosand_polygons.clear();
	}
	int maxretry = 15;
	int curretry = 0;
	vector<Coordinate>* cas = new vector<Coordinate>();

	double curx = 0;
	while (curx < max_x) {
		curretry = 0;
		Coordinate last = cas->size() == 0 ? Coordinate(-inte, 0) : cas->at(cas->size() - 1);
		while (curretry < maxretry) {
			auto randvalue = rand() / double(RAND_MAX) * 2 - 1;
			double y = randvalue * amp;
			CoordinateArraySequence* templine = new CoordinateArraySequence();
			templine->add(last);
			templine->add(Coordinate(curx, y));
			if (unionpt.get() != NULL) {
				auto templinestring = factory->createLineString(templine);
				auto cross = templinestring->crosses(unionpt.get());
				factory->destroyGeometry(templinestring);
				if (!cross) {
					cas->push_back(Coordinate(curx, y));
					break;
				}
			}
			else {
				cas->push_back(Coordinate(curx, y));
				break;
			}
			curretry++;
		}
		if (curretry == maxretry) {
			auto randvalue = rand() / double(RAND_MAX) * 2 - 1;
			double y = randvalue * amp;
			cas->push_back(Coordinate(curx, y));
		}
		curx += inte;
	}
	return unique_ptr<vector<Coordinate>>(cas);
}

wellsfilter::wellsfilter(const vector<int>& ignoreids, const vector<facie>& fs)
	:ignoreids(new vector<int>(ignoreids)), facies(new vector<facie>(fs))
{

}

const vector<int>* wellsfilter::get_ignoreids() const
{
	return this->ignoreids.get();
}

const vector<facie>* wellsfilter::get_facie_filter() const
{
	return this->facies.get();
}
