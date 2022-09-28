#include "section.h"
#include "bezier.h"
#include "mathhelp.h"
#include "geos/geos.h"

ch_section::ch_section(double max_width, double max_thick, double ay, double ch_z, double wid_pos)
{
	this->max_width = max_width;
	this->max_thick = max_thick;
	this->ay = ay;
	this->ch_z = ch_z;
	this->wid_pos = wid_pos;
}

facie ch_section::get_face() const
{
	return facie::channel_sand;
}

bool ch_section::z_is_in(double z) const
{
	if (ay <= 0) {
		return false;
	}
	if (z <= ch_z) {
		auto thick = get_thick();
		return ch_z - z < thick;
	}
	return false;
}

bool ch_section::is_in_maxthick(const double& tole) const
{
	return abs(this->wid_pos - (this->ay * this->max_width)) <= tole;
}

double ch_section::cross_lenth(const double& top, const double& bot) const
{
	auto _top = min(top, ch_z);
	auto _bot = max(bot, ch_z - get_thick());
	return max(0.0, _top - _bot);
}

double ch_section::get_thick() const
{
	return this->get_thick_by_pos(this->wid_pos);
}

unique_ptr<CoordinateArraySequence> ch_section::get_bezier_lines(const double& step) const
{
	auto instep = step;
	if (instep > this->max_width / 5.0) {
		instep = this->max_width / 5.0;
	}
	CoordinateArraySequence* cas = new CoordinateArraySequence();
	
	for (double dis = 0; dis < this->max_width; dis += instep)
	{
		auto y = ch_z - this->get_thick_by_pos(dis);
		cas->add(Coordinate(dis - wid_pos, y));
	}
	cas->add(Coordinate(max_width - wid_pos, ch_z));
	cas->add(Coordinate(-wid_pos, ch_z));
	return unique_ptr<CoordinateArraySequence>(cas);
}

double ch_section::get_thick_by_pos(const double& pos) const
{
	if (ay <= 0.5)
	{
		float by = (float)(-log(2) / log(ay));
		float thick = (float)(4 * max_thick * pow(pos / max_width, by) * abs(1 - pow(pos / max_width, by)));//ºñ¶È ÕýÖµ
		return thick;
	}
	else
	{
		float cy = (float)(-log(2) / log(1 - ay));
		float thick = (float)(4 * max_thick * pow((1 - pos / max_width), cy) * abs(1 - pow((1 - pos / max_width), cy)));
		return thick;
	}
}

levee_section::levee_section(double levee_width, double levee_b_thick, double levee_c_thick, double wid_pos, double ch_z)
{
	cur = NULL;
	if (levee_width > 0 && levee_b_thick > 0 && levee_c_thick > 0) {
		vector<Coordinate> ctr = { Coordinate(0,0),Coordinate(levee_width * 0.05,levee_b_thick * 0.56),Coordinate(levee_width * 0.19,levee_b_thick),
		Coordinate(levee_width * 0.31,levee_b_thick * 0.9),Coordinate(levee_width * 0.5,levee_b_thick * 0.61),
		Coordinate(levee_width * 0.75,levee_b_thick * 0.3),Coordinate(levee_width,0) };
		this->cur = new curveinterp(ctr);
		this->levee_c_thick = levee_c_thick;
		this->levee_b_thick = levee_b_thick;
		this->levee_width = levee_width;
		this->ch_z = ch_z;
		this->wid_pos = wid_pos;
	}
	else {
		throw exception(" ");
	}
}

levee_section::~levee_section()
{
	if (cur != NULL) {
		delete cur;
		cur = NULL;
	}
}

facie levee_section::get_face() const
{
	return facie::levee;
}

bool levee_section::z_is_in(double z) const
{
	if (wid_pos > levee_width || wid_pos < 0) {
		return false;
	}
	else {
		if (z - ch_z > levee_b_thick || ch_z - z > levee_c_thick) {
			return false;
		}
		else {
			auto yb = cur->get_dot(wid_pos).y + ch_z;
			auto yc = ch_z - ((wid_pos > 0.8 * levee_width) ? levee_c_thick - ((wid_pos - 0.8 * levee_width) * 5 * levee_c_thick / levee_width) : levee_c_thick);
			if (z >= yc && z <= yb) {
				return true;
			}
			else {
				return false;
			}
		}
	}
}

facie crevasse_section::get_face() const
{
	return facie::crevasse;
}

bool crevasse_section::z_is_in(double z) const
{
	return false;
}

default_section::default_section()
{
}

facie default_section::get_face() const
{
	return facie::_default;
}

bool default_section::z_is_in(double z) const
{
	return false;
}
