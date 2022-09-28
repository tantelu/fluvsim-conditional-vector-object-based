#include "basedefine.h"

rect3d::rect3d(const double& xmin, const double& xmax, const double& ymin, const double& ymax, const double& zmin, const double& zmax)
{
	if (xmax < xmin || ymax < ymin || zmax < zmin) {
		throw std::out_of_range("The scope must be legal");
	}
	this->xmin = xmin;
	this->xmax = xmax;
	this->ymin = ymin;
	this->ymax = ymax;
	this->zmin = zmin;
	this->zmax = zmax;
}

rect3d::rect3d(const rect3d& other) {
	this->xmin = other.xmin;
	this->xmax = other.xmax;
	this->ymin = other.ymin;
	this->ymax = other.ymax;
	this->zmin = other.zmin;
	this->zmax = other.zmax;
}

double rect3d::get_xmin() const { return xmin; }
double rect3d::get_xmax()const { return xmax; }
double rect3d::get_ymin() const { return ymin; }
double rect3d::get_ymax()const { return ymax; }
double rect3d::get_zmin()const { return zmin; }
double rect3d::get_zmax()const { return zmax; }

channelzone::channelzone(const double& angel, const Coordinate& startdot, const double& mainlen)
{
	this->start_dot = startdot;
	this->angel = angel;
	this->mainlen = mainlen;
}

channelzone::channelzone(const channelzone& zone)
{
	this->start_dot = zone.start_dot;
	this->angel = zone.angel;
	this->mainlen = zone.mainlen;
}

double channelzone::get_angel() const
{
	return angel;
}

const Coordinate& channelzone::get_start_dot() const
{
	return start_dot;
}

double channelzone::get_mainlen() const
{
	return mainlen;
}

modelpar::modelpar(const double& xy_tol, const double& z_tol, const double& ctrl_dis, const int& undefine)
{
	this->xy_tol = xy_tol;
	this->z_tol = z_tol;
	this->ctrl_dis = ctrl_dis;
	this->undefinevalue = undefine;
}

modelpar::modelpar(const modelpar& other)
{
	this->xy_tol = other.xy_tol;
	this->z_tol = other.z_tol;
	this->ctrl_dis = other.ctrl_dis;
	this->undefinevalue = other.undefinevalue;
}

double modelpar::get_xytol() const
{
	return xy_tol;
}

double modelpar::get_ztol() const
{
	return z_tol;
}

double modelpar::get_ctrl_dis() const
{
	return ctrl_dis;
}

void modelpar::set_ctrl_dis(const double& dis)
{
	this->ctrl_dis = dis;
}

channelpar::channelpar(const channelzone& chzone, const modelpar& modelpar, const double& width, const double& thick, const double& wave, const double& amp) :
	zone(chzone), m_par(modelpar)
{
	this->width = width;
	this->thick = thick;
	this->amp = amp;
	this->wave = wave;
	this->levee_b_ratio = 0;
	this->levee_c_ratio = 0;
	this->levee_w_ratio = 0;
	this->m_par.set_ctrl_dis(min(modelpar.get_ctrl_dis(), zone.get_mainlen() / 5));
}

channelpar::channelpar(const channelpar& chpar) :m_par(chpar.m_par), zone(chpar.zone)
{
	this->width = chpar.width;
	this->thick = chpar.thick;
	this->amp = chpar.amp;
	this->wave = chpar.wave;
	this->levee_b_ratio = chpar.levee_b_ratio;
	this->levee_c_ratio = chpar.levee_c_ratio;
	this->levee_w_ratio = chpar.levee_w_ratio;
}

double channelpar::get_width() const
{
	return width;
}

double channelpar::get_thick() const
{
	return thick;
}

double channelpar::get_wave() const
{
	return wave;
}

double channelpar::get_amp() const
{
	return amp;
}

double channelpar::get_levee_w_ratio() const
{
	return levee_w_ratio;
}

double channelpar::get_levee_b_ratio() const
{
	return levee_b_ratio;
}

double channelpar::get_levee_c_ratio() const
{
	return levee_c_ratio;
}

const channelzone& channelpar::get_zone() const
{
	return zone;
}

const modelpar& channelpar::get_modelpar() const
{
	return m_par;
}

well3d::well3d(const int& wellid, const Coordinate& pos, const facie& facie, const double& top, const double& bot)
{
	this->wellid = wellid;
	this->pos = pos;
	this->f = facie;
	this->top = top;
	this->bot = bot;
}

well3d::well3d(const well3d& other)
{
	this->wellid = other.wellid;
	this->pos = other.pos;
	this->f = other.f;
	this->top = other.top;
	this->bot = other.bot;
}

const Coordinate& well3d::get_pos() const
{
	return pos;
}

double well3d::get_top() const
{
	return top;
}

double well3d::get_bot() const
{
	return bot;
}

double well3d::get_thick() const
{
	return top - bot;
}

int well3d::get_id() const
{
	return wellid;
}

facie well3d::get_facie() const
{
	return f;
}

well3d_sequence::well3d_sequence(const vector<well3d>* wells)
{
	vector<well3d>* inner = new vector<well3d>();
	inner->reserve(wells->size());
	inner->insert(inner->begin(), wells->begin(), wells->end());
	sort(inner->begin(), inner->end(), [](const well3d& a, const well3d& b) {return a.get_id() < b.get_id(); });
	if (inner->size() > 20) {
		this->mid = (inner->begin()->get_id() + inner->rbegin()->get_id()) / 2;
	}
	else {
		this->mid = 0;
	}
	inwells = unique_ptr<vector<well3d>>(inner);
}

size_t well3d_sequence::get_size() const
{
	return inwells.get()->size();
}

const well3d& well3d_sequence::get_at(size_t index) const
{
	return inwells.get()->at(index);
}

const well3d& well3d_sequence::get_with_id(const int& wellid) const
{
	auto wells = inwells.get();
	if (wellid <= mid) {

		for (auto it = wells->begin(); it != wells->end(); it++)
		{
			if (it->get_id() == wellid) {
				return *it;
			}
		}
	}
	else {
		for (auto it = wells->rbegin(); it != wells->rend(); it++)
		{
			if (it->get_id() == wellid) {
				return *it;
			}
		}
	}
}

bool well3d_sequence::is_empty() const
{
	return inwells.get()->empty();
}

grid3d::grid3d(const modelpar& m_par, const rect3d& range, const Coordinate& size)
	:inrange(range), insize(size), pt_par(new modelpar(m_par))
{
	xcount = (int)ceil((inrange.get_xmax() - inrange.get_xmin()) / (insize.x));
	ycount = (int)ceil((inrange.get_ymax() - inrange.get_ymin()) / (insize.y));
	zcount = (int)ceil((inrange.get_zmax() - inrange.get_zmin()) / (insize.z));
	auto len = get_gridunits_count();
	auto vs = new int[len];
	int undefine = pt_par.get()->get_undefinevalue();
	for (size_t i = 0; i < len; i++)
	{
		*(vs + i) = undefine;
	}
	facies = unique_ptr<int[]>(vs);
}

double grid3d::get_ratio(int f) const
{
	auto fs = facies.get();
	auto count = get_gridunits_count();
	double acc = 0;
	for (size_t i = 0; i < count; i++)
	{
		if (*(fs + i) == f) {
			acc++;
		}
	}
	return acc / count;
}

size_t grid3d::get_gridunits_count() const
{
	return xcount * ycount * zcount;
}

void grid3d::add_unit(const Coordinate& pos3d, int data)
{
	int x = (int)((pos3d.x - this->inrange.get_xmin()) / this->insize.x);
	int y = (int)((pos3d.y - this->inrange.get_ymin()) / this->insize.y);
	int z = (int)((pos3d.z - this->inrange.get_zmin()) / this->insize.z);
	if (x >= 0 && x < xcount && y >= 0 && z >= 0 && y < ycount && z < zcount)
	{
		int index = z * xcount * ycount + y * xcount + x;
		*(facies.get() + index) = data;
	}
}

gridpoint::gridpoint(double p_x, double p_y, double p_z, void* p_data) :p(p_x, p_y, p_z), data(p_data)
{
}

gridpoint::gridpoint(const Coordinate& p_p, void* p_data) : p(p_p), data(p_data)
{
}

grid2d::grid2d(const Envelope& range, const Coordinate& size) : inrange(range), insize(size)
{
	xcount = (int)ceil((inrange.getMaxX() - inrange.getMinX()) / (insize.x));
	ycount = (int)ceil((inrange.getMaxY() - inrange.getMinY()) / (insize.y));
	auto len = get_gridunits_count();
	facies.reserve(len);
	for (size_t i = 0; i < len; i++)
	{
		facies.push_back(0);
	}
}

size_t grid2d::get_gridunits_count() const
{
	return xcount * ycount;
}

void grid2d::add_unit(const Coordinate& pos, int data)
{
	int x = (int)((pos.x - this->inrange.getMinX()) / this->insize.x);
	int y = (int)((pos.y - this->inrange.getMinY()) / this->insize.y);
	if (x >= 0 && x < xcount && y >= 0 && y < ycount)
	{
		int index = y * xcount + x;
		facies[index] = data;
	}
}

void grid2d::add_unit(const size_t& index, int data)
{
	facies[index] = data;
}
