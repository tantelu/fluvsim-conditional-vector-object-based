#pragma once
#include "basedefine.h"
#include "curveinterp.h"

using namespace tk;
/// <summary>
/// 切片
/// </summary>
class section {
public:
	virtual bool z_is_in(double z) const = 0;

	virtual facie get_face() const = 0;
};

/// <summary>
/// 河道切片
/// </summary>
class ch_section :public section
{
public:
	/// <summary>
	/// 最大宽度
	/// </summary>
	double max_width;
	/// <summary>
	/// 最大厚度
	/// </summary>
	double max_thick;
	/// <summary>
	/// 最厚点位置，由曲率计算得到（曲率圆在中线右侧 则ay<0.5）
	/// </summary>
	double ay;
	/// <summary>
	/// 切片河面高层
	/// </summary>
	double ch_z;
	/// <summary>
	/// 宽度位置 最左侧为0 最大为max_with  左右和ay的值对应 
	/// </summary>
	double wid_pos;

	ch_section(double max_width, double max_thick, double ay, double ch_z, double wid_pos);

	facie get_face() const override;

	bool z_is_in(double z) const override;

	bool is_in_maxthick(const double& tole) const;

	double get_z() const { return ch_z; }

	double cross_lenth(const double& top, const double& bot) const;

	double get_thick() const;

	unique_ptr<CoordinateArraySequence> get_bezier_lines(const double& step) const;
private:
	double get_thick_by_pos(const double& pos) const;
};

class levee_section :public section
{
private:
	curveinterp* cur;
	double levee_width;
	double levee_c_thick;
	double levee_b_thick;
	double wid_pos;
	double ch_z;
public:

	/// <summary>
	/// 构造
	/// </summary>
	/// <param name="levee_width"></param>
	/// <param name="levee_b_thick"></param>
	/// <param name="levee_c_thick"></param>
	/// <param name="wid_pos">靠近河道为0 最大为天然堤宽度</param>
	/// <param name="ch_z">河道顶面高</param>
	levee_section(double levee_width, double levee_b_thick, double levee_c_thick, double wid_pos, double ch_z);

	~levee_section();

	facie get_face() const override;

	bool z_is_in(double z) const override;
};

class crevasse_section :public section
{
public:

	facie get_face() const override;

	bool z_is_in(double z) const override;
};

class default_section :public section
{
public:
	default_section();

	facie get_face() const override;

	bool z_is_in(double z) const override;
};

