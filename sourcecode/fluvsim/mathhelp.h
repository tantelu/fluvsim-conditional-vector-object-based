#pragma once
#include "basedefine.h"
#include <geos\operation\distance\DistanceOp.h>
#include <geos\linearref\LocationIndexedLine.h>
#include <geos/algorithm/Distance.h>
#include <geos/algorithm/LineIntersector.h>

using namespace geos::algorithm;
using namespace geos::operation::distance;
using namespace geos::linearref;

class geo_math {
public:
	static double random(double min,double max) {
		double number = min + static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / (max - min)));
		return number;
	}

	static Coordinate rotate(const Coordinate& pos, const double& angel) {
		auto x = pos.x * cos(angel) + pos.y * sin(angel);
		auto y = pos.y * cos(angel) - pos.x * sin(angel);
		return Coordinate(x, y);
	}

	static Coordinate get_original_cor(const Coordinate& channel_start, const double& channel_angel) {
		double x = -channel_start.x * cos(channel_angel) - channel_start.y * sin(channel_angel);
		double y = channel_start.x * sin(channel_angel) - channel_start.y * cos(channel_angel);
		return Coordinate(x, y);
	}

	static Coordinate cor_translate(const Coordinate& pos, const Coordinate& cor_pos, const double& angel) {
		auto x = (pos.x - cor_pos.x) * cos(angel) + (pos.y - cor_pos.y) * sin(angel);
		auto y = -(pos.x - cor_pos.x) * sin(angel) + (pos.y - cor_pos.y) * cos(angel);
		return Coordinate(x, y);
	}

	static double cross(const Coordinate& dot1, const Coordinate& dot2) {
		return dot1.x * dot2.y - dot2.x * dot1.y;
	}

	static Coordinate add(const Coordinate& dot1, const Coordinate& dot2) {
		return Coordinate(dot1.x + dot2.x, dot1.y + dot2.y);
	}

	static Coordinate sub(const Coordinate& dot1, const Coordinate& dot2) {
		return Coordinate(dot1.x - dot2.x, dot1.y - dot2.y);
	}

	static bool normal(Coordinate& normal)
	{
		auto dis = sqrt(normal.x * normal.x + (normal.y * normal.y));
		if (dis >= 0.000001)
		{
			normal.x = normal.x / dis;
			normal.y = normal.y / dis;
			return true;
		}
		return false;
	}

	static bool normal(const Coordinate& start, const Coordinate& end, Coordinate& normal)
	{
		auto dis = sqrt((end.x - start.x) * (end.x - start.x) + (end.y - start.y) * (end.y - start.y));
		if (dis >= 0.000001)
		{
			normal.x = (end.x - start.x) / dis;
			normal.y = (end.y - start.y) / dis;
			return true;
		}
		return false;
	}

	static Coordinate normal(const Coordinate& start, const Coordinate& end)
	{
		auto dis = sqrt((end.x - start.x) * (end.x - start.x) + (end.y - start.y) * (end.y - start.y));
		if (dis >= 0.000001)
		{
			return Coordinate((end.x - start.x) / dis, (end.y - start.y) / dis);
		}
		return Coordinate::getNull();
	}

	static unique_ptr<CoordinateArraySequence> intersection(const Coordinate& start, const Coordinate& end, const LineString* line) {
		CoordinateArraySequence* cas = new CoordinateArraySequence();
		auto dots = line->getCoordinatesRO();
		LineIntersector li;
		for (size_t i = 0; i < dots->size() - 1; i++)
		{
			li.computeIntersection(start, end, dots->getAt(i), dots->getAt(i + 1));
			if (li.hasIntersection()) {
				cas->add(li.getIntersection(0));
			}
		}
		return unique_ptr<CoordinateArraySequence>(cas);
	}

	static double point_min_pos_to_line(const Coordinate& point, const LineString* line, size_t& lineindex, Coordinate& cross)
	{
		Coordinate pt(point.x, point.y);
		Point* p = GeometryFactory::getDefaultInstance()->createPoint(pt);
		DistanceOp op(p, line);
		auto dis = op.distance();
		auto dots = op.nearestPoints();
		cross = dots.get()->getAt(1);
		LocationIndexedLine loc(line);
		auto pos = loc.indexOf(cross);
		lineindex = pos.getSegmentIndex();
		GeometryFactory::getDefaultInstance()->destroyGeometry(p);
		return dis;
	}

	static double point_min_pos_to_line(const Coordinate& point, const LineString* line)
	{
		Coordinate pt(point.x, point.y);
		Point* p = GeometryFactory::getDefaultInstance()->createPoint(pt);
		DistanceOp op(p, line);
		auto dis = op.distance();
		GeometryFactory::getDefaultInstance()->destroyGeometry(p);
		return dis;
	}

	static void inpolygon_dots(const geos::geom::Polygon* poly, const vector<double>& sortedx, double oriX, double y, double step, vector<Coordinate>& indots)
	{
		auto factory = GeometryFactory::getDefaultInstance();
		for (int i = 0; i < sortedx.size() - 1; i++)
		{
			Point* p = factory->createPoint(Coordinate((sortedx[i + 1] + sortedx[i]) / 2, y));
			if (poly->contains(p))
			{
				double start = sortedx[i] - oriX;
				double end = sortedx[i + 1] - oriX;
				auto xmin = oriX + step * ceil(sortedx[i] - oriX) / step;
				auto xmax = oriX + step * floor(sortedx[i + 1] - oriX) / step;
				for (double x = xmin; x <= xmax; x += step)
				{
					indots.emplace_back(x, y);
				}
			}
			factory->destroyGeometry(p);
		}
	}

	static bool calculate_control_point(
		const std::vector<Coordinate>& rawPointVector,
		std::vector<Coordinate>& firstControlPointVector,
		std::vector<Coordinate>& secondControlPointVector)
	{
		if (rawPointVector.size() < 2)
		{
			return false;
		}
		std::size_t nPointSize = rawPointVector.size() - 1;
		const Coordinate* pRawPoint = rawPointVector.data();
		if (1 == nPointSize)
		{
			// 3P1 = 2P0 + P3
			firstControlPointVector.resize(1);
			firstControlPointVector[0].x = (2 * pRawPoint[0].x + pRawPoint[1].x) / 3;
			firstControlPointVector[0].y = (2 * pRawPoint[0].y + pRawPoint[1].y) / 3;
			// P2 = 2P1 ¨C P0
			secondControlPointVector.resize(1);
			secondControlPointVector[0].x = 2 * firstControlPointVector[0].x - pRawPoint[0].x;
			secondControlPointVector[0].y = 2 * firstControlPointVector[0].y - pRawPoint[0].y;
			return true;
		}

		std::vector<double> rhs(nPointSize);
		double* pTmp = rhs.data();
		for (std::size_t i = 1; i < nPointSize - 1; ++i)
		{
			pTmp[i] = 4 * pRawPoint[i].x + 2 * pRawPoint[i + 1].x;
		}
		pTmp[0] = pRawPoint[0].x + 2 * pRawPoint[1].x;
		pTmp[nPointSize - 1] = (8 * pRawPoint[nPointSize - 1].x + pRawPoint[nPointSize].x) / 2.0;

		std::vector<double> x;
		get_first_ctrl_points(rhs, x);

		for (std::size_t i = 1; i < nPointSize - 1; ++i)
			pTmp[i] = 4 * pRawPoint[i].y + 2 * pRawPoint[i + 1].y;
		pTmp[0] = pRawPoint[0].y + 2 * pRawPoint[1].y;
		pTmp[nPointSize - 1] = (8 * pRawPoint[nPointSize - 1].y + pRawPoint[nPointSize].y) / 2.0;

		std::vector<double> y;
		get_first_ctrl_points(rhs, y);
		double* pX = x.data();
		double* pY = y.data();

		firstControlPointVector.resize(nPointSize);
		secondControlPointVector.resize(nPointSize);
		Coordinate* pFirstPoints = firstControlPointVector.data();
		Coordinate* pSecondPoints = secondControlPointVector.data();
		for (std::size_t i = 0; i < nPointSize; ++i)
		{
			// Second control point
			pFirstPoints[i].x = x[i];
			pFirstPoints[i].y = y[i];
			if (i < nPointSize - 1)
			{
				pSecondPoints[i].x = 2 * pRawPoint[i + 1].x - pX[i + 1];
				pSecondPoints[i].y = 2 * pRawPoint[i + 1].y - pY[i + 1];
			}
			else
			{
				pSecondPoints[i].x = (pRawPoint[nPointSize].x + pX[nPointSize - 1]) / 2;
				pSecondPoints[i].y = (pRawPoint[nPointSize].y + pY[nPointSize - 1]) / 2;
			}
		}
		return true;
	}

	static vector<int> smallestSufficientTeam(const vector<int>& req_skills, const vector<vector<int>>& people) {
		unordered_map<int, int> skill2idx;
		int skillCnt = req_skills.size();
		for (int i = 0; i < skillCnt; i++) {
			skill2idx[req_skills[i]] = i;
		}
		int peopleNum = people.size();
		vector<long long> people2skill(peopleNum);
		vector<bool> ignore(peopleNum, false);
		for (int i = 0; i < peopleNum; i++) {
			long long skill = 0;
			for (int sk : people[i]) skill += ((long long)1 << skill2idx[sk]);
			people2skill[i] = skill;
		}
		vector<vector<int>> skill2people(skillCnt);
		for (int i = 0; i < peopleNum; i++) {
			for (int j = 0; j < peopleNum; j++) {
				if (i == j) continue;
				if (people2skill[i] != people2skill[j] && (people2skill[i] | people2skill[j]) == people2skill[j]) {
					ignore[i] = true;
					break;
				}
			}
			if (!ignore[i]) {
				for (int j = 0; j < skillCnt; j++) {
					if (people2skill[i] & ((long long)1 << j)) {
						skill2people[j].emplace_back(i);
					}
				}
			}
		}
		auto targetState = (1 << skillCnt) - 1;
		int minGroup = skillCnt + 1;
		vector<int> group;
		vector<int> ans;
		solve(skill2people, people2skill, group, 0, 0, targetState, skillCnt, minGroup, ans);
		return ans;
	}

private:
	static std::vector<double> get_first_ctrl_points(
		const std::vector<double>& rhs, std::vector<double>& x)
	{
		std::size_t n = rhs.size();
		x.resize(n);
		std::vector<double> tmp(n);

		double b = 2.0;
		x[0] = rhs[0] / b;
		for (std::size_t i = 1; i < n; ++i) // Decomposition and forward substitution.
		{
			tmp[i] = 1 / b;
			b = (i < n - 1 ? 4.0 : 3.5) - tmp[i];
			x[i] = (rhs[i] - x[i - 1]) / b;
		}
		for (std::size_t i = 1; i < n; ++i)
			x[n - i - 1] -= tmp[n - i] * x[n - i]; // Backsubstitution.
		return x;
	}

	static void solve(vector<vector<int>>& skill2people, vector<long long>& people2skill, vector<int>& group, int state, int step, int& targetState, int& skillCnt, int& minGroup, vector<int>& ans) {
		if (state == targetState) {
			ans = group;
			minGroup = step;
		}
		else {
			if (step + 1 >= minGroup) return;
			for (int i = 0; i < skillCnt; i++) {
				if (!(state & (1 << i))) {
					for (int p : skill2people[i]) {
						group.push_back(p);
						solve(skill2people, people2skill, group, (state | people2skill[p]), step + 1, targetState, skillCnt, minGroup, ans);
						group.pop_back();
					}
					break;
				}
			}
		}
	}
};
