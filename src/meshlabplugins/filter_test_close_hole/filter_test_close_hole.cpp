/****************************************************************************
* MeshLab                                                           o o     *
* A versatile mesh processing toolbox                             o     o   *
*                                                                _   O  _   *
* Copyright(C) 2005                                                \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/

#include "filter_test_close_hole.h"
#include <QtGui>
#include <tuple>
#include <cmath>
#include <numeric>
#include <map>

#include <vcg/complex/algorithms/hole.h>
#include <vcg/space/color4.h>

//#define _N_DEBUG1
#ifdef _N_DEBUG1
#define N_LOG1(...) qDebug(__VA_ARGS__)
#else
#define N_LOG1(...)
#endif

#define _N_DEBUG_FIND_VERT
#ifdef _N_DEBUG_FIND_VERT
#define N_LOG_FIND_VERT(...) qDebug(__VA_ARGS__)
#else
#define N_LOG_FIND_VERT(...)
#endif

#define _N_DEBUG_FILL_VERT
#ifdef _N_DEBUG_FILL_VERT
#define N_LOG_FILL_VERT(...) qDebug(__VA_ARGS__)
#else
#define N_LOG_FILL_VERT(...)
#endif

using namespace std;
using namespace vcg;

// our extra functions

// functions declaration
const char* pointToString(Point3m p);
int solveQuadraticEquation(float a, float b, float c,float &x1, float &x2);
Point3m findHoleCenterPoint(CMeshO& cm, std::vector<int> hole);

// functions implementation
float distance2Points(Point3m p1, Point3m p2) 
{
	return sqrt(pow(p2.X() - p1.X(), 2) + pow(p2.Y() - p1.Y(), 2) + pow(p2.Z() - p1.Z(), 2));
}

float calcAvgDistance(std::vector<float> v_distance)
{
    return v_distance.size() == 0 ? 0 : std::accumulate(v_distance.begin(), v_distance.end(), decltype(v_distance)::value_type(0)) / v_distance.size();
}

Point3m calcAvgPoint(std::vector<int> vIndex, CMeshO& cm)
{
	assert(vIndex.size() != 0);

	Point3m rs(0, 0, 0);
	for (int index: vIndex)
	{
        rs += cm.vert[index].P();
	}
	
	return rs / vIndex.size();
}

/**
 *  Find next point to fill in the hole by create new Isosceles triangle
 *        * R
 *       *|*         RM is perpendicular with AB
 *      * | *        M is the midpoint of AB
 *     *  |  *
 *  A * * * * * B
 *        M
*/
Point3m findFilledVertByIsosceles(Point3m p1, Point3m p2, float filling_point_distance, Point3m hole_center, float edge_distance)
{
	N_LOG_FIND_VERT("\n\nStart log find vertext to fill");
	N_LOG_FIND_VERT("First point %s", pointToString(p1));
	N_LOG_FIND_VERT("Second point %s", pointToString(p2));
	N_LOG_FIND_VERT("Hole center point %s", pointToString(hole_center));
	N_LOG_FIND_VERT("border edge distance %f, filling_point_distance %f", edge_distance, filling_point_distance);
	// Find the height of Isosceles
	float h = sqrt(pow(filling_point_distance, 2) - pow(edge_distance / 2, 2));
	// Find the midpoint of border edge
	Point3m pM = (p1 + p2) / 2;
	N_LOG_FIND_VERT("Midpoint of border %s, isosceles triangle height %f", pointToString(pM), h); 
	// Find coord vector AB 
	Point3m vAB = p2 - p1;
	N_LOG_FIND_VERT("Vector border edge %s", pointToString(vAB)); 
	// Find the new z of filling point
	float z = (p1.Z() + p2.Z()) / 2;
	// We have a system of equations with M(xm, ym, zm) and R(xr, yr, zr)
	// vMR . vAB(a, b, c) = 0 => (xr - xm)*a + (yr - ym)*b + (zr- zm)*c = 0
	// MR = h => (xr-xm)^2 + (yr-ym)^2 + (zr-zm)^2 = h
	// dx = xr - xm, dy = yr - ym, dz = zr - zm => dr*a + dy*b + dz*c = 0 => dx = - (dy*b + dz*c) / a
	// => ((dy*b + dz*c) / a )^2 + dy^2 + dz^2 = h^2
	// => (dy*b + dz*c)^2 + a^2*dy^2 + a^2*dz^x - a^2*h^2 = 0
	// => (b + a^2)*dy^x
	float dx1 = 0, dx2 = 0;
	float dy1 = 0, dy2 = 0;
	float dz = z - pM.Z();

	// with a.x^2 + b.x + c = 0
	float a = pow(vAB.X(), 2) + pow(vAB.Y(), 2);
	float b = 2 * pow(vAB.X(), 2) * vAB.Y() * vAB.Z() * dz;
	float c = (pow(vAB.X(), 2) * dz + pow(vAB.Z(), 2)) * dz - pow(vAB.X(), 2) * pow(h, 2);
	N_LOG_FIND_VERT("Params to solve equation %f %f %f", a, b, c); 

	solveQuadraticEquation(a, b, c, dy1, dy2);

	dx1 = - (dy1 * vAB.Y() + dz * vAB.Z()) / vAB.X();
	dx2 = - (dy2 * vAB.Y() + dz * vAB.Z()) / vAB.X();
	N_LOG_FIND_VERT("Coord of vector MR dy1 %f, dy2 %f, dx1 %f, dx2 %f", dy1, dy2, dx1, dx2);

	// N_LOG_FIND_VERT("Test cross product vAB and vMR with dy1: %f", vAB.X() * dx1 + vAB.Y() * dy1 + vAB.Z() * dz);
	// N_LOG_FIND_VERT("Test cross product vAB and vMR with dy2: %f", vAB.X() * dx2 + vAB.Y() * dy2 + vAB.Z() * dz);
	// N_LOG_FIND_VERT("Test quadratic equation dy1: %f", a * pow(dy1, 2) + b * dy1 + dz);
	// N_LOG_FIND_VERT("Test quadratic equation dy2: %f", a * pow(dy2, 2) + b * dy2 + dz);
	// N_LOG_FIND_VERT("Test length of MR and height with dy1: %f", sqrt(pow(dx1, 2) + pow(dy1, 2) + pow(dz, 2)));
	// N_LOG_FIND_VERT("Test length of MR and height with dy2: %f", sqrt(pow(dx2, 2) + pow(dy2, 2) + pow(dz, 2)));

	// we have 2 points R1(x1, y1, z), R2(x2, y2, z)
	float x1 = dx1 + pM.X();
	float y1 = dy1 + pM.Y();
	float x2 = dx2 + pM.X();
	float y2 = dy2 + pM.Y();
	Point3m pR1(x1, y1, z);
	Point3m pR2(x2, y2, z);
	N_LOG_FIND_VERT("Midpoint of border %s, isosceles triangle height %f", pointToString(pM), h); 
	N_LOG_FIND_VERT("pR1  %s", pointToString(pR1)); 
	N_LOG_FIND_VERT("pR2  %s", pointToString(pR2)); 

	// N_LOG_FIND_VERT("Test length of MR and height with pR1: %f", distance2Points(pM, pR1));
	// N_LOG_FIND_VERT("Test length of MR and height with pR2: %f", distance2Points(pM, pR2));

	// get the length of R1C and R2C with C is the hole_center
	float dc1 = distance2Points(pR1, hole_center);
	float dc2 = distance2Points(pR2, hole_center);
	N_LOG_FIND_VERT("Result of fill point %s bc dc1 = %f, dc2 = %f", pointToString(dc1 > dc2 ? pR2 : pR1), dc1, dc2); 

	N_LOG_FIND_VERT("End log find vertext to fill \n\n");
	// get the point closer to center
	return dc1 > dc2 ? pR2 : pR1;
}

int solveQuadraticEquation(float a, float b, float c, float &x1, float &x2){
    float delta = b*b - 4*a*c;
    if(delta<0){
        x1=x2=0.0;
        return 0;
    }
    else if(delta==0){
        x1 = x2 = -b/(2*a);
        return 1;
    }
    else{
        delta = sqrt(delta);
        x1 = (-b + delta) / (2*a);
        x2 = (-b - delta) / (2*a);
        return 2;
    }
}

const char* pointToString(Point3m p)
{
	string s = "(x, y, z) - (" + to_string(p.X()) + ", " + to_string(p.Y()) + ", " + to_string(p.Z()) + ")";
	return s.c_str();
}

void fillHoleByIsoscelesTriangle(CMeshO& cm, std::vector<int> hole, std::vector<float> vDistance, float nextPointDistance)
{
	if (hole.size() == 0) 
	{
		return;
	}

    Point3m centerPoint = findHoleCenterPoint(cm, hole);
	// map<int, vector<int>> fillPointMap;
	int prevFillIndex = -1;
	for (int i = 0; i < hole.size(); ++i)
	{
        int nextIndex = i == hole.size() - 1 ? 0 : i + 1;
		N_LOG_FILL_VERT("index of vert %i %i", hole[i], hole[nextIndex]);
		int firstIndex = hole[i];
		int secondIndex = hole[nextIndex];
		auto firstVert = cm.vert[hole[i]];
		auto secondVert = cm.vert[hole[nextIndex]];
		N_LOG_FILL_VERT("firstVert %s", pointToString(firstVert.P()));
		N_LOG_FILL_VERT("secondVert %s", pointToString(secondVert.P()));

        Point3m fillPoint = findFilledVertByIsosceles(firstVert.P(), secondVert.P(), nextPointDistance, centerPoint, vDistance[i]);
		N_LOG_FILL_VERT("filling point %s", pointToString(fillPoint));

		// check last filled point in map to group near point
		if (prevFillIndex != -1)
		{
			Point3m prevFillPoint = cm.vert[prevFillIndex].P();
			if (distance2Points(prevFillPoint, fillPoint) < nextPointDistance)
			{
				// make average current and last filled point, should be avg of all component points
                cm.vert[prevFillIndex].P() = (prevFillPoint + fillPoint) / 2;
				vcg::tri::Allocator<CMeshO>::AddFace(cm, firstIndex, secondIndex, prevFillIndex);
				continue;
			}
		}
		

		// add point to mesh
        CMeshO::VertexIterator vi = vcg::tri::Allocator<CMeshO>::AddVertices(cm, 1);
        vi->P() = fillPoint;
		prevFillIndex = vi->Index();
        N_LOG_FILL_VERT("mesh vert size %i, first index %i, second index %i, new index %i", cm.vert.size(), firstIndex, secondIndex, vi->Index());
        vcg::tri::Allocator<CMeshO>::AddFace(cm, firstIndex, secondIndex, vi->Index());
	}

	return;
}

Point3m findHoleCenterPoint(CMeshO& cm, std::vector<int> hole)
{
	assert(hole.size() > 0);

	Point3m centerPoint(0, 0, 0);
	for (int idx: hole)
	{
		centerPoint += cm.vert[idx].P();
	}

	return centerPoint / hole.size();
}

float calcAvgHoleEdge(CMeshO& cm, std::vector<int> hole)
{
	if (hole.size() == 0) 
	{
		return 0;
	}

    float d = 0;
    Point3m prevPoint = cm.vert[hole.back()].P();
    for (int index: hole)
	{
		Point3m currPoint = cm.vert[index].P();
		d += distance2Points(prevPoint, currPoint);

		prevPoint = currPoint;
	}

	return d / hole.size();
}

void fillHoleByCenter(CMeshO& cm, std::vector<int> hole, float extra, float ratio)
{
    if (hole.size() <= 3) {
        return;
    }

	Point3m centerP(0, 0, 0);
	float maxZ = cm.vert[hole[0]].P().Z();
	for (int idx: hole)
	{
		centerP += cm.vert[idx].P();
		qDebug("point x, y, z, index %i (%f, %f, %f)", cm.vert[idx].Index(), cm.vert[idx].P().X(), cm.vert[idx].P().Y(), cm.vert[idx].P().Z());
		if (cm.vert[idx].P().Z() > maxZ) {
			maxZ = cm.vert[idx].P().Z();
		}
	}
	centerP /= hole.size();
	qDebug("hole center point x, y, z (%f, %f, %f) \n\n", centerP.X(), centerP.Y(), centerP.Z());

	CMeshO::VertexIterator vi = vcg::tri::Allocator<CMeshO>::AddVertices(cm, 1);
	// centerP.Z() = centerP.Z() * ratio;
	centerP.Z() = maxZ;
	vi->P() = centerP;
	vi->C() = vcg::Color4b(0, 255, 255, 255);
	// CMeshO::VertexIterator vi = vcg::tri::Allocator<CMeshO>::AddVertex(cm, centerP);

	int prevI = -1;
	int firstI = -1;
	for (int idx: hole)
	{
		if (prevI == -1)
		{
			prevI = idx;
			firstI = idx;
			continue;
		}

		qDebug("new mesh point 1 x, y, z, index %i (%f, %f, %f)", cm.vert[prevI].P().X(), cm.vert[prevI].P().Y(), cm.vert[prevI].P().Z(), cm.vert[prevI].Index());
		qDebug("new mesh point 2 x, y, z, index %i (%f, %f, %f)", cm.vert[idx].P().X(), cm.vert[idx].P().Y(), cm.vert[idx].P().Z(), cm.vert[idx].Index());
		qDebug("new mesh point 3 x, y, z, index %i (%f, %f, %f)", vi->P().X(), vi->P().Y(), vi->P().Z(), vi->Index());

		// CMeshO::FaceIterator fi = vcg::tri::Allocator<CMeshO>::AddFaces(cm, 1);
		// fi->V(0)=prevP;
		// fi->V(1)=v;
		// fi->V(2)=&*vi;
		vcg::tri::Allocator<CMeshO>::AddFace(cm, prevI, idx, vi->Index());

		prevI = idx;
	}
	if (firstI != -1) {
		vcg::tri::Allocator<CMeshO>::AddFace(cm, firstI, prevI, vi->Index());
	}

	return;
}

bool checkHoleSize(CMeshO& cm, std::vector<int> hole, float threshold, Point3m centerPoint)
{
	float totalDistance = 0;
	for (int i = 0; i < hole.size(); i++)
	{
		int index = hole[i];
		float dCenter = distance2Points(cm.vert[index].P(), centerPoint);
		totalDistance += dCenter;
	}

	return totalDistance / hole.size() < threshold;
}

void fillHoleRingByRing(CMeshO& cm, std::vector<int> hole, float threshold, bool stepByStep, std::vector<float> vRatio, float avgZRatio)
{
	if (hole.size() == 0) 
	{
		return;
	}

    Point3m centerPoint = findHoleCenterPoint(cm, hole);
    centerPoint.Z() = centerPoint.Z() * avgZRatio;
    // centerPoint.Z() = centerPoint.Z() * pow(1.073515, (distance2Points(cm.vert[hole[0]].P(), centerPoint) / threshold));
    // centerPoint.Z() = centerPoint.Z() * testRatio;
	// float startAvgEdge = calcAvgHoleEdge(cm, hole);

	while (true)
	{
		float avgEdge = calcAvgHoleEdge(cm, hole);
		// if (avgEdge < threshold)
		if (checkHoleSize(cm, hole, threshold, centerPoint) || avgEdge < threshold)
		{
			fillHoleByCenter(cm, hole, avgEdge, 1);
			break;
		}
		std::vector<int> reducedHole;

		int prevIndex = hole[hole.size() - 1];
		int prevFilledIndex = -1;
		int firstFilledIndex = -1;
		for (int i = 0; i < hole.size(); i++)
		{
			int index = hole[i];
			float dCenter = distance2Points(cm.vert[index].P(), centerPoint);
			float k = avgEdge / dCenter;
			float ratio = vRatio[i];
			qDebug("distance to center %f and ratio %f", dCenter, k);
			assert(dCenter == dCenter);
            Point3m fillPoint = ((centerPoint - cm.vert[index].P()) * k) + cm.vert[index].P();

			// check last filled point in map to group near point
			// if (prevFilledIndex != -1)
			// {
			// 	Point3m prevFillPoint = cm.vert[prevFilledIndex].P();
			// 	if (distance2Points(prevFillPoint, fillPoint) < avgEdge/2)
			// 	{
			// 		// make average current and last filled point, should be avg of all component points
			// 		// cm.vert[prevFilledIndex].P() = (prevFillPoint + fillPoint) / 2;
			// 		vcg::tri::Allocator<CMeshO>::AddFace(cm, prevIndex, index, prevFilledIndex);

			// 		prevIndex = index;
			// 		reducedHole.push_back(index);

			// 		continue;
			// 	}
			// }

			// add point to mesh
			CMeshO::VertexIterator vi = vcg::tri::Allocator<CMeshO>::AddVertices(cm, 1);
			qDebug("LogRatio During fill Border Vertex index %i next index %i coord (x, y, z): (%f, %f, %f) ratio %f \n", index, vi->Index(),
				cm.vert[index].P().X(), cm.vert[index].P().Y(), cm.vert[index].P().Z(), ratio);
			// fillPoint.Z() = cm.vert[index].P().Z() * 1.073515;
			vi->P() = fillPoint;
			vi->C() = vcg::Color4b(0, 255, 255, 255);
			vcg::tri::Allocator<CMeshO>::AddFace(cm, prevIndex, index, vi->Index());
			if (i > 0)
			{
            	vcg::tri::Allocator<CMeshO>::AddFace(cm, prevFilledIndex, prevIndex, vi->Index());
			} else {
				firstFilledIndex = vi->Index();
			}

			reducedHole.push_back(vi->Index());

			prevIndex = index;
			prevFilledIndex = vi->Index();
		}
		if (firstFilledIndex != prevFilledIndex) {
        	vcg::tri::Allocator<CMeshO>::AddFace(cm, prevFilledIndex, prevIndex, firstFilledIndex);
		}

		hole = reducedHole;

		if (stepByStep) {
			break;
		}
	}
	
	return;
}

// end extra functions

/**
 * @brief Constructor usually performs only two simple tasks of filling the two lists
 *  - typeList: with all the possible id of the filtering actions
 *  - actionList with the corresponding actions. If you want to add icons to your filtering actions you can do here by construction the QActions accordingly
 */
FilterFillHolePlugin::FilterFillHolePlugin()
{
	typeList = {
        FP_TEST_DP_CLOSE_HOLE,
		FP_TEST_CLOSE_HOLE
    };

	QCoreApplication *app = QCoreApplication::instance();
    for (ActionIDType tt : types())
    {
        QAction *act = new QAction(filterName(tt), this);
        actionList.push_back(act);

        if (app != nullptr)
        {
            if (tt == FP_TEST_CLOSE_HOLE)
            {
                //				act->setShortcut(QKeySequence ("Ctrl+Del"));
                act->setIcon(QIcon(":/images/fill_hole.svg"));
                //				act->setPriority(QAction::HighPriority);
            }
        }
    }
}

QString FilterFillHolePlugin::pluginName() const
{
	return "FilterCloseHole";
}

FilterFillHolePlugin::FilterClass FilterFillHolePlugin::getClass(const QAction * a) const
{
	switch (ID(a))
	{
	case FP_TEST_DP_CLOSE_HOLE: return FilterPlugin::Remeshing;
	case FP_TEST_CLOSE_HOLE: return FilterPlugin::Remeshing;

	default: assert(0); return FilterPlugin::Generic;
	}

	return FilterPlugin::Generic;
}

/**
 * @brief ST() must return the very short string describing each filtering action
 * (this string is used also to define the menu entry)
 * @param filterId: the id of the filter
 * @return the name of the filter
 */
QString FilterFillHolePlugin::filterName(ActionIDType filterId) const
{
	switch (filterId)
	{
	case FP_TEST_DP_CLOSE_HOLE:
		return "surface fill hole";
	case FP_TEST_CLOSE_HOLE:
		return "find mesh hole";
	default:
		return "";
	}
}

/**
 * @brief // Info() must return the longer string describing each filtering action
 * (this string is used in the About plugin dialog)
 * @param filterId: the id of the filter
 * @return an info string of the filter
 */
QString FilterFillHolePlugin::filterInfo(ActionIDType filterId) const
{
	switch (filterId)
	{
	case FP_TEST_DP_CLOSE_HOLE:
		return "fill hole algo.";
	case FP_TEST_CLOSE_HOLE:
		return "find hole algo.";
	default:
		return "Unknown Filter";
	}
}

/**
 * @brief FilterFillHolePlugin::getPreConditions
 * @return
 */
int FilterFillHolePlugin::getPreConditions(const QAction *action) const
{
    switch (ID(action))
	{
	case FP_TEST_DP_CLOSE_HOLE:
		return MeshModel::MM_FACENUMBER;
	case FP_TEST_CLOSE_HOLE:
		return MeshModel::MM_FACENUMBER;
	}

	return MeshModel::MM_NONE;
}

/**
 * @brief FilterFillHolePlugin::postCondition
 * @return
 */
int FilterFillHolePlugin::postCondition(const QAction *action) const
{
	switch (ID(action))
	{
	case FP_TEST_DP_CLOSE_HOLE:
		return MeshModel::MM_GEOMETRY_AND_TOPOLOGY_CHANGE;
	case FP_TEST_CLOSE_HOLE:
		return MeshModel::MM_GEOMETRY_AND_TOPOLOGY_CHANGE;
	}

	return MeshModel::MM_ALL;
}

/**
 * @brief This function define the needed parameters for each filter. Return true if the filter has some parameters
 * it is called every time, so you can set the default value of parameters according to the mesh
 * For each parameter you need to define,
 * - the name of the parameter,
 * - the default value
 * - the string shown in the dialog
 * - a possibly long string describing the meaning of that parameter (shown as a popup help in the dialog)
 * @param action
 * @param m
 * @param parlst
 */
void FilterFillHolePlugin::initParameterList(const QAction *action, MeshModel &m, RichParameterList &parlst)
{
	float maxVal;
	QStringList curvCalcMethods;
	QStringList curvColorMethods;
	QStringList loopWeightLst;

	switch (ID(action))
	{
	case FP_TEST_DP_CLOSE_HOLE:
		parlst.addParam(RichInt ("MaxHoleSize",(int)30,"Max size to be closed ","The size is expressed as number of edges composing the hole boundary"));
		parlst.addParam(RichBool("Selected",m.cm.sfn>0,"Close holes with selected faces","Only the holes with at least one of the boundary faces selected are closed"));
		parlst.addParam(RichBool("NewFaceSelected",true,"Select the newly created faces","After closing a hole the faces that have been created are left selected. Any previous selection is lost. Useful for example for smoothing the newly created holes."));
		parlst.addParam(RichBool("SelfIntersection",true,"Prevent creation of selfIntersecting faces","When closing an holes it tries to prevent the creation of faces that intersect faces adjacent to the boundary of the hole. It is an heuristic, non intersetcting hole filling can be NP-complete."));
		break;
	case FP_TEST_CLOSE_HOLE:
	{
		QStringList algoType;
		algoType.push_back("Ring by ring");
		algoType.push_back("Center point");
		algoType.push_back("Isosceles");
		algoType.push_back("Get hole information");
		parlst.addParam(RichEnum("algo", 0, algoType, tr("Algorithm type"), tr("Choose the algorithm to close hole")));
        parlst.addParam(RichFloat("threshold", 0, "Threshold", "Set a threshold > 0 for filling hole step by step"));
        parlst.addParam(RichFloat("ratio", 0, "Ratio", "User test ratio for Z"));
        parlst.addParam(RichInt("maxholesize", int(0), "Max hole size", "Size of a hole is the number of boundary face of that hole"));
		parlst.addParam(RichBool("enableBorderColor", false, "Change color of border face and border vertex"));
	}
		break;
	default:
		break;
	}
}

/**
 * @brief The Real Core Function doing the actual mesh processing.
 * @param action
 * @param md: an object containing all the meshes and rasters of MeshLab
 * @param par: the set of parameters of each filter
 * @param cb: callback object to tell MeshLab the percentage of execution of the filter
 * @return true if the filter has been applied correctly, false otherwise
 */
std::map<std::string, QVariant> FilterFillHolePlugin::applyFilter(
    const QAction *action,
    const RichParameterList &par,
    MeshDocument &md,
    unsigned int& /*postConditionMask*/,
    vcg::CallBackPos *cb)
{
	MeshModel &m = *(md.mm());
	switch (ID(action))
	{
	case FP_TEST_DP_CLOSE_HOLE:
	{
		m.updateDataMask(MeshModel::MM_FACEFACETOPO);
		if (tri::Clean<CMeshO>::CountNonManifoldEdgeFF(m.cm) > 0)
		{
			throw MLException("Mesh has some not 2-manifold edges, filter requires edge manifoldness");
		}

		size_t OriginalSize = m.cm.face.size();
		int MaxHoleSize = par.getInt("MaxHoleSize");
		bool SelectedFlag = par.getBool("Selected");
		bool SelfIntersectionFlag = par.getBool("SelfIntersection");
		bool NewFaceSelectedFlag = par.getBool("NewFaceSelected");
		int holeCnt;
		if (SelfIntersectionFlag)
			holeCnt = tri::Hole<CMeshO>::EarCuttingIntersectionFill<tri::SelfIntersectionEar<CMeshO>>(m.cm, MaxHoleSize, SelectedFlag, cb);
		else
			holeCnt = tri::Hole<CMeshO>::EarCuttingFill<vcg::tri::MinimumWeightEar<CMeshO>>(m.cm, MaxHoleSize, SelectedFlag, cb);
		log("Closed %i holes and added %i new faces",holeCnt,m.cm.fn-OriginalSize);
		assert(tri::Clean<CMeshO>::IsFFAdjacencyConsistent(m.cm));
		m.UpdateBoxAndNormals();

		// hole filling filter does not correctly update the border flags (but the topology is still ok!)
		if (NewFaceSelectedFlag)
		{
			tri::UpdateSelection<CMeshO>::FaceClear(m.cm);
			for (size_t i = OriginalSize; i < m.cm.face.size(); ++i)
				if (!m.cm.face[i].IsD())
					m.cm.face[i].SetS();
		}
		break;
	}
	case FP_TEST_CLOSE_HOLE:
	{
		m.updateDataMask(MeshModel::MM_FACEFACETOPO);
		m.updateDataMask(MeshModel::MM_VERTCOLOR);
		m.updateDataMask(MeshModel::MM_FACECOLOR);
		CMeshO & cm = m.cm;
		if (tri::Clean<CMeshO>::CountNonManifoldEdgeFF(m.cm) > 0)
		{
			throw MLException("Mesh has some not 2-manifold edges, filter requires edge manifoldness");
		}

		std::vector<tri::Hole<CMeshO>::Info> vinfo;
		bool Selected = false; // TODO: put in param
		std::vector<std::vector<CVertexO*>> vholeV;
		int expandedVBit = CVertexO::NewBitFlag();
		// int borderVBit = CVertexO::NewBitFlag();
		int borderSharedBit = 15;
        std::vector<std::vector<int>> vholeI;
        std::vector<std::tuple<std::vector<int>, std::vector<float>, std::vector<float>, float>> vholeInfo;

		std::vector<std::vector<int>> vHoleFaceIndex;

        // vcg::tri::Hole<CMeshO>::GetInfo(cm, Selected, vinfo);

		tri::UpdateFlags<CMeshO>::FaceClearV(cm);

        for(tri::Hole<CMeshO>::FaceIterator fi = cm.face.begin(); fi!=cm.face.end(); ++fi)
        {
			if(!(*fi).IsD())
			{
				if(Selected && !(*fi).IsS())
				{
					//if I have to consider only the selected triangles e
					//what I'm considering isn't marking it and moving on
					(*fi).SetV();
				}
				else
				{
						for(int j =0; j<3 ; ++j)
						{
							if( face::IsBorder(*fi,j) && !(*fi).IsV() )
							{//Found a board face not yet visited.
                                (*fi).SetV();
								tri::Hole<CMeshO>::PosType sp(&*fi, j, (*fi).V(j));
								tri::Hole<CMeshO>::PosType fp=sp;
								int holesize=0;

                                // (*fi).C() = vcg::Color4b(255, 0, 255, 255);
								// CVertexO* v1 = (*fi).V(j);
								// (*v1).C() = vcg::Color4b(255, 0, 255, 255);
								// CVertexO* v2 = (*fi).V((j + 1) % 3);
								// (*v2).C() = vcg::Color4b(255, 0, 255, 255);

								qDebug("Hole %i detected \n", vinfo.size() + 1);
								std::vector<CVertexO*> vBorderVertex;
								std::vector<int> vBorderIndex;
								std::vector<float> vDistanceVert;
								std::vector<float> vRatio;
								float totalZ = 0;
								int countZ = 0;
								float totalExpandZ = 0;
								int countExpandZ = 0;
								Point3m prevPoint = sp.v->P();
								bool hasPrevP = true;
								std::vector<int> vFaceIndex;					
								float ratio = 0;
								float totalRatio = 0;
								int countRatio = 0;

								tri::Hole<CMeshO>::Box3Type hbox;
								hbox.Add(sp.v->cP());
								//printf("Looping %i : (face %i edge %i) \n", VHI.size(),sp.f-&*m.face.begin(),sp.z);
								sp.f->SetV();
								do
								{
									sp.f->SetV();
									hbox.Add(sp.v->cP());
									++holesize;
									sp.NextB();
									sp.f->SetV();

									// Set corlor of border face, vertex
                                    // sp.f->C().SetHSVColor(0, 1.0f, 1.0f);
									// Set corlor of extended boundary
									sp.v->SetUserBit(borderSharedBit);
                                    CVertexO* expandedFP = sp.f->V2(sp.z);

									totalZ += sp.v->P().Z();
									++countZ;
									if (!expandedFP->IsUserBit(borderSharedBit)) {
										// expandedFP->C() = vcg::Color4b(255, 255, 0, 255);
										totalExpandZ += expandedFP->P().Z();
										++countExpandZ;
									}

									qDebug("LogRatio Border Vertex index %i coord (x, y, z): (%f, %f, %f) index %d \n", 
										sp.v->Index(), sp.v->P().X(), sp.v->P().Y(), sp.v->P().Z(), sp.v->Index());
									qDebug("LogRatio Extended Vertex index %i coord (x, y, z): (%f, %f, %f) index %d \n", 
										expandedFP->Index(), expandedFP->P().X(), expandedFP->P().Y(), expandedFP->P().Z(), expandedFP->Index());
									ratio =  sp.v->P().Z() / expandedFP->P().Z();
									qDebug("LogRatio Ratio %f \n", ratio);
									++countRatio;
									totalRatio += ratio;

									if (ratio < 1) {
										// sp.v->C() = vcg::Color4b(255, 0, 255, 255);
										// expandedFP->C() = vcg::Color4b(255, 255, 0, 255);
									}

									int vIndex = sp.v->Index();

                                	// qDebug("Border Vertex index %i coord (x, y, z): (%f, %f, %f) \n", vIndex, sp.v->P().X(), sp.v->P().Y(), sp.v->P().Z());
									vBorderVertex.push_back(sp.v);
									vBorderIndex.push_back(sp.v->Index());
									vRatio.push_back(ratio);
									vFaceIndex.push_back(sp.f->Index());

									if (!hasPrevP) 
									{
										hasPrevP = true;
										prevPoint = sp.v->P();
                                    } else {
                                        float d = distance2Points(sp.v->P(), prevPoint);
										assert(d);
                                        vDistanceVert.push_back(d);
										// set prevP after calc distance
										prevPoint = sp.v->P();
									}

									assert(sp.IsBorder());
								}while(sp != fp);

								// vDistanceVert.push_back(distance2Points(sp.v->P(), prevPoint));
								log("Ratio %f", totalRatio / countRatio);

								qDebug("End hole point log");
								vholeV.push_back(vBorderVertex);
								vholeI.push_back(vBorderIndex);
								float averageZRatio = (totalZ/countZ) / (totalExpandZ/countExpandZ);
								if (averageZRatio != averageZRatio) {
									averageZRatio = 1;
								}
								// assert(averageZRatio == averageZRatio);
                                vholeInfo.push_back(std::make_tuple(vBorderIndex, vDistanceVert, vRatio, averageZRatio));
								vHoleFaceIndex.push_back(vFaceIndex);

								//I recovered the information on the whole hole
                                vinfo.push_back( tri::Hole<CMeshO>::Info(sp,holesize,hbox) );
							}
						}//for on the edges of the triangle
				}//S & !S
			}//!IsD()
		}//for principale!!!

		bool enableBorderColor = par.getBool("enableBorderColor");
		float threshold = par.getFloat("threshold");
		float userRatio = par.getFloat("ratio");
        float maxHoleSize = par.getInt("maxholesize");
		if (enableBorderColor)
		{
            for (std::tuple<std::vector<int>, std::vector<float>, std::vector<float>, float> hole: vholeInfo)
			{
				std::vector<int> vVertIndex;
				std::vector<float> vDistance;
				std::vector<float> vRatio;
                float avgZRatio;
                tie(vVertIndex, vDistance, vRatio, avgZRatio) = hole;

				if (maxHoleSize > 0 && vVertIndex.size() > maxHoleSize) {
					continue;
				}

				for (int i = 0; i < vVertIndex.size(); i++)
				{
                    int vi = vVertIndex[i];
					cm.vert[vi].C() = vcg::Color4b(255, 0, 255, 255);
				}
			}

			for (std::vector<int> faceIdx: vHoleFaceIndex)
			{
				if (faceIdx.size() > 30) {
					continue;
				}

				for (int fi: faceIdx)
				{
					cm.face[fi].C() = vcg::Color4b(255, 0, 0, 255);
					// cm.face[fi].C().SetHSVColor(0, 1.0f, 1.0f);
					// break;
				}
			}
		}

        
		bool stepByStep = true;
		switch(par.getEnum("algo"))
		{
			case 0:
			{
				// 1. calc average edge as threshold
				// 2. add new points which in the line to center
				// 3. when distance to center < threshold fill by center
                for (std::tuple<std::vector<int>, std::vector<float>, std::vector<float>, float> hole: vholeInfo)
				{
					std::vector<int> vVertIndex;
					std::vector<float> vDistance;
					std::vector<float> vRatio;
					float avgZRatio;
					tie(vVertIndex, vDistance, vRatio, avgZRatio) = hole;
                    if (maxHoleSize > 0 && vVertIndex.size() > maxHoleSize) {
						continue;
					}
					float avgDistance = calcAvgDistance(vDistance);
					if (threshold <= 0) {
                        threshold = avgDistance/2;
						stepByStep = false;
					}
					if (userRatio > 0) {
						avgZRatio = userRatio;
					}
					qDebug("start one hole filling with threshold %f", threshold);

					fillHoleRingByRing(cm, vVertIndex, threshold, stepByStep, vRatio, avgZRatio);

					qDebug("End one hole filling");
				}
			}
				break;
			case 1: 
			{
				for (std::vector<int> hole: vholeI) 
				{
					if (maxHoleSize > 0 && hole.size() > maxHoleSize) {
						continue;
					}

                    fillHoleByCenter(cm, hole, 0, 1);
				}
			} // end case 1
				break;
			case 2: 
			{
				// 1. calculate average distance
				// 2. add points by new average distance
				// 3. modify and improve
                for (std::tuple<std::vector<int>, std::vector<float>, std::vector<float>, float> hole: vholeInfo)
				{
					std::vector<int> vVertIndex;
					std::vector<float> vDistance;
                    std::vector<float> vRatio;
					float avgZRatio;
                    tie(vVertIndex, vDistance, vRatio, avgZRatio) = hole;
                    if (maxHoleSize > 0 && vVertIndex.size() > maxHoleSize) {
						continue;
					}
					float avgDistance = calcAvgDistance(vDistance);
					qDebug("start one hole filling with distance %f", avgDistance);

					fillHoleByIsoscelesTriangle(cm, vVertIndex, vDistance, avgDistance);

					qDebug("End one hole filling");
				}

			}
				break;
			case 3:
			{
				int count = 0;
				for (std::tuple<std::vector<int>, std::vector<float>, std::vector<float>, float> hole: vholeInfo)
				{
					std::vector<int> vVertIndex;
					std::vector<float> vDistance;
					std::vector<float> vRatio;
					float avgZRatio;
					tie(vVertIndex, vDistance, vRatio, avgZRatio) = hole;
                    if (maxHoleSize > 0 && vVertIndex.size() > maxHoleSize) {
						continue;
					}
					float avgDistance = calcAvgDistance(vDistance);

					++count;
					log("Hole number %i size %i average edge %f avg ratio %f", count, vVertIndex.size(), avgDistance, avgZRatio);
				}
			} break;
			default:
				assert(0);
		} // end switch

        log("Found %i holes",vinfo.size());
		break;
	}

	default:
        wrongActionCalled(action);
	}
	return std::map<std::string, QVariant>();
}

MESHLAB_PLUGIN_NAME_EXPORTER(FilterFillHolePlugin)
