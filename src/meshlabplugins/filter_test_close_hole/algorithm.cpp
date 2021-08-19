#include "filter_test_close_hole.h"
#include "algorithm.h"

#include <cmath>
#include <numeric>
#include <map>

#include <QtGui>

#define _N_DEBUG_FILL_VERT
#ifdef _N_DEBUG_FILL_VERT
#define N_LOG_FILL_VERT(...) qDebug() << __VA_ARGS__
#else
#define N_LOG_FILL_VERT(...)
#endif

// #define _N_DEBUG_FILL_CENTER
#ifdef _N_DEBUG_FILL_CENTER
#define N_LOG_FILL_CENTER(...) qDebug(__VA_ARGS__)
#else
#define N_LOG_FILL_CENTER(...)
#endif


// functions implementation
float distance2Points(Point3m p1, Point3m p2) 
{
	return sqrt(pow(p2.X() - p1.X(), 2) + pow(p2.Y() - p1.Y(), 2) + pow(p2.Z() - p1.Z(), 2));
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

const char* pointToString(Point3m p)
{
    std::string s = "(x, y, z) - (" + std::to_string(p.X()) + ", " + std::to_string(p.Y()) + ", " + std::to_string(p.Z()) + ")";
	return s.c_str();
}

QString PointToQString(Point3m p) {

	return QString("(x, y, z) - (%1, %2, %3)").arg(QString::number(p.X()), QString::number(p.Y()), QString::number(p.Z()));
}

float CalcAvgHoleEdge(CMeshO& cm, std::vector<int> hole)
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

Point3m CalcHoleCenter(CMeshO& cm, std::vector<int> hole) {
	Point3m centerP(0, 0, 0);
	for (int idx: hole)
	{
		centerP += cm.vert[idx].P();
	}
	centerP /= hole.size();

	return centerP;
}

int chooseEdgeNearAvg(float d1, float d2, float avgd) {
    assert(d1 > 0 && d2 > 0 && avgd > 0);
	float diff1 = abs(d1 - avgd);
	float diff2 = abs(d2 - avgd);
	
	return diff1 > diff2 ? 1 : 2;
}

float calcAvgDistanceToCenter(CMeshO& cm, std::vector<int> hole, Point3m centerPoint)
{
	float totalDistance = 0;
	for (int i = 0; i < hole.size(); i++)
	{
		int index = hole[i];
		float dCenter = distance2Points(cm.vert[index].P(), centerPoint);
		totalDistance += dCenter;
	}

	return totalDistance / hole.size();
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

Point3m calcFillingPoint(Point3m boundaryPoint, Point3m centerPoint, float avgEdge, float ratio) {
	float dCenter = distance2Points(boundaryPoint, centerPoint);
	float k = avgEdge / dCenter;

	N_LOG_FILL_VERT(QString("distance to center %1 and ratio %2").arg(QString::number(dCenter), QString::number(k)));
	assert(dCenter == dCenter);
	Point3m fillPoint = ((centerPoint - boundaryPoint) * k) + boundaryPoint;

	return fillPoint;
}

bool checkCurrPointEdgeOk(Point3m curFill, Point3m curBoundary, Point3m prevFill, Point3m prevBoundary) {
	float dCurr = distance2Points(curFill, prevBoundary);
	float dPrev = distance2Points(prevFill, curBoundary);

	return dCurr < dPrev;
}

bool checkNewPrevPointDistance(Point3m newFill, Point3m prevFill, float threshold) {
	float d = distance2Points(newFill, prevFill);

	return d > threshold;
}

float calcCenterZChange(CMeshO& cm, Point3m center, float avgEdge, std::vector<int> hole, std::vector<float> vZChange) {
	// 1. each boundary - ratio, check distance to center d, factor to center = d/avgEdge, rounded factor
	// 2. only get min rounded factor, total factor, count factor
	// 3. average of 2 is result
	int count = 0;
	float totalZChange = 0;
	int minFactor = 59999;
	for (int i = 0; i < hole.size(); i++) {
		int bindex = hole[i];
		Point3m bvertex = cm.vert[bindex].P();
		float zChange = vZChange[i];

		float d = distance2Points(bvertex, center);
		float factor = d / avgEdge;
		int roundFactor = round(factor);

		if (roundFactor < minFactor) {
			minFactor = roundFactor;
			totalZChange = zChange;
			count = 1;
		} else {
			totalZChange += zChange;
			++count;
		}
	}

    return (totalZChange / count) * minFactor;
}

Point3m CalcCenterChange(CMeshO& cm, Point3m center, float avgEdge, std::vector<int> hole, std::vector<Point3m> expHole) {

	int count = 0;
	Point3m totalChange(0, 0, 0);
	int minFactor = 59999;
	for (int i = 0; i < hole.size(); i++) {
		int bindex = hole[i];
		Point3m bvertex = cm.vert[bindex].P();
		Point3m cChange = expHole[i];

		float d = distance2Points(bvertex, center);
		float factor = d / avgEdge;
		int roundFactor = round(factor);

		if (roundFactor < minFactor) {
			minFactor = roundFactor;
			totalChange = cChange;
			count = 1;
		} else {
			totalChange += cChange;
			++count;
		}
	}

    return (totalChange / count) * minFactor;
}

int findStepToCenter(Point3m center, Point3m boundary, float avgEdge, float & dRatio) {
	// avgEdge *= 0.866;
	float d = distance2Points(center, boundary);
	dRatio = d / avgEdge;
	int step = round(dRatio);
	
	return step;
}

int findMaxStepToCenter(CMeshO& cm, Point3m center, float avgEdge, std::vector<int> hole) {
	int maxStep = 0;
	for (int index: hole) {
		Point3m boundary = cm.vert[index].P();
		float dRatio;
		int step = findStepToCenter(center, boundary, avgEdge, dRatio);
		if (step > maxStep) {
			maxStep = step;
		}
	}

	return maxStep;
}

std::vector<int> rearrangeHole(CMeshO& cm, Point3m center, float avgEdge, std::vector<int> hole) {
	int k = findMaxStepToCenter(cm, center, avgEdge, hole);
	int maxStepIndex = 0;
	std::vector<int> rearrange;

	for (int i = 0; i < hole.size(); ++i) {
        int index = hole[i];
		Point3m boundary = cm.vert[index].P();
		float dRatio;
		int step = findStepToCenter(center, boundary, avgEdge, dRatio);
		if (step == k) {
			maxStepIndex = i;
			break;
		}
	}
	for (int i = maxStepIndex + 1; i < hole.size(); ++i) {
		rearrange.push_back(hole[i]);
	}
	for (int i = 0; i <= maxStepIndex; ++i) {
		rearrange.push_back(hole[i]);
	}

	return rearrange;
}

std::vector<int> reduceHoleByConnectNearby(CMeshO& cm, std::vector<int> hole, float avgEdge, Point3m center) {
	std::vector<int> cHole = hole;
	bool complete = false;
	// std::vector<int> reducedHole;

	do {
		if (cHole.size() <= 3) {
			break;
		}

		complete = true;
		// reducedHole.clear();
		for (int i = 0; i < cHole.size() - 2; i++) {
			int firstIndex = cHole[i];
			int secondIndex = cHole[i+1];
			int thirdIndex = cHole[i+2];
			if (i == cHole.size() - 1) {
				secondIndex = cHole[0];
				thirdIndex = cHole[1];
			}
			Point3m firstPoint = cm.vert[firstIndex].P();
			Point3m secondPoint = cm.vert[secondIndex].P();
			Point3m thirdPoint = cm.vert[thirdIndex].P();
			Point3m mid13Point = (firstPoint + thirdPoint) / 2;
			float dSecond = distance2Points(secondPoint, center);
			float dMid13 = distance2Points(mid13Point, center);
			float d13 = distance2Points(firstPoint, thirdPoint);

			float d21 = distance2Points(secondPoint, firstPoint);
			float d23 = distance2Points(secondPoint, thirdPoint);

			if (d13 < avgEdge && dSecond > avgEdge && dMid13 < dSecond && (d13 < d21*1.2 || d13 < d23*1.2)) {
				// add face
				vcg::tri::Allocator<CMeshO>::AddFace(cm, firstIndex, secondIndex, thirdIndex);
				cm.face.back().C() = vcg::Color4b::Red;
				// NOTE: mutate the hole
                auto it = cHole.begin() + i + 1;
				cHole.erase(it);
				complete = false;
				// reducedHole.push_back(secondIndex);
				break;
			} 
			
			// reducedHole.push_back(firstIndex);
		}
	} while(!complete);

	// return reducedHole;
	return cHole;
}

float DistanceVector(Point3m p) {
	return sqrt(pow(p.X(), 2) + pow(p.Y(), 2) + pow(p.Z(), 2));
}

float CalcCenterZChangeUsingExpVertex(CMeshO& cm, Point3m center, float avgEdge, std::vector<int> hole, std::vector<int> expHole) {
	// 1. each boundary - ratio, check distance to center d, factor to center = d/avgEdge, rounded factor
	// 2. only get min rounded factor, total factor, count factor
	// 3. average of 2 is result
	int count = 0;
	float totalZChange = 0;
	int minFactor = 59999;
	for (int i = 0; i < hole.size(); i++) {
		int bindex = hole[i];
		Point3m bvertex = cm.vert[bindex].P();
		int expIdx = expHole[i];
		Point3m evertex = cm.vert[expIdx].P();

		float zChange = bvertex.Z() - evertex.Z();

		float d = distance2Points(bvertex, center);
		float factor = d / avgEdge;
		int roundFactor = round(factor);

		if (roundFactor < minFactor) {
			minFactor = roundFactor;
			totalZChange = zChange;
			count = 1;
		} else {
			totalZChange += zChange;
			++count;
		}
	}

    return (totalZChange / count) * minFactor;
}

void FillHoleByCenter(CMeshO& cm, std::vector<int> hole, float extra, float ratio)
{
    if (hole.size() <= 3) {
        return;
    }

	Point3m centerP(0, 0, 0);
	float maxZ = cm.vert[hole[0]].P().Z();
	for (int idx: hole)
	{
		centerP += cm.vert[idx].P();
		N_LOG_FILL_CENTER("point x, y, z, index %i (%f, %f, %f)", cm.vert[idx].Index(), cm.vert[idx].P().X(), cm.vert[idx].P().Y(), cm.vert[idx].P().Z());
		if (cm.vert[idx].P().Z() > maxZ) {
			maxZ = cm.vert[idx].P().Z();
		}
	}
	centerP /= hole.size();
	N_LOG_FILL_CENTER("hole center point x, y, z (%f, %f, %f) \n\n", centerP.X(), centerP.Y(), centerP.Z());

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

		N_LOG_FILL_CENTER("new mesh point 1 x, y, z, index %i (%f, %f, %f)", cm.vert[prevI].P().X(), cm.vert[prevI].P().Y(), cm.vert[prevI].P().Z(), cm.vert[prevI].Index());
		N_LOG_FILL_CENTER("new mesh point 2 x, y, z, index %i (%f, %f, %f)", cm.vert[idx].P().X(), cm.vert[idx].P().Y(), cm.vert[idx].P().Z(), cm.vert[idx].Index());
		N_LOG_FILL_CENTER("new mesh point 3 x, y, z, index %i (%f, %f, %f)", vi->P().X(), vi->P().Y(), vi->P().Z(), vi->Index());

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

void FillHoleByCenterRefined(CMeshO& cm, std::vector<int> hole, float extra, float ratio, Point3m centerP) {
	if (hole.size() < 3) {
        return;
    }

	if (hole.size() == 3) {
		int firstIdx = cm.vert[hole[0]].Index();
		int secondIdx = cm.vert[hole[1]].Index();
		int thirdIdx = cm.vert[hole[2]].Index();
		vcg::tri::Allocator<CMeshO>::AddFace(cm, firstIdx, secondIdx, thirdIdx);
		return;
	}

	// Point3m centerP(0, 0, 0);
	// float maxZ = cm.vert[hole[0]].P().Z();
	// for (int idx: hole)
	// {
	// 	centerP += cm.vert[idx].P();
	// 	qDebug("point x, y, z, index %i (%f, %f, %f)", cm.vert[idx].Index(), cm.vert[idx].P().X(), cm.vert[idx].P().Y(), cm.vert[idx].P().Z());
	// 	if (cm.vert[idx].P().Z() > maxZ) {
	// 		maxZ = cm.vert[idx].P().Z();
	// 	}
	// }
	// centerP /= hole.size();
	qDebug("hole center point x, y, z (%f, %f, %f) \n\n", centerP.X(), centerP.Y(), centerP.Z());

	CMeshO::VertexIterator vi = vcg::tri::Allocator<CMeshO>::AddVertices(cm, 1);
	int centerIdx = vi->Index();
	// centerP.Z() = centerP.Z() * ratio;
	// centerP.Z() = maxZ;
	vi->P() = centerP;
	vi->C() = vcg::Color4b(0, 255, 255, 255);
	// CMeshO::VertexIterator vi = vcg::tri::Allocator<CMeshO>::AddVertex(cm, centerP);

	int prevI = -1;
	int firstI = -1;
	for (int i = 0; i < hole.size() - 1; i++)
	{
		int idx = hole[i];

		int firstIdx = hole[i];
		int secondIdx = hole[(i + 1) % hole.size()];
		int thirdIdx = hole[(i + 2) % hole.size()];

		float d13 = distance2Points(cm.vert[firstIdx].P(), cm.vert[thirdIdx].P());
		float d1c = distance2Points(cm.vert[firstIdx].P(), centerP);
		float d2c = distance2Points(cm.vert[secondIdx].P(), centerP);

		float d21 = distance2Points(cm.vert[secondIdx].P(), cm.vert[firstIdx].P());
		float d23 = distance2Points(cm.vert[secondIdx].P(), cm.vert[thirdIdx].P());

		if (d13 < d2c && d2c > d1c && (d13 < d21*1.2 || d13 < d23*1.2)) {
			// add face 13c, 123; i++
			vcg::tri::Allocator<CMeshO>::AddFace(cm, firstIdx, thirdIdx, centerIdx);
			cm.face.back().C() = vcg::Color4b::Gray;
			vcg::tri::Allocator<CMeshO>::AddFace(cm, firstIdx, secondIdx, thirdIdx);
			cm.face.back().C() = vcg::Color4b::Gray;

			++i;
			
			if (i == hole.size() - 2) {
				// add face 3c & first
				vcg::tri::Allocator<CMeshO>::AddFace(cm, hole[0], thirdIdx, centerIdx);
				cm.face.back().C() = vcg::Color4b::Gray;
			}
		} else {
			// add face 12c; continue
			vcg::tri::Allocator<CMeshO>::AddFace(cm, firstIdx, secondIdx, centerIdx);
			cm.face.back().C() = vcg::Color4b::Gray;
			if (i == hole.size() - 2) {
				// add face 23c
				vcg::tri::Allocator<CMeshO>::AddFace(cm, secondIdx, thirdIdx, centerIdx);
				cm.face.back().C() = vcg::Color4b::Gray;
			}
		}
	}

	return;
}

void FillHoleRingByRingRefined(CMeshO& cm, std::vector<int> hole, float startAvgEdge, Point3m holeCenter, bool stepByStep, std::vector<float> vZChange, float adjustRatio)
{
	// if (hole.size() == 0) 
	// {
	// 	return;
	// }

    // Point3m centerPoint = findHoleCenterPoint(cm, hole);
	Point3m centerPoint = holeCenter;

	// float avgCenterDistance = calcAvgDistanceToCenter(cm, hole, centerPoint);
	// float startAvgEdge = CalcAvgHoleEdge(cm, hole);
	// float factor = avgCenterDistance / startAvgEdge;

	// float centerZChange = calcCenterZChange(cm, centerPoint, startAvgEdge, hole, vZChange);
	// assert(centerZChange);
    // centerPoint.Z() = centerPoint.Z() + centerZChange;

	while (true)
	{
		// Point3m centerPoint = findHoleCenterPoint(cm, hole);
		// float centerZChange = calcCenterZChange(cm, centerPoint, startAvgEdge, hole, vRatio);
    	// centerPoint.Z() = centerPoint.Z() + centerZChange;
		if (hole.size() < 3) {
			break;
		}

		hole = reduceHoleByConnectNearby(cm, hole, startAvgEdge, centerPoint);
		hole = rearrangeHole(cm, centerPoint, startAvgEdge, hole);

		float avgEdge = CalcAvgHoleEdge(cm, hole);
		
		// init check step to center
		int maxStepCenter = findMaxStepToCenter(cm, centerPoint, avgEdge, hole);
		
		// if (avgEdge < threshold)
		if (checkHoleSize(cm, hole, startAvgEdge, centerPoint) || maxStepCenter <= 1)
		{
			// fillHoleByCenter(cm, hole, avgEdge, 1);
            FillHoleByCenterRefined(cm, hole, avgEdge, 1, centerPoint);
			break;
		}

		std::vector<int> reducedHole;

		// insert first point by using last boundary point
		Point3m firstFillPoint = calcFillingPoint(cm.vert[hole.back()].P(), centerPoint, avgEdge, 1);
		CMeshO::VertexIterator firstFillV = vcg::tri::Allocator<CMeshO>::AddVertices(cm, 1);
		N_LOG_FILL_VERT(QString("LogFirstFill During fill Border Vertex index %1 next index %2 coord (x, y, z): (%3, %4, %5) \n").arg(
			QString::number(hole.back()), QString::number(firstFillV->Index()), QString::number(cm.vert[hole.back()].P().X()), QString::number(cm.vert[hole.back()].P().Y()), 
			QString::number(cm.vert[hole.back()].P().Z())
			));
		firstFillV->P() = firstFillPoint;
		firstFillV->C() = vcg::Color4b(0, 255, 255, 255);
		int firstFillIndex = firstFillV->Index();

		int prevIndex = hole.back();
		int prevFilledIndex = firstFillV->Index();
		bool skipPrev = false;
		for (int i = 0; i < hole.size(); i++)
		{
			int index = hole[i];
			// float ratio = vRatio[i]; // TODO: remove
			// if (cm.vert[index].P().Z() > centerPoint.Z()) {
			// 	cm.vert[index].C() = vcg::Color4b::Red;
			// 	qDebug("border z %f center z %f", cm.vert[index].P().Z(), centerPoint.Z());
			// }

			Point3m fillPoint = calcFillingPoint(cm.vert[index].P(), centerPoint, avgEdge, 1);
			int fillIndex = -1;

			// check step to center
			float dRatio;
			int stepCenter = findStepToCenter(centerPoint, cm.vert[index].P(), avgEdge, dRatio);
			assert(stepCenter <= maxStepCenter);
			if (stepCenter < maxStepCenter) {
				// reducedHole.push_back(index);
				if (skipPrev) {
					prevFilledIndex = prevIndex;
				} else {
					vcg::tri::Allocator<CMeshO>::AddFace(cm, prevIndex, index, prevFilledIndex);
					cm.face.back().C() = vcg::Color4b::Yellow;
				}
				prevIndex = index;
				skipPrev = true;

				if (index == hole.back()) {
					vcg::tri::Allocator<CMeshO>::AddFace(cm, firstFillIndex, index, prevFilledIndex);
					cm.face.back().C() = vcg::Color4b::Yellow;
					reducedHole.push_back(firstFillIndex);
				} else {
					reducedHole.push_back(index);
				}

				continue;
			}

			// fill new point
			CMeshO::VertexIterator vi;
			if (index != hole.back()) {
				// add point to mesh
				vi = vcg::tri::Allocator<CMeshO>::AddVertices(cm, 1);
				N_LOG_FILL_VERT(QString("LogFillPoint During fill Border Vertex step %1 index %2 next index %3 coord (x, y, z): (%4, %5, %6) \n").arg(
					QString::number(stepCenter), QString::number(index), QString::number(vi->Index()),
					QString::number(cm.vert[index].P().X()), QString::number(cm.vert[index].P().Y()), QString::number(cm.vert[index].P().Z())
				));
				// fillPoint.Z() = cm.vert[index].P().Z() * 1.073515;
				// if (stepCenter >= 2) {
					// fillPoint.Z() += centerZChange * (dRatio - 1) * adjustRatio;
					// qDebug("LogRatio During fill Border Vertex index %f \n", centerZChange * (dRatio - 1) * adjustRatio);
				// }
				// float czChange = fillPoint.Z() - cm.vert[index].P().Z();
				// fillPoint.Z() += czChange * adjustRatio;
				vi->P() = fillPoint;
				vi->C() = vcg::Color4b(0, 255, 255, 255);
				fillIndex = vi->Index();

				// check distance with new and prev fill points
				// if (!checkNewPrevPointDistance(fillPoint, cm.vert[prevFilledIndex].P(), threshold)) {
				// 	vcg::tri::Allocator<CMeshO>::AddFace(cm, prevFilledIndex, index, prevIndex);
				// 	prevIndex = index;
				// 	cm.vert[index].C() = vcg::Color4b(255, 255, 0, 255);
				// 	continue;
				// }
			} else {
				fillIndex = firstFillIndex;
			}
			
			if (checkCurrPointEdgeOk(fillPoint, cm.vert[index].P(), cm.vert[prevFilledIndex].P(), cm.vert[prevIndex].P())) {
				vcg::tri::Allocator<CMeshO>::AddFace(cm, prevIndex, index, fillIndex);
				cm.face.back().C() = vcg::Color4b::Green;
				vcg::tri::Allocator<CMeshO>::AddFace(cm, prevIndex, prevFilledIndex, fillIndex);
				cm.face.back().C() = vcg::Color4b::Green;
			} else {
				vcg::tri::Allocator<CMeshO>::AddFace(cm, prevIndex, index, prevFilledIndex);
				cm.face.back().C() = vcg::Color4b::Magenta;
				vcg::tri::Allocator<CMeshO>::AddFace(cm, prevFilledIndex, index, fillIndex);
				cm.face.back().C() = vcg::Color4b::Magenta;
			}

			if (skipPrev) {
				reducedHole.pop_back();
				skipPrev = false;
			}

			reducedHole.push_back(fillIndex);

			prevIndex = index;
			prevFilledIndex = fillIndex;
		}

		hole = reducedHole;

		if (stepByStep) {
			break;
		}
	}
	
	return;
}
