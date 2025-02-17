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

#include "filter_test_filling_hole.h"
#include "algorithm.h"

#include <QtGui>
#include <QMessageBox>
#include <QGridLayout>
#include <QSpacerItem>

#include <tuple>
#include <cmath>

#include <vcg/complex/algorithms/hole.h>
#include <vcg/space/color4.h>


// #define _N_DEBUG_FIND_VERT
#ifdef _N_DEBUG_FIND_VERT
#define N_LOG_FIND_VERT(...) qDebug(__VA_ARGS__)
#else
#define N_LOG_FIND_VERT(...)
#endif

#define _N_DEBUG_RUN_ALGO
#ifdef _N_DEBUG_RUN_ALGO
#define N_LOG_RUN_ALGO(...) qDebug(__VA_ARGS__)
#else
#define N_LOG_RUN_ALGO(...)
#endif


using namespace std;
using namespace vcg;

// our logic
// rotate mesh
void Freeze(MeshModel *m)
{
	tri::UpdatePosition<CMeshO>::Matrix(m->cm, m->cm.Tr,true);
	tri::UpdateBounding<CMeshO>::Box(m->cm);
	m->cm.shot.ApplyRigidTransformation(m->cm.Tr);
	m->cm.Tr.SetIdentity();
}

void ApplyTransform(MeshDocument &md, const Matrix44m &tr, bool toAllFlag, bool freeze,
					bool invertFlag=false, bool composeFlage=true)
{
	if(toAllFlag)
	{
		MeshModel   *m=NULL;
		while ((m=md.nextVisibleMesh(m)))
		{
			if(invertFlag) m->cm.Tr = Inverse(m->cm.Tr);
			if(composeFlage) m->cm.Tr = tr * m->cm.Tr;
			else m->cm.Tr=tr;
			if(freeze) Freeze(m);
		}

		for (int i = 0; i < md.rasterList.size(); i++)
			if (md.rasterList[0]->visible)
				md.rasterList[i]->shot.ApplyRigidTransformation(tr);
	}
	else
	{
		MeshModel   *m=md.mm();
		if(invertFlag) m->cm.Tr = Inverse(m->cm.Tr);
		if(composeFlage) m->cm.Tr = tr * m->cm.Tr;
		else m->cm.Tr=tr;
		if(freeze) Freeze(md.mm());
	}
}

Matrix44m rotateHoleCenter(MeshDocument &md, Point3m holeCenter) {
	Matrix44m trRot, trTran, trTranInv, transfM;
	// Point3m axis, tranVec;

	Point3m tranVec(0, 0, 0); // suppose holeCenter is P(x, y, z) tranVec is vector OP 
	Point3m zAxis(0, 0, 1);
	Point3m tranAxis = holeCenter ^ zAxis; // tranAxis is cross product of OP and Oz axis

	Scalarm angleRad= Angle(zAxis, holeCenter);
	Scalarm angleDeg = angleRad * (180.0 / 3.141592653589793238463);

	trRot.SetRotateDeg(angleDeg,tranAxis);
	trTran.SetTranslate(tranVec);
	trTranInv.SetTranslate(-tranVec);
	transfM = trTran*trRot*trTranInv;

	ApplyTransform(md, transfM, false, true);
	return transfM;
}

void rotateInverse(MeshDocument &md, Matrix44m mt) {
	Matrix44m imt = Inverse(mt);

	ApplyTransform(md, imt, false, true);
}

std::vector<int> getExpandedVertexIndex(CVertexO* vp) {
	std::vector<CVertexO*> vExpVert;
    vcg::face::VVStarVF<CFaceO>(vp, vExpVert);
	std::vector<int> vExpVertIdx;
	for (CVertexO* v : vExpVert) {
		vExpVertIdx.push_back(v->Index());
	}

	return vExpVertIdx;
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
                act->setIcon(QIcon(":/images/3d-icon-1.png"));
                //				act->setPriority(QAction::HighPriority);
            }
        }
    }
}

QString FilterFillHolePlugin::pluginName() const
{
	return "FilterTestCloseHole";
}

FilterFillHolePlugin::FilterClass FilterFillHolePlugin::getClass(const QAction * a) const
{
	switch (ID(a))
	{
	case FP_TEST_DP_CLOSE_HOLE: return FilterPlugin::NTest;
	case FP_TEST_CLOSE_HOLE: return FilterPlugin::NTest;

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
		return "Test ML close hole";
	case FP_TEST_CLOSE_HOLE:
		return "Test filling hole on mesh";
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
		return "Test close hole with Meshlab algorithm";
	case FP_TEST_CLOSE_HOLE:
		return "Our proposed method for filling hole on mesh";
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
		return MeshModel::MM_GEOMETRY_AND_TOPOLOGY_CHANGE + MeshModel::MM_TRANSFMATRIX;
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
		algoType.push_back("Fill large hole ring by ring");
		algoType.push_back("Fill small hole by centroid point");
		algoType.push_back("Rotate hole centroid to Oz");
		algoType.push_back("Get hole information");
		algoType.push_back("Test no rotate");
		algoType.push_back("Test slope");
		parlst.addParam(RichEnum("algo", 0, algoType, tr("Algorithm type"), tr("Choose the algorithm to close hole")));
		// parlst.addParam(RichPoint3f("hole_centroid", Point3f(0,0,0), "Centroid point", "Centroid point of a hole for running step by step"));
        // parlst.addParam(RichFloat("expect_edge_length", 0, "Boundary edge length", "If value = 0 then expected edge length is the average hole edges"));
        // parlst.addParam(RichFloat("threshold_ratio", 1.0, "Threshold ratio", "Threshold = Threhold ratio * Expected edge length"));
		parlst.addParam(RichDynamicFloat("adjust_ratio", 0.0f, 0.0f, 1.0f, "Adjustment z-coordinate ratio", "Adjustment ratio to adjust z-coordinate when filling new points"));
        parlst.addParam(RichInt("max_hole_size", int(0), "Max hole size", "Size of a hole is the number of boundary face of that hole"));
		// optional flags
		parlst.addParam(RichBool("selected", m.cm.sfn>0, "Apply algorithm with selected faces", "Only the holes with at least one of the boundary faces selected are applied"));
		// parlst.addParam(RichBool("one_ring", false, "Apply algorithm to fill one ring only", "Fill the hole one ring to the center"));
        parlst.addParam(RichBool("prevent_rotate", false, "Prevent rotate the hole", "Prevent rotate the mesh to move hole centroid to z-axis"));
		parlst.addParam(RichBool("enable_border_color", false, "Change color of boundary faces and boundary vertexes"));
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
		m.updateDataMask(MeshModel::MM_FACECOLOR);
		m.updateDataMask(MeshModel::MM_FACEFACETOPO);
		m.updateDataMask(MeshModel::MM_VERTCOLOR);
		m.updateDataMask(MeshModel::MM_VERTFACETOPO);

		CMeshO & cm = m.cm;
		if (tri::Clean<CMeshO>::CountNonManifoldEdgeFF(m.cm) > 0)
		{
			throw MLException("Mesh has some not 2-manifold edges, filter requires edge manifoldness");
		}

		// get parameters
		bool enableBorderColor = par.getBool("enable_border_color");
		// Point3m ucentroid = par.getPoint3m("hole_centroid");
		Point3m ucentroid(0, 0, 0);
		// float expectedEdgeLength = par.getFloat("expect_edge_length");
		float expectedEdgeLength = 0;
		// float thresholdRatio = par.getFloat("threshold_ratio");
		float thresholdRatio = 1;
		float adjustRatio = par.getFloat("adjust_ratio");
		// float adjustRatio = 0;
        float maxHoleSize = par.getInt("max_hole_size");
		bool selected = par.getBool("selected"); 
		// bool isOneRing = par.getBool("one_ring");
		bool isOneRing = false;
		bool preventRotate = par.getBool("prevent_rotate");

		std::vector<tri::Hole<CMeshO>::Info> vinfo;
		int borderVBit = CVertexO::NewBitFlag();
		int expandedVBit = CVertexO::NewBitFlag();
        std::vector<std::vector<int>> vholeIdx;
		std::vector<HoleVertInfo> vHoleVertInfo;
		std::vector<std::vector<HoleVertData>> vHoleData;
		std::vector<std::vector<int>> vHoleFaceIndex;

		// clear Visited flag for all face in mesh
		tri::UpdateFlags<CMeshO>::FaceClearV(cm);

        for(tri::Hole<CMeshO>::FaceIterator fi = cm.face.begin(); fi!=cm.face.end(); ++fi)
        {
			if((*fi).IsD())
			{
				continue;
			}

			if(selected && !(*fi).IsS())
			{
				//if I have to consider only the selected triangles e
				//what I'm considering isn't marking it and moving on
				(*fi).SetV();
				continue;
			}

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

					N_LOG_FIND_VERT("Hole %i detected \n", vinfo.size() + 1);
					std::vector<int> vBorderIndex;
					std::vector<float> vDistanceVert;
					Point3m prevPoint = sp.v->P();
					bool hasPrevP = true;
					std::vector<int> vFaceIndex;					
					float ratio = 0;
					HoleVertInfo vertInfo;
					std::vector<HoleVertData> vHoleVertData;

					tri::Hole<CMeshO>::Box3Type hbox;
					hbox.Add(sp.v->cP());

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
						sp.v->SetUserBit(borderVBit);
						CVertexO* expandedFP = sp.f->V2(sp.z);

						if (!expandedFP->IsUserBit(expandedVBit)) {
							// expandedFP->C() = vcg::Color4b(255, 255, 0, 255);
						}

						N_LOG_FIND_VERT("Border Vertex index %i coord (x, y, z): (%f, %f, %f) index %d \n", 
							sp.v->Index(), sp.v->P().X(), sp.v->P().Y(), sp.v->P().Z(), sp.v->Index());
						N_LOG_FIND_VERT("Extended Vertex index %i coord (x, y, z): (%f, %f, %f) index %d \n", 
							expandedFP->Index(), expandedFP->P().X(), expandedFP->P().Y(), expandedFP->P().Z(), expandedFP->Index());

						float zChange = sp.v->P().Z() - expandedFP->P().Z();

						int vIndex = sp.v->Index();

						// qDebug("Border Vertex index %i coord (x, y, z): (%f, %f, %f) \n", vIndex, sp.v->P().X(), sp.v->P().Y(), sp.v->P().Z());
						vBorderIndex.push_back(sp.v->Index());
						vFaceIndex.push_back(sp.f->Index());
						vertInfo.vHoleVertIndex.push_back(sp.v->Index());
						vertInfo.vExpHoleVertIndex.push_back(expandedFP->Index());
						// vertInfo.vZChange.push_back(zChange);
						std::vector<int> vExpVertIdx = getExpandedVertexIndex(sp.v);
						vertInfo.vSetExpVertIndex.push_back(vExpVertIdx);
						HoleVertData vData = { vIndex, vExpVertIdx, expandedFP->P() };
						vHoleVertData.push_back(vData);

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

					qDebug("End hole point log");
					vholeIdx.push_back(vBorderIndex);
					vHoleVertInfo.push_back(vertInfo);
					vHoleData.push_back(vHoleVertData);
					vHoleFaceIndex.push_back(vFaceIndex);

					//I recovered the information on the whole hole
					vinfo.push_back( tri::Hole<CMeshO>::Info(sp,holesize,hbox) );
				}
			}//end for edges of the triangle

		}//end for faces of the mesh!!!

		if (enableBorderColor)
		{
			for (std::vector<int> hole: vholeIdx) 
			{
				if (maxHoleSize > 0 && hole.size() > maxHoleSize) {
					continue;
				}

				for (int i = 0; i < hole.size(); i++)
				{
                    int vi = hole[i];
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

		// bool stop = false;
		switch(par.getEnum("algo"))
		{
			case 0: {
				for (HoleVertInfo hole: vHoleVertInfo) {

					std::vector<int> vVertIndex = hole.vHoleVertIndex;

                    if (maxHoleSize > 0 && vVertIndex.size() > maxHoleSize) {
						continue;
					}

					Point3m holeCenter = CalcHoleCenter(cm, vVertIndex);
					log("Center hole before rotate (x, y, z) = (%f, %f, %f)", holeCenter.X(), holeCenter.Y(), holeCenter.Z());
					// rotate the hole center to z-axis
					Matrix44m transMt = rotateHoleCenter(md, holeCenter);

					holeCenter = CalcHoleCenter(cm, vVertIndex); // re-calc after rotate
					log("Center hole after rotate (x, y, z) = (%f, %f, %f)", holeCenter.X(), holeCenter.Y(), holeCenter.Z());

					// compute z-change
					std::vector<float> vZChange2;
					for (int i = 0; i < vVertIndex.size(); i++) {
						int bIdx = vVertIndex[i];

						Point3m avgExpPoint(0, 0, 0);
						int count = 0;
						std::vector<int> vExpVertIdx = hole.vSetExpVertIndex[i];
						if (vExpVertIdx.size() < 2) {
							assert(false);
                        } else if (vExpVertIdx.size() == 2) {
                            avgExpPoint += (cm.vert[vExpVertIdx[0]].P() + cm.vert[vExpVertIdx[1]].P());
							count = 2;
						} else {
							for (int expIdx: vExpVertIdx) {
								// cm.vert[expIdx].C() = vcg::Color4b::Red;
								if (cm.vert[expIdx].IsUserBit(borderVBit)) {
									continue;
								}

								avgExpPoint += cm.vert[expIdx].P();
								++count;
							}
						}
                        // assert(count);

                        if (count == 0) {
                            for (int expIdx: vExpVertIdx) {

                                avgExpPoint += cm.vert[expIdx].P();
                                ++count;
                            }
                        }
						// if (count == 0) {
						// 	cm.vert[bIdx].C() = vcg::Color4b::Red;
						// 	stop = true;
						// 	break;
						// }
						avgExpPoint /= count;
						float zChange = cm.vert[bIdx].P().Z() - avgExpPoint.Z();
						// assert(zChange == zChange);
						vZChange2.push_back(zChange);
						hole.vExpPoint.push_back(avgExpPoint);
						// Point3m bPoint = cm.vert[bIdx].P();
						// float dxy = sqrt( pow(avgExpPoint.X() - bPoint.X(), 2) + pow(avgExpPoint.Y() - bPoint.Y(), 2) );
						// float sl = zChange / dxy;
						// hole.vSlope.push_back(sl);
					}
                    // assert(hole.vZChange.size() == hole.vSlope.size());
					// hole.vZChange = hole.vSlope;
					// hole.vZChange = vZChange2;

					float edgeLength = expectedEdgeLength;
					if (edgeLength <= 0) {
						edgeLength = CalcAvgHoleEdge(cm, vVertIndex);
					}
					// float holeThreshold = edgeLength * thresholdRatio;
					// if (holeThreshold <= 0) {
                    //     holeThreshold = edgeLength;
					// } 
					bool stepByStep = isOneRing;
                    // qDebug("start one hole filling with threshold %f", holeThreshold);
					// if (stop) {
					// 	break;
					// }

					// float centerZChange = calcCenterZChange(cm, holeCenter, edgeLength, vVertIndex, hole.vZChange);
					float centerZChange = CalcCenterZChangeUsingAvgExpVertex(cm, holeCenter, edgeLength, vVertIndex, hole.vExpHoleVertIndex);
					holeCenter.Z() = holeCenter.Z() + centerZChange;
					log("Center hole after change z (x, y, z) = (%f, %f, %f) and dz %f", holeCenter.X(), holeCenter.Y(), holeCenter.Z(), centerZChange);

                    FillHoleRingByRingRefined(cm, vVertIndex, edgeLength, holeCenter, stepByStep, centerZChange, adjustRatio);

					// revert the rotation
					rotateInverse(md, transMt);
					vcg::tri::UpdateNormal<CMeshO>::PerVertexNormalizedPerFaceNormalized(m.cm);
				}

			} break;

			case 1: 
			{
				for (std::vector<int> hole: vholeIdx) 
				{
					if (maxHoleSize > 0 && hole.size() > maxHoleSize) {
						continue;
					}

					// rotate the hole center to z-axis
					Point3m holeCenter = CalcHoleCenter(cm, hole);
					rotateHoleCenter(md, holeCenter);

					// fill hole
                    FillHoleByCenter(cm, hole, 0, 1);
				}
			} break;

			case 2:
			{
				for (std::vector<int> hole: vholeIdx) 
				{
					if (maxHoleSize > 0 && hole.size() > maxHoleSize) {
						continue;
					}

					// rotate the hole center to z-axis
					Point3m holeCenter = CalcHoleCenter(cm, hole);
					rotateHoleCenter(md, holeCenter);
				}
			} break;

			case 3:
			{
				QMessageBox msgBox;
				msgBox.setText("Holes information");
				QString text = "";

				int count = 0;
				for (std::vector<int> hole: vholeIdx) 
				{
					if (maxHoleSize > 0 && hole.size() > maxHoleSize) {
						continue;
					}

					float avgDistance = CalcAvgHoleEdge(cm, hole);
                    Point3m holeCenter = CalcHoleCenter(cm, hole);

					++count;
					log("Hole number %i size %i average edge %f centroid point C(%f, %f, %f)", count, hole.size(), avgDistance, holeCenter.X(), holeCenter.Y(), holeCenter.Z());
					text.append(QString("Hole number %1 size %2 average edge %3 centroid point C(%4, %5, %6) \n").arg(
                        QString::number(count), QString::number(hole.size()), QString::number(avgDistance),
						QString::number(holeCenter.X()), QString::number(holeCenter.Y()), QString::number(holeCenter.Z())
						));
				}

				msgBox.setInformativeText(text);
				QSpacerItem* horizontalSpacer = new QSpacerItem(500, 0, QSizePolicy::Minimum, QSizePolicy::Expanding);
    			QGridLayout* layout = (QGridLayout*)msgBox.layout();
    			layout->addItem(horizontalSpacer, layout->rowCount(), 0, 1, layout->columnCount());
				msgBox.exec();
			} break;
			
			case 4:
			{
				for (HoleVertInfo hole: vHoleVertInfo) {

					std::vector<int> vVertIndex = hole.vHoleVertIndex;

                    if (maxHoleSize > 0 && vVertIndex.size() > maxHoleSize) {
						continue;
					}

					Point3m holeCenter = CalcHoleCenter(cm, vVertIndex);

					// compute z-change
					std::vector<Point3m> vChange;
					for (int i = 0; i < vVertIndex.size(); i++) {
						int bIdx = vVertIndex[i];

						Point3m avgExpPoint(0, 0, 0);
						int count = 0;
						std::vector<int> vExpVertIdx = hole.vSetExpVertIndex[i];
						if (vExpVertIdx.size() < 2) {
							assert(false);
                        } else if (vExpVertIdx.size() == 2) {
                            avgExpPoint += (cm.vert[vExpVertIdx[0]].P() + cm.vert[vExpVertIdx[1]].P());
							count = 2;
						} else {
							for (int expIdx: vExpVertIdx) {
								// cm.vert[expIdx].C() = vcg::Color4b::Red;
								if (cm.vert[expIdx].IsUserBit(borderVBit)) {
									continue;
								}

								avgExpPoint += cm.vert[expIdx].P();
								++count;
							}
						}
                        // assert(count);

                        if (count == 0) {
                            for (int expIdx: vExpVertIdx) {

                                avgExpPoint += cm.vert[expIdx].P();
                                ++count;
                            }
                        }
						// if (count == 0) {
						// 	cm.vert[bIdx].C() = vcg::Color4b::Red;
						// 	stop = true;
						// 	break;
						// }
						avgExpPoint /= count;
						Point3m cChange = cm.vert[bIdx].P() - avgExpPoint;
						// assert(zChange == zChange);
						// vZChange2.push_back(zChange);
						hole.vExpPoint.push_back(avgExpPoint);
						// Point3m bPoint = cm.vert[bIdx].P();
						// float dxy = sqrt( pow(avgExpPoint.X() - bPoint.X(), 2) + pow(avgExpPoint.Y() - bPoint.Y(), 2) );
						// float sl = zChange / dxy;
						// hole.vSlope.push_back(sl);
					}
                    // assert(hole.vZChange.size() == hole.vSlope.size());
					// hole.vZChange = hole.vSlope;
					// hole.vZChange = vZChange2;

					float edgeLength = expectedEdgeLength;
					if (edgeLength <= 0) {
						edgeLength = CalcAvgHoleEdge(cm, vVertIndex);
					}
					// float holeThreshold = edgeLength * thresholdRatio;
					// if (holeThreshold <= 0) {
                    //     holeThreshold = edgeLength;
					// } 
					bool stepByStep = isOneRing;
                    // qDebug("start one hole filling with threshold %f", holeThreshold);
					// if (stop) {
					// 	break;
					// }

					Point3m centerChange = CalcCenterChange(cm, holeCenter, edgeLength, vVertIndex, hole.vExpPoint);
					holeCenter = holeCenter + centerChange;

                    FillHoleRingByRingRefined(cm, vVertIndex, edgeLength, holeCenter, stepByStep, 0.0, adjustRatio);

					qDebug("End one hole filling");
				}
			} break;

			case 5:
			{
				for (HoleVertInfo hole: vHoleVertInfo) {

					std::vector<int> vVertIndex = hole.vHoleVertIndex;

                    if (maxHoleSize > 0 && vVertIndex.size() > maxHoleSize) {
						continue;
					}

					Point3m holeCenter = CalcHoleCenter(cm, vVertIndex);

					// compute slope
					// compute z-change
					// compute hole center
					// fill hole

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
