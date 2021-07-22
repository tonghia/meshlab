#include "filter_test_anything.h"

#include <vector>
#include <iostream>
#include <ctime>
#include <stdlib.h>
#include <QtGui>
#define _YES_I_WANT_TO_USE_DANGEROUS_STUFF
#include <vcg/math/old_lin_algebra.h>
#include <vcg/space/colorspace.h>
#include <wrap/io_trimesh/io_mask.h>
#include <vcg/math/histogram.h>
#include <vcg/complex/algorithms/update/color.h>

using namespace vcg;

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

TestAnythingPlugin::TestAnythingPlugin()
{
	typeList = {
		TEST_ANYTHING,
		TEST_FILL_COLOR,
		TEST_FILL_FACE_COLOR,
		TEST_ROTATE
	};

	for(ActionIDType tt : types()) {
			QAction* act = new QAction(filterName(tt), this);
			actionList.push_back(act);
	}
}

QString TestAnythingPlugin::pluginName() const
{
	return "TestAnythingPlugin";
}

QString TestAnythingPlugin::filterName(ActionIDType filter) const
{
	switch (filter)
	{
	case TEST_ANYTHING:
		return QString("Test Anything");
	case TEST_FILL_COLOR:
		return QString("Test Fill Color");
	case TEST_FILL_FACE_COLOR:
		return QString("Test Fill Face Color");
	case TEST_ROTATE:
		return QString("Test Rotate Mesh");
	}
	return QString("Unknown filter");
}

QString TestAnythingPlugin::filterInfo(ActionIDType filterId) const
{
	switch (filterId)
	{
	case TEST_ANYTHING:
		return tr("Test anything");
	case TEST_FILL_COLOR:
		return tr("Test fill color in code not from Qt Color picker parameter");
	case TEST_FILL_FACE_COLOR:
		return tr("Test fill face color");
	case TEST_ROTATE:
		return tr("Test rotate mesh");
	}
	return QString("Unknown Filter");
}

TestAnythingPlugin::FilterClass TestAnythingPlugin::getClass(const QAction *action) const
{
	switch (ID(action))
	{
	case TEST_ANYTHING:
		return FilterClass(FilterPlugin::PointSet + FilterPlugin::VertexColoring);
	case TEST_FILL_COLOR:
		return FilterClass(FilterPlugin::NTest);
	case TEST_FILL_FACE_COLOR:
		return FilterClass(FilterPlugin::NTest);
	case TEST_ROTATE:
		return FilterClass(FilterPlugin::NTest);
	default:
		assert(0);
	}
	return FilterPlugin::Generic;
}

void TestAnythingPlugin::initParameterList(const QAction *action, MeshModel &m, RichParameterList &parlst)
{
	switch (ID(action))
	{
        case TEST_ANYTHING: {
				// QColor color1 = QColor(0, 0, 0, 255);
				// parlst.addParam(RichColor("color1", color1, "Color:", "Sets the color to apply to vertices."));
				// parlst.addParam(RichBool("onSelected", false, "Only on selection", "If checked, only affects selected vertices"));
				break;
        }
        case TEST_FILL_COLOR: {
			break;
		}
		case TEST_FILL_FACE_COLOR:
		{
			break;
		}
		case TEST_ROTATE:
		{
			QStringList rotMethod;
			rotMethod.push_back("X axis");
			rotMethod.push_back("Y axis");
			rotMethod.push_back("Z axis");
			rotMethod.push_back("custom axis");
			parlst.addParam(RichEnum("rotAxis", 0, rotMethod, tr("Rotation on:"), tr("Choose a method")));
			QStringList rotCenter;
			rotCenter.push_back("origin");
			rotCenter.push_back("barycenter");
			rotCenter.push_back("custom point");
			parlst.addParam(RichEnum("rotCenter", 0, rotCenter, tr("Center of rotation:"), tr("Choose a method")));
			parlst.addParam(RichDynamicFloat("angle",0,-360,360,"Rotation Angle","Angle of rotation (in <b>degree</b>). If snapping is enabled this value is rounded according to the snap value"));
			parlst.addParam(RichPoint3f("customAxis",Point3f(0,0,0),"Custom axis","This rotation axis is used only if the 'custom axis' option is chosen."));
			parlst.addParam(RichPoint3f("customCenter",Point3f(0,0,0),"Custom center","This rotation center is used only if the 'custom point' option is chosen."));
			parlst.addParam(RichBool("snapFlag", false, "Snap angle", "If selected, before starting the filter will remove any unreferenced vertex (for which curvature values are not defined)"));
			parlst.addParam(RichFloat("snapAngle",30,"Snapping Value","This value is used to snap the rotation angle (i.e. if the snapping value is 30, 227 becomes 210)."));
			parlst.addParam(RichBool ("Freeze",true,"Freeze Matrix","The transformation is explicitly applied, and the vertex coordinates are actually changed"));
			parlst.addParam(RichBool ("allLayers",false,"Apply to all visible Layers","If selected the filter will be applied to all visible mesh layers"));

			break;
		}
	default:
		assert(0);
	}
}

// /////////////////////////////////////////////
// Main
// /////////////////////////////////////////////
std::map<std::string, QVariant> TestAnythingPlugin::applyFilter(
		const QAction *action,
		const RichParameterList &par,
		MeshDocument &md,
		unsigned int & /*postConditionMask*/,
		vcg::CallBackPos * /*cb*/)
{
	std::map<std::string, QVariant> outputValues;

	MeshModel *m = md.mm();
	CMeshO & cm = m->cm;

	switch(ID(action)) {
		case TEST_FILL_COLOR: {
			// First try to clone from filter_colorproc CP_FILLING
			// vcg::tri::UpdateColor<CMeshO>::PerVertexConstant(meshModel->cm, vcg::Color4b(vcg::Color4b::Red), false);

			// Second try to bring their code here
			//     bool selected = false;
			//     Color4b vs=Color4b::Red;
			//	for(auto vi=cm.vert.begin();vi!=cm.vert.end();++vi) {
			//		if(!(*vi).IsD()){
			//			if(!selected || (*vi).IsS())
			//			{
			//				(*vi).C() = Color4b(255, 255, 0, 255);
			//			}
			//		}
			//	}

			// Finally use our code
			int meshVertexNumber = cm.VN();
			for (int i = 1; i < meshVertexNumber; i++)
			{
						cm.vert[i].C() = vcg::Color4b(vcg::Color4b::Red);
						//  cm.vert[i].C() = Color4b(255, 255, 0, 255);
			}

			m->updateDataMask(MeshModel::MM_VERTCOLOR);
			break;
		}
		case TEST_FILL_FACE_COLOR:
		{
			m->updateDataMask(MeshModel::MM_FACECOLOR);

			// every parser variables is related to face attributes.
			for(CMeshO::FaceIterator fi = cm.face.begin(); fi != cm.face.end(); ++fi)
			{
				// set new color for this iteration
				(*fi).C() = Color4b(255, 0, 255, 255);
			}

			break;
		}
		case TEST_ROTATE:
		{
			Matrix44m trRot, trTran, trTranInv, transfM;
			Point3m axis, tranVec;

			switch(par.getEnum("rotAxis"))
			{
			case 0: axis=Point3m(1,0,0); break;
			case 1: axis=Point3m(0,1,0);break;
			case 2: axis=Point3m(0,0,1);break;
			case 3: axis=par.getPoint3m("customAxis");break;
			}
			switch(par.getEnum("rotCenter"))
			{
			case 0: tranVec=Point3m(0,0,0); break;
			case 1: tranVec= cm.Tr * cm.bbox.Center(); break;
			case 2: tranVec=par.getPoint3m("customCenter");break;
			}

			Scalarm angleDeg= par.getDynamicFloat("angle");
			Scalarm snapAngle = par.getFloat("snapAngle");
			if(par.getBool("snapFlag"))
			{
				angleDeg = floor(angleDeg / snapAngle)*snapAngle;
				//par.setValue("angle", DynamicFloatValue(angleDeg));
			}

			trRot.SetRotateDeg(angleDeg,axis);
			trTran.SetTranslate(tranVec);
			trTranInv.SetTranslate(-tranVec);
			transfM = trTran*trRot*trTranInv;

			ApplyTransform(md,transfM,par.getBool("allLayers"),par.getBool("Freeze"));

			break;
		}
	}

	return outputValues;
}

QString TestAnythingPlugin::filterScriptFunctionName(ActionIDType filterID)
{
	switch (filterID)
	{
	case TEST_ANYTHING:
		return QString("TestAnything");
	case TEST_FILL_COLOR:
		return QString("TestFillColor");
	case TEST_FILL_FACE_COLOR:
		return QString("TestFillFaceColor");
	case TEST_ROTATE:
		return QString("TestRotate");
	default:
		assert(0);
	}
	return QString();
}

int TestAnythingPlugin::getRequirements(const QAction *action)
{
	switch (ID(action))
	{
	case TEST_ANYTHING:
		return MeshModel::MM_VERTCOORD;
    case TEST_FILL_COLOR:
        return MeshModel::MM_NONE;
	case TEST_FILL_FACE_COLOR:
		return MeshModel::MM_NONE;
	case TEST_ROTATE:
		return MeshModel::MM_NONE;
	default:
		return 0;
	}
}

int TestAnythingPlugin::postCondition(const QAction *action) const
{
	switch (ID(action))
	{
	case TEST_ANYTHING:
		return MeshModel::MM_VERTCOLOR;
	case TEST_FILL_COLOR:
		return MeshModel::MM_NONE;
	case TEST_FILL_FACE_COLOR:
		return MeshModel::MM_NONE;
	case TEST_ROTATE:
		return MeshModel::MM_TRANSFMATRIX;
	default:
		return MeshModel::MM_UNKNOWN;
	}
}

int TestAnythingPlugin::getPreConditions(QAction *action) const
{
	switch (ID(action))
	{
	case TEST_ANYTHING:
		return MeshModel::MM_NONE;
	case TEST_FILL_COLOR:
		return MeshModel::MM_NONE;
	case TEST_FILL_FACE_COLOR:
		return MeshModel::MM_NONE;
	case TEST_ROTATE:
		return MeshModel::MM_NONE;
	default:
		return 0;
	}
}

FilterPlugin::FilterArity TestAnythingPlugin::filterArity(const QAction *filter) const
{
	switch(ID(filter))
	{
		case TEST_ANYTHING:
			return FilterPlugin::VARIABLE;
		case TEST_FILL_COLOR:
			return FilterPlugin::SINGLE_MESH;
		case TEST_FILL_FACE_COLOR:
			return FilterPlugin::SINGLE_MESH;
		case TEST_ROTATE:
			return FilterPlugin::SINGLE_MESH;
	}

	return FilterPlugin::VARIABLE;
}

MESHLAB_PLUGIN_NAME_EXPORTER(TestAnythingPlugin)
