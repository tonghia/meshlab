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

TestAnythingPlugin::TestAnythingPlugin()
{
	typeList = {
		TEST_ANYTHING,
		TEST_FILL_COLOR,
		TEST_FILL_FACE_COLOR
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
	}

	return FilterPlugin::VARIABLE;
}

MESHLAB_PLUGIN_NAME_EXPORTER(TestAnythingPlugin)
