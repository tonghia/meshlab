#ifndef FILTER_TESTBEZIER_H
#define FILTER_TESTBEZIER_H

#include <QObject>
#include <common/plugins/interfaces/filter_plugin.h>

// Octrees
#include <vcg/space/index/octree.h>
#include <vcg/space/index/octree_template.h>
#include <wrap/utils.h>

// Points-Matrixes
#include <vcg/math/matrix33.h>
#include <vcg/space/deprecated_point3.h>
#include <vcg/space/box3.h>

using namespace std ;
using namespace vcg ;
using namespace tri ;


// /////////////////////// TestFilterPlugin

class TestBezierPlugin : public QObject, public FilterPlugin
{
    Q_OBJECT
		MESHLAB_PLUGIN_IID_EXPORTER(FILTER_PLUGIN_IID)
    Q_INTERFACES(FilterPlugin)

public:
       /* naming convention :
                 - FP -> Filter Plugin
                 - name of the plugin separated by _
       */
  enum {
    TEST_BEZIER
  } ;
	TestBezierPlugin();

	QString pluginName() const;
	virtual QString filterInfo(ActionIDType filter) const;
	virtual QString filterName(ActionIDType filter) const;

	void initParameterList(
		const QAction* action,
		MeshModel &m,
		RichParameterList &parlst);

	QString filterScriptFunctionName(ActionIDType) ;
	std::map<std::string, QVariant> applyFilter(
		const QAction *action,
		const RichParameterList & /*parameters*/,
		MeshDocument &md,
		unsigned int& postConditionMask,
		vcg::CallBackPos * cb) ;
	int getRequirements(const QAction*);
	int getPreConditions(QAction*) const;
	int postCondition(const QAction* ) const;
	FilterClass getClass(const QAction *a) const;
	FilterArity filterArity(const QAction* filter) const;

private:
    // Octree class

    typedef Octree<CVertexO, float > MyMeshOctree;

    // MyPoint3

    typedef Point3<float> MyPoint3 ;
    typedef Matrix33<float> MyMatrix33 ;

    void print_coords(CMeshO::CoordType) ;
    void print_matrix33(MyMatrix33) ;
};

#endif // FILTER_TESTBEZIER_H
