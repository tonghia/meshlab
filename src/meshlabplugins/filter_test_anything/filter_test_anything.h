#ifndef filter_test_anything_H
#define filter_test_anything_H

#include <QObject>
#include <common/plugins/interfaces/filter_plugin.h>


// /////////////////////// TestFilterPlugin

class TestAnythingPlugin : public QObject, public FilterPlugin
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
		TEST_ANYTHING,
		TEST_FILL_COLOR,
		TEST_FILL_FACE_COLOR
	};

	TestAnythingPlugin();

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

};

#endif // filter_test_anything_H
