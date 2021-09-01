#ifndef FILTERFILLHOLE_PLUGIN_H
#define FILTERFILLHOLE_PLUGIN_H

#include <common/plugins/interfaces/filter_plugin.h>

class FilterFillHolePlugin : public QObject, public FilterPlugin
{
	Q_OBJECT
	MESHLAB_PLUGIN_IID_EXPORTER(FILTER_PLUGIN_IID)
	Q_INTERFACES(FilterPlugin)

public:
	enum { FP_TEST_DP_CLOSE_HOLE, FP_TEST_CLOSE_HOLE } ;

	FilterFillHolePlugin();

	QString pluginName() const;
	QString filterName(ActionIDType filter) const;
	QString filterInfo(ActionIDType filter) const;

	FilterClass getClass(const QAction*) const;
	void initParameterList(const QAction*, MeshModel &/*m*/, RichParameterList & /*parent*/);
	std::map<std::string, QVariant> applyFilter(
			const QAction* action,
			const RichParameterList & parameters,
			MeshDocument &md,
			unsigned int& postConditionMask,
			vcg::CallBackPos * cb);
	int postCondition(const QAction *action) const;
	int getPreConditions(const QAction *action) const;
	FilterArity filterArity(const QAction *) const {return SINGLE_MESH;}

};

#endif
