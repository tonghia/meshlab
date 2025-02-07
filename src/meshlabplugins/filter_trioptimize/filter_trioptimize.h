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

#ifndef TRIOPTIMIZEFILTERSPLUGIN_H
#define TRIOPTIMIZEFILTERSPLUGIN_H

#include <QObject>
#include <common/plugins/interfaces/filter_plugin.h>

class TriOptimizePlugin : public QObject, public FilterPlugin
{
	Q_OBJECT
	MESHLAB_PLUGIN_IID_EXPORTER(FILTER_PLUGIN_IID)
	Q_INTERFACES(FilterPlugin)

public:
	enum { 
		// mesh improvement by edge flipping
		FP_CURVATURE_EDGE_FLIP,
		FP_PLANAR_EDGE_FLIP,
		// Laplacian smooth that do not moves vertices far from the surface
		FP_NEAR_LAPLACIAN_SMOOTH
	};

	TriOptimizePlugin();
	
	QString pluginName() const;
	QString filterName(ActionIDType filter) const;
	QString filterInfo(ActionIDType filter) const;
	void initParameterList(const QAction*, MeshModel &/*m*/, RichParameterList & /*parent*/);
	std::map<std::string, QVariant> applyFilter(
			const QAction* action,
			const RichParameterList & parameters,
			MeshDocument &md,
			unsigned int& postConditionMask,
			vcg::CallBackPos * cb);
	int getRequirements(const QAction*);
	FilterClass getClass(const QAction *) const;
	int postCondition(const QAction* ) const;
	FilterArity filterArity(const QAction *) const {return SINGLE_MESH;}

};

#endif
