#ifndef EXTRADECORATEPLUGIN_H
#define EXTRADECORATEPLUGIN_H

#include <common/plugins/interfaces/decorate_plugin.h>
#include <common/ml_shared_data_context/ml_shared_data_context.h>

// #include <wrap/gui/coordinateframe.h>

class DecorateMeshPlugin : public QObject, public DecoratePlugin
{
    Q_OBJECT
	MESHLAB_PLUGIN_IID_EXPORTER(DECORATE_PLUGIN_IID)
	Q_INTERFACES(DecoratePlugin)
	QString decorationName(ActionIDType filter) const;
	QString decorationInfo(ActionIDType filter) const;
	QString pluginName() const;

    vcg::Color4b textColor;

    enum {
		DP_SHOW_HOLE
	};

public:
    DecorateMeshPlugin()
    {
        typeList <<
            DP_SHOW_HOLE;

        ActionIDType tt;
        foreach(tt, types())
        {
            actionList << new QAction(decorationName(tt), this);
        }
        QAction *ap;
        foreach(ap, actionList) {
            ap->setCheckable(true);
        }
    }

    void initGlobalParameterList(const QAction* /*format*/, RichParameterList& /*globalparam*/);

	bool startDecorate(const QAction *, MeshDocument &, const RichParameterList *, GLArea *);
	bool startDecorate(const QAction *, MeshModel &, const RichParameterList *, GLArea *);
	void decorateMesh(const QAction *, MeshModel &, const RichParameterList *, GLArea *, QPainter *, GLLogStream &);
	void decorateDoc(const QAction *, MeshDocument &, const RichParameterList *, GLArea *, QPainter *, GLLogStream &);
	void endDecorate(const QAction *, MeshModel &, const RichParameterList *, GLArea *);
	// void endDecorate(const QAction *, MeshDocument &, const RichParameterList *, GLArea *);

	bool isDecorationApplicable(const QAction *, const MeshModel&, QString&) const;

	int getDecorationClass(const QAction *) const;

	void DrawVertLabel(std::vector<CVertexO*> &vv, QPainter *gla);

    inline QString MeshHoleNormalLength() const { return  "MeshLab::Decoration::Mesh::Hole::NormalLength" ; }
	inline QString MeshHoleNormalWidth() const { return  "MeshLab::Decoration::Mesh::Hole::NormalWidth" ; }
	inline QString MeshHoleNormalVertColor() const { return  "MeshLab::Decoration::Mesh::Hole::NormalVertColor" ; }
	inline QString MeshHoleNormalFaceColor() const { return  "MeshLab::Decoration::Mesh::Hole::NormalFaceColor" ; }
	inline QString MeshHoleNormalSelection() const { return  "MeshLab::Decoration::Mesh::Hole::NormalSelection"; }
	inline QString MeshHoleNormalVertFlag() const { return  "MeshLab::Decoration::Mesh::Hole::NormalVertFlag" ; }
	inline QString MeshHoleNormalFaceFlag() const { return  "MeshLab::Decoration::Mesh::Hole::NormalFaceFlag" ; }
	
};

#endif
