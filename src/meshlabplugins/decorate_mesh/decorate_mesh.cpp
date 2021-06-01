#include "decorate_mesh.h"
// #include <wrap/gl/addons.h>
// #include <vcg/complex/algorithms/stat.h>
// #include <vcg/complex/algorithms/bitquad_support.h>
#include <common/GLExtensionsManager.h>
#include <meshlab/glarea.h>
#include <wrap/qt/checkGLError.h>
#include <wrap/qt/gl_label.h>
// #include <QGLShader>
#include <meshlab/glarea_setting.h>
#include <wrap/gl/gl_type_name.h>
//
#include <vcg/complex/algorithms/hole.h>
#include <vcg/space/color4.h>
//
using namespace vcg;
using namespace std;


QString DecorateMeshPlugin::decorationInfo(ActionIDType filter) const
{
	switch(filter)
	{
        case DP_SHOW_HOLE: return tr("Draw hole normals");
	}
	assert(0);
	return QString();
}

QString DecorateMeshPlugin::pluginName() const
{
	return "DecorateMesh";
}

QString DecorateMeshPlugin::decorationName(ActionIDType filter) const
{
	switch(filter)
	{
    case DP_SHOW_HOLE: return QString("Show Mesh Hole");
	}
	return QString("error!");
}

void DecorateMeshPlugin::decorateDoc(const QAction* a, MeshDocument &md, const RichParameterList *rm, GLArea *gla, QPainter *painter, GLLogStream &/*_log*/)
{
	QFont qf;
	
	switch (ID(a))
	{
		
	} // end switch
}

void DecorateMeshPlugin::decorateMesh(const QAction* a, MeshModel &m, const RichParameterList *rm, GLArea *gla, QPainter *painter, GLLogStream &_log)
{
	this->setLog(&_log);
	QFont qf;
	
	textColor = rm->getColor4b( GLAreaSetting::textColorParam());
	
	glPushMatrix();
	glMultMatrix(m.cm.Tr);
	switch (ID(a))
	{
		
	case DP_SHOW_HOLE:
	{
		std::vector<CVertexO*> vBoundaryVertex;

		glPushAttrib(GL_ENABLE_BIT );
		Scalarm NormalLen=rm->getFloat(MeshHoleNormalLength());
		Scalarm NormalWid = rm->getFloat(MeshHoleNormalWidth());
		vcg::Color4b VertNormalColor = rm->getColor4b(MeshHoleNormalVertColor());
		vcg::Color4b FaceNormalColor = rm->getColor4b(MeshHoleNormalFaceColor());
		bool showselection = rm->getBool(MeshHoleNormalSelection());
		
		Scalarm LineLen = m.cm.bbox.Diag()*NormalLen;
		
		//query line width range
		GLfloat widthRange[2];
		widthRange[0] = 1.0f; widthRange[1] = 1.0f;
		glGetFloatv(GL_ALIASED_LINE_WIDTH_RANGE, widthRange);
		
		//set line width according to the width range
		NormalWid = (NormalWid < widthRange[0]) ?  widthRange[0] :  NormalWid;
		NormalWid = (NormalWid > widthRange[1]) ?  widthRange[1] :  NormalWid;
		
		//store current linewidth and set new line width
		GLfloat lineWidthtmp[1];
		glGetFloatv(GL_LINE_WIDTH, lineWidthtmp);
		glLineWidth(NormalWid);
		
		glDisable(GL_LIGHTING);
		glDisable(GL_TEXTURE_2D);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
		glBegin(GL_LINES);

        if(rm->getBool(MeshHoleNormalVertFlag())) // vert Normals
		{
			glColor(VertNormalColor);
			// for(CMeshO::VertexIterator vi=m.cm.vert.begin();vi!=m.cm.vert.end();++vi) if(!(*vi).IsD())
			// {
			// 	if ((!showselection) || (showselection && vi->IsS()))
			// 	{
			// 		glVertex((*vi).P());
			// 		glVertex((*vi).P() + (*vi).N()*LineLen);
			// 	}
			// }
			// typename std::vector<tri::Hole<CMeshO>::Info>::iterator ith;

            // for(auto vi = vvi.begin(); vi!= vvi.end(); ++vi)
			// {
			// 	if ((!showselection) || (showselection && vi->IsS()))
			// 	{
			// 		glVertex((*vi).P());
			// 		glVertex((*vi).P() + (*vi).N()*LineLen);
			// 	}
			// }
		}
        if(rm->getBool(MeshHoleNormalFaceFlag())) // face Normals
		{
			glColor(FaceNormalColor);
            // for(auto fi=m.cm.face.begin();fi!=m.cm.face.end();++fi) if(!(*fi).IsD())
			// {
			// 	if ((!showselection) || (showselection && fi->IsS()))
			// 	{
			// 		Point3m b = Barycenter(*fi);
			// 		glVertex(b);
			// 		glVertex(b + (*fi).N()*LineLen);
			// 	}
			// }
		}

		// Start Find hole
		m.updateDataMask(MeshModel::MM_FACEFACETOPO);
		m.updateDataMask(MeshModel::MM_VERTCOLOR);
		m.updateDataMask(MeshModel::MM_FACECOLOR);
		CMeshO & cm = m.cm;

		vector<tri::Hole<CMeshO>::Info> vinfo;
		bool Selected = showselection;

		// vector<CMeshO::FaceIterator> vfi;
		// vector<CMeshO::VertexIterator> vvi;

		tri::UpdateFlags<CMeshO>::FaceClearV(cm);
        for(CMeshO::FaceIterator fi = cm.face.begin(); fi!=cm.face.end(); ++fi)
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

								tri::Hole<CMeshO>::Box3Type hbox;
								hbox.Add(sp.v->cP());
								//printf("Looping %i : (face %i edge %i) \n", VHI.size(),sp.f-&*m.face.begin(),sp.z);
                                qDebug("Looping %i : (face %i edge %i) \n", vinfo.size(),sp.f-&*cm.face.begin(),sp.z);
								sp.f->SetV();
								do
								{
									sp.f->SetV();
									hbox.Add(sp.v->cP());
									++holesize;
									sp.NextB();
									sp.f->SetV();

									// Set corlor of border face, vertex
									sp.f->C().SetHSVColor(0, 1.0f, 1.0f);
									sp.v->C() = vcg::Color4b(255, 0, 255, 255);

									// vfi.push_back(fi);
                                    // vvi.push_back(sp.v);
									vBoundaryVertex.push_back(sp.v);

                                    // glVertex((*sp.v).P());
                                    // glVertex((*sp.v).P() + (*sp.v).N()*LineLen);

									assert(sp.IsBorder());
								}while(sp != fp);

								//I recovered the information on the whole hole
                                // vinfo.push_back( tri::Hole<CMeshO>::Info(sp,holesize,hbox) );
							}
						}//for on the edges of the triangle
				}//S & !S
			}//!IsD()
		}//for principale!!!
		// End Find hole

		for(CVertexO* vi: vBoundaryVertex)
		{
			glVertex((*vi).P());
            glVertex((*vi).P() + (*vi).N()*LineLen);
		}
		
		glEnd();
		//restore previous line width
		glLineWidth(lineWidthtmp[0]);
		glPopAttrib();

        DrawVertLabel(vBoundaryVertex, painter);

	} break;
		
	} // end switch;
	glPopMatrix();
}


void DecorateMeshPlugin::DrawVertLabel(std::vector<CVertexO*> &vv,QPainter *painter)
{
	glPushAttrib(GL_LIGHTING_BIT | GL_CURRENT_BIT | GL_DEPTH_BUFFER_BIT );
	glDepthFunc(GL_ALWAYS);
	glDisable(GL_LIGHTING);
	glColor3f(.4f,.4f,.4f);
    for(CVertexO* vi: vv)
	{
		if(!(*vi).IsD())
			glLabel::render(painter, (*vi).P(),tr("%1").arg(vi->Index()),glLabel::Mode(textColor));
	}
	glPopAttrib();
}

int DecorateMeshPlugin::getDecorationClass(const QAction *action) const
{
	switch(ID(action))
	{
    case DP_SHOW_HOLE: return DecorateMeshPlugin::PerMesh;
	}
	assert (0);
	return 0;
}

bool DecorateMeshPlugin::isDecorationApplicable(const QAction* action, const MeshModel& m, QString &ErrorMessage) const
{
	
	return true;
}


bool DecorateMeshPlugin::startDecorate(const QAction * action, MeshDocument &, const RichParameterList *, GLArea *)
{
	switch(ID(action))
	{
		
	}
	return true;
}


void DecorateMeshPlugin::endDecorate(const QAction * action, MeshModel &m, const RichParameterList *, GLArea *)
{
	switch(ID(action))
	{

	default: break;
	}
}

bool DecorateMeshPlugin::startDecorate(const QAction * action, MeshModel &m, const RichParameterList *rm, GLArea *gla)
{
	switch(ID(action))
	{
	
	}
	return true;
}

void DecorateMeshPlugin::initGlobalParameterList(const QAction* action, RichParameterList &parset)
{
	
	switch(ID(action))
    {

    case DP_SHOW_HOLE:
	{
		parset.addParam(RichFloat(MeshHoleNormalLength(),0.05f,"Vector Length","The length of the normal expressed as a percentage of the bbox of the mesh"));
		parset.addParam(RichFloat(MeshHoleNormalWidth(), 1.0f,"Normal Width","The width of the normal expressed in pixels"));
		parset.addParam(RichColor(MeshHoleNormalVertColor(),QColor(102, 102, 255, 153),QString("Curr Vert Normal Color"),QString("Current Vert Normal Color")));
		parset.addParam(RichColor(MeshHoleNormalFaceColor(),QColor(102, 102, 255, 153),QString("Curr Face Normal Color"),QString("Current Face Normal Color")));
		parset.addParam(RichBool(MeshHoleNormalVertFlag(),true,"Per Vertex",""));
		parset.addParam(RichBool(MeshHoleNormalFaceFlag(),true,"Per Face",""));
		parset.addParam(RichBool(MeshHoleNormalSelection(), false, "Show Selected", ""));
    } break;

	}
	
}


MESHLAB_PLUGIN_NAME_EXPORTER(DecorateMeshPlugin)
