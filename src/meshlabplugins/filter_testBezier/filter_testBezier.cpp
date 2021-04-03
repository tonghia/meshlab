#include "filter_testBezier.h"

#include <vector>
#include <iostream>
#include <ctime>
#include <stdlib.h>
#include <QtGui>
#define _YES_I_WANT_TO_USE_DANGEROUS_STUFF
#include <vcg/math/old_lin_algebra.h>
#include <vcg/space/color4.h>
#include <wrap/io_trimesh/io_mask.h>
#include <vcg/math/histogram.h>

#define MAX 999

#define k_max 1

#define EPSILON 0.000001f

#define norm_p 1
#define ex_bound 2
#define in_bound 3
#define ins_point 4

#define same_line 11 // any number

// DEBUG EXTRACT EXTERNAL BOUNDARY
//#define _M_DEBUG
#ifdef _M_DEBUG
#define M_LOG(...) qDebug(__VA_ARGS__)
#else
#define M_LOG(...)
#endif

// DEBUG EXTRACT INTERNAL BOUNDARY
//#define _M_DEBUG2
#ifdef _M_DEBUG2
#define M_LOG2(...) qDebug(__VA_ARGS__)
#else
#define M_LOG2(...)
#endif

// DEBUG FILLING THE HOLE
//#define _M_DEBUG3
#ifdef _M_DEBUG3
#define M_LOG3(...) qDebug(__VA_ARGS__)
#else
#define M_LOG3(...)
#endif

// DEBUG TANGENT PLANE (WHEN FILLING THE HOLE)
//#define _M_DEBUG4
#ifdef _M_DEBUG4
#define M_LOG4(...) qDebug(__VA_ARGS__)
#else
#define M_LOG4(...)
#endif

// Body

TestBezierPlugin::TestBezierPlugin()
{
	typeList = {
		TEST_BEZIER
	};

    QCoreApplication* app = QCoreApplication::instance();
    for(ActionIDType tt : types()) {
        QAction* act = new QAction(filterName(tt), this);
        actionList.push_back(act);

        if (app != nullptr) {
            if(tt==TEST_BEZIER){
//				act->setShortcut(QKeySequence ("Ctrl+Del"));
                act->setIcon(QIcon(":/images/TestBezier.png"));
//				act->setPriority(QAction::HighPriority);
            }
        }
    }
}

QString TestBezierPlugin::pluginName() const
{
	return "TestBezierPlugin";
}

QString TestBezierPlugin::filterName(ActionIDType filter) const
{
	switch (filter)
	{
	case TEST_BEZIER:
		return QString("Test Bezier points");
	}
	return QString("Unknown filter");
}

QString TestBezierPlugin::filterInfo(ActionIDType filterId) const
{
	switch (filterId)
	{
	case TEST_BEZIER:
		return tr("Create points and apply Bezier curves");
	}
	return QString("Unknown Filter");
}

TestBezierPlugin::FilterClass TestBezierPlugin::getClass(const QAction *action) const
{
	switch (ID(action))
	{
	case TEST_BEZIER:
		return FilterClass(FilterPlugin::PointSet + FilterPlugin::VertexColoring);
	default:
		assert(0);
	}
	return FilterPlugin::Generic;
}

void TestBezierPlugin::initParameterList(const QAction *action, MeshModel &m, RichParameterList &parlst)
{
	switch (ID(action))
	{
	case TEST_BEZIER:
		parlst.addParam(RichBool("case1", true, "External Boundary Extraction", "Case 1"));
		parlst.addParam(RichBool("case2", true, "Internal Boundary Extraction", "Case 2"));
		parlst.addParam(RichBool("case3", true, "Filling the Holes", "Case 3"));
		parlst.addParam(RichBool("case4", true, "Refine points of the Holes", "Case 4"));
		//parlst.addParam(new RichInt("x", 3, "X", "x variable"));
		//parlst.addParam(new RichInt("y", 4, "Y", "y variable"));
		//parlst.addParam(new RichInt("k", 2, "K", "k variable"));
		break;
	default:
		assert(0);
	}
}

//========================================================
// Calculate Bezier

struct Point
{
	float x;
	float y;
	float z;
};
Point Points[MAX][MAX];

// the level of detail of the surface
unsigned int LOD = 20;

// debug by print out value x,y,z of point
void M_debug_point(CVertexO a)
{
	M_LOG4("Point at (%d,%d): %f, %f, %f", (int)a.P().X(), (int)a.P().Y(), a.P().X(), a.P().Y(), a.P().Z());
}

Point CalculateU(float t, int row)
{

	// the final point
	Point p;

	// the t value inverted
	float it = 1.0f - t;

	// calculate blending functions
	float b0 = t * t * t;
	float b1 = 3 * t * t * it;
	float b2 = 3 * t * it * it;
	float b3 = it * it * it;

	// sum the effects of the Points and their respective blending functions
	p.x = b0 * Points[row][0].x +
				b1 * Points[row][1].x +
				b2 * Points[row][2].x +
				b3 * Points[row][3].x;

	p.y = b0 * Points[row][0].y +
				b1 * Points[row][1].y +
				b2 * Points[row][2].y +
				b3 * Points[row][3].y;

	p.z = b0 * Points[row][0].z +
				b1 * Points[row][1].z +
				b2 * Points[row][2].z +
				b3 * Points[row][3].z;

	return p;
}

Point CalculateV(float t, Point *pnts)
{
	Point p;

	// the t value inverted
	float it = 1.0f - t;

	// calculate blending functions
	float b0 = t * t * t;
	float b1 = 3 * t * t * it;
	float b2 = 3 * t * it * it;
	float b3 = it * it * it;

	// sum the effects of the Points and their respective blending functions
	p.x = b0 * pnts[0].x +
				b1 * pnts[1].x +
				b2 * pnts[2].x +
				b3 * pnts[3].x;

	p.y = b0 * pnts[0].y +
				b1 * pnts[1].y +
				b2 * pnts[2].y +
				b3 * pnts[3].y;

	p.z = b0 * pnts[0].z +
				b1 * pnts[1].z +
				b2 * pnts[2].z +
				b3 * pnts[3].z;

	return p;
}

Point Calculate(float u, float v)
{

	// first of all we will need to evaluate 4 curves in the u
	// direction. The points from those will be stored in this
	// temporary array
	Point temp[4];

	// calculate each point on our final v curve
	temp[0] = CalculateU(u, 0);
	temp[1] = CalculateU(u, 1);
	temp[2] = CalculateU(u, 2);
	temp[3] = CalculateU(u, 3);

	// having got 4 points, we can use it as a bezier curve
	// to calculate the v direction. This should give us our
	// final point
	//
	return CalculateV(v, temp);
}

// check 4 neighborhood
bool miss_1_of_4_neigh(CVertexO **vertices, int x, int y, int extend_y, int max_x, int max_y, int min_x, int min_y, int k)
{
	//return 4 neightborhood of point(x,y) are axist or not

	if (x == min_x || y == min_y)
	{
		M_LOG("x,y(%d,%d) at external boundary --> Select this point", x, y);
		return true;
	}

	if (x == max_x || y == max_y)
	{
		M_LOG("x,y(%d,%d) at external boundary --> Select this point", x, y);
		return true;
	}

	if (vertices[(x + k + 1) * extend_y + (y + k)] == NULL || // Check right == NULL
			vertices[(x + k - 1) * extend_y + (y + k)] == NULL || // Check left == NULL
			vertices[(x + k) * extend_y + (y + k - 1)] == NULL || // Check bottom == NULL
			vertices[(x + k) * extend_y + (y + k + 1)] == NULL)		// Check top == NULL
	{
		M_LOG("x,y(%d,%d) miss at least 1 neightborhood --> Select this point", x, y);
		return true;
	}
	else
	{
		M_LOG("x,y(%d,%d) have 4 neightborhood --> Reject this point", x, y);
		return false;
	}
}

// check 8 neighborhoods exits internal boundary or not
bool check_interal_8_neigh(int *vertices, int x, int y, int extend_y, int k)
{

	bool i = false;
	if (vertices[(x + k + 1) * extend_y + (y + k)] == in_bound)
	{
		i = true;
	}
	else if (vertices[(x + k - 1) * extend_y + (y + k)] == in_bound)
	{
		i = true;
	}
	else if (vertices[(x + k) * extend_y + (y + k - 1)] == in_bound)
	{
		i = true;
	}
	else if (vertices[(x + k) * extend_y + (y + k + 1)] == in_bound)
	{
		i = true;
	}

	else if (vertices[(x + k - 1) * extend_y + (y + k - 1)] == in_bound)
	{
		i = true;
	}
	else if (vertices[(x + k - 1) * extend_y + (y + k + 1)] == in_bound)
	{
		i = true;
	}
	else if (vertices[(x + k + 1) * extend_y + (y + k - 1)] == in_bound)
	{
		i = true;
	}
	else if (vertices[(x + k + 1) * extend_y + (y + k + 1)] == in_bound)
	{
		i = true;
	}

	return i;
}

// check 8 neighborhoods have at least one normal point or not
bool no_norm_8_neigh(int *vertices, int x, int y, int extend_y, int k)
{

	bool i = true;
	if (vertices[(x + k + 1) * extend_y + (y + k)] == norm_p)
	{
		i = false;
	}
	else if (vertices[(x + k - 1) * extend_y + (y + k)] == norm_p)
	{
		i = false;
	}
	else if (vertices[(x + k) * extend_y + (y + k - 1)] == norm_p)
	{
		i = false;
	}
	else if (vertices[(x + k) * extend_y + (y + k + 1)] == norm_p)
	{
		i = false;
	}

	else if (vertices[(x + k - 1) * extend_y + (y + k - 1)] == norm_p)
	{
		i = false;
	}
	else if (vertices[(x + k - 1) * extend_y + (y + k + 1)] == norm_p)
	{
		i = false;
	}
	else if (vertices[(x + k + 1) * extend_y + (y + k - 1)] == norm_p)
	{
		i = false;
	}
	else if (vertices[(x + k + 1) * extend_y + (y + k + 1)] == norm_p)
	{
		i = false;
	}

	return i;
}

// check 8 neighbordhood
void find_first_of_8_neigh(CVertexO **vertices, int x_curr, int y_curr, int &x, int &y, int max_y, int k)
{
	if (vertices[max_y * (x_curr + k) + (y_curr + k - 1)] != NULL)
	{
		x = x_curr;
		y = y_curr - 1;
	}
	else if (vertices[max_y * (x_curr + k - 1) + (y_curr + k - 1)] != NULL)
	{
		x = x_curr - 1;
		y = y_curr - 1;
	}
	else if (vertices[max_y * (x_curr + k - 1) + (y_curr + k)] != NULL)
	{
		x = x_curr - 1;
		y = y_curr;
	}
	else if (vertices[max_y * (x_curr + k - 1) + (y_curr + k + 1)] != NULL)
	{
		x = x_curr - 1;
		y = y_curr + 1;
	}
	else if (vertices[max_y * (x_curr + k) + (y_curr + k + 1)] != NULL)
	{
		x = x_curr;
		y = y_curr + 1;
	}
	else if (vertices[max_y * (x_curr + k + 1) + (y_curr + k + 1)] != NULL)
	{
		x = x_curr + 1;
		y = y_curr + 1;
	}
	else if (vertices[max_y * (x_curr + k + 1) + (y_curr + k)] != NULL)
	{
		x = x_curr + 1;
		y = y_curr;
	}
	else if (vertices[max_y * (x_curr + k + 1) + (y_curr + k - 1)] != NULL)
	{
		x = x_curr + 1;
		y = y_curr - 1;
	}
}

// return index of one dimention array
int Cal_index_k(int x, int y, int max_y, int k)
{
	return (x + k) * max_y + (y + k);
}

// find the point match between 2 line
//----------------------------------------------------
void find_match(int x1, int y1,
								int x2, int y2,
								int x3, int y3,
								int x4, int y4,
								float &x_result, float &y_result,
								bool &status_intersection)
{
	// return x_result,y_result which is the intersection created by line (1-2) and (3-4)

	//for easy imagination, we alway use (x1,y1) is the p_curr
	float m1, c1, m2, c2;
	//float x1, y1, x2, y2;
	float dx1, dy1, dx2, dy2;

	float intersection_X, intersection_Y;

	M_LOG("Program to find the intersecting point of two lines:");

	M_LOG("the first line is (p_curr,p_prev) ");
	M_LOG("x1,y1(%d,%d) \t x2,y2(%d,%d) ", x1, y1, x2, y2);
	M_LOG("x3,y3(%d,%d) \t x4,y4(%d,%d) ", x3, y3, x4, y4);

	dx1 = x2 - x1;
	dy1 = y2 - y1;
	m1 = dy1 / dx1; // if x1 == x2 therefore, the linear is: X = x1 = x2
	// y = mx + c
	// intercept c = y - mx
	c1 = y1 - m1 * x1; // which is same as y2 - slope * x2

	dx2 = x4 - x3;
	dy2 = y4 - y3;
	m2 = dy2 / dx2;		 // if x3 == x4 therefore, the linear is: X = x3 = x4
	c2 = y4 - m2 * x4; // which is same as y2 - slope * x2

	if (dx1 == 0 && dx2 == 0)
	{
		//because 2 line created by (p_curr, p_prev) and side of square are never the same
		status_intersection = false;
		M_LOG("Parallel with y axis ");
		M_LOG("-------------");
	}
	else if (dy1 == 0 && dy2 == 0)
	{
		//because 2 line created by (p_curr, p_prev) and side of square are never the same
		status_intersection = false;
		M_LOG("Parallel with x axis ");
		M_LOG("-------------");
	}
	else if (dx1 == 0)
	{
		M_LOG("Equation of line1: X = %f ", (float)x1);
		M_LOG("Equation of line2: Y = %f*X %c %f ", m2, ((c2 < 0) ? ' ' : '+'), c2);

		intersection_X = (float)x1;
		intersection_Y = m2 * intersection_X + c2;
		M_LOG("Intersecting Point: = %f, %f ", intersection_X, intersection_Y);
		x_result = intersection_X;
		y_result = intersection_Y;
		status_intersection = true;
		M_LOG("-------------");
	}
	else if (dx2 == 0)
	{
		M_LOG("Equation of line1: Y = %f*X %c %f ", m1, ((c1 < 0) ? ' ' : '+'), c1);
		M_LOG("Equation of line2: X = %f ", (float)x3);

		intersection_X = (float)x3;
		intersection_Y = m1 * intersection_X + c1;
		M_LOG("Intersecting Point: = %f, %f ", intersection_X, intersection_Y);
		x_result = intersection_X;
		y_result = intersection_Y;
		status_intersection = true;
		M_LOG("-------------");
	}
	else if ((m1 - m2) == 0)
	{
		status_intersection = false;
		M_LOG("No Intersection between the lines");
	}
	else
	{
		M_LOG("Equation of line1: Y = %f*X %c %f ", m1, ((c1 < 0) ? ' ' : '+'), c1);
		//std::cout << m1 << "X " << ((c1 < 0) ? ' ' : '+') << c1 << "";
		M_LOG("Equation of line2: Y = %f*X %c %f ", m2, ((c2 < 0) ? ' ' : '+'), c2);
		//std::cout << m2 << "X " << ((c2 < 0) ? ' ' : '+') << c2 << "";

		intersection_X = (c2 - c1) / (m1 - m2);
		intersection_Y = m2 * intersection_X + c2;
		M_LOG("Intersecting Point: = %f, %f ", intersection_X, intersection_Y);
		x_result = intersection_X;
		y_result = intersection_Y;
		status_intersection = true;
		M_LOG("-------------");
	}
}

//find distance between 2 point
float distance(float x1, float y1, float x2, float y2)
{
	return sqrt(pow((x1 - x2), 2) + pow((y1 - y2), 2));
}

//find distance between 2 point in 3D
double distance3D(CVertexO p1, CVertexO p2)
{

	return sqrt(pow((p1.P().X() - p2.P().X()), 2) + pow((p1.P().Y() - p2.P().Y()), 2) + pow((p1.P().Z() - p2.P().Z()), 2));
}

// compare 2 float
bool Compare2Float(float A, float B)
{
	float diff = A - B;
	//M_LOG("diff: %f", diff);
	//M_LOG("EPSILON: %f", EPSILON);
	return (diff < EPSILON) && (-diff < EPSILON);
}

// compare 2 double
// return true when difference of 2 double < Epsilon
bool Compare2Double(double A, double B)
{
	double diff = A - B;
	//M_LOG("diff: %f", diff);
	//M_LOG("EPSILON: %f", EPSILON);
	return (diff < EPSILON) && (-diff < EPSILON);
}

//check whether if p_curr is the center point or not
bool p_curr_is_center(float x_center, float y_center, float x1, float y1, float x2, float y2)
{
	//return p_curr is the center point(1-2) or not
	float d1 = distance(x_center, y_center, x1, y1);
	//M_LOG("distance(x_center(%f), y_center(%f), x1(%f), y1(%f)): %f", x_center, y_center, x1, y1, d1);
	float d2 = distance(x_center, y_center, x2, y2);
	//M_LOG("distance(x_center(%f), y_center(%f), x2(%f), y2(%f)): %f", x_center, y_center, x2, y2, d2);
	float d3 = distance(x1, y1, x2, y2);
	//M_LOG("distance(x1(%f), y1(%f), x2(%f), y2(%f)): %f", x1, y1, x2, y2, d3);

	if (Compare2Float(d3, d1 + d2))
	{
		M_LOG("This point (%f,%f) is center (Point (%f,%f) & (%f,%f))", x_center, y_center, x1, y1, x2, y2);
		return true;
	}
	else
	{
		M_LOG("This point (%f,%f) is NOT center (Point (%f,%f) & (%f,%f))", x_center, y_center, x1, y1, x2, y2);
		return false;
	}
}

// find the intersection between RNk(p) and half_line(p_curr, p_prev)
void Comp_p_Prev_p_Curr(
		int x_prev, int y_prev,
		int x_curr, int y_curr,
		int max_y, int k_ins,
		float &x_result, float &y_result, int &line)
{
	// return x_result,y_result which is the intersection between RNk(p) and half_line(p_curr, p_prev)
	int x1 = x_curr - k_ins;
	int y1 = y_curr - k_ins;

	int x2 = x_curr - k_ins;
	int y2 = y_curr + k_ins;

	int x3 = x_curr + k_ins;
	int y3 = y_curr + k_ins;

	int x4 = x_curr + k_ins;
	int y4 = y_curr - k_ins;

	//find 4 intersection between (p_curr, p_prev) and sides of square RN(k)
	//***check on the x_line below (Point 1 & 2)***
	bool status_intersection;
	M_LOG("Start find intersection of (p_curr, p_prev) and 1-2");
	find_match(x_curr, y_curr, x_prev, y_prev, x1, y1, x2, y2, x_result, y_result, status_intersection);
	if (status_intersection)
	{
		bool check = p_curr_is_center(x_curr, y_curr, x_prev, y_prev, x_result, y_result);
		if (check)
		{
			M_LOG("p_curr(%d,%d) is center p_prev(%d,%d) and p_result(%f,%f)", x_curr, y_curr, x_prev, y_prev, x_result, y_result);
			M_LOG("->REJECT point(%f,%f)  ", x_result, y_result);
		}
		if ((int)x_result == x1 && (y2 >= y_result) && (y_result >= y1) && !check)
		{
			line = 12;
			M_LOG("This point(%f,%f) belong (Point 1(%d,%d) & 2(%d,%d))", x_result, y_result, x1, y1, x2, y2);
			return;
		}
	}

	//***check on the x_line below (Point 2 & 3)***
	M_LOG("Start find intersection of (p_curr, p_prev) and 2-3");
	find_match(x_curr, y_curr, x_prev, y_prev, x2, y2, x3, y3, x_result, y_result, status_intersection);
	if (status_intersection)
	{
		bool check = p_curr_is_center(x_curr, y_curr, x_prev, y_prev, x_result, y_result);
		if (check)
		{
			M_LOG("p_curr(%d,%d) is center p_prev(%d,%d) and p_result(%f,%f)", x_curr, y_curr, x_prev, y_prev, x_result, y_result);
			M_LOG("->REJECT point(%f,%f)  ", x_result, y_result);
		}
		if ((int)y_result == y2 && (x3 >= x_result) && (x_result >= x2) && !check)
		{
			line = 23;
			M_LOG("This point(%f,%f) belong (Point 2(%d,%d) & 3(%d,%d))", x_result, y_result, x2, y2, x3, y3);
			return;
		}
	}

	//***check on the x_line below (Point 3 & 4)***
	M_LOG("Start find intersection of (p_curr, p_prev) and 3-4");
	find_match(x_curr, y_curr, x_prev, y_prev, x3, y3, x4, y4, x_result, y_result, status_intersection);
	if (status_intersection)
	{
		bool check = p_curr_is_center(x_curr, y_curr, x_prev, y_prev, x_result, y_result);
		if (check)
		{
			M_LOG("p_curr(%d,%d) is center p_prev(%d,%d) and p_result(%f,%f)", x_curr, y_curr, x_prev, y_prev, x_result, y_result);
			M_LOG("->REJECT point(%f,%f)  ", x_result, y_result);
		}
		if ((int)x_result == x3 && (y3 >= y_result) && (y_result >= y4) && !check)
		{
			line = 34;
			M_LOG("This point(%f,%f) belong (Point 3(%d,%d) & 4(%d,%d))", x_result, y_result, x3, y3, x4, y4);
			return;
		}
	}

	//***check on the x_line below (Point 4 & 1)***
	M_LOG("Start find intersection of (p_curr, p_prev) and 4-1");
	find_match(x_curr, y_curr, x_prev, y_prev, x4, y4, x1, y1, x_result, y_result, status_intersection);
	if (status_intersection)
	{
		bool check = p_curr_is_center(x_curr, y_curr, x_prev, y_prev, x_result, y_result);
		if (check)
		{
			M_LOG("p_curr(%d,%d) is center p_prev(%d,%d) and p_result(%f,%f)", x_curr, y_curr, x_prev, y_prev, x_result, y_result);
			M_LOG("->REJECT point(%f,%f)  ", x_result, y_result);
		}
		if ((int)y_result == y4 && (x4 >= x_result) && (x_result >= x1) && !check)
		{
			line = 41;
			M_LOG("This point(%f,%f) belong (Point 4(%d,%d) & 1(%d,%d))", x_result, y_result, x4, y4, x1, y1);
		}
	}
}

// Move next point on the RNk(p)
bool move_next_RNk(int &x_RNk, int &y_RNk, int x_curr, int y_curr, int k)
{
	//check RNK(p) is on which sides of square
	//return next position

	// 4 coners
	int x1 = x_curr - k;
	int y1 = y_curr - k;

	int x2 = x_curr - k;
	int y2 = y_curr + k;

	int x3 = x_curr + k;
	int y3 = y_curr + k;

	int x4 = x_curr + k;
	int y4 = y_curr - k;

	M_LOG("\t Prepare to move to next point on RNk");
	M_LOG("Check x_RNk,y_RNk (%d, %d) on which side of square", x_RNk, y_RNk);

	//check on side 1-2
	M_LOG("Check x_RNk,y_RNk (%d, %d) on line 1-2", x_RNk, y_RNk);
	if (p_curr_is_center(x_RNk, y_RNk, x1, y1, x2, y2) && (y_RNk < y2))
	{
		M_LOG("End check (%d, %d) and move to next point (%d, %d)", x_RNk, y_RNk, x_RNk, y_RNk + 1);
		y_RNk++;
		return true;
	}
	//check on side 2-3
	M_LOG("Check x_RNk,y_RNk (%d, %d) on line 2-3", x_RNk, y_RNk);
	if (p_curr_is_center(x_RNk, y_RNk, x2, y2, x3, y3) && (x_RNk < x3))
	{
		M_LOG("End check (%d, %d) and move to next point (%d, %d)", x_RNk, y_RNk, x_RNk + 1, y_RNk);
		x_RNk++;
		return true;
	}
	//check on side 3-4
	M_LOG("Check x_RNk,y_RNk (%d, %d) on line 3-4", x_RNk, y_RNk);
	if (p_curr_is_center(x_RNk, y_RNk, x3, y3, x4, y4) && (y_RNk > y4))
	{
		M_LOG("End check (%d, %d) and move to next point (%d, %d)", x_RNk, y_RNk, x_RNk, y_RNk - 1);
		y_RNk--;
		return true;
	}
	//check on side 4-1
	M_LOG("Check x_RNk,y_RNk (%d, %d) on line 4-1", x_RNk, y_RNk);
	if (p_curr_is_center(x_RNk, y_RNk, x4, y4, x1, y1) && (x_RNk > x1))
	{
		M_LOG("End check (%d, %d) and move to next point (%d, %d)", x_RNk, y_RNk, x_RNk - 1, y_RNk);
		x_RNk--;
		return true;
	}

	qDebug("ERROR at move_next_RNk() ");
	qDebug("ERROR at (%d,%d)", x_curr, y_curr);
	return false;
}

// check float the same with int
bool check_float(float f)
{
	if (f - (int)f == 0)
		return true;
	else
		return false;
}

int round_num(float a)
{
	if (a >= 0)
		return (int)(a + 0.5);
	else
		return (int)(a - 0.5);
}

void round_by_line(float x, float y, int line, int &x_new, int &y_new)
{
	if (line == 12)
	{
		x_new = (int)x;
		y_new = (int)y;
	}
	else if (line == 23)
	{
		x_new = (int)x;
		y_new = (int)y;
	}
	else if (line == 34)
	{
		x_new = (int)x;
		if ((y - (int)y) > EPSILON)
			y_new = (int)y + 1;
		else
			y_new = (int)y;
	}
	else if (line == 41)
	{
		if ((x - (int)x) > EPSILON)
			x_new = (int)x + 1;
		else
			x_new = (int)x;
		y_new = (int)y;
	}
}

// check whether a point has checked or not
bool not_exits_point(CVertexO *vertices, int total, int x, int y)
{
	M_LOG("start check number of total checked points");
	if (total == 0)
		return true;
	M_LOG("start check array checked points");
	for (int i = 0; i < total; i++)
	{
		//M_LOG("Old points: i(%d) x(%f)-y(%f)", i, vertices[i].P().X(), vertices[i].P().Y());
		if ((vertices[i].P().X() == x) && (vertices[i].P().Y() == y))
		{
			return false;
		}
	}
	M_LOG("not exits");
	return true;
}

// set a point into array
void set_point(CVertexO *vertices, int &index, int x, int y)
{
	CVertexO a;
	a.P() = vcg::Point3f(x, y, 0.0);
	vertices[index] = a;
	index++;
}

// set first point for inter hole
void first_inter_point(int &min_x, int &min_y, int total_inter_points, vector<CVertexO> points_list)
{
	min_x = points_list[0].P().X();
	min_y = points_list[0].P().Y();
	//M_LOG2("debug min point, min_x (%d) - min_y(%d)", min_x, min_y);

	for (int i = 0; i < total_inter_points; i++)
	{
		if (points_list[i].P().X() <= min_x)
		{
			min_x = points_list[i].P().X();
			min_y = points_list[i].P().Y(); // set temp y
		}
	}

	// find min y when equal x
	for (int i = 0; i < total_inter_points; i++)
	{
		if (points_list[i].P().X() == min_x)
		{
			if (points_list[i].P().Y() <= min_y)
            {
                min_y = points_list[i].P().Y();
				//M_LOG2("debug inside min point, min_x (%d) - min_y(%d)", min_x, min_y);
			}
		}
	}
}

// erase one iterior in vector list
bool erase(vector<CVertexO> &myNumbers_in, CVertexO number_in)
{
	bool found = false;
	for (unsigned int i = 0; i < myNumbers_in.size(); i++)
	{
		if ((myNumbers_in[i].P().X() == number_in.P().X()) &&
				(myNumbers_in[i].P().Y() == number_in.P().Y()))
		{
			myNumbers_in.erase(myNumbers_in.begin() + i);
			M_LOG2("delete: i(%d) at x(%f) y(%f)", i, number_in.P().X(), number_in.P().Y());
			found = true;
			break;
		}
	}
	if (!found)
	{
		M_LOG2("No point in list to delete: x(%f) y(%f)", number_in.P().X(), number_in.P().Y());
	}
	return found;
}

// Bubble Sort Function for Descending Order, abc = 1: x; abc = 0 y
void BubbleSort(vector<CVertexO> &point, int abc)
{
	int i, j, flag = 1; // set flag to 1 to start first pass
	CVertexO temp;			// holding variable
	int numLength = point.size();
	for (i = 1; (i <= numLength) && flag; i++)
	{
		flag = 0;
		for (j = 0; j < (numLength - 1); j++)
		{
			if (abc == 1)
			{
				if (point[j + 1].P().X() < point[j].P().X()) // ascending order simply changes to <
				{
					temp = point[j]; // swap elements
					point[j] = point[j + 1];
					point[j + 1] = temp;
					flag = 1; // indicates that a swap occurred.
				}
			}
			if (abc == 0)
			{
				if (point[j + 1].P().Y() < point[j].P().Y()) // ascending order simply changes to <
				{
					temp = point[j]; // swap elements
					point[j] = point[j + 1];
					point[j + 1] = temp;
					flag = 1; // indicates that a swap occurred.
				}
			}
		}
	}
	return; //arrays are passed to functions by address; nothing is returned
}

// find the interpolate of Bezier Curve through 4 points
int interpolate(CVertexO p0, double u, CVertexO p1, double v, CVertexO p2, CVertexO p3, CVertexO &c1, CVertexO &c2)
{
	double a = 0.0, b = 0.0, c = 0.0, d = 0.0, det = 0.0;
	CVertexO q1, q2;

	if ((u <= 0.0) || (u >= 1.0) || (v <= 0.0) || (v >= 1.0) || (u >= v))
	{
		return 0; /* failure */
		M_LOG3("Error at calculate u,v");
	}

	a = 3 * (1 - u) * (1 - u) * u;
	b = 3 * (1 - u) * u * u;
	c = 3 * (1 - v) * (1 - v) * v;
	d = 3 * (1 - v) * v * v;
	//M_LOG3("a: %f", a);
	//M_LOG3("b: %f", b);
	//M_LOG3("c: %f", c);
	//M_LOG3("d: %f", d);
	det = a * d - b * c;
	//M_LOG3("det: %f", det);
	/* unnecessary, but just in case... */
	if (det == 0.0)
	{
		return 0; /* failure */
		M_LOG3("Error at calculate det = a*d - b*c");
	}

	q1.P().X() = p1.P().X() - ((1 - u) * (1 - u) * (1 - u) * p0.P().X() + u * u * u * p3.P().X());
	q1.P().Y() = p1.P().Y() - ((1 - u) * (1 - u) * (1 - u) * p0.P().Y() + u * u * u * p3.P().Y());
	q1.P().Z() = p1.P().Z() - ((1 - u) * (1 - u) * (1 - u) * p0.P().Z() + u * u * u * p3.P().Z());

	q2.P().X() = p2.P().X() - ((1 - v) * (1 - v) * (1 - v) * p0.P().X() + v * v * v * p3.P().X());
	q2.P().Y() = p2.P().Y() - ((1 - v) * (1 - v) * (1 - v) * p0.P().Y() + v * v * v * p3.P().Y());
	q2.P().Z() = p2.P().Z() - ((1 - v) * (1 - v) * (1 - v) * p0.P().Z() + v * v * v * p3.P().Z());

	c1.P().X() = d * q1.P().X() - b * q2.P().X();
	c1.P().Y() = d * q1.P().Y() - b * q2.P().Y();
	c1.P().Z() = d * q1.P().Z() - b * q2.P().Z();
	c1.P().X() /= det;
	c1.P().Y() /= det;
	c1.P().Z() /= det;

	c2.P().X() = (-c) * q1.P().X() + a * q2.P().X();
	c2.P().Y() = (-c) * q1.P().Y() + a * q2.P().Y();
	c2.P().Z() = (-c) * q1.P().Z() + a * q2.P().Z();
	c2.P().X() /= det;
	c2.P().Y() /= det;
	c2.P().Z() /= det;

	return 1; /* success */
}

// find a point on Bezier curve
CVertexO drawBezier_curve(CVertexO A, CVertexO B, CVertexO C, CVertexO D, double t)
{
	CVertexO P;
	P.P().X() = pow((1 - t), 3) * A.P().X() + 3 * t * pow((1 - t), 2) * B.P().X() + 3 * (1 - t) * pow(t, 2) * C.P().X() + pow(t, 3) * D.P().X();
	P.P().Y() = pow((1 - t), 3) * A.P().Y() + 3 * t * pow((1 - t), 2) * B.P().Y() + 3 * (1 - t) * pow(t, 2) * C.P().Y() + pow(t, 3) * D.P().Y();
	P.P().Z() = pow((1 - t), 3) * A.P().Z() + 3 * t * pow((1 - t), 2) * B.P().Z() + 3 * (1 - t) * pow(t, 2) * C.P().Z() + pow(t, 3) * D.P().Z();

	return P;
}

// solve cubic equation
// return case 1: only 1 real root
// return case 2: 3 real root
// return case 3: 1 real root, 2 unreal
int solve_cubic(double a, double b, double c, double d, double &t1, double &t2, double &t3)
{
	double x1, x2, x3;

	//qDebug("\n a*x^3+b*x^2+c*x+d = 0\n");
	/*
	a = 1;
	b = 6;
	c = 12;
	d = 8;
	*/
	//M_LOG3("%fx^3 + %fx^2 + %fx + %f = 0", a, b, c, d);

	double f = (3 * c / a - b * b / (a * a)) / 3;																			//qDebug("f: %f\n", f);
	double g = (2 * b * b * b / (a * a * a) - 9 * b * c / (a * a) + 27 * d / a) / 27; //qDebug("g: %f\n", g);
	double discriminant = g * g / 4 + f * f * f / 27;																	//qDebug("discriminant: %f\n", discriminant);

	if ((f == 0) & (g == 0) & (discriminant == 0))
	{ //When All 3 Roots Are Real and Equal
		M_LOG3("All 3 Roots Are Real and Equal");
		double tt = d / a;
		if (tt < 0)
		{
			tt = -1 * tt;
			x1 = -1 * pow(tt, 1.0 / 3) * (-1);
		}
		else
		{
			x1 = pow(tt, 1.0 / 3) * (-1);
		}
		M_LOG3("x1 = x2 = x3 = %f", x1);
		t1 = x1;
		t2 = x1;
		t3 = x1;
		return 1;
	}
	if (discriminant <= 0)
	{ // ALL 3 Roots Are Real
		M_LOG3("ALL 3 Roots Are Real");
		double i = sqrt((g * g / 4) - discriminant); //qDebug("i: %f\n", i);
		double j = pow(i, 1.0 / 3);									 //qDebug("j: %f\n", j);
		double k = acos(-(g / (2 * i)));						 //qDebug("k: %f\n", k);
		double l = -1 * j;													 //qDebug("l: %f\n", l);
		double m = cos(k / 3);											 //qDebug("m: %f\n", m);
		double n = sqrt(3) * sin(k / 3);						 //qDebug("n: %f\n", n);
		double p = (b / (3 * a)) * (-1);						 //qDebug("p: %f\n", p);
		x1 = 2 * j * cos(k / 3) - (b / (3 * a));
		M_LOG3("x1 = %f", x1);
		x2 = l * (m + n) + p;
		M_LOG3("x2 = %f", x2);
		x3 = l * (m - n) + p;
		M_LOG3("x3 = %f", x3);
		t1 = x1;
		t2 = x2;
		t3 = x3;
		return 2;
	}
	if (discriminant > 0)
	{ //When Only 1 Root Is Real
		M_LOG3("Only 1 Root Is Real");
		double r = -(g / 2) + sqrt(discriminant); //qDebug("r: %f\n", r);
		double s;
		if (r < 0)
		{
			s = -1 * pow(r, 1.0 / 3); //qDebug("s: %f\n", s);
		}
		else
		{
			s = pow(r, 1.0 / 3); //qDebug("s: %f\n", s);
		}
		double t = -(g / 2) - sqrt(discriminant); //qDebug("t: %f\n", t);
		double u;
		if (t < 0)
		{
			t = -1 * t;
			u = -1 * pow(t, 1.0 / 3); //qDebug("u: %f\n", u);
		}
		else
		{
			u = pow(t, 1.0 / 3); //qDebug("u: %f\n", u);
		}
		x1 = s + u - (b / (3 * a));
		M_LOG3("x1 = %f", x1);
		double x2_1 = -1 * (s + u) / 2 - (b / (3 * a));
		double x2_2 = (s - u) * sqrt(3) / 2;
		M_LOG3("x2 = %f + i*%f", x2_1, x2_2);
		M_LOG3("x3 = %f - i*%f", x2_1, x2_2);
		t1 = x1;
		t2 = x2_1;
		t3 = x2_2;
		return 3;
	}
	qDebug("f: %f", f);
	qDebug("g: %f", g);
	qDebug("discriminant: %f", discriminant);
	return 0;
}

// solve quadratic equation
// return case 1: 2 real root
// return case 2: 1 real root
// return case 3: 2 unreal
int solve_quadratic(double a, double b, double c, double &t1, double &t2)
{

	double x1, x2, determinant, realPart, imaginaryPart;
	//cout << "Enter coefficients a, b and c: ";
	//cin >> a >> b >> c;
	determinant = b * b - 4 * a * c;

	if (determinant > 0)
	{
		x1 = (-b + sqrt(determinant)) / (2 * a);
		x2 = (-b - sqrt(determinant)) / (2 * a);
		//cout << "Roots are real and different." << endl;
		//cout << "x1 = " << x1 << endl;
		//cout << "x2 = " << x2 << endl;
		t1 = x1;
		t2 = x2;
		return 1;
	}

	else if (determinant == 0)
	{
		//cout << "Roots are real and same." << endl;
		x1 = (-b + sqrt(determinant)) / (2 * a);
		//cout << "x1 = x2 =" << x1 << endl;
		t1 = x1;
		return 2;
	}

	else
	{
		realPart = -b / (2 * a);
		imaginaryPart = sqrt(-determinant) / (2 * a);
		//cout << "Roots are complex and different." << endl;
		//cout << "x1 = " << realPart << "+" << imaginaryPart << "i" << endl;
		//cout << "x2 = " << realPart << "-" << imaginaryPart << "i" << endl;
		qDebug("ERROR at solve quadratic equation");
	}

	return 0;
}

// check difference of 4 point less than a double or not
bool z_4_points_less_than_x(CVertexO A, CVertexO B, CVertexO C, CVertexO D, CVertexO E, double x)
{
	double array[5];
	array[0] = A.P().Z();
	array[1] = B.P().Z();
	array[2] = C.P().Z();
	array[3] = D.P().Z();
	array[4] = E.P().Z();
	double min = array[0];
	double max = array[0];

	for (int i = 0; i < 5; i++)
	{
		min = (min < array[i]) ? min : array[i];
		max = (max > array[i]) ? max : array[i];
	}

	//qDebug("delta: %f", (max - min));

	if ((max - min) <= x)
	{
		return true;
	}
	else
	{
		return false;
	}
}

//vector in math
struct Vector
{
	float x;
	float y;
	float z;
};

// Calculate a Vector from 2 points
Vector Cal_Vector_AB(CVertexO B, CVertexO A)
{
	Vector V;
	V.x = A.P().X() - B.P().X();
	V.y = A.P().Y() - B.P().Y();
	V.z = A.P().Z() - B.P().Z();
	return V;
}

// Calculate a plus of 2 vectors
Vector Vector_plus_Vector(Vector A, Vector B)
{
	Vector v;
	v.x = A.x + B.x;
	v.y = A.y + B.y;
	v.z = A.z + B.z;
	return v;
}

// Calculate cross product of 2 vectors AxB
Vector Cross_A_B(Vector A, Vector B)
{
	Vector v;
	v.x = A.y * B.z - A.z * B.y;
	v.y = A.z * B.x - A.x * B.z;
	v.z = A.x * B.y - A.y * B.x;
	return v;
}

// Find Normal vector at 8-connectivity (using 8 neighbour of each)
Vector Vec_Normal_8(CVertexO **points_list, int max_y, int k, CVertexO a)
{

	float x = a.P().X();
	float y = a.P().Y();
	bool c1 = false, c2 = false, c3 = false, c4 = false, c5 = false, c6 = false, c7 = false, c8 = false;
	CVertexO a1, a2, a3, a4, a5, a6, a7, a8;
	Vector b1, b2, b3, b4, b5, b6, b7, b8;
	Vector b12, b23, b34, b45, b56, b67, b78, b81;
	Vector sum = {0, 0, 0};
	// check each point of 8
	if (points_list[Cal_index_k(x - k, y - k, max_y, k)] != NULL)
	{
		c1 = true;
		a1 = *points_list[Cal_index_k(x - k, y - k, max_y, k)];
	}
	if (points_list[Cal_index_k(x - k, y, max_y, k)] != NULL)
	{
		c2 = true;
		a2 = *points_list[Cal_index_k(x - k, y, max_y, k)];
	}
	if (points_list[Cal_index_k(x - k, y + k, max_y, k)] != NULL)
	{
		c3 = true;
		a3 = *points_list[Cal_index_k(x - k, y + k, max_y, k)];
	}
	if (points_list[Cal_index_k(x, y + k, max_y, k)] != NULL)
	{
		c4 = true;
		a4 = *points_list[Cal_index_k(x, y + k, max_y, k)];
	}
	if (points_list[Cal_index_k(x + k, y + k, max_y, k)] != NULL)
	{
		c5 = true;
		a5 = *points_list[Cal_index_k(x + k, y + k, max_y, k)];
	}
	if (points_list[Cal_index_k(x + k, y, max_y, k)] != NULL)
	{
		c6 = true;
		a6 = *points_list[Cal_index_k(x + k, y, max_y, k)];
	}
	if (points_list[Cal_index_k(x + k, y - k, max_y, k)] != NULL)
	{
		c7 = true;
		a7 = *points_list[Cal_index_k(x + k, y - k, max_y, k)];
	}
	if (points_list[Cal_index_k(x, y - k, max_y, k)] != NULL)
	{
		c8 = true;
		a8 = *points_list[Cal_index_k(x, y - k, max_y, k)];
	}

	M_LOG4("point a: %f,%f,%f", a.P().X(), a.P().Y(), a.P().Z());

	if (c1 == true && c2 == true)
	{
		b1 = Cal_Vector_AB(a, a1);
		b2 = Cal_Vector_AB(a, a2);
		b12 = Cross_A_B(b1, b2);
		sum = Vector_plus_Vector(sum, b12);
		//M_LOG4("point a1: %f,%f,%f", a1.P().X(),a1.P().Y(),a1.P().Z());
		//M_LOG4("point a2: %f,%f,%f", a2.P().X(),a2.P().Y(),a2.P().Z());
		//M_LOG4("vector b1(a1-a): %f,%f,%f", b1.x,b1.y,b1.z);
		//M_LOG4("vector b2(a2-a): %f,%f,%f", b2.x,b2.y,b2.z);
		//M_LOG4("vector b12(b1xb2): %f,%f,%f", b12.x,b12.y,b12.z);
		//M_LOG4("vector sum: %f,%f,%f", sum.x,sum.y,sum.z);
	}
	if (c2 == true && c3 == true)
	{
		b2 = Cal_Vector_AB(a, a2);
		b3 = Cal_Vector_AB(a, a3);
		b23 = Cross_A_B(b2, b3);
		sum = Vector_plus_Vector(sum, b23);
		//M_LOG4("point a2: %f,%f,%f", a2.P().X(),a2.P().Y(),a2.P().Z());
		//M_LOG4("point a3: %f,%f,%f", a3.P().X(),a3.P().Y(),a3.P().Z());
		//M_LOG4("vector b2(a2-a): %f,%f,%f", b2.x,b2.y,b2.z);
		//M_LOG4("vector b3(a3-a): %f,%f,%f", b3.x,b3.y,b3.z);
		//M_LOG4("vector b23(b2xb3): %f,%f,%f", b23.x,b23.y,b23.z);
		//M_LOG4("vector sum: %f,%f,%f", sum.x,sum.y,sum.z);
	}
	if (c3 == true && c4 == true)
	{
		b3 = Cal_Vector_AB(a, a3);
		b4 = Cal_Vector_AB(a, a4);
		b34 = Cross_A_B(b3, b4);
		sum = Vector_plus_Vector(sum, b34);
		//M_LOG4("point a3: %f,%f,%f", a3.P().X(),a3.P().Y(),a3.P().Z());
		//M_LOG4("point a4: %f,%f,%f", a4.P().X(),a4.P().Y(),a4.P().Z());
		//M_LOG4("vector b3(a3-a): %f,%f,%f", b3.x,b3.y,b3.z);
		//M_LOG4("vector b4(a4-a): %f,%f,%f", b4.x,b4.y,b4.z);
		//M_LOG4("vector b34(b3xb4): %f,%f,%f", b34.x,b34.y,b34.z);
		//M_LOG4("vector sum: %f,%f,%f", sum.x,sum.y,sum.z);
	}
	if (c4 == true && c5 == true)
	{
		b4 = Cal_Vector_AB(a, a4);
		b5 = Cal_Vector_AB(a, a5);
		b45 = Cross_A_B(b4, b5);
		sum = Vector_plus_Vector(sum, b45);
		//M_LOG4("point a4: %f,%f,%f", a4.P().X(),a4.P().Y(),a4.P().Z());
		//M_LOG4("point a5: %f,%f,%f", a5.P().X(),a5.P().Y(),a5.P().Z());
		//M_LOG4("vector b4(a4-a): %f,%f,%f", b4.x,b4.y,b4.z);
		//M_LOG4("vector b5(a5-a): %f,%f,%f", b5.x,b5.y,b5.z);
		//M_LOG4("vector b45(b4xb5): %f,%f,%f", b45.x,b45.y,b45.z);
		//M_LOG4("vector sum: %f,%f,%f", sum.x,sum.y,sum.z);
	}
	if (c5 == true && c6 == true)
	{
		b5 = Cal_Vector_AB(a, a5);
		b6 = Cal_Vector_AB(a, a6);
		b56 = Cross_A_B(b5, b6);
		sum = Vector_plus_Vector(sum, b56);
		//M_LOG4("point a5: %f,%f,%f", a5.P().X(),a5.P().Y(),a5.P().Z());
		//M_LOG4("point a6: %f,%f,%f", a6.P().X(),a6.P().Y(),a6.P().Z());
		//M_LOG4("vector b5(a5-a): %f,%f,%f", b5.x,b5.y,b5.z);
		//M_LOG4("vector b6(a6-a): %f,%f,%f", b6.x,b6.y,b6.z);
		//M_LOG4("vector b56(b5xb6): %f,%f,%f", b56.x,b56.y,b56.z);
		//M_LOG4("vector sum: %f,%f,%f", sum.x,sum.y,sum.z);
	}
	if (c6 == true && c7 == true)
	{
		b6 = Cal_Vector_AB(a, a6);
		b7 = Cal_Vector_AB(a, a7);
		b67 = Cross_A_B(b6, b7);
		sum = Vector_plus_Vector(sum, b67);
		//M_LOG4("point a6: %f,%f,%f", a6.P().X(),a6.P().Y(),a6.P().Z());
		//M_LOG4("point a7: %f,%f,%f", a7.P().X(),a7.P().Y(),a7.P().Z());
		//M_LOG4("vector b6(a6-a): %f,%f,%f", b6.x,b6.y,b6.z);
		//M_LOG4("vector b7(a7-a): %f,%f,%f", b7.x,b7.y,b7.z);
		//M_LOG4("vector b67(b6xb7): %f,%f,%f", b67.x,b67.y,b67.z);
		//M_LOG4("vector sum: %f,%f,%f", sum.x,sum.y,sum.z);
	}
	if (c7 == true && c8 == true)
	{
		b7 = Cal_Vector_AB(a, a7);
		b8 = Cal_Vector_AB(a, a8);
		b78 = Cross_A_B(b7, b8);
		sum = Vector_plus_Vector(sum, b78);
		//M_LOG4("point a7: %f,%f,%f", a7.P().X(),a7.P().Y(),a7.P().Z());
		//M_LOG4("point a8: %f,%f,%f", a8.P().X(),a8.P().Y(),a8.P().Z());
		//M_LOG4("vector b7(a7-a): %f,%f,%f", b7.x,b7.y,b7.z);
		//M_LOG4("vector b8(a8-a): %f,%f,%f", b8.x,b8.y,b8.z);
		//M_LOG4("vector b78(b7xb8): %f,%f,%f", b78.x,b78.y,b78.z);
		//M_LOG4("vector sum: %f,%f,%f", sum.x,sum.y,sum.z);
	}
	if (c1 == true && c8 == true)
	{
		b1 = Cal_Vector_AB(a, a8);
		b8 = Cal_Vector_AB(a, a1);
		b81 = Cross_A_B(b8, b1);
		sum = Vector_plus_Vector(sum, b81);
		//M_LOG4("point a8: %f,%f,%f", a8.P().X(),a8.P().Y(),a8.P().Z());
		//M_LOG4("point a1: %f,%f,%f", a1.P().X(),a1.P().Y(),a1.P().Z());
		//M_LOG4("vector b8(a8-a): %f,%f,%f", b8.x,b8.y,b8.z);
		//M_LOG4("vector b1(a1-a): %f,%f,%f", b1.x,b1.y,b1.z);
		//M_LOG4("vector b81(b8xb1): %f,%f,%f", b81.x,b81.y,b81.z);
		//M_LOG4("vector sum: %f,%f,%f", sum.x,sum.y,sum.z);
	}

	//Calculate unit vector
	//float sqr = (float) sqrt(sum.x*sum.x + sum.y*sum.y + sum.z*sum.z);
	//sum.x = sum.x/sqr;
	//sum.y = sum.y/sqr;
	//sum.z = sum.z/sqr;
	//M_LOG4("sqr: %f", sqr);
	M_LOG4("Unit Vector sum: %f,%f,%f \n", sum.x, sum.y, sum.z);

	return sum;
}

// Find z of CVertexO a from 1 dot product
float find_z(Vector normal_8, CVertexO b, CVertexO a)
{
	float z_of_a = -1 * (normal_8.x * (a.P().X() - b.P().X()) + normal_8.y * (a.P().Y() - b.P().Y())) / normal_8.z + b.P().Z();
	return z_of_a;
}

// Refine 8 z of an inserted point (from dot product with normal vectors of 8 connectivity)
float z_refine(CVertexO **points_list, int max_y, int k, CVertexO a)
{
	int x = (int)a.P().X();
	int y = (int)a.P().Y();
	float z_ave = 0;
	CVertexO a1, a2, a3, a4, a5, a6, a7, a8;
	bool c1 = false, c2 = false, c3 = false, c4 = false, c5 = false, c6 = false, c7 = false, c8 = false;
	Vector b1, b2, b3, b4, b5, b6, b7, b8;
	int num_z = 0;

	if (points_list[Cal_index_k(x - k, y - k, max_y, k)] != NULL)
	{
		c1 = true;
		a1 = *points_list[Cal_index_k(x - k, y - k, max_y, k)];
		b1 = Vec_Normal_8(points_list, max_y, k, a1);
		z_ave += find_z(b1, a1, a);
		num_z++;
		M_LOG4("z_ave 1: %f \n", find_z(b1, a1, a));
	}
	if (points_list[Cal_index_k(x - k, y, max_y, k)] != NULL)
	{
		c2 = true;
		a2 = *points_list[Cal_index_k(x - k, y, max_y, k)];
		b2 = Vec_Normal_8(points_list, max_y, k, a2);
		z_ave += find_z(b2, a2, a);
		num_z++;
		M_LOG4("z_ave 2: %f \n", find_z(b2, a2, a));
	}
	if (points_list[Cal_index_k(x - k, y + k, max_y, k)] != NULL)
	{
		c3 = true;
		a3 = *points_list[Cal_index_k(x - k, y + k, max_y, k)];
		b3 = Vec_Normal_8(points_list, max_y, k, a3);
		z_ave += find_z(b3, a3, a);
		num_z++;
		M_LOG4("z_ave 3: %f \n", find_z(b3, a3, a));
	}
	if (points_list[Cal_index_k(x, y + k, max_y, k)] != NULL)
	{
		c4 = true;
		a4 = *points_list[Cal_index_k(x, y + k, max_y, k)];
		b4 = Vec_Normal_8(points_list, max_y, k, a4);
		z_ave += find_z(b4, a4, a);
		num_z++;
		M_LOG4("z_ave 4: %f \n", find_z(b4, a4, a));
	}
	if (points_list[Cal_index_k(x + k, y + k, max_y, k)] != NULL)
	{
		c5 = true;
		a5 = *points_list[Cal_index_k(x + k, y + k, max_y, k)];
		b5 = Vec_Normal_8(points_list, max_y, k, a5);
		z_ave += find_z(b5, a5, a);
		num_z++;
		M_LOG4("z_ave 5: %f \n", find_z(b5, a5, a));
	}
	if (points_list[Cal_index_k(x + k, y, max_y, k)] != NULL)
	{
		c6 = true;
		a6 = *points_list[Cal_index_k(x + k, y, max_y, k)];
		b6 = Vec_Normal_8(points_list, max_y, k, a6);
		z_ave += find_z(b6, a6, a);
		num_z++;
		M_LOG4("z_ave 6: %f \n", find_z(b6, a6, a));
	}
	if (points_list[Cal_index_k(x + k, y - k, max_y, k)] != NULL)
	{
		c7 = true;
		a7 = *points_list[Cal_index_k(x + k, y - k, max_y, k)];
		b7 = Vec_Normal_8(points_list, max_y, k, a7);
		z_ave += find_z(b7, a7, a);
		num_z++;
		M_LOG4("z_ave 7: %f \n", find_z(b7, a7, a));
	}
	if (points_list[Cal_index_k(x, y - k, max_y, k)] != NULL)
	{
		c8 = true;
		a8 = *points_list[Cal_index_k(x, y - k, max_y, k)];
		b8 = Vec_Normal_8(points_list, max_y, k, a8);
		z_ave += find_z(b8, a8, a);
		num_z++;
		M_LOG4("z_ave 8: %f \n", find_z(b8, a8, a));
	}
	z_ave = z_ave / num_z;
	M_LOG4("z_ave: %f \n", z_ave);
	return z_ave;
}

// Move to next position on 8 connect
bool move_next_p_8_con(int &pre_x, int &pre_y, int curr_x, int curr_y)
{
	if ((curr_x - pre_x) == 1 && (curr_y - pre_y) == 1)
	{
		pre_y++;
		return true;
	}
	if ((curr_x - pre_x) == 1 && (curr_y - pre_y) == 0)
	{
		pre_y++;
		return true;
	}
	if ((curr_x - pre_x) == 1 && (curr_y - pre_y) == -1)
	{
		pre_x++;
		return true;
	}
	if ((curr_x - pre_x) == 0 && (curr_y - pre_y) == -1)
	{
		pre_x++;
		return true;
	}
	if ((curr_x - pre_x) == -1 && (curr_y - pre_y) == -1)
	{
		pre_y--;
		return true;
	}
	if ((curr_x - pre_x) == -1 && (curr_y - pre_y) == 0)
	{
		pre_y--;
		return true;
	}
	if ((curr_x - pre_x) == -1 && (curr_y - pre_y) == 1)
	{
		pre_x--;
		return true;
	}
	if ((curr_x - pre_x) == 0 && (curr_y - pre_y) == 1)
	{
		pre_x--;
		return true;
	}
	return false;
}

// /////////////////////////////////////////////
// Main
// /////////////////////////////////////////////
std::map<std::string, QVariant> TestBezierPlugin::applyFilter(
		const QAction *action,
		const RichParameterList &par,
		MeshDocument &md,
		unsigned int & /*postConditionMask*/,
		vcg::CallBackPos * /*cb*/)
{
	std::map<std::string, QVariant> outputValues;

	//int k = par.getInt("k");
	bool case1 = par.getBool("case1");
	bool case2 = par.getBool("case2");
	bool case3 = par.getBool("case3");
	bool case4 = par.getBool("case4");

	if (md.mm() == NULL)
		return outputValues;

	MeshModel &m = *(md.mm());
	//m.cm.vert.RadiusEnabled = true;

	/*
	for (int i = 0; i < 10; i++)
	{
	//m.cm.vert[i].C().SetRGB(255, 0, 0);

	CVertexO vertex;
	vertex.P() = m.cm.vert[i].P() * 2;
	vertex.C().SetRGB(0, 255 ,0) ;
	m.cm.vert.push_back(vertex);
	}
        */

	//convert 1-dimention to 2-dimention:
	int vertex_numbers = m.cm.vn;
	M_LOG("\n\n Number of point in the file: %d", vertex_numbers); //count the number of point in the file

	//------clone all mesh to tracking_mesh and find max index of x & y

	int max_index_x = 0, max_index_y = 0;
	int min_index_x = (int)m.cm.vert[0].P().X();
	int min_index_y = (int)m.cm.vert[0].P().Y();
	M_LOG("First value at index (0,0) x:%d --- y:%d ", min_index_x, min_index_y);

	int sta_index_x = 0, sta_index_y = 0;

	for (int i = 0; i < vertex_numbers; i++)
	{
		//find max
		if ((int)m.cm.vert[i].P().X() > max_index_x)
			max_index_x = (int)m.cm.vert[i].P().X();
		if ((int)m.cm.vert[i].P().Y() > max_index_y)
			max_index_y = (int)m.cm.vert[i].P().Y();
		//find min
		if ((int)m.cm.vert[i].P().X() < min_index_x)
			min_index_x = (int)m.cm.vert[i].P().X();
		if ((int)m.cm.vert[i].P().Y() < min_index_y)
			min_index_y = (int)m.cm.vert[i].P().Y();
	}

	//find start
	sta_index_x = min_index_x;
	sta_index_y = max_index_y;
	for (int i = 0; i < vertex_numbers; i++)
	{
		if ((int)m.cm.vert[i].P().X() == min_index_x)
		{
			if ((int)m.cm.vert[i].P().Y() <= sta_index_y)
			{
				sta_index_y = (int)m.cm.vert[i].P().Y();
			}
		}
	}

	M_LOG("Min value at x:%d --- y:%d ", min_index_x, min_index_y);
	M_LOG("Max index at x:%d --- y:%d ", max_index_x, max_index_y);

	//----create a clone for tracking expanded with k----
	int expanded_x = max_index_x + 2 * k_max;
	int expanded_y = max_index_y + 2 * k_max;
	int max_size = expanded_x * expanded_y;
	int *track_mesh = new int[max_size];
	memset(track_mesh, NULL, max_size * sizeof(int));

	//---start to clone---
	for (int i = 0; i < vertex_numbers; i++)
	{
		int x = (int)m.cm.vert[i].P().X();
		int y = (int)m.cm.vert[i].P().Y();
		track_mesh[Cal_index_k(x, y, expanded_y, k_max)] = norm_p;
	}

	//---clone to another array to compute---
	CVertexO **xy = new CVertexO *[max_size];
	memset(xy, NULL, max_size * sizeof(CVertexO *));
	for (int i = 0; i < vertex_numbers; i++)
	{
		int a = (int)m.cm.vert[i].P().X();
		int b = (int)m.cm.vert[i].P().Y();
		xy[Cal_index_k(a, b, expanded_y, k_max)] = &m.cm.vert[i];
	}

	//---create list of points was tracked---
	CVertexO *tracked_points = new CVertexO[vertex_numbers];
	int tracked_index = 0;
	int total_tracked = 0;

	//-----------------------START EXTRACT EXTERNAL BOUNDARY-----------------------
	if (case1)
	{
		bool stop = false;
		bool find;

		int k = 1;

		// Setup the first point
		CVertexO p_First, p_Curr, p_Prev, p_Prev_k;
		int x_First = sta_index_x, y_First = sta_index_y;
		int x_Curr = sta_index_x, y_Curr = sta_index_y;

		M_LOG("Value at p_First x(%d):%d --- y(%d):%d ",
					x_Curr, (int)xy[Cal_index_k(x_Curr, y_Curr, expanded_y, k_max)]->P().X(),
					y_Curr, (int)xy[Cal_index_k(x_Curr, y_Curr, expanded_y, k_max)]->P().Y());

		//p_First.ImportData(*xy[start_index]);
		//track_mesh[start_index] = ex_bound;
		//xy[start_index]->C().SetRGB(255, 0, 0);
		//Log("Value at x:%d --- y:%d ", (int)p_First.P().X(), (int)p_First.P().Y());

		// Setup the previous point
		int x_Prev, y_Prev;
		//find_first_of_8_neigh(xy, min_index_x, min_index_y, x_Prev, y_Prev, expanded_y, k);
		//**************debug begin from previous point*****************
		x_Prev = x_Curr;
		y_Prev = y_Curr - 1;
		/////////////////////////////////////////////////////////////

		//**************debug begin from this point*****************
		/*
		x_Curr = 854;
		y_Curr = 205;
		x_Prev = 852;
		y_Prev = 207;
		track_mesh[Cal_index_k(x_Curr, y_Curr, expanded_y, k_max)] = ex_bound;
		xy[Cal_index_k(x_Curr, y_Curr, expanded_y, k_max)]->C().SetRGB(255, 0, 0);
		track_mesh[Cal_index_k(x_Prev, y_Prev, expanded_y, k_max)] = ex_bound;
		xy[Cal_index_k(x_Prev, y_Prev, expanded_y, k_max)]->C().SetRGB(255, 0, 0);
		*/
		/////////////////////////////////////////////////////////////

		float x_Prev_k, y_Prev_k;
		int x_pre_k, y_pre_k;
		//Comp_p_Prev_p_Curr(x_Prev, y_Prev, x_Curr, y_Curr, expanded_y, k, x_Prev_k, y_Prev_k);
		int stop_loop = 0;
		int stop_loop2 = 0;

		int x_ex_1st_track, y_ex_1st_track;
		bool ex_1st_track = false;

		int x_pre_k_1st, y_pre_k_1st;
		bool k_1st_track = false; // used for cycle, avoid loop forever

		bool pass_1st_prev = true;

		while (!stop)
		{

			// Setup the p_Prev_k (which created by (p_curr,p_prev) and RNk(p)) and create new after loop
			M_LOG("\n\n\t\t\t *** NEW PROCESS ***");
			M_LOG("*** p_prev (%d,%d) ***", x_Prev, y_Prev);
			M_LOG("*** p_curr (%d,%d) ***", x_Curr, y_Curr);
			int line_k;
			M_LOG("\t ** Begin to find to intersection (p_Prev_k) between line (p_curr,p_prev) and square (RNk) **");
			Comp_p_Prev_p_Curr(x_Prev, y_Prev, x_Curr, y_Curr, expanded_y, k, x_Prev_k, y_Prev_k, line_k);
			M_LOG("\t ***** x_Prev_k, y_Prev_k (%f,%f) *****", x_Prev_k, y_Prev_k);
			find = false;
			round_by_line(x_Prev_k, y_Prev_k, line_k, x_pre_k, y_pre_k);

			// tracking first prev_k
			if (k_1st_track)
			{
				x_pre_k_1st = x_pre_k;
				y_pre_k_1st = y_pre_k;
				k_1st_track = true;
			}

			//while (q = p_Prev_k; (any q belong RNk(p) && !find);) {
			while (!find)
			{

				move_next_RNk(x_pre_k, y_pre_k, x_Curr, y_Curr, k);
				M_LOG("\n\t\t\t * Move to next point RNk(p) x_Prev_k, y_Prev_k (%d,%d) *", x_pre_k, y_pre_k);

				// Case 1: check point from i to k
				for (int i = k; i <= k; i++)
				{
					//q_i = intersection between RNi(p) and discrete segment(p, q).
					float x_qi, y_qi;

					M_LOG("\n\t Find intersection between RNi(p)-[i=%d] and discrete segment p_prev_k,p_curr(%d,%d),(%d,%d).",
								i, x_pre_k, y_pre_k, x_Curr, y_Curr);

					int line_i;
					Comp_p_Prev_p_Curr(x_pre_k, y_pre_k, x_Curr, y_Curr, expanded_y, i, x_qi, y_qi, line_i);
					int x_i, y_i;
					round_by_line(x_qi, y_qi, line_i, x_i, y_i);

					M_LOG("\t\t\t ** (i=%d)The intersection between RNi(p) and discrete segment(p, q) is: (%f,%f) **",
								i, x_qi, y_qi);
					M_LOG("\t\t\t ** (i=%d)The intersection between RNi(p) and discrete segment(p, q) is: (%d,%d) **",
								i, x_i, y_i);

					M_LOG("Min index (%d,%d); Max index (%d,%d)", min_index_x, min_index_y, max_index_x, max_index_y);

					// check if out of range => skip
					if ((min_index_x <= x_i) && (x_i <= max_index_x) && (min_index_y <= y_i) && (y_i <= max_index_y))
					{
						// skip some interator when meet tracked points
						//if (k > 1 &&
						//pass_1st_prev &&
						//(x_i == x_Prev) &&
						//(y_i == y_Prev))
						if ((k > 1) &&
								pass_1st_prev)
						{
							if ((track_mesh[Cal_index_k(x_i, y_i, expanded_y, k_max)] == ex_bound) ||
									((x_i == x_Prev) && (y_i == y_Prev)))
							{
								M_LOG("REJECT (%d,%d) (k > 1) when check this point be tracked", x_i, y_i);
								break;
							}
							else
								pass_1st_prev = false;
						}

						// IF CASE of case 4: tracking first external boundary point (belong i)
						if (!ex_1st_track &&
								(track_mesh[Cal_index_k(x_i, y_i, expanded_y, k_max)] == ex_bound))
						{
							x_ex_1st_track = x_i;
							y_ex_1st_track = y_i;
							ex_1st_track = true;
							M_LOG("\t SET (%d,%d) is the first tracked point meet", x_ex_1st_track, y_ex_1st_track);
						}

						// no point -> exit this interator
						if (xy[Cal_index_k(x_i, y_i, expanded_y, k_max)] == NULL)
						{
							M_LOG("At (%d,%d), there is no point", x_i, y_i);
							continue;
						}

						// pass all requirement => choose this point
						if ((track_mesh[Cal_index_k(x_i, y_i, expanded_y, k_max)] != ex_bound) && (xy[Cal_index_k(x_i, y_i, expanded_y, k_max)] != NULL))
						{
							if (miss_1_of_4_neigh(xy, x_i, y_i, expanded_y,
																		max_index_x, max_index_y,
																		min_index_x, min_index_y, k_max))
							{
								M_LOG("\t\t\t *** Found the next point (%d,%d) ***", x_i, y_i);
								//p_Prev = p;
								x_Prev = x_Curr;
								y_Prev = y_Curr;
								M_LOG("\t\t\t *** New point p_prev (%d,%d) ***", x_Prev, y_Prev);
								//p_Next = q_i;
								x_Curr = x_i;
								y_Curr = y_i;
								M_LOG("\t\t\t *** New point p_curr (%d,%d) ***", x_Curr, y_Curr);

								// set info
								set_point(tracked_points, total_tracked, x_i, y_i);
								ex_1st_track = false;
								tracked_index = total_tracked;
								track_mesh[Cal_index_k(x_i, y_i, expanded_y, k_max)] = ex_bound;
								xy[Cal_index_k(x_i, y_i, expanded_y, k_max)]->C().SetHSVColor(0, 1.0f, 1.0f); // red color

								find = true;
								stop_loop++;
								stop_loop2 = 0;
								k_1st_track = true;
								pass_1st_prev = true;
								k = 1;

								break;
							}
						}
					}
					else
					{
						M_LOG("At (%d,%d), out of range", x_i, y_i);
						pass_1st_prev = false; // set this when out of range
						break;
					}
				}

				// Case 2: If not found any point -> set first tracked point is next point
				if ((x_pre_k == x_pre_k_1st) &&
						(y_pre_k == y_pre_k_1st) &&
						ex_1st_track &&
						k == k_max)
				{
					M_LOG("\t\t\t *** Set the tracked point to new point (%d,%d) ***", x_ex_1st_track, y_ex_1st_track);
					//p_Prev = p;
					x_Prev = x_Curr;
					y_Prev = y_Curr;
					M_LOG("\t\t\t *** New point p_prev (%d,%d) ***", x_Prev, y_Prev);
					//p_Next = q_i;
					x_Curr = x_ex_1st_track;
					y_Curr = y_ex_1st_track;
					M_LOG("\t\t\t *** New point p_curr (%d,%d) ***", x_Curr, y_Curr);

					// set info
					set_point(tracked_points, total_tracked, x_ex_1st_track, y_ex_1st_track);
					ex_1st_track = false;
					tracked_index = total_tracked;

					find = true;
					stop_loop++;
					stop_loop2 = 0;
					k_1st_track = true;
					pass_1st_prev = true;
					k = 1;
				}
				else
				{
					M_LOG("\t\t\t *** ex_1st_track (%d,%d) ***", x_ex_1st_track, y_ex_1st_track);
					M_LOG("\t\t\t *** pre_k_1st (%d,%d) ***", x_pre_k_1st, y_pre_k_1st);
					M_LOG("\t\t\t *** pre_k (%d,%d) ***", x_pre_k, y_pre_k);
					M_LOG("\t\t\t *** ex_1st_track (%d) ***", ex_1st_track);
				}

				// increasing k when not found any point
				if ((x_pre_k == x_pre_k_1st) &&
						(y_pre_k == y_pre_k_1st) &&
						k < k_max)
				{
					k++;
					M_LOG("\t\t\t *** INCREASING K (%d) ***", k);
					stop_loop2 = 0;
					M_LOG("\t\t\t *** RESET STOP LOOP (%d) ***", stop_loop2);
					break;
				}

				//qDebug("total_tracked: %d", total_tracked);

				int max_loop_2 = 4 * (2 * k + 1) + 2;
				stop_loop2++;
				M_LOG("(%d) STOP LOOP HERE, (max: %d)", stop_loop2, max_loop_2);
				if (stop_loop2 == max_loop_2)
				{
					find = true;
					stop = true;
				}

				//if (stop_loop == 500) { find = true; stop = true; }
			}
			// stop of while(!find)

			//if (p.P() == p_First.P()) {
			if (x_Curr == x_First && y_Curr == y_First)
			{
				M_LOG("\t\t\t *** END ***");
				stop = true;
			}
		} // stop of while(!stop)

		//-----------------------START EXTRACT INTERNAL BOUNDARY-----------------------
		// CHECKBOX 1&2
		if (case2)
		{
			vector<vector<CVertexO>> hole_list;
			vector<CVertexO> points_list;
			int total_inter_points = 0;
			int min_x = 0, min_y = 0;

			// find the inter points
			for (int i = 0; i < vertex_numbers; i++)
			{
				int x = (int)m.cm.vert[i].P().X();
				int y = (int)m.cm.vert[i].P().Y();

				if ((track_mesh[Cal_index_k(x, y, expanded_y, k_max)] != ex_bound) && miss_1_of_4_neigh(xy, x, y, expanded_y,
																																																max_index_x, max_index_y,
																																																min_index_x, min_index_y, k_max))
				{
					track_mesh[Cal_index_k(x, y, expanded_y, k_max)] = in_bound;
					// debug
                    xy[Cal_index_k(x, y, expanded_y, k_max)]->C().SetHSVColor(0.6f, 1.0f, 0.392f); // Yellow & Brown color
					CVertexO temp = *xy[Cal_index_k(x, y, expanded_y, k_max)];
					points_list.push_back(temp);
					total_inter_points++;
				}
			}

			for (unsigned int i = 0; i < points_list.size(); i++)
			{
				M_LOG2("ALL DATA: i(%d): x(%f) y(%f)", i, points_list[i].P().X(), points_list[i].P().Y());
			}

			// kk is number of holes
			int kk = 0;

			// start to find inter points of holes
			while (total_inter_points > 0)
			{

				first_inter_point(min_x, min_y, total_inter_points, points_list);
				M_LOG2("min_x (%d), min_y (%d)", min_x, min_y);

				int x_1st_inter = min_x;
				int y_1st_inter = min_y;

				vector<CVertexO> hole;

				int x_curr_inter = min_x + 0;
				int y_curr_inter = min_y + 0;
				M_LOG2("Start point x_curr_inter (%d), y_curr_inter (%d)", x_curr_inter, y_curr_inter);

				x_1st_inter = x_curr_inter;
				y_1st_inter = y_curr_inter;

				int x_pre_inter = x_curr_inter + 0;
				int y_pre_inter = y_curr_inter + 1;

				//while (xy[Cal_index_k(x_pre_inter, y_pre_inter, expanded_y, k_max)] == NULL)
				//{
				//move_next_RNk(x_pre_inter, y_pre_inter, x_curr_inter, y_curr_inter, 1);
				//}
				int temp = 0;
				int x_pre_inter_1st = 0;
				int y_pre_inter_1st = 0;
				while (true)
				{
					move_next_RNk(x_pre_inter, y_pre_inter, x_curr_inter, y_curr_inter, 1);

					if (xy[Cal_index_k(x_pre_inter, y_pre_inter, expanded_y, k_max)] == NULL)
					{
						// set first point in hole
						x_pre_inter_1st = x_pre_inter;
						y_pre_inter_1st = y_pre_inter;
						int x_temp = x_pre_inter;
						int y_temp = y_pre_inter;
						move_next_RNk(x_temp, y_temp, x_curr_inter, y_curr_inter, 1);
						if (track_mesh[Cal_index_k(x_temp, y_temp, expanded_y, k_max)] == in_bound)
							break;
					}
					temp++;

					if (temp == 8)
					{
						x_pre_inter = x_pre_inter_1st;
						y_pre_inter = y_pre_inter_1st;
						break;
					}
				}
				M_LOG2("Start point x_pre_inter (%d), y_pre_inter (%d)", x_pre_inter, y_pre_inter);

				bool stop_inter = false;
				bool found_1_point = false;
				M_LOG2("------------index of hole: %d------------", kk);

				while (!stop_inter)
				{
					// move pre point to next point on RNk ring
					move_next_RNk(x_pre_inter, y_pre_inter, x_curr_inter, y_curr_inter, 1);
					M_LOG2("* Move to next point RNk(p) x_pre_inter (%d), y_pre_inter (%d) *", x_pre_inter, y_pre_inter);

					if (track_mesh[Cal_index_k(x_pre_inter, y_pre_inter, expanded_y, k_max)] >= norm_p)
					{
						if (track_mesh[Cal_index_k(x_pre_inter, y_pre_inter, expanded_y, k_max)] == ex_bound)
						{
							// if meet external point -> check whether it have 1 neighbourd is internal or no normal in neighbourds
							if (check_interal_8_neigh(track_mesh, x_pre_inter, y_pre_inter, expanded_y, k_max) ||
									no_norm_8_neigh(track_mesh, x_pre_inter, y_pre_inter, expanded_y, k_max))
							{
								//track_mesh[Cal_index_k(x_pre_inter, y_pre_inter, expanded_y, k_max)] = in_bound; // set ex to in
								xy[Cal_index_k(x_pre_inter, y_pre_inter, expanded_y, k_max)]->C().SetHSVColor(0.3f, 1.0f, 1.0f); // orange color
								M_LOG2("\t\t\t * Found External point: x_pre_inter (%d), y_pre_inter (%d) *", x_pre_inter, y_pre_inter);
								int x_temp = x_pre_inter;
								int y_temp = y_pre_inter;

								x_pre_inter = x_curr_inter;
								y_pre_inter = y_curr_inter;

								x_curr_inter = x_temp;
								y_curr_inter = y_temp;

								found_1_point = true;
								CVertexO temp = *xy[Cal_index_k(x_temp, y_temp, expanded_y, k_max)];
								//erase(points_list, temp);
								hole.push_back(temp);
								//total_inter_points--;
							}
						}
						else if (track_mesh[Cal_index_k(x_pre_inter, y_pre_inter, expanded_y, k_max)] == in_bound)
						{
							xy[Cal_index_k(x_pre_inter, y_pre_inter, expanded_y, k_max)]->C().SetHSVColor(3.0f, 1.0f, 1.0f); // violet color
							M_LOG2("\t\t\t * Found Interior point: x_pre_inter (%d), y_pre_inter (%d) *", x_pre_inter, y_pre_inter);
							int x_temp = x_pre_inter;
							int y_temp = y_pre_inter;

							x_pre_inter = x_curr_inter;
							y_pre_inter = y_curr_inter;

							x_curr_inter = x_temp;
							y_curr_inter = y_temp;

							found_1_point = true;
							CVertexO temp = *xy[Cal_index_k(x_temp, y_temp, expanded_y, k_max)];
							bool found = erase(points_list, temp);
							if (found)
								total_inter_points--;
							hole.push_back(temp);
						}
					}
					else
					{
						M_LOG2("* No point: x_pre_inter (%d), y_pre_inter (%d) *", x_pre_inter, y_pre_inter);
					}

					if (found_1_point && (x_curr_inter == x_1st_inter) && (y_curr_inter == y_1st_inter))
					{
						stop_inter = true;
						hole_list.push_back(hole); // copy a hole to list holes
						kk++;
					}

					// debug when travelling number of point, set by kk
					/*
					if (kk == 300) {
					stop_inter = true;
					//found_1_point = true;
					total_inter_points = 0;
					}
					else kk++;
					*/
				} // stop while (!stop_inter)

				// debug after a number of holes
				/*
				if (kk == 3) {
				stop_inter = true;
				found_1_point = true;
				total_inter_points = 0;
				}*/
				//else kk++;

				//qDebug("number of holes: %d", kk);

				M_LOG2("\t\t\t------------DONE ONE HOLE (m: %d)------------", kk);
			} // stop total_inter_points

			for (unsigned int i = 0; i < hole_list.size(); i++)
			{
				//M_LOG3("HOLE i(%d)", i);
				for (unsigned int j = 0; j < hole_list.at(i).size(); j++)
				{
					//M_LOG3("* X (%f), Y (%f) *", hole_list.at(i).at(j).P().X(), hole_list.at(i).at(j).P().Y());
				}
			}

			//-----------------------START FILLING HOLES-----------------------
			// CHECKBOX 1&2&3
			if (case3)
			{
				M_LOG3("--------------------------------------------------------------------------------------------------------------");
				M_LOG3("\t\t\t START FILLING THE HOLE");

				vector<CVertexO> new_points_list_row;
				vector<CVertexO> new_points_list_col;

				// sorting list of point to matrix
				for (unsigned int i = 0; i < hole_list.size(); i++)
				//for (unsigned int i = 9; i < 10; i++)
				{
					M_LOG3("HOLE i(%d)", i);
					vector<vector<CVertexO>> one_hole_row; // for row
					vector<vector<CVertexO>> one_hole_col; // for col
					vector<CVertexO> one_row;
					vector<CVertexO> one_col;

					int x = hole_list.at(i).at(0).P().X();
					int y = hole_list.at(i).at(0).P().Y();
					int points = hole_list.at(i).size();

					// find min of x,y of a hole
					for (int j = 0; j < points; j++)
					{
						if (hole_list.at(i).at(j).P().Y() <= y)
						{
							y = hole_list.at(i).at(j).P().Y();
						}
					}
					for (int j = 0; j < points; j++)
					{
						if (hole_list.at(i).at(j).P().X() <= x)
						{
							x = hole_list.at(i).at(j).P().X();
						}
					}
					int x_min = x;
					int y_min = y;
					M_LOG3("* x_min (%d), y_min (%d) *", x_min, y_min);

					x = hole_list.at(i).at(0).P().X();
					y = hole_list.at(i).at(0).P().Y();

					// find max of x,y of a hole
					for (int j = 0; j < points; j++)
					{
						if (hole_list.at(i).at(j).P().Y() >= y)
						{
							y = hole_list.at(i).at(j).P().Y();
						}
					}
					for (int j = 0; j < points; j++)
					{
						if (hole_list.at(i).at(j).P().X() >= x)
						{
							x = hole_list.at(i).at(j).P().X();
						}
					}
					int x_max = x;
					int y_max = y;
					M_LOG3("* x_max (%d), y_max (%d) *", x_max, y_max);

					// copy a row (the same Y) to one_row
					// find left+right of begin and end row
					for (int a = y_min; a <= y_max; a++)
					{
						for (int j = 0; j < points; j++)
						{
							if (a == hole_list.at(i).at(j).P().Y())
							{
								//one_row.push_back(*xy[Cal_index_k(hole_list.at(i).at(j).P().X(), a, expanded_y, k_max)]);
								one_row.push_back(hole_list.at(i).at(j));
							}
						}
						BubbleSort(one_row, 1);
						M_LOG3("------------");
						M_LOG3("min X: %f, Y: %f", one_row.begin()->P().X(), one_row.begin()->P().Y());
						M_LOG3("max X: %f, Y: %f", one_row.at(one_row.size() - 1).P().X(), one_row.at(one_row.size() - 1).P().Y());
						CVertexO begin, end;
						vector<CVertexO> one_row_new;
						for (unsigned int b = 0; b < one_row.size(); b++)
						{
							M_LOG3("X(%f), Y(%f)", one_row.at(b).P().X(), one_row.at(b).P().Y());
						}
						if (one_row.size() > 1)
						{
							for (unsigned int b = 0; b < one_row.size(); b++)
							{
								M_LOG3("curr X: %f", one_row.at(b).P().X());
								if (b > 0)
								{
									int begin_x = one_row.at(b - 1).P().X();
									int end_x = one_row.at(b).P().X();
									int end_left = one_row.at(b).P().X() - 1;
									if (((end_x - begin_x) > 1) &&
											(xy[Cal_index_k(end_left, one_row.at(b).P().Y(), expanded_y, k_max)] == NULL))
									{
										begin = one_row.at(b - 1);
										end = one_row.at(b);
										one_row_new.push_back(begin);
										one_row_new.push_back(end);
										one_hole_row.push_back(one_row_new);
										one_row_new.clear();
									}
								}
							}
							//one_hole_row.push_back(one_row);
							//one_row.clear();
						}
						one_row.clear();
					}
					M_LOG3("------------Created 1 Hole by ROW------------");

					////////////////////////////////////////////////////////
					// copy a col (the same X) to one_col
					// find left+right of begin and end col
					///*
					for (int a = x_min; a <= x_max; a++)
					{
						for (int j = 0; j < points; j++)
						{
							if (a == hole_list.at(i).at(j).P().X())
							{
								//one_col.push_back(*xy[Cal_index_k(hole_list.at(i).at(j).P().X(), a, expanded_y, k_max)]);
								one_col.push_back(hole_list.at(i).at(j));
							}
						}
						BubbleSort(one_col, 0);
						M_LOG3("------------");
						M_LOG3("min X: %f, Y: %f", one_col.begin()->P().X(), one_col.begin()->P().Y());
						M_LOG3("max X: %f, Y: %f", one_col.at(one_col.size() - 1).P().X(), one_col.at(one_col.size() - 1).P().Y());
						CVertexO begin, end;
						vector<CVertexO> one_col_new;
						for (unsigned int b = 0; b < one_col.size(); b++)
						{
							M_LOG3("X(%f), Y(%f)", one_col.at(b).P().X(), one_col.at(b).P().Y());
						}
						if (one_col.size() > 1)
						{
							for (unsigned int b = 0; b < one_col.size(); b++)
							{
								M_LOG3("curr Y: %f", one_col.at(b).P().Y());
								if (b > 0)
								{
									int begin_y = one_col.at(b - 1).P().Y();
									int end_y = one_col.at(b).P().Y();
									int end_left = one_col.at(b).P().Y() - 1;
									if (((end_y - begin_y) > 1) &&
											(xy[Cal_index_k(one_col.at(b).P().X(), end_left, expanded_y, k_max)] == NULL))
									{
										begin = one_col.at(b - 1);
										end = one_col.at(b);
										one_col_new.push_back(begin);
										one_col_new.push_back(end);
										one_hole_col.push_back(one_col_new);
										one_col_new.clear();
									}
								}
							}
							//one_hole_row.push_back(one_row);
							//one_row.clear();
						}
						one_col.clear();
					}
					M_LOG3("------------Created 1 Hole by COL------------");
					//*/

					for (unsigned int i1 = 0; i1 < one_hole_row.size(); i1++)
					{
						M_LOG3("------ Number of rows: i1(%d)", i1);
						for (unsigned int j = 0; j < one_hole_row.at(i1).size(); j++)
						{
							M_LOG3("* X (%f), Y (%f) *", one_hole_row.at(i1).at(j).P().X(), one_hole_row.at(i1).at(j).P().Y());
						}
					}
					////////////////////////////////////////////////////////
					for (unsigned int i1 = 0; i1 < one_hole_col.size(); i1++)
					{
						M_LOG3("------ Number of col: i1(%d)", i1);
						for (unsigned int j = 0; j < one_hole_col.at(i1).size(); j++)
						{
							M_LOG3("* X (%f), Y (%f) *", one_hole_col.at(i1).at(j).P().X(), one_hole_col.at(i1).at(j).P().Y());
						}
					}

					// create new points on row
					for (unsigned int i2 = 0; i2 < one_hole_row.size(); i2++)
					{
						M_LOG3("\t- At row(i): (%d) -", i2);
						int check_size = one_hole_row.at(i2).begin()->P().X() - one_hole_row.at(i2).at(one_hole_row.at(i2).size() - 1).P().X() + 1;
						if (one_hole_row.at(i2).size() == 1)
							M_LOG3("Row have only one point -> SKIP");
						else if ((one_hole_row.at(i2).size() - check_size) == 0)
						{
							M_LOG3("No missing point -> SKIP");
						}
						else
						{
							int x_temp, y_temp;
							// setup 2 point on Curve
							//CVertexO p1 = one_hole_row.at(i).at(0);
							//CVertexO p2 = *one_hole_row.at(i).end();
							CVertexO p1 = *one_hole_row.at(i2).begin();
							CVertexO p2 = one_hole_row.at(i2).at(one_hole_row.at(i2).size() - 1);
							M_LOG3("p1 (x,y,z) : (%f,%f,%f)", p1.P().X(), p1.P().Y(), p1.P().Z());
							M_LOG3("p2 (x,y,z) : (%f,%f,%f)", p2.P().X(), p2.P().Y(), p2.P().Z());

							CVertexO p0, p3;

							//setup start point
							x_temp = p1.P().X() - 1;
							y_temp = p1.P().Y();
							if (xy[Cal_index_k(x_temp, y_temp, expanded_y, k_max)] != NULL)
								p0 = *xy[Cal_index_k(x_temp, y_temp, expanded_y, k_max)];
							else
							{
								p0.P().X() = p1.P().X() - 1;
								p0.P().Y() = p1.P().Y();
								p0.P().Z() = p1.P().Z();
							}
							M_LOG3("p0 (x,y,z) : (%f,%f,%f)", p0.P().X(), p0.P().Y(), p0.P().Z());

							//setup end point
							x_temp = p2.P().X() + 1;
							y_temp = p2.P().Y();
							if (xy[Cal_index_k(x_temp, y_temp, expanded_y, k_max)] != NULL)
								p3 = *xy[Cal_index_k(x_temp, y_temp, expanded_y, k_max)];
							else
							{
								p3.P().X() = p2.P().X() + 1;
								p3.P().Y() = p2.P().Y();
								p3.P().Z() = p2.P().Z();
							}
							M_LOG3("p3 (x,y,z) : (%f,%f,%f)", p3.P().X(), p3.P().Y(), p3.P().Z());

							///*
							//double u = (p1.P().X() - p0.P().X()) / (p3.P().X() - p0.P().X());
							double u = (distance3D(p0, p1)) / (distance3D(p0, p1) + distance3D(p1, p2) + distance3D(p2, p3));
							M_LOG3("u:%f", u);
							//double v = (p2.P().X() - p0.P().X()) / (p3.P().X() - p0.P().X());
							double v = (distance3D(p0, p1) + distance3D(p1, p2)) / (distance3D(p0, p1) + distance3D(p1, p2) + distance3D(p2, p3));
							M_LOG3("v:%f", v);

							CVertexO c1, c2;

							int result = interpolate(p0, u, p1, v, p2, p3, c1, c2);

							M_LOG3("c1 (x,y,z) : (%f,%f,%f)", c1.P().X(), c1.P().Y(), c1.P().Z());
							M_LOG3("c2 (x,y,z) : (%f,%f,%f)", c2.P().X(), c2.P().Y(), c2.P().Z());

							if (result == 1)
							{
								M_LOG3("The calculated control points are successed");

								for (int j = p1.P().X() + 1; j < p2.P().X(); j++)
								{
									//double t = (j - p0.P().X()) / (p3.P().X() - p0.P().X());
									double t;
									double a = -1 * p0.P().X() + 3 * c1.P().X() - 3 * c2.P().X() + p3.P().X();
									double b = 3 * p0.P().X() - 6 * c1.P().X() + 3 * c2.P().X();
									double c = -3 * p0.P().X() + 3 * c1.P().X();
									double d = p0.P().X() - j;
									double t1 = 0, t2 = 0, t3 = 0;
									M_LOG3("\nStart to solve Equation for rows");
									M_LOG3("%fx^3 + %fx^2 + %fx + %f = 0", a, b, c, d);
									if (Compare2Double(a, 0.0) && Compare2Double(b, 0.0))
									{
										t = -1 * d / c;
										M_LOG3("ROW ****");
									}
									else if (Compare2Double(a, 0.0) && !Compare2Double(b, 0.0))
									{
										int case_quadratic = solve_quadratic(b, c, d, t1, t2);
										if (case_quadratic == 2)
										{
											t = t1;
										}
										if (case_quadratic == 1)
										{
											bool case_1 = false, case_2 = false;
											if ((u < t1) && (t1 < v))
											{
												t = t1;
												case_1 = true;
											}
											if ((u < t2) && (t2 < v))
											{
												t = t2;
												case_2 = true;
											}
											if (case_1 && case_2)
												qDebug("ERROR at case_1&2");
										}
										if (case_quadratic == 0)
											qDebug("ERROR when solve quadratic function at ROW");
										M_LOG3("t:%f", t);
										M_LOG3("u:%f < t:%f < v:%f", u, t, v);
									}
									else if (!Compare2Double(a, 0.0))
									{
										int case_cubic = solve_cubic(a, b, c, d, t1, t2, t3);
										if ((case_cubic == 1) || (case_cubic == 3))
										{
											t = t1;
										}
										if (case_cubic == 2)
										{
											bool case_1 = false, case_2 = false, case_3 = false;
											if ((u < t1) && (t1 < v))
											{
												t = t1;
												case_1 = true;
											}
											if ((u < t2) && (t2 < v))
											{
												t = t2;
												case_2 = true;
											}
											if ((u < t3) && (t3 < v))
											{
												t = t3;
												case_3 = true;
											}
											//if (case_1 && case_2) qDebug("ERROR at case_1&2");
											//if (case_2 && case_3) qDebug("ERROR at case_2&3");
											//if (case_3 && case_1) qDebug("ERROR at case_3&1");
										}
										if (case_cubic == 0)
											qDebug("ERROR when solve cubic function at ROW");
										M_LOG3("t:%f", t);
										M_LOG3("u:%f < t:%f < v:%f", u, t, v);
									}

									CVertexO new_point = drawBezier_curve(p0, c1, c2, p3, t);
									// check curve are straingle line or not
									M_LOG3("ROW:");
									//if (z_4_points_less_than_x(p0, p1, p2, p3, new_point, 0.0)) {
									if (z_4_points_less_than_x(p0, p1, p2, p3, new_point, EPSILON))
									{
										new_point.SetUserBit(same_line);
										M_LOG3("This point on the straingle line: %f, %f, %f", new_point.P().X(), new_point.P().Y(), new_point.P().Z());
									}
									//else {
									//M_LOG3("Check straingle line: %f, %f, %f, %f, %f", p0.P().Z(), p1.P().Z(), p2.P().Z(), p3.P().Z(), new_point.P().Z());
									//}
									//new_point.C().SetRGB(0, 255, 0);
									M_LOG3("new_point.P().X():%f", new_point.P().X());
									M_LOG3("new_point.P().Y():%f", new_point.P().Y());
									M_LOG3("new_point.P().Z():%f", new_point.P().Z());
									new_points_list_row.push_back(new_point);
								}

								M_LOG3("The new points have created and added");
							}
							else
							{
								M_LOG3("Unable to calculate Bezier control points.");
							}
							//*/
						}
					}

					////////////////////////////////////////////////////////
					// create new points on col
					///*
					for (unsigned int i2 = 0; i2 < one_hole_col.size(); i2++)
					{
						M_LOG3("\t- At col(i): (%d) -", i2);
						int check_size = one_hole_col.at(i2).begin()->P().X() - one_hole_col.at(i2).at(one_hole_col.at(i2).size() - 1).P().X() + 1;
						if (one_hole_col.at(i2).size() == 1)
							M_LOG3("Col have only one point -> SKIP");
						else if ((one_hole_col.at(i2).size() - check_size) == 0)
						{
							M_LOG3("No missing point -> SKIP");
						}
						else
						{
							int x_temp, y_temp;
							// setup 2 point on Curve
							//CVertexO p1 = one_hole_col.at(i).at(0);
							//CVertexO p2 = *one_hole_col.at(i).end();
							CVertexO p1 = *one_hole_col.at(i2).begin();
							CVertexO p2 = one_hole_col.at(i2).at(one_hole_col.at(i2).size() - 1);
							M_LOG3("p1 (x,y,z) : (%f,%f,%f)", p1.P().X(), p1.P().Y(), p1.P().Z());
							M_LOG3("p2 (x,y,z) : (%f,%f,%f)", p2.P().X(), p2.P().Y(), p2.P().Z());

							CVertexO p0, p3;

							//setup start point
							x_temp = p1.P().X();
							y_temp = p1.P().Y() - 1;
							if (xy[Cal_index_k(x_temp, y_temp, expanded_y, k_max)] != NULL)
								p0 = *xy[Cal_index_k(x_temp, y_temp, expanded_y, k_max)];
							else
							{
								p0.P().X() = p1.P().X();
								p0.P().Y() = p1.P().Y() - 1;
								p0.P().Z() = p1.P().Z();
							}
							M_LOG3("p0 (x,y,z) : (%f,%f,%f)", p0.P().X(), p0.P().Y(), p0.P().Z());

							//setup end point
							x_temp = p2.P().X();
							y_temp = p2.P().Y() + 1;
							if (xy[Cal_index_k(x_temp, y_temp, expanded_y, k_max)] != NULL)
								p3 = *xy[Cal_index_k(x_temp, y_temp, expanded_y, k_max)];
							else
							{
								p3.P().X() = p2.P().X();
								p3.P().Y() = p2.P().Y() + 1;
								p3.P().Z() = p2.P().Z();
							}
							M_LOG3("p3 (x,y,z) : (%f,%f,%f)", p3.P().X(), p3.P().Y(), p3.P().Z());

							//double u = (p1.P().X() - p0.P().X()) / (p3.P().X() - p0.P().X());
							double u = (distance3D(p0, p1)) / (distance3D(p0, p1) + distance3D(p1, p2) + distance3D(p2, p3));
							M_LOG3("u:%f", u);
							//double v = (p2.P().X() - p0.P().X()) / (p3.P().X() - p0.P().X());
							double v = (distance3D(p0, p1) + distance3D(p1, p2)) / (distance3D(p0, p1) + distance3D(p1, p2) + distance3D(p2, p3));
							M_LOG3("v:%f", v);

							CVertexO c1, c2;

							int result = interpolate(p0, u, p1, v, p2, p3, c1, c2);

							M_LOG3("c1 (x,y,z) : (%f,%f,%f)", c1.P().X(), c1.P().Y(), c1.P().Z());
							M_LOG3("c2 (x,y,z) : (%f,%f,%f)", c2.P().X(), c2.P().Y(), c2.P().Z());

							if (result == 1)
							{
								M_LOG3("The calculated control points are successed");

								for (int j = p1.P().Y() + 1; j < p2.P().Y(); j++)
								{
									//double t = (j - p0.P().Y()) / (p3.P().Y() - p0.P().Y());
									double t;
									double a = -1 * p0.P().Y() + 3 * c1.P().Y() - 3 * c2.P().Y() + p3.P().Y();
									double b = 3 * p0.P().Y() - 6 * c1.P().Y() + 3 * c2.P().Y();
									double c = -3 * p0.P().Y() + 3 * c1.P().Y();
									double d = p0.P().Y() - j;
									double t1 = 0, t2 = 0, t3 = 0;
									M_LOG3("\nStart to solve Cubic Equation for cols");
									M_LOG3("%fx^3 + %fx^2 + %fx + %f = 0", a, b, c, d);
									if (Compare2Double(a, 0.0) && Compare2Double(b, 0.0))
									{
										t = -1 * d / c;
										M_LOG3("COL ****");
									}
									else if (Compare2Double(a, 0.0) && !Compare2Double(b, 0.0))
									{
										int case_quadratic = solve_quadratic(b, c, d, t1, t2);
										if (case_quadratic == 2)
										{
											t = t1;
										}
										if (case_quadratic == 1)
										{
											bool case_1 = false, case_2 = false;
											if ((u < t1) && (t1 < v))
											{
												t = t1;
												case_1 = true;
											}
											if ((u < t2) && (t2 < v))
											{
												t = t2;
												case_2 = true;
											}
											if (case_1 && case_2)
												qDebug("ERROR at case_1&2");
										}
										if (case_quadratic == 0)
											qDebug("ERROR when solve quadratic function at COL");
										M_LOG3("t:%f", t);
										M_LOG3("u:%f < t:%f < v:%f", u, t, v);
									}
									else if (!Compare2Double(a, 0.0))
									{
										int case_cubic = solve_cubic(a, b, c, d, t1, t2, t3);
										if ((case_cubic == 1) || (case_cubic == 3))
										{
											t = t1;
										}
										if (case_cubic == 2)
										{
											bool case_1 = false, case_2 = false, case_3 = false;
											if ((u < t1) && (t1 < v))
											{
												t = t1;
												case_1 = true;
											}
											if ((u < t2) && (t2 < v))
											{
												t = t2;
												case_2 = true;
											}
											if ((u < t3) && (t3 < v))
											{
												t = t3;
												case_3 = true;
											}
											//if (case_1 && case_2) qDebug("ERROR at case_1&2");
											//if (case_2 && case_3) qDebug("ERROR at case_2&3");
											//if (case_3 && case_1) qDebug("ERROR at case_3&1");
										}
										if (case_cubic == 0)
											qDebug("ERROR when solve cubic function at COL");
										M_LOG3("t:%f", t);
										M_LOG3("u:%f < t:%f < v:%f", u, t, v);
									}

									CVertexO new_point = drawBezier_curve(p0, c1, c2, p3, t);
									// check curve are straingle line or not
									M_LOG3("COL:");
									//if (z_4_points_less_than_x(p0, p1, p2, p3, new_point, 0.0)) {
									if (z_4_points_less_than_x(p0, p1, p2, p3, new_point, EPSILON))
									{
										new_point.SetUserBit(same_line);
										M_LOG3("This point on the straingle line: %f, %f, %f", new_point.P().X(), new_point.P().Y(), new_point.P().Z());
									}
									//else {
									//M_LOG3("Check straingle line: %f, %f, %f, %f, %f", p0.P().Z(), p1.P().Z(), p2.P().Z(), p3.P().Z(), new_point.P().Z());
									//}
									//new_point.C().SetRGB(0, 255, 0);
									M_LOG3("new_point.P().X():%f", new_point.P().X());
									M_LOG3("new_point.P().Y():%f", new_point.P().Y());
									M_LOG3("new_point.P().Z():%f", new_point.P().Z());
									new_points_list_col.push_back(new_point);
								}

								M_LOG3("The new points have created and added");
							}
							else
							{
								M_LOG3("Unable to calculate Bezier control points.");
							}
						}
					}
					//*/

					M_LOG3("\t\t\t------------DONE ONE HOLE (created new points)------------");
				}

				// get the everage points
				M_LOG3("\t\t\t------------START TO CALCULATE AVERAGE POINT AT THE SAME X,Y------------");
				vector<CVertexO> all_new_points_list;
				vector<CVertexO> all_new_points_list_dup;
				///*
				for (unsigned int i4 = 0; i4 < new_points_list_row.size(); i4++)
				{
					for (unsigned int j4 = 0; j4 < new_points_list_col.size(); j4++)
					{
						int x1 = round_num(new_points_list_row.at(i4).P().X());
						int y1 = round_num(new_points_list_row.at(i4).P().Y());
						int x2 = round_num(new_points_list_col.at(j4).P().X());
						int y2 = round_num(new_points_list_col.at(j4).P().Y());
						if ((x1 == x2) && (y1 == y2))
						{
							CVertexO new_point;
							new_point.P().X() = x1;
							new_point.P().Y() = y1;
							if (new_points_list_row.at(i4).IsUserBit(same_line) && !new_points_list_col.at(j4).IsUserBit(same_line))
							{
								new_point.P().Z() = new_points_list_row.at(i4).P().Z();
								M_LOG3("Choose this point of straingle line of ROW: %f, %f, %f", new_point.P().X(), new_point.P().Y(), new_point.P().Z());
							}
							else if (!new_points_list_row.at(i4).IsUserBit(same_line) && new_points_list_col.at(j4).IsUserBit(same_line))
							{
								new_point.P().Z() = new_points_list_col.at(j4).P().Z();
								M_LOG3("Choose this point of straingle line of COL: %f, %f, %f", new_point.P().X(), new_point.P().Y(), new_point.P().Z());
							}
							else
							{
								new_point.P().Z() = (new_points_list_row.at(i4).P().Z() + new_points_list_col.at(j4).P().Z()) / 2.0;
							}

							new_point.C().SetHSVColor(1.2f, 1.0f, 1.0f);
							all_new_points_list.push_back(new_point);
							all_new_points_list_dup.push_back(new_point);
						}
					}
				}

				// add new points here
				if (!case4)
				{
					for (unsigned int i5 = 0; i5 < all_new_points_list.size(); i5++)
					{
						m.cm.vert.push_back(all_new_points_list.at(i5));
					}
					m.cm.vn = m.cm.vn + all_new_points_list.size();
				}

				//*/
				// debug
				/*
                                for (unsigned int i4 = 0; i4 < new_points_list_row.size(); i4++)
				{
                                int x1 = round_num(new_points_list_row.at(i4).P().X());
                                int y1 = round_num(new_points_list_row.at(i4).P().Y());

				CVertexO new_point;
				new_point.P().X() = x1;
				new_point.P().Y() = y1;
                                new_point.P().Z() = new_points_list_row.at(i4).P().Z();
				new_point.C().SetRGB(0, 255, 0); // green color
                                all_new_points_list.push_back(new_point);
				}*/

				// CHECKBOX 1&2&3&4
				if (case4)
				{
					// REFINE ALL NEW POINTS THAT WE FOUND
					vector<CVertexO> all_new_points_list_refine;

					//---clone to another array to compute & refine---
					CVertexO **xy_dup = new CVertexO *[max_size];
					memset(xy_dup, NULL, max_size * sizeof(CVertexO *));
					for (int i = 0; i < vertex_numbers; i++)
					{
						int a = (int)m.cm.vert[i].P().X();
						int b = (int)m.cm.vert[i].P().Y();
						xy_dup[Cal_index_k(a, b, expanded_y, k_max)] = &m.cm.vert[i];
					}

					//add new points to duplicate mesh & find minimum index
					int curr_insert_p_x = all_new_points_list.at(0).P().X();
					int curr_insert_p_y = all_new_points_list.at(0).P().Y();
					CVertexO new_refine_point;
					for (unsigned int i = 0; i < all_new_points_list.size(); i++)
					{
						int a = all_new_points_list.at(i).P().X();
						int b = all_new_points_list.at(i).P().Y();
						xy_dup[Cal_index_k(a, b, expanded_y, k_max)] = &all_new_points_list.at(i);
						xy_dup[Cal_index_k(a, b, expanded_y, k_max)]->SetUserBit(ins_point);
						if (curr_insert_p_x >= a)
						{
							if (curr_insert_p_y >= b)
							{
								curr_insert_p_x = a;
								curr_insert_p_y = b;
							}
						}
						xy_dup[Cal_index_k(a, b, expanded_y, k_max)] = &all_new_points_list.at(i);
						//M_debug_point(all_new_points_list.at(i));
					}

					M_LOG4("\n --------- BEGINNING REFINE: --------- ");

					// Move next point in 8 connectivity of new points & refine z of normal vector
					int pre_insert_p_x = curr_insert_p_x - 1;
					int pre_insert_p_y = curr_insert_p_y - 1;
					for (unsigned int i = 0; i < all_new_points_list.size(); i++)
					{
						//refine new inserted points
						M_LOG4("\n BEGINNING POINT TO REFINE: %d,%d ", curr_insert_p_x, curr_insert_p_y);
						//all_new_points_list_refine.at(i) = all_new_points_list.at(i);
						new_refine_point.P().X() = curr_insert_p_x;
						new_refine_point.P().Y() = curr_insert_p_y;
						float z_abc = z_refine(xy_dup, expanded_y, k_max,
																	 *xy_dup[Cal_index_k(curr_insert_p_x, curr_insert_p_y, expanded_y, k_max)]);
						new_refine_point.P().Z() = z_abc;
						if (z_abc != xy_dup[Cal_index_k(curr_insert_p_x, curr_insert_p_y, expanded_y, k_max)]->P().Z())
						{
							M_LOG4("****** DIFFERENT ******");
							M_debug_point(*xy_dup[Cal_index_k(curr_insert_p_x, curr_insert_p_y, expanded_y, k_max)]);
							M_debug_point(new_refine_point);
							M_LOG4("****** DIFFERENT ****** \n");
						}

						///////////////////////////
						//Debug new point
                        new_refine_point.C().SetHSVColor(0.6f, 1.0f, 1.0f);
						//M_LOG4("Refine new point:");
						//M_debug_point(*xy_dup[Cal_index_k(curr_insert_p_x, curr_insert_p_y, expanded_y, k_max)]);
						//M_debug_point(new_refine_point);
						all_new_points_list_refine.push_back(new_refine_point);
						/////////////////////////////

						// clear mark insert_point of current
						xy_dup[Cal_index_k(curr_insert_p_x, curr_insert_p_y, expanded_y, k_max)]->ClearUserBit(ins_point);
						M_LOG4("Clear Userbit ins_point:");
						M_debug_point(*xy_dup[Cal_index_k(curr_insert_p_x, curr_insert_p_y, expanded_y, k_max)]);

						// remove the un_adjusted point in current list
						M_LOG4("all_new_points_list_dup.size(): %d", all_new_points_list_dup.size());
						for (unsigned int m = 0; m < all_new_points_list_dup.size(); m++)
						{
							if (((int)all_new_points_list_dup.at(m).P().X() == curr_insert_p_x) &&
									((int)all_new_points_list_dup.at(m).P().Y() == curr_insert_p_y))
							{
								M_LOG4("DELETE POINT:");
								M_debug_point(all_new_points_list_dup.at(m));
								erase(all_new_points_list_dup, all_new_points_list_dup.at(m));
								M_LOG4("all_new_points_list_dup.size(): %d", all_new_points_list_dup.size());
								break;
							}
						}

						// list all points
						M_LOG4("\n LIST REMAIN POINTS:");
						M_LOG4("all_new_points_list_dup.size(): %d", all_new_points_list_dup.size());
						for (unsigned int m = 0; m < all_new_points_list_dup.size(); m++)
						{
							if (m < 3)
							{
								M_debug_point(all_new_points_list_dup.at(m));
							}
						}

						// move next point in list of new points by travel in 8 connect
						M_LOG4("\n Change pre x,y");
						M_LOG4("pre_x, pre_y: %d,%d", pre_insert_p_x, pre_insert_p_y);
						int j = 0;
						for (j = 0; j < 8; j++)
						{
							move_next_p_8_con(pre_insert_p_x, pre_insert_p_y, curr_insert_p_x, curr_insert_p_y);
							if (xy_dup[Cal_index_k(pre_insert_p_x, pre_insert_p_y, expanded_y, k_max)] != NULL)
							{
								if (xy_dup[Cal_index_k(pre_insert_p_x, pre_insert_p_y, expanded_y, k_max)]->IsUserBit(ins_point))
								{
									M_LOG4("pre_x, pre_y move to: %d,%d", pre_insert_p_x, pre_insert_p_y);
									break;
								}
							}
						}

						//Move to next point is same hole or not
						if (!all_new_points_list_dup.empty())
						{
							//moving in the same hole
							if (j < 8)
							{
								M_LOG4("This move is belong same hole. SWAP current & pre");
								int temp_x = curr_insert_p_x;
								int temp_y = curr_insert_p_y;
								curr_insert_p_x = pre_insert_p_x;
								curr_insert_p_y = pre_insert_p_y;
								pre_insert_p_x = temp_x;
								pre_insert_p_y = temp_y;
							}
							//move to next point in another hole
							else
							{
								// find the smallest index of new current list
								M_LOG4("End point. Move to new point of another hole");
								curr_insert_p_x = all_new_points_list_dup.at(0).P().X();
								curr_insert_p_y = all_new_points_list_dup.at(0).P().Y();
								for (unsigned int m = 0; m < all_new_points_list_dup.size(); m++)
								{
									int a = all_new_points_list_dup.at(m).P().X();
									int b = all_new_points_list_dup.at(m).P().Y();
									if (curr_insert_p_x >= a)
									{
										if (curr_insert_p_y >= b)
										{
											curr_insert_p_x = a;
											curr_insert_p_y = b;
										}
									}
								}
								pre_insert_p_x = curr_insert_p_x - 1;
								pre_insert_p_y = curr_insert_p_y - 1;
							}
							M_LOG4("PRE:");
							M_LOG4("pre_x, pre_y move to: %d,%d", pre_insert_p_x, pre_insert_p_y);
							M_LOG4("CURR:");
							M_LOG4("curr_x, curr_y move to: %d,%d", curr_insert_p_x, curr_insert_p_y);
						}
					}

					// add new points here
					for (unsigned int i5 = 0; i5 < all_new_points_list_refine.size(); i5++)
					{
						m.cm.vert.push_back(all_new_points_list_refine.at(i5));
					}
					m.cm.vn = m.cm.vn + all_new_points_list_refine.size();
				}
			}
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////// debug
	/*
	int index_x = par.getInt("x");
	int index_y = par.getInt("y");
	//convert 1-dimention to 2-dimention:
	vertex_numbers = m.cm.vn;
	M_LOG("\n\n Number of point in the file: %d", vertex_numbers);     //count the number of point in the file

	//------clone all mesh to tracking_mesh and find max index of x & y

	max_index_x = 0, max_index_y = 0;
	min_index_x = (int)m.cm.vert[0].P().X();
	min_index_y = (int)m.cm.vert[0].P().Y();
	sta_index_x = 0, sta_index_y = 0;
	//Log("value at x:%d --- y:%d ", min_index_x, min_index_y);

	for (int i = 0; i < vertex_numbers; i++)
	{
	//find max
	if ((int)m.cm.vert[i].P().X() > max_index_x) max_index_x = (int)m.cm.vert[i].P().X();
	if ((int)m.cm.vert[i].P().Y() > max_index_y) max_index_y = (int)m.cm.vert[i].P().Y();
	//find min
	if ((int)m.cm.vert[i].P().X() < min_index_x) min_index_x = (int)m.cm.vert[i].P().X();
	if ((int)m.cm.vert[i].P().Y() < min_index_y) min_index_y = (int)m.cm.vert[i].P().Y();
	}

	//find start
	sta_index_x = min_index_x;
	sta_index_y = max_index_y;
	for (int i = 0; i < vertex_numbers; i++)
	{
	if ((int)m.cm.vert[i].P().X() == min_index_x)
	{
	if ((int)m.cm.vert[i].P().Y() <= sta_index_y)
	{
	sta_index_y = (int)m.cm.vert[i].P().Y();
	}
	}
	}

        //M_LOG("value at x:%d --- y:%d ", min_index_x, min_index_y);
	M_LOG("Max index at x:%d --- y:%d ", max_index_x, max_index_y);

	//----create a clone for tracking expanded with k----
	expanded_x = max_index_x + 2 * k_max;
	expanded_y = max_index_y + 2 * k_max;
	max_size = expanded_x * expanded_y;

	//---clone to another array to compute---
	CVertexO** xy2 = new CVertexO*[max_size];
	memset(xy2, NULL, max_size * sizeof(CVertexO*));
	for (int i = 0; i < vertex_numbers; i++) {
	int a = (int)m.cm.vert[i].P().X();
	int b = (int)m.cm.vert[i].P().Y();
	xy2[Cal_index_k(a, b, expanded_y, k_max)] = &m.cm.vert[i];
	}
	if (xy2[expanded_y * (index_x + k_max) + (index_y + k_max)] != NULL) {
	xy2[expanded_y * (index_x + k_max) + (index_y + k_max)]->C().SetRGB(0, 0, 255); // Blue color
        //M_LOG3("new_points_list_row[0].P().X():%f", new_points_list_row[0].P().X());
        //M_LOG3("new_points_list_row[0].P().Y():%f", new_points_list_row[0].P().Y());
        //M_LOG3("new_points_list_row[0].P().Z():%f", new_points_list_row[0].P().Z());
        //m.cm.vert.push_back(new_points_list_row[0]);
	}
	else {
	M_LOG("There is no record at x(%d) --- y(%d)", index_x, index_y);
	}

	delete[] xy2;
*/
	//end of debug

	delete[] xy;
	md.mm()->updateDataMask(MeshModel::MM_VERTCOLOR);
	return outputValues;
}

QString TestBezierPlugin::filterScriptFunctionName(ActionIDType filterID)
{
	switch (filterID)
	{
	case TEST_BEZIER:
		return QString("TestBezier");
	default:
		assert(0);
	}
	return QString();
}

int TestBezierPlugin::getRequirements(const QAction *action)
{
	switch (ID(action))
	{
	case TEST_BEZIER:
		return MeshModel::MM_VERTCOORD;
	default:
		return 0;
	}
}

int TestBezierPlugin::postCondition(const QAction *action) const
{
	switch (ID(action))
	{
	case TEST_BEZIER:
		return MeshModel::MM_VERTCOLOR;
	default:
		return MeshModel::MM_UNKNOWN;
	}
}

int TestBezierPlugin::getPreConditions(QAction *action) const
{
	switch (ID(action))
	{
	case TEST_BEZIER:
		return MeshModel::MM_VERTCOORD;
	default:
		return 0;
	}
}

FilterPlugin::FilterArity TestBezierPlugin::filterArity(const QAction *filter) const
{
	// switch(ID(filter))
	// {

	// }
	return FilterPlugin::NONE;
}

MESHLAB_PLUGIN_NAME_EXPORTER(TestBezierPlugin)
