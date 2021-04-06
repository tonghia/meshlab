#include "filter_bezier2.h"

#include <vector>
#include <iostream>
#include <tuple>
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
#define _M_DEBUG
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

// DEBUG EXTRA
#define _M_DEBUG_E
#ifdef _M_DEBUG_E
#define N_LOG_E(...) qDebug(__VA_ARGS__)
#else
#define N_LOG_E(...)
#endif

// Body

FilterBezier2Plugin::FilterBezier2Plugin()
{
    typeList = {
        TEST_BEZIER2
    };

    QCoreApplication *app = QCoreApplication::instance();
    for (ActionIDType tt : types())
    {
        QAction *act = new QAction(filterName(tt), this);
        actionList.push_back(act);

        if (app != nullptr)
        {
            if (tt == TEST_BEZIER2)
            {
                //				act->setShortcut(QKeySequence ("Ctrl+Del"));
                act->setIcon(QIcon(":/images/3d-icon-1.png"));
                //				act->setPriority(QAction::HighPriority);
            }
        }
    }
}

QString FilterBezier2Plugin::pluginName() const
{
    return "FilterBezier2Plugin";
}

QString FilterBezier2Plugin::filterName(ActionIDType filter) const
{
    switch (filter)
    {
    case TEST_BEZIER2:
        return QString("Test Bezier 2 points");
    }
    return QString("Unknown filter");
}

QString FilterBezier2Plugin::filterInfo(ActionIDType filterId) const
{
    switch (filterId)
    {
    case TEST_BEZIER2:
        return tr("Create points and apply Bezier curves");
    }
    return QString("Unknown Filter");
}

FilterBezier2Plugin::FilterClass FilterBezier2Plugin::getClass(const QAction *action) const
{
    switch (ID(action))
    {
    case TEST_BEZIER2:
        return FilterClass(FilterPlugin::PointSet + FilterPlugin::VertexColoring);
    default:
        assert(0);
    }
    return FilterPlugin::Generic;
}

void FilterBezier2Plugin::initParameterList(const QAction *action, MeshModel &m, RichParameterList &parlst)
{
    switch (ID(action))
    {
    case TEST_BEZIER2:
        parlst.addParam(RichBool("step1", true, "External Boundary Extraction", "Case 1"));
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

// return index of one dimention array
int Cal_index_k(int x, int y, int max_y, int k)
{
    return (x + k) * max_y + (y + k);
}

// functions for find external boundary

//find distance between 2 point
float distance(float x1, float y1, float x2, float y2)
{
    return sqrt(pow((x1 - x2), 2) + pow((y1 - y2), 2));
}

// compare 2 float
bool Compare2Float(float A, float B)
{
    float diff = A - B;
    return (diff < EPSILON) && (-diff < EPSILON);
}

//check whether if p_curr(x_center, y_center) is the center point of p1(x1, y1) and p2(x2, y2) or not
bool isPCurrCenter(float x_center, float y_center, float x1, float y1, float x2, float y2)
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

enum N_LINE
{
    LINE_0 = 0,
    LINE_12 = 12,
    LINE_23 = 23,
    LINE_34 = 34,
    LINE_41 = 41
};

/**
 * find match point of line ((x1, y1), (x2, y2)) and line ((x3, y3), (x4, y4))
 */
tuple<float, float, bool> findMatch(int x1, int y1, int x2, int y2, int x3, int y3, int x4, int y4) {
    float intersectionX, intersectionY;

    M_LOG("Find the intersecting point of 2 lines");
    // y = m * x + c
    // with 2 point (x1, y1), (x2, y2) we have a system of linear equations
    // y1 = m * x1 + c; y2 = m * x2 + c;
    // => y2 - y1 = m * (x2 - x1) + c
    // => m = (x2 - x1) / (y2 - y1); c = y1 = m * x1
    // Note: check cases y = c or x = a

    float dx21 = x2 - x1;
    float dy21 = y2 - y1;
    float m21 = dy21 / dx21;
    float c21 = y1 - m21 * x1;

    float dx43 = x4 - x3;
    float dy43 = y4 - y3;
    float m43 = dy43 / dx43;
    float c43 = y3 - m43 * x3;

    M_LOG("Program to find the intersecting point of two lines:");
    M_LOG("the first line is (p_curr,p_prev) ");
    N_LOG_E("x1,y1(%d,%d) \t x2,y2(%d,%d) \t", x1, y1, x2, y2);
    M_LOG("dx21 %f \t dy21 %f \t m21 %f \t c21 %f", dx21, dy21, m21, c21);
    N_LOG_E("x3,y3(%d,%d) \t x4,y4(%d,%d) ", x3, y3, x4, y4);
    M_LOG("dx43 %f \t dy43 %f \t m43 %f \t c43 %f ", x3, y3, x4, y4);

    if (dx21 == 0 && dx43 == 0)
    {
        M_LOG("Parallel with y axis ");
        M_LOG("-------------");
        return {0, 0, false};
    }
    else if (dy21 == 0 && dy43 == 0)
    {
        M_LOG("Parallel with x axis ");
        M_LOG("-------------");
        return {0, 0, false};
    }
    else if (dx21 == 0)
    {
        // line 21: x = x1
        // line 43: y = m2 * x + c43
        intersectionX = x1;
        intersectionY = m43 * intersectionX + c43;
        M_LOG("Intersecting Point: = %f, %f ", intersectionX, intersectionY);
        M_LOG("-------------");
        return {intersectionX, intersectionY, true};
    }
    else if (dx43 == 0)
    {
        // line 21: y = m1 * x + c21
        // line 43: x = x3
        intersectionX = x3;
        intersectionY = m21 * intersectionX + c21;
        M_LOG("Intersecting Point: = %f, %f ", intersectionX, intersectionY);
        M_LOG("-------------");
        return {intersectionX, intersectionY, true};
    }
    else if ((m21 - m43) == 0)
    {
        M_LOG("2 Lines are paralled, no intersection");
        return {0, 0, false};
    }
    else
    {
        // line 1: y = m21 * x + c21
        // line 2: y = m43 * x + c43
        intersectionX = (c21 - c43) / (m21 - m43);
        intersectionY = m21 * intersectionX + c21;
        M_LOG("Intersecting Point: = %f, %f ", intersectionX, intersectionY);
        M_LOG("-------------");
        return {intersectionX, intersectionY, true};
    }
}

/** find the intersection between RNk(p) and half_line(p_curr, p_prev)
 *
 * @param[in] x_prev
 * @param[in] y_prev
 * @param[in] x_curr
 * @param[in] y_curr
 * @param[in] max_y
 * @param[in] k_ins
 *
 * @result tuple<x_result, y_result, line>
 *
 * 2---------3
 * | \  |  / |
 * | -  C  - |
 * | /  |  \ |
 * 1---------4
 * we have 4 lines: 12, 23, 34, 41
 * from point C(x_curr, y_curr) we calculate point 1, 2, 3, 4 with unit distance k_ins
 * from point P(x_prev, y_prev) we have half line CP
 * find the intersection between half line CP and lines 12, 23, 34, 41
 */
tuple<float, float, N_LINE> Comp_p_Prev_p_Curr(
    int x_prev, int y_prev,
    int x_curr, int y_curr,
    int max_y, int k_ins)
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

    float x_result, y_result;
    bool status_intersection;

    //find 4 intersection between (p_curr, p_prev) and sides of square RN(k)

    //***check on the x_line below (Point 1 & 2)***
    M_LOG("Start find intersection of (p_curr, p_prev) and 1-2");
    tie(x_result, y_result, status_intersection) = findMatch(x_curr, y_curr, x_prev, y_prev, x1, y1, x2, y2);
    if (status_intersection)
    {
        bool check = isPCurrCenter(x_curr, y_curr, x_prev, y_prev, x_result, y_result);
        if (check)
        {
            M_LOG("p_curr(%d,%d) is center p_prev(%d,%d) and p_result(%f,%f)", x_curr, y_curr, x_prev, y_prev, x_result, y_result);
            M_LOG("->REJECT point(%f,%f)  ", x_result, y_result);
        }
        else if ((int)x_result == x1 && (y2 >= y_result) && (y_result >= y1))
        {
            M_LOG("This point(%f,%f) belong (Point 1(%d,%d) & 2(%d,%d))", x_result, y_result, x1, y1, x2, y2);
            return {x_result, y_result, LINE_12};
        }
    }

    //***check on the x_line below (Point 2 & 3)***
    M_LOG("Start find intersection of (p_curr, p_prev) and 2-3");
    tie(x_result, y_result, status_intersection) = findMatch(x_curr, y_curr, x_prev, y_prev, x2, y2, x3, y3);
    if (status_intersection)
    {
        bool check = isPCurrCenter(x_curr, y_curr, x_prev, y_prev, x_result, y_result);
        if (check)
        {
            M_LOG("p_curr(%d,%d) is center p_prev(%d,%d) and p_result(%f,%f)", x_curr, y_curr, x_prev, y_prev, x_result, y_result);
            M_LOG("->REJECT point(%f,%f)  ", x_result, y_result);
        }
        else if ((int)y_result == y2 && (x3 >= x_result) && (x_result >= x2))
        {
            M_LOG("This point(%f,%f) belong (Point 2(%d,%d) & 3(%d,%d))", x_result, y_result, x2, y2, x3, y3);
            return {x_result, y_result, LINE_23};
        }
    }

    //***check on the x_line below (Point 3 & 4)***
    M_LOG("Start find intersection of (p_curr, p_prev) and 3-4");
    tie(x_result, y_result, status_intersection) = findMatch(x_curr, y_curr, x_prev, y_prev, x3, y3, x4, y4);
    if (status_intersection)
    {
        bool check = isPCurrCenter(x_curr, y_curr, x_prev, y_prev, x_result, y_result);
        if (check)
        {
            M_LOG("p_curr(%d,%d) is center p_prev(%d,%d) and p_result(%f,%f)", x_curr, y_curr, x_prev, y_prev, x_result, y_result);
            M_LOG("->REJECT point(%f,%f)  ", x_result, y_result);
        }
        else if ((int)x_result == x3 && (y3 >= y_result) && (y_result >= y4))
        {
            M_LOG("This point(%f,%f) belong (Point 3(%d,%d) & 4(%d,%d))", x_result, y_result, x3, y3, x4, y4);
            return {x_result, y_result, LINE_34};
        }
    }

    //***check on the x_line below (Point 4 & 1)***
    M_LOG("Start find intersection of (p_curr, p_prev) and 4-1");
    tie(x_result, y_result, status_intersection) = findMatch(x_curr, y_curr, x_prev, y_prev, x4, y4, x1, y1);
    if (status_intersection)
    {
        bool check = isPCurrCenter(x_curr, y_curr, x_prev, y_prev, x_result, y_result);
        if (check)
        {
            M_LOG("p_curr(%d,%d) is center p_prev(%d,%d) and p_result(%f,%f)", x_curr, y_curr, x_prev, y_prev, x_result, y_result);
            M_LOG("->REJECT point(%f,%f)  ", x_result, y_result);
        }
        else if ((int)y_result == y4 && (x4 >= x_result) && (x_result >= x1))
        {
            M_LOG("This point(%f,%f) belong (Point 4(%d,%d) & 1(%d,%d))", x_result, y_result, x4, y4, x1, y1);
            return {x_result, y_result, LINE_41};
        }
    }

    assert(0);
}

/** round Coordinate of a point by move it to the given line
 *
 * @param[in] x x-coordinate value
 * @param[in] y y-coordinate value
 * @param[in] line
 *
 * @result tuple<x_result, y_result>
 *
 * 2---------3
 * | \  |  / |
 * | -  C  - |
 * | /  |  \ |
 * 1---------4
 * line 12, 23: round the coordinate value by converting float to int
 * line 34: adjust y-coordinate value
 * line 41: adjust x-coordinate value
 */
tuple<int, int> roundByLine(float x, float y, N_LINE line) {
    if (line == LINE_12)
    {
        return {(int)x, (int)y};
    }
    else if (line == LINE_23)
    {
        return {(int)x, (int)y};
    }
    else if (line == LINE_34)
    {
        int newY = ((y - (int)y) > EPSILON) ? ((int)y + 1) : (int)y;
        return {(int)x, newY};
    }
    else if (line == LINE_41)
    {
        int newX = ((x - (int)x) > EPSILON) ? ((int)x + 1) : (int)x;
        return {newX, (int)y};
    }

    assert(0);
}

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

//find distance between 2 point in 3D
double distance3D(CVertexO p1, CVertexO p2)
{

    return sqrt(pow((p1.P().X() - p2.P().X()), 2) + pow((p1.P().Y() - p2.P().Y()), 2) + pow((p1.P().Z() - p2.P().Z()), 2));
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
    if (isPCurrCenter(x_RNk, y_RNk, x1, y1, x2, y2) && (y_RNk < y2))
    {
        M_LOG("End check (%d, %d) and move to next point (%d, %d)", x_RNk, y_RNk, x_RNk, y_RNk + 1);
        y_RNk++;
        return true;
    }
    //check on side 2-3
    M_LOG("Check x_RNk,y_RNk (%d, %d) on line 2-3", x_RNk, y_RNk);
    if (isPCurrCenter(x_RNk, y_RNk, x2, y2, x3, y3) && (x_RNk < x3))
    {
        M_LOG("End check (%d, %d) and move to next point (%d, %d)", x_RNk, y_RNk, x_RNk + 1, y_RNk);
        x_RNk++;
        return true;
    }
    //check on side 3-4
    M_LOG("Check x_RNk,y_RNk (%d, %d) on line 3-4", x_RNk, y_RNk);
    if (isPCurrCenter(x_RNk, y_RNk, x3, y3, x4, y4) && (y_RNk > y4))
    {
        M_LOG("End check (%d, %d) and move to next point (%d, %d)", x_RNk, y_RNk, x_RNk, y_RNk - 1);
        y_RNk--;
        return true;
    }
    //check on side 4-1
    M_LOG("Check x_RNk,y_RNk (%d, %d) on line 4-1", x_RNk, y_RNk);
    if (isPCurrCenter(x_RNk, y_RNk, x4, y4, x1, y1) && (x_RNk > x1))
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
std::map<std::string, QVariant> FilterBezier2Plugin::applyFilter(
    const QAction *action,
    const RichParameterList &par,
    MeshDocument &md,
    unsigned int & /*postConditionMask*/,
    vcg::CallBackPos * /*cb*/)
{
    std::map<std::string, QVariant> outputValues;

    MeshModel *m = md.mm();
    CMeshO & cm = m->cm;
    int numberSelectedFaces = cm.sfn;
    int numberSelectedVertices = cm.svn;
    int numberPolygonalVertices = cm.pvn;
    int numberPolygonalFaces = cm.pfn;
    N_LOG_E("number of selected faces: %d /n selected vertices: %d /n polygonal vertices: %d /n polygonal faces: %d",
            numberSelectedFaces,
            numberSelectedVertices,
            numberPolygonalVertices,
            numberPolygonalFaces);

    int meshVertexNumber = cm.VN();
    M_LOG("\n\n Number of vertex of the mesh: %d", meshVertexNumber);

    // init max, min of X, Y of all points by first point of mesh
    auto p0 = cm.vert[0].P();
    int maxX = (int)p0.X();
    int minX = (int)p0.X();
    int maxY = (int)p0.Y();
    int minY = (int)p0.Y();
    M_LOG("First value at index (0,0) x:%d --- y:%d ", minX, minY);
    // init start point X, Y
    int startX = (int)p0.X();
    int startY = (int)p0.Y();
    // find max, min of X and Y value of all points in the current mesh
    for (int i = 1; i < meshVertexNumber; i++)
    {
        auto point = cm.vert[i].P();
        if ((int)point.X() > maxX)
        {
            maxX = (int)point.X();
        }
        if ((int)point.Y() > maxY)
        {
            maxY = (int)point.Y();
        }
        if ((int)point.X() < minX)
        {
            minX = (int)point.X();
            // with all points at min X, find a point with min Y
            startX = minX;
            if ((int)point.Y() <= startY)
            {
                startY = (int)point.Y();
            }
        }
        if ((int)point.Y() < minY)
        {
            minY = (int)point.Y();
        }
    }
    N_LOG_E("maxX: %d /n minX: %d /n maxY: %d /n minY: %d /n startX: %d /n startY: %d",
            maxX, minX, maxY, minY, startX, startY);
    M_LOG("Min value at x:%d --- y:%d ", minX, minY);
    M_LOG("Max index at x:%d --- y:%d ", maxX, maxY);

    bool step1 = par.getBool("step1");
    // bool case2 = par.getBool("case2");
    // bool case3 = par.getBool("case3");
    // bool case4 = par.getBool("case4");


    //----create a clone for tracking expanded with k----
    int expanded_x = maxX + 2 * k_max;
    const int expanded_y = maxY + 2 * k_max;
    int max_size = expanded_x * expanded_y;
    int *track_mesh = new int[max_size];
    memset(track_mesh, NULL, max_size * sizeof(int));

    //---start to clone---
    for (int i = 0; i < meshVertexNumber; i++)
    {
        int x = (int)cm.vert[i].P().X();
        int y = (int)cm.vert[i].P().Y();
        track_mesh[Cal_index_k(x, y, expanded_y, k_max)] = norm_p;
    }

    //---clone to another array to compute---
    CVertexO **xy = new CVertexO *[max_size];
    memset(xy, NULL, max_size * sizeof(CVertexO *));
    for (int i = 0; i < meshVertexNumber; i++)
    {
        int a = (int)cm.vert[i].P().X();
        int b = (int)cm.vert[i].P().Y();
        xy[Cal_index_k(a, b, expanded_y, k_max)] = &cm.vert[i];
    }

    //---create list of points was tracked---
    CVertexO *tracked_points = new CVertexO[meshVertexNumber];
    int tracked_index = 0;
    int total_tracked = 0;

    // extract external boundary
    if (step1)
    {
        bool stop = false;
        bool find = false;
        int k = 1;

        // Setup first, current point
        int x_First = startX, y_First = startY;
        int x_Curr = startX, y_Curr = startY;

        // Setup the previous point
        int x_Prev = x_Curr;
        int y_Prev = y_Curr - 1;

        float x_Prev_k, y_Prev_k;
        int x_pre_k, y_pre_k;

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
            N_LINE line_k;
            M_LOG("\t ** Begin to find to intersection (p_Prev_k) between line (p_curr,p_prev) and square (RNk) **");
            tie(x_Prev_k, y_Prev_k, line_k) = Comp_p_Prev_p_Curr(x_Prev, y_Prev, x_Curr, y_Curr, expanded_y, k);
            M_LOG("\t ***** x_Prev_k, y_Prev_k (%f,%f) *****", x_Prev_k, y_Prev_k);
            find = false;
            tie(x_pre_k, y_pre_k) = roundByLine(x_Prev_k, y_Prev_k, line_k);

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
                for (int i = k; i <=k; i++)
                {
                    //q_i = intersection between RNi(p) and discrete segment(p, q)
                    float x_qi, y_qi;

                    M_LOG("\n\t Find intersection between RNi(p)-[i=%d] and discrete segment p_prev_k, p_curr(%d,%d),(%d,%d).",
                          i, x_pre_k, y_pre_k, x_Curr, y_Curr);

                    N_LINE line_i;
                    tie(x_qi, y_qi, line_i) = Comp_p_Prev_p_Curr(x_pre_k, y_pre_k, x_Curr, y_Curr, expanded_y, i);
                    const auto [x_i, y_i] = roundByLine(x_qi, y_qi, line_i);

                    M_LOG("\t\t\t ** (i=%d)The intersection between RNi(p) and discrete segment(p, q) is: (%f,%f) **",
                          i, x_qi, y_qi);
                    M_LOG("\t\t\t ** (i=%d)The intersection between RNi(p) and discrete segment(p, q) is: (%d,%d) **",
                          i, x_i, y_i);

                    M_LOG("Min index (%d,%d); Max index (%d,%d)", minX, minY, maxX, maxY);

                    // check if out of range => skip
                    if ((minX <= x_i) && (x_i <= maxX) && (minY <= y_i) && (y_i <= maxY))
                    {
                        const int index_k = Cal_index_k(x_i, y_i, expanded_y, k_max);
                        //skip some iterator when meet tracked points
                        if ((k > 1) && pass_1st_prev)
                        {
                            if ((track_mesh[index_k] == ex_bound) ||
                                ((x_i == x_Prev) && (y_i == y_Prev)))
                            {
                                M_LOG("REJECT (%d, %d) (k > 1) when check this point be tracked", x_i, y_i);
                                break;
                            }
                            else
                            {
                                pass_1st_prev = false;
                            }
                        }

                        // IF CASE of case 4: tracking first external boundary point (belong i)
                        if (!ex_1st_track &&
                            (track_mesh[index_k] == ex_bound))
                        {
                            x_ex_1st_track = x_i;
                            y_ex_1st_track = y_i;
                            ex_1st_track = true;
                            M_LOG("\t SET (%d, %d) is the first tracked point meet", x_ex_1st_track, y_ex_1st_track);
                        }

                        // no point -> exit this iterator
                        if (xy[index_k] == NULL)
                        {
                            M_LOG("At (%d, %d), there is no point", x_i, y_i);
                            continue;
                        }

                        // pass all requirement => choose this point
                        if ((track_mesh[index_k] != ex_bound) &&
                            (xy[index_k] != NULL))
                        {
                            if (miss_1_of_4_neigh(xy, x_i, y_i, expanded_y, maxX, maxY, minX, minY, k_max))
                            {
                                M_LOG("\t\t\t *** Found the next point (%d, %d) ***", x_i, y_i);
                                // p_Prev = p;
                                x_Prev = x_Curr;
                                y_Prev = y_Curr;
                                M_LOG("\t\t\t **** New point p_prev (%d, %d) ***", x_Prev, y_Prev);
                                // p_Next = q_i;
                                x_Curr = x_i;
                                y_Curr = y_i;
                                M_LOG("\t\t\t *** New point p_curr (%d, %d) *** asdf", x_Curr, y_Curr);

                                // set info
                                set_point(tracked_points, total_tracked, x_i, y_i);
                                M_LOG("\t\t\t a11111");

                                ex_1st_track = false;
                                M_LOG("\t\t\t a22222");

                                tracked_index = total_tracked;
                                M_LOG("\t\t\t a3");

                                track_mesh[index_k] = ex_bound;
                                M_LOG("\t\t\t a4");

                                xy[index_k]->C().SetHSVColor(0, 1.0f, 1.0f); // red color

                                M_LOG("\t\t\t a5");

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
                        M_LOG("At (%d, %d), out of range", x_i, y_i);
                        pass_1st_prev = false; // set this when out of range
                        break;
                    }
                }

                // Case 2: If not found any point -> set first tracked point is next point
                if ((x_pre_k == x_pre_k_1st) &&
                    (y_pre_k == y_pre_k_1st) &&
                    ex_1st_track && k == k_max)
                {
                    M_LOG("\t\t\t *** Set the tracked point to new point (%d, %d) ***", x_ex_1st_track, y_ex_1st_track);
                    //p_Prev = p;
                    x_Prev = x_Curr;
                    y_Prev = y_Curr;
                    M_LOG("\t\t\t *** New point p_prev (%d, %d) ***", x_Prev, y_Prev);
                    //p_Next = q_i;
                    x_Curr = x_ex_1st_track;
                    y_Curr = y_ex_1st_track;
                    M_LOG("\t\t\t *** New point p_curr (%d, %d) ***", x_Curr, y_Curr);

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
                    qDebug("\t\t\t *** ex_1st_track (%d,%d) ***", x_ex_1st_track, y_ex_1st_track);
                    qDebug("\t\t\t *** pre_k_1st (%d,%d) ***", x_pre_k_1st, y_pre_k_1st);
                    qDebug("\t\t\t *** pre_k (%d,%d) ***", x_pre_k, y_pre_k);
                    qDebug("\t\t\t *** ex_1st_track (%d) ***", ex_1st_track);
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
                // M_LOG("total_tracked: %d", total_tracked);

                int max_loop_2 = 4 * (2 * k + 1) + 2;
                stop_loop2++;
                M_LOG("(%d) STOP LOOP HERE, (max: %d)", stop_loop2, max_loop_2);
                if (stop_loop2 == max_loop_2)
                {
                    find = true;
                    stop = true;
                }
            }
            // stop of while(!find)
            // if (p.P() == p_First.P())
            if (x_Curr == x_First && y_Curr == y_First)
            {
                M_LOG("\t\t\t *** END ***");
                stop = true;
            }
        }
        // end while (!stop)
    }
    // end if (step1)


    md.mm()->updateDataMask(MeshModel::MM_VERTCOLOR);
    return outputValues;
}

QString FilterBezier2Plugin::filterScriptFunctionName(ActionIDType filterID)
{
    switch (filterID)
    {
    case TEST_BEZIER2:
        return QString("TestBezier");
    default:
        assert(0);
    }
    return QString();
}

int FilterBezier2Plugin::getRequirements(const QAction *action)
{
    switch (ID(action))
    {
    case TEST_BEZIER2:
        return MeshModel::MM_VERTCOORD;
    default:
        return 0;
    }
}

int FilterBezier2Plugin::postCondition(const QAction *action) const
{
    switch (ID(action))
    {
    case TEST_BEZIER2:
        return MeshModel::MM_VERTCOLOR;
    default:
        return MeshModel::MM_UNKNOWN;
    }
}

int FilterBezier2Plugin::getPreConditions(QAction *action) const
{
    switch (ID(action))
    {
    case TEST_BEZIER2:
        return MeshModel::MM_VERTCOORD;
    default:
        return 0;
    }
}

FilterPlugin::FilterArity FilterBezier2Plugin::filterArity(const QAction *filter) const
{
    // switch(ID(filter))
    // {

    // }
    return FilterPlugin::NONE;
}

MESHLAB_PLUGIN_NAME_EXPORTER(FilterBezier2Plugin)
