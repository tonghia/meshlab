float distance2Points(Point3m p1, Point3m p2) ;
float calcAvgDistance(std::vector<float> v_distance);
Point3m calcAvgPoint(std::vector<int> vIndex, CMeshO& cm);
Point3m findFilledVertByIsosceles(Point3m p1, Point3m p2, float filling_point_distance, Point3m hole_center, float edge_distance);
int solveQuadraticEquation(float a, float b, float c, float &x1, float &x2);
const char* pointToString(Point3m p);
QString PointToQString(Point3m p);
Point3m findHoleCenterPoint(CMeshO& cm, std::vector<int> hole);
float CalcAvgHoleEdge(CMeshO& cm, std::vector<int> hole);
Point3m calcHoleCenter(CMeshO& cm, std::vector<int> hole);
float calcAvgDistanceToCenter(CMeshO& cm, std::vector<int> hole, Point3m centerPoint);
bool checkHoleSize(CMeshO& cm, std::vector<int> hole, float threshold, Point3m centerPoint);
Point3m calcFillingPoint(Point3m boundaryPoint, Point3m centerPoint, float avgEdge, float ratio);
bool checkCurrPointEdgeOk(Point3m curFill, Point3m curBoundary, Point3m prevFill, Point3m prevBoundary);
bool checkNewPrevPointDistance(Point3m newFill, Point3m prevFill, float threshold);
float calcCenterRatio(CMeshO& cm, Point3m center, float avgEdge, std::vector<int> hole, std::vector<float> vratio);
int findStepToCenter(Point3m center, Point3m boundary, float avgEdge);
int findMaxStepToCenter(CMeshO& cm, Point3m center, float avgEdge, std::vector<int> hole);
std::vector<int> rearrangeHole(CMeshO& cm, Point3m center, float avgEdge, std::vector<int> hole);
std::vector<int> reduceHoleByConnectNearby(CMeshO& cm, std::vector<int> hole, float avgEdge, Point3m center);

// Fill hole algorithm
void fillHoleByCenter(CMeshO& cm, std::vector<int> hole, float extra, float ratio);
void fillHoleByCenterRefined(CMeshO& cm, std::vector<int> hole, float extra, float ratio);
void fillHoleByIsoscelesTriangle(CMeshO& cm, std::vector<int> hole, std::vector<float> vDistance, float nextPointDistance);
void fillHoleRingByRing(CMeshO& cm, std::vector<int> hole, float threshold, bool stepByStep, std::vector<float> vRatio, float avgZRatio);
void fillHoleRingByRingRefined(CMeshO& cm, std::vector<int> hole, float threshold, bool stepByStep, std::vector<float> vRatio, float avgZRatio);
