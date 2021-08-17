struct HoleVertInfo {
	std::vector<int> vHoleVertIndex;
	std::vector<std::vector<int>> vSetExpVertIndex;
	std::vector<Point3m> vExpPoint;
	std::vector<float> vZChange;
	std::vector<float> vSlope;
};

struct HoleVertData {
    int holeVertIndex;
	std::vector<int> setExpVertIndex;
	Point3m expPoint;
	float zChange;
	float slope;
};

float distance2Points(Point3m p1, Point3m p2);
Point3m calcAvgPoint(std::vector<int> vIndex, CMeshO& cm);
const char* pointToString(Point3m p);
float calcAvgDistanceToCenter(CMeshO& cm, std::vector<int> hole, Point3m centerPoint);
bool checkHoleSize(CMeshO& cm, std::vector<int> hole, float threshold, Point3m centerPoint);
Point3m calcFillingPoint(Point3m boundaryPoint, Point3m centerPoint, float avgEdge, float ratio);
bool checkCurrPointEdgeOk(Point3m curFill, Point3m curBoundary, Point3m prevFill, Point3m prevBoundary);
bool checkNewPrevPointDistance(Point3m newFill, Point3m prevFill, float threshold);
float calcCenterRatio(CMeshO& cm, Point3m center, float avgEdge, std::vector<int> hole, std::vector<float> vratio);
int findStepToCenter(Point3m center, Point3m boundary, float avgEdge, float&);
int findMaxStepToCenter(CMeshO& cm, Point3m center, float avgEdge, std::vector<int> hole);
std::vector<int> rearrangeHole(CMeshO& cm, Point3m center, float avgEdge, std::vector<int> hole);
std::vector<int> reduceHoleByConnectNearby(CMeshO& cm, std::vector<int> hole, float avgEdge, Point3m center);
float calcCenterZChange(CMeshO& cm, Point3m center, float avgEdge, std::vector<int> hole, std::vector<float> vZChange);

Point3m CalcHoleCenter(CMeshO& cm, std::vector<int> hole);
QString PointToQString(Point3m p);
float CalcAvgHoleEdge(CMeshO& cm, std::vector<int> hole);

// Fill hole algorithm
void FillHoleByCenter(CMeshO& cm, std::vector<int> hole, float extra, float ratio);
void FillHoleByCenterRefined(CMeshO& cm, std::vector<int> hole, float extra, float ratio);
void FillHoleRingByRingRefined(CMeshO& cm, std::vector<int> hole, float startAvgEdge, Point3m holeCenter, bool stepByStep, std::vector<float> vZChange, float adjustRatio);
