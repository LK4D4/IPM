#ifndef GRAPH_H_INCLUDED
#define GRAPH_H_INCLUDED
#include "funcs.h"
/*struct Vertices {
	CoordXYZ grid;
	int TriangleNumber;
};*/
struct Rib { //данные о ребрах
	CoordXYZ a, b; //координаты вершин, упорядочены по модулю
	int TriangleNumber; //номер треугольника
};
struct SimplifiedRib {
	int RibNumber;
	int TriangleNumber;
};
vector<SimplifiedRib> CreateInitialStructure(Surfaces& s, PointsXYZ &m);
bool CheckSide(vector<SimplifiedRib>);
bool CheckGraphCoherency(vector<SimplifiedRib>, int);
#endif
