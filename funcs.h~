// Copyright (C) M.Iakobovski, E.Kornilina, A.Morozov, 2010, v.1

#ifndef FUNCS_H_INCLUDED
#define FUNCS_H_INCLUDED

#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <math.h>

#define Copyright(name,x,y,z) \
static char *name = "\r\n" #x " - Copyright (C) M.Iakobovski, E.Kornilina, A.Morozov, " #y " V. " #z " ";

#include "param.h"

#define kva(a) ((a)*(a))

using std::cout;
using std::endl;
using std::cin;
using std::vector;

/*struct Coord{
	double x,y,z;
};*/
typedef vector <double> CoordXYZ;
typedef vector <CoordXYZ> PointsXYZ;

struct SurfaceCoord{
	int a_index, b_index;
	double delta;
};
typedef vector <SurfaceCoord> SurfacePoints;
typedef vector <SurfacePoints> Surfaces;

typedef vector <int> Tetras;

struct TriangleParams{
	double area;
	CoordXYZ center;
	CoordXYZ normal;
};

//typedef vector <double> CoordXYZ;

#define kNumberOfCoord 3
#define kNumberOfTetVertexes 4
//const int kNumberOfSurfaceVertexes = 3;
const int kNumberOfProcTetra = 1000000;
//параметры поверхности
const int SurfaceType = 0; //0 - эллипсоид; 1 - параллелепипед; 2 - конус
//эллипсоид:

//параллелепипед:
const double p_x0 =  5.35;
const double p_y0 =  5.23;
const double p_z0 =  5.26;
const double p_x1 =  11.42;
const double p_y1 =  11.12;
const double p_z1 =  11.57;

//конус:
const double k_x0 = 5.11467342;
const double k_y0 = 8.;
const double k_z0 = 8.;
const double k_x1 = 11.007743285;
const double k_y1 = k_y0;
const double k_z1 = k_z0;
const double k_a0 = 3;
const double k_b0 = 3;
const double k_a1 = 3;
const double k_b1 = 3;

/*bool ellipse(Coord a);

bool parallel(Coord a);

double mydistance(Coord a, Coord b);

bool isIn(Coord a);

bool isIntersect(Coord n, Coord m);

Coord intersectSearch(Coord a, Coord b);
*/

double keisu(CoordXYZ a);
SurfacePoints FindSurface(Tetras a, PointsXYZ& m);
TriangleParams GetTriangleParams(SurfacePoints &s, PointsXYZ &m);
void CreateWff(char *name, Surfaces &s, PointsXYZ &m);
double CalcIntegral(Surfaces &s, PointsXYZ &m);
int RoundControl(SurfacePoints &s, PointsXYZ &m);
Surfaces DivideSquares(SurfacePoints &sur, PointsXYZ &m);
CoordXYZ GetCoordSurf(SurfaceCoord c, PointsXYZ &m);
void PutSurface(char *name, Surfaces &s, PointsXYZ &m);
int tetra_to_wff(vector<Tetras> &tetrahedron, PointsXYZ &mesh);

#endif // FUNCS_H_INCLUDED
