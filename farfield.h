// Copyright (C) M.Iakobovski, E.Kornilina, A.Morozov, 2010, v.1

#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>
#include <limits>
#include <fstream>
#include "funcs.h"

using std::cout;
using std::endl;
using std::cin;
using std::vector;
using std::numeric_limits;
using std::ifstream;
using std::ofstream;

const double kSpeedOfSound = 1; //331.46; // c
const double PI = 3.14159;


typedef size_t CoordNumber;
//typedef vector <CoordNumber> PointsNumbers; //для локальных объектов координата - unsigned номер
//									 //в общем массиве

double GetVectorNorma(CoordXYZ &x);
double GetScalarProduct(CoordXYZ &a,CoordXYZ &b);
CoordXYZ GetVectorProduct(CoordXYZ &a,CoordXYZ &b);

// --------------------------------------------

#define kNumberOfCoord 3
#define kNumberOfTetVertexes 4
const int kNumberOfSurfaceVertexes = 3;
//const long kNumberOfAllVertexes = 10;

#define kva(a) ((a)*(a))
#define delta_func(i,j) ((i)==(j)?(1):(0))


// --------------------------------------------


class ControlSurface{
public:
	ControlSurface(double r_sphere_ = 0.1, int k1_lat_=100, int k2_lon_=2):r_sphere(r_sphere_)
	{
		current_cell_number = -1;
		index_lon = -1;
		index_lat = -1;
		k1_lat=k1_lat_;
        k2_lon=k2_lon_;
	}
private:
	double cell_square;
	CoordXYZ cell_normal;
	CoordXYZ cell_mass_centre;
	size_t number_of_cells;
	int current_cell_number;
	double GetSquare(double phi1, double phi2);
	//size_t k1_lat, k2_lon;
	int index_lat, index_lon;
	int k1_lat;
	int k2_lon;

public:
	double r_sphere;
	void SetMeshParams(size_t k1, size_t k2){
		k1_lat = k1;
		k2_lon = k2;
	}
	int GetNextCellParams(	double &cell_square,
							CoordXYZ &cell_normal,
							CoordXYZ &cell_mass_centre);
	void RewindCells (); // перемотать в начало
};

class DistantField{
protected:
	CoordXYZ w;
	CoordXYZ control_point, position_of_source;

	CoordXYZ ModifyPositionByFlow(CoordXYZ x0);
	double ModifyTimeByPosition(CoordXYZ x, double t);
    CoordXYZ ModifyNormalByFlow(CoordXYZ n0, CoordXYZ R0, CoordXYZ R);
	int k1_lat;
	int k2_lon;


public:
	DistantField(){}
	DistantField(CoordXYZ &w_, double omega_, CoordXYZ &control_point_,
		CoordXYZ & position_of_source_, double r_sphere_ = 1., int k1_lat_=100, int k2_lon_=2)
	{
		w.resize(3);
		w = w_;
		control_point.resize(3);
		control_point = control_point_;
		position_of_source.resize(3);
		position_of_source = position_of_source_;
		k1_lat=k1_lat_;
        k2_lon=k2_lon_;
        omega=omega_;
        r_sphere=r_sphere_;
	}
	double CountIntegral(double time_in_control_point);
	double CountBigIntegral(double time_in_control_point);
	protected:
		double r_sphere, omega;

protected:
	double Get_tc(CoordXYZ R0);
};

//------------------------------------------------

class SourceOfNoise : public DistantField{
public:
	SourceOfNoise(){}
	SourceOfNoise(double omega_, CoordXYZ &position_of_source_, CoordXYZ w_)
	{
		position_of_source.resize(3);
		position_of_source = position_of_source_;
		w.resize(3);
		w = w_;
        omega=omega_;
	}
	double omega;
	CoordXYZ position_of_source;
	CoordXYZ GetUDerivative(CoordXYZ &position_of_target, double time_in_target);
	CoordXYZ GetU(CoordXYZ &position_of_target, double time_in_target);
	double GetRhoDerivative(CoordXYZ &position_of_target, double time_in_target);
	double GetRho(CoordXYZ &position_of_target, double time_in_target);

    double GetRhoA(CoordXYZ &x0, double t );

};


CoordXYZ minus2(CoordXYZ a, CoordXYZ b);
CoordXYZ plus2(CoordXYZ a, CoordXYZ b);
CoordXYZ smult2(double a, CoordXYZ b);

void PutTrgToWff(CoordXYZ xtrg[3]);

double det4(CoordXYZ &r0,CoordXYZ &r1,CoordXYZ &r2,CoordXYZ &r3);

//double CalcIntegralOfFullTrg(Coord xtrg[3], double T, SourceOfNoise &source_of_noise, CoordXYZ &control_point);
double CalcIntegralOfFullTrg(SurfacePoints p, PointsXYZ &m, double T, SourceOfNoise &source_of_noise, CoordXYZ &control_point, int &err);

//double CalcIntegralOfSmallTrg(CoordXYZ xtrg[3], double time, BasePoints &bp, CoordXYZ &control_point);

double GetVectorNorma1(CoordXYZ x);
double ModifyTimeByPosition(CoordXYZ x,double t);

int GetCellParams(  CoordXYZ r[3], // массив координат вершин треугольника
                    double &cell_square,
                    CoordXYZ &cell_normal,
					CoordXYZ &cell_mass_centre);

void test_of_basepoint(void);
void putvect(CoordXYZ v);

void CheckTrgSurface(Surfaces &s, PointsXYZ &m);


// ------------------------------------------------

class BasePoints
{
    #define kbase 4
    CoordXYZ rpyr[kbase]; // координаты опорных вершин
    int ip[4]; // номера опорных вершин

    double a_line_form[kbase][kbase]; // коэффициенты линейной формы
                                       // f=GetVal(r,f)

    int solv4eq_ort(int i_of_non_zero_ort);

public:

    BasePoints(){}
    BasePoints(CoordXYZ rpyr[kbase], int ip[4], int &err);
    double GetVal(CoordXYZ r, int ip[4], CoordXYZ f);
};

// -----------------------------------------------

double CalcIntegralOfSmallTrg(CoordXYZ xtrg[3], double T, SourceOfNoise &source_of_noise, CoordXYZ &control_point);
double CalcIntegralOfSmallTrgBp(CoordXYZ xtrg[3], double time, BasePoints &bp, CoordXYZ &control_point);

