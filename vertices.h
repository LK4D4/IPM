// Copyright (C) M.Iakobovski, E.Kornilina, A.Morozov, 2010, v.1
#define nOfFuncs 5 // rho rhoder uder[3]
#define nOfsTimeStrips 30

class Vertex
{
public:
    Vertex(){}

    int itime; // число слоёв отставания по времени
    double a[nOfFuncs]; // весовые множители


};

class Vertices
{
    double timestep; // шаг дискретизации по времени

    int n; // общее число вершин

    Vertex **vert; // собственно данные [номер вершины][номер слоя]

    int imintime;
    int imaxtime;
    int max_strips_per_vertex;

    int incorrect;

public:

    Vertices(){incorrect=0;imintime=2; imaxtime=1; max_strips_per_vertex=0;}

int Init(int n_, double timestep_);

double Svertka(Surfaces &s, PointsXYZ &m, SourceOfNoise &sn, double time);

int StripOfTime(double time);
int AddData(int iv, int i_func_num, double value_of_func, double time);
int PutData(char *name);
};

