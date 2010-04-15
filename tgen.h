// Copyright (C) M.Iakobovski, 2010, v.1
#include "./funcs.h"
int meshgen(char * name, int n);
void tgen(std::ofstream &fout, CoordXYZ &p1, CoordXYZ &p2, CoordXYZ &p3, int &k, int n);

class Edge
{
public:

    int i1,i2; // номера вершин - концов ребра
    int ic;    // номер вершины на середине ребра, изначально -1
    int ie1,ie2; // номера рёбер - потомков ребра, изначально -1

    Edge()
    {
        i1=-1;
        i2=-1;
        ic=-1;
        ie1=-1;
        ie2=-1;
    }

    Edge(int i1_, int i2_)
    {
        i1=i1_;
        i2=i2_;
        ic=-1;
        ie1=-1;
        ie2=-1;
    }

    void put(void)
    {
        printf("i1=%d i2=%d ic=%d ie1=%d ie2=%d\n",i1,i2,ic,ie1,ie2);
    }
};

typedef vector <int> iTrg;
typedef vector <iTrg> iTrgs;
typedef vector <Edge> Edges;

class RefTriangles
{
    double R1,R2; // радиусы сфер размещения узлов
    int nsteps; // число шагов измельчения

    Edges edges;
    PointsXYZ points;
    iTrgs itrgs;

public:


    RefTriangles()
    {
        R1=.98;
        R2=1.02;
//        edges.resize(1000000);
    }

    int meshgen(char * name, int n);
    void tgen(CoordXYZ &p1, CoordXYZ &p2, CoordXYZ &p3, int &k, int level);

    int meshgent(char * name, int n);
    void tgent(int ie1, int ie2, int ie3, int level); // измельчение одного треугольника

    int AddEdge(int i1, int i2); // возвращает номер ребра
    int AddPoint(CoordXYZ p); // возвращает номер точки
    void AddiTrg(int i1, int i2, int i3);
    int WriteMsh(char *name);

    int GetNumOfEdge(int ie, int iv1, int iv2);

    int is_the_same(Edge &e1, Edge &e2);
    void ModyfyEdge(int ie, int i1, int i2, CoordXYZ q,    int &ic, int &ie1, int &ie2);

};
