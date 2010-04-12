// Copyright (C) M.Iakobovski, 2010, v.1

#include "funcs.h"
#include "farfield.h"
#include "vertices.h"
#include "tgen.h"

int RefTriangles::meshgent(char * name, int n)
{
    CoordXYZ p[4];

    for(int i=0;i<4;i++)
        p[i].resize(3);

    p[0][0]=0;
    p[0][1]=0;
    p[0][2]=R1;

    p[1][0]=0;
    p[1][1]=R1*2.*sqrt(2.)/3.;
    p[1][2]=-R1/3.;

    p[2][0]=p[1][1]*sqrt(3.)/2.;
    p[2][1]=-p[1][1]/2.;
    p[2][2]=-R1/3.;

    p[3][0]=-p[1][1]*sqrt(3.)/2.;
    p[3][1]=-p[1][1]/2.;
    p[3][2]=-R1/3.;

    int i0=AddPoint(p[0]);
    int i1=AddPoint(p[1]);
    int i2=AddPoint(p[2]);
    int i3=AddPoint(p[3]);

    int ie01=AddEdge(i0,i1);
    int ie02=AddEdge(i0,i2);
    int ie03=AddEdge(i0,i3);
    int ie12=AddEdge(i1,i2);
    int ie13=AddEdge(i1,i3);
    int ie23=AddEdge(i2,i3);

    tgent(ie01,ie02,ie12,n);
    tgent(ie02,ie03,ie23,n);
    tgent(ie03,ie01,ie13,n);
    tgent(ie12,ie23,ie13,n);

    WriteMsh(name);

    return 0;
}

void RefTriangles::ModyfyEdge(int ie, int i1, int i2, CoordXYZ q,    int &ic, int &ie1, int &ie2)
{
    ic=edges[ie].ic;
    if(ic>=0)
    {
        ie1=GetNumOfEdge(ie,i1,ic);
        ie2=GetNumOfEdge(ie,ic,i2);
    }
    else
    {
        ic=AddPoint(q);
        ie1=AddEdge(i1,ic);
        ie2=AddEdge(ic,i2);
        edges[ie].ic=ic;
        edges[ie].ie1=ie1;
        edges[ie].ie2=ie2;
    }
}

void RefTriangles::tgent(int ie1, int ie2, int ie3, int level) // измельчение одного треугольника
{
    int i1,i2,i3; //  номера вершин измельчаемого треугольника

    i1=edges[ie1].i1;
    i2=edges[ie1].i2;
// ie1: (i1 i2)

    i3=edges[ie2].i1;
    if(i3==i1 || i3==i2)
        i3=edges[ie2].i2;

    if(i2==edges[ie3].i1 || i2==edges[ie3].i2)
    {
        i1=edges[ie1].i2;
        i2=edges[ie1].i1;
    }
// ie2: (i2 i3)
// ie3: (i1 i3)


//  определены номера трёх точек - вершин измельчаемого треугольника

//cout << i1 << " " << i2 << " " << i3 << endl;

    if(!level)
    {
        AddiTrg(i1,i2,i3);
        return;
    }

    CoordXYZ q(3),q1(3),q2(3),q3(3);

    double  d;
    q=plus2(points[i1],points[i2]);  d=GetVectorNorma1(q); q1=smult2(R1/d,q);
    q=plus2(points[i2],points[i3]);  d=GetVectorNorma1(q); q2=smult2(R1/d,q);
    q=plus2(points[i3],points[i1]);  d=GetVectorNorma1(q); q3=smult2(R1/d,q);

    //  добавим, при необходимости, новые точки и рёбра

    int i12,i23,i31;
    int ie11,ie12,ie21,ie22,ie31,ie32,ie1_2,ie2_3,ie3_1;

    ModyfyEdge(ie1,i1,i2,q1,   i12, ie11,ie12);
    ModyfyEdge(ie2,i2,i3,q2,   i23, ie21,ie22);
    ModyfyEdge(ie3,i3,i1,q3,   i31, ie31,ie32);

    //  номера новых рёбер

    ie1_2=AddEdge(i31,i23);
    ie2_3=AddEdge(i12,i31);
    ie3_1=AddEdge(i23,i12);

    tgent(ie1_2,ie2_3,ie3_1,level-1);
    tgent(ie11,ie32,ie2_3,level-1);
    tgent(ie12,ie21,ie3_1,level-1);
    tgent(ie22,ie31,ie1_2,level-1);
}

int RefTriangles::GetNumOfEdge(int ie, int iv1, int iv2)
{
    Edge e(iv1,iv2);
    size_t k;

    if(is_the_same(e,edges[edges[ie].ie1]))
        k=edges[ie].ie1;
    else
        k=edges[ie].ie2;

    if(k<0 || k>edges.size())
        {
            cout << " ERROR 5 : GetNumOfEdge k=" << k << " ie=" << ie << " " << iv1 << " " << iv2 << endl;
            edges[ie].put();
            exit(1);
        }

    return k;
}

int RefTriangles::is_the_same(Edge &e1, Edge &e2)
{
    if(
     ((e1.i1==e2.i1) && (e1.i2==e2.i2))
     ||
     ((e1.i1==e2.i2) && (e1.i2==e2.i1))
      )
      return 1;

    return 0;
}

int RefTriangles::AddEdge(int i1, int i2) // возвращает номер ребра
{
    int n=edges.size(),k=0;
    Edge e(i1,i2);

//cout << n << " <> " << e.i1 << " " << e.i2 << endl;

    for(k=0;k<n;k++)
        if(is_the_same(e,edges[k]))
            return k;

    edges.push_back(e);
    return n;
}

int RefTriangles::AddPoint(CoordXYZ p) // возвращает номер точки
{
    points.push_back(p);
    return points.size()-1;
}

void RefTriangles::AddiTrg(int i1, int i2, int i3)
{
    iTrg t(3);
    t[0]=i1;
    t[1]=i2;
    t[2]=i3;
    itrgs.push_back(t);
}

int RefTriangles::WriteMsh(char *name)
{
    size_t i,kpyr=0;
    double d;

    int ntrg=itrgs.size();  // базовые треугольники
    int npts=points.size(); // базовые точки

    int npts_res=2*npts+ntrg;  // число получившихся точек

    FILE *fcoord,*fpyr;

    // Сформируем дополнительные точки
    CoordXYZ q(3);

    for(i=0;i<(size_t)npts;i++)
    {
        q=smult2(R2/R1,points[i]);
        AddPoint(q);
    }

    for(i=0;i<(size_t)ntrg;i++)
    {
        q=plus2(plus2(points[itrgs[i][0]],points[itrgs[i][1]]),points[itrgs[i][2]]);
        d=GetVectorNorma1(q);
        q=smult2(R2/d,q);
        AddPoint(q);
    }

// для сетки
{
    fcoord=fopen("Tmesh\\coordinate.msh","w");

    for(i=0;i<points.size();i++)
        fprintf(fcoord,"%g %g %g\n",points[i][0],points[i][1],points[i][2]);

    fclose(fcoord);

    fpyr=fopen("Tmesh\\tetrahedron.msh","w");

    for(i=0;i<itrgs.size();i++)
    {
        int i1=itrgs[i][0];
        int i2=itrgs[i][1];
        int i3=itrgs[i][2];
        int i7=i+2*npts;

        int j1,j2,j1up,j2up;

        i1++;i2++;i3++;i7++;

//-------------

        fprintf(fpyr,"%d %d %d %d\n",i1,i2,i3,i7); kpyr++;

//-------------
if(1)
{
        if(i1<i2)   { j1=i1; j2=i2; }
        else        { j1=i2; j2=i1; }

        j1up = j1+npts;  j2up = j2+npts;

        fprintf(fpyr,"%d %d %d %d\n",j1,j2,j2up,i7); kpyr++;
        fprintf(fpyr,"%d %d %d %d\n",j1,j2up,j1up,i7); kpyr++;
//-------------

        if(i2<i3)   { j1=i2; j2=i3; }
        else        { j1=i3; j2=i2; }

        j1up = j1+npts;  j2up = j2+npts;

        fprintf(fpyr,"%d %d %d %d\n",j1,j2,j2up,i7); kpyr++;
        fprintf(fpyr,"%d %d %d %d\n",j1,j2up,j1up,i7); kpyr++;
//-------------

        if(i3<i1)   { j1=i3; j2=i1; }
        else        { j1=i1; j2=i3; }

        j1up = j1+npts;  j2up = j2+npts;

        fprintf(fpyr,"%d %d %d %d\n",j1,j2,j2up,i7); kpyr++;
        fprintf(fpyr,"%d %d %d %d\n",j1,j2up,j1up,i7); kpyr++;
}
//-------------
    }

    fclose(fpyr);

    // вывод общей информации

    std::ofstream foutmesh("Tmesh\\mesh.msh");
    foutmesh << npts_res << endl;
    foutmesh << kpyr << endl;
    foutmesh << 0 << endl;
    foutmesh.close();
}

// для визуализации
if(0)
{
    FILE *fout=fopen("r.wff","w");

    fprintf(fout,"# points %d \n",npts_res);

    for(i=0;i<points.size();i++)
        fprintf(fout,"v %g %g %g\n",points[i][0],points[i][1],points[i][2]);

    for(i=0;i<itrgs.size();i++)
    {
        int i1=itrgs[i][0];
        int i2=itrgs[i][1];
        int i3=itrgs[i][2];
        int i7=i+2*npts;

        int j1,j2,j1up,j2up;

        i1++;i2++;i3++;i7++;

//-------------

//      fprintf(fpyr,"%d %d %d %d\n",i1,i2,i3,i7);
        fprintf(fout,"f %d %d %d\n",i1,i2,i3);

//-------------

        if(i1<i2)   { j1=i1; j2=i2; }
        else        { j1=i2; j2=i1; }

        j1up = j1+npts;  j2up = j2+npts;

//        fprintf(fpyr,"f %d %d %d %d\n",j1,j2,j2up,i7);
//        fprintf(fpyr,"f %d %d %d %d\n",j1,j2up,j1up,i7);

        fprintf(fout,"f %d %d %d\n",j1,j2,j2up);
        fprintf(fout,"f %d %d %d\n",j1,j2up,i7);

//-------------

        if(i2<i3)   { j1=i2; j2=i3; }
        else        { j1=i3; j2=i2; }

        j1up = j1+npts;  j2up = j2+npts;

//        fprintf(fpyr,"%d %d %d %d\n",j1,j2,j2up,i7);
//        fprintf(fpyr,"%d %d %d %d\n",j1,j2up,j1up,i7);

        fprintf(fout,"f %d %d %d\n",j1,j2,j2up);
        fprintf(fout,"f %d %d %d\n",j1,j2up,i7);
//-------------

        if(i3<i1)   { j1=i3; j2=i1; }
        else        { j1=i1; j2=i3; }

        j1up = j1+npts;  j2up = j2+npts;

//        fprintf(fpyr,"%d %d %d %d\n",j1,j2,j2up,i7);
//        fprintf(fpyr,"%d %d %d %d\n",j1,j2up,j1up,i7);

        fprintf(fout,"f %d %d %d\n",j1,j2,j2up);
        fprintf(fout,"f %d %d %d\n",j1,j2up,i7);
//-------------
    }
    fclose(fout);
}


    return 0;
}



