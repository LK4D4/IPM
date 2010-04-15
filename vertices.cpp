// Copyright (C) M.Iakobovski, E.Kornilina, A.Morozov, 2010, v.1

#include "funcs.h"
#include "farfield.h"
#include "vertices.h"

int Vertices::Init(int n_, double timestep_)
{
//    return 0;
     n=n_;
     timestep=timestep_;

     vert=new Vertex * [n];

     for(int i=0;i<n;i++)
        vert[i]=NULL;

    return 0;
}

double Vertices::Svertka(Surfaces &s, PointsXYZ &m, SourceOfNoise &sn, double time)
{
    double sum=0.;

    for(int iv=0;iv<n;iv++)
      if(vert[iv])
      {
        for(int i=0;i<nOfsTimeStrips && vert[iv][i].itime;i++)
           {
                double t=time+vert[iv][i].itime*timestep;

                CoordXYZ u_der = sn.GetUDerivative(m[iv], t);
                double rho_der = sn.GetRhoDerivative(m[iv], t);
                double rho = sn.GetRho(m[iv], t);

                sum+=(rho*vert[iv][i].a[0]);
                sum+=(rho_der*vert[iv][i].a[1]);

                for(int j=0;j<3;j++)
                        sum+=(u_der[j]*vert[iv][i].a[j+2]);
           }
      }

    return sum;
}

int Vertices::AddData(int iv, int i_func_num, double value_of_func, double time)
{
    int i;

    if(incorrect)
        return 2;

//if(!i_func_num)  cout << "AddData " << iv << " " << time << " : ";

    if(!vert[iv])
    {
        vert[iv]=new Vertex[nOfsTimeStrips];
        for(i=0;i<nOfsTimeStrips;i++)
        {
            vert[iv][i].itime=0;

            for(int j=0;j<nOfFuncs;j++)
                vert[iv][i].a[j]=0;
        }
    }

    int it=StripOfTime(time);

    if(imintime>imaxtime)
        imintime=imaxtime=it;

    if(imintime>it)imintime=it;
    if(imaxtime<it)imaxtime=it;

    for(i=0;(i<nOfsTimeStrips) &&
            vert[iv][i].itime &&
            (it!=vert[iv][i].itime);
            i++);

//    if(!i_func_num)       cout << i << " " << it << endl;

    if(i==nOfsTimeStrips)
    {
        PutData("err.txt");
        cout << endl << "ERROR 4 : Vertices::AddData needed more then " << nOfsTimeStrips <<
            " Strips of time for vertex " << iv << endl;

        incorrect=1;
        return 1;
        exit(0);
    }


    vert[iv][i].itime=it;

    vert[iv][i].a[i_func_num]+=value_of_func;

    if(max_strips_per_vertex<i) max_strips_per_vertex=i;

    return 0;
}

int Vertices::StripOfTime(double time)
{
//    return 0;
    return (int)(time/timestep)-1;
}

int Vertices::PutData(char *name)
{
//    return 0;
	std::ofstream fout(name);
    int it;

    if(incorrect)
    {
        cout << "Data incorrect " << endl;
        return 2;
    }

    cout << name << endl;

    fout << timestep << endl;
    fout << imintime << endl;
    fout << imaxtime << endl;
    fout << max_strips_per_vertex << endl;
    fout << nOfFuncs << endl;

    for(it=imintime;it<=imaxtime;it++)
    {
        int count=0;

        for(int q=0;q<2;q++)
        {
            if(q)
            {
                if(count)
                    fout << it << " " << count << endl;
                else
                    break;
            }

            for(int iv=0;iv<n;iv++)
            {
                if(vert[iv])
                {
                for(int i=0;i<nOfsTimeStrips && vert[iv][i].itime;i++)
                    {
                        if(vert[iv][i].itime==it)
                        {
                            if(!q)
                                count++;
                            else
                            {
                                fout << "  " << iv << endl;
                                for(int k=0;k<nOfFuncs;k++)
                                    fout << "       " << vert[iv][i].a[k] << endl ;
//                                    fout << vert[iv][i].a[k] << " " ;
  //                              fout << endl;
                            }
                        }
                    }
                }
            }
        }
    }

    fout << endl;

    cout << "imintime= " << imintime << " imaxtime= " << imaxtime << " max_strips_per_vertex= " << max_strips_per_vertex <<endl;
    return 0;
}

