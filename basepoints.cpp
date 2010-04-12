// Copyright (C) M.Iakobovski, E.Kornilina, A.Morozov, 2010, v.1

#include "funcs.h"
#include "farfield.h"
#include "vertices.h"

BasePoints::BasePoints(CoordXYZ rpyr_[kbase], int ip_[4],int &err)
{
    err=0;

    for (int i = 0; i < kbase; i++)
    {
        rpyr[i]=rpyr_[i];
        ip[i]=ip_[i];
    }

    for(int i=0;i<kbase;i++)
    {
        // надо решить четыре системы уравнений для 4х ортов значения функции

        err=solv4eq_ort(i);

        if(err)
            break;
    }
}

double BasePoints::GetVal(CoordXYZ r, int ip_[4], CoordXYZ f_)
{
    double g=0.;

if(0)
    printf("( %d %d %d %d : %g %g %g %g : %g %g %g\n",
        ip_[0],ip_[1],ip_[2],ip_[3],
        f_[0]   ,f_[1]   ,f_[2]   ,f_[3],
        r[0]   ,r[1]   ,r[2]
//        r[0][0],r[0][1],r[0][2]
        );
    for (int i1 = 0; i1 < kbase; i1++)
    {
        int i;
        for(i=0;i<kbase && ip[i]!=ip_[i1];i++);

        if(i==kbase)
        {
            cout << " ERROR 6 : unpresent vertex " << ip_[i1] << " in base points of pyr " << endl;
            exit(1);
        }

//        if(i!=i1) cout << endl << i1 <<" mda " << i << endl;

        double q=a_line_form[i][kbase-1];

        for (int j = 0; j < kbase-1; j++)
            q+=(a_line_form[i][j]*r[j]);

        g+=(q*f_[i1]);
    }

    return g;
}

//===============================================

int BasePoints::solv4eq_ort(int i_of_non_zero_ort)
{
    int i,j;

    CoordXYZ r[kbase];
    CoordXYZ b(kbase);

//    cout << "r_[i]" << endl;    for(i=0;i<kbase;i++)    {        putvect(r_[i]);        cout << endl;    }    cout << " " << endl;

    for(i=0;i<kbase;i++)
        b[i]=(i==i_of_non_zero_ort);

    for(i=0;i<kbase;i++)
        r[i].resize(4);

#if 0
    for(i=0;i<kbase;i++)
    {
        for(int j=0;j<3;j++)
            r[i][j]=rpyr[i][j];
    }

    for(i=0;i<kbase;i++)
        r[i][3]=1;
#else
    for(i=0;i<kbase-1;i++)
    {
        for(int j=0;j<kbase;j++)
            r[i][j]=rpyr[j][i];
    }

    for(i=0;i<kbase;i++)
        r[3][i]=1;
#endif

//    cout << "b" << endl;        putvect(b);        cout << endl;    cout << " " << endl;
//    cout << endl << "r[i]" << endl;    for(i=0;i<kbase;i++)    {        putvect(r[i]);        cout << endl;    }    cout << " " << endl;

    double D=det4(r[0],r[1],r[2],r[3]); // параметры - столбцы

    if(fabs(D)<1e-14)
        {
            // !!! эта ситуация означает, что пирамида имеет нулевой объём
            cout << "ERROR 3 small D= " << D << " - zero volume of base pyramid" << endl;


            for(i=0;i<kbase;i++)
            {
                for(j=0;j<kbase;j++)
                    cout << r[i][j] << "\t" ;
                cout << endl;
            }

            return 1;
        }

    a_line_form[i_of_non_zero_ort][0]=det4(   b,r[1],r[2],r[3])/D;
    a_line_form[i_of_non_zero_ort][1]=det4(r[0],   b,r[2],r[3])/D;
    a_line_form[i_of_non_zero_ort][2]=det4(r[0],r[1],   b,r[3])/D;
    a_line_form[i_of_non_zero_ort][3]=det4(r[0],r[1],r[2],   b)/D;

//    cout << endl << "D=" << D << "  a_line: i_of_non_zero_ort=" << i_of_non_zero_ort << endl;    for(i=0;i<4;i++)        cout << a_line[i] << " " ;    cout << endl;    cout << endl << "======================" << endl << endl;

    return 0;
}


//----------------------------

void test_of_basepoint(void)
{
    int i,j,ret;
    SourceOfNoise source_of_noise;
    CoordXYZ control_point;
    CoordXYZ rp[4];
    CoordXYZ r(3);
    int ip[4]={0,1,2,3};
    int ip1[4]={3,1,2,0};

    for(i=0;i<4;i++)
        rp[i].resize(3);

// (0,0,0)
    for(j=0;j<3;j++)
        rp[0][j]=-1;

        rp[0][2]=-7;

// (1,0,0)
// (0,2,0)
// (0,0,3)
    for(i=1;i<4;i++)
        for(j=0;j<3;j++)
            rp[i][j]=(i==(j+1));

    BasePoints bp(rp,ip,ret);

    cout << "ret=" << ret << endl;

    CoordXYZ f(4);
    f[0]=4;
    f[1]=2;
    f[2]=3;
    f[3]=1;

    for(i=0;i<4;i++)
    {
        putvect(rp[i]); cout << " :=: " << bp.GetVal(rp[i],ip,f) << endl;
    }
    cout << endl;

    r[0]=.0;
    r[1]=.5;
    r[2]=.5;

    putvect(r); cout << "  =  " << bp.GetVal(r,ip1,f) << endl;

/*
    cout << endl << "f( " ;
    putvect(r);
    cout << " : " ;
    putvect(f);
    cout << " ) = " << bp.GetVal(r,f) << endl;
*/
}

//-----------------------------------

double det3(CoordXYZ &x, CoordXYZ &y, CoordXYZ &z){
	return x[0]*y[1]*z[2]-x[0]*y[2]*z[1]-x[1]*y[0]*z[2]+x[1]*y[2]*z[0]+x[2]*y[0]*z[1]-x[2]*y[1]*z[0];
}

double det4(CoordXYZ &r0,CoordXYZ &r1,CoordXYZ &r2,CoordXYZ &r3)
{
   return
    -r0[3]*det3(r1,r2,r3)
    +r1[3]*det3(r0,r2,r3)
    -r2[3]*det3(r0,r1,r3)
    +r3[3]*det3(r0,r1,r2);
}

