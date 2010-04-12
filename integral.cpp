// Copyright (C) M.Iakobovski, E.Kornilina, A.Morozov, 2010, v.1

#include "funcs.h"
#include "farfield.h"
#include "vertices.h"

CoordXYZ IntN(3);
double S_Surf;
//extern std::ofstream fout1;
//extern int kp;

int const isbp=1;

Vertices vrt;

double CalcIntegral(Surfaces &s, PointsXYZ &m){
//return;
	int ns=s.size(); // число треугольников
    CoordXYZ x[3];

    for (int i = 0; i < 3; i++)
        x[i].resize(3);

//============= init ===================

    vrt.Init(m.size(),p_DeltaT);

	CoordXYZ w(3);
	CoordXYZ control_point(3);
	CoordXYZ position_of_source(3);

	w[0] = 0;
	w[1] = 0;
	w[2] = 0;

	position_of_source[0] = e_x0;
	position_of_source[1] = e_y0;
	position_of_source[2] = e_z0;

	control_point[0] = p_control_x;
	control_point[1] = p_control_y;
	control_point[2] = p_control_z;

	double omega = p_omega;

	SourceOfNoise source_of_noise(omega, position_of_source,w);//,control_point);

//========== end of init ===============


    double Sintegral;

//  for(double z=90.;z<100.;z+=(100.-90.)/200.)
    double T=0.;
//    for(;T<110.;T+=(110.-90.)/230.)
    {
    for (int i = 0; i < 3; ++i)
        IntN[i]=0;
        S_Surf=0.;

//     control_point[2]=z;

       Sintegral=0;

       for(int j=0;j<ns;j++)
       {
            int err;

            Sintegral+=CalcIntegralOfFullTrg(s[j],m,T,source_of_noise, control_point,err);

            if(err)
            {
                cout << " ERROR 3_ ret=" << err << " j= " << j << endl;
                return 0.;
            }
       }

        {
         static int ugu=0;
         if(!ugu)   cout << " S_Surf=" << S_Surf << " (" <<  S_Surf/(4.*M_PI*kva(p_Rsphere)) << ") IntN " << IntN[0] << " " << IntN[1] << " " << IntN[2] << endl<< endl;
         ugu=1;
        }

//        cout << control_point[2] << " " << Sintegral << " " << source_of_noise.GetRho(control_point,T) << endl;
       if(!isbp)
          cout << T << "\t" << source_of_noise.GetRho(control_point,T) << " " << Sintegral << endl;

    }

       if(isbp)
       {
         vrt.PutData("koeff.txt");

         // следует выполнить свёртку!

         int nk=150,ik;
         double dt=(p_T2-p_T1)/nk;

         for(ik=0,T=p_T1;ik<=nk;ik++,T+=dt)
            {
            Sintegral=vrt.Svertka(s,m,source_of_noise,T);
            cout << T  << "\t" << source_of_noise.GetRho(control_point,T) << " " << Sintegral << endl;
            }
       }

    cout << endl;
    cout  << " S_Surf=" << S_Surf << " (" <<  S_Surf/(4.*M_PI*kva(p_Rsphere)) << ") Int_N " << IntN[0] << " " << IntN[1] << " " << IntN[2] << endl;

    return Sintegral;
}

int intime(double t){
//    return 1;
    if(t>=0 && t<= 1. )
        return 1;
    return 0;
}

int add_to_ip(int *ip, int &k, int ind);
int add_to_ip(int *ip, int &k, int ind)
{
    if(k>=4) return 1;

    for(int i=0;i<k;i++)
        if(ip[i]==ind)
            return 0;

    ip[k++]=ind;

    return 0;
}

int ip[4];

double CalcIntegralOfFullTrg(SurfacePoints p, PointsXYZ &m, double time, SourceOfNoise &source_of_noise, CoordXYZ &control_point, int &err)
{
    int i,j,n=4;
    double rho=0;
    CoordXYZ strg[3],h1,h2;
    CoordXYZ xtrg[3];

    if(p.size()!=3) { cout << " ERROR 8_1 " << p.size() << "!=3\n"; exit(0);}

    for (int i = 0; i < 3; i++)
        xtrg[i] = GetCoordSurf(p[i],m);

//    cout << " CalcIntegralOfFullTrg started " << endl;

    CoordXYZ pyr_points[4];
    for (int i = 0; i < kbase; i++)
        pyr_points[i].resize(3);

    // сформируем пирамиду
    int k=0;

    for (int i = 0; i < 3; i++)
    {
        add_to_ip(ip,k,p[i].a_index);
        add_to_ip(ip,k,p[i].b_index);
    }

    if(k!=4)
    {
        cout << "ERROR 6 : In base tetrahtdron only " << k << " nodes <4.";
        err=2;
        return 0;
    }

    //   cout << ip[0] << " " << ip[1] << " " << ip[2] << " " << ip[3] << " " << endl;
    for (int i = 0; i < 4; i++)
       pyr_points[i]=m[ip[i]];

    BasePoints bp(pyr_points,ip,err);

    if(err)
        return err;

    for (int i = 0; i < 3; i++)
        strg[i].resize(3);

/*
 к этому моменту есть :
    - координаты вершин базовой пирамиды
    - номера вершин базовой пирамиды
    - функция вычисления значения поля в точке r
      при заданных значениях f в узлах пирамиды
*/

    for(i=0;i<n;i++)
    for(j=0;j<n-i;j++)
    {
        strg[0]=
        plus2(
            plus2(xtrg[0],
                    smult2((double)i/(double)n,minus2(xtrg[1],xtrg[0]))),
                        smult2((double)j/(double)n,minus2(xtrg[2],xtrg[0])));
        strg[1]=
        plus2(
            plus2(xtrg[0],
                    smult2((double)(i+1)/(double)n,minus2(xtrg[1],xtrg[0]))),
                        smult2((double)j/(double)n,minus2(xtrg[2],xtrg[0])));
        strg[2]=
        plus2(
            plus2(xtrg[0],
                    smult2((double)i/(double)n,minus2(xtrg[1],xtrg[0]))),
                        smult2((double)(j+1)/(double)n,minus2(xtrg[2],xtrg[0])));

        if(!isbp) rho+=CalcIntegralOfSmallTrg(strg, time, source_of_noise, control_point);
        else      rho+=CalcIntegralOfSmallTrgBp(strg, time, bp, control_point);

        if(j<n-i-1)
        {
            strg[0]=
            plus2(
                plus2(xtrg[0],
                        smult2((double)(i+1)/(double)n,minus2(xtrg[1],xtrg[0]))),
                            smult2((double)(j+1)/(double)n,minus2(xtrg[2],xtrg[0])));
            strg[2]=
            plus2(
                plus2(xtrg[0],
                        smult2((double)(i+1)/(double)n,minus2(xtrg[1],xtrg[0]))),
                            smult2((double)j/(double)n,minus2(xtrg[2],xtrg[0])));
            strg[1]=
            plus2(
                plus2(xtrg[0],
                        smult2((double)i/(double)n,minus2(xtrg[1],xtrg[0]))),
                            smult2((double)(j+1)/(double)n,minus2(xtrg[2],xtrg[0])));

        if(!isbp) rho+=CalcIntegralOfSmallTrg(strg, time, source_of_noise, control_point);
        else      rho+=CalcIntegralOfSmallTrgBp(strg, time, bp, control_point);
        }
    }


    return rho;
}


double CalcIntegralOfSmallTrgBp(CoordXYZ xtrg[3], double time, BasePoints &bp, CoordXYZ &control_point)
{
	double cell_square;
	CoordXYZ cell_normal(3), cell_mass_centre(3);

    if (GetCellParams(xtrg, cell_square, cell_normal,cell_mass_centre)){
             return  0;
    }

    CoordXYZ R(3);

    for (int i = 0; i < 3; ++i)
        R[i] = control_point[i] - cell_mass_centre[i];

    double norma_R = GetVectorNorma1(R);
    if (!norma_R){
        cout << " Error 1 : norma_r=0" << endl;
        exit(1);
    }

 	double SUM = 0;

    double t = ModifyTimeByPosition(R, time);

    CoordXYZ f(kbase);
 	for(int i=0;i<kbase;i++)
        f[i]=0;

// 	double rez[kbase];

    // найдём отклик от каждой вершины
    // этот фрагмент следует вычислять для каждого узла пирамиды
    // и для каждой функции отдельно

    double rho=0, rho_der=0, u_der[3];

    for(int nonZeroNode=0; nonZeroNode<4; nonZeroNode++) // по каждой вершине
    {
        f[nonZeroNode]=1;
        double rez=bp.GetVal(cell_mass_centre,ip,f); // величина отклика от вершины, без учета функциональной зависимости
        f[nonZeroNode]=0;

        for(int nonZeroFunc=0; nonZeroFunc<5; nonZeroFunc++)
        {
            /*
            Эти функции равны 0 или 1 в зависимости от nonZeroFunc
            CoordXYZ u_der = bp.GetUDerivative(cell_mass_centre, t);
            double rho_der = bp.GetRhoDerivative(cell_mass_centre, t);
            double rho = bp.GetRho(cell_mass_centre, t);
            */

            switch(nonZeroFunc)
            {
                case 0: rho=1; rho_der=0; u_der[0]=0; u_der[1]=0; u_der[2]=0; break;
                case 1: rho=0; rho_der=1; u_der[0]=0; u_der[1]=0; u_der[2]=0; break;
                case 2: rho=0; rho_der=0; u_der[0]=1; u_der[1]=0; u_der[2]=0; break;
                case 3: rho=0; rho_der=0; u_der[0]=0; u_der[1]=1; u_der[2]=0; break;
                case 4: rho=0; rho_der=0; u_der[0]=0; u_der[1]=0; u_der[2]=1; break;
            }

            SUM=0.;

            for (int i = 0; i < 3; ++i)
            {
                SUM +=
                    1./(4*PI)*
                    (
                    u_der[i]/norma_R
                    +
                    R[i]/kva(norma_R)*(rho_der + rho/norma_R)
                    )
                    *cell_normal[i]
                    *cell_square;
            }

            int err=vrt.AddData(ip[nonZeroNode],nonZeroFunc,SUM*rez,t);

            if(err)
                exit(0);
        }
    }

    // Расчет инвариантов
	S_Surf+=cell_square;
	for (int i = 0; i < 3; ++i)
        IntN[i]+=(cell_square*cell_normal[i]);

    return SUM;
}

double CalcIntegralOfSmallTrg(CoordXYZ xtrg[3], double time, SourceOfNoise &source_of_noise, CoordXYZ &control_point)
{
// return 1;
//    PutTrgToWff(xtrg);

	double cell_square;
	CoordXYZ cell_normal(3), cell_mass_centre(3);

    if (GetCellParams(xtrg, cell_square, cell_normal,cell_mass_centre)){
 //           cout << " Error # 01 : norma_r=0" << endl;
            return  0;
    }

    CoordXYZ R(3);

    for (int i = 0; i < 3; ++i)
        R[i] = control_point[i] - cell_mass_centre[i];

    double norma_R = GetVectorNorma1(R);
    if (!norma_R){
        cout << " Error # 01 : norma_r=0" << endl;
        exit(1);
    }

    double t = ModifyTimeByPosition(R, time);

	CoordXYZ u_der = source_of_noise.GetUDerivative(cell_mass_centre, t);
	double rho_der = source_of_noise.GetRhoDerivative(cell_mass_centre, t);
	double rho = source_of_noise.GetRho(cell_mass_centre, t);

//    cout << norma_R << " " << rho << " " << rho_der << " " << u_der[0] << " ! " << u_der[1] << " " << u_der[2] << endl;
//    cout << cell_square << " ! " << cell_normal[0] << " " << cell_normal[1] << " " << cell_normal[2] << endl;

 	double SUM = 0;
	for (int i = 0; i < 3; ++i)
	{
		SUM +=1./(4*PI)*
            (
            u_der[i]/norma_R
            +
            R[i]/kva(norma_R)*(rho_der + rho/norma_R)
            )
            *cell_normal[i]
            *cell_square;

        IntN[i]+=(cell_square*cell_normal[i]);
    }

//cout << "!! IntN[" << i << "]=" <<  IntN[i] << "\t" << cell_square << "\t" << cell_normal[i] << endl;

	S_Surf+=cell_square;

//        cout << " cs " << cell_square*cell_normal[0] << " " << cell_square*cell_normal[1] << " " << cell_square*cell_normal[2] << endl;

    return SUM;
}



//=============================================
#define ModifyPositionByFlow(x) (x)

CoordXYZ SourceOfNoise::GetUDerivative(CoordXYZ &x0, double t){
	CoordXYZ x(3);
	for (int i = 0; i < 3; ++i)
		x[i] = x0[i] - position_of_source[i];
	x = ModifyPositionByFlow(x);
	CoordXYZ null(3);
	for (int i = 0; i < 3; ++i)
		null[i] = 0;
	double norma_x = GetVectorNorma1(x);
	if (!norma_x)
		return null;

    if(!intime(t-norma_x))
    		return null;

	CoordXYZ u_der(3);
	for (int i = 0; i < 3; ++i)
		u_der[i] = x[i] / kva(norma_x) * (  omega*cos(omega*(t-norma_x))
						+ 1./norma_x * sin(omega*(t-norma_x)) );


//	for (int i = 0; i < 3; ++i)		u_der[i] = 0.;

	return u_der;
}

double SourceOfNoise::GetRhoDerivative(CoordXYZ &x0, double t){
	CoordXYZ x(3),x1(3);
	for (int i = 0; i < 3; ++i)
		x1[i] = x0[i] - position_of_source[i];
	x = ModifyPositionByFlow(x1);
	double norma_x = GetVectorNorma1(x);
	if (!norma_x)
        {
        cout << " ERROR: GetRhoDerivative : small norma_x=" << norma_x << endl;

//        for (int i = 0; i < 3; ++i) cout <<

        exit(0);
//		return numeric_limits<double>::max();
        }

    if(!intime(t-norma_x))
    		return 0;

	return omega/norma_x * cos(omega*(t-norma_x));
}

double SourceOfNoise::GetRho(CoordXYZ &x0, double t){
	CoordXYZ x(3);
	for (int i = 0; i < 3; ++i)
		x[i] = x0[i] - position_of_source[i];
	x = ModifyPositionByFlow(x);
	double norma_x = GetVectorNorma1(x);

	if (!norma_x)
        {
        cout << " ERROR: GetRho : small norma_x=" << norma_x << endl;
        exit(0);
//		return numeric_limits<double>::max();
        }

    if(!intime(t-norma_x))		return 0;

	return 1./norma_x * sin(omega*(t-norma_x));
}



//=============================================

CoordXYZ GetVectorProduct(CoordXYZ &a,CoordXYZ &b)
{
    CoordXYZ r(3);

    r[0]= a[1]*b[2]-a[2]*b[1];
    r[1]=-a[0]*b[2]+a[2]*b[0];
    r[2]= a[0]*b[1]-a[1]*b[0];

    return r;
}

double GetVectorNorma1(CoordXYZ x){
	double norma = 0.;
	for (int i = 0; i < 3; ++i){
		norma += kva(x[i]);
	}
	return sqrt(norma);
}

double ModifyTimeByPosition(CoordXYZ x,double t){
	double norma_x = GetVectorNorma1(x);
#define kSpeedOfSound 1
	return t - norma_x/kSpeedOfSound;
}

int GetCellParams(  CoordXYZ r[3], // массив координат вершин треугольника
                    double &cell_square,
                    CoordXYZ &cell_normal,
					CoordXYZ &cell_mass_centre)
{
	CoordXYZ v1(3),v2(3);

    for(int i=0;i<3;i++)
    {
        cell_mass_centre[i] = (r[0][i]+r[1][i]+r[2][i])/3.;
        v1[i]=r[1][i]-r[0][i];
        v2[i]=r[2][i]-r[0][i];
    }

    cell_normal = GetVectorProduct(v1,v2);
    cell_square = GetVectorNorma1(cell_normal);

    if (cell_square < 1e-10)
        return 1;

    for(int i=0;i<3;i++)
        cell_normal[i]/=cell_square;

    cell_square/=2.;

	return 0;
}

CoordXYZ minus2(CoordXYZ a, CoordXYZ b)
{
    CoordXYZ c(3);
    for(int i=0;i<3;i++)
        c[i]=a[i]-b[i];

    return c;
}

CoordXYZ plus2(CoordXYZ a, CoordXYZ b)
{
    CoordXYZ c(3);
    for(int i=0;i<3;i++)
        c[i]=a[i]+b[i];

    return c;
}

CoordXYZ smult2(double a, CoordXYZ b)
{
    CoordXYZ c(3);
    for(int i=0;i<3;i++)
        c[i]=a*b[i];

    return c;
}

void putvect(CoordXYZ v)
{
    for(size_t i=0;i<v.size();i++)
        cout << " " << v[i];
}


