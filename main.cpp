// Copyright (C) M.Iakobovski, E.Kornilina, A.Morozov, 2010, v.1

//#pragma once

#define _main_parm_

#include "funcs.h"
#include "farfield.h"
#include "vertices.h"
#include "tgen.h"
#include "string.h"
#include "graph.h"

Copyright(ncpr,FarFieldK,2010,1);

#define meshdir1 "cubicsmall"   // IndexOfFirstPoint 1
#define meshdir2 "MESH_1800110" // IndexOfFirstPoint 1
#define meshdir3 "SETKA1"        // IndexOfFirstPoint 0
#define meshdir4 "mesh2"        // IndexOfFirstPoint 1
#define meshdir5 "Tmesh"        // IndexOfFirstPoint 1
#define IndexOfFirstPoint 1

//#define R_R 10.

int main(int argc, char *argv[])
{
    if(argc<6)
        {
        cout << ncpr << endl << endl;
        cout << "USAGE: FarFieldKgen <Mesh dir> DeltaT Xsph Ysph zSph Rsphere "
                "[ xContr yContr zContr [ T1 T2 [Omega]]]" << endl;
        exit(0);
        }

    strncpy(p_outdir,"./",10000);

    strncpy(p_meshdir,argv[1],10000);

    sscanf(argv[2],"%le",&p_DeltaT);
    sscanf(argv[3],"%le",&p_xc);
    sscanf(argv[4],"%le",&p_yc);
    sscanf(argv[5],"%le",&p_zc);
    sscanf(argv[6],"%le",&p_Rsphere);

    p_control_x=p_xc;
    p_control_y=p_yc;
    p_control_z=p_zc+90.5;

    if(argc>=10)
    {
    sscanf(argv[7],"%le",&p_control_x);
    sscanf(argv[8],"%le",&p_control_y);
    sscanf(argv[9],"%le",&p_control_z);
    }

    p_T1=90.;
    p_T2=110.;

    if(argc>=12)
    {
    sscanf(argv[10],"%le",&p_T1);
    sscanf(argv[11],"%le",&p_T2);
    }

    p_omega=M_PI;
    if(argc>=13)
        sscanf(argv[12],"%le",&p_omega);

    printf("%s DeltaT=%g Rc=(%g,%g,%g) Rs=%g omega=%g\n",
        p_meshdir,p_DeltaT,p_xc,p_yc,p_zc,p_Rsphere,p_omega);

    printf("Rcontrol=(%g,%g,%g) T1=%g T2=%g\n",
        p_control_x,p_control_y,p_control_z,p_T1,p_T2);

    e_a=p_Rsphere;
    e_b=p_Rsphere;
    e_c=p_Rsphere;

    e_x0=p_xc;
    e_y0=p_yc;
    e_z0=p_zc;

    if(!strstr(p_meshdir,"Tmesh"))
    {
        RefTriangles RT;
        RT.meshgent("tg.wff",3);
    }

//    test_of_basepoint(); exit(0);

	int numc, numt; //количество узлов и тетраэдров

    #define NN 10000
    char name[NN+1];

    strncpy(name,p_meshdir,NN);
    strncat(name,"/mesh.msh",NN);
    std::ifstream meshf(name);
	if(!meshf) { cout << "Unable to find file " << name << "\n"; exit(1); }
	meshf >> numc >> numt;
	meshf.close();

    cout << "\n numc=" << numc << " numt=" << numt << endl ;

    strncpy(name,p_meshdir,NN);
    strncat(name,"/coordinate.msh",NN);
	std::ifstream coords(name);
	if(!coords) { cout << "Unable to find file " << name << "\n"; exit(1); }

    int km=0;
	CoordXYZ tempmesh(kNumberOfCoord);
	PointsXYZ mesh(numc,tempmesh); //вектор со всеми узлами решетки
	for(int j = 0; j < numc; j++){
		for(int i = 0; i < kNumberOfCoord; ++i)
		{
		    coords >> mesh[j][i];
		}
        if(((++km)%100000)==0)        cout << "coord " << km << " get ok" << endl;
	}
	coords.close(); //закрываем файловый поток
    cout << "coord " << km << " get ok" << endl;

    strncpy(name,p_meshdir,NN);
    strncat(name,"/tetrahedron.msh",NN);
	std::ifstream points(name);
	if(!points) { cout << "Unable to find file " << name << "\n"; exit(1); } //проверяем наличие файла

    Surfaces SurfaceList;
    vector<Tetras> IntersectTetrahedron;

    int ko=0,kp=0;

    while(!points.eof()) //считываем тетраэдры до конца файла
    {
        if(ko>=numt)break;

        Tetras temptetra(kNumberOfTetVertexes);

        //читаем по kNumberOfProcTetra тетраэдров
        vector<Tetras> tetrahedron; //<=1000000 тетраэдров
        tetrahedron.reserve(kNumberOfProcTetra);

        int k = 0;
        while(!points.eof() && k < kNumberOfProcTetra)
        {
            for(int i =0; i < kNumberOfTetVertexes; ++i)
            {
                points >> temptetra[i];
                temptetra[i]-=IndexOfFirstPoint;

            }

            tetrahedron.push_back(temptetra);
            k++;

            if(((++ko)%100000)==0) cout << "tetrahedron " << ko << endl;

            if(ko>=numt)break;
        }

        cout << "tetrahedron " << ko << endl;
        tetra_to_wff(tetrahedron,mesh);

        for(int i = 0; i<k; ++i) //обрабатываем считанные тетраэдры
        {
 			SurfacePoints sur = FindSurface(tetrahedron[i], mesh);

			if(sur.size()==4)
			{
				Surfaces temp = DivideSquares(sur,mesh);

				if(GetTriangleParams(temp[0],mesh).area)
					SurfaceList.push_back(temp[0]);

				if(GetTriangleParams(temp[1],mesh).area)
					SurfaceList.push_back(temp[1]);
			}
			else if(sur.size()==3)
				if(GetTriangleParams(sur,mesh).area)
					SurfaceList.push_back(sur);

        if(((++kp)%100000)==0)        cout << "tetrahedron " << kp << " processed ok. SurfaceList: " << SurfaceList.size() << endl;
        }
	}

    cout << "tetrahedron " << ko  << endl;

	double RealArea = 0;
	for(size_t i =0; i < SurfaceList.size(); ++i){
        if(RoundControl(SurfaceList[i],mesh))
            {
               cout << " ERROR 8_4 i=" << i << " SurfaceList.size()=" << SurfaceList.size() <<"\n" ;
               exit(1);
            }
            RealArea+=GetTriangleParams(SurfaceList[i],mesh).area;
	}

    points.close();
	double area = 4*M_PI*kva(e_a);
	PutSurface(strcat(p_outdir,"_surf_trangles.txt"),SurfaceList,mesh);

//	CheckTrgSurface(SurfaceList,mesh);
//	CheckTrgSurface(SurfaceList,mesh);

//	PutSurface("surf_0.txt",SurfaceList,mesh);

	CreateWff(strcat(p_outdir,"_surf_trangles.wff"),SurfaceList,mesh);

	CalcIntegral(SurfaceList,mesh);
	cout << "Площадь интерполированной поверхности: " << RealArea << endl << "Площадь реальной поверхности: " << area << endl << "Отличаются на: " << (area - RealArea)/area*100 << "%" << endl; //сравниваем площадь поверхности с площадью построенной поверхности
	vector<SimplifiedRib> ribs = CreateInitialStructure(SurfaceList,mesh);
	if(CheckSide(ribs)) cout << "All ribs conjugate to two sides! OK" <<  endl;
	else cout << "Have some troubles with ribs! Check data!" << endl;
	if(CheckGraphCoherency(ribs,SurfaceList.size())) cout << "Graph is coherent" << endl;
	else cout << "Graph is not coherent" << endl;
	return 0;
}
