#include "funcs.h"
#include "farfield.h"

bool ellipse(Coord a){ //функция эллипсоида
    if((kva((a[0]-e_x0)/(e_a)) + kva((a[1]-e_y0)/(e_b)) + kva((a[2]-e_z0)/(e_c)))<1) return true;
    return false;
}

bool parallel(Coord a){ //функция параллелепипеда
		if((a[0] < p_x1 && a[0] > p_x0) && (a[1] < p_y1 && a[1] > p_y0) && (a[1] < p_z1 && a[1] > p_z0)) return true;
		return false;
}

bool cone(Coord a){ //функция конуса
	if(a[0] <= k_x0 || a[0] >= k_x1) return false;
	double temp_a = (k_a1-k_a0)*(a[0] - k_x0)/(k_x1 - k_x0) + k_a0; double temp_b = (k_b1-k_b0)*(a[0] - k_x0)/(k_x1 - k_x0) + k_b0;
	if(kva(a[1]-k_y0/temp_a) + kva(a[2]-k_z0/temp_b) < 1) return true;
	return false;
}

double mydistance(Coord a, Coord b){ //расстояние между точками
		return sqrt(kva(a[0]-b[0])+kva(a[1]-b[1])+kva(a[2]-b[2]));
}

bool isIn(Coord a){ //определяет внутри или снаружи поверхности точка
	switch (SurfaceType){
		case 0: return ellipse(a);//внутри эллипсоида
		case 1: return parallel(a); //внутри параллелепипеда
		case 2: return cone(a); //внутри конуса
	}
	return false;
}

bool isIntersect(Coord n, Coord m){ //функция возвращает true, если ребро пересекает поверхность и false, если нет
	if(isIn(m)!= isIn(n)) return true;
	return false;
}

Coord intersectSearch(Coord a, Coord b){ //ищем пересечение ребра с поверхностью
	Coord k;
	if(mydistance(a,b) < 0.000001) return a;
	Coord c(3);
	c[0] = (a[0]+b[0])/2; c[1] = (a[1]+b[1])/2; c[2] = (a[2]+b[2])/2;
	if(isIn(c) != isIn(a)) k = intersectSearch(a,c);
	else k = intersectSearch(c,b);
	return k;
}

double determinant(Coord x, Coord y, Coord z){
	return x[0]*y[1]*z[2]-x[0]*y[2]*z[1]-x[1]*y[0]*z[2]+x[1]*y[2]*z[0]+x[2]*y[0]*z[1]-x[2]*y[1]*z[0];
}

/*поиск поверхности, возвращает вектор длины 0, если тетраэдр не пересекает поверхность,
длины 3, если поверхность треугольник, 4, если четырехугольник. В векторе структуры типа SurfaceCoord.
delta вычисляется от меньшей по номеру вершины к большей.*/
SurfacePoints FindSurface(Tetras a, Points& mesh){
	SurfaceCoord temp;
	vector<Coord> tempc; //координаты вершин поверхности
	SurfacePoints ret; //поверхность
//проверяем разные ребра тетраэдра на пересечения
	if(isIntersect(mesh[a[0]],mesh[a[1]])){
		if(a[0] < a[1]){
			Coord p = intersectSearch(mesh[a[0]],mesh[a[1]]);
			temp.a_index=a[0]; temp.b_index=a[1];
			temp.delta = (mydistance(mesh[a[0]],p)/mydistance(mesh[a[0]],mesh[a[1]]));
			ret.push_back(temp);

		}
		else {
			Coord p = intersectSearch(mesh[a[1]],mesh[a[0]]);
			temp.a_index=a[1]; temp.b_index=a[0];
			temp.delta = (mydistance(mesh[a[1]],p)/mydistance(mesh[a[0]],mesh[a[1]]));
			ret.push_back(temp);
		}
	}

	if(isIntersect(mesh[a[0]],mesh[a[2]])){
		if(a[0] < a[2]){
			Coord p = intersectSearch(mesh[a[0]],mesh[a[2]]);
			temp.a_index=a[0]; temp.b_index=a[2];
			temp.delta = (mydistance(mesh[a[0]],p)/mydistance(mesh[a[0]],mesh[a[2]]));
			ret.push_back(temp);
		}
		else {
			Coord p = intersectSearch(mesh[a[2]],mesh[a[0]]);
			temp.a_index=a[2]; temp.b_index=a[0];
			temp.delta = (mydistance(mesh[a[2]],p)/mydistance(mesh[a[0]],mesh[a[2]]));
			ret.push_back(temp);
		}
	}

	if(isIntersect(mesh[a[0]],mesh[a[3]])){
		if(a[0] < a[3]){
			Coord p = intersectSearch(mesh[a[0]],mesh[a[3]]);
			temp.a_index=a[0]; temp.b_index=a[3];
			temp.delta = (mydistance(mesh[a[0]],p)/mydistance(mesh[a[0]],mesh[a[3]]));
			ret.push_back(temp);
		}
		else {
			Coord p = intersectSearch(mesh[a[3]],mesh[a[0]]);
			temp.a_index=a[3]; temp.b_index=a[0];
			temp.delta = (mydistance(mesh[a[3]],p)/mydistance(mesh[a[0]],mesh[a[3]]));
			ret.push_back(temp);
		}
	}

	if(isIntersect(mesh[a[1]],mesh[a[2]])){
		if(a[1] < a[2]){
			Coord p = intersectSearch(mesh[a[1]],mesh[a[2]]);
			temp.a_index=a[1]; temp.b_index=a[2];
			temp.delta = (mydistance(mesh[a[1]],p)/mydistance(mesh[a[1]],mesh[a[2]]));
			ret.push_back(temp);
		}
		else {
			Coord p = intersectSearch(mesh[a[2]],mesh[a[1]]);
			temp.a_index=a[2]; temp.b_index=a[1];
			temp.delta = (mydistance(mesh[a[2]],p)/mydistance(mesh[a[1]],mesh[a[2]]));
			ret.push_back(temp);
		}
	}

	if(isIntersect(mesh[a[1]],mesh[a[3]])){
		if(a[1] < a[3]){
			Coord p = intersectSearch(mesh[a[1]],mesh[a[3]]);
			temp.a_index=a[1]; temp.b_index=a[3];
			temp.delta = (mydistance(mesh[a[1]],p)/mydistance(mesh[a[1]],mesh[a[3]]));
			ret.push_back(temp);
		}
		else {
			Coord p = intersectSearch(mesh[a[3]],mesh[a[1]]);
			temp.a_index=a[3]; temp.b_index=a[1];
			temp.delta = (mydistance(mesh[a[3]],p)/mydistance(mesh[a[1]],mesh[a[3]]));
			ret.push_back(temp);
		}
	}

	if(isIntersect(mesh[a[2]],mesh[a[3]])){
		if(a[2] < a[3]){
			Coord p = intersectSearch(mesh[a[2]],mesh[a[3]]);
			temp.a_index=a[2]; temp.b_index=a[3];
			temp.delta = (mydistance(mesh[a[2]],p)/mydistance(mesh[a[2]],mesh[a[3]]));
			ret.push_back(temp);
		}
		else {
			Coord p = intersectSearch(mesh[a[3]],mesh[a[2]]);
			temp.a_index=a[3]; temp.b_index=a[2];
			temp.delta = (mydistance(mesh[a[3]],p)/mydistance(mesh[a[2]],mesh[a[3]]));
			ret.push_back(temp);
		}
	}
	return ret;
}

double V_pyr_of_sph(Coord p, Coord q, Coord g){ 
	CoordXYZ a(3),b(3),c(3); 

	a[0]=p[0]-e_x0; a[1]=p[1]-e_y0; a[2]=p[2]-e_z0; 
	b[0]=q[0]-e_x0; b[1]=q[1]-e_y0; b[2]=q[2]-e_z0; 
	c[0]=g[0]-e_x0; c[1]=g[1]-e_y0; c[2]=g[2]-e_z0; 

	return 
	( 
	(a[0])*(b[1])*(c[2])- 
	(a[0])*(b[2])*(c[1])- 
	(a[1])*(b[0])*(c[2])+ 
	(a[1])*(b[2])*(c[0])+ 
	(a[2])*(b[0])*(c[1])- 
	(a[2])*(b[1])*(c[0]) 
	) 
	/6.;
} 

void CheckTrgSurface(Surfaces &SurfaceList, Points &mesh){ 
	int err=0; 
	int ns=SurfaceList.size(); 
	Coord x[3]; 

	for (int i = 0; i < 3; i++) 
		x[i].resize(3); 

	for(int j=0;j<ns;j++){
		for (int i = 0; i < SurfaceList[j].size(); i++){
			int ip1=SurfaceList[j][i].a_index; 
			int ip2=SurfaceList[j][i].b_index; 
			double delta=SurfaceList[j][i].delta; 

			for(int l=0;l<3;l++) 
			x[i][l]=mesh[ip1][l]*(1-delta)+mesh[ip2][l]*delta; 
		} 

		double v1=V_pyr_of_sph(x[0],x[1],x[2]); 
		double v2=V_pyr_of_sph(x[0],x[2],x[1]); 

		if(v1<0) { 
			cout << j << " Error v1 = " << v1 << " < 0 ; v2 = " << v2 << endl;
			for(int k = 0; k < 3; ++k){
				for(int z = 0; z < 3; ++z){
					cout << x[k][z] << " ";
				}
				cout << endl;
			}
			err++; 

			int t; 
			double d; 
			t=SurfaceList[j][2].a_index; SurfaceList[j][2].a_index=SurfaceList[j][1].a_index; SurfaceList[j][1].a_index=t; 
			t=SurfaceList[j][2].b_index; SurfaceList[j][2].b_index=SurfaceList[j][1].b_index; SurfaceList[j][1].b_index=t; 
			d=SurfaceList[j][2].delta; SurfaceList[j][2].delta=SurfaceList[j][1].delta; SurfaceList[j][1].delta=d; 
		} 
	} 

	cout << "Check: number of incorrect faces = " << err << endl ; 
} 

Coord GetCoordSurf(SurfaceCoord c, Points &m){ //выдает координаты точек, заданных структурой SurfaceCoord
	Coord t(3);
	
	double lambda = c.delta;
	t[0] = (m[c.a_index][0]*(1-lambda)+lambda*m[c.b_index][0]);
	t[1] = (m[c.a_index][1]*(1-lambda)+lambda*m[c.b_index][1]);
	t[2] = (m[c.a_index][2]*(1-lambda)+lambda*m[c.b_index][2]);
	
	return t;
}
double GetTriangleSquare(SurfacePoints &s, Points &m){ //считает площадь поверхности, заданной вектором треугольников и четырехугольников
	// Heron's formula
	Points triangle(3);
	triangle[0] = GetCoordSurf(s[0],m);
	triangle[1] = GetCoordSurf(s[1],m);
	triangle[2] = GetCoordSurf(s[2],m);
	double a, b, c, p;
	a = mydistance(triangle[0],triangle[1]);
	b = mydistance(triangle[1],triangle[2]);
	c = mydistance(triangle[0],triangle[2]);
	p = (a + b + c) / 2;
	return sqrt(p*(p-a)*(p-b)*(p-c));
}

void CreateWff(Surfaces &s, Points &m){

	int ns=s.size();
	cout << "SurfaceList.size= " << ns << endl;
	std::ofstream fout("res_w0.wff");

	int k=1,l;
	double V=0;

	for(int j=0;j<ns;j++)
	{
		Coord x[3];
		for (int i = 0; i < 3; i++)
			x[i].resize(3);
		for (int i = 0; i < s[j].size(); i++)
		{
			fout << "v" ;
			for(l=0;l<3;l++)
			{
				x[i] = GetCoordSurf(s[j][i],m);
				fout << " " << x[i][l] ;
			}
			fout << endl;
		}

		int ind[4]={0,1,2,3};
		
		fout << "f " << k+ind[0] << " " << k+ind[1] << " " << k+ind[2] << endl;

		k+=s[j].size(); // номер следующей точки
		}
	fout.close();
}

Surfaces DivideSquares(SurfacePoints &sur, Points &mesh){
	Surfaces ret(2);
	ret[0].resize(3); ret[1].resize(3);
	Coord vex[4];
	
	if(sur[0].a_index != sur[1].a_index && sur[0].b_index != sur[1].b_index && sur[0].a_index != sur[1].b_index && sur[0].b_index != sur[1].a_index){
		SurfaceCoord t;
		t = sur[3];
		sur[3] = sur[1];
		sur[1] = t;
	}
	else if(sur[0].a_index != sur[2].a_index && sur[0].b_index != sur[2].b_index && sur[0].a_index != sur[2].b_index && sur[0].b_index != sur[2].a_index){
		SurfaceCoord t;
		t = sur[3];
		sur[3] = sur[2];
		sur[2] = t;
	}
	
	for(int i = 0; i < 4; ++i){
		vex[i] = GetCoordSurf(sur[i],mesh);
	}
	if(mydistance(vex[0],vex[3]) < mydistance(vex[1],vex[2])){
		ret[0][0].a_index=sur[0].a_index; ret[0][0].b_index=sur[0].b_index; ret[0][0].delta=sur[0].delta;
		ret[0][1].a_index=sur[1].a_index; ret[0][1].b_index=sur[1].b_index; ret[0][1].delta=sur[1].delta;
		ret[0][2].a_index=sur[2].a_index; ret[0][2].b_index=sur[2].b_index; ret[0][2].delta=sur[2].delta;
		ret[1][0].a_index=sur[1].a_index; ret[1][0].b_index=sur[1].b_index; ret[1][0].delta=sur[1].delta;
		ret[1][1].a_index=sur[2].a_index; ret[1][1].b_index=sur[2].b_index; ret[1][1].delta=sur[2].delta;
		ret[1][2].a_index=sur[3].a_index; ret[1][2].b_index=sur[3].b_index; ret[1][2].delta=sur[3].delta;
	}
	else{
		ret[0][0].a_index=sur[1].a_index; ret[0][0].b_index=sur[1].b_index; ret[0][0].delta=sur[1].delta;
		ret[0][1].a_index=sur[0].a_index; ret[0][1].b_index=sur[0].b_index; ret[0][1].delta=sur[0].delta;
		ret[0][2].a_index=sur[3].a_index; ret[0][2].b_index=sur[3].b_index; ret[0][2].delta=sur[3].delta;
		ret[1][0].a_index=sur[0].a_index; ret[1][0].b_index=sur[0].b_index; ret[1][0].delta=sur[0].delta;
		ret[1][1].a_index=sur[3].a_index; ret[1][1].b_index=sur[3].b_index; ret[1][1].delta=sur[3].delta;
		ret[1][2].a_index=sur[2].a_index; ret[1][2].b_index=sur[2].b_index; ret[1][2].delta=sur[2].delta;
	}
	
	return ret;
}

void RoundControl(SurfacePoints &s, Points &mesh){
	Coord vec1(3), vec2(3), vec3(3);
	Coord vex[3];
	bool detector;
	for(int i = 0; i < 3; ++i){
		vex[i] = GetCoordSurf(s[i],mesh);
	}
	
	vec1[0] = vex[1][0] - vex[0][0]; vec1[1] = vex[1][1] - vex[0][1]; vec1[2] = vex[1][2] - vex[0][2];
	vec2[0] = vex[2][0] - vex[0][0]; vec2[1] = vex[2][1] - vex[0][1]; vec2[2] = vex[2][2] - vex[0][2];
	if(s[0].delta == 0){
		vec3[0] = mesh[s[0].b_index][0] - vex[0][0];
		vec3[1] = mesh[s[0].b_index][1] - vex[0][1];
		vec3[2] = mesh[s[0].b_index][2] - vex[0][2];
		if(!isIn(mesh[s[0].b_index])){
			if(determinant(vec3,vec1,vec2)<0) {
				SurfaceCoord t;
				t = s[1];
				s[1] = s[2];
				s[2] = t;
			}
		}
		else{
			if(determinant(vec3,vec1,vec2)>0) {
				SurfaceCoord t;
				t = s[1];
				s[1] = s[2];
				s[2] = t;
			}
		}
	}
	else if(s[0].delta == 1){
		vec3[0] = mesh[s[0].a_index][0] - vex[0][0];
		vec3[1] = mesh[s[0].a_index][1] - vex[0][1];
		vec3[2] = mesh[s[0].a_index][2] - vex[0][2];
		if(!isIn(mesh[s[0].a_index])){
			if(determinant(vec3,vec1,vec2)<0) {
				SurfaceCoord t;
				t = s[1];
				s[1] = s[2];
				s[2] = t;
			}
		}
		else{
			if(determinant(vec3,vec1,vec2)>0) {
				SurfaceCoord t;
				t = s[1];
				s[1] = s[2];
				s[2] = t;
			}
		}
	}
	else{
		if(!isIn(mesh[s[0].a_index])){
			vec3[0] = mesh[s[0].a_index][0] - vex[0][0];
			vec3[1] = mesh[s[0].a_index][1] - vex[0][1];
			vec3[2] = mesh[s[0].a_index][2] - vex[0][2];
		}
		else{
			vec3[0] = mesh[s[0].b_index][0] - vex[0][0];
			vec3[1] = mesh[s[0].b_index][1] - vex[0][1];
			vec3[2] = mesh[s[0].b_index][2] - vex[0][2];
		}
		if(determinant(vec3,vec1,vec2)==0){
			vec1[0] = vex[0][0] - e_x0; vec1[1] = vex[0][1] - e_y0; vec1[2] = vex[0][2] - e_z0;
			vec2[0] = vex[1][0] - e_x0; vec2[1] = vex[1][1] - e_y0; vec2[2] = vex[1][2] - e_z0;
			vec3[0] = vex[2][0] - e_x0; vec3[1] = vex[2][1] - e_y0; vec3[2] = vex[2][2] - e_z0;
			if(determinant(vec3,vec1,vec2)<0){
				SurfaceCoord t;
				t = s[1];
				s[1] = s[2];
				s[2] = t;
			}
		}
		
		else if(determinant(vec3,vec1,vec2)<0){
			SurfaceCoord t;
			t = s[1];
			s[1] = s[2];
			s[2] = t;
		}
	}
// 	if(determinant(vec3,vec1,vec2)==0){ 
// 		cout << "Wuba\n";
// 		cout << vec1[0] << " " << vec1[1] << " " << vec1[2] << endl;
// 		cout << vec2[0] << " " << vec2[1] << " " << vec2[2] << endl;
// 		cout << vec3[0] << " " << vec3[1] << " " << vec3[2] << endl;
// 	}
}