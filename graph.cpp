#include "graph.h"
#include <list>
#include <stack>

using std::list;
using std::stack;

void qsortmy1(vector<Rib> &f, int lb, int rb, int coordnum){
     if(lb >= rb) return;

    Rib tmp; 
    int l, r, m;
     l = lb; r = rb; m = (lb + rb) >> 1;

    while(l < r){
         while((f[l].a[coordnum] <= f[m].a[coordnum]) && (l < m)) l++;
         while((f[r].a[coordnum] >= f[m].a[coordnum]) && (r > m)) r--;
         tmp = f[l]; f[l] = f[r]; f[r] = tmp;
         if(l == m) m = r;
         else if(r == m) m = l;
    }
     qsortmy1(f, lb, m - 1, coordnum);
     qsortmy1(f, m + 1, rb, coordnum);
}

void qsortmy2(vector<Rib> &f, int lb, int rb, int coordnum){
     if(lb >= rb) return;

    Rib tmp; 
    int l, r, m;
     l = lb; r = rb; m = (lb + rb) >> 1;

    while(l < r){
         while((f[l].b[coordnum] <= f[m].b[coordnum]) && (l < m)) l++;
         while((f[r].b[coordnum] >= f[m].b[coordnum]) && (r > m)) r--;
         tmp = f[l]; f[l] = f[r]; f[r] = tmp;
         if(l == m) m = r;
         else if(r == m) m = l;
    }
     qsortmy2(f, lb, m - 1,coordnum);
     qsortmy2(f, m + 1, rb,coordnum);
}


void qsorttriangle(vector<SimplifiedRib> &a, int lb, int rb){
     if(lb >= rb) return;

    SimplifiedRib tmp; 
    int l, r, m;
     l = lb; r = rb; m = (lb + rb) >> 1;

    while(l < r){
         while((a[l].TriangleNumber <= a[m].TriangleNumber) && (l < m)) l++;
         while((a[r].TriangleNumber >= a[m].TriangleNumber) && (r > m)) r--;
         tmp = a[l]; a[l] = a[r]; a[r] = tmp;
         if(l == m) m = r;
         else if(r == m) m = l;
    }
     qsorttriangle(a, lb, m - 1);
     qsorttriangle(a, m + 1, rb);
}

/*void qsorttriangletest(vector<Rib> &a, int lb, int rb){
     if(lb >= rb) return;

    Rib tmp; 
    int l, r, m;
     l = lb; r = rb; m = (lb + rb) >> 1;

    while(l < r){
         while((a[l].TriangleNumber <= a[m].TriangleNumber) && (l < m)) l++;
         while((a[r].TriangleNumber >= a[m].TriangleNumber) && (r > m)) r--;
         tmp = a[l]; a[l] = a[r]; a[r] = tmp;
         if(l == m) m = r;
         else if(r == m) m = l;
    }
     qsorttriangletest(a, lb, m - 1);
     qsorttriangletest(a, m + 1, rb);
}*/

vector<SimplifiedRib> CreateInitialStructure(Surfaces& s, PointsXYZ &m){
	vector<Rib> k;
	k.resize(s.size()*3);
	//cout << s.size() << " " << k.size() << endl;
	for(int i = 0; i < s.size(); i++){ // заполняем k ребрами с номерами треугольников
		k[i*3].TriangleNumber = k[i*3+1].TriangleNumber = k[i*3+2].TriangleNumber = i;
		k[i*3].a = k[i*3+1].a = GetCoordSurf(s[i][0],m);
		k[i*3].b = k[i*3+2].a = GetCoordSurf(s[i][1],m);
		k[i*3+1].b = k[i*3+2].b = GetCoordSurf(s[i][2],m);
	}
	
	for(int i = 0; i < k.size(); i++){ //упорядочиваем вершины каждого ребра по величине модуля
		CoordXYZ t;
	/*	cout << i << endl;
		cout << "a: " << k[i].a[0] << " " << k[i].a[1] << " "<< k[i].a[2] << " ";
		cout << "b: " << k[i].b[0] << " " << k[i].b[1] << " "<< k[i].b[2] << endl;
		if(!((i+1)%3)) cout << k[i].TriangleNumber << endl;*/
		if(keisu(k[i].a) > keisu(k[i].b)) {
			t = k[i].a;
			k[i].a = k[i].b;
			k[i].b = t;
		}
	}
	
	qsortmy1(k,0,k.size()-1,0);
	int lb = 0;
	int rb = 0;
	for(int i = 1; i < k.size(); i++){
		
		if(k[i].a[0] != k[i-1].a[0]){
			rb = i-1;
			qsortmy1(k,lb,rb,1);
			lb = i;
		}
		else if(i == k.size()-1){
			rb = i;
			qsortmy1(k,lb,rb,1);
		}
	}
	
	lb = rb = 0;
	for(int i = 1; i < k.size(); i++){
		if(k[i].a[0] != k[i-1].a[0] || k[i].a[1] != k[i-1].a[1]){
			rb = i-1;
			qsortmy1(k,lb,rb,2);
			lb = i;
		}
		else if(i == k.size()-1){
			rb = i;
			qsortmy1(k,lb,rb,2);
		}
	}

	lb = rb = 0;
	for(int i = 1; i < k.size(); i++){
		if(k[i].a[0] != k[i-1].a[0] || k[i].a[1] != k[i-1].a[1] || k[i].a[2] != k[i-1].a[2]){
			rb = i-1;
			qsortmy2(k,lb,rb,0);
			lb = i;
		}
		else if(i == k.size()-1){
			rb = i;
			qsortmy2(k,lb,rb,0);
		}
	}
	
	lb = rb = 0;
	for(int i = 1; i < k.size(); i++){
		if(k[i].a[0] != k[i-1].a[0] || k[i].a[1] != k[i-1].a[1] || k[i].a[2] != k[i-1].a[2] || k[i].b[0] != k[i-1].b[0]){
			rb = i-1;
			qsortmy2(k,lb,rb,1);
			lb = i;
		}
		else if(i == k.size()-1){
			rb = i;
			qsortmy2(k,lb,rb,1);
		}
	}
	
	lb = rb = 0;
	for(int i = 1; i < k.size(); i++){
		if(k[i].a[0] != k[i-1].a[0] || k[i].a[1] != k[i-1].a[1] || k[i].a[2] != k[i-1].a[2] || k[i].b[0] != k[i-1].b[0] || k[i].b[1] != k[i-1].b[1]){
			rb = i-1;
			qsortmy2(k,lb,rb,2);
			lb = i;
		}
		else if(i == k.size()-1){
			rb = i;
			qsortmy2(k,lb,rb,2);
		}
	}

	//qsorttriangletest(k,0,k.size()-1);
	/*for(int i = 0; i < k.size();i++){
		if(k[i].TriangleNumber == 3779){
		cout << i << " " << k[i].TriangleNumber << endl;
		cout << "a: " << k[i].a[0] << " " << k[i].a[1] << " "<< k[i].a[2] << " ";
		cout << "b: " << k[i].b[0] << " " << k[i].b[1] << " "<< k[i].b[2] << endl;
		}
	}*/
	vector<SimplifiedRib> ret;
	ret.resize(k.size());
	ret[0].RibNumber = 0; ret[0].TriangleNumber = k[0].TriangleNumber;
	int z = 0;
	for(int i = 1; i < k.size(); i++){
		if(k[i].a != k[i-1].a || k[i].b != k[i-1].b) z++;
		ret[i].RibNumber = z; ret[i].TriangleNumber = k[i].TriangleNumber;
	}
	//qsorttriangle(ret,0,ret.size()-1);
	/*for(int i = 0; i < ret.size(); i++){
		if(ret[i].RibNumber == 6805){
		cout << i << " ";
		cout << "RibNumber: " << ret[i].RibNumber << " ";
		cout << "TriangleNumber: " << ret[i].TriangleNumber << endl;
		}
	}*/
	return ret;
}

bool CheckSide(vector<SimplifiedRib> r){
	int k = 1;
	bool flag = true;
	int z = 0;
	for(int i = 1; i < r.size(); i++){
		if(r[i].RibNumber != r[i-1].RibNumber) {
			if(k!=2) { cout << "One side conjuration! Triangle: " << r[i-1].TriangleNumber << endl; flag = false; z++; }
			k = 0;
		}
		k++;
	}
	cout << "Number of ribs conjucted on one side: " << z << endl;
	return flag;
}

bool CheckGraphCoherency(vector<SimplifiedRib> r, int n){
	vector< list<int> > graph(n);

	for(int i = 0; i < r.size()-1; i++){
		if(r[i].RibNumber == r[i+1].RibNumber){
			graph[r[i+1].TriangleNumber].push_back(r[i].TriangleNumber);
			graph[r[i].TriangleNumber].push_back(r[i+1].TriangleNumber);
		}
	}
	list<int>::iterator j;
	/*for(int i = 0; i < graph.size(); i++){
		cout << i << " ";
		for(j = graph[i].begin(); j != graph[i].end(); ++j){
			cout << *j << " ";
		}
		cout << endl;
		if(graph[i].size() != 3) cout << "Warning: " << i << endl;
	}*/
	stack<int> Q;
	vector<bool> isChecked(n,false);
	isChecked[0] = true; Q.push(0);
	
	while(Q.size()){
		int t = Q.top();
		Q.pop();
		for(j = graph[t].begin(); j != graph[t].end(); ++j){
			if(!isChecked[*j]) { isChecked[*j] = true; Q.push(*j); }
		}
	}
	int s = 0;
	for(int i = 0; i < isChecked.size(); i++){
		if(isChecked[i]) s++;
		else cout << "Triangle with troubles: " << i << endl;
	}
	cout << "Number of triangles not in main component: " << n - s << endl;
	if(s == n) return true;
	return false;
}
