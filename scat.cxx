#include <cmath>
#include <fstream>
#include <iostream>
#include <cstdlib>

using namespace std;

void path_init (double* const l, const double lambda0, const int N) {
  
  for (int i = 0; i<N; i++) l[i] = - lambda0 * log(double(1.0*rand()/RAND_MAX));
  
}

int main() {
  
  srand(time(NULL));
  const int N = 100000;
  const double d = 1;
  const int Gn = 100;
  const double v0 = 0.01;
  const double v1 = 1;
  const double Ra = 0.01;
  const double Rs = 0.3;
  double lambda0 = v0 / (Ra + Rs);
  double lambda1 = v1 / (Ra + Rs);
  
  double* l = new double[N];
  path_init ( l, lambda0, N);
  
  double t[Gn];
  double t2[Gn];
  double t_temp[Gn];
  for (int i = 0; i<Gn; i++) {
    t[i] = 0;
    t2[i] = 0;
    t_temp[i] = 0;
  }
  
  int n1;		// Nummer des Gitterpunktes, wo sich das Teilchen aufhält
  int n2;		// Nummer des Gitterpunktes, wo sich das Teilchen landet
  double r1;		// Restposition im Gitterpunkt
  double r2;		// Restposition im Zielgitterpunkt
  double w;
  double s;
  
  int p = 0;		// Teilchenindex
  double v = v0;	// Teilchengeschwindigkeit
  double x = 0;		// Teilchenort
  while (p<N) {
    if ( ((l[p]+x)>Gn*d) || ((l[p]+x)<0) ) {
      n1 = int(x/d);
      r1 = x - n1*d;
      if (l[p]>=0){ 
	t[n1] += (d-r1)/v;
	for (int i = n1+1; i<Gn; i++) t[i] += d/v;
      }
      else {
	t[n1] += r1/v;
	for (int i = n1-1; i>=0; i--) t[i] += d/v;
      }
      p++;
      v = v0;
      x = 0;
    }
    
    else {
      n1 = int(x/d);
      n2 = int( (l[p]+x) / d );
      r1 = x - n1*d;
      r2 = l[p] + x - n2*d;
      if (n2>n1) {
	t[n1] += (d-r1)/v;
	for (int i = n1+1; i<n2; i++) t[i] += d/v;
	t[n2] += r2/v;
      }
      else {
	if (n2<n1) {
	  t[n1] += r1/v;
	  for (int i = n1-1; i>n2; i--) t[i] += d/v;
	  t[n2] += (d-r2)/v;
	}
	else
	  t[n1] += abs(r1-r2) / v;
      }
      
      w = double(1.0*rand()/RAND_MAX);
      if ( w <= Ra / (Ra+Rs) ) {
	p++;
	v = v0;
	x = 0;
      }
      else {
	x += l[p];
	s = 2.0*double(1.0*rand()/RAND_MAX) - 1.0;
	s = (s > 0) - (s < 0);
	l[p] = - lambda1 * log(double(1.0*rand()/RAND_MAX)) * s;
	v = v1;
      }
      
    }
    
  }
  
  double t_max = t[0];
  for (int i = 0; i<Gn; i++) t[i] /= N;
  
  string name = "out";
  ofstream out(name.c_str());
  for (int i = 0; i<Gn; i++)
    out << i << "\t" << t[i] << endl;
  out.close();
  
  delete[] l;
  return 0;
}