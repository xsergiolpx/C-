/*
 *		Sergio Ballesteros
 */

#include <stdio.h>
#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <string>

using namespace std;
void diferencias(vector<double>& t, vector<double>& x);
void tridiagonal(vector<double>& c, vector<double>& b, vector<double>& a, vector<double>& d, vector<double>& x);
void mostrar(vector<double>& u, vector<double>& v);

double t0 = 1, tf = 3;
double x0 = 1, xf = -2;
double h = 0.2;

int main(int argc, char** argv)
{
	int tam = (tf - t0) / h + 1;
	vector<double> x(tam - 2);
	vector<double> t(tam);
	diferencias(t, x);
	mostrar(t, x);
	return 0;
}

void mostrar(vector<double>& u, vector<double>& v)
{
	int espacio = 12;
	int precision = 4;
	cout  << right << setw(espacio) << "x" << right << setw(espacio) << "u" << endl;
	cout << endl;
	for(int i = 0; i < v.size(); i++) {
		cout  << right << setw(espacio) << setprecision(precision) << u[i + 1] << right << setw(espacio) << setprecision(precision) << v[i] << endl;
	}
}

void diferencias(vector<double>& t, vector<double>& x)
{
	for(int i = 0; i < t.size(); i++) {
		t[i] = t0 + h * (i);
	}
	int tamDiag = t.size() - 2;
	vector<double> c(tamDiag); //diagonal superior
	vector<double> b(tamDiag); //diagonal principal
	vector<double> a(tamDiag); //diagonal inferior
	vector<double> d(tamDiag); //terminos independientes
	for(int i = 0; i < tamDiag; i++) {
		c[i] = 1;
		b[i] = -(2 + h * h * (1 - t[i + 1] / 5));
		a[i] = 1;
		d[i] = h * h * t[i + 1];
	}
	d[0] = d[0] - x0;
	d[tamDiag - 1] = d[tamDiag - 1] - xf;
	tridiagonal(c, b, a, d, x);
}

void tridiagonal(vector<double>& c, vector<double>& b, vector<double>& a, vector<double>& d, vector<double>& x)
{
	int N = b.size();
	vector<double> c_star(N);
	vector<double> d_star(N);
	c_star[0] = c[0] / b[0];
	d_star[0] = d[0] / b[0];
	for (int i = 1; i < N; i++) {
		double m = 1.0 / (b[i] - a[i] * c_star[i - 1]);
		c_star[i] = c[i] * m;
		d_star[i] = (d[i] - a[i] * d_star[i - 1]) * m;
	}
	x[N - 1] = d_star[N - 1];
	for (int i = N - 1; i-- > 0; ) {
		x[i] = d_star[i] - c_star[i] * x[i + 1];
	}
}
