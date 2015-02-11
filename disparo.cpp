//Runge-Kutta Felhberg funciona bien, pero queda llamar a las funciones con punteros de funciones

#include <stdio.h>
#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <string>

using namespace std;

//Variables globales
double h=0.2; //precision
int espacio=16; //espaciado enre columnas
double tf=3;

//Prototipos
void runge(vector<double> &t, vector<double> &x, vector<double> &y);
void runge2(vector<double> &t, vector<double> &x, vector<double> &y);
double f(double t, double x, double y);
double g(double t, double x, double y);
double g2(double t, double x, double y);
void disparo(vector<double> &t, vector<double> &x, vector<double> &y, double S2, double Rd);
void disparo2(vector<double> &t, vector<double> &x, vector<double> &y, double S2, double Rd);

//Mostrar valores
void mostrar(vector<double> &t, vector<double> &x, vector<double> &y)
{
	for (unsigned i=0; i<t.size()-1; i++) {
		cout  << right << setw(espacio) << t[i] << right << setw(espacio) << x[i] << right << setw(espacio) << y[i] << endl;
	}
	cout << "--------------------------------------------------------------------------------------------------------\n";
}

//Inicio del programa
int main(int argc, char **argv)
{
	int tam=1; //tamano de los vectores iniciales
	vector<double> y (tam); //derivada de orden 5
	vector<double> x (tam); //solucion de orden 5
	vector<double> t (tam); //variable independiente

	//Condiciones de contorno
	x[0]=2;
	y[0]=-1.5;
	t[0]=1;
	disparo(t,x,y,-3,-1);

	//Caso no lineal
	x[0]=1;
	y[0]=-1.5;
	t[0]=1;
	disparo2(t,x,y,-3,-2);

	cout << "\n\nPresiona enter para salir ";
	cin.ignore();
	return 0;
}

//Planteamiento de la ecuacion diferencia como sistema de ecuaciones
double f(double t, double x, double y)
{
	return (y);
}

double g(double t, double x, double y)
{
	return (t+(1-t/5.0)*x);
}

double g2(double t, double x, double y)
{
	return(t+(1-t/5)*y*x);
}

//Runge-kutta-Felberg
void runge2(vector<double> &t, vector<double> &x, vector<double> &y)
{
	double k1,k2,k3,k4,k5,k6,m1,m2,m3,m4,m5,m6;
	signed i = 0;
	x.resize(1);
	y.resize(1);
	t.resize(1);
	while(t[i-1]<tf) {
		//Redimensionamos los vectores para que quepa un elemento mas
		x.resize(x.size()+1);
		y.resize(y.size()+1);
		t.resize(t.size()+1);

		k1 = h*f(t[i],x[i],y[i]);
		m1 = h*g2(t[i],x[i],y[i]);
		k2 = h*f(t[i]+h/4.0,x[i]+k1/4.0,y[i]+m1/4.0);
		m2 = h*g2(t[i]+h/4.0,x[i]+k1/4.0,y[i]+m1/4.0);
		k3 = h*f(t[i]+3.0*h/8.0,x[i]+3.0*k1/32.0+9.0*k2/32.0,y[i]+3.0*m1/32.0+9.0*m2/32.0);
		m3 = h*g2(t[i]+3.0*h/8.0,x[i]+3.0*k1/32.0+9.0*k2/32.0,y[i]+3.0*m1/32.0+9.0*m2/32.0);
		k4 = h*f(t[i]+12.0*h/13.0,x[i]+1932.0*k1/2197.0-7200.0*k2/2197.0+7296.0*k3/2197.0,y[i]+1932.0*m1/2197.0-7200.0*m2/2197.0+7296.0*m3/2197.0);
		m4 = h*g2(t[i]+12.0*h/13.0,x[i]+1932.0*k1/2197.0-7200.0*k2/2197.0+7296.0*k3/2197.0,y[i]+1932.0*m1/2197.0-7200.0*m2/2197.0+7296.0*m3/2197.0);
		k5 = h*f(t[i]+h,x[i]+439.0*k1/216.0-8.0*k2+3680.0*k3/513.0-845.0*k4/4104.0,y[i]+439.0*m1/216.0-8.0*m2+3680.0*m3/513.0-845.0*m4/4104.0);
		m5 = h*g2(t[i]+h,x[i]+439.0*k1/216.0-8.0*k2+3680.0*k3/513.0-845.0*k4/4104.0,y[i]+439.0*m1/216.0-8.0*m2+3680.0*m3/513.0-845.0*m4/4104.0);
		k6 = h*f(t[i]+h/2.0,x[i]-8.0*k1/27.0+2.0*k2-3544.0*k3/2565.0+1859.0*k4/4104.0-11*k5/40,y[i]-8.0*m1/27.0+2.0*m2-3544.0*m3/2565.0+1859.0*m4/4104.0-11*m5/40);
		m6 = h*g2(t[i]+h/2.0,x[i]-8.0*k1/27.0+2.0*k2-3544.0*k3/2565.0+1859.0*k4/4104.0-11*k5/40,y[i]-8.0*m1/27.0+2.0*m2-3544.0*m3/2565.0+1859.0*m4/4104.0-11*m5/40);
		x[i+1] = x[i] + (16.0*k1/135.0+6656.0*k3/12825.0+28561.0*k4/56430.0-9.0*k5/50.0+2*k6/55.0);
		y[i+1] = y[i] + (16.0*m1/135.0+6656.0*m3/12825.0+28561.0*m4/56430.0-9.0*m5/50.0+2*m6/55.0);

		//Damos un paso en el tiempo
		t[i+1]=t[0]+(i+1)*h;
		//Aumentamos contador
		i += 1;
	}
}

void runge(vector<double> &t, vector<double> &x, vector<double> &y)
{
	double k1,k2,k3,k4,k5,k6,m1,m2,m3,m4,m5,m6;
	signed i = 0;
	x.resize(1);
	y.resize(1);
	t.resize(1);
	while(t[i-1]<tf) {
		//Redimensionamos los vectores para que quepa un elemento mas
		x.resize(x.size()+1);
		y.resize(y.size()+1);
		t.resize(t.size()+1);
		k1 = h*f(t[i],x[i],y[i]);
		m1 = h*g(t[i],x[i],y[i]);
		k2 = h*f(t[i]+h/4.0,x[i]+k1/4.0,y[i]+m1/4.0);
		m2 = h*g(t[i]+h/4.0,x[i]+k1/4.0,y[i]+m1/4.0);
		k3 = h*f(t[i]+3.0*h/8.0,x[i]+3.0*k1/32.0+9.0*k2/32.0,y[i]+3.0*m1/32.0+9.0*m2/32.0);
		m3 = h*g(t[i]+3.0*h/8.0,x[i]+3.0*k1/32.0+9.0*k2/32.0,y[i]+3.0*m1/32.0+9.0*m2/32.0);
		k4 = h*f(t[i]+12.0*h/13.0,x[i]+1932.0*k1/2197.0-7200.0*k2/2197.0+7296.0*k3/2197.0,y[i]+1932.0*m1/2197.0-7200.0*m2/2197.0+7296.0*m3/2197.0);
		m4 = h*g(t[i]+12.0*h/13.0,x[i]+1932.0*k1/2197.0-7200.0*k2/2197.0+7296.0*k3/2197.0,y[i]+1932.0*m1/2197.0-7200.0*m2/2197.0+7296.0*m3/2197.0);
		k5 = h*f(t[i]+h,x[i]+439.0*k1/216.0-8.0*k2+3680.0*k3/513.0-845.0*k4/4104.0,y[i]+439.0*m1/216.0-8.0*m2+3680.0*m3/513.0-845.0*m4/4104.0);
		m5 = h*g(t[i]+h,x[i]+439.0*k1/216.0-8.0*k2+3680.0*k3/513.0-845.0*k4/4104.0,y[i]+439.0*m1/216.0-8.0*m2+3680.0*m3/513.0-845.0*m4/4104.0);
		k6 = h*f(t[i]+h/2.0,x[i]-8.0*k1/27.0+2.0*k2-3544.0*k3/2565.0+1859.0*k4/4104.0-11*k5/40,y[i]-8.0*m1/27.0+2.0*m2-3544.0*m3/2565.0+1859.0*m4/4104.0-11*m5/40);
		m6 = h*g(t[i]+h/2.0,x[i]-8.0*k1/27.0+2.0*k2-3544.0*k3/2565.0+1859.0*k4/4104.0-11*k5/40,y[i]-8.0*m1/27.0+2.0*m2-3544.0*m3/2565.0+1859.0*m4/4104.0-11*m5/40);
		x[i+1] = x[i] + (16.0*k1/135.0+6656.0*k3/12825.0+28561.0*k4/56430.0-9.0*k5/50.0+2*k6/55.0);
		y[i+1] = y[i] + (16.0*m1/135.0+6656.0*m3/12825.0+28561.0*m4/56430.0-9.0*m5/50.0+2*m6/55.0);

		//Damos un paso en el tiempo
		t[i+1]=t[0]+(i+1)*h;
		//Aumentamos contador
		i += 1;
	}
}

void disparo(vector<double> &t, vector<double> &x, vector<double> &y, double S2, double Rd)
{
	double S1=y[0];
	runge(t,x,y);
	double R1=x[x.size()-2];
	y[0]=S2;
	runge(t,x,y);
	double R2=x[x.size()-2];
	y[0]=S2+(Rd-R2)*(S2-S1)/(R2-R1);
	do {
		S1=S2;
		S2=y[0];
		R1=R2;
		runge(t,x,y);
		R2=x[x.size()-2];
		y[0]=S2+(Rd-R2)*(S2-S1)/(R2-R1);
	} while(abs(abs(Rd)-abs(R2))>0.000001);
	cout << "\n\nCaso lineal:    Pendiente S = " << y[0] << "\n\n";
	mostrar(t,x,y);
}


void disparo2(vector<double> &t, vector<double> &x, vector<double> &y, double S2, double Rd)
{

	double S1=y[0];
	runge2(t,x,y);
	double R1=x[x.size()-2];
	y[0]=S2;
	runge2(t,x,y);
	double R2=x[x.size()-2];
	y[0]=S2+(Rd-R2)*(S2-S1)/(R2-R1);
	do {
		S1=S2;
		S2=y[0];
		R1=R2;
		runge2(t,x,y);
		R2=x[x.size()-2];
		y[0]=S2+(Rd-R2)*(S2-S1)/(R2-R1);

	} while(abs(abs(Rd)-abs(R2))>0.000001);
	cout << "\n\nCaso no lineal:    Pendiente S = " << y[0] << "\n\n";
	mostrar(t,x,y);
}
