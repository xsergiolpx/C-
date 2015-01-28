#include <stdio.h>
#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
using namespace std;

//Variables globales
double h=0.1;
int tam=100; //tamano de los vectores
int espacio=16; //espaciado enre columnas

void runge(vector<double> &t, vector<double> &x, vector<double> &y, vector<double> &y4);
double fx(double t, double x, double y);
double fy(double t, double x, double y);

int main(int argc, char **argv)
{
	vector<double> y (tam); //solucion de orden 5
	vector<double> x (tam);
        vector<double> t (tam);
	vector<double> y4 (tam); //solucion de orden 4
	//Condiciones de contorno
	y[0]=1;
	x[0]=0;
	y4[0]=y[0];
        t[0]=0;
	cout << endl;
	runge(t,x,y,y4);
	cout << "\n\nPresiona enter para salir ";
	cin.ignore();
	return 0;
}

//Valor de la derivada x''=-6x'+7x
double fx(double t, double x, double y){
	return (y);
}

double fy(double t, double x, double y){
	return (-6*y+7*x);
}

//Runge-kutta-Felberg
void runge(vector<double> &t, vector<double> &x, vector<double> &y, vector<double> &y4){
	double kx1, kx2, kx3, kx4, kx5, kx6;
        double ky1, ky2, ky3, ky4, ky5, ky6;
	double error = 0;
	//cout << "################## Metodo de Runge-Kutta- Fehlberg ##################\n\n";
	//cout << "donde y5 es la solucion de 5º orden, y y4 la de 4º orden\n\n";
	//cout << setw(4) << right << "x" << right << setw(espacio) << "y5" << right << setw(espacio) << "y4" << right << setw(espacio) << "error\n\n";
	for(signed i=0;i<tam-1;i++){
		kx1=h*fx(t[i],x[i],y[i]);
                ky1=h*fy(t[i],x[i],y[i]);
		kx2=h*fx(t[i]+h/2,x[i]+kx1/2,y[i]+ky1/2);
                ky2=h*fy(t[i]+h/2,x[i]+kx1/2,y[i]+ky1/2);
		kx3=h*fx(t[i]+3/8*h,x[i]+3*kx1/32+9*kx2/32,y[i]+3*ky1/32+9*ky2/32);
                ky3=h*fy(t[i]+3/8*h,x[i]+3*kx1/32+9*kx2/32,y[i]+3*ky1/32+9*ky2/32);
		kx4=h*fx(t[i]+12*h/13,x[i]+1932*kx1/2197-7200*kx2/2197+7296*kx3/2197,y[i]+1932*ky1/2197-7200*ky2/2197+7296*ky3/2197);
                ky4=h*fy(t[i]+12*h/13,x[i]+1932*kx1/2197-7200*kx2/2197+7296*kx3/2197,y[i]+1932*ky1/2197-7200*ky2/2197+7296*ky3/2197);
                kx5=h*fx(t[i]+h,x[i]+439*kx1/216-8*kx2+3680*kx3/513-845*kx4/4104,y[i]+439*ky1/216-8*ky2+3680*ky3/513-845*ky4/4104);
		ky5=h*fy(t[i]+h,x[i]+439*kx1/216-8*kx2+3680*kx3/513-845*kx4/4104,y[i]+439*ky1/216-8*ky2+3680*ky3/513-845*ky4/4104);
                kx6=h*fx(t[i]+h/2,x[i]-8*kx1/27+2*kx2-3544*kx3/2565+1859*kx4/4104-11*kx5/40,y[i]-8*ky1/27+2*ky2-3544*ky3/2565+1859*ky4/4104-11*ky5/40);
		ky6=h*fy(t[i]+h/2,x[i]-8*kx1/27+2*kx2-3544*kx3/2565+1859*kx4/4104-11*kx5/40,y[i]-8*ky1/27+2*ky2-3544*ky3/2565+1859*ky4/4104-11*ky5/40);
                
                //y4[i+1]=y[i]+25*k1/216+1408*k3/2565+2197*k4/4104-k5/5;
		y[i+1]=y[i]+16*ky1/135+6656*ky3/12825+28561*ky4/56430-9*ky5/50+2*ky6/55;
                x[i+1]=x[i]+16*kx1/135+6656*kx3/12825+28561*kx4/56430-9*kx5/50+2*kx6/55;
		t[i+1]=t[0]+(i+1)*h;
                cout << "i  = " << i << " ##  x =  " << x[i]<< endl;
                //error += abs(abs(x[i])-abs(y4[i]));
                //Descomentar la siguiente linea para ver todos los valores
                //cout << setw(4) << right<< x[i] << right << setw(espacio) << y[i] << right << setw(espacio) << y4[i]<< right << setw(espacio) << error << endl;
	}
        //cout << setw(4) << right<< x[tam-2] << right << setw(espacio) << y[tam-2] << right << setw(espacio) << y4[tam-2]<< right << setw(espacio) << error << endl;
	
}
