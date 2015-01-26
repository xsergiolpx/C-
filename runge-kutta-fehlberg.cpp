#include <stdio.h>
#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
using namespace std;
//Variables globales
double h=0.1;
int tam=32; //tamano de los vectores
int espacio=16; //espaciado enre columnas

void runge(vector<double> &x, vector<double> &y, vector<double> &y4);
double f(double x, double y);

int main(int argc, char **argv)
{
	vector<double> y (tam); //solucion de orden 5
	vector<double> x (tam);
	vector<double> y4 (tam); //solucion de orden 4
	//Condiciones de contorno
	y[0]=-1;
	y4[0]=-1;
	x[0]=0;
	cout << endl;
	runge(x,y,y4);
	cout << "\n\nPresiona enter para salir ";
	cin.ignore();
	return 0;
}

//Valor de la derivada
double f(double x, double y){
	return (-2*x-y);
}

//Runge-kutta-Felberg
void runge(vector<double> &x, vector<double> &y, vector<double> &y4){
	double k1, k2, k3, k4, k5, k6;
	double error = 0;
	cout << "################## Metodo de Runge-Kutta- Fehlberg ##################\n\n";
	cout << "donde y5 es la solucion de 5º orden, y y4 la de 4º orden\n\n";
	cout << setw(4) << right << "x" << right << setw(espacio) << "y5" << right << setw(espacio) << "y4" << right << setw(espacio) << "error\n\n";
	for(signed i=0;i<tam-1;i++){
		k1=h*f(x[i],y[i]);
		k2=h*f(x[i]+h/4,y[i]+k1/4);
		k3=h*f(x[i]+3*h/8,y[i]+3*k1/32+9*k2/32);
		k4=h*f(x[i]+12*h/13,y[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197);
		k5=h*f(x[i]+h,y[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104);
		k6=h*f(x[i]+h/2,y[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40);
		y4[i+1]=y[i]+25*k1/216+1408*k3/2565+2197*k4/4104-k5/5;
		y[i+1]=y[i]+16*k1/135+6656*k3/12825+28561*k4/56430-9*k5/50+2*k6/55;
		x[i+1]=x[0]+(i+1)*h;
                error += abs(abs(y[i])-abs(y4[i]));
                //Descomentar la siguiente linea para ver todos los valores
                //cout << setw(4) << right<< x[i] << right << setw(espacio) << y[i] << right << setw(espacio) << y4[i]<< right << setw(espacio) << error << endl;
	}
        cout << setw(4) << right<< x[tam-2] << right << setw(espacio) << y[tam-2] << right << setw(espacio) << y4[tam-2]<< right << setw(espacio) << error << endl;
	
}
