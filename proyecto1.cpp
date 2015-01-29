#include <stdio.h>
#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <fstream>

using namespace std;

//Variables globales
double h=0.001;
int tam=2; //tamano de los vectores
int espacio=16; //espaciado enre columnas

//Parametros del problema
double Cd = 0.47;
double Pf = 1.29;
double A = 3.14*0.1*0.1;
double g = 9.81;
double m = 1; //Cambiar

void runge(vector<double> &t, vector<double> &x, vector<double> &y, vector<double> &y4, vector<double> &x4);
double fx(double t, double x, double y);
double fy(double t, double x, double y);

int main(int argc, char **argv)
{
	vector<double> y (tam); //derivada de orden 5
	vector<double> x (tam); //solucion de orden 5
	vector<double> t (tam); //variable independiente
	vector<double> y4 (tam); //derivada de orden 4
	vector<double> x4 (tam); //solucion de orden 4
	//Condiciones de contorno
	y[0]=0;
	x[0]=-1000;
	y4[0]=y[0];
	x4[0]=x[0];
	t[0]=0;
	cout << endl;
	runge(t,x,y,y4,x4);
	cout << "\n\nPresiona enter para salir ";
	cin.ignore();
	return 0;
}

//Valor de la derivada
double fx(double t, double x, double y){
	return (y);
}

double fy(double t, double x, double y){
	return (g-1/(2*m)*Cd*Pf*A*y*y);
}

//Runge-kutta-Felberg
void runge(vector<double> &t, vector<double> &x, vector<double> &y, vector<double> &y4, vector<double> &x4){
	double kx1, kx2, kx3, kx4, kx5, kx6;
	double ky1, ky2, ky3, ky4, ky5, ky6;
	double errorx = 0;
	double errory = 0;
	signed i = 0;
	cout << "################## Caida de un objeto ##################\n\n";
	cout << "Con: x5, x4 la posicion del objeto (soluciones de 5º y 4º orden)\n\n";
	cout << "Con: y5, y4 la velocidad del objeto (soluciones de 5º y 4º orden)\n\n";
	cout  << right << setw(espacio) << "t" << right << setw(espacio) << "x5" << right << setw(espacio) << "x4" << right << setw(espacio) << "y5"
	<< right << setw(espacio) << "y4" << right << setw(espacio) << "errorx" << right << setw(espacio) << "errory\n\n";
	
	//Comienza la salida de datos
	ofstream datos;
	datos.open ("datos.txt");
	while(x[i]<0){
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
		
		//Soluciones
        	y4[i+1]=y[i]+25*ky1/216+1408*ky3/2565+2197*ky4/4104-ky5/5;
		x4[i+1]=x[i]+25*kx1/216+1408*kx3/2565+2197*kx4/4104-kx5/5;
		y[i+1]=y[i]+16*ky1/135+6656*ky3/12825+28561*ky4/56430-9*ky5/50+2*ky6/55;
		x[i+1]=x[i]+16*kx1/135+6656*kx3/12825+28561*kx4/56430-9*kx5/50+2*kx6/55;
		//Errores
        	errorx += abs(abs(x[i])-abs(x4[i]));
        	errory += abs(abs(y[i])-abs(y4[i]));
		//Archivo de texto
		datos << t[i] << "	" << x[i] << "	" << y[i] << endl;
        	cout  << right << setw(espacio) << t[i] << right << setw(espacio) << x[i] << right << setw(espacio) << x4[i] << right << setw(espacio) << y[i] << right << setw(espacio) << y4[i] << right << setw(espacio) << errorx << right << setw(espacio) << errory<< endl;
		//Damos un paso
		t[i+1]=t[0]+(i+1)*h;
		//Redimensionamos los vectores para que quepa un elemento mas
		x.resize(x.size()+1);
		x4.resize(x4.size()+1);
		y.resize(y.size()+1);
		y4.resize(y4.size()+1);
		t.resize(t.size()+1);
		//Aumentamos contador
		i += 1;
	}
	datos.close();
}

  
