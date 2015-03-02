
#include <stdio.h>
#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <string>

using namespace std;

//Variables globales
double h=0.05; //precision
int espacio=16; //espaciado enre columnas
const double pi = atan(1)*4;

//Parametros del problema
double Pf = 1.29;
double g = 9.81;

//Prototipos
//Pasar como variables tiempo, posicion, velocidad, posicion y velocidad para RKF4, masa, nombre del archivo de salida.
void runge(vector<double> &t, vector<double> &x, vector<double> &y, vector<double> &y4, vector<double> &x4, double m, double Cd, double A,string nombre);
//Sistema de ecuaciones
double fx(double t, double x, double y);
double fy(double t, double x, double y, double m);
//Calcular la masa con la densidad y radio de la esfera
double masa(double p, double r);

//Inicio del programa
int main(int argc, char **argv)
{	
	int tam=1; //tamano de los vectores iniciales
	vector<double> y (tam); //derivada de orden 5
	vector<double> x (tam); //solucion de orden 5
	vector<double> t (tam); //variable independiente
	vector<double> y4 (tam); //derivada de orden 4
	vector<double> x4 (tam); //solucion de orden 4

	//Condiciones de contorno
	y[0]=0; //velocidad inicial
	x[0]=-1000; //comienza en 1000 metros de altura, pero se tomo como el eje negativo
	y4[0]=y[0]; //igual para RKF4
	x4[0]=x[0];
	t[0]=0; //empieza en t=0
	cout << endl;
	
	//Se calculan todos los casos
	cout << "Valores finales antes de llegar al suelo en unidades del S.I\n\n t = tiempo\n x5/x4 = posicion final con RKF5/4\n y5/y4 = velocidad final con RKF5/4\n errorx/y = error acumulativo en la posicion y velocidad entre RKF 4 y 5\n\n";
	
	//Esfera maciza de hierro
	string nombre="esfera-hierro-maciza";
	double r=0.1; //radio de la esfera
	double p=7874; //densidad de la esfera
	double m, m1, m2; //masa
	double Cd = 0.47; //coeficiente de arrastre
	double A=pi*r*r; //area transversal
	m1=masa(p,r);
	runge(t,x,y,y4,x4,m1,Cd,A,nombre);

	//Esfera hueca de hierro
	nombre="esfera-hierro-hueca";
	double r2=0.09; //radio interno
	m2=masa(p,r)-masa(p,r2);
	runge(t,x,y,y4,x4,m2,Cd,A,nombre);

	//Esfera de madera de pino
	nombre="esfera-madera-pino";
	p=500;
	m = masa(p,r);
	runge(t,x,y,y4,x4,m,Cd,A,nombre);

	//Balon
	nombre="balon-futbol";
	m = 0.4;
	runge(t,x,y,y4,x4,m,Cd,A,nombre);
	
	//Bola maciza con paracaidas
	nombre="esfera-hierro-maciza-con-paracaidas";
	Cd = 0.8;
	m = 2+m1;
	A = pi*1.5*1.5;
	runge(t,x,y,y4,x4,m,Cd,A,nombre);
	
	//Bola hueca con paracaidas
	nombre="esfera-hierro-hueca-con-paracaidas";
	m = 2+m2;
	A = pi*1.5*1.5;
	runge(t,x,y,y4,x4,m,Cd,A,nombre);
	
	cout << "\n\nPresiona enter para salir ";
	cin.ignore();
	return 0;
}

//Planteamiento de la ecuacion diferencia como sistema de ecuaciones
double fx(double t, double x, double y)
{
	return (y);
}

double fy(double t, double x, double y, double m, double Cd, double A)
{
	return (g-1/(2*m)*Cd*Pf*A*y*y);
}

//Calculo de la masa
double masa(double p, double r)
{
	return(p*4/3*pi*r*r*r);
}

//Runge-kutta-Felberg
void runge(vector<double> &t, vector<double> &x, vector<double> &y, vector<double> &y4, vector<double> &x4, double m, double Cd, double A, string nombre)
{
	double kx1, kx2, kx3, kx4, kx5, kx6;
	double ky1, ky2, ky3, ky4, ky5, ky6;
	double errorx = 0;
	double errory = 0;
	signed i = 0;
	cout << "# Caso:  " << string(nombre) << "     Masa:  " << m<< endl;
	//cout << "Con: x5, x4 la posicion del objeto (soluciones de 5º y 4º orden)\n\n";
	//cout << "Con: y5, y4 la velocidad del objeto (soluciones de 5º y 4º orden)\n\n";
	cout  << right << setw(espacio) << "t" << right << setw(espacio) << "x5" << right << setw(espacio) << "x4" << right << setw(espacio) << "y5"
	      << right << setw(espacio) << "y4" << right << setw(espacio) << "errorx" << right << setw(espacio) << "errory\n\n";

	//Comienza la salida de datos
	ofstream datos;
	datos.open (string(nombre+".txt").c_str());
	//Inicio del calculo de RKF 4 y 5. Uso del while para salir cuando llegue al suelo el objeto
	while(x[i]<0) {
		//Redimensionamos los vectores para que quepa un elemento mas
		x.resize(x.size()+1);
		x4.resize(x4.size()+1);
		y.resize(y.size()+1);
		y4.resize(y4.size()+1);
		t.resize(t.size()+1);

		//Calculamos las k
		kx1=h*fx(t[i],x[i],y[i]);
		ky1=h*fy(t[i],x[i],y[i],m,Cd,A);
		kx2=h*fx(t[i]+h/2,x[i]+kx1/2,y[i]+ky1/2);
		ky2=h*fy(t[i]+h/2,x[i]+kx1/2,y[i]+ky1/2,m,Cd,A);
		kx3=h*fx(t[i]+3/8*h,x[i]+3*kx1/32+9*kx2/32,y[i]+3*ky1/32+9*ky2/32);
		ky3=h*fy(t[i]+3/8*h,x[i]+3*kx1/32+9*kx2/32,y[i]+3*ky1/32+9*ky2/32,m,Cd,A);
		kx4=h*fx(t[i]+12*h/13,x[i]+1932*kx1/2197-7200*kx2/2197+7296*kx3/2197,y[i]+1932*ky1/2197-7200*ky2/2197+7296*ky3/2197);
		ky4=h*fy(t[i]+12*h/13,x[i]+1932*kx1/2197-7200*kx2/2197+7296*kx3/2197,y[i]+1932*ky1/2197-7200*ky2/2197+7296*ky3/2197,m,Cd,A);
		kx5=h*fx(t[i]+h,x[i]+439*kx1/216-8*kx2+3680*kx3/513-845*kx4/4104,y[i]+439*ky1/216-8*ky2+3680*ky3/513-845*ky4/4104);
		ky5=h*fy(t[i]+h,x[i]+439*kx1/216-8*kx2+3680*kx3/513-845*kx4/4104,y[i]+439*ky1/216-8*ky2+3680*ky3/513-845*ky4/4104,m,Cd,A);
		kx6=h*fx(t[i]+h/2,x[i]-8*kx1/27+2*kx2-3544*kx3/2565+1859*kx4/4104-11*kx5/40,y[i]-8*ky1/27+2*ky2-3544*ky3/2565+1859*ky4/4104-11*ky5/40);
		ky6=h*fy(t[i]+h/2,x[i]-8*kx1/27+2*kx2-3544*kx3/2565+1859*kx4/4104-11*kx5/40,y[i]-8*ky1/27+2*ky2-3544*ky3/2565+1859*ky4/4104-11*ky5/40,m,Cd,A);

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
		//Descomentar la siguiente linea para ver todos los datos
		//cout  << right << setw(espacio) << t[i] << right << setw(espacio) << x[i] << right << setw(espacio) << x4[i] << right << setw(espacio) << y[i] << right << setw(espacio) << y4[i] << right << setw(espacio) << errorx << right << setw(espacio) << errory<< endl;
		
		//Damos un paso en el tiempo
		t[i+1]=t[0]+(i+1)*h;
		//Aumentamos contador
		i += 1;
	}
	datos.close();
	i=i-1; //Le restamos el ultimo valor a i para que sea el valor antes de llegar al suelo
	cout  << right << setw(espacio) << t[i] << right << setw(espacio) << x[i] << right << setw(espacio) << x4[i] << right << setw(espacio) << y[i] << right << setw(espacio) << y4[i] << right << setw(espacio) << errorx << right << setw(espacio) << errory<< "\n";
	cout << "-----------------------------------------------------------------------------------------------------------------------\n\n\n";
}
