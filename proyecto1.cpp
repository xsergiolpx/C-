#include <stdio.h>
#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <string>
#include <ctime>

using namespace std;

//Variables de la placa
double largo = 5;
double ancho = 9;
double h = 1;

//Constantes del problema
double k = 0.16;
double d = 0.5;
double Q = 0.6;
double H = 0.073;
double ur = 25;

//Elementos
int sizeX = largo / h + 1;
int sizeY = ancho / h + 1;
double tolerancia = 0.000001;

//Prototipos
void mostrar(double **M); //Muestra la matriz de temperaturas de la placa y escribe los datos en un fichero
double modulo(double **M); //Calcula el modulo de la matriz
void mostrarvector(vector<double> &v); //Muestra un vector
void tridiagonal(vector<double>& c, vector<double>& b, vector<double>& a, vector<double>& d, vector<double>& x); //Resolucion del sistema de ecuaciones tridiagonal
void calcularMatriz(double **M, vector<double>& X); //Calcula la matriz de las temperaturas a partir del vector proporcionado

int main(int argc, char **argv)
{
	//Tiempo de calculo
	const clock_t begin_time = clock();
	//Comprobamos si el valor de h vale
	if(floor((largo / h) != ((largo ) / h) || floor((ancho ) / h) != (ancho ) / h))  {
		cout << "El valor de h = " << h << " no da un numero de pasos entero, elije otro valor";
		return 0;
	}
//Crear matriz de temperaturas
	double** matriz = new double*[sizeX];
	for(int i = 0; i < sizeX; ++i) {
		matriz[i] = new double[sizeY];
	}

	//Condiciones de contorno
	for (int i = 0; i < sizeX; i++) {
		matriz[i][0] = 20;
		matriz[i][sizeY - 1] = 20;
	}

	//Vector de soluciones y vectores tridiagonal
	vector<double> X (0);
	vector<double> a (0);
	vector<double> b (0);
	vector<double> c (0);

	//Vectores de la matriz de coeficientes
	for(int i = 0; i < (sizeY - 2) * 2 + (sizeX - 2) * (sizeY - 2); i++) {
		a.push_back(0);
		b.push_back(1);
		c.push_back(0);
		X.push_back(0);
	}
	//Comenzamos a resolver el sistema
	double moduloAnterior = 0;
	double moduloNuevo = 0;
	int bucles = 0;
	do {
		calcularMatriz(matriz, X);
		moduloAnterior = modulo(matriz);

		//Vector independiente
		vector<double> B (0);
		for(int j = 1; j < sizeY - 1; j++) {
			B.push_back(1 / (2. * h * H / k + 4.) * (matriz[0][j - 1] + matriz[0][j + 1] + 2. * matriz[1][j] + (50. * h * H * d - h * h * Q) / (k * d)));
		}

		for (int i = 1; i < sizeX - 1; i++) {
			for(int j = 1; j < sizeY - 1; j++) {
				B.push_back(1. / 4. * (matriz[i - 1][j] + matriz[i][j - 1] + matriz[i][j + 1] + matriz[i + 1][j] - h * h * Q / (k * d)));
			}
		}

		for(int j = 1; j < sizeY - 1; j++) {
			B.push_back(1. / 4.*(2.*matriz[sizeX - 2][j] + matriz[sizeX - 1][j - 1] + matriz[sizeX - 1][j + 1] + 30.*h - h * h * Q / (k * d)));
		}

		tridiagonal(c, b, a, B, X);
		calcularMatriz(matriz, X);
		moduloNuevo = modulo(matriz);
		bucles += 1;
	} while(abs(moduloNuevo - moduloAnterior) > tolerancia); //Condicion de salida
	mostrar(matriz);
	//mostrarvector(X);
	cout << "Sistema de " << (sizeY - 2) * 2 + (sizeX - 2) * (sizeY - 2) << " ecuaciones" << endl;
	cout << "Numero de iteraciones: " << bucles;
	cout << endl << "Tiempo de ejecución: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << " segundos\n";
	return 0;
}

//Funciones auxiliares

void mostrar(double **M)
{
	cout << "\nDistribucion de temperaturas en la placa:\n\n";
	ofstream datos;
	datos.open ("datos.txt");
	for (int i = 0; i < sizeX; ++i) {
		for (int j = 0; j < sizeY; j++) {
			cout << setw(10) << right << M[i][j] << " ";
			datos << M[i][j] << "	";
		}
		cout << endl;
		datos << endl;
	}
	datos.close();
	cout << endl;
	cout << endl;

	//Ver el punto mas frio
	double puntoFrio = M[0][0];
	double puntoCaliente = M[0][0];
	int x = 0, x2 = 0;
	int y = 0, y2 = 0;
	for (int i = 0; i < sizeX; ++i) {
		for (int j = 0; j < sizeY ; j++) {
			if(M[i][j] < puntoFrio) {
				puntoFrio = M[i][j];
				x = i;
				y = j;
			}
			if(M[i][j] > puntoCaliente) {
				puntoCaliente = M[i][j];
				x2 = i;
				y2 = j;
			}
		}
	}
	cout << "Punto más frío: " << puntoFrio << endl;
	cout << "Coordenadas (i, j): " << "(" << x << ", " << y << ") y (" << x << ", " << sizeY - 1 - y << ") \n\n";
	cout << "Punto más caliente: " << puntoCaliente << endl;
	cout << "Coordenadas (i, j): " << "(" << x2 << ", " << y2 << ") y (" << x2 << ", " << sizeY - 1 - y2 << ") \n\n";

}

double modulo(double **M)
{
	double suma = 0;
	for(int i = 0; i < sizeX; i++) {
		for(int j = 0; j < sizeY; j++) {
			suma += M[i][j] * M[i][j];
		}
	}
	return sqrt(suma);
}

void mostrarvector(vector<double> &v)
{
	for(unsigned i = 0; i < v.size(); i++) {
		cout << " i=" << i << ": " << v[i] << " \n";
	}
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

void calcularMatriz(double **matriz, vector<double>& X )
{

	//Calcular matriz
	int ultimaJ = 0;
	for(int j = 1; j < sizeY - 1; j++) {
		matriz[0][j] = X[j - 1];
		ultimaJ = j;
	}
	int ultimaJ2 = 0 + ultimaJ;
	for(int i = 1; i < sizeX - 1; i++) {
		for(int j = 1; j < sizeY - 1; j++) {
			matriz[i][j] = X[(j - 1 + ultimaJ * i)];
			ultimaJ2 = (j + ultimaJ * i);
		}
	}
	for(int j = 1; j < sizeY - 1; j++) {
		matriz[sizeX - 1][j] = X[j - 1 + ultimaJ2];
	}
}
