/*
 *
 * En este programa se incluyen los métodos: explicito, Crank-Nicolson y Tetha.
 *
 */

#include <stdio.h>
#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <string>

using namespace std;

// Prototipos
void mostrar(double** M, string nombre); //Muestra y guarda una matriz
void mostrarvector(vector<double> &v);
void tridiagonal(vector<double>& c, vector<double>& b, vector<double>& a, vector<double>& d, vector<double>& x);
void condiciones(double** M); //Si se le da una matriz, la devuelve con las condiciones de contorno e iniciales del problema
double modulo(double **M); //Calcula el modulo de una matriz

// Variables de la placa
double k = 0.13;
double cc = 0.11;
double p = 7.8;
double dx = 0.25;
double r = 0.5;
double tetha = 2. / 3.;
double tamX = 2.;
double dt = cc * p * dx * dx * r / (k);
double tolerancia = 0.01; // Tolerancia en la temperatura opcional para dejar de calcular

// Tamaño de la matriz
int sizeX = tamX / dx + 1; //Tamaño en el eje x
int sizeT = 15; //Intervalos de tiempo

int main(int argc, char** argv)
{
	//Comprobamos si el dx es valido
	if(floor(tamX / dx) != (tamX / dx)) {
		cout << "El valor de dx = " << dx << " no da un numero de pasos entero, elije otro valor";
		return 0;
	}

	// Crear matriz de temperaturas u para el metodo explicito
	double** u = new double* [sizeX];
	for(int i = 0; i < sizeX; ++i) {
		u[i] = new double[sizeT];
	}

	condiciones(u);
	//Fin de crear la matriz

	//Empezamos el método explícito
	for(int j = 1; j < sizeT ; j++) {
		for(int i = 1; i < sizeX - 1; i++) {
			u[i][j] = r * (u[i - 1][j - 1] + u[i + 1][j - 1]) + (1.0 - 2.0 * r) * u[i][j - 1];
		}
	}
	cout << "\n####################################   MÉTODO EXPLÍCITO   ####################################\n";
	string nombre = "explicito";
	mostrar(u, nombre);

	/*
	 *	Método de Cranck-Nicolson
	 *
	 */

	// Creamos otra matriz de temperaturas para no pisar la del método explicito
	double** v = new double* [sizeX];
	for(int i = 0; i < sizeX; ++i) {
		v[i] = new double[sizeT];
	}
	condiciones(v);
	// Fin de crear la matriz de temperaturas

	//Vector de soluciones y vectores tridiagonal
	vector<double> X (0); // vector de soluciones
	vector<double> a (0); //diagonal inferior
	vector<double> b (0); //diagonal principal
	vector<double> c (0); //diagonal superior
	vector<double> d (0); //terminos independientes

	//Se definen los valores de los vectores según la bibliografia: http://www.ugr.es/~lorente/APUNTESCN/capitulo5.pdf (página 6)
	//Vectores inciciales
	for(int i = 0; i < sizeX - 2; i++ ) {
		a.push_back(-r / 2);
		b.push_back(1. + r);
		c.push_back(-r / 2);
		X.push_back(0);
	}

	//Comienzo del calculo de cada elemento por el método de Crank-Nicolson
	for (int j = 0; j < sizeT - 1; j++) {
		d.resize(0);
		// Se determina el valor de cada elemento de d
		d.push_back((1 - r) * v[1][j] + (r / 2.) * v[2][j]);
		for(int i = 1; i < sizeX - 2; i++) {
			d.push_back((r / 2.) * v[i][j] + (1 - r) * v[i + 1][j] + (r / 2.) * v[i + 2][j]);
		}
		d.push_back((r / 2.) * v[sizeX - 3][j] + (1 - r) * v[sizeX - 2][j]);
		// Fin de determinar elementos de d
		tridiagonal(c, b, a, d, X);
		for(int i = 0; i < sizeX - 2; i++) {
			v[i + 1][j + 1] = X[i]; // Se inscribe la solución en la matriz de temperaturas
		}
	}
	cout << "\n####################################   MÉTODO DE CRANK-NICOLSON   ####################################\n";
	nombre = "crank-nicholson";
	mostrar(v, nombre);

	/*
	 *	Método Tetha
	 *
	 */

	// Creamos otra matriz de temperaturas para no pisar la del método explicito
	double** w = new double* [sizeX];
	for(int i = 0; i < sizeX; ++i) {
		w[i] = new double[sizeT];
	}
	condiciones(w);
	// Fin de crear la matriz de temperaturas

	//Vector de soluciones y vectores tridiagonal
	X.resize (0); // vector de soluciones
	a.resize (0); //diagonal inferior
	b.resize (0); //diagonal principal
	c.resize (0); //diagonal superior
	d.resize (0); //terminos independientes

	//Se definen los valores de los vectores según la bibliografia: http://www.ugr.es/~lorente/APUNTESCN/capitulo5.pdf (página 6)
	//Vectores inciciales
	for(int i = 0; i < sizeX - 2; i++ ) {
		a.push_back(-r * tetha);
		b.push_back(1. + r * tetha * 2.);
		c.push_back(-r * tetha);
		X.push_back(0);
	}

	//Condiciones de salida
	double moduloAnterior = 0;
	double moduloNuevo = 0;
	//Comienzo del calculo de cada elemento por el método de Crank-Nicolson
	for (int j = 0; j < sizeT - 1; j++) {
		d.resize(0);
		// Se determina el valor de cada elemento de d
		d.push_back((1 - r * (1 - tetha) * 2) * w[1][j] + (r * (1. - tetha)) * w[2][j]);
		for(int i = 1; i < sizeX - 2; i++) {
			d.push_back((r * (1. - tetha)) * w[i][j] + (1 - r * (1 - tetha) * 2) * w[i + 1][j] + (r * (1. - tetha)) * w[i + 2][j]);
		}
		d.push_back((r * (1. - tetha)) * w[sizeX - 3][j] + (1 - r * (1 - tetha) * 2) * w[sizeX - 2][j]);
		// Fin de determinar elementos de d
		moduloAnterior = modulo(w);
		tridiagonal(c, b, a, d, X);
		for(int i = 0; i < sizeX - 2; i++) {
			w[i + 1][j + 1] = X[i]; // Se inscribe la solución en la matriz de temperaturas
		}
		moduloNuevo = modulo(w);
		/*if(abs(moduloNuevo - moduloAnterior) < tolerancia) {
			break;		//Descomentar para usar la tolerancia como condicion de salida
		}*/
	}
	cout << "\n####################################   MÉTODO THETA   ####################################\n";
	nombre = "tetha";
	mostrar(w, nombre);
	cout << "Con:\n\nr = " << r << "\ntheta = "  << tetha << "\n\n";
	return 0;
}

void mostrar(double** M, string nombre)
{
	cout << "\nTemperaturas en el tiempo:\n\n";
	ofstream datos;
	datos.open (string(nombre + ".txt").c_str());
	cout << setw(15) << "x = ";
	for (int j = 0; j < sizeX; j++) {
		cout << setw(11) << right << setprecision(4) << dx*j;
	}
	cout << endl << endl;
	for(int i = 0; i < sizeT; ++i) {
		cout << setw(10) << right << "t = " << setw(6) << right << dt*i ;
		for(int j = 0; j < sizeX; j++) {
			cout << setw(10) << right << M[j][i] << " ";
			datos << M[j][i] << "	";
		}
		cout << endl;
		datos << endl;
	}
	datos.close();
	cout << endl;
	cout << endl;
}

void condiciones(double** M)
{
	// Condiciones iniciales del problema
	for(int i = 0; i < sizeX; i++) {
		double xx = i * dx;
		if(i * dx < tamX / 2) {
			M[i][0] = 100.0 * xx;
		} else {
			M[i][0] = 200.0 - 100.0 * xx;
		}
	}

	// Condiciones de contorno a cero grados
	for(int j = 0; j < sizeT; j++) {
		M[0][j] = 0.;
		M[sizeX - 1][j] = 0.;
	}
}

void tridiagonal(vector<double>& c, vector<double>& b, vector<double>& a, vector<double>& d, vector<double>& x)
{
	/*
	 * c es la diagonal superior
	 * b es la diagonal principal
	 * a es la diagonal inferior
	 * d es el vector independiente
	 */

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

void mostrarvector(vector<double> &v)
{
	for(unsigned i = 0; i < v.size(); i++) {
		cout << " i=" << i << ": " << v[i] << " \n";
	}
}

double modulo(double **M)
{
	double suma = 0;
	for(int i = 0; i < sizeX; i++) {
		for(int j = 0; j < sizeT; j++) {
			suma += M[i][j] * M[i][j];
		}
	}
	return sqrt(suma);
}
