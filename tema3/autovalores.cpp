/*
 *Se utiliza la librería Eigen para calcular las matrices inversas, multiplicaciones de matriz por vector y demás
 *
 */

#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

void diag(const Ref<Matrix3d> m, Ref<Vector3d> v, double& autovalor)
{
	double normaAntes, normaDespues;
	v << 1, 1, 1;
	do {
		normaAntes = v.norm();
		// Multiplicamos la matriz por el vector:
		v = m * v;
		// Vemos el valor máximo absoluto del vector
		autovalor = v.cwiseAbs().maxCoeff();
		if(v.maxCoeff() != autovalor) {
			autovalor = -autovalor;
		}
		// Lo normalizamos
		v = v / autovalor;
		normaDespues = v.norm();
	} while(abs(normaAntes - normaDespues) > 0.0000001);
}

int main()
{
	Matrix3d m;
	 m << 6, 2, 5, 2, 2, 3, 5, 3, 6;
	// m << 3, -1, 0, -2, 4, -3, 0, -1, 1;
	// m << 4, -1, 1, 1, 1, 1, -2, 0, -6;
	Vector3d v;
	double autovalor1 = 0; // Calculemos el primer autovalor
	diag(m, v, autovalor1); //Realizamos el proceso de encontrar el autovalor
	cout << "--> Autovalor:\n\n " << autovalor1;
	cout << "\n\n--> Autovector: \n\n" << v << "\n\n###############\n\n";

	// Calculemos el segundo autovalor
	Matrix3d minvertida = m.inverse(); //Usamos la matriz invertida
	Vector3d u;
	double autovalor2 = 0;
	diag(minvertida, u, autovalor2);
	autovalor2 = 1. / autovalor2; // Inversa del autovalor
	cout << "--> Autovalor:\n\n " << autovalor2;
	cout << "\n\n--> Autovector: \n\n" << u << "\n\n###############\n\n";

	// Calculemos el tercer autovalor usando la traza
	double autovalor3 = m.trace() - (autovalor1 + autovalor2);
	cout << "--> Autovalor:\n\n " << autovalor3;
	// Desplazamos la matriz
	Matrix3d mdesplazada;
	mdesplazada = m;
	for(int i = 0; i < 3; i++) {
		mdesplazada(i, i) = m(i, i) - autovalor3;
	}
	Vector3d w;
	Matrix3d mdesplazadaInvertida = mdesplazada.inverse(); //Invertimos la matriz desplazada
	diag(mdesplazadaInvertida, w, autovalor3);
	cout << "\n\n--> Autovector: \n\n" << w << "\n\n###############";
}
