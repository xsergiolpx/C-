#include <stdio.h>
#include <iostream>
#include <vector>
#include <iomanip>

using namespace std; 
float h=0.1;
double x0=0; //punto inicial

double f(double x, double y){
	return (-2*x-y);
}

void euler(vector<double> &x, vector<double> &y, vector<double> &dy){
	//Metodo de Euler normal
	cout <<  "################      Metodo de Euler      ################\n\n";
	cout << "    Xn      Yn      Y'n      hY'n\n\n";
	int espacio=10;
		for (unsigned i=0;i<y.size()-1;i++){
		x[i]=x0+h*i;
		dy[i]=f(x[i],y[i]);
		y[i+1]=y[i]+h*dy[i];
		//Salida por pantalla
		cout << setprecision(5) << right << setw(4) << x[i] << right << setw(espacio) << y[i] << right << setw(espacio) << dy[i]<< right << setw(espacio) << h*dy[i] << endl;
	}
}

void eulerMejorado(vector<double> &x2, vector<double> &y2, vector<double> &dy2,vector<double> &y2Euler,vector<double> &dyAv){
	//Mostrar las columnas en pantalla
	int espacio=11;
	cout << "\n################      Metodo de Euler mejorado      ################\n\n";
	cout <<  right << setw(4) << "Xn"<< right << setw(espacio) <<
	"Yn"<< right << setw(espacio) <<  "hY'n" << right << setw(espacio) <<
	"Yn+1"<< right << setw(espacio) << "hY'n+1,p"<< right << setw(espacio) <<
	"hY'av"<< right << setw(espacio) << "Yn+1,c\n\n";
	
	for (unsigned i=0;i<y2.size()-1;i++){
		x2[i]=x0+h*i;
		x2[i+1]=x0+h*(i+1); // Punto siguiente para calcular Euler normal
		dy2[i]=f(x2[i],y2[i]);
		y2Euler[i+1]=y2[i]+h*dy2[i]; //Valor de la funcion en el punto siguiente con Euler normal
		dy2[i+1]=-2*x2[i+1]-y2Euler[i+1]; //Valor de la derivada en el punto siguiente
		dyAv[i]=(dy2[i]+dy2[i+1])/2; //Valor medio
		y2[i+1]=y2[i]+h*dyAv[i];
		//Salida por pantalla
		cout << setprecision(4) << right << setw(4) << x2[i] << right  << setw(espacio) << y2[i]
		<< right << setw(espacio) << h*dy2[i] << right << setw(espacio) << y2[i+1]
		<< right << setw(espacio) << h*dy2[i+1] << right << setw(espacio) << h*dyAv[i] 
		<< right << setw(espacio) << y2[i+1] << endl;
	}
}

int main(int argc, char **argv)
{	
	//resolver y'=-2x-y;   y(0)=-1
	int tam=6; //Tamaño del vector
	//Vectores
	vector<double> y (tam);
	vector<double> dy (tam);
	vector<double> x (tam);
	y[0]=-1; //Condicion de contorno
	euler(x,y,dy);

	
	//Euler mejorado
	tam=7; //Mostramos un valor mas
	vector<double> y2 (tam);
	vector<double> dy2 (tam);
	vector<double> x2 (tam);
	vector<double> y2Euler (tam);
	vector<double> dyAv (tam);
	y2[0]=-1;
	
	eulerMejorado(x2,y2,dy2,y2Euler,dyAv);
	
	cout << "\n\nPresiona enter para continuar\n\n";
	cin.ignore();
	return 0;
}
