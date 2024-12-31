#include <iostream>
#include <fstream>
#include <iomanip>
#include "Wavelet_Analysis.h"

int main()
{
	int N = 512;

	//Здесь менять номер этапа
	int Stage = 1;

	std::vector<std::complex<double>> Z(N), Koef_Psi, Koef_Fi;

	std::vector<std::complex<double>> P, Q, Z_Rec;

	//Параметры а и b из варианта
	for (int n = 0; n < N; n++)
	{
		if (n >= 128 && n <= 255) Z[n].real(sin(fabs(pow(n - 128, 1.86)) / 128.0));
		if (n >= 384 && n <= 447) Z[n].real(sin(fabs(pow(n - 128, 2.3)) / 128.0));
	}

	//Здесь менять базис. Complex_Shannon - Шеннон, Dobeshi - Добеши D6
	Com_Methods::Wavelet_Analysis Wavelet_Test(N, Com_Methods::Wavelet_Analysis::Basis_Type::Complex_Shannon);

	Wavelet_Test.Analysis_Phase(Stage, Z, Koef_Psi, Koef_Fi);

	Wavelet_Test.Synthesis_Phase(Stage, Koef_Psi, Koef_Fi, P, Q, Z_Rec);

	std::ofstream file("Recsignal.txt");
	for (int i = 0; i < N; i++)
	{
		file << Z_Rec[i].real() << "\n";
	}
	file.close();

	std::ofstream file1("Psi.txt");
	std::ofstream file2("Fi.txt");
	for (int i = 0; i < Koef_Fi.size(); i++){
		file1 << Koef_Psi[i].real() << "\n";
		file2 << Koef_Fi[i].real() << "\n";
	}
	file1.close();
	file2.close();

	std::vector<std::complex<double>> E(N);
	for (int i = 0; i < N; i++){
		E[i].real(fabs(Z[i].real() - Z_Rec[i].real()));
		E[i].imag(fabs(Z[i].imag() - Z_Rec[i].imag()));
	}

	std::ofstream file3("Err.txt");
	for (int i = 0; i < N; i++){
		file3 << E[i].real() << "\n";
	}
	file3.close();

	std::ofstream file4("Z.txt");
	for (int i = 0; i < N; i++){
		file4 << Z[i].real() << "\n";
	}
	file4.close();

}