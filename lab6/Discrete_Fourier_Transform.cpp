#include "Descrete_Fourier_Transform.h"

// прямое преобразование
// Data - входные данные, Result - массив результата
std::vector<std::complex<double>> DFT(const std::vector<std::complex<double>> &Data)
{
	int N = Data.size();
	std::vector<std::complex<double>> Result(N);
	std::complex<double> Exp;
	for (int m = 0; m < N; ++m)
	{
		Result[m].real(0);
		Result[m].imag(0);
		for (int n = 0; n < N; ++n)
		{
			Exp.real(cos(-2.0 * PI * m * n / N));
			Exp.imag(sin(-2.0 * PI * m * n / N));
			Result[m] += Data[n] * Exp;
		}
	}
	return Result;
}

std::vector<std::complex<double>> IDFT(const std::vector<std::complex<double>> &Data)
{
	int N = Data.size();
	std::vector<std::complex<double>> Result(N);
	std::complex<double> Exp;
	for (int m = 0; m < N; ++m)
	{
		Result[m].real(0);
		Result[m].imag(0);
		for (int n = 0; n < N; ++n)
		{
			Exp.real(cos(2.0 * PI * m * n / N));
			Exp.imag(sin(2.0 * PI * m * n / N));
			Result[m] += Data[n] * Exp;
		}
		Result[m] /= double(N);
	}
	return Result;
}

// быстрое преобразование (длина вектора - чётное число)
// Data - входные данные, Result - массив результата
std::vector<std::complex<double>> FFT(const std::vector<std::complex<double>> &Data)
{
	int N = Data.size(), M = N / 2;
	std::vector<std::complex<double>> Result(N);

	std::complex<double> Exp, U, V;

	for (int m = 0; m < M; m++)
	{
		U.real(0.0);
		U.imag(0.0);
		V.real(0.0);
		V.imag(0.0);
		for (int n = 0; n < M; n++)
		{
			Exp.real(cos(-2.0 * PI * m * n / M));
			Exp.imag(sin(-2.0 * PI * m * n / M));
			U += Data[2 * n] * Exp;
			V += Data[2 * n + 1] * Exp;
		}

		Exp.real(cos(-2.0 * PI * m / N));
		Exp.imag(sin(-2.0 * PI * m / N));
		Result[m] = U + Exp * V;
		Result[m + M] = U - Exp * V;
	}
	return Result;
}

//------------------------------------------------------------------------------------------

// обратное быстрое преобразование (длина вектора чётная)
// Data - входные данные, Result - массив результата
std::vector<std::complex<double>> IFFT(const std::vector<std::complex<double>> &Data)
{
	int N = Data.size();
	auto Result = FFT(Data);
	std::complex<double> Val;
	for (int i = 1; i <= N / 2; i++)
	{
		Val = Result[i];
		Result[i] = Result[N - i] / double(N);
		Result[N - i] = Val / double(N);
	}
	Result[0] /= double(N);
	return Result;
}