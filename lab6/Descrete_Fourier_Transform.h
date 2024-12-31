#pragma once
#include <iostream>
#include <iomanip>
#include <vector>
#include <complex>

#define PI 3.1415926535897932

std::vector<std::complex<double>> FFT(const std::vector<std::complex<double>>& Data);
std::vector<std::complex<double>> IFFT(const std::vector<std::complex<double>>& Data);
std::vector<std::complex<double>> DFT(const std::vector<std::complex<double>>& Data);
std::vector<std::complex<double>> IDFT(const std::vector<std::complex<double>>& Data);