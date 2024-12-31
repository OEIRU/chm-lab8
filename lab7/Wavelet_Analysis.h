#pragma once

#ifndef Wavelet_Analysis_H
#define Wavelet_Analysis_H

#include <vector>
#include <complex>
#include "Operators.h"
#include "Discrete_Fourier_Transformation.h"

namespace Com_Methods
{
	class Wavelet_Analysis
	{
	private:
		//��������� � ����������� ��������
		std::vector<std::complex<double>> U, V;
		//������� f � g, ������ ������� ��������� p-������� �������-������
		std::vector<std::vector<std::complex<double>>> f, g;

	public:
		//���� �������
		enum class Basis_Type
		{
			//����������� ����� �������
			Complex_Shannon = 1,
			Dobeshi = 2,
		};

		//�����������: ���������� �������� ���������� � ������������ ���������
		//Count_Data - ����� �������� � ������� �������, Type - ��� �������-������
		Wavelet_Analysis(int Count_Data, Basis_Type Type);

	private:
		//����� ��� ���������� ������� �������-��������
		//Stages - ���������� ������
		void Wavelet_Filters_System(int Stages);
	
		//����� ��� ���������� �������-������
		//Stage - ����� �����
		//Psi - ��������������� ����������
		//Fi  - �������������� ����������
		void Wavelet_Basis(int Stage,
						   std::vector<std::vector<std::complex<double>>>& Psi,
						   std::vector<std::vector<std::complex<double>>>& Fi);

	public:
		//���� ������� �������
		//Stage - ����� �����
		//Data  - ������ �������
		//Koef_Psi - ������������ ���������� �� ��������������� ����������
		//Koef_Fi  - ������������ ���������� �� �������������� ����������
		void Analysis_Phase(int Stage, 
							const std::vector<std::complex<double>>& Data,
							std::vector<std::complex<double>>& Koef_Psi,
							std::vector<std::complex<double>>& Koef_Fi);

		//���� ������� (�������������� �������)
		//Stage - ����� �����
		//Koef_Psi  - ������������ ���������� �� ��������������� ����������
		//Koef_Fi   - ������������ ���������� �� �������������� ����������
		//P � Q - ������� ������������� �������� ������� �� ��������� {Psi} � {Fi} ����� Stage
		//Recovery  - ��������������� ������ �� ����� Stage
		void Synthesis_Phase(int Stage,
								const std::vector<std::complex<double>>& Koef_Psi,
								const std::vector<std::complex<double>>& Koef_Fi,
								std::vector<std::complex<double>>& P,
								std::vector<std::complex<double>>& Q,
								std::vector<std::complex<double>>& Recovery);
	};
}

#endif
