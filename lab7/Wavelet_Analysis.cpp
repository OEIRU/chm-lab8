#include "Wavelet_Analysis.h"

namespace Com_Methods
{
	//�����������: ���������� �������� ���������� � ������������ ���������
	//Count_Data - ����� �������� � ������� �������, Type - ��� �������-������
	Wavelet_Analysis::Wavelet_Analysis(int Count_Data, Basis_Type Type)
	{
		//����� �������� � ������� �������
		int N = Count_Data;
		U.clear(); U.resize(N);
		V.clear(); V.resize(N);

		switch (Type)
		{
			//����������� ����� �������
			case Basis_Type::Complex_Shannon:
			{
				if (N % 4 != 0) throw std::exception();
					
				std::complex<double> Val;
					
				U[0] = V[0] = 1.0 / sqrt(2.0);
					
				for (int i = 1; i < N; i++)
				{
					Val.real(sqrt(2.0) / N * cos(PI * i / N) * sin(PI * i / 2.0) / sin(PI * i / N));
					Val.imag(-sqrt(2.0) / N * sin(PI * i / N) * sin(PI * i / 2.0) / sin(PI * i / N));
					U[i] = Val;
					V[i] = pow(-1, i) * Val;
				}
				break;
			}

			case Basis_Type::Dobeshi:
			{
				double a = 1 - sqrt(10), b = 1 + sqrt(10), c = sqrt(5 + 2 * sqrt(10));

				// U[0] = (sqrt(2) / 32) * (b + c);
				// U[1] = (sqrt(2) / 32) * (2 * a + 3 * b + 3 * c);
				// U[2] = (sqrt(2) / 32) * (6 * a + 4 * b + 2 * c);
				// U[3] = (sqrt(2) / 32) * (6 * a + 4 * b - 2 * c);
				// U[4] = (sqrt(2) / 32) * (2 * a + 3 * b - 3 * c);
				// U[5] = (sqrt(2) / 32) * (b - c);

				U[0] = 0.47046721;
				U[1] = 1.14111692;
				U[2] = 0.650365;
				U[3] = -0.19093442;
				U[4] = -0.12083221;
				U[5] = 0.0498175;

				for (int i = 6; i < N; i++)
					U[i] = 0;

				V[0] = -U[1];
				V[1] = U[0];

				for (int i = 2; i < N - 4; i++)
					V[i] = 0;
				
				V[N - 4] = -U[5];
				V[N - 3] = U[4];
				V[N - 2] = -U[3];
				V[N - 1] = U[2];
			}
		}
	}
	
	//--------------------------------------------------------------------------------------------------------------------

	//���������� ������� �������-��������
	//Stages - ���������� ������ 
	void Wavelet_Analysis::Wavelet_Filters_System(int Stages)
	{
		//��������� �������
		Operators Operator;

		//����� �������� � ������� �������
		int N = U.size();

		//�������-�������
		std::vector<std::vector<std::complex<double>>> U_Filters(Stages), V_Filters(Stages);

		//������� 1-�����
		U_Filters[0] = U;
		V_Filters[0] = V;

		//��������� �������
		for (int i = 1; i < Stages; i++)
		{
			//����� ��������� � �������� i-�����
			int Count_Elements = N / int(pow(2, i));

			U_Filters[i].resize(Count_Elements);
			V_Filters[i].resize(Count_Elements);

			//���������� �������� i-�����
			for (int n = 0; n < Count_Elements; n++)
			{
				int Max_Index = int(pow(2, i));
				for (int k = 0; k < Max_Index; k++)
				{
					U_Filters[i][n] += U_Filters[0][n + k * N / Max_Index];
					V_Filters[i][n] += V_Filters[0][n + k * N / Max_Index];
				}
			}
		}

		//������ DFT ��� �������� ������ "*"
		Com_Methods::Discrete_Fourier_Transformation DFT;
		
		//������� ��� �������� �������������� ���������� ���������� ������
		std::vector<std::complex<double>> U_u, U_v;

		//��������� ������� f � g
		f.resize(Stages); g.resize(Stages);
		f[0] = V_Filters[0];
		g[0] = U_Filters[0];
		for (int i = 1; i < Stages; i++)
		{
			//���������� ������� U_i
			Operator.Upsampling_Operator(i, U_Filters[i], U_u);
			//���������� ������� V_i
			Operator.Upsampling_Operator(i, V_Filters[i], U_v);
			//������
			DFT.Convolution(g[i - 1], U_v, f[i]);
			DFT.Convolution(g[i - 1], U_u, g[i]);
		}
	}
	
	//--------------------------------------------------------------------------------------------------------------------

	//���������� �������-������
	//Stage - ����� �����
	//Psi - ��������������� ����������
	//Fi  - �������������� ����������
	void Wavelet_Analysis::Wavelet_Basis(int Stage,
										 std::vector<std::vector<std::complex<double>>>& Psi,
										 std::vector<std::vector<std::complex<double>>>& Fi)
	{
		//��������� ������ � ���������� ������������
		Operators Operator;

		//����� ��������
		int Count_Data = U.size();

		//����� ��������� � �������-�����e Stage-����� 
		int Count_Bas_Elements = Count_Data / int(pow(2, Stage));

		//�������� �� ������������� ������� �������� Stage-�����
		if (g.size() < Stage) Wavelet_Filters_System(Stage + 1);

		//��������� �������� �������-������
		Psi.resize(Count_Bas_Elements); Fi.resize(Count_Bas_Elements);
			
		for (int i = 0; i < Count_Bas_Elements; i++)
		{
			int Index = int(pow(2, Stage)) * i;
			Operator.Cyclic_Shift(Index, f[Stage - 1], Psi[i]);
			Operator.Cyclic_Shift(Index, g[Stage - 1], Fi[i]);
		}
	}

	//--------------------------------------------------------------------------------------------------------------------

	//���� ������� �������
	//Stage - ����� �����
	//Data  - ������ �������
	//Koef_Psi - ������������ ���������� �� ��������������� ����������
	//Koef_Fi  - ������������ ���������� �� �������������� ����������
	void  Wavelet_Analysis::Analysis_Phase(int Stage,
										   const std::vector<std::complex<double>>& Data,
										   std::vector<std::complex<double>>& Koef_Psi,
										   std::vector<std::complex<double>>& Koef_Fi)
	{
		//�������� ���������� ������������
		Operators Operator;

		//�������-�����
		std::vector<std::vector<std::complex<double>>> Psi, Fi;

		//�������� �������-������
		Wavelet_Basis(Stage, Psi, Fi);

		//���������� �������� ���������
		int Count_Bas_Elements = Psi.size();

		//���� �������
		for (int Bas_Element_Index = 0; Bas_Element_Index < Count_Bas_Elements; Bas_Element_Index++)
		{
			//���� �������: ��������� ������������ ������� ������� �� �������� �������
			Koef_Psi.push_back(Operator.Dot_Product(Data, Psi[Bas_Element_Index]));
			Koef_Fi.push_back (Operator.Dot_Product(Data, Fi [Bas_Element_Index]));
		}
	}

	//--------------------------------------------------------------------------------------------------------------------

	//���� ������� (�������������� �������)
	//Stage - ����� �����
	//Koef_Psi  - ������������ ���������� �� ��������������� ����������
	//Koef_Fi   - ������������ ���������� �� �������������� ����������
	//P � Q - ������� ������������� �������� ������� �� ��������� {Psi} � {Fi} ����� Stage
	//Recovery  - ��������������� ������ �� ����� Stage
	void  Wavelet_Analysis::Synthesis_Phase(int Stage,
											const std::vector<std::complex<double>>& Koef_Psi,
											const std::vector<std::complex<double>>& Koef_Fi,
											std::vector<std::complex<double>>& P,
											std::vector<std::complex<double>>& Q,
											std::vector<std::complex<double>>& Recovery)
	{
		//�������-�����
		std::vector<std::vector<std::complex<double>>> Psi, Fi;

		//�������� �������-������
		Wavelet_Basis(Stage, Psi, Fi);

		//���������� �������� ���������
		int Count_Bas_Elements = Psi.size();

		//���������� �������� �������
		int Count_Data = U.size();

		//���� �������: ���������� ��������� ������ ���� ���������� ������� �� ������ 
		for (int Data_Index = 0; Data_Index < Count_Data; Data_Index++)
		{
			//���������� ��� ���������� ���� P � Q
			std::complex<double> P_Stage(0, 0), Q_Stage(0, 0);

			for (int Bas_Element_Index = 0; Bas_Element_Index < Count_Bas_Elements; Bas_Element_Index++)
			{
				//��������� �������������� �� Stage-�����
				P_Stage += Koef_Fi[Bas_Element_Index]  * Fi[Bas_Element_Index][Data_Index];
				//�������������� ���������� � �������������� �� (Stage + 1)-�����
				Q_Stage += Koef_Psi[Bas_Element_Index] * Psi[Bas_Element_Index][Data_Index];

			}
			P.push_back(P_Stage);
			Q.push_back(Q_Stage);
			Recovery.push_back(P[Data_Index] + Q[Data_Index]);
		}
	}
}
