// csv读取.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <fstream>
#include <stdio.h>
#include <iostream>
#include<string>
#include<string.h>
#include<stdlib.h>
#include <vector>
#include <math.h>
using namespace std;
#define pi 3.1415926535 

const int y = 5;																//一元y-1次方程

int dft(char filename[], char filename2[], char filename3[], char filename4[], char filename5[], char filename6[], char filename7[], char filename8[], char filename9[], char filename10[], char filename11[], char filename12[])	//傅里叶转换
{

	/**************************************************************************************
	模块说明:
	读出所有数据vector<string> date1;过滤掉逗号压入vector<string> date2;
	*****************************************************************************************/

	fstream file;                                     //【1】声明一个文件输入输出流对象
	file.open(filename, std::ios::in);
	string strtmp;
	int date1num = 0, date2num = 0;
	vector<string> date1;
	vector<string> date2;
	char b[10000];
	while (getline(file, strtmp))
	{
		date1num++;
		//std::cout << strtmp << std::endl;
		date1.push_back(strtmp);
		strcpy(b, strtmp.c_str());
		char*p;
		p = strtok(b, ",");
		while (p != NULL)
		{
			date2num++;
			//std::cout << p << std::endl;
			date2.push_back(p);
			p = strtok(NULL, ",");
		}
	}
	file.close();


	int rows, cols;
	cout << date1num << endl;
	rows = date1.size();
	cout << "zonghangshu " << rows << endl;
	cout << date2num << endl;
	cout << "zongshujushu " << date2.size() << endl;
	cols = date2.size() / date1.size();
	cout << "zonglieshu " << cols << endl;


	/***************************************************************************************
	模块说明:
	将读入的CSV表格(字符串)格式-----转化为-----可以用于计算的double型表格
	****************************************************************************************/


	vector<double> date;
	for (int i = 0; i < date2num;)
	{
		double g = atof(date2[i].c_str());
		date.push_back(g);
		i++;
	}
	cout << "zhuanhuanhoudezongshuju" << date.size() << endl;




	/**************************************************************************************
	模块说明:
	对每一列的数据进行傅里叶转换求得实部R和虚部I
	*****************************************************************************************/
	ofstream oFile2, oFile3;
	//打开要输出的文件
	oFile2.open(filename2, ios::out | ios::trunc);
	oFile3.open(filename3, ios::out | ios::trunc);

	//int t = cols / 2;
	vector<double> R;
	vector<double> I;
	vector<double> R0;
	
	for (int j = 0; j < cols; j++)  /*第j列*/
	{

		for (int k = 0; k < rows; k++)
		{
			double Re = 0.0, Im = 0.0;
			for (int z = 0; z < rows; z++)  /*第z行*/
			{
				int n = cols * z + j;
				Re += (date[n])* cos(2 * pi*k*z / rows) / rows;
				Im += -(date[n])* sin(2 * pi*k*z / rows) / rows;

			}
			oFile2 << Re << ",";
			oFile3 << Im << ",";


			if (k == 40)
			{
				R.push_back(Re);
				I.push_back(Im);
				cout << "k=40" << "," << Re << "," << Im << endl;

			}
			if (k == 0)
			{
				R0.push_back(Re);
				
				cout << "k=0" << "," << Re << endl;

			}

		}
		oFile2 << endl;
		oFile3 << endl;

	}
	oFile2.close();
	oFile3.close();


	ofstream oFile8;
	oFile8.open(filename8, ios::out | ios::trunc);

	vector<double> ak;
	vector<double> bk;
	vector<double> ak0;
	oFile8 << " ak" << "," << " bk " <<","<<"ak0"<< endl;
	for (int i = 0; i < cols; i++)
	{
		ak.push_back( 2*R[i]);
		bk.push_back(2*I[i]);
		ak0.push_back(2*R0[i]);
		oFile8 << R[i] << "," << I[i] <<","<<R0[i]<< endl;
	}



	/***************************************************************************************
	模块说明:
	高斯消元法曲线拟合
	****************************************************************************************/



	/***************************************************************************************
	模块说明:
	求解  实部   方程的系数xR[n]
	****************************************************************************************/


	
	vector<double> DbR;												//实部
	vector<double> bbR;												//实部

	ofstream oFile4;												//打开拟合曲线的导数和函数值要输出的文件
	oFile4.open(filename4, ios::out | ios::trunc);

	double rr[105];
	for (int i = 0; i < 105; i++)
	{
		rr[i] =-0.024923 + i*0.049846 / 104;
		cout << rr[i] << ",";
	}
	cout << endl;
	
	for (int t = 0; t < 0.5; t++)
	{
		double xR[y];												//实部曲线拟合方程系数
		double s[y];												//虚部右侧
		double m[y][y];												//虚部左侧
		//double rr = (t * 0.2 + 0.04)*0.025;							//曲线拟合时的起点的r值
		for (int i = 0; i < y; i++)									//求解m[][] 行列式左侧 
		{
			for (int j = 0; j < y; j++)
			{

				double H = 0;
				for (int z = 10; z < 50; z++)
				{
					H += pow(rr[z], i + j);
				}
				m[i][j] = H;
				cout << "m[" << i << "][" << j << "]" << "=" << m[i][j]<<",";
			}
			cout << endl;
		}
	
		for (int i = 0; i < y; i++)								// 实部行列式右侧
		{
			double H = 0;
			for (int z = 10; z < 50; z++)
			{
				H += ak[z] * pow(rr[z], i);
			}
			s[i] = H;
			cout << "s[" << i << "]" << "=" << s[i] << ",";
		}

		double c[y];									 //存储初等行变换的系数，用于行的相减
		int i;
		for (int k = 0; k < y - 1; k++)					//消元的整个过程如下，总共n-1次消元过程
		{
			for (i = k + 1; i < y; i++)					//求出第K次初等行变换的系数
				c[i] = m[i][k] / m[k][k];

			for (i = k + 1; i < y; i++)					//第K次的消元计算
			{
				for (int j = 0; j < y; j++)
				{
					m[i][j] = m[i][j] - c[i] * m[k][j];
				}
				s[i] = s[i] - c[i] * s[k];
			}
		}

		xR[y - 1] = s[y - 1] / m[y - 1][y - 1];							//先计算出最后一个未知数

		for (int i = y - 2; i >= 0; i--)									//求出每个未知数的值
		{
			double sum = 0;
			for (int j = i + 1; j < y; j++)
			{
				sum += m[i][j] * xR[j];
			}
			xR[i] = (s[i] - sum) / m[i][i];				//xR[i]方程系数
		}

		cout << " the solution of the equations is:" << endl;
		cout << endl;
		for (i = 0; i < y; i++)
			cout << "shibuxishu a" << i << "=" << xR[i] << endl;

		for (int z = 10; z < 50; z++)			//实部
		{
			double DB = 0, b1 = 0;
			for (i = 1; i < y; i++)
			{
				DB += pow(i, 2)* xR[i] * pow(rr[z], i - 2);
			}
			//cout << "DbR" << r << "=" << DB << endl;
			DbR.push_back(DB);

			for (i = 0; i < y; i++)
			{
				b1 += xR[i] * pow(rr[z], i);
			}

			bbR.push_back(b1);
			//cout << "bbR" << r << "=" << b1 << endl;
			oFile4 << rr[z] << "," << ak[z] << "," << DB << "," << b1 << endl;
		}
	}
	oFile4.close();

	/***************************************************************************************
	模块说明:
	求解虚部方程的系数xR[n],DbR[],bbR[]
	****************************************************************************************/

	vector<double> DbI;												//虚部
	vector<double> bbI;												//虚部

	ofstream oFile5;												//打开拟合曲线的导数和函数值要输出的文件
	oFile5.open(filename5, ios::out | ios::trunc);

	for (int t = 0; t < 0.5; t++)
	{
		double xI[y];												//虚部曲线拟合方程系数
	double s[y];												//虚部右侧
	double m[y][y];												//虚部左侧
	//double rr = t * 0.047962 / 108 * 5 - 0.022649;			//曲线拟合时的起点的r值
	for (int i = 0; i < y; i++)									//求解m[][] 行列式左侧 
	{
		for (int j = 0; j < y; j++)
		{
			double H = 0;
			for (int z = 10; z < 50; z++)
			{
				H += pow(rr[z], i + j);
			}
			m[i][j] = H;
			cout << "m[" << i << "][" << j << "]" << "=" << m[i][j]<<",";
		}
		cout << endl;
	}
	
	for (int i = 0; i < y; i++)									// 虚部行列式右侧
	{
		double H = 0;
		for (int z = 10; z < 50; z++)
		{
			H += bk[z] * pow(rr[z], i);
		}
		s[i] = H;
		cout << "s[" << i << "]"  << "=" << s[i]<<",";
	}

	double c[y];											 //存储初等行变换的系数，用于行的相减
	int i;
	for (int k = 0; k < y - 1; k++)							//消元的整个过程如下，总共n-1次消元过程
	{
		for (i = k + 1; i < y; i++)							//求出第K次初等行变换的系数
			c[i] = m[i][k] / m[k][k];

		for (i = k + 1; i < y; i++)							//第K次的消元计算
		{
			for (int j = 0; j < y; j++)
			{
				m[i][j] = m[i][j] - c[i] * m[k][j];
			}
			s[i] = s[i] - c[i] * s[k];
		}
	}

	xI[y - 1] = s[y - 1] / m[y - 1][y - 1];							//先计算出最后一个未知数

	for (int i = y - 2; i >= 0; i--)									//求出每个未知数的值
	{
		double sum = 0;
		for (int j = i + 1; j < y; j++)
		{
			sum += m[i][j] * xI[j];
		}
		xI[i] = (s[i] - sum) / m[i][i];							//xI[i]方程系数
	}

	cout << " the solution of the equations is:" << endl;
	cout << endl;
	for (i = 0; i < y; i++)
		cout << "xubuxishu a" << i << "=" << xI[i] << endl;


	for (int z = 10; z < 50; z++)									//虚部
	{
		double DB = 0, b1 = 0;
		for (i = 1; i < y; i++)
		{
			DB += pow(i, 2)* xI[i] * pow(rr[z], i - 2);
		}

		DbI.push_back(DB);

		for (i = 0; i < y; i++)
		{
			b1 += xI[i] * pow(rr[z], i);
		}

		bbI.push_back(b1);

		oFile5 <<rr[z] << "," <<bk[z]<<","<< DB << "," << b1 << endl;
	}
}
	oFile5.close();


	/***************************************************************************************
	模块说明:
	`费用方程求解
	****************************************************************************************/

	ofstream File6;												//打开要输出的文件
	File6.open(filename6, ios::out | ios::trunc);

	File6 << "," << 1 << endl;
	File6 << "nProfiles" << "," << "nChannels" << "," << "MinSampTimeValue[msec]" << "," << "ChannelDistance[mm]" << "," << "VelResolution[mm/s]" << "," << "start[mm]" << "," << " GainStart" << "," << " GainEnd" << "," << "nTDX" << endl;;
	File6 << "1024" << "," << "30" << "," << "0" << "," << "0.69" << "," << "24" << "," << "0" << "," << "0" << "," << "20" << "," << "1" << endl;

	ofstream File10;												//打开要输出的文件
	File10.open(filename10, ios::out | ios::trunc);
	File10 << "Pr" << "," << "F" << endl;
	ofstream File11;												//打开要输出的文件
	File11.open(filename11, ios::out | ios::trunc);
	File11 << "Pi" << "," << "F" << endl;
	ofstream File12;												//打开要输出的文件
	File12.open(filename12, ios::out | ios::trunc);
	File12 << "v" << "," << "F" << endl;

	double  w0 =2* pi;

	double F , W, U, H;
	for (int v = -1000; v <= 500; v += 15)
	{
		F = 1000;
		for (double B = -500; B< 200; B += 15)
		{

		for (double A = -200; A< 1000; A += 15)
		{
		

			double Re1 = 0.0, Im1 = 0.0, C = 0.0;
			for (int i = 0; i < 40; i++)
			{
				
				double D = 0;
				Re1 = -A * pow(10, -3) +w0 * bbI[i] - v * DbR[i] * pow(10, -6);
				Im1 = -B * pow(10, -3) - w0 * bbR[i] -v * DbI[i] * pow(10, -6);
				D = Re1 * Re1 + Im1 * Im1;
				C += D;
			}

			if (C < F)
			{
				F = C;
				W = v;
				U = B;
				H = A;
			}
			
		}

		}
		File12 << v << "," << F << endl;
	}
	
	for (double B = -500; B< 500; B += 15)
	{
		F = 1000;
		for (int v = -1000; v <= 500; v += 15)
		{

			for (double A = -200; A< 1000; A += 15)
			{


				double Re1 = 0.0, Im1 = 0.0, C = 0.0;
				for (int i = 0; i < 40; i++)
				{

					double D = 0;
					Re1 = -A * pow(10, -3) + w0 * bbI[i] - v * DbR[i] * pow(10, -6);
					Im1 = -B * pow(10, -3) - w0 * bbR[i] - v * DbI[i] * pow(10, -6);
					D = Re1 * Re1 + Im1 * Im1;
					C += D;
				}

				if (C < F)
				{
					F = C;
					W = v;
					U = B;
					H = A;
				}

			}

		}
		File11 << B << "," << F << endl;
	}

	
	for (double A = -1000; A < 1000; A += 15)
	{
		F = 1000;
		for (double B = -500; B < 500; B += 15)
		{
			for (int v = -1000; v <= 500; v += 15)
			{


				double Re1 = 0.0, Im1 = 0.0, C = 0.0;
				for (int i = 0; i < 40; i++)
				{

					double D = 0;
					Re1 = -A * pow(10, -3) + w0 * bbI[i] - v * DbR[i] * pow(10, -6);
					Im1 = -B * pow(10, -3) - w0 * bbR[i] - v * DbI[i] * pow(10, -6);
					D = Re1 * Re1 + Im1 * Im1;
					C += D;
				}

				if (C < F)
				{
					F = C;
					W = v;
					U = B;
					H = A;
				}

			}

		}
		File10 << A << "," << F << endl;
	}
	//cout << "v=" << W << "," << "B=" << U << "," << "A=" << H << endl;
	/*
	double F1 = 10000, W1, U1, H1;
	for (double v = W - 200; v < W + 200; v += 5)
	{
		for (double B = U - 300; B < U + 200; B += 5)
		{
		
		for (double A = H - 400; A < H + 400; A += 5)
		{
		
			double Re1 = 0.0, Im1 = 0.0, C = 0.0;
			for (int i =0; i < 40; i++)
			{
				
				double D = 0;
				Re1 = -A * pow(10, -3) +w0 * bbI[i] - v * DbR[i] * pow(10, -6);
				Im1 = -B * pow(10, -3) -w0 * bbR[i] - v * DbI[i] * pow(10, -6);
				D = Re1 * Re1 + Im1 * Im1;
				C += D;
			}

			if (C < F1)
			{
				F1 = C;
				W1 = v;
				U1 = B;
				H1 = A;
			}

		}

		}
	}
	cout << "v=" << W1 << "," << "B=" << U1 << "," << "A=" << H1 << endl;
	double F2 = 10000, W2, U2, H2;
	for (double v = W1 - 30; v < W1 + 30; v += 1)
	{
		for (double B = U1 - 30; B < U1 + 30; B += 1)
		{
	
		for (double A = H1 - 30; A < H1 + 30; A += 1)
		{

			double Re1 = 0.0, Im1 = 0.0, C = 0.0;
			for (int i = 0; i < 40; i++)
			{

				double D = 0;
				Re1 = -A * pow(10, -3) + w0 * bbI[i] - v * DbR[i]* pow(10, -6);
				Im1 = -B * pow(10, -3) - w0 *bbR[i] - v * DbI[i] * pow(10, -6);
				D = Re1 * Re1 + Im1 * Im1;
				C += D;
			}

			if (C < F2)
			{
				F2 = C;
				W2 = v;
				U2 = B;
				H2 = A;
			}

		}

		}
	}
	cout << "v=" << W2 << "," << "B=" << U2 << "," << "A=" << H2 << endl;

	for (double A = H2 - 10; A < H2 + 10; A += 0.1)
	{
		File6 << "," << A;
	}
	File6 << endl;

	double F3 = 10000, W3, U3, H3;
	for (double v = W2 - 10; v < W2 + 5; v += 0.1)
	{
		File6 << v << ",";
		for (double B = U2 - 5; B < U2 + 5; B += 0.1)
		{
	
		for (double A = H2 - 10; A < H2 + 10; A += 0.1)
		{

			double Re1 = 0.0, Im1 = 0.0, C = 0.0;
			for (int i = 0; i < 40; i++)
			{

				double D = 0;
				Re1 = -A * pow(10, -3) + w0 * bbI[i] - v * DbR[i] * pow(10, -6);
				Im1 = -B * pow(10, -3) - w0 * bbR[i]  - v * DbI[i] * pow(10, -6);
				D = Re1 * Re1 + Im1 * Im1;
				C += D;
			}
			File6 << C << ",";
			if (C < F3)
			{
				F3 = C;
				W3 = v;
				U3 = B;
				H3 = A;
			}

		}

		}
		File6 << endl;
	}

	cout << "min=" << F3 << endl;
	cout << "v=" << W3 << "," << "B=" << U3 << "," << "A=" << H3 << endl;
	*/
	File6.close();
	File10.close();
	File11.close();
	File12.close();
	return 0;
}

int main()
{
	char filename[] = "C://vel data//vel_f_00_.csv";     //theo dates
	char filename3[] = "C://vel data//im.csv";				//fft im
	char filename2[] = "C://vel data//re.csv";				//fft im
	char filename4[] = "C://vel data//nihe ak.csv";			//re ak(r)
	char filename5[] = "C://vel data//nihe bk.csv";			//im bk(r)
	char filename6[] = "C://vel data//cost function.csv";    //cost function result
	char filename7[] = "C://vel data//nihe r0.0125.csv";
	char filename8[] = "C://vel data//akbk.csv";				//fft re
	char filename9[] = "C://vel data//nihe ak0.csv";				//fft re
	char filename10[] = "C://vel data//cost function_Pr.csv";    //cost function result
	char filename11[] = "C://vel data//cost function_Pi.csv";    //cost function result
	char filename12[] = "C://vel data//cost function_v.csv";    //cost function result

	dft(filename, filename2, filename3, filename4, filename5, filename6, filename7, filename8, filename9, filename10, filename11, filename12);

	return 0;
}
