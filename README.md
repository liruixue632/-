
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

　const int y = 5;																//y-1次方程式

int dft(char filename[], char filename2[], char filename3[], char filename4[], char filename5[], char filename6[], char filename7[], char filename8[])	
{

/**************************************************************************************
 パート说明:
  ファイルのデータをvector<string> date1に読み込む；コンマを取り除いてvector<string> date2に読み込む
  縦軸はtime(s)、横軸はパイプ半径の距離r(mm) 
*****************************************************************************************/


	fstream file;                                     
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
	rows = date1.size();                                   //データの行数
	cout << "rows_data number" << rows << endl;
	cout << date2num << endl;
	cout << "zongshujushu " << date2.size() << endl;
	cols = date2.size() / date1.size();　　　　　　　　　　　//データの列数
	cout << "cols_date number " << cols << endl;


	/***************************************************************************************
	パート说明:
      　　　　　　 String型　→　double変換
	****************************************************************************************/


	vector<double> date;
	for (int i = 0; i < date2num;)
	{
		double g = atof(date2[i].c_str());
		date.push_back(g);
		i++;
	}
	cout << "double_data number" << date.size() << endl;




　　/**************************************************************************************
　　　　 パート说明:
     　　　　　　　　各列の double型のデータにフーリエ変換をして、実部R、虚部Iと主波数を得る。
 　*****************************************************************************************/
  
	ofstream oFile2, oFile3;　　　　　　　　　　　　　　　　
	oFile2.open(filename2, ios::out | ios::trunc);　　　　//FFTフーリエ変換　実部Rファイル
	oFile3.open(filename3, ios::out | ios::trunc);　　　　//FFTフーリエ変換　虚部Iファイル　

	//int t = cols / 2;
	vector<double> R;　　　　　　　　　　　　　　　　　　　　//実部R
	vector<double> I;　　　　　　　　　　　　　　　　　　　　//虚部I
	vector<double> Rr;　　　　　　　　　　　　　　　　　　　//主波数列の実部Rr
	vector<double> Ii;　　　　　　　　　　　　　　　　　　　//主波数列の虚部Ii
	for (int j = 0; j < cols;j++)  /*第j列*/
	{

		for (int k = 0; k < rows;k++)
		{
			double Re = 0.0, Im = 0.0;
			for (int z = 0; z < rows;z++)  /*第z行*/
			{
				int n = cols * z + j;
				Re += (date[n])* cos(2 * pi*k*z / rows) / rows;     //FFTフーリエ変換
				Im += -(date[n])* sin(2 * pi*k*z / rows) / rows;　　//FFTフーリエ変換
				
			}
			oFile2 << Re << ",";　　　　　　　　　　　　　
			oFile3 << Im << ",";


			if (k == 10)　　　　　　　　　　　　　　　　　　　　//K==10 データの主波数
			{
				R.push_back(Re);
				I.push_back(Im);
				cout << "k=10" << "," << Re << "," << Im << endl;

			}
			
		}
		oFile2 << endl;
		oFile3 << endl;
		
	}
	oFile2.close();
	oFile3.close();

	ofstream oFile8;
	oFile8.open(filename8, ios::out | ios::trunc);

	vector<double> ak;　　　　　　　　　　　　　　　　　　　　　//主波数の実部Ｒ→ak
	vector<double> bk;                                       //主波数の虚部Ｉ→bk
	oFile8 << " ak" << "," << " bk " << endl;
	for (int i = 0; i < cols; i++)
	{
		ak.push_back(R[i]);
		bk.push_back(I[i]);
		oFile8 << R[i] << "," << I[i] << endl;
	}


/***************************************************************************************
　パート说明:
方程式的係数xR[n]の求める　　主波数実部akをパイプ半径の関数  ak=xR[0]+xR[1]*r+xR[2]*r*r+xR[3]*r*r*r+xR[4]*r*r*r*r
****************************************************************************************/


	int t;
	vector<double> DbR;												//実部
	vector<double> bbR;												//実部

	ofstream oFile4;												
	oFile4.open(filename4, ios::out | ios::trunc);

	for (t = 0; t < 3; t++)
	{
		double xR[y];												
		double s[y];											
		double m[y][y];												
		double rr =  t * 0.2 / 29;							
		for (int i = 0; i < y; i++)									 
		{
			for (int j = 0; j < y; j++)
			{

				m[i][j] = pow(rr, j);
			}
			rr += 2 * 0.02 / 29;
		}
		for (int i = 0; i < y; i++)								// 实部行列式右侧
		{
			s[i] = ak[2 * i + 10 * t];
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

		rr = t * 0.2 / 29;
		double q = rr + 0.2 / 29;
		double cr = 0.02 / 29;
		for (double r = rr; r < q;r+=cr)			//实部
		{
			double DB = 0, b1 = 0;
			for (i = 1; i < y; i++)
			{
				DB += pow(i, 2)* xR[i] * pow(r, i - 2);
			}
			//cout << "DbR" << r << "=" << DB << endl;
			DbR.push_back(DB);

			for (i = 0; i < y; i++)
			{
				b1 += xR[i] * pow(r, i);
			}

			bbR.push_back(b1);
			//cout << "bbR" << r << "=" << b1 << endl;
			oFile4 << r << "," << DB << "," << b1 << endl;
		}
	}
	oFile4.close();

/***************************************************************************************
 パート说明:
 方程式的係数xI[n]の求める　　主波数虚部bkをパイプ半径の関数  bk=xI[0]+xI[1]*r+xI[2]*r*r+xI[3]*r*r*r+xI[4]*r*r*r*r
 ****************************************************************************************/

	vector<double> DbI;												
	vector<double> bbI;												

	ofstream oFile5;												
	oFile5.open(filename5, ios::out | ios::trunc);

	for (t = 0; t < 3; t++)
	{
		double xI[y];												
		double s[y];												
		double m[y][y];												
		double rr = t * 0.2 / 29;							
		for (int i = 0; i < y; i++)								 
		{
			for (int j = 0; j < y; j++)
			{

				m[i][j] = pow(rr, j);
			}
			rr += 2 * 0.02 / 29;
		}
		for (int i = 0; i < y; i++)
		{
			s[i] = bk[2 * i + 10 * t];
		}

		double c[y];											 
		int i;
		for (int k = 0; k < y - 1; k++)							
		{
			for (i = k + 1; i < y; i++)							
				c[i] = m[i][k] / m[k][k];

			for (i = k + 1; i < y; i++)							
			{
				for (int j = 0; j < y; j++)
				{
					m[i][j] = m[i][j] - c[i] * m[k][j];
				}
				s[i] = s[i] - c[i] * s[k];
			}
		}

		xI[y - 1] = s[y - 1] / m[y - 1][y - 1];							
		for (int i = y - 2; i >= 0; i--)								
		{
			double sum = 0;
			for (int j = i + 1; j < y; j++)
			{
				sum += m[i][j] * xI[j];
			}
			xI[i] = (s[i] - sum) / m[i][i];							
		}

		cout << " the solution of the equations is:" << endl;
		cout << endl;
		for (i = 0; i < y; i++)
			cout << "xubuxishu a" << i << "=" << xI[i] << endl;

		rr = t * 0.2 / 29;
		double q = rr + 0.2 / 29;
		double cr = 0.02 / 29;
		for (double r = rr; r < q; r += cr)									
		{
			double DB = 0, b1 = 0;
			for (i = 1; i < y; i++)
			{
				DB += pow(i, 2)* xI[i] * pow(r, i - 2);
			}
		
			DbI.push_back(DB);

			for (i = 0; i < y; i++)
			{
				b1 += xI[i] * pow(r, i);
			}

			bbI.push_back(b1);
		
			oFile5 << r << "," << DB << "," << b1 << endl;
		}
	}
	oFile5.close();

　/***************************************************************************************
	パート说明:
			cost function の求める
****************************************************************************************/

	ofstream File6;												
	File6.open(filename6, ios::out | ios::trunc);

	File6 << "," << 1 << endl;
	File6 << "nProfiles" << "," << "nChannels" << "," << "MinSampTimeValue[msec]" << "," << "ChannelDistance[mm]" << "," << "VelResolution[mm/s]" << "," << "start[mm]" << "," << " GainStart" << "," << " GainEnd" << "," << "nTDX" << endl;;
	File6 << "1024" << "," << "30" << "," << "0" << "," << "0.69" << "," << "24" << "," << "0" << "," << "0" << "," << "20" << "," << "1" << endl;　　//ファイルの横軸
	for (int A = 0; A < 800; A++)　　　　　　
	{
		File6 << "," << A;　　　　　　　　　　　//ファイルの縦軸
	}
	File6 << endl;
	double  w0 = pi;

	double F = 10000, W, U, H;
	for (int v = 0; v<= 1000; v += 15)             //v:粘性
	{
		for (double B= -500; B< 2000; B += 15)　　　　//B:圧力の虚部
		{
			for (double A = -1500;A< 1000; A+= 15)      //A:圧力の実部
			{
				double Re1 = 0.0, Im1 = 0.0, C = 0.0;
				for (int i = 1; i < 30; i++)
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

	}
	cout << "v=" << W << "," << "B=" << U << "," << "A=" << H<< endl;
	double F1 = 10000, W1, U1, H1;
	for (double v = W - 200; v < W + 200; v+= 5)
	{
		for (double B= U - 300; B < U+ 200;B += 5)
		{
			for (double A =H- 400; A < H+ 400;A += 5)
			{
				double Re1 = 0.0, Im1 = 0.0, C = 0.0;
				for (int i = 1; i < 30; i++)
				{
					
					double D = 0;
					Re1 = -A * pow(10, -3) + w0 * bbI[i] - v * DbR[i] * pow(10, -6);
					Im1 = -B * pow(10, -3) - w0 * bbR[i] - v * DbI[i] * pow(10, -6);
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
		for (double B =U1 - 30; B < U1 + 30; B+= 1)
		{
			for (double A = H1 - 30; A < H1 + 30; A += 1)
			{		
				double Re1 = 0.0, Im1 = 0.0, C = 0.0;
				for (int i = 1; i < 30; i++)
				{
					double D = 0;
					Re1 = -A * pow(10, -3) + w0 * bbI[i] - v * DbR[i] * pow(10, -6);
					Im1 = -B * pow(10, -3) - w0 * bbR[i] - v * DbI[i] * pow(10, -6);
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
	double F3 = 10000, W3, U3, H3;
	for (double v = W2 - 10; v < W2 + 5; v += 0.1)
	{
		for (double B= U2 - 5;B < U2 + 5; B += 0.1)
		{
			for (double A = H2 - 10;A < H2 + 10; A += 0.1)
			{
				
				double Re1 = 0.0, Im1 = 0.0, C = 0.0;
				for (int i = 1; i < 30; i++)
				
					double D = 0;
					Re1 = -A * pow(10, -3) + w0 * bbI[i] - v * DbR[i] * pow(10, -6);
					Im1 = -B * pow(10, -3) - w0 * bbR[i] - v * DbI[i] * pow(10, -6);
					D = Re1 * Re1 + Im1 * Im1;
					C += D;
				}

				if (C < F3)
				{
					F3 = C;
					W3 = v;
					U3 = B;
					H3 = A;
				}

			}

		}
	}

	File6 << endl;
	cout << "min=" << F3 << endl;
	cout << "v=" << W3 << "," << "B=" << U3 << "," << "A=" << H3 << endl;
	File6.close();

	return 0;
}

int main()
{
	char filename[] = "F://theo data 100_157_0_.csv";     //理論データ、粘性v=100、圧力Pの実部A＝314、虚部＝０(輸入)
	char filename3[] = "F://im.csv";		      //FFT im(r)　(輸出)
	char filename2[] = "F://re.csv";		      //FFT Re(r)　　(輸出)
	char filename4[] = "F://nihe ak.csv";		      //主波数Re ak(r)　　(輸出)　
	char filename5[] = "F://nihe bk.csv";	              //主波数im bk(r)　　(輸出)
	char filename6[] = "F://cost function.csv";           //cost function result　　(輸出)
	char filename8[] = "F://akbk.csv";		    
	
	dft(filename, filename2, filename3, filename4, filename5, filename6, filename7, filename8);

	return 0;
	
}
