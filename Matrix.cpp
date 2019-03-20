#include "Matrix.h"

// решение матрицы

Matrix::Matrix()
{
	ifstream fin("input.txt");
	fin >> maxiter >> eps;
	fin >> cond1 >> cond2 >> cond3 >> cond4 >> cond5 >> cond6 >> cond7 >> cond8;
	//fin >> nx >>  ny; // nx = n3 ny=n5
	fin >> x0 >> x1 >> x2 >> x3 >> y0 >> y1 >> y2;
	//fin >> n1 >> n2 >> nx >> n4 >> ny;
//	fin >> g1 >> g2 >> g3 >> g4 >> g5 >> g6 >> g7 >> g8;	
	fin >> kx1 >> kx2 >> kx3 >> ky1 >> ky2; //
	fin >> qx1 >> qx2 >> qx3 >> qy1 >> qy2;; //
	//n1 = n2 = n4 = 5;
	//n =(nx + 1) * (ny + 1);
	n = (kx1 + kx2 + kx3 - 2)*(ky1 + ky2 - 1);
	mn = (kx1 + kx2 + kx3 - 2);
	n1 = kx1 - 1;
	n2 = kx1 + kx2 - 2;
	n4 = ky1 - 1;
	nx = kx1 + kx2 + kx3 - 3;
	ny = ky1 + ky2 - 2;
 	w = 1;
	di.resize(n);
	l1.resize(n - 1);
	u1.resize(n - 1);
	l2.resize(n - mn - 1);
	u2.resize(n - mn - 1);
	F.resize(n);
	nev.resize(n);
	xT.resize(n);
	pog.resize(n);
	xP.resize(n);
	x.resize(n);
	for (int i = 0; i < nev.size(); i++)
	{
		nev[i] = 100;
	}	
	for (int i = 0; i < u2.size(); i++)
	{
		u2[i] = 0;
	}
	for (int i = 0; i < u1.size(); i++)
	{
		u1[i] = 0;
	}
	for (int i = 0; i < di.size(); i++)
	{
		di[i] = 0;
	}
	for (int i = 0; i < l1.size(); i++)
	{
		l1[i] = 0;
	}
	for (int i = 0; i < l2.size(); i++)
	{
		l2[i] = 0;
	}
	for (int i = 0; i < F.size(); i++)
	{
		F[i] = 0;
	}
	for (int i = 0; i < xT.size(); i++)
	{
		xT[i] = 0;
		x[i] = 0;
		xP[i] = 0;
	}
	
}

//умножение строки матрицы на вектор
real Matrix::multVV( int i, vector<real> xP)
{
	
	real res = 0;

		if (i > 0)
		{
			res += l1[i - 1] * xP[i - 1];
			if (i > mn)
			{
				res += l2[i - mn - 1] * xP[i - mn - 1];
				/*if (i > m + 2)
				{
					res += l3[i - m  - 3] * xP[i - m - 3];
				}*/
			}
		}
	
		res += di[i] * xP[i]; ///////////////////////////////////

		if (i < n - 1)
		{
			res += u1[i] * xP[i + 1];
			if (i < n - mn - 2)
			{
				res += u2[i] * xP[i + mn + 1];
				/*if (i < n - m - 3)
				{
					res += u3[i] * xP[i + m + 3];
				}*/
			}
		}
	return res;
}

// вычисление нормы в эвклидовом порстранстве
real Matrix::norm(vector<real> x)
{
	real norma = 0;
	for (int i = 0; i < n; i++)
	{
		norma += pow(x[i], 2);
	}
	norma = sqrt(norma);
	return norma;
}

void  Matrix::print()
{
	ofstream answer("answer.txt");
	answer.precision(20);
	answer << " u= pow(x,4) + pow(y,4)" << endl;
	answer << " iter=" << iter << endl;
//	answer << " obusl=" << obusl << endl;
	answer << " node " <<"\t xcoord " << "\t ycoord " << "\t u " <<"\tu* - u" << endl;
	answer.scientific;
	int ix, iy;
	
	for (int i = 0; i < n; i++)
	{
		ix = i%(nx+1);
		iy = i % (ny + 1);
		answer <<  i << "\t  " << xl[ix] << "\t  " << yl[iy] << "\t  " << x[i] << "\t  " << pog[i] << endl;
	}
	
	
}

void Matrix::run()
{  
	//if (metod == 1) jacobi(x);
	//if (metod == 2) zeidel(x);
	////if (metod == 3) zeidel(x);
	//for (int i = 0; i < n; i++)
	//{
	//	pog[i] = x[i] - xT[i];
	//}
	//pogr = norm(pog) / normxT;
	//obusl = pogr / nv;
	//print();

	Form_cross_point();
	Make_cell();
	jacobi(x);
	for (int i = 0; i < n; i++)
	{
		pog[i] = x[i] - xT[i];
	}
	normxT = norm(xT);
	pogr = norm(pog) / normxT;
	obusl = pogr / nv;
	print();
}

// метод якоби
void Matrix::jacobi( vector<real> &x) /////////////////////////////////////// выход по е
{
	//m = nx;
	normF = norm(F);
	int exit = 0;
	for (int k = 1; k < maxiter && exit != 1; k++)
	{
		iter = k;
				for (int i = 0; i < n; i++)
				{
					x[i] = xP[i] + (w/di[i])*(F[i] - multVV(i, xP));
				}
			xP = x;//xP.swap(x);
			cout << "iter:" << k;
			double multmw;
			for (int i = 0; i < n; i++)
				{
				multmw = multVV(i, xP);
				nev[i] = F[i] - multmw;		
				}
			 nv = norm(nev) / normF;
			cout << "      nev:" <<  nv << "\r";
			if ((nv) < eps)
			{
				exit = 1;
				cout << "EXIT" << endl;
			}
	}
	
}
// метод гаусса 
void Matrix::zeidel(vector<real> &x) ///////
{
	int exit = 0;
	for (int k = 1; k < maxiter&& exit != 1; k++)
	{
		iter = k;
			for (int i = 0; i < n; i++)
			{
				x[i] = x[i] +( w / di[i])*(F[i] - multVV(i, x));
			}
			cout << "k:" << k;
			for (int i = 0; i < n; i++)
			{
				nev[i] = F[i] - multVV(i, x);
			}
			 nv = norm(nev) / normF;

			cout << "    nev:" << nv << "\r";
			if ((nv) < eps)
			{
				exit = 1;
				cout << "\n    EXIT" << endl;
			}
	}
}
