#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <conio.h>
#include <math.h>
#include <vector>
#include <fstream>
#include <locale>
#include <iostream>
using namespace std;

typedef double real;
typedef double dubl;

class Matrix
{
public:
	int cond1, cond2, cond3, cond4, cond5, cond6, cond7, cond8;//  ��� ��������
	int n1, n2,  n4, nx, ny, mn; //  ������ ����� ������
	// _________________ ny
	//|                 |
	//|     _______     | n4
	//|    |       |    |
	//|    |       |    | 
	//|    |       |    |
	//|    |       |    |
	//|____|       |____|
	//     n1     n2    nx
	double x0, x1, x2, x3, y0, y1, y2;
	int  kx1, kx2, kx3, ky1, ky2; ;
	double  qx1, qx2, qx3, qy1, qy2;// // ��� �� � � �� �, ����. ��������, �������� ������� �������
	double *xl, *yl, **xy, *hx, *hy; // ���������� �����
	int n, m, maxiter; // ������ �������, ����������, �������� ��������, ������ �����
	vector<real> di, l1, l2, l3, u1, u2, u3, F, xP, x, nev, xT, pog; // �������, ������ ������ �����, �����������, � �������, , � ������
	real eps, normF, normxT; // �����������, ����� ������� F
	real w, obusl, nv, pogr; // �������� ����������, ����� ���������������
	int metod, iter; // ����� ������, ����� ����������� ��������
	Matrix();
	void run();
	void print();
	real multVV(int i, vector<real> xP);
	void jacobi( vector<real> &x);
	void zeidel( vector<real> &x);
	real norm(vector<real> x);
	//������ ����������� 
	double First_Grad(double x, double y, int index);
	//���������� �������� �� ������� 
	void Calc_matrix_board(int x, int y, int cond);
	//���������� �������� � ���������� �������
	void Calc_matrix_empty(int x, int y);
	//���������� �������� ������ �������
	void Calc_matrix_in(int x, int y);
	//������� ��� ������
	double Lambda(double x, double y);
	//������� ��� �����
	double Sigma(double x, double y);
	//-(����� ������ ������� �����������)
	double Function(double x, double y);
	//������� u
	double Your_function(double x, double y);
	//��������� ����������� ������ ����
	int Get_a_global_number(int x, int y);
	//����������� ������� �����
	void Form_cross_point();
	//����������  ����-����������� �������
	void Make_cell();
};


