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
	int cond1, cond2, cond3, cond4, cond5, cond6, cond7, cond8;//  тип краевого
	int n1, n2,  n4, nx, ny, mn; //  номера узлов границ
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
	double  qx1, qx2, qx3, qy1, qy2;// // шаг по х и по у, коэф. разрядки, значение краевых условий
	double *xl, *yl, **xy, *hx, *hy; // координаты узлов
	int n, m, maxiter; // размер матрицы, промежуток, максимум итераций, размер блока
	vector<real> di, l1, l2, l3, u1, u2, u3, F, xP, x, nev, xT, pog; // матрица, вектор правой части, приближение, х текущее, , х точное
	real eps, normF, normxT; // погрешность, норма вектора F
	real w, obusl, nv, pogr; // параметр релаксации, число обусловленности
	int metod, iter; // номер метода, колво выполненных итераций
	Matrix();
	void run();
	void print();
	real multVV(int i, vector<real> xP);
	void jacobi( vector<real> &x);
	void zeidel( vector<real> &x);
	real norm(vector<real> x);
	//первая производная 
	double First_Grad(double x, double y, int index);
	//вычисление значений на границе 
	void Calc_matrix_board(int x, int y, int cond);
	//вычисление значений в вырезанной области
	void Calc_matrix_empty(int x, int y);
	//вычисление значений внутри области
	void Calc_matrix_in(int x, int y);
	//функция для лямбды
	double Lambda(double x, double y);
	//функция для сигмы
	double Sigma(double x, double y);
	//-(сумма вторых частных производных)
	double Function(double x, double y);
	//функция u
	double Your_function(double x, double y);
	//получение глобального номера узла
	int Get_a_global_number(int x, int y);
	//заполнениие массива узлов
	void Form_cross_point();
	//заполнение  пяти-диагонально матрицы
	void Make_cell();
};


