#include"Matrix.h"
// построение матрицы

double Matrix::Lambda(double x, double y)
{
	return 1.0;
}

double Matrix::Sigma(double x, double y)
{
	return 0;
}

double Matrix::Function(double x, double y) //-(dx^2 + dy^2)
{
	return -1*(12*x*x + 12*y*y);
}

double Matrix::Your_function(double x, double y)  // искомая функция u
{
	return pow(x,4) + pow(y,4);
}

double Matrix::First_Grad(double x, double y, int index) {//первая производная
	if (index == 1) return 4 * pow(x, 3); /////x
	if (index == 2) return 4 * pow(y, 3); ////y
}

void Matrix::Form_cross_point()
{
	
	xl = new double[nx];
	yl = new double[ny];
	hx = new double[nx];
	hy = new double[ny];
	xl[0] = x0;
	yl[0] = y0;
	//double h = (x1 - x0)*(1 - qx1) / (1 - pow(qx1, kx1));
	double h = (x1 - x0) / (kx1 - 1);
	
	for (int i = 1; i < kx1; i++)
	{
		xl[i] = xl[i - 1] + h;
		h *= qx1;
		
	}
	//h = (x2 - x1)*(1 - qx2) / (1 - pow(qx2, kx2));
	h = (x2 - x1) / (kx2 - 1);
	
	for (int i = kx1; i < n2; i++)
	{
		xl[i] = xl[i - 1] + h;
		h *= qx2;
		
	}
	//h = (x3 - x2)*(1 - qx3) / (1 - pow(qx3, kx3));
	h = (x3 - x2) / (kx3 - 1);

	for (int i = n2; i <= nx; i++)
	{
		xl[i] = xl[i - 1] + h;
		h *= qx3;
		
	}
	//h = (y1 - y0)*(1 - qy1) / (1 - pow(qy1, ky1));
	h = (y1 - y0) / (ky1 - 1);


	for (int i = 1; i < ky1; i++)
	{
		yl[i] = yl[i - 1] + h;
		h *= qy1;
		
	} 
	//h = (y2 - y1)*(1 - qy2) / (1 - pow(qy2, ky2));
	h = (y2 - y1) / (ky2 - 1);

	
	for (int i = ky1; i <= ny; i++)
	{
		yl[i] = yl[i - 1] + h;
		h *= qy1;
		
	}
	
}

void Matrix::Calc_matrix_board(int x, int y, int cond) { 
	int m = Get_a_global_number(x, y);
	if (cond == 2) {
		if (x == n1 || x == nx) {

			di[m] = (Lambda(xl[x], yl[y]) / (xl[x] - xl[x - 1]));//hx
			l1[m - 1] = (-Lambda(xl[x], yl[y]) / (xl[x] - xl[x - 1]));//hx
			F[m] = First_Grad(xl[x], yl[y], 1);
		}
		else {
			if (x == 0 || x == n2) {
				di[m] = (Lambda(xl[x], yl[y]) / (xl[x + 1] - xl[x])); //hx
				u1[m] = (-Lambda(xl[x], yl[y]) / (xl[x + 1] - xl[x]));//hx
				F[m] = -First_Grad(xl[x], yl[y], 1);//////----------------
			}
			else {
				if (y == ny) {
					di[m] = (Lambda(xl[x], yl[y]) / (yl[x] - yl[x - 1]));//hy
					l2[m - (nx + 1)] = (-Lambda(xl[x], yl[y]) / (yl[x] - yl[x - 1]));//hy
					F[m] = First_Grad(xl[x], yl[y], 2);
				}
				else {
					if (y == 0 || y == n4) {
						di[m] = (Lambda(xl[x], yl[y]) / (yl[x + 1] - yl[x]));//hy
						u2[m] = (-Lambda(xl[x], yl[y]) / (yl[x + 1] - yl[x]));//hy
						F[m] = -First_Grad(xl[x], yl[y], 2);
					}
					else {
						cout << "not complete";
						cin.get();
					}
				}
			}
		}
		 ///=======const=======
		xT[m] = Your_function(xl[x], yl[y]);
	}
	if(cond == 1)// 1 краевые
	{ 
	di[m] = 1;
	F[m] = Your_function(xl[x], yl[y]); //Function(xl[x], yl[y]);======const=====
	xT[m] = Your_function(xl[x], yl[y]);//Function(xl[x], yl[y])
	}	

}

void Matrix::Calc_matrix_empty(int x, int y) {
	int m = Get_a_global_number(x, y);
	di[m] = 1;
	F[m] = 0;
	xT[m] = 0;
}

void Matrix::Calc_matrix_in(int x, int y) //////////////////////////////
{
	int m = Get_a_global_number(x, y);
	di[m] = 2 * Lambda(xl[x], yl[y]) * (1 / pow((xl[x + 1] - xl[x - 1])/2,2) + 1 / (pow((xl[x + 1] - xl[x - 1]) / 2, 2))) + Sigma(x,y);
	l1[m - 1] = (-1) / (pow((xl[x + 1] - xl[x - 1]) / 2, 2));
	l2[m - (nx + 1)] = (-1) / (pow((yl[x + 1] - yl[x - 1]) / 2, 2));
	u1[m] = (-1) / (pow((xl[x + 1] - xl[x - 1]) / 2, 2));
	u2[m] = (-1) / (pow((yl[x + 1] - yl[x - 1]) / 2, 2));
	F[m] = Function(xl[x], yl[y]); ///
	xT[m] = Your_function(xl[x], yl[y]);///////////

}

void Matrix::Make_cell() {

	int x, y;
	//первая граница
	for (x = 0, y = 0; x < n1; x++) Calc_matrix_board(x, y, cond1);
	//вторая граница
	for (x = n1, y = 0; y < n4; y++) Calc_matrix_board(x, y, cond2);
	//третья граница
	for (x = n1, y = n4; x < n2 ; x++)
		Calc_matrix_board(x, y, cond3);
	//четвертая граница
	for (x = n2, y = 0; y < n4; y++) Calc_matrix_board(x, y, cond4);
	//пятая граница
	for (x = n2, y = 0; x < nx; x++) Calc_matrix_board(x, y, cond5);
	//шестая граница
	for (x = nx, y = 1; y < ny; y++)
		Calc_matrix_board(x, y, cond6);
	//седьмая граница
	for (x = 1, y = ny; x < nx; x++) 
		Calc_matrix_board(x, y, cond7);
	//восьмая граница
	for (x = 0, y = 1; y < ny; y++)
		Calc_matrix_board(x, y, cond8);
	////для пустоты
	for (x = n1 + 1; x < n2; x++)
		for (y = 0; y < n4 ; y++)
			Calc_matrix_empty(x, y);

	//для внутренностей
	for (x = 1 ;  x < n1  ; x++)
		for (y = 1; y < ny ; y++)
			Calc_matrix_in(x,y);
	for (x = n2 + 1; x < nx ; x++)
		for (y = 1; y < ny ; y++)
			Calc_matrix_in(x, y);
	for (y = n4 + 1; y < ny; y++)
		for (x = n1 ; x <= n2 ; x++)		
			Calc_matrix_in(x, y);

	if (cond1 == 2 && cond8 == 2)
	{
		x = 0;
		y = 0;
		m = Get_a_global_number(x, y);
		di[m] = (Lambda(xl[x], yl[y]) / (xl[x + 1] - xl[x])); //=-la
		u1[m] = (-Lambda(xl[x], yl[y]) / (xl[x + 1] - xl[x]));
		F[m] = -First_Grad(xl[x], yl[y], 1);
		xT[m] = Your_function(xl[x], yl[y]);
	}
	else Calc_matrix_board(0, 0, 1);

	if (cond7 == 2 && cond8 == 2) {
		x = 0;
		y = ny;
		m = Get_a_global_number(x, y);
		di[m] = (Lambda(xl[x], yl[y]) / (xl[x + 1] - xl[x]));//= (-Lambda(xl[x], yl[y]) / hx);
		u1[m] = (-Lambda(xl[x], yl[y]) / (xl[x + 1] - xl[x]));
		F[m] = -First_Grad(xl[x], yl[y], 1);
		xT[m] = Your_function(xl[x], yl[y]);
	}
	else Calc_matrix_board(0, ny, 1);

	if (cond1 == 2 && cond2 == 2){
	x = n1;
	y = 0;
	m = Get_a_global_number(x, y);
	di[m] = (Lambda(xl[x], yl[y]) / (xl[x] - xl[x - 1]));
	l1[m - 1] = (-Lambda(xl[x], yl[y]) / ((xl[x] - xl[x - 1])));
	F[m] = First_Grad(xl[x], yl[y], 1);
	xT[m] = Your_function(xl[x], yl[y]);
	}
	else Calc_matrix_board(n1, 0, 1);

	if (cond2 == 2 && cond3 == 2) {
		x = n1;
		y = n4;
		m = Get_a_global_number(x, y);
		di[m] = (Lambda(xl[x], yl[y]) / (xl[x] - xl[x - 1]));
		l1[m - 1] = (-Lambda(xl[x], yl[y]) / (xl[x] - xl[x - 1]));
		F[m] = First_Grad(xl[x], yl[y], 1);
		xT[m] = Your_function(xl[x], yl[y]);
	}
	else Calc_matrix_board(n1, n4, 1);

	if (cond4 == 2 && cond5 == 2) {
	x = n2;
	y = 0;
	m = Get_a_global_number(x, y);
		di[m] = (Lambda(xl[x], yl[y]) / (xl[x + 1] - xl[x])); //=-la
		u1[m] = (-Lambda(xl[x], yl[y]) / (xl[x + 1] - xl[x]));
		F[m] = -First_Grad(xl[x], yl[y], 1);
	xT[m] = Your_function(xl[x], yl[y]);
	}
		else Calc_matrix_board(n2, 0, 1);

	if (cond3 == 2 && cond4 == 2) {
	x = n2;
	y = n4;
	m = Get_a_global_number(x, y);
	di[m] = (Lambda(xl[x], yl[y]) / (xl[x + 1] - xl[x])); //=-la
	u1[m] = (-Lambda(xl[x], yl[y]) / (xl[x + 1] - xl[x]));
	F[m] = -First_Grad(xl[x], yl[y], 1);
	xT[m] = Your_function(xl[x], yl[y]);
	}
		else Calc_matrix_board(n2, n4, 1);

	if (cond5 == 2 && cond6 == 2) {
		x = nx;
		y = 0;
		m = Get_a_global_number(x, y);
		di[m] = (Lambda(xl[x], yl[y]) / (xl[x] - xl[x - 1]));
		l1[m - 1] = (-Lambda(xl[x], yl[y]) / (xl[x] - xl[x - 1]));
		F[m] = First_Grad(xl[x], yl[y], 1);
		xT[m] = Your_function(xl[x], yl[y]);
	}
		else Calc_matrix_board(nx, 0, 1);

	if (cond6 == 2 && cond7 == 2) {
		x = nx;
		y = ny;
		m = Get_a_global_number(x, y);
		di[m] = (Lambda(xl[x], yl[y]) / (xl[x] - xl[x - 1]));
		l1[m - 1] = (-Lambda(xl[x], yl[y]) / (xl[x] - xl[x - 1]));
		F[m] = First_Grad(xl[x], yl[y], 1);
		xT[m] = Your_function(xl[x], yl[y]);
	}
		else Calc_matrix_board(nx, ny, 1);
	
}

int Matrix::Get_a_global_number(int x, int y)
{
	int number = y * (nx + 1) + x;
	return number;
}
