#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;
double Steppolinom(double,int,double**);
void Matrix_inversion(double**,int);
void Matrix_multiply(double**,double**,double**,int,int,int);

void Polinom(double* polinom,double* y_values,int e,int s){
ofstream fpolinom;
double** secondmat = new double* [s + 1];
for (int i = 0; i < s + 1; i++)
  secondmat[i] = new double[1];	///матрица решений
for (int i = 0; i < s + 1; i++)
  secondmat[i][0] = 0;

for (int i = 0; i < s + 1; i++){
  for (int j = 0; j < e; j++)
    secondmat[i][0] += polinom[j] * pow(y_values[j], i);
}


/// Создадим матрицу сумм произведений х в разных степенях|| САМАЯ БОЛЬШАЯ МАТРИЦА В УРАВНЕНИИ
double** firstmat = new double* [s + 1];
for (int i = 0; i < s + 1; i++)
  firstmat[i] = new double[s + 1];

for (int i = 0; i < s + 1; i++){
  for (int j = 0; j < s + 1; j++){
    if (i == 0 && j == 0)
      firstmat[i][j] = e;
    else {
      for (int k = 0; k < e; k++)
        firstmat[i][j] += pow(y_values[k], i + j);
    }
    }
}

///Создадим матрицу коэффициентов
double** koef = new double* [s + 1];
for (int i = 0; i < s + 1; i++)
  koef[i] = new double[1];
for (int i = 0; i < s + 1; i++)
  koef[i][0] = 0.0;

Matrix_inversion(firstmat, s + 1);

Matrix_multiply(firstmat, secondmat, koef, s + 1, s + 1, 1);

fpolinom.open("Polinom.dat");

int k = y_values[e - 1];

for (double i = 0; i < k; i += 0.1){
   fpolinom <<  Steppolinom(i, s, koef) << " " <<  i << endl;
  }
  fpolinom.close();
  system("python3 pol.py");
}

double Steppolinom(double x, int s, double** koef)
{
	double y = 0.0;
	for (int i = 0; i <= s; i++)  y += koef[i][0] * pow(x, i);

	return y;
}

//ОБРАТНАЯ МАТРИЦА
void Matrix_inversion(double** Inverse_Matrix, int N)
{
  double temp;
  double** UnitMatrix = new double* [N];

  for (int i = 0; i < N; i++)
		UnitMatrix[i] = new double[N];

  for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
        //Создание единичной матрицы (главная диогональ 1)
            UnitMatrix[i][j] = 0.0;
            if (i == j)
                UnitMatrix[i][j] = 1.0;
        }
    }
    //-------------

    for (int k = 0; k < N; k++) {
        temp = Inverse_Matrix[k][k];
        for (int j = 0; j < N; j++){
            Inverse_Matrix[k][j] /= temp;
            UnitMatrix[k][j] /= temp;
        }
        for (int i = k + 1; i < N; i++){
            temp =  Inverse_Matrix[i][k];
            for (int j = 0; j < N; j++){
                Inverse_Matrix[i][j] -= Inverse_Matrix[k][j] * temp;
                UnitMatrix[i][j] -=  UnitMatrix[k][j] * temp;
            }
        }
    }

    for (int k = N - 1; k > 0; k--){
        for (int i = k - 1; i >= 0; i--){
            temp = Inverse_Matrix[i][k];
            for (int j = 0; j < N; j++){
                Inverse_Matrix[i][j] -= Inverse_Matrix[k][j] * temp;
                UnitMatrix[i][j] -= UnitMatrix[k][j] * temp;
            }
        }
    }

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			  Inverse_Matrix[i][j] =  UnitMatrix[i][j];

	for (int i = 0; i < N; i++)
		delete[]  UnitMatrix[i];

	delete[]  UnitMatrix;

}
///Перемножение матриц с d - строки и столбцы первой матрицы, f - столбец второй

void Matrix_multiply(double** A, double** B, double** C, int c, int d, int f)
{
	for (int i = 0; i < c; i++)
	{
		for (int j = 0; j < f; j++)
		{
			for (int r = 0; r < d; r++)
				C[i][j] += A[i][r] * B[r][j];
		}
	}
}



double a_cof(double u, int n) {
 double temp = u;
 for (int i = 1; i < n; i++) temp *= (u - i);
 return temp;
}

double a_cof2(double u, int n) {
 double temp = u;
 for (int i = 1; i < n; i++) temp *= (u + i);
 return temp;
}

double fact(int n) {
 double f = 1.;
 for (int i = 1; i <= n; i++) f *= i;
 return f;
}


double Newton_two(double* x, double* y, short p, double _x) {
 //(назад)
 int n = p;
 double yy[n][n];
 for (int i = 0; i < n; i++) yy[i][0]=x[i];
 for (int i = 1; i < n; i++) {
  	for (int j = n - 1; j >= i; j--)
			yy[j][i] = yy[j][i - 1] - yy[j - 1][i - 1];
 }
  double sum = yy[n-1][0];
  double u = (_x - y[n-1]) / (y[1] - y[0]);
  for (int i = 1; i < n; i++) {
   sum += (a_cof2(u, i) * yy[n-1][i]) / fact(i);
  }
 return sum;
}


double Newton_one(double* x, double* y, short n, double _x) {

 //(вперёд)
 double yy[n][n];
 for (int i = 0; i < n; i++){ yy[i][0] = x[i];  cout << yy[i][0]  << " "; }
 cout << endl;
 for (int i = 1; i < n; i++){
  for (int j = 0; j < n - i; j++){
		yy[j][i] = yy[j + 1][i - 1] - yy[j][i - 1];
    cout << yy[j][i]  << " ";
    }
 }
  double sum = yy[0][0];
  double u = (_x - y[0]) / (y[1] - y[0]);
  for (int i = 1; i < n; i++) {
   sum += (a_cof(u, i) * yy[0][i]) / fact(i);
  }
  return sum;
}




double Lagrange(double* x, double* y, short n, double _x)
{
	double result = 0.0;
	for (short i = 0; i < n; i++)
	{
		double P = 1.0;
		for (short j = 0; j < n; j++)
			if (j != i)
				P *= (_x - y[j])/ (y[i] - y[j]);
    if(x[i] != 999.9)
		  result += P * x[i];
	}
  if(result < 100){
	   return result;
   }
}




void Select(double** data, int row,int col){
    ofstream flagrange,fnewton_one,fnewton_two,fpolinom,origlag,orignew1,orignew2,origpol;
    int r,c,e;
    system("clear");
    cout << endl;
    cout << "Select data(row number,column number of elements. 'row col elements'): "; cin >> r >> c >> e;
    double* x_values = new double [e];
    double* y_values = new double [e];
    double* year = new double [e];
    double* yy_values = new double [row];
    double* polinom = new double[row];
     for (int i = 0; i < e; i++){
            x_values[i] =  data[i+r][c+1];
            y_values[i] = i;
            year[i] = data[i][0];
            cout << x_values[i] << " | " << year[i] << endl;
          }
    for(int i = 0; i < row; i++){
      polinom[i] = data[i+r][c+1];
      yy_values[i] = i;
    }
        cout << string(1, ' ');
        cout << endl;
        int select = 0;
        cin.clear();
        cout << "Select Method: " << endl
             << "1. Lagrange." << endl
             << "2. Newton_one. " << endl
             << "3. Newton_two." << endl
             << "4. Polinom." << endl;
        cin >> select;
        double di = 0.1;
        int s = 0;
          switch (select){
            case 1: flagrange.open("Lagrange.dat"); for(double i = 0, j = 0; i < e; i += di, j++)  flagrange << Lagrange(x_values,y_values,e,i) << " " << i << endl;
            flagrange.close(); system("python3 lag.py"); origlag.open("orignlag.dat"); for(int i = 0; i < e; i++) origlag << x_values[i] << " " << y_values[i] << endl; origlag.close(); break;
            case 2: fnewton_one.open("Newton_one.dat"); for(double i = 0; i < e; i += di) fnewton_one << Newton_one(x_values,y_values,e,i) << " " << i << endl;
            fnewton_one.close(); system("python3 new1.py");orignew1.open("orignew1.dat"); for(int i = 0; i < e; i++) orignew1 << x_values[i] << " " << y_values[i] << endl; orignew1.close(); break;
            case 3: fnewton_two.open("Newton_two.dat"); for(double i = 0; i < e; i += di) fnewton_two << Newton_two(x_values,y_values,e,i) << " " << i << endl;
            fnewton_two.close(); system("python3 new2.py"); orignew2.open("orignew2.dat"); for(int i = 0; i < e; i++)  orignew2 << x_values[i] << " " << y_values[i] << endl; orignew2.close(); break;
            case 4: cout << "Degree: "; cin >> s; Polinom(polinom,yy_values,row,s); origpol.open("origpol.dat"); for(int i = 0; i < row; i++) origpol << polinom[i] << " " << yy_values[i] << endl; origpol.close(); break;
            case 0: break;
          }
    }


void delete_array2d(double** A,int size){
  for (int i=0; i < size; i++)
    delete A [i];
  delete [] A;
}



int main()
{
    fstream fdata;
    double **matr,**temp;
    char *str = new char [1024];
    int size = 0;
    double t = 0;
    int row = 0,col = 0;
    string filename = "data.txt";

    fdata.open(filename,ios::in);
    if(fdata){
      while(!fdata.eof()){
        fdata >> t;
        size++;
      }
    }
    fdata.close();

    fdata.open(filename,ios::in);
    if(fdata){
      while(!fdata.eof()){
        fdata.getline(str, 1024, '\n');
          row++;
        }
      row--;
      col = size/(row - 1);
    }
    fdata.close();

    matr = new double* [size];

    for(int i = 0; i < size; i++){
      matr[i] = new double[size];
    }


    fdata.open(filename,ios::in);
    if(fdata){
      for (int i = 0; i < row;i++){
        for(int j = 0; j < col; j++){
          if(matr[i][j] != 999.9)
            fdata >> matr[i][j];
        }
      }
    }
    fdata.close();

    Select(matr,row,col);

    delete_array2d(matr,size);
    return 0;
}
