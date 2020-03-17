#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

int N = 20; // количество промежуточных точек
int K = 40; // дополнительное разбиение отрезков для построения графика

double* Jordan(double** a, double** b, int N)
{

    double* x = new double[N];

    double R;
    for (int k = 0; k < N; k++)
    {
        for (int i = 0; i < N; i++)
        {
            if (i != k)
            {
                R = a[i][k] / a[k][k];
                for (int j = k; j < N; j++)
                {
                    a[i][j] -= a[k][j] * R;
                }
                b[0][i] -= b[0][k] * R;
            }
        }
    }
    for (int i = 0; i < N; i++)
    {
        x[i] = b[0][i] / a[i][i];
    }

    return x;
}

double* CubicInterpolation(double* src_x, double* src_y, double* new_x) { // src_x - начальные х, src_y - начальные у, new_x - новые х для нахождения у

    double* h = new double[N - 1];
    double** gamma = new double* [1];
    *gamma = new double[N - 2];
    double** matrix = new double* [N - 2];
    double* c = new double[N - 2]; // коэффициенты с
    double* d = new double[N - 1]; // коэффициенты d
    double* b = new double[N - 1]; // коэффициенты b
    double* new_y = new double[(N - 1) * K + N]; // массив новых y
    double* ci = new double[N - 1]; // массив для хранения коэффициентов с
    double* a = new double[N - 1]; // коэффициенты а

    for (int i = 0; i < N - 2; i++)
    {
        matrix[i] = new double[N - 2];
    }

    for (int i = 0; i < N - 1; i++) // вычисление h
    {
        h[i] = src_x[i + 1] - src_x[i];
    }

    for (int i = 0; i < N - 2; i++)  // вычисление gamma
    {
        gamma[0][i] = 3 * ((src_y[i + 2] - src_y[i + 1]) / (h[i + 1]) - (src_y[i + 1] - src_y[i]) / (h[i]));
    }

    for (int i = 0; i < N - 2; i++)
    {
        for (int j = 0; j < N - 2; j++)
        {
            matrix[i][j] = 0;
        }
    }

    for (int i = 0; i < N - 2; i++) // построение матрицы, из которой будут находиться коэффициенты с
    {
        matrix[i][i] = 2 * (h[i] + h[i + 1]);
        if (i == 0)
        {
            matrix[i][i + 1] = h[i + 1];
        }
        else
        {
            if (i + 1 == N - 2)
            {
                matrix[i][i - 1] = h[i];
            }
            else
            {
                matrix[i][i + 1] = h[i + 1];
                matrix[i][i - 1] = h[i];
            }
        }
    }

    c = Jordan(matrix, gamma, N - 2); // Решение матрицы методом Жордана и получение коэффициента

    for (int i = 0; i < N - 1; i++)// переписывание коэффициентов с в другой массив с добавлением первого нулевого элемента
    {
        if (i == 0)
        {
            ci[0] = 0;
        }
        else
        {
            ci[i] = c[i - 1];
        }
    }

    for (int i = 0; i < N - 2; i++)// вычисление коэффициентов d
    {
        d[i] = ((ci[i + 1] - ci[i]) / (3.0 * h[i]));
    }
    d[N - 2] = (-1 * ci[N - 2]) / (3.0 * h[N - 2]);

    for (int i = 0; i < N - 2; i++) // вычисление коэффициентов b
    {
        b[i] = ((src_y[i + 1] - src_y[i]) / h[i]) - (((ci[i + 1] + 2 * ci[i]) / 3.0)) * h[i];
    }
    b[N - 2] = ((src_y[N - 1] - src_y[N - 2]) / h[N - 2]) - (2.0 / 3.0) * (ci[N - 2] * h[N - 2]);

    for (int i = 0; i < N - 1; i++) // вычисление коэффициентов а
    {
        a[i] = src_y[i];
    }

    for (int k = 0; k < N - 1; k++)
    {
        for (int i = 1 + k * (K + 1); i < K + 2 + k * (K + 1); i++)
        {
            new_y[i] = a[k] + (b[k] * (new_x[i] - src_x[k])) + (ci[k] * pow((new_x[i] - src_x[k]), 2)) + (d[k] * pow((new_x[i] - src_x[k]), 3)); // Кубический сплайн
        }
    }
    return new_y;
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




    int a = 0;
    int b = 0;
    double h;
    std::string path1 = "TestFile1.txt";

    double *x = new double[N];
    double *y = new double[N];

    double* x1 = new double[(N - 1) * K + N];

    h = (b - a) / (double)(N - 1);

    for (int i = 0; i < N; i++)
    {
        cout << matr[i+1][2] << endl;
        x[i] = i;
        y[i] = matr[i+1][2];
    }

    for (int i = 0; i < (N - 1) * K + N; i++)
    {
        x1[i] = (a + (i * 4.0)) / ((N - 1) * (K + 1));
    }

    std::ofstream outfile(path1);

    double* y1 = CubicInterpolation(x, y, x1);

    int k = 1;
    for (int i = 1; i < (N - 1) * K + N; i++)
    {
            outfile << x1[i] << " " << y1[i] << "\n";
    }

    outfile.close();
    return 0;
}
