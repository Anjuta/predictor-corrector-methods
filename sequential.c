/*
Последовательная реализация предикторно-корректорного метода
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


// Вычисление результатов функций в системе
void func(double *x, double *result) {
  result[0] = x[0] * x[1] - x[0] + 1 / (x[6] * x[6] + x[7] * x[7] + 1);
  result[1] = -x[0] * x[1] + x[1] + 1 / (x[7] * x[7] + x[8] * x[8] + 1);
  result[2] = x[2] * x[3] - x[2] + 1 / (x[8] * x[8] + x[9] * x[9] + 1);
  result[3] = -x[2] * x[3] + x[3] + 1 / (x[9] * x[9] + x[0] * x[0] + 1);
  result[4] = x[4] * x[5] - x[4] + 1 / (x[0] * x[0] + x[1] * x[1] + 1);
  result[5] = -x[4] * x[5] + x[5] + 1 / (x[1] * x[1] + x[2] * x[2] + 1);
  result[6] = x[6] * x[7] - x[6] + 1 / (x[2] * x[2] + x[3] * x[3] + 1);
  result[7] = -x[6] * x[7] + x[7] + 1 / (x[3] * x[3] + x[4] * x[4] + 1);
  result[8] = x[8] * x[9] - x[8] + 1 / (x[4] * x[4] + x[5] * x[5] + 1);
  result[9] = -x[8] * x[9] + x[9] + 1 / (x[5] * x[5] + x[6] * x[6] + 1);
}


// Метод Эйлера для разгона
void e_method(double *x, double *new_x, double h, int N) {
  double  func_x[N];

  func(x, func_x);
  for (int i = 0; i < N; i++)
    new_x[i] = x[i] + h * func_x[i];
}


void predictor(double  *xnc, double  *xn1c, double  *result, double h, int N) {
  double  func_x[N], func_x1[N];

  func(xnc, func_x);
  func(xn1c, func_x1);
  for (int i = 0; i < N; i++)
    result[i] = xnc[i] + h / 2 * (3 * func_x[i] - func_x1[i]);
}


void corrector(double  *xnc, double  *xn1p, double  *result, double h, int N) {
  double  func_x[N], func_xp[N];

  func(xnc, func_x);
  func(xn1p, func_xp);
  for (int i = 0; i < N; i++)
    result[i] = xnc[i] + h / 2 * (func_x[i] + func_xp[i]);
}


int main() {
  int N = 10;  // Количество уравнений в системе
  double h = 0.000001;  // Шаг для метода
  double x[N], xn1[N], xn2[N], xnc[N], xnc1[N], xnc2[N];

  // Заполняем начальные значения для вектора x
  for (int i = 0; i < N; i++)
    xnc[i] = 1 + 0.1 * i;

  // Предварительное вычисление для разгона
  e_method(xnc, xnc1, h, N);

  for (double t = h; t <= 100; t += h) {
    predictor(xnc1, xnc, xn2, h, N);
    corrector(xnc1, xn2, xnc2, h, N);
    memcpy(xnc, xnc1, N * sizeof(double));
    memcpy(xnc1, xnc2, N * sizeof(double));
  }

  for(int i = 0; i < N; i++)
    printf("x[%i] = %lf\n", i, xnc1[i]);

  return 0;
}
