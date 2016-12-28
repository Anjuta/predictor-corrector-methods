/*
Параллельная реализация предикторно-корректорного метода
*/
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


// Вычисление значений функций в системе
void func(double  *x, double  *result) {
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
void e_method(double  *x, double  *new_x, double h, int N) {
  double  func_x[N];

  func(x, func_x);
  for (int i = 0; i < N; i++)
    new_x[i] = x[i] + h * func_x[i];
}


void predictor(double  *xnp, double  *xn1c, double  *result, double h, int N) {
  double  func_x[N];
  func(xnp, func_x);
  for (int i = 0; i < N; i++)
    result[i] = xn1c[i] + 2 * h * func_x[i];
}


void corrector(double  *xn1c, double  *xnp, double  *result, double h, int N) {
  double  func_x[N], func_xc[N];
  func(xnp, func_x);
  func(xn1c, func_xc);
  for (int i = 0; i < N; i++)
    result[i] = xn1c[i] + h / 2 * (func_x[i] + func_xc[i]);
}


int main() {
  int N = 10;
  double h = 0.000001; // h = 0.00001;
  double  x[N], xn1[N], xn2[N], xnc[N], xnc1[N];

  // Заполняем начальные значения для вектора x
  for(int i = 0; i < N; i++)
    x[i] = 1 + 0.1 * i;

  // Предварительное вычисление y_n+1^p для разгона
  e_method(x, xn1, h, N);

  #pragma omp parallel num_threads (2)
  {
    int thread_num = omp_get_thread_num();
    for(double t = h; t <= 100; t += h) {
      if (thread_num == 0)
        predictor(xn1, xnc, xn2, h, N);
      else
        corrector(xnc, xn1, xnc1, h, N);
      // Ждем окончания всех вычислений
      #pragma omp barrier
      #pragma omp master
      {
        memcpy(xnc, xnc1, N * sizeof(double));
        memcpy(xn1, xn2, N * sizeof(double));
      }
    }
  }

  for(int i = 0; i < N; i++)
    printf("x[%i] = %lf\n", i, xnc1[i]);

  return 0;
}
