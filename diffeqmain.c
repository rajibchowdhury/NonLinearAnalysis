#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "diffeqheader.h"

void printParameters(parameters* parameter)
{
  printf("Gamma = %f, C = %f, Omega = %f, Omega_0 = %f, A = %f\n", parameter->gamma, parameter->c, parameter->omega, parameter->omega_0, parameter->a);
  printf("Initial position = %f, initial velocity = %f, Range = %f, Time Step = %f\n", parameter->X_0, parameter->X_DOT_0, parameter->range, parameter->delta_t);
}

void printArray(double** array, double len)
{
  int i;

  for (i = 0; i < len; i++)
  {
    printf("%f\t%f\n", array[0][i], array[2][i]);
  }
}


parameters* setdiffeq3(double gamma, double c, double omega, double omega_0, double a, double x_0, double x_dot_0, double range, double delta_t)
{
  parameters* info = malloc(sizeof(parameters));

  info->gamma = gamma;
  info->c = c;
  info->omega = omega;
  info->omega_0 = omega_0;
  info->a = a;

  info->X_0 = x_0;
  info->X_DOT_0 = x_dot_0;

  info->range = range;
  info->delta_t = delta_t;

  return info;
}

int main(void)
{
  int i, j;
  double **array1, **array2, **array3, **array4;
  double ***array;
  parameters* information;

  // Change parameters here
  double gamma = 0.7;
  double c = 0.0;
  double omega = 0.1;
  double omega_0 = 0.0;
  double a = 0.4;

  double x_0 = 0.0;
  double x_dot_0 = 0.0;

  double range = 2000;
  double delta_t = 0.001;

  int len = (int) range / delta_t;

  information = setdiffeq3(gamma, c, omega, omega_0, a, x_0, x_dot_0, range, delta_t);

  double precision;
  double c1, c2;
  int range2;

  c1 = 0;
  // Error when c2 isnt integer?
  c2 = c1 + 4.0;
  precision = 0.001;
  range2 = (int) (c2 - c1) / precision;

  double** i_v = ivcharacteristic(information, c1, c2, precision);

  for (i = 0; i < range2; i++)
  {
    printf("%f\t%f\n", i_v[0][i], i_v[1][i]);
  }
  // printf("%f\n", averageVelocity(practicefunction(information), information));

  return 0;
}
