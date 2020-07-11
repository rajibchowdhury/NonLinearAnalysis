#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "diffeqheader.h"

double func(parameters* info, double t, double x)
{
  double val;

  double gamma = info->gamma;
  double c = info->c;

  val = (double) (c - sin(x)) / gamma;
  //printf("val = %f\n", val);
  return val;
}

double rg2(parameters* info, double t, double x)
{
  double k1, k2, k3, k4;
  double sol;

  double delta_t = info->delta_t;

  k1 = delta_t * func(info, t, x);
  //printf("k1 = %f\n", k1);
  k2 = delta_t * func(info, t + (0.5*delta_t), x + (0.5*k1));

  k3 = delta_t * func(info, t + (0.5*delta_t), x + (0.5*k2));

  k4 = delta_t * func(info, t + delta_t, x + k3);

  sol = (double) ((k1 + (2*k2) + (2*k3) + k4) / 6);
  //printf("k = %f\n", sol);
  return sol;
}



double** rungekutta2(parameters* info)
{
  int i, j;

  // Initial conditions
  double X_0 = info->X_0;
  double X_DOT_0 = info->X_DOT_0;

  // Parameters
  double gamma = info->gamma;
  double c = info->c;
  double omega = info->omega;
  double a = info->a;

  // Time Steps
  double delta_t = info->delta_t;
  double range = info->range;
  int len = (int) range / delta_t;

  // Allocate memory for vector
  double **array = malloc(sizeof(double*) * 3);
  for (i = 0; i < 3; i++)
  {
    array[i] = malloc(sizeof(double) * len);
  }

  // Initialization of initial conditions
  array[0][0] = 0.0;
  array[1][0] = X_0;
  array[2][0] = X_DOT_0;

  // Runge Kutta Routine
  for (i = 1; i < len; i++)
  {
    array[0][i] = array[0][i-1] + delta_t;
    array[1][i] = array[1][i-1] + rg2(info, array[0][i-1], array[1][i-1]);
    //array[2][i] = array[2][i-1] + rgx3(info, array[0][i-1], array[1][i-1], array[2][i-1]);
  }

  return array;
}
