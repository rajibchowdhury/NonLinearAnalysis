#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "diffeqheader.h"

// Linear Oscillator
double getxdotdot(parameters* info, double time, double x, double xdot)
{
  double gamma = info->gamma;
  double omega = info->omega;
  double omega_0 = info->omega_0;
  double c = info->c;

  return (c*sin(omega*time)) - (gamma * xdot) - (omega_0 * omega_0 * x);
}

double** diffeq1(parameters* info)
{
  int i, j;

  double gamma = info->gamma;
  double omega = info->omega;
  double omega_0 = info->omega_0;
  double c = info->c;

  // Initial Conditions
  double X_0 = info->X_0;
  double X_DOT_0 = info->X_DOT_0;

  // Time Step and Time
  double delta_t = info->delta_t;
  double range = info->range; //s
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

  double xdot;
  double xdotdot;

  for (i = 1; i < len; i++)
  {
    array[0][i] = array[0][i-1] + delta_t;
    array[1][i] = array[1][i-1] + (delta_t * array[2][i-1]);
    array[2][i] = array[2][i-1] + (delta_t * getxdotdot(info, array[0][i-1], array[1][i-1], array[2][i-1]));
  }

  return array;
}
