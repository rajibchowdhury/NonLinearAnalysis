#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "diffeqheader.h"

// Euler's Method on nonlinear oscillator with periodic driving
// Performance, passing struct vs double?
double getxdoubledot3(parameters* info, double t, double x, double xdot)
{
  // Parameters
  double gamma = info->gamma;
  double c = info->c;
  double omega = info->omega;
  double a = info->a;

  // Change equation here
  return c + (a*sin(omega*t)) - (gamma*xdot) - sin(x);
}


// Performs time iteration
double** diffeq3(parameters* info)
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

  double xdotdot;

  for (j = 1; j < len; j++)
  {
    array[0][j] = array[0][j-1] + delta_t;
    array[1][j] = array[1][j-1] + (delta_t * array[2][j-1]);

    xdotdot = getxdoubledot3(info, array[0][j-1], array[1][j-1], array[2][j-1]);
    array[2][j] = array[2][j-1] + (delta_t * xdotdot);

  }

  return array;
}
