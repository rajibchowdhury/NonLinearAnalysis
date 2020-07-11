#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Nonlinear oscillator
// Time independent
double getxdotdot(double x, double xdot)
{
  double gamma = 1.0;
  double c = 0.5;

  return c - (gamma * xdot) - sin(x);
}

int main(void)
{
  double **array = malloc(sizeof(double*) * 3);
  double range = 50;
  double delta_t = 0.025;
  int len = (int) range / delta_t;
  int i;
  int j;

  // Initial Conditions
  double X_0 = 20;
  double X_DOT_0 = 0;

  // constants
  double gamma = 1.0;
  double c = 1.0;

  for (i = 0; i < 3; i++)
  {
    array[i] = malloc(sizeof(double) * len);
  }

  // Time Array
  array[0][0] = 0;
  for (j = 1; j < len; j++)
  {
    array[0][j] = array[0][j-1] + delta_t;
  }

  array[0][0] = 0.0;
  array[1][0] = X_0;
  array[2][0] = X_DOT_0;

  double xdot;
  double xdotdot;

  for (j = 1; j < len; j++)
  {
    xdot = array[2][j-1];
    array[1][j] = array[1][j-1] + (xdot * delta_t);

    xdotdot = getxdotdot(array[1][j-1], array[2][j-1]);
    array[2][j] = array[2][j-1] + (xdotdot * delta_t);
  }

  printf("Nonlinear oscillator (gamma = %f; c = %f; range = %f, delta t = %f)\n", gamma, c, range, delta_t);
  printf("Time\tx\tx'\t\n");
  for (j = 1; j < len; j++)
  {
    printf("%f\t%f\t%f\t\n", array[0][j], array[1][j], array[2][j]);
  }

  return 0;
}
