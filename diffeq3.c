#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Nonlinear oscillator with driving force
// Time dependent
double getxdoubledot(double t, double x, double xdot)
{
  // Parameters
  double gamma = 1.1;
  double c = 0.0;
  double omega = 1.0;
  double a = 1.5;

  return c + (a*sin(omega*t)) - (gamma*xdot) - sin(x);
}

int main(void)
{
  int i, j;

  // Initial conditions
  double X_0 = 20.0;
  double X_DOT_0 = 0.0;
  //double X_DOUBLE_DOT_0 = 0;

  // Parameters
  double gamma = 1.1;
  double c = 0.0;
  double omega = 1.0;
  double a = 1.5;

  // Time Steps
  double delta_t = 0.05;
  double range = 50; //s
  int len = (int) range / delta_t;

  // Allocate memory for vector
  double **array = malloc(sizeof(double*) * 3);
  for (i = 0; i < 3; i++)
  {
    array[i] = malloc(sizeof(double) * len);
  }

  // Initialization of initial conditions
  array[0][0] = X_0;
  array[1][0] = X_DOT_0;
  array[2][0] = 0.0;

  double xdotdot;

  for (j = 1; j < len; j++)
  {
    array[0][j] = array[0][j-1] + (delta_t * array[1][j-1]);

    xdotdot = getxdoubledot(array[2][j-1], array[0][j-1], array[1][j-1]);
    array[1][j] = array[1][j-1] + (delta_t * xdotdot);

    array[2][j] = array[2][j-1] + delta_t;
  }

  printf("Non linear oscillator with periodic driving (gamma = %f; c = %f; a = %f, omega = %f, initial position = %f, initial velocity = %f)\n", gamma, c, a, omega, X_0, X_DOT_0);
  printf("Time\tx\tx'\t\n");
  for (j = 0; j < len; j++)
  {
    printf("%f\t%f\t%f\t\n", array[2][j], array[0][j], array[1][j]);
  }

}
