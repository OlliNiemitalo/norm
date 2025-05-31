#include <stdio.h>
#include <math.h>
#include "keyboard.h"
#include "opti.hpp"
#include <limits>

class NormProblem : public Opti::Problem {
private:
  int numParams;
  int numSamples;
  double *min;
  double *max;
  double *x;
  double *y;
  double error_multiplier;

public:

  // numParams must be odd
  NormProblem(int numParams, int numSamples, double startX = 0.001, double endX = 1.0, double error_multiplier = 1.01, double* candidate = NULL) : numParams(numParams), numSamples(numSamples), error_multiplier(error_multiplier) {
    min = new double[numParams];
    max = new double[numParams];
    x = new double[numSamples];
    if (candidate != NULL) {
      for (int i = 0; i < numParams; i++) {
        min[i] = candidate[i]-fabs(candidate[i])*(1.0/65536);
        max[i] = candidate[i]+fabs(candidate[i])*(1.0/65536);
      }
    } else {
      for (int i = 0; i < numParams; i++) {
        min[i] = -0.5;
        max[i] = 0.5;
      }
    }
    for (int i = 0; i < numSamples; i++) {
      // x[i] = startX + (endX-startX)*i/(numSamples-1);  // Uniform sampling
      x[i] = startX + (endX-startX)*(0.5 - 0.5*cos(M_PI*i/(numSamples-1)));  // Like Chebyshev nodes but including end points, to make MSE-optimimal similar to maxabs-error-optimal
    }
  }

  double *getMin() {
    return min;
  }

  double *getMax() {
    return max;
  }

  void print(double *params) {
    printf("Printout:\n");
    for (int i = 0; i < numParams; i += 3) {
      printf("(");
      for (int j = 0; j < 3; j++) {
        printf("%.20f", params[i+j]);
        if (j < 2) {
          printf(", ");
        }
      }
      printf("),\n");
    }
    printf("\n");
    for (int i = 0; i < numParams; i += 3) {
      for (int j = 0; j < 3; j++) {
        printf(j < 2? "%.20f x^%d + " : "%.20f x^%d", params[i+j], j*2+1);
      }
      printf("\n");
    }
    printf("\n");
  }

  double costFunction(double *params, double compare) {
    double maxAbsErr = 0.0;
    params[0] = fabs(params[0]);
    for (int j = 0; j*3 < numParams; j++) {
      params[j*3] = params[0];
    }
    for (int i = 0; i < numSamples; i++) {
      double y = x[i];
      double y_plus_error = x[i];
      for (int j = 0; j < numParams/3; j++) {
        double orig_y = y;
        double orig_y_plus_error = y_plus_error;
        y = params[j*3]*y + params[j*3+1]*(y*(y*y)) + params[j*3+2]*(y*(y*y)*(y*y));
        y_plus_error = params[j*3]*y_plus_error + params[j*3+1]*(y_plus_error*(y_plus_error*y_plus_error)) + params[j*3+2]*(y_plus_error*(y_plus_error*y_plus_error)*(y_plus_error*y_plus_error));
        if (y_plus_error < y) {
          std::swap(y, y_plus_error);
        }
        y_plus_error *= error_multiplier;
      }
      double absErr = std::max(fabs(y_plus_error - 1.0), fabs(y - 1.0));
      if (absErr > maxAbsErr) {
        maxAbsErr = absErr;
        if (maxAbsErr > compare) {
          return compare;
        }
      }
    }
    return maxAbsErr;
  } 

  int getNumDimensions() {
    return numParams;
  }

  ~NormProblem() {
    delete[] min;
    delete[] max;
    delete[] x;
  }
};

int main() {
  INITKEYBOARD;
  NormProblem problem(3*5, 65537, 0.001, 1.0, 1.01);  // Would also use cushion=0.029158505 but cushion is not implemented
  Opti::DERecombinator deRecombinator(0.999, 0.76);
  Opti::DE optimizer(&problem, 1000, &deRecombinator);
  for(int t = 0;; t++) {
    double bestcost = optimizer.evolve();
    if (!(t % 10000)) {
      printf("gen=%d, bestcost=%.20f, average=%.20f\n", t, bestcost, optimizer.averageCost());
      if (kbhit()) {
        printf("Parameter vector printout:\n");
        problem.print(optimizer.best());
        printf("Best cost %f\n", problem.costFunction(optimizer.best(), std::numeric_limits<double>::max()));
        if (getch() == 27) {
          break;
        }
        getch();
      }
    }
  }
  DEINITKEYBOARD;
  return 0;
}

// Compile with:
// g++ -g -O0 optimize.cpp opti.cpp
// g++ optimize.cpp opti.cpp -ffast-math -march=native -O3 -Wno-unused-result