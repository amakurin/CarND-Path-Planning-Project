#ifndef UTILS_H
#define UTILS_H

#include <math.h>


namespace path_planning
{

constexpr double pi() { return M_PI; }

double deg2rad(double x) { return x * pi() / 180; }

double rad2deg(double x) { return x * 180 / pi(); }

double distance(double x1, double y1, double x2, double y2)
{
  double x_diff = (x2-x1);  
  double y_diff = (y2-y1);  
  return sqrt(x_diff * x_diff + y_diff * y_diff);
}

}

#endif /* UTILS_H */