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

template <typename T> int sign(T val) {
    return (T(0) < val) - (val < T(0));
}

double mph2ms(double mph) { return mph*0.44704; }

double norm_angle (double angle_radians) {
  if (fabs(angle_radians) > path_planning::pi()) {
      double pi2 = (path_planning::pi() * 2.);
      double correction = (floor((angle_radians + path_planning::pi()) / pi2) * pi2);
      angle_radians -= correction;
  }
  return angle_radians;
}

}

#endif /* UTILS_H */