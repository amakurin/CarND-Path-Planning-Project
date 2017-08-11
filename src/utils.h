#ifndef UTILS_H
#define UTILS_H

#include <math.h>
#include <vector>
#include <chrono>

namespace path_planning
{

typedef std::chrono::time_point<std::chrono::high_resolution_clock> instant;  

instant now(){
	return std::chrono::high_resolution_clock::now();
}

double get_duration(instant t1, instant t2){
	std::chrono::duration<double> t_diff = t2-t1;
	return t_diff.count();
}

constexpr double pi() { return M_PI; }

bool between(double test, double x0, double x1){
  return (test > x0 && test < x1);
}

bool not_zero(double test){
  return fabs(test)>0.00001;
}


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

std::vector<double> discrete_timeline(double dt, int steps){
  std::vector<double> timeline;
  for (int i=0; i<=steps; i++){
  	timeline.push_back(i*dt);
  }
  return timeline;
} 


double get_lane_center(int lane_index, double lane_width){
	// ugly workaround of false "Out of the lane" of simulator at few waypoints 
	if (lane_index == 2) return 9.8;
	return (lane_index + 0.5)* lane_width; 
}

int get_lane_index(double frenet_d, double lane_width, int lanes){
  int index = -1;
  while ((index<lanes) && (frenet_d > (index+1) * lane_width))
    ++index;

  if (index == lanes)
    index = -1;

  return index;
}

std::ostream& operator << (std::ostream& os, const std::vector<double>& state) {
    os << "{";
    for (int i=0; i<state.size(); i++){
    	os << state[i];
    	if (i == state.size()-1) os << "} ";
    	else os<<", ";
    }
    return os;
}



}

#endif /* UTILS_H */