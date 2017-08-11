#ifndef MAP_H
#define MAP_H

#include <algorithm>
#include <fstream>
#include <vector>
#include "utils.h"
#include "spline.h"

#include "matplotlibcpp.h"

using namespace std;
namespace plt = matplotlibcpp;

namespace path_planning
{

struct Map
{
  double max_s_ = 6945.554;
  double lane_width_ = 4;
  int lanes_ = 3;

  tk::spline x_s_;
  tk::spline y_s_;

  vector<double> maps_x_;
  vector<double> maps_y_;
  vector<double> maps_s_;
  vector<double> maps_d_x_;
  vector<double> maps_d_y_;

  Map(string map_file){
    ifstream in_map_(map_file.c_str(), ifstream::in);

    string line;
    while (getline(in_map_, line)) 
    {
      istringstream iss(line);
      double x;
      double y;
      float s;
      float d_x;
      float d_y;
      iss >> x;
      iss >> y;
      iss >> s;
      iss >> d_x;
      iss >> d_y;
      maps_x_.push_back(x);
      maps_y_.push_back(y);
      maps_s_.push_back(s);
      maps_d_x_.push_back(d_x);
      maps_d_y_.push_back(d_y);
    }
    int map_size = maps_s_.size();

    vector<double> ss_vals(maps_s_), xs_vals(maps_x_), ys_vals(maps_y_);
    ss_vals.push_back(max_s_);
    xs_vals.push_back(maps_x_[0]);
    ys_vals.push_back(maps_y_[0]);

    x_s_.set_points(ss_vals, xs_vals);  
    y_s_.set_points(ss_vals, ys_vals);  
  }

  void Plot(){
    vector<double> xs, ys;
    double prev_x = 0;
    double prev_y = 0;
    for (int i=0; i< int(max_s_/0.02); i++){
      double cur_s = 0.02*i;
      double cur_x = x_s_(cur_s);
      double cur_y = y_s_(cur_s);
      xs.push_back(cur_x);
      ys.push_back(cur_y);
      prev_x = cur_x;
      prev_y = cur_y;
    }

    xs.push_back(x_s_(max_s_));
    ys.push_back(y_s_(max_s_));

    plt::subplot(2, 2, 1);
    plt::title("Red = Waypoints, Black=x(s)+y(s)");
    //plt::axis("equal");
    plt::xlim(0, 2500);
    plt::plot(maps_x_, maps_y_, "r."
             , xs, ys, "k-");

    plt::subplot(2, 2, 2);
    plt::title("y(s)");
    //plt::axis("equal");
    plt::plot(maps_s_, maps_y_, "b-");

    plt::subplot(2, 2, 3);
    plt::title("s(x)");
    plt::xlim(0, 2500);
    //plt::axis("equal");
    plt::plot(maps_x_, maps_s_, "b-");

    plt::show();

  }

  int ClosestWaypoint(double x, double y)
  {
    double closest_len = 100000; //large number
    int closest_waypoint = 0;

    for(int i = 0; i < maps_x_.size(); i++)
    {
      double map_x = maps_x_[i];
      double map_y = maps_y_[i];
      double dist = distance(x, y, map_x, map_y);
      if(dist < closest_len)
      {
        closest_len = dist;
        closest_waypoint = i;
      }
    }
    return closest_waypoint;
  }

  double FrenetSCycle(double s){
    return s - int(s / max_s_) * max_s_;
  }

  vector<double> ToFrenet(double x, double y){
    int wp1 = ClosestWaypoint(x, y);

    double frenet_s = maps_s_[wp1];
    
    double ds = 0.001;
   
    double norm_distance = 1000;

    while (fabs(norm_distance) > 0.01){
      double x0 = x_s_(frenet_s);
      double y0 = y_s_(frenet_s);
      double s1 = FrenetSCycle(frenet_s+ds);
      double dx = x_s_(s1)-x0;
      double dy = y_s_(s1)-y0;

      double slope = atan2(dy, dx);

      norm_distance = (y - y0) * sin(slope) + (x - x0) * cos(slope);
    
      frenet_s = FrenetSCycle(frenet_s + norm_distance);
    }

    double x0 = x_s_(frenet_s);
    double y0 = y_s_(frenet_s);
    
    double s1 = FrenetSCycle(frenet_s+ds);
    double dx = x_s_(s1)-x0;
    double dy = y_s_(s1)-y0;

    double slope = atan2(dy, dx);

    double frenet_d = (y - y0) * cos(slope) - (x - x0) * sin(slope);

    return {frenet_s, -frenet_d};
  }

  vector<double> ToCartesian(double s, double d){
    double frenet_s = FrenetSCycle(s);
    double x = x_s_(frenet_s);
    double y = y_s_(frenet_s);

    double ds = 0.01;
    double s1 = FrenetSCycle(frenet_s+ds);
    double dx = x_s_(s1)-x;
    double dy = y_s_(s1)-y;

    double norm = atan2(dx, dy);

    x += d*cos(norm);
    y -= d*sin(norm);

    return {x, y};
  }

  double Curvature(double s){
    double frenet_s = FrenetSCycle(s);
    vector<double> dx = x_s_.derivative(frenet_s, 2);
    vector<double> dy = y_s_.derivative(frenet_s, 2);
    
    double square_sum = dx[0]*dx[0] + dy[0]*dy[0];
    double zero_thre = 1e-5; 
    if (square_sum < zero_thre) square_sum = zero_thre;

    double curvature = fabs(dx[0]*dy[1] - dy[0]*dx[1])/sqrt(square_sum*square_sum*square_sum);

    return curvature;
  }

  vector<double> ToFrenetVelocity(double vx, double vy, double frenet_s){

    double x0 = x_s_(frenet_s);
    double y0 = y_s_(frenet_s);
    double ds = 0.001;
    double s1 = FrenetSCycle(frenet_s+ds);
    double dx = x_s_(s1)-x0;
    double dy = y_s_(s1)-y0;

    double slope = atan2(dy, dx);
    
    double cos_slope = cos(slope);
    double sin_slope = sin(slope);

    double frenet_vs = vy * sin_slope + vx * cos_slope;
    double frenet_vd = vy * cos_slope - vx * sin_slope;
    return {frenet_vs, -frenet_vd};
  }


  int LaneIndex(double d){
    return get_lane_index(d, lane_width_, lanes_);
  }

  double LaneCenter(int lane_index){
    return get_lane_center(lane_index, lane_width_); 
  }

  double SubstractS(double minuend, double subtrahend){
    double s_1 = FrenetSCycle(minuend);
    double s_2 = FrenetSCycle(subtrahend);
    double dist = s_1 - s_2;
    double abs_dist = fabs(dist);
    if (abs_dist > max_s_ * 0.5) {
      double sign = 1;
      if (dist > 0) sign = -1;
      dist = sign * (max_s_ - abs_dist);
    }
    return dist;
  }

};

} /* map */




#endif /* MAP_H */