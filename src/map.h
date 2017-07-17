#ifndef MAP_H
#define MAP_H

#include <fstream>
#include <vector>
#include "utils.h"

using namespace std;

namespace path_planning
{

class Map
{
private:

public:
  vector<double> maps_x_;
  vector<double> maps_y_;
  vector<double> maps_s_;
  vector<double> maps_d_x_;
  vector<double> maps_d_y_;

  Map();

  void ReadFrom(string map_file);
  int ClosestWaypoint(double x, double y);
  int NextWaypoint(double x, double y, double theta);
  vector<double> ToFrenet(double x, double y, double theta);
  vector<double> ToCartesian(double s, double d);
};

Map::Map()
{

}

void Map::ReadFrom(string map_file){
  // Load up map values for waypoint's x,y,s and d normalized normal vectors

  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

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
}

int Map::ClosestWaypoint(double x, double y)
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

int Map::NextWaypoint(double x, double y, double theta)
{
  int closest_waypoint = ClosestWaypoint(x, y);

  double map_x = maps_x_[closest_waypoint];
  double map_y = maps_y_[closest_waypoint];

  double heading = atan2((map_y - y), (map_x - x));
  double angle = abs(theta - heading);

  if(angle > pi()/4)
  {
    closest_waypoint++;
  }

  return closest_waypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> Map::ToFrenet(double x, double y, double theta)
{
  int next_wp = NextWaypoint(x, y, theta);

  int prev_wp;
  prev_wp = next_wp-1;
  if(next_wp == 0)
  {
    prev_wp  = maps_x_.size()-1;
  }

  double n_x = maps_x_[next_wp] - maps_x_[prev_wp];
  double n_y = maps_y_[next_wp] - maps_y_[prev_wp];
  double x_x = x - maps_x_[prev_wp];
  double x_y = y - maps_y_[prev_wp];

  // find the projection of x onto n
  double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
  double proj_x = proj_norm*n_x;
  double proj_y = proj_norm*n_y;

  double frenet_d = distance(x_x,x_y,proj_x,proj_y);

  //see if d value is positive or negative by comparing it to a center point

  double center_x = 1000-maps_x_[prev_wp];
  double center_y = 2000-maps_y_[prev_wp];
  double centerToPos = distance(center_x,center_y,x_x,x_y);
  double centerToRef = distance(center_x,center_y,proj_x,proj_y);

  if(centerToPos <= centerToRef)
  {
    frenet_d *= -1;
  }

  // calculate s value
  double frenet_s = 0;
  for(int i = 0; i < prev_wp; i++)
  {
    frenet_s += distance(maps_x_[i], maps_y_[i], maps_x_[i+1], maps_y_[i+1]);
  }

  frenet_s += distance(0, 0, proj_x, proj_y);

  return {frenet_s, frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> Map::ToCartesian(double s, double d)
{
  int prev_wp = -1;

  while(s > maps_s_[prev_wp+1] && (prev_wp < (int)(maps_s_.size()-1) ))
  {
    prev_wp++;
  }

  int wp2 = (prev_wp+1) % maps_x_.size();

  double heading = atan2((maps_y_[wp2] - maps_y_[prev_wp])
                        ,(maps_x_[wp2] - maps_x_[prev_wp]));
  // the x,y,s along the segment
  double seg_s = (s-maps_s_[prev_wp]);

  double seg_x = maps_x_[prev_wp] + seg_s * cos(heading);
  double seg_y = maps_y_[prev_wp] + seg_s * sin(heading);

  double perp_heading = heading-pi()/2;

  double x = seg_x + d*cos(perp_heading);
  double y = seg_y + d*sin(perp_heading);

  return {x,y};

}

}


#endif /* MAP_H */