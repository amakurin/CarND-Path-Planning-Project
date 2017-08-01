#ifndef MAP_H
#define MAP_H

#include <fstream>
#include <vector>
#include "utils.h"
#include "spline.h"

#include "matplotlibcpp.h"

using namespace std;
namespace plt = matplotlibcpp;

namespace path_planning
{

class Map
{
  public:
  struct WaypointLinkage {
    double alpha12;
    tk::spline spline;
    double ratio;
  };
private:
  double seg_fract_;
  double max_s_ = 6945.554;
public:
  vector<double> maps_x_;
  vector<double> maps_y_;
  vector<double> maps_s_;
  vector<double> maps_d_x_;
  vector<double> maps_d_y_;
  vector<WaypointLinkage> maps_linkage_;

  Map(string map_file, double curvature_segment_fraction = 0.1);
  int ClosestWaypoint(double x, double y);
  int NextWaypoint(double x, double y, double theta);
  vector<double> ToFrenetOld(double x, double y, double theta);
  vector<double> ToFrenet(double x, double y, double theta);
  vector<double> ToCartesian(double s, double d);
};

Map::Map(string map_file, double curvature_segment_fraction)
{
  seg_fract_ = curvature_segment_fraction;

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

  double tot_s = 0;
  for (int i = 0; i < map_size; i++){
    int wp0 = i-1;
    if (wp0 < 0) wp0 = map_size - 1;
    int wp1 = i;
    int wp2 = (i+1) % map_size;
    int wp3 = (i+2) % map_size;

    double s1 = distance(maps_x_[wp0],maps_y_[wp0],maps_x_[wp1],maps_y_[wp1]);
    double s2 = distance(maps_x_[wp2],maps_y_[wp2],maps_x_[wp1],maps_y_[wp1]);
    double s3 = distance(maps_x_[wp2],maps_y_[wp2],maps_x_[wp3],maps_y_[wp3]);

    double s01 = seg_fract_*s1;
    double s12 = seg_fract_*s2;

    double alpha01 = atan2((maps_y_[wp1] - maps_y_[wp0])
                     ,(maps_x_[wp1] - maps_x_[wp0]));

    double alpha12 = atan2((maps_y_[wp2] - maps_y_[wp1])
                     ,(maps_x_[wp2] - maps_x_[wp1]));

    double alpha23 = atan2((maps_y_[wp3] - maps_y_[wp2])
                     ,(maps_x_[wp3] - maps_x_[wp2]));

    double alpha = norm_angle(alpha12 - alpha01);

    double alpha3 = norm_angle(alpha23 - alpha12);

    tk::spline spline;
    
    vector<double> x, y;

    double cos_alpha = cos(alpha);
    double sin_alpha = sin(alpha);
    if (i>0){
      double x_spline = -s1*seg_fract_*1.1*cos(alpha);
      x.push_back(x_spline);
      x_spline += s1;
      double y_spline = maps_linkage_[i-1].spline(x_spline);
      y.push_back(s1 * sin_alpha - x_spline * sin_alpha + y_spline * cos_alpha);
      
       x_spline = -s1*seg_fract_*cos(alpha);
      x.push_back(x_spline);
      x_spline += s1;
       y_spline = maps_linkage_[i-1].spline(x_spline);
      y.push_back(s1 * sin_alpha - x_spline * sin_alpha + y_spline * cos_alpha);
    }else{
      x.push_back(-s1*cos_alpha);
      y.push_back(s1*sin_alpha);
    }
    x.push_back(0);
    y.push_back(0);
    x.push_back(s2);
    y.push_back(0);
    x.push_back(s2+s3*cos(alpha3));
    y.push_back(s3*sin(alpha3));

    spline.set_points(x, y);  

    double delta = 0.05;
    int cnt = (int) ceil(s2/delta);
    double arcl = 0;
    for (int k = 0; k<cnt; k++){
      double x1 = k*delta;
      double y1 = spline(x1);
      double x2 = min(x1+delta,s2);
      double y2 = spline(x2);
      double dx = x2-x1;
      double dy = y2-y1;
      arcl += sqrt(dx*dx + dy*dy);
    }
    //cout<<"i="<<i<<"\n"
    //<< "s2="<<s2
    //<<"\tarl="<<arcl
    //<<"\tdiff="<<arcl-s2
    //<<"\tmaps_s="<<maps_s_[i]<<"\ttot_s="<<tot_s
    //<<"\n";
    maps_s_[i] = tot_s;
    tot_s+=arcl;

    WaypointLinkage linkage;
    linkage.alpha12 = alpha12;
    linkage.spline = spline;
    linkage.ratio = s2/arcl;

    maps_linkage_.push_back(linkage);


    if (i <0){//6
      int steps = 250;

      std::vector<double> x_vals0(steps, 0.0);
      std::vector<double> y_vals0(steps, 0.0);
      std::vector<double> x_vals1(steps, 0.0);
      std::vector<double> y_vals1(steps, 0.0);
      std::vector<double> x_vals2;
      std::vector<double> y_vals2;

      for (int k=0; k<steps; k++){
        double incr = k*0.5-s1*cos(alpha); 
        x_vals0[k] = incr;
        y_vals0[k] = spline(x_vals0[k]);
        x_vals1[k] = incr;
        if (incr < s2) x_vals2.push_back(incr);
        if (incr < 0){
          y_vals1[k] = -incr*sin(alpha);
          double prev_x = s1+x_vals2[k]*cos(alpha);
          double prev_y = maps_linkage_[i-1].spline(prev_x);
          y_vals2.push_back(s1*sin(alpha) -prev_x * sin(alpha) + prev_y * cos(alpha));
        }
        else if (incr > s2){
          y_vals1[k] = (incr-s2)*sin(alpha3);
        }
        else 
        {
          y_vals1[k] = 0;
          double prev_x = s1+x_vals2[k]*cos(alpha);
          double prev_y = maps_linkage_[i-1].spline(prev_x);
          y_vals2.push_back(s1*sin(alpha)-prev_x * sin(alpha) + prev_y * cos(alpha));
        }
      }

      cout<< "WP "<<i<<"\t alpha="<<alpha <<'\n'
      <<"a01="<<alpha01<<"\ta12="<<alpha12<<'\n'
      <<"s1="<<s1<<"\ts2="<<s2<<'\n'
      <<"spline(x[1])="<<spline(x[1])<<"\tspline(x[2])="<<spline(x[2])<<'\n'
      <<"s1 + s12*cos(alpha)="<<x[3]<<"\ts12*sin(alpha)="<<y[3]<<"\tspline="<<spline(x[3])<<'\n'
      <<"----\n";
      
      plt::subplot(1, 1, 1);
      plt::title("y(x) Red - spline, Blue - linear, Green - prev spline");
      plt::axis("equal");
      plt::plot(x_vals0, y_vals0, "r-"
              , x_vals1, y_vals1, "b-"
              , x_vals2, y_vals2, "g-"
              );
      plt::show();
      //break;
    }
  }
  max_s_ = tot_s;
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
  double angle = abs(norm_angle(theta - heading));

  if(angle > pi()/4)
  {
    closest_waypoint++;
  }

  return closest_waypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> Map::ToFrenetOld(double x, double y, double theta)
{
  int next_wp = NextWaypoint(x, y, theta);

  int wp1;
  wp1 = next_wp-1;
  if(next_wp == 0)
  {
    wp1  = maps_x_.size()-1;
  }

  double n_x = maps_x_[next_wp] - maps_x_[wp1];
  double n_y = maps_y_[next_wp] - maps_y_[wp1];
  double x_x = x - maps_x_[wp1];
  double x_y = y - maps_y_[wp1];

  // find the projection of x onto n
  double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
  double proj_x = proj_norm*n_x;
  double proj_y = proj_norm*n_y;

  double frenet_d = distance(x_x,x_y,proj_x,proj_y);

  //see if d value is positive or negative by comparing it to a center point

  double center_x = 1000-maps_x_[wp1];
  double center_y = 2000-maps_y_[wp1];
  double centerToPos = distance(center_x,center_y,x_x,x_y);
  double centerToRef = distance(center_x,center_y,proj_x,proj_y);

  if(centerToPos <= centerToRef)
  {
    frenet_d *= -1;
  }

  // calculate s value
  double frenet_s = 0;
  for(int i = 0; i < wp1; i++)
  {
    frenet_s += distance(maps_x_[i], maps_y_[i], maps_x_[i+1], maps_y_[i+1]);
  }

  frenet_s += distance(0, 0, proj_x, proj_y);

  return {frenet_s, frenet_d};

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> Map::ToFrenet(double x, double y, double theta)
{
  int wp2 = NextWaypoint(x, y, theta);

  int wp1;
  wp1 = wp2-1;
  if (wp2 == 0) wp1 = maps_x_.size()-1;
  
  Map::WaypointLinkage link = maps_linkage_[wp1];

  double x_diff = x - maps_x_[wp1];
  double y_diff = y - maps_y_[wp1];

  double sin_alpha12 = sin(link.alpha12);
  double cos_alpha12 = cos(link.alpha12);

  double seg_x = y_diff * sin_alpha12 + x_diff * cos_alpha12;
  double seg_y = y_diff * cos_alpha12 - x_diff * sin_alpha12;

  double dx = 0.001;
  
  double x0 = seg_x;
  double y0 = link.spline(x0);
  double norm_distance = 100;
  while (fabs(norm_distance)>dx/10){
    double dy = link.spline(x0+dx) - y0;

    double m = dy/dx;
    double b = y0 + x0/m;

    norm_distance = (seg_x+m*seg_y-m*b)/sqrt(1+m*m);

    //cout<< "wp1="<<wp1
    //<< "\nx0="<<x0<<"\ty0="<<y0
    //<< "\nnd="<<norm_distance
    //<<"\n";
    x0 += norm_distance*cos(atan2(dy,dx));
    y0 = link.spline(x0);
  }

  double frenet_s = maps_s_[wp1] + x0;
  double frenet_d = distance(seg_x, seg_y, x0, y0)*sign(seg_y);

  return {frenet_s, -frenet_d};
}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> Map::ToCartesian(double s, double d)
{
  double s_ = s - floor(s/max_s_)*max_s_ + 0.00001;

  int wp1 = -1;

  while(s_ > maps_s_[wp1+1] && (wp1 < (int)(maps_s_.size()-1) ))
  {
    wp1++;
  }

  double seg_s = s_ - maps_s_[wp1];

  WaypointLinkage link = maps_linkage_[wp1];

  double seg_x = seg_s;
  double seg_y = link.spline(seg_x);
  double delta = 0.001;
  double spline_norm = atan2(link.spline(seg_x+delta) - seg_y, delta) - pi()/2;
  
  seg_x += d*cos(spline_norm);
  seg_y += d*sin(spline_norm);

  double sin_head = sin(link.alpha12);
  double cos_head = cos(link.alpha12);

  double x = maps_x_[wp1] + seg_x * cos_head - seg_y * sin_head;
  double y = maps_y_[wp1] + seg_x * sin_head + seg_y * cos_head;

  //cout<< wp1 << "----\n"
  //  <<"x="<<x <<"\ty="<<y<<'\n';
  return {x,y}; 
}

}

#endif /* MAP_H */