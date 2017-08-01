
#include <uWS/uWS.h>

#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "utils.h"
#include "map.h"
#include "vehicle.h"
#include "trajectory_generator.h"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

using namespace std;
using namespace path_planning;
// for convenience
using json = nlohmann::json;

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) 
  {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) 
  {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

int main() {
  uWS::Hub h;

  Map map("../data/highway_map.csv");
 
  double v_max = mph2ms(45);
  double a_max = 10; 
  double j_max = 10;

  auto sdpoint = map.ToFrenet(909.48, 1128.67, 0);

  double s1 = sdpoint[0];  
  double sv1 = 0;  
  double sa1 = 0;  

  int time_horizon = 200;
  double s2 = s1+10;  
  double sv2 = v_max;  
  double sa2 = 0;  

  double dt = 0.02;
  
  vector<vector<double>> s_traj;
  for (int i = 0; i < 1; i++){

    vector<double> is = {s1, sv1, sa1};
    vector<double> fs = {s2, sv2, sa2};
    vector<double> lims = {a_max, j_max};
    s_traj = generate_trajectory_1Dq(is, fs, lims, dt, time_horizon);

    vector<double> s_vals;
    vector<double> sd_vals;
    vector<double> sdd_vals;
    vector<double> sddd_vals;

    vector<double> t_vals;
    for (int i =0; i < s_traj.size(); i++){
      cout<<i<<"\tt= "<<i*dt<<"\ts= "<<s_traj[i][0]<<"\tv= "<<s_traj[i][1]
      <<"\ta= "<<s_traj[i][2]<<"\tj= "<<s_traj[i][3]<<"\n";
      s_vals.push_back(s_traj[i][0]);
      sd_vals.push_back(s_traj[i][1]);
      sdd_vals.push_back(s_traj[i][2]);
      sddd_vals.push_back(s_traj[i][3]);

      t_vals.push_back(i*dt);
    }

    plt::subplot(1, 1, 1);
    plt::title("y(x) Red - traj, Green - v, Blue - a");
    //plt::axis("equal");
    plt::plot(t_vals, s_vals, "r-"
            , t_vals, sd_vals, "g-"
            , t_vals, sdd_vals, "b-"
            , t_vals, sddd_vals, "m-");

//    plt::subplot(2, 1, 2);
//    //plt::axis("equal");
//    plt::plot(t_vals, d_vals, "r-"
//            , t_vals, dd_vals, "g-"
//            , t_vals, ddd_vals, "b-");
    plt::show();
  }
  return 0;

  Vehicle car;

  h.onMessage([&map, &car](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length, uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
          // Main car's localization Data
          double car_x = j[1]["x"];
          double car_y = j[1]["y"];
          double car_s = j[1]["s"];
          double car_d = j[1]["d"];
          double car_yaw = j[1]["yaw"];
          double car_speed = j[1]["speed"];

          // Previous path data given to the Planner
          auto previous_path_x = j[1]["previous_path_x"];
          auto previous_path_y = j[1]["previous_path_y"];
          // Previous path's end s and d values 
          double end_path_s = j[1]["end_path_s"];
          double end_path_d = j[1]["end_path_d"];

          // Sensor Fusion Data, a list of all other cars on the same side of the road.
          auto sensor_fusion = j[1]["sensor_fusion"];

          json msgJson;

          car_yaw = deg2rad(car_yaw);

          int next_wp = map.NextWaypoint(car_x, car_y, car_yaw);
          double next_x = map.maps_x_[next_wp] + map.maps_d_x_[next_wp]*6;
          double next_y = map.maps_y_[next_wp] + map.maps_d_y_[next_wp]*6;
          double next_s = map.maps_s_[next_wp];

          double dist = distance(car_x, car_y, next_x, next_y); 

          vector<double> acc = car.Update(mph2ms(car_speed));

          cout //<< "Current XY:\t"<< car_x<<",\t" << car_y
          //<< "\nNext WP:\t"<<next_x<<",\t" << next_y
          //<< "\nDistance:\t"<< dist
          //<< "\nCurrent SD:\t"<< car_s<<",\t" << car_d
          //<< "\nNext WP:\t"<<next_s
          << "\nCurrent Heading:\t"<<car_yaw
          << "\nCurrent SPEED:\t"<<car_speed
          //<< "\nCurrent ACCEL:\t"<<acc[1]
          //<< "\nCurrent LAT:\t"<<acc[0]
          //<< "\nPREV PATH SIZE:\t"<<previous_path_x.size()
          <<"\n----\n";

          auto sdpoint = map.ToFrenet(car_x, car_y, car_yaw);
          auto point = map.ToCartesian(sdpoint[0], sdpoint[1]);
          //cout << "CALC XY:\t"<< point[0]<<",\t" << point[1]
          //<<"\n----\n";


          double speed_lim_mph = 45.0;
          double speed_lim_ms = mph2ms(speed_lim_mph);

          int horizon = 50;
          int n = 10;

          //cout << "prev_size="<<previous_path_x.size()<< "\ttraj_size="<<car.last_trajectory_s_.size()<<'\n';
          int diff = car.last_trajectory_s_.size() - previous_path_x.size();

          if (diff>0){
            auto begin = car.last_trajectory_s_.begin();
            car.last_trajectory_s_.erase(begin, begin + diff);
            begin = car.last_trajectory_d_.begin();
            car.last_trajectory_d_.erase(begin, begin + diff);
          }
          
          if (n < car.last_trajectory_s_.size()){
            car.last_trajectory_s_.erase(car.last_trajectory_s_.begin()+n, car.last_trajectory_s_.end());
            car.last_trajectory_d_.erase(car.last_trajectory_d_.begin()+n, car.last_trajectory_d_.end());
          }

          int traj_size = car.last_trajectory_s_.size();
          vector<double> init_s = {sdpoint[0], 0, 0};
          vector<double> init_d = {sdpoint[1], 0, 0};
          if (traj_size > 0){
            init_s = *(car.last_trajectory_s_.end()-1);
            init_d = *(car.last_trajectory_d_.end()-1);
          }

          double dt = 0.02;

          int time_horizon = horizon - traj_size;
          vector<double> fs = {init_s[0]+speed_lim_ms*dt*time_horizon, speed_lim_ms, 0};
          vector<vector<double>> s_traj = generate_trajectory_1D(init_s, fs, dt*time_horizon, dt);

          vector<double> fd = {6, 0, 0};
          vector<vector<double>> d_traj = generate_trajectory_1D(init_d, fd, dt*time_horizon, dt);

          for(int i = 1; i < s_traj.size(); i++)
          {
            car.last_trajectory_s_.push_back(s_traj[i]);
            car.last_trajectory_d_.push_back(d_traj[i]);
            //cout<<"GEN i="<<i<<"\ts="<<s_traj[i][0]<<"\td="<<d_traj[i][0]<<"\n";
          }
          //for(int i = 0; i < horizon - traj_size; i++)
          //{
          //  vector<double> s = {init_s[0] + 0.41 * (i + 1), 0, 0};
          //  car.last_trajectory_s_.push_back(s);
          //}

          vector<double> next_x_vals;
          vector<double> next_y_vals;

          for (int i = 0; i < car.last_trajectory_s_.size(); i++){
            auto point = map.ToCartesian(car.last_trajectory_s_[i][0], car.last_trajectory_d_[i][0]);
            //cout<<"TRAJ XY i="<<i<<"\tx="<<point[0]<<"\ty="<<point[1]<<"\n";
            next_x_vals.push_back(point[0]);
            next_y_vals.push_back(point[1]);
          }

          // TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;

          auto msg = "42[\"control\","+ msgJson.dump()+"]";

          //this_thread::sleep_for(chrono::milliseconds(1000));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    // no need to close here Each close raises onDisconnection -> segmentation fault
    //ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}






