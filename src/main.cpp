
#include <uWS/uWS.h>

#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "utils.h"
#include "map.h"
#include "trajectory_generator.h"
#include "path_planner.h"
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

//  map.Plot();


  PathPlanner path_planner(map);

// ************************  
// PLOTTING FOR THE REPORT
// ************************  
//  vector<pair<vector<vector<double>>, vector<Car>>> confs =
//  {
//  {
//  {{780.419, 18.512, 1.51745} , {6.25211, 0.596295, 0.849652} }
//  , {
//    Car({737.051, 17.3417, 0} , {1.99625, 0, 0} )
//    , Car({727.131, 20.8778, 0} , {5.98311, 0, 0} )
//    , Car({758.303, 15.2982, 0} , {10.0916, 0, 0} )
//    , Car({881.687, 17.311, 0} , {5.96817, 0, 0} )
//    , Car({755.25, 19.5138, 0} , {6.10506, 0, 0} )
//    , Car({798.226, 17.0967, 0} , {5.92148, 0, 0} )
//    , Car({789.087, 14.8169, 0} , {1.86633, 0, 0} )
//    , Car({734.343, 16.3393, 0} , {9.94021, 0, 0} )
//    , Car({777.251, 14.83, 0} , {10.0242, 0, 0} )
//    , Car({764.724, 16.0519, 0} , {2.11897, 0, 0} )
//    , Car({681.92, 19.5318, 0} , {2.08397, 0, 0} )
//    , Car({674.879, 19.3906, 0} , {6.04891, 0, 0} )
//  }
//  }
//  };
//
//  for (auto conf:confs){
//    path_planner.PlotConfiguration(conf.first[0] , conf.first[1], conf.second);
//
//  }
//
//  return 0;


  h.onMessage([&map, &path_planner](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length, uWS::OpCode opCode) {
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
          double car_yaw = j[1]["yaw"];
          double car_speed = j[1]["speed"];

          // Previous path data given to the Planner
          vector<double> previous_path_x = j[1]["previous_path_x"];
          vector<double> previous_path_y = j[1]["previous_path_y"];

          // Sensor Fusion Data, a list of all other cars on the same side of the road.
          auto sensor_fusion = j[1]["sensor_fusion"];

          // Generate trajectory
          discrete_path_2D path = path_planner.UpdatePath(
            {car_x, car_y, deg2rad(car_yaw), mph2ms(car_speed)}
            , {previous_path_x, previous_path_y}
            , sensor_fusion);

          json msgJson;

          // Define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
          msgJson["next_x"] = path.first;
          msgJson["next_y"] = path.second;

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






