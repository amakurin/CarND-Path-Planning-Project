#ifndef PATH_PLANNER_H
#define PATH_PLANNER_H
#include <vector>
#include <map>
#include <utility>
#include <chrono>
#include <algorithm>
#include "trajectory_generator.h"
#include "map.h"
#include "util.h"
#include "matplotlibcpp.h"
#include "aabb.h"

namespace plt = matplotlibcpp;
using namespace std;
using namespace path_planning;

namespace path_planning
{

// Neighboring car data container
struct Car{
  // state in S
  vector<double> s_state_;
  // state in D
  vector<double> d_state_;
  // car lane index
  int lane_index_ = -1;
  // distance to me
  double s_distance_ = 1000;
  // predicted trajectory
  continuous_path_2D prediction_;
  // AABB tree build by trajectory
  AABBTree aabb_;

  // flag for plotting
  bool is_referenced_ = 0;

  Car(vector<double> s_state, vector<double> d_state){
    s_state_ = s_state;
    d_state_ = d_state;
  }
};


struct PathPlanner{
  // timing options
  double planning_horizon_ = 5;
  double dt_ = 0.02;
  int planning_time_steps_ = planning_horizon_/dt_;
  int derivation_steps_ = 5;
  
  //kinematic limits
  double a_max_ = 10;
  double j_max_ = 10;
  double v_max_ = mph2ms(46);

  //cars dimensions, used in collision checks
  double cars_length_ = 8;
  double cars_width_ = 2.5;

  //safety options
  double safety_distance_ = cars_length_+2;
  double safety_time_ = 0.5;

  //current path
  vector<PathVariation1D> current_path_var_;
  //used steps of current path
  int used_steps_ = 0;
  //s_shift of current path (for local coordinates transformation)
  double s_shift_ = 0;

  Map map_;

  PathPlanner(Map road_map): map_(road_map){
  }

  double CalcSafetyRegion(double s_velocity){
    return safety_distance_ + s_velocity * safety_time_;
  }

  int GetPastSteps(discrete_path_2D path_remainder){
    int result = 0;
    if (current_path_var_.size()>0)
      result = planning_time_steps_ - path_remainder.first.size();
    return result; 
  }

  // wraps sensor_fusion data to Car containers, predicts states of cars in current_time_step
  vector<Car> PredictCars(
    const vector<vector<double>> & sensor_fusion
    , int current_time_step){
    vector<Car> cars;

    double current_time = dt_ * current_time_step;
    for (int i = 0; i < sensor_fusion.size(); ++i)
    {
      auto car_info = sensor_fusion[i];
      double x = car_info[1];
      double y = car_info[2];
      double vx = car_info[3];
      double vy = car_info[4];
      // transform to Frenet Frame (cant use provided S and D)
      auto car_sdpoint = map_.ToFrenet(x, y);
      
      double s = car_sdpoint[0];
      double d = car_sdpoint[1];

      // Calculate velocity in Frenet Frame
      vector<double> frenet_v = map_.ToFrenetVelocity(vx, vy, s);

      // Init State
      vector<double> s_state = {s, frenet_v[0], 0};
      vector<double> d_state = {d, 0, 0};

      // Predict state to current_time_step
      s_state[0] =  map_.FrenetSCycle(s_state[0] + (s_state[1]  + 0.5 * s_state[2] * current_time) * current_time);
      
      s_state[1] += s_state[2] * current_time;

      d_state[0] += (frenet_v[1]  + 0.5 * d_state[2] * current_time) * current_time;
      d_state[1] += d_state[2] * current_time;

      Car car(s_state, d_state);
      cars.push_back(car);
    }
    return cars;
  }

  // Filters cars in s_radius around of s_state, predicts their trajectories and builds AABB trees
  vector<Car> GetWatchableCars(
    vector<double> & s_state
    , vector<double> & d_state
    , vector<Car> & cars
    , double s_radius){
    vector<Car> watchable_cars;
    for (auto car : cars){
      car.lane_index_ = map_.LaneIndex(car.d_state_[0]);
      car.s_distance_ = map_.SubstractS(car.s_state_[0], s_state[0]);
      if (car.lane_index_ != -1 && fabs(car.s_distance_)< s_radius){
        car.prediction_ = {keep_state_trajectory_1D(car.s_state_)
                          , keep_state_trajectory_1D(car.d_state_)};
        auto s_path = discretize_path(car.prediction_.first, dt_, planning_time_steps_);
        auto d_path = discretize_path(car.prediction_.second, dt_, planning_time_steps_);
        car.aabb_ = AABBTree::Build(s_path, d_path, cars_length_, cars_width_);

        watchable_cars.push_back(car);
      }
    }
    return watchable_cars;
  }

  // Finds closest leading car index for lane, or -1 if not found
  int ClosestCarAhead(int lane_index, const vector<Car> & watchable_cars){
    int closest_car_index = -1;
    double min_distance = 1000;
    for (int i = 0; i < watchable_cars.size(); i++){
      Car car = watchable_cars[i];
      if (car.lane_index_ == lane_index 
        && car.s_distance_ > 0
        && car.s_distance_ < min_distance){
        min_distance = car.s_distance_;
        closest_car_index = i;
      }
    }
    return closest_car_index;
  }

  // copies 1D pathes for lane from source to target
  void CopyLaneVariations(int lane, vector<PathVariation1D> & source, vector<PathVariation1D> & target){
    for (int i=0; i< source.size(); i++){
      auto var = source[i];
      var.lane_index_ = lane;
      target.push_back(var);
    }
  }

  // Generates S(t) and D(t) sets of 1D trajectories dependong on 
  void BuildVariations(
  // my s state s, ds, dds
   const vector<double>& s_state
  // my d state d, dd, ddd
  , const vector<double>& d_state
  // neighboring cars states
  , vector<Car> & watchable_cars
  // resulting S(t) set
  , vector<PathVariation1D> & s_variations
  // resulting D(t) set
  , vector<PathVariation1D> & d_variations){

    int my_lane = map_.LaneIndex(d_state[0]);
    double safety_region = CalcSafetyRegion(s_state[1]);

    // safety S trajectories to not collide with leading car
    vector<PathVariation1D> my_lane_common_s_vars;

    // inspect my lane and surrounding lanes
    vector<int> lanes = {my_lane};
    if (my_lane > 0) lanes.push_back(my_lane - 1);
    if (my_lane < map_.lanes_ - 1) lanes.push_back(my_lane + 1);

    for (auto lane : lanes){
      bool has_leading = false;
      // if has leading car
      int closest_car_ahead_index = ClosestCarAhead(lane, watchable_cars);
      if (closest_car_ahead_index > -1){
        auto closest_car = (watchable_cars.begin()+closest_car_ahead_index);
        // close enough to reference it
        if (closest_car->s_distance_ <  1.5 * safety_region){
          closest_car->is_referenced_ = true;
          has_leading = true; 

          vector<PathVariation1D> vars;
          // if it is my lane      
          if (lane == my_lane){
            // reference it
            vars = variate_1D(s_state, {ReferenceState::FollowState(closest_car->s_state_, safety_distance_, safety_time_)});  
            // and remember if too close
            if (closest_car->s_distance_ <  1.5 * safety_distance_){
              CopyLaneVariations(lane, vars, my_lane_common_s_vars);
            }
          }  
          // if it is not my line
          else{
            // reference it if is not to close
            if (closest_car->s_distance_ > safety_distance_){
              vars = variate_1D(s_state, {ReferenceState::FollowState(closest_car->s_state_, safety_distance_, safety_time_)});
            }
          } 
          // copy variations to main S set
          CopyLaneVariations(lane, vars, s_variations);
        }
      }

      // if doesn't have leading
      if (!has_leading){
        // reference max_velocity, zero_acceleration
        auto vars = variate_1D(s_state, {ReferenceState::KeepVelocityState(v_max_)});
        CopyLaneVariations(lane, vars, s_variations);
      }
      if (lane != my_lane){
        // add safety options if exitst
        CopyLaneVariations(lane, my_lane_common_s_vars, s_variations);
      }

      // D variations for lane
      vector<PathVariation1D> lane_d_variations = variate_1D(d_state, {ReferenceState::KeepLaneState(lane, map_.lane_width_, map_.lanes_)});
      CopyLaneVariations(lane, lane_d_variations, d_variations);

    }
  }

  // combines feasible of S set with each feasible of D set
  void CombineFeasible(
    vector<PathVariation1D> & s_variations
    , vector<PathVariation1D> & d_variations
    , vector<PathVariation2D> & result){

    for (int j=0; j< d_variations.size(); j++){
      if (d_variations[j].IsFeasible(j_max_, a_max_, v_max_)){
        for (int i=0; i< s_variations.size(); i++){         
          if (s_variations[i].lane_index_ == d_variations[j].lane_index_){
            if (s_variations[i].IsFeasible(j_max_, a_max_, v_max_)){
              result.push_back(PathVariation2D(&s_variations[i], &d_variations[j]));
            }
          }
        }
      }
    }
  }

  // checks for collision of each proposed 2D path with each os neighboring cars
  // used only for plotting
  void CheckCollisions(
    vector<PathVariation2D> & vars_2D
    , const vector<Car> & watchable_cars
    , int my_lane){

    for (auto car : watchable_cars){
      if (car.lane_index_ != my_lane || car.s_distance_ > 0){
        for (int i=0; i< vars_2D.size(); i++){
          auto s_discrete = vars_2D[i].s_variation_->DiscretePath(dt_, planning_time_steps_);
          auto d_discrete = vars_2D[i].d_variation_->DiscretePath(dt_, planning_time_steps_);
          int index = car.aabb_.CollisionIndex(s_discrete, d_discrete, cars_length_, cars_width_); 
          if (index != -1){
            vars_2D[i].is_collided_ = true;
          }
        }
      }
    }
  }

  // finds index of first, collision free path
  int OptimalCollisionFree(
    vector<PathVariation2D> & vars_2D
    , const vector<double> & s_state
    , const vector<double> & d_state
    , const vector<Car> & watchable_cars
    , int my_lane){

    // sort by cost
    PathVariation2DComparator comparator(s_state, d_state, planning_horizon_, v_max_);
    sort(vars_2D.begin(), vars_2D.end(), comparator);

    // check collisions for each path starting from most optimal 
    int found = -1; 
    for (int i=0; i< vars_2D.size(); i++){
      for (auto car : watchable_cars){
        if (car.lane_index_ != my_lane || car.s_distance_ > 0){
          auto s_discrete = vars_2D[i].s_variation_->DiscretePath(dt_, planning_time_steps_);
          auto d_discrete = vars_2D[i].d_variation_->DiscretePath(dt_, planning_time_steps_);
          int index = car.aabb_.CollisionIndex(s_discrete, d_discrete, cars_length_, cars_width_); 
          if (index != -1){
            vars_2D[i].is_collided_ = true;
          }
        }
      }
      if (!vars_2D[i].is_collided_){
        found = i;
        break;
      }
    }
    return found;
  }

// Main path generation routine
  discrete_path_2D UpdatePath(
    // current car state from simulator: x, y, yaw, speed
    vector<double> car_state
    // path remainder from simulator
    , discrete_path_2D path_remainder
    // sensor fusion from simulator
    , vector<vector<double>> sensor_fusion
    ){
//OutputFusion(sensor_fusion);

    // derive path to compensate latency 
    discrete_path_2D path = derive_discrete_path(path_remainder, derivation_steps_);

    // calculate instant of initial state
    int derived_steps = path.first.size();
    used_steps_ += GetPastSteps(path_remainder);
    int state_step = derived_steps + used_steps_;
    double state_instant = dt_ * state_step;

    // find initial state 
    vector<double> s_state;
    vector<double> d_state;
    double s_original = 0;
    // from current path, if exitst
    if (current_path_var_.size()>0) {
      s_state = current_path_var_[0].EvalAllAt(state_instant);
      d_state = current_path_var_[1].EvalAllAt(state_instant);
      s_original = s_state[0] + s_shift_;
    }
    // or from current car state 
    else {
      auto sdpoint = map_.ToFrenet(car_state[0], car_state[1]); 
      s_state = {sdpoint[0], car_state[3], 0};
      d_state = {sdpoint[1], 0, 0};
      s_shift_ = s_state[0];
      s_original = s_state[0];
    }
    
    // find neighboring cars states
    auto cars = PredictCars(sensor_fusion, derived_steps);

    // transform all to local coordinates
    for (int i = 0; i < cars.size(); i++){
      cars[i].s_state_[0] = map_.SubstractS(cars[i].s_state_[0], s_original);
    }
    s_state[0] = 0;

//    OutputConfiguration(s_state, d_state, cars);

    // calculate my carrent lane
    int my_lane = map_.LaneIndex(d_state[0]);

    // calculate radius to look for reference cars
    double s_radius = v_max_ * planning_horizon_;

    // collect cars with build trajectories and AABB trees
    vector<Car> watchable_cars = GetWatchableCars(s_state, d_state, cars, s_radius);

    // propose 1D variations of path
    vector<PathVariation1D> s_variations;
    vector<PathVariation1D> d_variations;
    BuildVariations(s_state, d_state, watchable_cars, s_variations, d_variations);

    // select feasible 1D pathes and combine them to 2D 
    vector<PathVariation2D> vars_2D;
    CombineFeasible(s_variations, d_variations, vars_2D);
    
    // find new optimal collision free trajectory 
    int optimal_path_index = OptimalCollisionFree(vars_2D, s_state, d_state, watchable_cars, my_lane);

    // if found - set as current
    if (current_path_var_.size()==0 || optimal_path_index != -1){
      PathVariation2D  optimal_variation = vars_2D[optimal_path_index];
      current_path_var_ = 
        {*(optimal_variation.s_variation_), *(optimal_variation.d_variation_)};
      used_steps_ = -derived_steps;
      s_shift_ = s_original;
    }

    // discretize current path
    int start_step = derived_steps + used_steps_;
    auto s_path = current_path_var_[0].DiscretePath(dt_, start_step + planning_time_steps_);
    auto d_path = current_path_var_[1].DiscretePath(dt_, start_step + planning_time_steps_);

    // fill in resulting path 
    for (int i = start_step+1; i < s_path.size(); i++){
      auto point = map_.ToCartesian(s_shift_ + s_path[i], d_path[i]);
      path.first.push_back(point[0]);
      path.second.push_back(point[1]);
      // up to planning horizon
      if (path.first.size() >= planning_time_steps_) break;
    }
    return path;
  }

  //********************
  // PLOT UTILS
  //********************
  void OutputFusion(vector<vector<double>> sensor_fusion){
    cout<< "Fusion:\n";
  
    for (int i=0; i< sensor_fusion.size(); i++){
      auto car_info = sensor_fusion[i];
      double cr_x = car_info[1];
      double cr_y = car_info[2];
      double cr_vx = car_info[3];
      double cr_vy = car_info[4];
      //auto car_sdpoint = map_.ToFrenet(cr_x, cr_y);
      cout<< "--- "<< car_info[0]
      <<"\tx="<<cr_x<<"\ty="<<cr_y<<"\tvx="<<cr_vx<<"\tvy="<<cr_vy
      <<"\nsreal="<<car_info[5]<<"\tdreal="<<car_info[6]
      //<<"\ns="<<car_sdpoint[0]<<"\td="<<car_sdpoint[1]
      <<"\n";
    }
  }


  void Plot1DVariations(const vector<PathVariation1D> & variations, double dt, int steps){
    auto timeline = discrete_timeline(dt, steps);
    for (auto variation : variations){
      auto path = variation.DiscretePath(dt, steps);
      if (variation.IsFeasible(j_max_, a_max_, v_max_))
        plt::plot(timeline, path, "k-");
      else
        plt::plot(timeline, path, "grey");
    }
  }

  void Plot2DLanes(){
    vector<double> lane_line_s = {0,100};
    for (int i = 0; i <= map_.lanes_; i++){
      vector<double> lane_line_d;
      for (int j = 0; j < lane_line_s.size(); j++){
        lane_line_d.push_back(i * map_.lane_width_);
      }
      plt::plot(lane_line_s, lane_line_d, "grey");
    }  
  }

  void Plot2DVariations(const vector<PathVariation2D> & variations, int optimal_index, double dt, int steps){
    for (int i =0; i<variations.size(); i++){
      auto s_path = variations[i].s_variation_->DiscretePath(dt, steps);
      auto d_path = variations[i].d_variation_->DiscretePath(dt, steps);
      if (variations[i].is_collided_)
        plt::plot(s_path, d_path, "r-");
      else{
        if (i != optimal_index) plt::plot(s_path, d_path, "k-");
      }
    }
    if (optimal_index != -1){
      auto s_path = variations[optimal_index].s_variation_->DiscretePath(dt, steps);
      auto d_path = variations[optimal_index].d_variation_->DiscretePath(dt, steps);
      plt::plot(s_path, d_path, "g-");
    }
  }

  void Plot2DCars(vector<Car> cars){
    for (int i = 0; i < cars.size(); i++){
      auto s_path = discretize_path(cars[i].prediction_.first, dt_, planning_time_steps_);
      auto d_path = discretize_path(cars[i].prediction_.second, dt_, planning_time_steps_);
      vector<double> s_ends = {*(s_path.begin()),*(s_path.end()-1)};
      vector<double> d_ends = {*(d_path.begin()),*(d_path.end()-1)};
      if (cars[i].is_referenced_){
        plt::plot(s_path, d_path, "c--");
        plt::plot(s_ends, d_ends, "c>");
      }
      else {
        plt::plot(s_path, d_path, "b--");
        plt::plot(s_ends, d_ends, "b>");
      } 
    }
  }

  void OutputConfiguration(
    // my s state: s, ds, dds
      vector<double> s_state
    // my d state: d, dd, ddd
      , vector<double> d_state
    // neighbors states
      , vector<Car> cars){

    cout << "Configuarion ::\n"
    << s_state << ", " << d_state <<"\n"
    << ", {\n";
    for (int i = 0; i < cars.size(); i++){
      cout << "  ";
      if (i>0) cout<<", ";
      cout<<"Car("<< cars[i].s_state_ << ", " << cars[i].d_state_ <<")\n";
    }
    cout << "}\n";
  }

// Plots graphs of given road situation
  void PlotConfiguration(
    // my s state: s, ds, dds
    vector<double> s_state
    // my d state: d, dd, ddd
    , vector<double> d_state
    // neighbors states
    , vector<Car> cars){

    OutputConfiguration(s_state, d_state, cars);

    for (int i = 0; i < cars.size(); i++){
      cars[i].s_state_[0] -= s_state[0];
    }
    s_state[0] = 0;

    OutputConfiguration(s_state, d_state, cars);
    
    int my_lane = map_.LaneIndex(d_state[0]);

    double s_radius = v_max_ * planning_horizon_;
    instant t0 = now();

    vector<Car> watchable_cars = GetWatchableCars(s_state, d_state, cars, s_radius);

    instant t1 = now();

    cout << "s_radius:" << s_radius << "\n";
    cout << "safety_region:" << CalcSafetyRegion(s_state[1]) << "\n";
    cout << "safety_distance_:" << safety_distance_ << "\n";

    vector<PathVariation1D> s_variations;
    vector<PathVariation1D> d_variations;
 
    BuildVariations(s_state, d_state, watchable_cars, s_variations, d_variations);

    instant t2 = now();

    vector<PathVariation2D> vars_2D;
    CombineFeasible(s_variations, d_variations, vars_2D);
    
    instant t3 = now();

    int optimal_path_index = OptimalCollisionFree(vars_2D, s_state, d_state, watchable_cars, my_lane);

    instant t4 = now();

    CheckCollisions(vars_2D, watchable_cars, my_lane);

    instant t5 = now();

    auto get_cars = get_duration(t0,t1);
    cout << "Get cars time:" << get_cars << "\n";
    auto get_1d = get_duration(t1,t2);
    cout << "Gen 1D:" << get_1d << "\n";
    auto get_feas = get_duration(t2,t3);
    cout << "Feasibility:" << get_feas << "\n";
    auto get_opt = get_duration(t3,t4);
    cout << "Opt:" << get_opt << "\n";
    cout << "Total:" << get_cars+get_1d+get_feas+get_opt << "\n";
    auto get_collis = get_duration(t4,t5);
    cout << "ALL Collisions:" << get_collis << "\n";

    cout << "---" << "\n";
    cout << "Cars:" << watchable_cars.size() << "\n";
    cout << "S vars:" << s_variations.size() << "\n";
    cout << "D vars:" << d_variations.size() << "\n";
    cout << "VARS:" << vars_2D.size() << "\n";

    plt::subplot(3, 1, 1);
    plt::title("S variations. Black - feasible, grey - infeasible");
    //plt::axis("equal");
    Plot1DVariations(s_variations, dt_, planning_time_steps_);
    cout << "plS:\n";

    plt::subplot(3, 1, 2);
    plt::title("D variations. Black - feasible, grey - infeasible");
    //plt::axis("equal");
    Plot1DVariations(d_variations, dt_, planning_time_steps_);
    cout << "plD:\n";

    plt::subplot(3, 1, 3);
    plt::title("S+D. Black - possible, Red - collided, Green - optimal.");
    //plt::axis("equal");
    Plot2DLanes();
    Plot2DVariations(vars_2D, optimal_path_index,  dt_, planning_time_steps_);
    Plot2DCars(watchable_cars);

    plt::show();

  }

};



}
#endif /* PATH_PLANNER_H */