#ifndef TRAJECTORY_GENERATOR_H
#define TRAJECTORY_GENERATOR_H

#include <algorithm>
#include <cassert>
#include <vector>
#include <utility>
#include "poly.h"

using namespace std;

namespace path_planning
{
typedef vector<double> discrete_path_1D;

typedef pair<discrete_path_1D, discrete_path_1D> discrete_path_2D;

typedef pair<polynomial, double> trajectory_segment;

typedef vector<trajectory_segment> continuous_path_1D;

typedef pair<continuous_path_1D, continuous_path_1D> continuous_path_2D;

constexpr double TIME_INF = 1000;

struct PathVariation1D{
private:
  // cached values
  double integral_jerk_ = -1;
  double max_jerk_ = -1;
  double max_a_ = -1;
  double max_v_ = -1;

  double discretization_step_ = -1;
  double discretization_steps_ = -1;
  vector<vector<double>> discrete_trajectory_;

public:

  int lane_index_ = -1;

  continuous_path_1D path_;

  vector<double> state_var_;

  PathVariation1D(polynomial poly, double T) {
    state_var_ = polyeval_v(poly, T);
    path_ = {{poly, T}, {keep_state_poly(state_var_), TIME_INF}};
  }

  polynomial GetPoly(){
    return path_[0].first;
  }

  polynomial GetContinuation(){
    return path_[1].first;
  }

  double GetT(){
    return path_[0].second;
  }

  double GetIntegralJerk(){
    if (integral_jerk_ < 0){
      integral_jerk_ = poly_integral_jerk(GetPoly(), 0, GetT());
    }
    return integral_jerk_;
  }

  double GetMaxJerk(){
    if (max_jerk_ < 0){
      max_jerk_ = poly_max_jerk(GetPoly(), 0, GetT());
    }
    return max_jerk_;
  }

  double GetMaxAcceleration(){
    if (max_a_ < 0){
      max_a_ = poly_max_acceleration(GetPoly(), 0, GetT());
    }
    return max_a_;
  }

  double GetMaxVelocity(double search_lim = -1){
    if (max_v_ < 0){
      int steps = 10;
      double dt = GetT()/steps;
      max_v_ = 0;
      for (int i = 1; i <= steps; i++){
        double velocity = fabs(EvalAllAt(i*dt)[1]);
        if (velocity > max_v_) max_v_ = velocity;
        if (max_v_ >= search_lim) break;
      }
    }
    return max_v_;
  }

  vector<double> EvalAllAt(double t){
    auto poly = GetPoly();
    double T = GetT();
    double start_time = 0;
    if (t > T){
      start_time = T;
      poly = GetContinuation();
    }
    return polyeval_v(poly, t - start_time);
  }

  double EvalPathAt(double t){
    auto poly = GetPoly();
    double T = GetT();
    double start_time = 0;
    if (t > T){
      start_time = T;
      poly = GetContinuation();
    }
    return polyeval(poly, t - start_time);
  }

  vector<vector<double>>* DiscreteTrajectory(double dt, int steps){
    if (fabs(discretization_step_-dt) > 1e-10 
      || discretization_steps_ != steps){

      discretization_step_ = dt;
      discretization_steps_ = steps;

      discrete_trajectory_.clear();
      for (int i=0; i <= steps; i++){
        discrete_trajectory_.push_back(EvalAllAt(i * dt));
      }
    }
    return &discrete_trajectory_;
  }

  vector<double> DiscretePath(double dt, int steps){
    auto trajectory = DiscreteTrajectory(dt, steps);
    vector<double> path;
    for (int i=0; i < trajectory->size(); i++){
      path.push_back((*trajectory)[i][0]);
    }
    return path;
  }

  bool IsFeasible(double max_jerk, double max_a, double max_v){
    return 
      GetMaxJerk() < max_jerk 
      && GetMaxAcceleration() < max_a
      && GetMaxVelocity(max_v) <= max_v;
  }

};

struct PathVariation2D{
  PathVariation1D* s_variation_;
  PathVariation1D* d_variation_;
  bool is_collided_ = false;
  PathVariation2D(PathVariation1D* s_variation, PathVariation1D* d_variation){
    s_variation_ = s_variation;
    d_variation_ = d_variation;
  }

double CostFunctional(
  vector<double> s_state
  , vector<double> d_state
  , double time_horizon, double v_max){

    vector<double> final_s_state = s_variation_->EvalAllAt(time_horizon);
    double final_d = d_variation_->EvalPathAt(time_horizon);

    double s_diff = final_s_state[0] - s_state[0];
    double d_diff = fabs(final_d - d_state[0]);

    double v_cost = 1 - fabs(final_s_state[1])/v_max;
    double s_cost = 1 - s_diff/(v_max*time_horizon);
    double d_cost = d_diff;
    return v_cost + s_cost + 0.007 * d_cost + 0.01*d_variation_->GetT();
  }
};

struct PathVariation2DComparator {
  vector<double> s_state_;
  vector<double> d_state_;
  double time_horizon_;
  double v_max_;

  PathVariation2DComparator(
    vector<double> s_state
    , vector<double> d_state
    , double time_horizon
    , double v_max){
    time_horizon_ = time_horizon;
    v_max_ = v_max;
    s_state_ = s_state;
    d_state_ = d_state;
  }

  bool operator() (PathVariation2D& path1, PathVariation2D& path2){ 
    double cost1 = path1.CostFunctional(s_state_, d_state_, time_horizon_, v_max_);
    double cost2 = path2.CostFunctional(s_state_, d_state_, time_horizon_, v_max_);
    return cost1 < cost2;
  }
};


struct ReferenceState{
  enum reference_type {STATIC, DYNAMIC};
  reference_type type_; 
  vector<double> state_;
  
  vector<double> parameter_variations_;

  double safety_distance_ = 0;
  double safety_time_ = 0;
  
  ReferenceState(){
    type_ = STATIC;
  }

  ReferenceState(vector<double> state, reference_type type = STATIC){
    state_ = state;
    type_ = type;
  }

  static void NegativeVariations(
    ReferenceState& ref_state
    , double max_negative_deviation
    , int steps){
    ref_state.parameter_variations_.clear();
    double dparam = -max_negative_deviation/steps;
    for (int i=0; i<steps; i++){
      ref_state.parameter_variations_.push_back(i*dparam);
    }
  }

  static ReferenceState KeepVelocityState(
    double velocity
    , double max_dev_ratio = 0.3
    , int var_steps = 3){
    ReferenceState ref_state;
    ref_state.state_ = {velocity, 0};
    NegativeVariations(ref_state, velocity*max_dev_ratio, var_steps);
    return ref_state;
  }

  static ReferenceState FollowState(
    vector<double> state
    , double safety_distance
    , double safety_time
    , double max_dev_ratio = 0.3
    , int var_steps = 3){
    ReferenceState ref_state(state, DYNAMIC);
    ref_state.safety_distance_ = safety_distance;
    ref_state.safety_time_ = safety_time;
    NegativeVariations(ref_state, state[0]*max_dev_ratio, var_steps);
    return ref_state;
  }
  
  static void LaneVariations(
    ReferenceState& ref_state
    , double lane_width
    , int lanes){
    
    ref_state.parameter_variations_.clear();
    for (int i=0; i<lanes; i++){
      int var_index = i; 
      ref_state.parameter_variations_.push_back(var_index * lane_width);
    }
  }

  static ReferenceState LanesState(
    double lane_width
    , int lanes){
    ReferenceState ref_state;

    ref_state.state_ = {get_lane_center(0, lane_width), 0, 0};
    LaneVariations(ref_state, lane_width, lanes);
    return ref_state;
  }


  static ReferenceState KeepLaneState(
    int lane_index
    , double lane_width
    , int lanes){
    ReferenceState ref_state;

    assert(lane_index!=-1);
    double lane_center = get_lane_center(lane_index, lane_width);

    ref_state.state_ = {lane_center, 0, 0};
    ref_state.parameter_variations_ = {0};
    return ref_state;
  }

  static ReferenceState KeepLaneState(
    double frenet_d
    , double lane_width
    , int lanes){
    int lane_index = get_lane_index(frenet_d, lane_width, lanes);
    return ReferenceState::KeepLaneState(lane_index, lane_width, lanes);
  }

  virtual vector<double> Variate(double parameter_variation, double time_variation){
    vector<double> state_variation(state_);
    if (type_ == STATIC){
      state_variation[0] += parameter_variation;
    }
    else{
      double safety_region = safety_distance_ + safety_time_ * state_[1];
      state_variation[0] += parameter_variation 
                  + (state_[1] + state_[2] * time_variation * 0.5) * time_variation
                  - safety_region
                  ;
      state_variation[1] += state_[2]*time_variation;
    }
    return state_variation;
  }

};

// Main variation function. Builds variations of trajectories by variating time and referenced parameter
vector<PathVariation1D> variate_1D(
    vector<double> state0
    , vector<ReferenceState> ref_states
    , double dt = 1, int time_steps = 5
    ){
    vector<PathVariation1D> result;
    for (auto ref_state : ref_states){
      for (double param_variation : ref_state.parameter_variations_){
        for (int j = 1; j <= time_steps; j++){
          double T = j * dt; 
          auto state_variation = ref_state.Variate(param_variation, T);
          auto poly = solve_poly(state0, state_variation, T);
          result.push_back(PathVariation1D(poly, T));
        }
      }
    }
    return result;
}

// builds path of constant velocity and acceleration
continuous_path_1D keep_state_trajectory_1D(
  vector<double> state0){
  continuous_path_1D path; 
  auto new_coeffs = keep_state_poly(state0);
  path.push_back({new_coeffs, 1000});
  return path;
}

// discretizes path, with dt on steps, using time shift
discrete_path_1D discretize_path(continuous_path_1D path
  , double dt, int steps, double time_shift = 0){

  auto traj_seg = path.begin();
  auto final_seg = path.end();
  --final_seg;
  discrete_path_1D result;

  double start_time = 0;
  while ((time_shift > traj_seg->second) && (traj_seg != final_seg)){
    start_time = traj_seg->second;
    ++traj_seg;           
  }
  
  for (int i=0; i < steps; i++){
      double t = time_shift + i * dt;
      if ((t > traj_seg->second) && (traj_seg != final_seg)){
          start_time = traj_seg->second;
          ++traj_seg;           
      }

      auto coeffs = traj_seg->first;
      result.push_back(path_planning::polyeval(coeffs, t - start_time));
  }
  return result;
}

// takes max_steps from source
discrete_path_2D derive_discrete_path(discrete_path_2D source, int max_steps){
  discrete_path_2D result;
  int steps = source.first.size();
  if (steps > max_steps) steps = max_steps;
  for (int i =0; i<steps; i++){
    result.first.push_back(source.first[i]);
    result.second.push_back(source.second[i]);
  }
  return result;
}


}


#endif /* TRAJECTORY_GENERATOR_H */