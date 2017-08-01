#ifndef TRAJECTORY_GENERATOR_H
#define TRAJECTORY_GENERATOR_H

#include <vector>
#include <utility>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"

using namespace std;

namespace path_planning
{

Eigen::VectorXd polyfit(const vector<double> &initial_state
                        ,const vector<double> &goal_state
                        ,double t){
    double a0 = initial_state[0];
    double a1 = initial_state[1];
    double a2 = initial_state[2]/2;

    double t2 = t*t;
    double t3 = t*t2;
    double t4 = t*t3;
    double t5 = t*t4;

    Eigen::MatrixXd A(3, 3);
    A <<    t3, t4, t5,
            3*t2, 4*t3, 5*t4,
            6*t, 12*t2, 20*t3;

    Eigen::VectorXd b(3);
    b <<    goal_state[0]-(a0+a1*t+a2*t2),
            goal_state[1]-(a1+2*a2*t),
            goal_state[2]-(2*a2);

    auto QR = A.householderQr();
    auto coeffs = QR.solve(b);

    Eigen::VectorXd result(6);
    result <<   a0,
                a1,
                a2,
                coeffs[0],
                coeffs[1],
                coeffs[2];
    return result;
}

double polyeval(const Eigen::VectorXd &coeffs, double x) {
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

vector<double> polyeval4(const Eigen::VectorXd &coeffs, double x) {
  double f = 0.0;
  double df = 0.0;
  double ddf = 0.0;
  double dddf = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    f += coeffs[i] * pow(x, i);
    if (i>0) df += i*coeffs[i] * pow(x, i-1);
    if (i>1) ddf += (i-1)*i*coeffs[i] * pow(x, i-2);
    if (i>2) dddf += (i-2)*(i-1)*i*coeffs[i] * pow(x, i-3);
  }
  return {f, df, ddf, dddf};
}

vector<vector<double>> generate_trajectory_1D(vector<double> &initial_state
                            , vector<double> &goal_state
                            , double duration, double dt){
    Eigen::VectorXd coeffs;
    coeffs = polyfit(initial_state, goal_state, duration);
    //cout << "Coeffs=";
    //for (int i = 0; i < coeffs.size(); i++) {
    //  cout<<" "<< coeffs[i]<<"*tˆ"<<i<<"+";
    //}
    //cout<<"\n----\n";
    vector<vector<double>> result;
    for (int i=0; i < int(duration/dt); i++){
        result.push_back(path_planning::polyeval4(coeffs, i*dt));
    }
    return result;
}

Eigen::VectorXd constant_speed_traj(double p, double s){
    Eigen::VectorXd result(6); 
        result << p,
        s,
        0,
        0,
        0,
        0;
    return result;
}

typedef pair<double, Eigen::VectorXd> traj_quintic;

vector<traj_quintic> sustained_acceleration_pulse(
    double p1, double s1, double p2, double s2
    , double a_max, double j_max, double dt_max, double D){

    double adt = a_max*dt_max;
    double adt2 = adt*dt_max;

    double sa = s1 + adt*0.5;
    double sb = s2 - adt*0.5;

    double dpi2 = 1/(pi()*pi());
    
    double D1 = adt2 * (0.25 - dpi2) + s1 * dt_max;
    double D2 = (sb * sb - sa * sa) * 0.5/a_max;
    double D3 = adt2 * (0.25 + dpi2) + sb * dt_max;

    cout<< "SUSTAINED: \ts1="<<s1<<"\ts2="<<s2<<"\tD="<<D<<"\tD1="<<D1<<"\tD2="<<D2<<"\tD3="<<D3<<"\n";
    if (D1+D2+D3-D > 0){
        cout<< "\tD1+D2+D3-D="<<D1+D2+D3-D<< "\n";
        sb = -adt + sqrt(adt*adt - 2*a_max*(adt2*0.5 + s1*dt_max - sa*sa*0.5/a_max - D));
        s2 = sb + adt * 0.5;
        D2 = (sb * sb - sa * sa) * 0.5/a_max;
        D3 = adt2 * (0.25 + dpi2) + sb * dt_max;
        cout<< "--\tsb="<<sb<< "\ts2="<<s2<<"\tD2="<<D2<<"\tD3="<<D3<<"\n";

        cout<< "\tD1+D2+D3-D="<<D1+D2+D3-D<< "\n";
    }

    vector<double> cp_1 = {p1, s1 ,0};
    vector<double> cp_a = {cp_1[0]+D1, sa, a_max};
    vector<double> cp_b = {cp_a[0]+D2, sb, a_max};
    vector<double> cp_2 = {cp_b[0]+D3, s2, 0};
    vector<double> cp_f = {p2, s2, 0};
    
    double t_1a = dt_max;
    double t_ab = (sb-sa)/a_max;
    double t_b2 = dt_max;
    double t_2f = (D-D1-D2-D3)/s2;
    
    auto coeffs_1a = polyfit(cp_1, cp_a, t_1a);
    auto coeffs_ab = polyfit(cp_a, cp_b, t_ab);
    auto coeffs_b2 = polyfit(cp_b, cp_2, t_b2);
    Eigen::VectorXd coeffs_2f;
    if (t_2f > 0.1)
        coeffs_2f = polyfit(cp_2, cp_f, t_2f);
    else 
        coeffs_2f = constant_speed_traj(cp_2[0],cp_2[1]);

    return {
        {t_1a, coeffs_1a}, 
        {t_1a+t_ab, coeffs_ab},
        {t_1a+t_ab+t_b2, coeffs_b2},
        {t_1a+t_ab+t_b2+t_2f, coeffs_2f}};
}

vector<traj_quintic> acceleration_pulse(
    double p1, double p2, double s1, double s2
    , double a_max, double j_max, double dt_max
    , double D, double D_min, double ds_min){

    double ds = s2-s1;
    
    double dt_ramp = dt_max;
    double a_peak = a_max;
    if (D < D_min){
        double pi_j = pi()/j_max;
        double s1pi_j = pi_j*s1;
        double Dpi_2j = pi_j*D*0.5;
        double f;
        do {
            double dt_ramp2 = dt_ramp*dt_ramp; 
            f = (dt_ramp2*dt_ramp + s1pi_j*dt_ramp - Dpi_2j);
            double delta = f/(2*dt_ramp2+s1pi_j);
            cout<< "---\tdt_ramp="<<dt_ramp<< "\tdelta="<<delta<< "\tf="<<f<<"\n";
            dt_ramp -= delta;
        }    
        while (fabs(f)>0.0001);
        
        a_peak = 2*j_max*dt_ramp/pi();
    }
    cout<< "\tD="<<D<< "\tD_min="<<D_min
    << "\tds="<<ds<< "\tds_min="<<ds_min
    << "\ta_peak="<<a_peak<< "\tdt_ramp="<<dt_ramp
    << "\ta_peak*dt_ramp="<<a_peak*dt_ramp<<"\n";
    if (ds < a_peak*dt_ramp){
        cout<<"ds < a_peak*dt_ramp\n";
        dt_ramp = sqrt(pi()*(s2-s1)*0.5/j_max);
        a_peak = 2*j_max*dt_ramp/pi();
    }

    double dpi2 = 1/(pi()*pi());

    double adt = a_peak*dt_ramp;
    double adt2 = adt*dt_ramp;

    double sa = s1 + adt*0.5;
    double sb = s1 + adt;

    double D1 = adt2 * (0.25 - dpi2) + s1 * dt_ramp;
    double D2 = adt2 * (0.25 + dpi2) + sa * dt_ramp;
     cout<< "\tD1="<<D1<< "\tD2="<<D2
    << "\tsa="<<sa<< "\tsb="<<sb
    << "\ts1="<<s1<< "\ts2="<<s2
    <<"\n";

    vector<double> cp_1 = {p1, s1 ,0};
    vector<double> cp_a = {cp_1[0]+D1, sa, a_peak};
    vector<double> cp_b = {cp_a[0]+D2, sb, 0};
    vector<double> cp_2 = {p2, s2, 0};
    
    double t_1a = dt_ramp;
    double t_ab = dt_ramp;
    double t_b2 = (D-D1-D2)/sb;
    
    auto coeffs_1a = polyfit(cp_1, cp_a, t_1a);
    auto coeffs_ab = polyfit(cp_a, cp_b, t_ab);
    Eigen::VectorXd coeffs_b2;
    cout<<"yep! t_b2="<<t_b2<<"\n";
    if (t_b2 > 0.1)
        coeffs_b2 = polyfit(cp_b, cp_2, t_b2);
    else 
        coeffs_b2 = constant_speed_traj(cp_b[0],cp_b[1]);

    return {
        {t_1a, coeffs_1a}, 
        {t_1a+t_ab, coeffs_ab},
        {t_1a+t_ab+t_b2, coeffs_b2}};
}

vector<vector<double>> generate_trajectory_1Dq(
    vector<double> &initial_state
    , vector<double> &goal_state
    , vector<double> &kinematic_limits
    ,double dt, int horizon){

    double a_max = kinematic_limits[0];
    double j_max = kinematic_limits[1];

    double dt_max = pi()*a_max*0.5/j_max;

    double p1 = initial_state[0];
    double s1 = initial_state[1];
    double a1 = initial_state[2];

    double p2 = goal_state[0];
    double s2 = goal_state[1];
    double a2 = goal_state[2];
    
    double D = p2-p1;
    double D_min = a_max*dt_max*dt_max + 2*s1*dt_max;
    double ds = s2-s1;
    double ds_min = a_max*dt_max;

    vector<traj_quintic> traj_quintics;
    if ((D > D_min) && (ds > ds_min))
        traj_quintics = sustained_acceleration_pulse(p1, s1, p2, s2, a_max, j_max, dt_max, D);
    else 
        traj_quintics = acceleration_pulse(p1, p2, s1, s2, a_max, j_max, dt_max, D, D_min, ds_min);

    //cout<< "Timepoints:"<<"\t"<<t_1a
    //<<"\t"<<t_1a+t_ab<<"\t"<<t_1a+t_ab+t_b2
    //<<"\t"<<t_1a+t_ab+t_b2+t_2f<< "\n";
 
    auto traj_quintic = traj_quintics.begin();
    auto final_quintic = traj_quintics.end();
    --final_quintic;
    double start_time = 0;
    vector<vector<double>> result;
    for (int i=0; i < horizon; i++){
        double t = i * dt;
        if ((t > traj_quintic->first) && (traj_quintic != final_quintic)){
            start_time = traj_quintic->first;
            ++traj_quintic;           
        }

        auto coeffs = traj_quintic->second;
        result.push_back(path_planning::polyeval4(coeffs, t - start_time));
    }

    return result;
}


}


#endif /* TRAJECTORY_GENERATOR_H */