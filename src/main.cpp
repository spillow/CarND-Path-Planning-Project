#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"
#include <random>

#include "Eigen-3.3/Eigen/Dense"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
    auto found_null = s.find("null");
    auto b1 = s.find_first_of("[");
    auto b2 = s.find_first_of("}");
    if (found_null != string::npos) {
        return "";
    }
    else if (b1 != string::npos && b2 != string::npos) {
        return s.substr(b1, b2 - b1 + 2);
    }
    return "";
}

double distance(double x1, double y1, double x2, double y2)
{
    return sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

    double closestLen = 100000; //large number
    int closestWaypoint = 0;

    for (int i = 0; i < maps_x.size(); i++)
    {
        double map_x = maps_x[i];
        double map_y = maps_y[i];
        double dist = distance(x, y, map_x, map_y);
        if (dist < closestLen)
        {
            closestLen = dist;
            closestWaypoint = i;
        }

    }

    return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

    int closestWaypoint = ClosestWaypoint(x, y, maps_x, maps_y);

    double map_x = maps_x[closestWaypoint];
    double map_y = maps_y[closestWaypoint];

    double heading = atan2((map_y - y), (map_x - x));

    double angle = abs(theta - heading);

    if (angle > pi() / 4)
    {
        closestWaypoint++;
    }

    return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(
    double x, double y, double theta,
    const vector<double> &maps_x, const vector<double> &maps_y)
{
    int next_wp = NextWaypoint(x, y, theta, maps_x, maps_y);

    int prev_wp;
    prev_wp = next_wp - 1;
    if (next_wp == 0)
    {
        prev_wp = maps_x.size() - 1;
    }

    double n_x = maps_x[next_wp] - maps_x[prev_wp];
    double n_y = maps_y[next_wp] - maps_y[prev_wp];
    double x_x = x - maps_x[prev_wp];
    double x_y = y - maps_y[prev_wp];

    // find the projection of x onto n
    double proj_norm = (x_x*n_x + x_y*n_y) / (n_x*n_x + n_y*n_y);
    double proj_x = proj_norm*n_x;
    double proj_y = proj_norm*n_y;

    double frenet_d = distance(x_x, x_y, proj_x, proj_y);

    //see if d value is positive or negative by comparing it to a center point

    double center_x = 1000 - maps_x[prev_wp];
    double center_y = 2000 - maps_y[prev_wp];
    double centerToPos = distance(center_x, center_y, x_x, x_y);
    double centerToRef = distance(center_x, center_y, proj_x, proj_y);

    if (centerToPos <= centerToRef)
    {
        frenet_d *= -1;
    }

    // calculate s value
    double frenet_s = 0;
    for (int i = 0; i < prev_wp; i++)
    {
        frenet_s += distance(maps_x[i], maps_y[i], maps_x[i + 1], maps_y[i + 1]);
    }

    frenet_s += distance(0, 0, proj_x, proj_y);

    return{ frenet_s,frenet_d };

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(
    double s, double d,
    const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y,
    const vector<double> &maps_dx, const vector<double> &maps_dy)
{
#if 1
    int prev_wp = -1;

    while (s > maps_s[prev_wp + 1] && (prev_wp < (int)(maps_s.size() - 1)))
    {
        prev_wp++;
    }

    int wp2 = (prev_wp + 1) % maps_x.size();
    int wp3 = (prev_wp + 2) % maps_x.size();

    auto fillpoints = [&](tk::spline &s, const vector<double> &mx, const vector<double> &my)
    {
        vector<double> xs = { mx[prev_wp], mx[wp2], mx[wp3] };
        vector<double> ys = { my[prev_wp], my[wp2], my[wp3] };

        s.set_points(xs, ys);
    };

    tk::spline spline_x_s;
    fillpoints(spline_x_s, maps_s, maps_x);

    tk::spline spline_y_s;
    fillpoints(spline_y_s, maps_s, maps_y);

    tk::spline spline_dx_s;
    fillpoints(spline_dx_s, maps_s, maps_dx);

    tk::spline spline_dy_s;
    fillpoints(spline_dy_s, maps_s, maps_dy);

    double x = spline_x_s(s) + d * spline_dx_s(s);
    double y = spline_y_s(s) + d * spline_dy_s(s);

    return { x,y };
#else
    int prev_wp = -1;

    while (s > maps_s[prev_wp + 1] && (prev_wp < (int)(maps_s.size() - 1)))
    {
        prev_wp++;
    }

    int wp2 = (prev_wp + 1) % maps_x.size();

    double heading = atan2((maps_y[wp2] - maps_y[prev_wp]), (maps_x[wp2] - maps_x[prev_wp]));
    // the x,y,s along the segment
    double seg_s = (s - maps_s[prev_wp]);

    double seg_x = maps_x[prev_wp] + seg_s*cos(heading);
    double seg_y = maps_y[prev_wp] + seg_s*sin(heading);

    double perp_heading = heading - pi() / 2;

    double x = seg_x + d*cos(perp_heading);
    double y = seg_y + d*sin(perp_heading);

    return{ x,y };
#endif
}

const double TIME_STEP = 0.02; // 20 ms
const unsigned BUFFER_SIZE = 50;

class Ego
{
public:
    double x;
    double y;
    double s;
    double d;
    double yaw;
    double speed;

    double ref_vel;

    unsigned get_lane() const
    {
        if (d < 4)
            return 0;
        else if (d < 8)
            return 1;
        else
            return 2;
    }
};

class Vehicle
{
public:
    double id;
    double x;
    double y;
    double vx;
    double vy;
    double s;
    double d;

    double speed() const
    {
        return sqrt(vx*vx + vy*vy);
    }

    double _s(double t = 0.0) const
    {
        return s + speed() * t;
    }

    double _sd(double t = 0.0) const
    {
        return speed();
    }

    double _sdd(double t = 0.0) const
    {
        return 0.0;
    }

    double _d(double t = 0.0) const
    {
        return d;
    }

    double _dd(double t = 0.0) const
    {
        return 0.0;
    }

    double _ddd(double t = 0.0) const
    {
        return 0.0;
    }

    vector<double> state_predict(double t = 0.0) const
    {
        return {
            _s(t) + _sd(t) * t + 0.5 * _sdd(t) * t * t,
            _sd(t) + _sdd(t) * t,
            _sdd(t),
            _d(t) + _dd(t) * t + 0.5 * _ddd(t) * t * t,
            _dd(t) + _ddd(t) * t,
            _ddd(t)
        };
    }

    unsigned get_lane() const
    {
        if (d < 4)
            return 0;
        else if (d < 8)
            return 1;
        else
            return 2;
    }
};

enum BehaviorState
{
    KEEP_LANE,
    PREPARE_CHANGE_LANE,
    CHANGE_LANE
};

class SensorData
{
public:
    Ego ego;
    std::vector<Vehicle> vehicles;
};

class StateInfo
{
public:
    std::vector<double> previous_path_x;
    std::vector<double> previous_path_y;
    double end_path_s;
    double end_path_d;
    std::function<vector<double>(double, double)> getXY;

    double end_path_time() const
    {
        return (double)previous_path_x.size() * TIME_STEP;
    }
};

class Context
{
public:
    StateInfo Info;
    SensorData Data;

    double end_path_time() const
    {
        return Info.end_path_time();
    }
};

SensorData getSensorData(
    const std::vector<std::vector<double>> &sensor_fusion,
    double car_x, double car_y, double car_s, double car_d,
    double car_yaw, double car_speed, double ref_vel)
{
    std::vector<Vehicle> vehicles;
    for (auto &V : sensor_fusion)
    {
        vehicles.push_back({ V[0], V[1], V[2], V[3], V[4], V[5], V[6] });
    }

    Ego ego{ car_x, car_y, car_s, car_d, car_yaw, car_speed, ref_vel };

    return{ ego, vehicles };
}

StateInfo getStateInfo(
    const std::vector<double> &previous_path_x,
    const std::vector<double> &previous_path_y,
    double end_path_s,
    double end_path_d,
    std::function<vector<double>(double, double)> getXY)
{
    return{
        previous_path_x, previous_path_y,
        end_path_s, end_path_d, getXY
    };
}

Context getContext(const SensorData &Data, const StateInfo &Info)
{
    return{ Info, Data };
}

bool open_lane(const Context &Ctx, unsigned lane, vector<Vehicle> &block_vehicles)
{
    auto &Vehicles = Ctx.Data.vehicles;
    double t = Ctx.end_path_time();

    bool ret = true;

    for (auto &V : Vehicles)
    {
        if (V.get_lane() != lane)
            continue;

        double dist_diff = abs(Ctx.Info.end_path_s - V._s(t));

        if (dist_diff < 7.0)
        {
            block_vehicles.push_back(V);
            ret = false;
        }

        if (dist_diff < 30.0)
        {
            if (V._s(t) > Ctx.Info.end_path_s)
            {
                if (V.speed() > Ctx.Data.ego.ref_vel)
                    continue;
            }
            else
            {
                if (Ctx.Data.ego.ref_vel > V.speed())
                    continue;
            }
            block_vehicles.push_back(V);
            ret = false;
        }
    }

    return ret;
}

bool look_for_lane_change(const Context &Ctx)
{
    auto &Vehicles = Ctx.Data.vehicles;
    double t = Ctx.end_path_time();
    double car_pred_s = Ctx.Info.end_path_s;

    for (auto &V : Vehicles)
    {
        float d = V.d;
        if (V.get_lane() == Ctx.Data.ego.get_lane())
        {
            double vehicle_pred_s = V._s(t);

            if (vehicle_pred_s > car_pred_s && (vehicle_pred_s - car_pred_s) < 30.0)
            {
                return true;
            }
        }
    }

    return false;
}

std::vector<Vehicle> cars_ahead(const Context &Ctx, double s_begin, double s_end, unsigned lane)
{
    auto &Vehicles = Ctx.Data.vehicles;
    double t = Ctx.end_path_time();
    double car_pred_s = Ctx.Info.end_path_s;

    std::vector<Vehicle> in_range;

    for (auto &V : Vehicles)
    {
        if (V.get_lane() == lane)
        {
            double vehicle_pred_s = V._s(t);

            if (vehicle_pred_s > s_begin && vehicle_pred_s < s_end)
                in_range.push_back(V);
        }
    }

    return in_range;
}

bool check_lane_opening(const Context &Ctx, unsigned &lane, vector<Vehicle> &blocking_vehicles)
{
    bool ret = false;

    if (lane == 0)
    {
        if (ret = open_lane(Ctx, 1, blocking_vehicles))
            lane = 1;
    }
    else if (lane == 1)
    {
        bool left = open_lane(Ctx, 0, blocking_vehicles);
        bool right = open_lane(Ctx, 2, blocking_vehicles);

        ret = left || right;

        if (left && right)
        {
            // pick the one with less traffic ahead
            double s_begin = Ctx.Info.end_path_s;
            double s_end = s_begin + 70;

            if (cars_ahead(Ctx, s_begin, s_end, 0).size() <=
                cars_ahead(Ctx, s_begin, s_end, 2).size())
            {
                lane = 0;
            }
            else
            {
                lane = 2;
            }
        }
        else if (left)
            lane = 0;
        else if (right)
            lane = 2;
    }
    else if (lane == 2)
    {
        if (ret = open_lane(Ctx, 1, blocking_vehicles))
            lane = 1;
    }
    else
    {
        assert(0 && "not in a lane?");
    }

    return ret;
}

const Vehicle *get_vehicle_in_front(const Context &Ctx, unsigned lane, double &dist_s)
{
    auto &Vehicles = Ctx.Data.vehicles;
    double t = Ctx.end_path_time();
    double car_pred_s = Ctx.Info.end_path_s;

    const Vehicle *Closest = nullptr;
    double min_s = +INFINITY;

    for (auto &V : Vehicles)
    {
        if (V.get_lane() == lane)
        {
            double vehicle_pred_s = V._s(t);

            if (vehicle_pred_s > car_pred_s)
            {
                if (vehicle_pred_s < min_s)
                {
                    dist_s = vehicle_pred_s - car_pred_s;
                    min_s = vehicle_pred_s;
                    Closest = &V;
                }
            }
        }
    }

    return Closest;
}

void fill_straight_path(
    const Context &Ctx, double ref_vel, unsigned lane, vector<double> &new_x, vector<double> &new_y)
{
    std::vector<double> ptsx;
    std::vector<double> ptsy;

    auto &Ego = Ctx.Data.ego;
    double car_x = Ego.x;
    double car_y = Ego.y;
    double car_yaw = Ego.yaw;
    double car_s = Ego.s;
    double end_path_s = Ctx.Info.end_path_s;
    auto getXY = Ctx.Info.getXY;

    auto &previous_path_x = Ctx.Info.previous_path_x;
    auto &previous_path_y = Ctx.Info.previous_path_y;
    unsigned prev_size = previous_path_x.size();

    double ref_x = car_x;
    double ref_y = car_y;
    double ref_yaw = deg2rad(car_yaw);

    double car_pred_s = (prev_size > 0) ? end_path_s : car_s;

    if (prev_size < 2)
    {
        double prev_car_x = car_x - cos(car_yaw);
        double prev_car_y = car_y - sin(car_yaw);

        ptsx.push_back(prev_car_x);
        ptsx.push_back(car_x);

        ptsy.push_back(prev_car_y);
        ptsy.push_back(car_y);
    }
    else
    {
        ref_x = previous_path_x[prev_size - 1];
        ref_y = previous_path_y[prev_size - 1];

        double ref_x_prev = previous_path_x[prev_size - 2];
        double ref_y_prev = previous_path_y[prev_size - 2];

        ref_yaw = atan2(ref_y - ref_y_prev, ref_x - ref_x_prev);

        ptsx.push_back(ref_x_prev);
        ptsx.push_back(ref_x);

        ptsy.push_back(ref_y_prev);
        ptsy.push_back(ref_y);
    }

    std::vector<double> next_wp0 = getXY(car_pred_s + 30, (2 + 4 * lane));
    std::vector<double> next_wp1 = getXY(car_pred_s + 60, (2 + 4 * lane));
    std::vector<double> next_wp2 = getXY(car_pred_s + 90, (2 + 4 * lane));

    ptsx.push_back(next_wp0[0]);
    ptsx.push_back(next_wp1[0]);
    ptsx.push_back(next_wp2[0]);

    ptsy.push_back(next_wp0[1]);
    ptsy.push_back(next_wp1[1]);
    ptsy.push_back(next_wp2[1]);

    for (int i = 0; i < ptsx.size(); i++)
    {
        // shift reference angle
        double shift_x = ptsx[i] - ref_x;
        double shift_y = ptsy[i] - ref_y;

        ptsx[i] = (shift_x * cos(-ref_yaw) - shift_y * sin(-ref_yaw));
        ptsy[i] = (shift_x * sin(-ref_yaw) + shift_y * cos(-ref_yaw));
    }

    tk::spline s;

    s.set_points(ptsx, ptsy);

    double target_x = 30.0;
    double target_y = s(target_x);
    double target_dist = sqrt(target_x*target_x + target_y*target_y);

    double x_add_on = 0;

    for (int i = 1; i < BUFFER_SIZE - prev_size; i++)
    {
        // N * TIME_STEP * v = d
        double N = (target_dist / (TIME_STEP * ref_vel / 2.24));
        double x_point = x_add_on + target_x / N;
        double y_point = s(x_point);

        x_add_on = x_point;

        double x_ref = x_point;
        double y_ref = y_point;

        x_point = x_ref * cos(ref_yaw) - y_ref * sin(ref_yaw);
        y_point = x_ref * sin(ref_yaw) + y_ref * cos(ref_yaw);

        x_point += ref_x;
        y_point += ref_y;

        new_x.push_back(x_point);
        new_y.push_back(y_point);
    }
}

vector<double> JMT(vector<double> start, vector <double> end, double T)
{
    /*
    Calculate the Jerk Minimizing Trajectory that connects the initial state
    to the final state in time T.

    INPUTS

    start - the vehicles start location given as a length three array
    corresponding to initial values of [s, s_dot, s_double_dot]

    end   - the desired end state for vehicle. Like "start" this is a
    length three array.

    T - The duration, in seconds, over which this maneuver should occur.

    OUTPUT
    an array of length 6, each value corresponding to a coefficent in the polynomial
    s(t) = a_0 + a_1 * t + a_2 * t**2 + a_3 * t**3 + a_4 * t**4 + a_5 * t**5

    EXAMPLE

    > JMT( [0, 10, 0], [10, 10, 0], 1)
    [0.0, 10.0, 0.0, 0.0, 0.0, 0.0]
    */

    double si = start[0];
    double sid = start[1];
    double sidd = start[2];

    double sf = end[0];
    double sfd = end[1];
    double sfdd = end[2];

    double a0 = si;
    double a1 = sid;
    double a2 = 0.5*sidd;

    MatrixXd lhs = MatrixXd(3, 3);
    VectorXd rhs = VectorXd(3);

    double T2 = T*T;
    double T3 = T2*T;
    double T4 = T3*T;
    double T5 = T4*T;

    lhs << T3, T4, T5,
        3 * T2, 4 * T3, 5 * T4,
        6 * T, 12 * T2, 20 * T3;

    rhs << sf - (si + sid*T + 0.5*sidd*T2),
        sfd - (sid + sidd*T),
        sfdd - sidd;

    VectorXd res = lhs.inverse() * rhs;

    return{ a0, a1, a2, res[0], res[1], res[2] };
}

vector<double> add_vectors(const vector<double> &a, const vector<double> &b)
{
    assert(a.size() == b.size());
    vector<double> newvec;
    for (unsigned i = 0; i < a.size(); i++)
    {
        newvec.push_back(a[i] + b[i]);
    }

    return newvec;
}

std::pair<vector<double>, vector<double>> perturb_goal(vector<double> goal_s, vector<double> goal_d)
{
    vector<double> SIGMA_S{ 10.0, 4.0, 2.0 };
    vector<double> SIGMA_D{ 1.0, 1.0, 1.0 };

    default_random_engine gen;

    vector<double> new_goal_s;
    for (unsigned i = 0; i < 3; i++)
    {
        std::normal_distribution<double> distr_s(goal_s[i], SIGMA_S[i]);
        new_goal_s.push_back(distr_s(gen));
    }

    vector<double> new_goal_d;
    for (unsigned i = 0; i < 3; i++)
    {
        std::normal_distribution<double> distr_d(goal_d[i], SIGMA_D[i]);
        new_goal_d.push_back(distr_d(gen));
    }

    return std::make_pair(new_goal_s, new_goal_d);
}

typedef std::tuple<vector<double>, vector<double>, double> Trajectory;

void print_vec(const vector<double> &vec)
{
    std::cout << "[ ";
    for (unsigned i = 0; i < vec.size(); i++)
    {
        std::cout << vec[i] << ",";
    }

    std::cout << " ]" << std::endl;
}


Trajectory findTrajectory(
    vector<double> start_s, vector<double> start_d, vector<double> goal_s, vector<double> goal_d,
    double T, const Context &Ctx)
{
    // TODO: experiment with time range to search
    const double timestep = 0.5;
    const unsigned N_SAMPLES = 10;

    double t = T - 4.0 * timestep;
    vector<std::tuple<vector<double>, vector<double>, double>> goals;
    while (t <= T + 4.0 * timestep)
    {
        for (unsigned i = 0; i < N_SAMPLES; i++)
        {
            //auto perturbed = perturb_goal(goal_s, goal_d);
            //goals.push_back(std::make_tuple(perturbed.first, perturbed.second, t));
            goals.push_back(std::make_tuple(goal_s, goal_d, t));
        }

        t += timestep;
    }

    vector<Trajectory> trajectories;
    for (auto &goal : goals)
    {
        vector<double> s_goal;
        vector<double> d_goal;
        double t;
        std::tie(s_goal, d_goal, t) = goal;

        auto s_coeffs = JMT(start_s, s_goal, t);
        auto d_coeffs = JMT(start_d, d_goal, t);

        trajectories.push_back(Trajectory(s_coeffs, d_coeffs, t));
    }

    return trajectories[0];
}

// Evaluate a polynomial.
double polyeval(const vector<double> &coeffs, double x) {
    double result = 0.0;
    for (int i = 0; i < coeffs.size(); i++) {
        result += coeffs[i] * pow(x, i);
    }
    return result;
}

void convertTrajectoryToXY(
    Trajectory Traj, const StateInfo &Info, vector<double> &new_x, vector<double> &new_y)
{
    const unsigned NUM_SAMPLES = 5;

    auto &s_coeffs = std::get<0>(Traj);
    auto &d_coeffs = std::get<1>(Traj);
    double duration = std::get<2>(Traj);

    const double step = duration / (double)NUM_SAMPLES;

    vector<double> xvals;
    vector<double> yvals;
    vector<double> tvals;

    unsigned prev_size = Info.previous_path_x.size();

    if (prev_size >= 2)
    {
        xvals.push_back(Info.previous_path_x[prev_size - 2]);
        yvals.push_back(Info.previous_path_y[prev_size - 2]);
        tvals.push_back(-1 * TIME_STEP);
    }

    for (unsigned i = 1; i <= NUM_SAMPLES; i++)
    {
        double curr_t = step * (double)i;
        tvals.push_back(curr_t);

        double sval = polyeval(s_coeffs, curr_t);
        double dval = polyeval(d_coeffs, curr_t);

        auto xy = Info.getXY(sval, dval);

        xvals.push_back(xy[0]);
        yvals.push_back(xy[1]);
    }

    tk::spline x_spline;
    x_spline.set_points(tvals, xvals);

    tk::spline y_spline;
    y_spline.set_points(tvals, yvals);

    for (int i = 1; i < BUFFER_SIZE - prev_size + 1; i++)
    {
        double curr_t = (double)i * TIME_STEP;

        new_x.push_back(x_spline(curr_t));
        new_y.push_back(y_spline(curr_t));
    }
}

vector<double> polyderiv(const vector<double> &coeffs)
{
    vector<double> newcoeffs;
    for (unsigned i = 1; i < coeffs.size(); i++)
    {
        newcoeffs.push_back(coeffs[i] * double(i));
    }

    return newcoeffs;
}

struct EndPathInfo
{
    double s = 0;
    double sdot = 0;
    double sdotdot = 0;

    double d = 0;
    double ddot = 0;
    double ddotdot = 0;
};

int main() {
    uWS::Hub h;

    // Load up map values for waypoint's x,y,s and d normalized normal vectors
    vector<double> map_waypoints_x;
    vector<double> map_waypoints_y;
    vector<double> map_waypoints_s;
    vector<double> map_waypoints_dx;
    vector<double> map_waypoints_dy;

    // Waypoint map to read from
    string map_file_ = "../data/highway_map.csv";
    // The max s value before wrapping around the track back to 0
    double max_s = 6945.554;

    ifstream in_map_(map_file_.c_str(), ifstream::in);

    string line;
    while (getline(in_map_, line)) {
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
        map_waypoints_x.push_back(x);
        map_waypoints_y.push_back(y);
        map_waypoints_s.push_back(s);
        map_waypoints_dx.push_back(d_x);
        map_waypoints_dy.push_back(d_y);
    }

    auto getXY = [&](double s, double d)
    {
        return ::getXY(s, d,
            map_waypoints_s, map_waypoints_x, map_waypoints_y, map_waypoints_dx, map_waypoints_dy);
    };

    double ref_vel = 0; // mph
    unsigned lane = 1;

    BehaviorState State = KEEP_LANE;
    EndPathInfo EInfo;
    EInfo.d = 2 + 4 * (double)lane;

    h.onMessage(
        [&map_waypoints_x, &map_waypoints_y, &map_waypoints_s, &map_waypoints_dx, &map_waypoints_dy,
        &ref_vel, &lane, &State, &getXY,&EInfo]
    (uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length, uWS::OpCode opCode) {
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
                    const double car_x = j[1]["x"];
                    const double car_y = j[1]["y"];
                    const double car_s = j[1]["s"];
                    const double car_d = j[1]["d"];
                    const double car_yaw = j[1]["yaw"];
                    const double car_speed = j[1]["speed"];

                    // Previous path data given to the Planner
                    const auto previous_path_x = j[1]["previous_path_x"];
                    const auto previous_path_y = j[1]["previous_path_y"];
                    // Previous path's end s and d values 
                    const double end_path_s = j[1]["end_path_s"];
                    const double end_path_d = j[1]["end_path_d"];

                    // Sensor Fusion Data, a list of all other cars on the same side of the road.
                    const std::vector<std::vector<double>> sensor_fusion = j[1]["sensor_fusion"];

                    const SensorData Data =
                        getSensorData(sensor_fusion, car_x, car_y, car_s,
                            car_d, car_yaw, car_speed, ref_vel);

                    const StateInfo Info =
                        getStateInfo(previous_path_x, previous_path_y, end_path_s, end_path_d, getXY);

                    const Context Ctx = getContext(Data, Info);

                    const unsigned prev_size = previous_path_x.size();

                    const double end_path_time = Ctx.Info.end_path_time();

                    if (prev_size == 0)
                        EInfo.s = car_s;

                    /////

                    vector<double> next_x_vals;
                    vector<double> next_y_vals;

                    for (int i = 0; i < prev_size; i++)
                    {
                        next_x_vals.push_back(previous_path_x[i]);
                        next_y_vals.push_back(previous_path_y[i]);
                    }

                    vector<double> new_x;
                    vector<double> new_y;

                    switch (State)
                    {
                    case KEEP_LANE:
                    {
                        auto &E = EInfo;

                        if (E.sdot < 20) 
                        {
                            E.sdot += 0.1;
                        }

                        auto Traj = Trajectory({ E.s, E.sdot }, { E.d }, 3.0 );
                        convertTrajectoryToXY(Traj, Info, new_x, new_y);

                        unsigned num_new_points = new_x.size();
                        double time_projected = (double)num_new_points * TIME_STEP;
                        // sf = s0 + sdot*t
                        E.s += E.sdot * time_projected;

                        break;
                        /*
                        std::cout << "State = KEEP_LANE, lane = " << lane << std::endl;
                        bool ShouldChange = look_for_lane_change(Ctx);
                        vector<Vehicle> blocking_vehicles;
                        if (ShouldChange)
                        {
                            State = PREPARE_CHANGE_LANE;
                            ref_vel -= 0.224;
                        }
                        else if ((lane == 0 || lane == 2) && check_lane_opening(Ctx, lane, blocking_vehicles))
                        {
                            State = CHANGE_LANE;
                        }
                        else
                        {
                            if (ref_vel < 49.8)
                                ref_vel += 0.224;
                        }

                        fill_straight_path(Ctx, ref_vel, lane, new_x, new_y);
                        break;
                        */
                    }
                    case PREPARE_CHANGE_LANE:
                    {
                        std::cout << "State = PREPARE_CHANGE_LANE, lane = " << lane << std::endl;
                        vector<Vehicle> blocking_vehicles;
                        bool LaneOpening = check_lane_opening(Ctx, lane, blocking_vehicles);
                        if (LaneOpening)
                        {
                            State = CHANGE_LANE;
                        }
                        else
                        {
                            double dist;
                            const Vehicle *InFront =
                                get_vehicle_in_front(Ctx, Ctx.Data.ego.get_lane(), dist);

                            if (InFront)
                            {
                                if (dist < 10.0)
                                    ref_vel -= 0.15;
                                else if (ref_vel > InFront->speed() && dist < 20.0)
                                    ref_vel -= 0.15;
                                else if (ref_vel < 49.5 && dist > 30.0)
                                    ref_vel += 0.224;
                            }
                            else
                            {
                                State = KEEP_LANE;
                            }
                        }
                        fill_straight_path(Ctx, ref_vel, lane, new_x, new_y);
                        break;
                    }
                    case CHANGE_LANE:
                    {
                        std::cout << "State = CHANGE_LANE, lane = " << lane << std::endl;
                        double curr_d = Ctx.Data.ego.d;
                        double mid_lane_d = 2 + 4 * lane;
                        double left = mid_lane_d - 0.5;
                        double right = mid_lane_d + 0.5;
                        if (curr_d > left && curr_d < right)
                            State = KEEP_LANE;

                        double dist;
                        const Vehicle *InFront = get_vehicle_in_front(Ctx, lane, dist);

                        if (InFront)
                        {
                            if (dist < 20.0)
                            {
                                if (ref_vel > InFront->speed())
                                {
                                    ref_vel -= 0.15;
                                    std::cout << "slowing..." << std::endl;
                                }
                            }
                            else if (ref_vel < 49.8)
                                ref_vel += 0.224;
                        }
                        else if (ref_vel < 49.8)
                            ref_vel += 0.224;

                        fill_straight_path(Ctx, ref_vel, lane, new_x, new_y);
                        break;
                    }
                    default:
                        assert(0 && "unknown state!");
                    }

                    for (int i = 0; i < new_x.size(); i++)
                    {
                        next_x_vals.push_back(new_x[i]);
                        next_y_vals.push_back(new_y[i]);
                    }

                    // TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
                    json msgJson;
                    msgJson["next_x"] = next_x_vals;
                    msgJson["next_y"] = next_y_vals;

                    auto msg = "42[\"control\"," + msgJson.dump() + "]";

                    //this_thread::sleep_for(chrono::milliseconds(1000));
                    ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);

                }
            }
            else {
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
        }
        else {
            // i guess this should be done more gracefully?
            res->end(nullptr, 0);
        }
    });

    h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
        std::cout << "Connected!!!" << std::endl;
    });

    h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
        char *message, size_t length) {
        ws.close();
        std::cout << "Disconnected" << std::endl;
    });

    int port = 4567;
    if (h.listen(port)) {
        std::cout << "Listening to port " << port << std::endl;
    }
    else {
        std::cerr << "Failed to listen to port" << std::endl;
        return -1;
    }
    h.run();
}
