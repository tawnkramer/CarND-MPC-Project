#include "json.hpp"
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"

/*
*   Completed with help from https://www.youtube.com/&v=bOQuhpz3YfU
*   Self-Driving Car Project Q&A | MPC Controller
*/

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
  auto b2 = s.rfind("}]");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

// Evaluate a polynomial.
double polyeval(Eigen::VectorXd coeffs, double x) {
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

// Evaluate a derivative of polynomial
double polyeval_deriv(Eigen::VectorXd coeffs, double x) {
  double result = 0.0;
  for (int i = 1; i < coeffs.size(); i++) {
    result += coeffs[i] * i * pow(x, i-1);
  }
  return result;
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
                        int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A(xvals.size(), order + 1);

  for (int i = 0; i < xvals.size(); i++) {
    A(i, 0) = 1.0;
  }

  for (int j = 0; j < xvals.size(); j++) {
    for (int i = 0; i < order; i++) {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
}

int main() {
  uWS::Hub h;

  // MPC is initialized here!
  int lookAheadIter = 15;
  double lookAheadDt = 0.05;

  MPC mpc(lookAheadIter, lookAheadDt);

  mpc.SetDesiredVel_MPH(50.0);
  mpc.SetControlLatency_Sec(0.1);
  mpc.SetAccelLimits(1, -1);
  double steer_limit_rad = deg2rad(25.0);
  mpc.SetSteeringLimitAngle_rad(steer_limit_rad);

  h.onMessage([&mpc](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    string sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
      string s = hasData(sdata);
      if (s != "") {
        auto j = json::parse(s);
        string event = j[0].get<string>();
        if (event == "telemetry") {
          // j[1] is the data JSON object
          vector<double> ptsx = j[1]["ptsx"];
          vector<double> ptsy = j[1]["ptsy"];
          double px = j[1]["x"];
          double py = j[1]["y"];
          double psi = j[1]["psi"];
          double v = j[1]["speed"];

          //Transform incoming positions and rotations of car path so that the car starts at the origin
          //and the path extends along the x axis. This helps with the fit function because we have one and only
          //one value for y for every x. And the starting rotation is also zero - straight down the x axis.
          //Will help MPC evaluation and cross track error.
          vector<double> pos_x(ptsx.size());
          vector<double> pos_y(ptsx.size());
          for (int i=0; i<ptsx.size(); i++) {
            double shift_x = ptsx[i] - px;
            double shift_y = ptsy[i] - py;

            pos_x[i] = shift_x * cos (-psi) - shift_y * sin (-psi);
            pos_y[i] = shift_x * sin (-psi) + shift_y * cos (-psi);
          }

          double* ptrx = &pos_x[0];
          Eigen::Map<Eigen::VectorXd> ptsx_transform(ptrx, pos_x.size());

          double* ptry = &pos_y[0];
          Eigen::Map<Eigen::VectorXd> ptsy_transform(ptry, pos_y.size());

          const int order_polynimial = 3;
          auto coeffs = polyfit(ptsx_transform, ptsy_transform, order_polynimial);

          // in car frame coords of the car is (0, 0) with psi == 0
          double cte = polyeval(coeffs, 0);

          //full epsi can be simplified because psi is zero and px is zero
          //epsi = psi - atan(coeff[1] + 2 * px * coeffs[2] + 3 coeffs[3] * pow(px, 2))
          //double epsi = atan(polyeval_deriv (coeffs, 0));
          double epsi = -atan(coeffs[1]);

          Eigen::VectorXd state(6);

          //x , y, theta are all zero because of above trasnformations
          state << 0, 0, 0, v, cte, epsi;

          auto vars = mpc.Solve(state, coeffs);

          double steer_value = vars [0];
          double throttle_value = vars [1];

          json msgJson;
          
          // NOTE: Remember to divide by deg2rad(25) before you send the steering value back.
          // Otherwise the values will be in between [-deg2rad(25), deg2rad(25] instead of [-1, 1].
          msgJson["steering_angle"] = steer_value;
          msgJson["throttle"] = throttle_value;

          //Display the MPC predicted trajectory 
          vector<double> mpc_x_vals;
          vector<double> mpc_y_vals;

          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Green line
          //skip the first two vars which are our solution steering and throttle
          for(int i = 2; i < vars.size(); i++)
          {
              if(i % 2 == 0)
              {
                mpc_x_vals.push_back(vars[i]);
              }
              else
              {
                mpc_y_vals.push_back(vars[i]);
              }
          }

          msgJson["mpc_x"] = mpc_x_vals;
          msgJson["mpc_y"] = mpc_y_vals;

          //Display the waypoints/reference line
          vector<double> next_x_vals;
          vector<double> next_y_vals;

          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Yellow line

          double num_points = mpc.m_LookAheadIter;

          for(int i = 1; i < num_points; i++)
          {
              double x_val = i * 5;
              next_x_vals.push_back(x_val);
              next_y_vals.push_back(polyeval(coeffs, x_val));
          }

          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;


          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
          //std::cout << msg << std::endl;
          // Latency
          // The purpose is to mimic real driving conditions where
          // the car does actuate the commands instantly.
          //
          // Feel free to play around with this value but should be to drive
          // around the track with 100ms latency.
          //
          // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE
          // SUBMITTING.
          this_thread::sleep_for(chrono::milliseconds(100));
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
    ws.close();
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
