#include "json.hpp"
#include <zmq.hpp>
#include <math.h>
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

void show_usage()
{
    printf("Usage: mpc_zmq [--port <int>] [--iter <int>] [--dt <float>] [--mph <float>] [--latency <float>] [--steer_limit <float>] [--debug] [--help] \n\n");
		printf("port:       os port to bind socket. range(1024-65535)                            default: 5555\n");
		printf("iter:       number of steps to forward simulate car path. range(1-100)           default: 15\n");
		printf("dt:         time between steps of simulated car path. range(0.0001 - 1.0)        default: 0.1 seconds\n");
		printf("mph:        desired speed in miles per hour (0.001, 100)                         default: 10 mph\n");
		printf("latency:    estimated delay in control before change applied (0.001, 100)        default: 0.1 seconds\n");
		printf("steer_limit: degree range in steering (10, 90)                                   default: 15 degrees\n");
    printf("debug:      emit additional output in json for path info.                        default: false\n\n");
		printf("Summary:\n");
		printf("This will use zero message queue to listen for json telemetry, then solve the given path for the optimal next steering and throttle controls.\n");
		printf("The json input telemetry should take the form { ptsx[], ptsy[], x, y, psi, speed }\n");
		printf("The json return takes the form { steering, throttle } \n");
}

int main(int argc, char** argv) {
  
		// MPC default parameters
		int port = 5555;		
		int lookAheadIter = 15;
		double lookAheadDt = 0.05;
		float mph = 50.0f;
		float latency = 0.1f;
		float steer_limit = 15.0f;
    bool debug = false;

		for(int i = 1; i < argc; i++)
		{
				char* arg = argv[i];
				
				if (strcmp(arg, "--port") == 0 && i < argc - 1)
				{
				    port = atoi(argv[i + 1]); i++;
				}
				else if (strcmp(arg, "--iter") == 0 && i < argc - 1)
				{
				    lookAheadIter = atoi(argv[i + 1]); i++;
				}
				else if (strcmp(arg, "--dt") == 0 && i < argc - 1)
				{
				    lookAheadDt = atof(argv[i + 1]); i++;
				}
				else if (strcmp(arg, "--mph") == 0 && i < argc - 1)
				{
				    mph = atof(argv[i + 1]); i++;
				}
				else if (strcmp(arg, "--latency") == 0 && i < argc - 1)
				{
				    latency = atof(argv[i + 1]); i++;
				}
				else if (strcmp(arg, "--steer_limit") == 0 && i < argc - 1)
				{
				    steer_limit = atof(argv[i + 1]); i++;
				}
				else if (strcmp(arg, "--debug") == 0)
				{
				    debug = true;
				}
				else if (strcmp(arg, "--help") == 0 || strcmp(arg, "-h") == 0)
				{
				    show_usage();
				    exit(1);
				}
				else
				{
						printf("Error>> Unrecognized parameter: %s\n", argv[i]);
				    show_usage();
				    exit(1);
				}
		
		}

	 //  Prepare our context and socket
    zmq::context_t context (1);
    zmq::socket_t socket (context, ZMQ_REP);
    char connect_str[256];
    sprintf(connect_str, "tcp://*:%d", port);
    socket.bind (connect_str);

    printf("listening for mpc requests on port: %d\n", port);

		MPC mpc(lookAheadIter, lookAheadDt);

		mpc.SetDesiredVel_MPH(mph);
		mpc.SetControlLatency_Sec(latency);
		mpc.SetAccelLimits(1, -1);
		double steer_limit_rad = deg2rad(steer_limit);
		mpc.SetSteeringLimitAngle_rad(steer_limit_rad);

	 	while (true)
		{
			zmq::message_t request;

			//  Wait for next request from client
			socket.recv (&request);

			char* data = (char*)request.data();
			size_t length = request.size();

      string sdata = string(data).substr(0, length);
      auto j = json::parse(sdata);
      //cout << j << std::endl;
      
        // j[0] is the data JSON object
        vector<double> ptsx = j["ptsx"];
        vector<double> ptsy = j["ptsy"];
        double px = j["x"];
        double py = j["y"];
        double psi = j["psi"];
        double v = j["speed"];

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

        //x , y, theta are all zero because of above transformations
        state << 0, 0, 0, v, cte, epsi;

        auto vars = mpc.Solve(state, coeffs);

        double steer_value = vars [0];
        double throttle_value = vars [1];

        json msgJson;
        
        // NOTE: Remember to divide by deg2rad(25) before you send the steering value back.
        // Otherwise the values will be in between [-deg2rad(25), deg2rad(25] instead of [-1, 1].
        msgJson["steering"] = steer_value;
        msgJson["throttle"] = throttle_value;


        if(debug)
        {
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
        }
       
        auto msg = msgJson.dump();

				//cout << msg << "\n";
 				
				//  Send reply back to client
        zmq::message_t reply (msg.size());
        memcpy (reply.data (), msg.c_str(), msg.size());
        socket.send (reply);
  }

	return 0;
}

