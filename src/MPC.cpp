#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  double ref_v;
  double dt;
  MPC& mpc;

  FG_eval(Eigen::VectorXd coeffs, 
    double desired_vel, 
    double latency,
    MPC& _mpc) : mpc(_mpc) 
  {
    //we are assuming a third order polynomial
    assert(coeffs.size() == 4);
    this->coeffs = coeffs; 
    this->ref_v = desired_vel;
    this->dt = latency;
  }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)

    // The cost is stored is the first element of `fg`.
    // Any additions to the cost should be added to `fg[0]`.
    fg[0] = 0;

    double ref_cte = 0;
    double ref_epsi = 0;

    // The part of the cost based on the reference state.
    for (int t = 0; t < mpc.m_LookAheadIter; t++) {
      fg[0] += 2000 * CppAD::pow(vars[mpc.cte_start + t] - ref_cte, 2);
      fg[0] += 2000 * CppAD::pow(vars[mpc.epsi_start + t] - ref_epsi, 2);
      fg[0] += CppAD::pow(vars[mpc.v_start + t] - ref_v, 2);
    }

    // Minimize the use of actuators.
    for (int t = 0; t < mpc.m_LookAheadIter - 1; t++) {
      fg[0] += 5 * CppAD::pow(vars[mpc.delta_start + t], 2);
      fg[0] += 5 * CppAD::pow(vars[mpc.a_start + t], 2);
    }

    // Minimize the value gap between sequential actuations.
    for (int t = 0; t < mpc.m_LookAheadIter - 2; t++) {
      fg[0] += 200 * CppAD::pow(vars[mpc.delta_start + t + 1] - vars[mpc.delta_start + t], 2);
      fg[0] += 10 * CppAD::pow(vars[mpc.a_start + t + 1] - vars[mpc.a_start + t], 2);
    }

    //
    // Setup Constraints
    //
    // NOTE: In this section you'll setup the model constraints.

    // Initial constraints
    //
    // We add 1 to each of the starting indices due to cost being located at
    // index 0 of `fg`.
    // This bumps up the position of all the other values.
    fg[1 + mpc.x_start] = vars[mpc.x_start];
    fg[1 + mpc.y_start] = vars[mpc.y_start];
    fg[1 + mpc.psi_start] = vars[mpc.psi_start];
    fg[1 + mpc.v_start] = vars[mpc.v_start];
    fg[1 + mpc.cte_start] = vars[mpc.cte_start];
    fg[1 + mpc.epsi_start] = vars[mpc.epsi_start];

    // The rest of the constraints
    for (int t = 1; t < mpc.m_LookAheadIter; t++) {
      // The state at time t+1 .
      AD<double> x1 = vars[mpc.x_start + t];
      AD<double> y1 = vars[mpc.y_start + t];
      AD<double> psi1 = vars[mpc.psi_start + t];
      AD<double> v1 = vars[mpc.v_start + t];
      AD<double> cte1 = vars[mpc.cte_start + t];
      AD<double> epsi1 = vars[mpc.epsi_start + t];

      // The state at time t.
      AD<double> x0 = vars[mpc.x_start + t - 1];
      AD<double> y0 = vars[mpc.y_start + t - 1];
      AD<double> psi0 = vars[mpc.psi_start + t - 1];
      AD<double> v0 = vars[mpc.v_start + t - 1];
      AD<double> cte0 = vars[mpc.cte_start + t - 1];
      AD<double> epsi0 = vars[mpc.epsi_start + t - 1];

      // Only consider the actuation at time t.
      AD<double> delta0 = vars[mpc.delta_start + t - 1];
      AD<double> a0 = vars[mpc.a_start + t - 1];

      //assumes a third order polynomial passed in..
      AD<double> f0 = coeffs[0] + coeffs[1] * x0  + coeffs[2] * x0 * x0  + coeffs[3] * x0 * x0 * x0;
      AD<double> psides0 = CppAD::atan(3 * coeffs[3] * x0 * x0 + 2 * coeffs[2] * x0 + coeffs[1]);

      // Here's `x` to get you started.
      // The idea here is to constraint this value to be 0.
      //
      // Recall the equations for the model:
      // x_[t+1] = x[t] + v[t] * cos(psi[t]) * dt
      // y_[t+1] = y[t] + v[t] * sin(psi[t]) * dt
      // psi_[t+1] = psi[t] + v[t] / Lf * delta[t] * dt
      // v_[t+1] = v[t] + a[t] * dt
      // cte[t+1] = f(x[t]) - y[t] + v[t] * sin(epsi[t]) * dt
      // epsi[t+1] = psi[t] - psides[t] + v[t] * delta[t] / Lf * dt
      fg[1 + mpc.x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[1 + mpc.y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
      fg[1 + mpc.psi_start + t] = psi1 - (psi0 + v0 * delta0 / mpc.Lf * dt);
      fg[1 + mpc.v_start + t] = v1 - (v0 + a0 * dt);
      fg[1 + mpc.cte_start + t] = cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
      fg[1 + mpc.epsi_start + t] = epsi1 - ((psi0 - psides0) - v0 * delta0 / mpc.Lf * dt);
    }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC(int lookaheadIter, double lookaheadDt) 
{
  SetLookAheadIter(lookaheadIter);
  SetLookAheadDt(lookaheadDt);

  //Reasonable defaults
  m_DesiredVel_MPH = 30;
  m_ControlLatency_Sec = 0.1;
  m_LookAheadIter = 10;
  m_LookAheadDt_Sec = 0.5;
  m_SteerLimit_Rad = 0.5;
  m_AccelUpperLimit = 1;
  m_AccelLowerLimit = -1;

  // The solver takes all the state variables and actuator
  // variables in a singular vector. Thus, we should to establish
  // when one variable starts and another ends to make our life easier.

  int N = lookaheadIter;

  x_start = 0;
  y_start = x_start + N;
  psi_start = y_start + N;
  v_start = psi_start + N;
  cte_start = v_start + N;
  epsi_start = cte_start + N;
  delta_start = epsi_start + N;
  a_start = delta_start + N - 1;
}

MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  double x = state[0];
  double y = state[1];
  double psi = state[2];
  double v = state[3];
  double cte = state[4];
  double epsi = state[5];

  const int num_vars_state = 6;

  // number of independent variables
  // N timesteps == N - 1 actuations
  size_t n_vars = m_LookAheadIter * num_vars_state + (m_LookAheadIter - 1) * 2;
  
  // Number of constraints
  size_t n_constraints = m_LookAheadIter * num_vars_state;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) 
  {
    vars[i] = 0;
  }

// Set the initial variable values
  vars[x_start] = x;
  vars[y_start] = y;
  vars[psi_start] = psi;
  vars[v_start] = v;
  vars[cte_start] = cte;
  vars[epsi_start] = epsi;


  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);

  // Set all non-actuators upper and lowerlimits
  // to the max negative and positive values.
  for (int i = 0; i < delta_start; i++) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }

  // The upper and lower limits of delta are set to -25 and 25
  // degrees (values in radians).
  for (int i = delta_start; i < a_start; i++) {
    vars_lowerbound[i] = -m_SteerLimit_Rad * Lf;
    vars_upperbound[i] = m_SteerLimit_Rad * Lf;
  }

  // Acceleration/decceleration upper and lower limits.
  for (int i = a_start; i < n_vars; i++) {
    vars_lowerbound[i] = m_AccelLowerLimit;
    vars_upperbound[i] = m_AccelUpperLimit;
  }

  // to handle delay before we actually may control the car
  // we need to freeze control commands with previous values
  // during latency delay
  double latency_ = m_ControlLatency_Sec;
  double dt = 0.1;
  int num_iter_delayed = int(std::ceil(latency_ / dt));
  for (int i=0; i < num_iter_delayed; i++) {
    vars_lowerbound [delta_start + i] = prev_actuations [0];
    vars_upperbound [delta_start + i] = prev_actuations [0];
    vars_lowerbound [a_start + i] = prev_actuations [1];
    vars_upperbound [a_start + i] = prev_actuations [1];
  }

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  //Tell the solver to set the initial state to our init values
  constraints_lowerbound[x_start] = x;
  constraints_lowerbound[y_start] = y;
  constraints_lowerbound[psi_start] = psi;
  constraints_lowerbound[v_start] = v;
  constraints_lowerbound[cte_start] = cte;
  constraints_lowerbound[epsi_start] = epsi;

  constraints_upperbound[x_start] = x;
  constraints_upperbound[y_start] = y;
  constraints_upperbound[psi_start] = psi;
  constraints_upperbound[v_start] = v;
  constraints_upperbound[cte_start] = cte;
  constraints_upperbound[epsi_start] = epsi;

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs, m_DesiredVel_MPH, m_ControlLatency_Sec, *this);

  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  //std::cout << "Cost " << cost << std::endl;

  vector<double> result;

  double steering = solution.x[delta_start + num_iter_delayed];
  double throttle = solution.x[a_start + num_iter_delayed];

  //remember our last command,s
  prev_actuations [0] = steering;
  prev_actuations [1] = throttle;

  result.push_back(steering / 0.436332); //our best guess as steering
  result.push_back(throttle); //our best guess as accelleration

  for(int i = 0; i < m_LookAheadIter-1; i++)
  {
      result.push_back(solution.x[x_start + i + 1]);
      result.push_back(solution.x[y_start + i + 1]);
  }

  return result;
}
