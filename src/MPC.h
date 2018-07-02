#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"

using namespace std;

class MPC {
 public:
  MPC(int lookaheadIter, double lookaheadDt);

  virtual ~MPC();

  void SetDesiredVel_MPH(double vel) { m_DesiredVel_MPH = vel; }
  void SetControlLatency_Sec(double sec) { m_ControlLatency_Sec = sec; }
  void SetLookAheadIter(int i) { m_LookAheadIter = i; }
  void SetLookAheadDt(double dt) { m_LookAheadDt_Sec = dt; }
  void SetSteeringLimitAngle_rad(double angle) { m_SteerLimit_Rad = angle; }
  void SetAccelLimits(double upperLimit, double lowerLimit)
  { 
    m_AccelUpperLimit = upperLimit;
    m_AccelLowerLimit = lowerLimit;
  }

  double  m_DesiredVel_MPH;
  double  m_ControlLatency_Sec;
  int     m_LookAheadIter;
  double  m_LookAheadDt_Sec;
  double  m_SteerLimit_Rad;
  double  m_AccelUpperLimit;
  double  m_AccelLowerLimit;

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuatotions.
  vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);

  double prev_actuations[2];

  size_t x_start;
  size_t y_start;
  size_t psi_start;
  size_t v_start;
  size_t cte_start;
  size_t epsi_start;
  size_t delta_start;
  size_t a_start;

  // Lf was tuned until the the radius formed by the simulating the model
  // presented in the classroom matched the previous radius.
  // This is the length from front to CoG that has a similar radius.
  double Lf = 2.67;

};

#endif /* MPC_H */
