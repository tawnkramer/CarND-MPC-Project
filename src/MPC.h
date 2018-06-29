#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"

using namespace std;

class MPC {
 public:
  MPC();

  virtual ~MPC();

  void SetDesiredVel_MPH(float vel) { m_DesiredVel_MPH = vel; }
  void SetControlLatency_Sec(float sec) { m_ControlLatency_Sec = sec; }

  float m_DesiredVel_MPH;
  float m_ControlLatency_Sec;

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuatotions.
  vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);
};

#endif /* MPC_H */
