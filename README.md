# CarND-Controls-MPC
Self-Driving Car Engineer Nanodegree Program

---
## Model

This project uses the simplified equations of motion to predict the vehicle dynamics through time. This forward simulation is in the context of an iterative solver that seeks to minimize an arbitrary cost function. This cost function combines factors to encourage reduced steering and acceleration changes to accomplish the smoothest, smallest steering change to maintain the course.

The equations of motion are:

```
x_[t+1] = x[t] + v[t] * cos(psi[t]) * dt
y_[t+1] = y[t] + v[t] * sin(psi[t]) * dt
psi_[t+1] = psi[t] + v[t] / Lf * delta[t] * dt
v_[t+1] = v[t] + a[t] * dt
cte[t+1] = f(x[t]) - y[t] + v[t] * sin(epsi[t]) * dt
epsi[t+1] = psi[t] - psides[t] + v[t] * delta[t] / Lf * dt
```

Where x, y are the 2D translation coordinates, psi is the steering angle, v is the velocity, cte is the cross track error, and espi is the speed error.
Lf is the length from front to CoG that has a similar radius, which was tuned until the the radius formed by the simulating the model presented in the classroom matched the previous radius.

## Predictive

The simulation is used over N iterations with an adhoc dt which represents the time delta over each interval. 10 iterations with a dt of 0.1 would look ahead 1 second of time. N was tried at 10, 12, 15, 20 iteratiions. N was chosen at 15 steps which gave a balance between too reactive to short term error and the stability given with longer term judgment. Performance is also a factor, as the more time steps simulated will delay the utlitmate steering decision. Dt was tried at 0.05, 0.1, and 0.2 seconds. Dt of 0.05 seconds was chosen emperically to give a fine grain simulation without sacrificing too much distance.


## Controller

The IPOPT optimization library was used with CPPAD auto differentation library to derive a driving solution which minimized oscillation and did not over react to momentary error as the PID can do. Some 100 milliseconds latency was added to the control loop to account for the mechanical delay in the gears and actuators of a real vehicle. To account for this delay, the simulation was contrained over the first few iterations to maintain the previous steering and acceleration command. The first available command after the simulated latency was used as the next vehicle command. Each time step discards the previous predictions and simulates the next N timesteps anew with the latest information.

## Results

![gif](mpc_results.mpc)

Yellow line is a 3rd order polynomial fit through the waypoint lines. The green line plots the N simulated time steps indicicating the results of the solver.

----


## Dependencies

* cmake >= 3.5
 * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1(mac, linux), 3.81(Windows)
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools]((https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)
* [uWebSockets](https://github.com/uWebSockets/uWebSockets)
  * Run either `install-mac.sh` or `install-ubuntu.sh`.
  * If you install from source, checkout to commit `e94b6e1`, i.e.
    ```
    git clone https://github.com/uWebSockets/uWebSockets
    cd uWebSockets
    git checkout e94b6e1
    ```
    Some function signatures have changed in v0.14.x. See [this PR](https://github.com/udacity/CarND-MPC-Project/pull/3) for more details.

* **Ipopt and CppAD:** Please refer to [this document](https://github.com/udacity/CarND-MPC-Project/blob/master/install_Ipopt_CppAD.md) for installation instructions.
* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page). This is already part of the repo so you shouldn't have to worry about it.
* Simulator. You can download these from the [releases tab](https://github.com/udacity/self-driving-car-sim/releases).
* Not a dependency but read the [DATA.md](./DATA.md) for a description of the data sent back from the simulator.


## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./mpc`.

