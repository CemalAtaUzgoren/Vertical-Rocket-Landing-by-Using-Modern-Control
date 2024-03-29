# Rocket Control Repository\
\
This repository contains MATLAB scripts and models for rocket control using Linear Quadratic Regulator (LQR) and Model Predictive Control (MPC) techniques.\
\
## LQR Control\
\
### Files\
- `LQR_control.asv`\
- `LQR_control.m`\
- `lqr_controller_simulation_for_rocket_landing.slx`\
- `Rocket.stl`\
- `slprj/sim/varcache/lqr_controller_simulation_for_rocket_landing/checksumOfCache.mat`\
- `slprj/sim/varcache/lqr_controller_simulation_for_rocket_landing/tmwinternal/simulink_cache.xml`\
- `slprj/sim/varcache/lqr_controller_simulation_for_rocket_landing/varInfo.mat`\
- `untitled.slx.autosave`\
\
### How to Use\
1. Run the `lqr_controller_simulation_for_rocket_landing.slx` model first.\
2. After running the simulation, execute the `LQR_control.m` script.\
\
### LQR (Linear Quadratic Regulator)\
LQR is a control strategy designed to minimize the quadratic cost of a linear system subject to control input constraints. In this context, it is applied to rocket control for landing.\
\
## MPC (Model Predictive Control)\
\
### Files\
- `matrix_riccati.asv`\
- `matrix_riccati.m`\
- `MPC_main_Rocket_landing.m`\
- `response.m`\
- `Rocket.stl`\
\
### How to Use\
1. Run the `MPC_main_Rocket_landing.m` script.\
\
### MPC (Model Predictive Control)\
MPC is an advanced control strategy that utilizes a predictive model of the system to optimize control inputs over a finite time horizon. This repository implements MPC for rocket landing control.\
\
## Oscillations and System Improvements (MPC Section)\
The simulations on the Mac platform may exhibit excessive oscillations. However, it is important to note that these oscillations can be addressed and improved through further developments. The challenges in resolving these oscillations are attributed to singularities in the system.\
\
Feel free to explore, modify, and use these scripts and models in your projects. If you encounter any issues or have suggestions for improvements, please create an issue or submit a pull request.\
\
Happy coding!\
