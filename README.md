# Charged Particle Simulator

The accurate prediction of relativistic particle trajectories in electromagnetic (EM) fields is vital for scientific and technological applications in astrophysics, plasma physics, and accelerator physics. This work presents a numerical study of relativistic particle motion in arbitrary EM fields, employing the Boris algorithm for its widespread use and effectiveness in solving equations of motion. By comparing computed results with analytical solutions, we validate the algorithm's accuracy and demonstrate its ability to accurately integrate the equations in a time-accurate and stable manner.



## How to run

- Run `create_project.sh`
- Input the desired project's name. A folder will be created with that name containing the cpp, gp and run.sh files
- Modify the `project_name.cpp` file to satisfy your simulation's needs. In particular, you will want to modify the initial and final time, the definitions
of the E and B fields and the initial conditions
- Run `run.sh` inside the terminal. You will now see the `output.dat` file and the `plot.png` file

