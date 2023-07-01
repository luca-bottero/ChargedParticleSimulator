g++ -fopenmp gyration.cpp ../ode_solvers.cpp --output gyration.out; ./gyration.out; gnuplot gyration.gp; gnuplot 3d_trajectory.gp
