#!/bin/bash

read -p "Enter the project name: " project_name

mkdir "$project_name"

cp template/template.cpp "$project_name/$project_name.cpp"

cp template/template.gp "$project_name/$project_name.gp"

echo "g++ -fopenmp $project_name.cpp ../ode_solvers.cpp --output $project_name.out; ./$project_name.out; gnuplot $project_name.gp" > "$project_name/run.sh"

chmod +x "$project_name/run.sh"

echo "Project '$project_name' created successfully."
