## Dependencies

* Python3
* OMP
* R

## How to run
Premise: The code for generating plots requires the python interpreter and two libraries to be installed, matplotlib, and pandas with the following commands:
pip3 install matplotlib
pip install pandas

When generating tables with R a package will be installed, known as kableExtra and webshot.

1.	Navigate to the project folder containing the makefile

2.	To clear previously obtained achievements and previous builds, enter the command
make clean

3.	(OMP+MPI) To generate the necessary directories and compile and linking the various source codes, enter the command
make all

4.	(OMP+MPI) To run the algorithm for making omp+mpi tests and producing results enter the command
make mpitest

5.	(OMP+MPI) To produce plots and tables enter the command
make mpi_plots_tables

6.	(OMP+CUDA) To compile and linking the various source codes, enter the command
make cudacompile

7.	(OMP+CUDA) To run the algorithm for making omp+cuda tests and producing results enter the command
make cudatest

8.	(OMP+CUDA) To produce plots and tables enter the command
make cuda_plots_tables

The results of the algorithms can be viewed in the "RB_Search_Report" folder, which has three subfolders depending on the size of the problem (the nodes of the tree).
The execution times of the algorithms and their average values can be viewed respectively in the "Informations" and "Results" folders, divided by optimization.
The results in graphical and tabular format can be viewed respectively in the "Plots" and "Tables" folders, divided by optimization.
To reduce the completion time, change the value of "iterations" in the makefile to the desired number (each optimization iterations lasts about 20 minutes).
