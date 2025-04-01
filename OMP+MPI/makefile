.PHONY: all clean test compile0 compile1 compile2 compile3 mpitest0 mpitest1 mpitest2 mpitest3 cudatest cudatest0 cudatest1 cudatest2 cudatest3 statistics table_plots

all: compile0 compile1 compile2 compile3 create_directories

clean:
	rm -rf ./Build/*
	rm -rf ./Informations
	rm -rf ./Results
	rm -rf ./Tables
	rm -rf ./Plots
	rm -rf ./RB_Search_Report/*

mpitest: mpitest0 mpitest1 mpitest2 mpitest3 mpistatistics

mpi_plots_tables: mpiplots mpitables

cudatest: cudatest0 cudatest1 cudatest2 cudatest3 cudastatistics

cuda_plots_tables: cudaplots cudatables

n_nodestotal := 10000 250000 3500000
n_nodes1 := 10000
n_nodes2 := 250000
n_nodes3 := 3500000
num_to_search1 := 322 3456
num_to_search2:= 2340 175430
num_to_search3 := 676330 3234567
omp_threads := 1 2 4 8 16
omp_threads_cuda := 1 2 4
mpi_process := 1 2 4 8
seed = 1
iterations := 1 2

create_directories:
	@echo "Creating directories..."
	@mkdir -p Informations
	@cd Informations && mkdir -p OMP_MPI && cd OMP_MPI && \
		mkdir -p opt0 && \
    mkdir -p opt1 && \
    mkdir -p opt2 && \
    mkdir -p opt3
	@cd Informations && mkdir -p OMP_CUDA && cd OMP_CUDA && \
		mkdir -p opt0 && \
    mkdir -p opt1 && \
    mkdir -p opt2 && \
    mkdir -p opt3
	@mkdir -p RB_Search_Report
	@cd RB_Search_Report && \
		for node_size in $(n_nodestotal); do \
			mkdir -p nodes$$node_size; \
		done;
	@mkdir -p Results
	@cd Results && mkdir -p OMP_MPI && cd OMP_MPI && \
		mkdir -p opt0 && \
    mkdir -p opt1 && \
    mkdir -p opt2 && \
    mkdir -p opt3
	@cd Results && mkdir -p OMP_CUDA && cd OMP_CUDA && \
		mkdir -p opt0 && \
    mkdir -p opt1 && \
    mkdir -p opt2 && \
    mkdir -p opt3
	@mkdir -p Tables
	@cd Tables && mkdir -p OMP_MPI && cd OMP_MPI && \
		mkdir -p opt0 && \
    mkdir -p opt1 && \
    mkdir -p opt2 && \
    mkdir -p opt3
	@cd Tables && mkdir -p OMP_CUDA && cd OMP_CUDA && \
		mkdir -p opt0 && \
    mkdir -p opt1 && \
    mkdir -p opt2 && \
    mkdir -p opt3
	@mkdir -p Plots
	@cd Plots && mkdir -p OMP_MPI && cd OMP_MPI && \
		mkdir -p opt0 && \
    mkdir -p opt1 && \
    mkdir -p opt2 && \
    mkdir -p opt3
	@cd Plots && mkdir -p OMP_CUDA && cd OMP_CUDA && \
		mkdir -p opt0 && \
    mkdir -p opt1 && \
    mkdir -p opt2 && \
    mkdir -p opt3
	@echo "Directories created!"

compile0:
	gcc -c -o ./Build/RB_Tree_Generator.o ./Source/RB_Tree_Generator.c -O0
	gcc -c -o ./Build/RB_Tree_Sequential.o ./Source/RB_Tree_Sequential.c -O0
	mpicc -c -o ./Build/RB_Tree_OMP_MPI.o ./Source/RB_Tree_OMP_MPI.c -fopenmp -O0

	gcc -o ./Build/RB_Tree_Sequential0.exe ./Build/RB_Tree_Sequential.o ./Build/RB_Tree_Generator.o -O0
	mpicc -o ./Build/RB_Tree_OMP_MPI0.exe ./Build/RB_Tree_OMP_MPI.o ./Build/RB_Tree_Generator.o -fopenmp -O0

compile1:
	gcc -c -o ./Build/RB_Tree_Generator.o ./Source/RB_Tree_Generator.c -O1
	gcc -c -o ./Build/RB_Tree_Sequential.o ./Source/RB_Tree_Sequential.c -O1
	mpicc -c -o ./Build/RB_Tree_OMP_MPI.o ./Source/RB_Tree_OMP_MPI.c -fopenmp -O1

	gcc -o ./Build/RB_Tree_Sequential1.exe ./Build/RB_Tree_Sequential.o ./Build/RB_Tree_Generator.o -O1
	mpicc -o ./Build/RB_Tree_OMP_MPI1.exe ./Build/RB_Tree_OMP_MPI.o ./Build/RB_Tree_Generator.o -fopenmp -O1

compile2:
	gcc -c -o ./Build/RB_Tree_Generator.o ./Source/RB_Tree_Generator.c -O2
	gcc -c -o ./Build/RB_Tree_Sequential.o ./Source/RB_Tree_Sequential.c -O2
	mpicc -c -o ./Build/RB_Tree_OMP_MPI.o ./Source/RB_Tree_OMP_MPI.c -fopenmp -O2

	gcc -o ./Build/RB_Tree_Sequential2.exe ./Build/RB_Tree_Sequential.o ./Build/RB_Tree_Generator.o -O2
	mpicc -o ./Build/RB_Tree_OMP_MPI2.exe ./Build/RB_Tree_OMP_MPI.o ./Build/RB_Tree_Generator.o -fopenmp -O2

compile3:
	gcc -c -o ./Build/RB_Tree_Generator.o ./Source/RB_Tree_Generator.c -O3
	gcc -c -o ./Build/RB_Tree_Sequential.o ./Source/RB_Tree_Sequential.c -O3
	mpicc -c -o ./Build/RB_Tree_OMP_MPI.o ./Source/RB_Tree_OMP_MPI.c -fopenmp -O3

	gcc -o ./Build/RB_Tree_Sequential3.exe ./Build/RB_Tree_Sequential.o ./Build/RB_Tree_Generator.o -O3
	mpicc -o ./Build/RB_Tree_OMP_MPI3.exe ./Build/RB_Tree_OMP_MPI.o ./Build/RB_Tree_Generator.o -fopenmp -O3

mpistatistics:
	gcc -c -o ./Build/Callingpy.o ./Source/Callingpy.c 
	gcc -o ./Build/Callingpy.exe ./Build/Callingpy.o
	./Build/Callingpy.exe 0

mpiplots:
	./Build/Callingpy.exe 2

mpitables:
	Rscript ./Source/maketables.R

mpitest0:
	@opt=0; \
	for node in $(n_nodes1); do \
		for it in $(iterations); do \
			for valuetosearch in $(num_to_search1); do \
				echo "----------------------------{OMP+MPI: $$node,$$opt,$$it,$$valuetosearch}----------------------------"; \
				./Build/RB_Tree_Sequential0.exe $$node $(seed) $$opt $$valuetosearch; \
				for omp in $(omp_threads); do \
					for mpi in $(mpi_process); do \
						mpiexec -n $$mpi ./Build/RB_Tree_OMP_MPI0.exe $$node $(seed) $$opt $$omp $$valuetosearch -fopenmp; \
					done; \
				done; \
			done; \
		done; \
	done; \
	for node in $(n_nodes2); do \
		for it in $(iterations); do \
			for valuetosearch in $(num_to_search2); do \
				echo "----------------------------{OMP+MPI: $$node,$$opt,$$it,$$valuetosearch}----------------------------"; \
				./Build/RB_Tree_Sequential0.exe $$node $(seed) $$opt $$valuetosearch; \
				for omp in $(omp_threads); do \
					for mpi in $(mpi_process); do \
						mpiexec -n $$mpi ./Build/RB_Tree_OMP_MPI0.exe $$node $(seed) $$opt $$omp $$valuetosearch -fopenmp; \
					done; \
				done; \
			done; \
		done; \
	done; \
	for node in $(n_nodes3); do \
		for it in $(iterations); do \
			for valuetosearch in $(num_to_search3); do \
				echo "----------------------------{OMP+MPI: $$node,$$opt,$$it,$$valuetosearch}----------------------------"; \
				./Build/RB_Tree_Sequential0.exe $$node $(seed) $$opt $$valuetosearch; \
				for omp in $(omp_threads); do \
					for mpi in $(mpi_process); do \
						mpiexec -n $$mpi ./Build/RB_Tree_OMP_MPI0.exe $$node $(seed) $$opt $$omp $$valuetosearch -fopenmp; \
					done; \
				done; \
			done; \
		done; \
	done;

mpitest1:
	@opt=1; \
	for node in $(n_nodes1); do \
		for it in $(iterations); do \
			for valuetosearch in $(num_to_search1); do \
				echo "----------------------------{OMP+MPI: $$node,$$opt,$$it,$$valuetosearch}----------------------------"; \
				./Build/RB_Tree_Sequential1.exe $$node $(seed) $$opt $$valuetosearch; \
				for omp in $(omp_threads); do \
					for mpi in $(mpi_process); do \
						mpiexec -n $$mpi ./Build/RB_Tree_OMP_MPI1.exe $$node $(seed) $$opt $$omp -fopenmp $$valuetosearch; \
					done; \
				done; \
			done; \
		done; \
	done; \
	for node in $(n_nodes2); do \
		for it in $(iterations); do \
			for valuetosearch in $(num_to_search2); do \
				echo "----------------------------{OMP+MPI: $$node,$$opt,$$it,$$valuetosearch}----------------------------"; \
				./Build/RB_Tree_Sequential1.exe $$node $(seed) $$opt $$valuetosearch; \
				for omp in $(omp_threads); do \
					for mpi in $(mpi_process); do \
						mpiexec -n $$mpi ./Build/RB_Tree_OMP_MPI1.exe $$node $(seed) $$opt $$omp $$valuetosearch -fopenmp; \
					done; \
				done; \
			done; \
		done; \
	done; \
	for node in $(n_nodes3); do \
		for it in $(iterations); do \
			for valuetosearch in $(num_to_search3); do \
				echo "----------------------------{OMP+MPI: $$node,$$opt,$$it,$$valuetosearch}----------------------------"; \
				./Build/RB_Tree_Sequential1.exe $$node $(seed) $$opt $$valuetosearch; \
				for omp in $(omp_threads); do \
					for mpi in $(mpi_process); do \
						mpiexec -n $$mpi ./Build/RB_Tree_OMP_MPI1.exe $$node $(seed) $$opt $$omp $$valuetosearch -fopenmp; \
					done; \
				done; \
			done; \
		done; \
	done;

mpitest2:
	@opt=2; \
	for node in $(n_nodes1); do \
		for it in $(iterations); do \
			for valuetosearch in $(num_to_search1); do \
				echo "----------------------------{OMP+MPI: $$node,$$opt,$$it,$$valuetosearch}----------------------------"; \
				./Build/RB_Tree_Sequential2.exe $$node $(seed) $$opt $$valuetosearch; \
				for omp in $(omp_threads); do \
					for mpi in $(mpi_process); do \
						mpiexec -n $$mpi ./Build/RB_Tree_OMP_MPI2.exe $$node $(seed) $$opt $$omp $$valuetosearch -fopenmp; \
					done; \
				done; \
			done; \
		done; \
	done; \
	for node in $(n_nodes2); do \
		for it in $(iterations); do \
			for valuetosearch in $(num_to_search2); do \
				echo "----------------------------{OMP+MPI: $$node,$$opt,$$it,$$valuetosearch}----------------------------"; \
				./Build/RB_Tree_Sequential2.exe $$node $(seed) $$opt $$valuetosearch; \
				for omp in $(omp_threads); do \
					for mpi in $(mpi_process); do \
						mpiexec -n $$mpi ./Build/RB_Tree_OMP_MPI2.exe $$node $(seed) $$opt $$omp $$valuetosearch -fopenmp; \
					done; \
				done; \
			done; \
		done; \
	done; \
	for node in $(n_nodes3); do \
		for it in $(iterations); do \
			for valuetosearch in $(num_to_search3); do \
				echo "----------------------------{OMP+MPI: $$node,$$opt,$$it,$$valuetosearch}----------------------------"; \
				./Build/RB_Tree_Sequential2.exe $$node $(seed) $$opt $$valuetosearch; \
				for omp in $(omp_threads); do \
					for mpi in $(mpi_process); do \
						mpiexec -n $$mpi ./Build/RB_Tree_OMP_MPI2.exe $$node $(seed) $$opt $$omp $$valuetosearch -fopenmp; \
					done; \
				done; \
			done; \
		done; \
	done;

mpitest3:
	@opt=3; \
	for node in $(n_nodes1); do \
		for it in $(iterations); do \
			for valuetosearch in $(num_to_search1); do \
				echo "----------------------------{OMP+MPI: $$node,$$opt,$$it,$$valuetosearch}----------------------------"; \
				./Build/RB_Tree_Sequential3.exe $$node $(seed) $$opt $$valuetosearch; \
				for omp in $(omp_threads); do \
					for mpi in $(mpi_process); do \
						mpiexec -n $$mpi ./Build/RB_Tree_OMP_MPI3.exe $$node $(seed) $$opt $$omp $$valuetosearch -fopenmp; \
					done; \
				done; \
			done; \
		done; \
	done; \
	for node in $(n_nodes2); do \
		for it in $(iterations); do \
			for valuetosearch in $(num_to_search2); do \
				echo "----------------------------{OMP+MPI: $$node,$$opt,$$it,$$valuetosearch}----------------------------"; \
				./Build/RB_Tree_Sequential3.exe $$node $(seed) $$opt $$valuetosearch; \
				for omp in $(omp_threads); do \
					for mpi in $(mpi_process); do \
						mpiexec -n $$mpi ./Build/RB_Tree_OMP_MPI3.exe $$node $(seed) $$opt $$omp $$valuetosearch -fopenmp; \
					done; \
				done; \
			done; \
		done; \
	done; \
	for node in $(n_nodes3); do \
		for it in $(iterations); do \
			for valuetosearch in $(num_to_search3); do \
				echo "----------------------------{OMP+MPI: $$node,$$opt,$$it,$$valuetosearch}----------------------------"; \
				./Build/RB_Tree_Sequential3.exe $$node $(seed) $$opt $$valuetosearch; \
				for omp in $(omp_threads); do \
					for mpi in $(mpi_process); do \
						mpiexec -n $$mpi ./Build/RB_Tree_OMP_MPI3.exe $$node $(seed) $$opt $$omp $$valuetosearch -fopenmp; \
					done; \
				done; \
			done; \
		done; \
	done;

cudacompile:
	gcc -c -o ./Build/RB_Tree_Generator.o ./Source/RB_Tree_Generator.c -O0
	gcc -c -o ./Build/RB_Tree_Sequential.o ./Source/RB_Tree_Sequential.c -O0
	gcc -o ./Build/RB_Tree_Sequential0.exe ./Build/RB_Tree_Sequential.o ./Build/RB_Tree_Generator.o -O0

	nvcc -Xcompiler -fopenmp -o ./Build/RB_Tree_OMP_CUDA0.exe ./Source/RB_Tree_OMP_CUDA.cu -O0
	nvcc -Xcompiler -fopenmp -o ./Build/RB_Tree_OMP_CUDA0L1.exe ./Source/RB_Tree_OMP_CUDA.cu -D L1_CACHE -O0

	nvcc -Xcompiler -fopenmp -o ./Build/RB_Tree_OMP_CUDA1.exe ./Source/RB_Tree_OMP_CUDA.cu -O1
	nvcc -Xcompiler -fopenmp -o ./Build/RB_Tree_OMP_CUDA1L1.exe ./Source/RB_Tree_OMP_CUDA.cu -D L1_CACHE -O1

	nvcc -Xcompiler -fopenmp -o ./Build/RB_Tree_OMP_CUDA2.exe ./Source/RB_Tree_OMP_CUDA.cu -O2
	nvcc -Xcompiler -fopenmp -o ./Build/RB_Tree_OMP_CUDA2L1.exe ./Source/RB_Tree_OMP_CUDA.cu -D L1_CACHE -O2

	nvcc -Xcompiler -fopenmp -o ./Build/RB_Tree_OMP_CUDA3.exe ./Source/RB_Tree_OMP_CUDA.cu -O3
	nvcc -Xcompiler -fopenmp -o ./Build/RB_Tree_OMP_CUDA3L1.exe ./Source/RB_Tree_OMP_CUDA.cu -D L1_CACHE -O3

cudatest0:
	@opt=0; \
	for node in $(n_nodes1); do \
		for it in $(iterations); do \
			for valuetosearch in $(num_to_search1); do \
				echo "----------------------------{CUDA: $$node,$$opt,$$it,$$valuetosearch}----------------------------"; \
				./Build/RB_Tree_Sequential0.exe $$node $(seed) $$opt $$valuetosearch; \
				for omp in $(omp_threads_cuda); do \
					./Build/RB_Tree_OMP_CUDA0.exe $$node $(seed) $$opt $$omp $$valuetosearch; \
					./Build/RB_Tree_OMP_CUDA0L1.exe $$node $(seed) $$opt $$omp $$valuetosearch; \
				done; \
			done; \
		done; \
	done; \
	for node in $(n_nodes2); do \
		for it in $(iterations); do \
			for valuetosearch in $(num_to_search2); do \
				echo "----------------------------{CUDA: $$node,$$opt,$$it,$$valuetosearch}----------------------------"; \
				./Build/RB_Tree_Sequential0.exe $$node $(seed) $$opt $$valuetosearch; \
				for omp in $(omp_threads_cuda); do \
					./Build/RB_Tree_OMP_CUDA0.exe $$node $(seed) $$opt $$omp $$valuetosearch; \
					./Build/RB_Tree_OMP_CUDA0L1.exe $$node $(seed) $$opt $$omp $$valuetosearch; \
				done; \
			done; \
		done; \
	done; \
	for node in $(n_nodes3); do \
		for it in $(iterations); do \
			for valuetosearch in $(num_to_search3); do \
				echo "----------------------------{CUDA: $$node,$$opt,$$it,$$valuetosearch}----------------------------"; \
				./Build/RB_Tree_Sequential0.exe $$node $(seed) $$opt $$valuetosearch; \
				for omp in $(omp_threads_cuda); do \
					./Build/RB_Tree_OMP_CUDA0.exe $$node $(seed) $$opt $$omp $$valuetosearch; \
					./Build/RB_Tree_OMP_CUDA0L1.exe $$node $(seed) $$opt $$omp $$valuetosearch; \
				done; \
			done; \
		done; \
	done;

cudatest1:
	@opt=1; \
	for node in $(n_nodes1); do \
		for it in $(iterations); do \
			for valuetosearch in $(num_to_search1); do \
				echo "----------------------------{CUDA: $$node,$$opt,$$it,$$valuetosearch}----------------------------"; \
				./Build/RB_Tree_Sequential1.exe $$node $(seed) $$opt $$valuetosearch; \
				for omp in $(omp_threads_cuda); do \
					./Build/RB_Tree_OMP_CUDA1.exe $$node $(seed) $$opt $$omp $$valuetosearch; \
					./Build/RB_Tree_OMP_CUDA1L1.exe $$node $(seed) $$opt $$omp $$valuetosearch; \
				done; \
			done; \
		done; \
	done; \
	for node in $(n_nodes2); do \
		for it in $(iterations); do \
			for valuetosearch in $(num_to_search2); do \
				echo "----------------------------{CUDA: $$node,$$opt,$$it,$$valuetosearch}----------------------------"; \
				./Build/RB_Tree_Sequential1.exe $$node $(seed) $$opt $$valuetosearch; \
				for omp in $(omp_threads_cuda); do \
					./Build/RB_Tree_OMP_CUDA1.exe $$node $(seed) $$opt $$omp $$valuetosearch; \
					./Build/RB_Tree_OMP_CUDA1L1.exe $$node $(seed) $$opt $$omp $$valuetosearch; \
				done; \
			done; \
		done; \
	done; \
	for node in $(n_nodes3); do \
		for it in $(iterations); do \
			for valuetosearch in $(num_to_search3); do \
				echo "----------------------------{CUDA: $$node,$$opt,$$it,$$valuetosearch}----------------------------"; \
				./Build/RB_Tree_Sequential1.exe $$node $(seed) $$opt $$valuetosearch; \
				for omp in $(omp_threads_cuda); do \
					./Build/RB_Tree_OMP_CUDA1.exe $$node $(seed) $$opt $$omp $$valuetosearch; \
					./Build/RB_Tree_OMP_CUDA1L1.exe $$node $(seed) $$opt $$omp $$valuetosearch; \
				done; \
			done; \
		done; \
	done;

cudatest2:
	@opt=2; \
	for node in $(n_nodes1); do \
		for it in $(iterations); do \
			for valuetosearch in $(num_to_search1); do \
				echo "----------------------------{CUDA: $$node,$$opt,$$it,$$valuetosearch}----------------------------"; \
				./Build/RB_Tree_Sequential2.exe $$node $(seed) $$opt $$valuetosearch; \
				for omp in $(omp_threads_cuda); do \
					./Build/RB_Tree_OMP_CUDA2.exe $$node $(seed) $$opt $$omp $$valuetosearch; \
					./Build/RB_Tree_OMP_CUDA2L1.exe $$node $(seed) $$opt $$omp $$valuetosearch; \
				done; \
			done; \
		done; \
	done; \
	for node in $(n_nodes2); do \
		for it in $(iterations); do \
			for valuetosearch in $(num_to_search2); do \
				echo "----------------------------{CUDA: $$node,$$opt,$$it,$$valuetosearch}----------------------------"; \
				./Build/RB_Tree_Sequential2.exe $$node $(seed) $$opt $$valuetosearch; \
				for omp in $(omp_threads_cuda); do \
					./Build/RB_Tree_OMP_CUDA2.exe $$node $(seed) $$opt $$omp $$valuetosearch; \
					./Build/RB_Tree_OMP_CUDA2L1.exe $$node $(seed) $$opt $$omp $$valuetosearch; \
				done; \
			done; \
		done; \
	done; \
	for node in $(n_nodes3); do \
		for it in $(iterations); do \
			for valuetosearch in $(num_to_search3); do \
				echo "----------------------------{CUDA: $$node,$$opt,$$it,$$valuetosearch}----------------------------"; \
				./Build/RB_Tree_Sequential2.exe $$node $(seed) $$opt $$valuetosearch; \
				for omp in $(omp_threads_cuda); do \
					./Build/RB_Tree_OMP_CUDA2.exe $$node $(seed) $$opt $$omp $$valuetosearch; \
					./Build/RB_Tree_OMP_CUDA2L1.exe $$node $(seed) $$opt $$omp $$valuetosearch; \
				done; \
			done; \
		done; \
	done;

cudatest3:
	@opt=3; \
	for node in $(n_nodes1); do \
		for it in $(iterations); do \
			for valuetosearch in $(num_to_search1); do \
				echo "----------------------------{CUDA: $$node,$$opt,$$it,$$valuetosearch}----------------------------"; \
				./Build/RB_Tree_Sequential3.exe $$node $(seed) $$opt $$valuetosearch; \
				for omp in $(omp_threads_cuda); do \
					./Build/RB_Tree_OMP_CUDA3.exe $$node $(seed) $$opt $$omp $$valuetosearch; \
					./Build/RB_Tree_OMP_CUDA3L1.exe $$node $(seed) $$opt $$omp $$valuetosearch; \
				done; \
			done; \
		done; \
	done; \
	for node in $(n_nodes2); do \
		for it in $(iterations); do \
			for valuetosearch in $(num_to_search2); do \
				echo "----------------------------{CUDA: $$node,$$opt,$$it,$$valuetosearch}----------------------------"; \
				./Build/RB_Tree_Sequential3.exe $$node $(seed) $$opt $$valuetosearch; \
				for omp in $(omp_threads_cuda); do \
					./Build/RB_Tree_OMP_CUDA3.exe $$node $(seed) $$opt $$omp $$valuetosearch; \
					./Build/RB_Tree_OMP_CUDA3L1.exe $$node $(seed) $$opt $$omp $$valuetosearch; \
				done; \
			done; \
		done; \
	done; \
	for node in $(n_nodes3); do \
		for it in $(iterations); do \
			for valuetosearch in $(num_to_search3); do \
				echo "----------------------------{CUDA: $$node,$$opt,$$it,$$valuetosearch}----------------------------"; \
				./Build/RB_Tree_Sequential3.exe $$node $(seed) $$opt $$valuetosearch; \
				for omp in $(omp_threads_cuda); do \
					./Build/RB_Tree_OMP_CUDA3.exe $$node $(seed) $$opt $$omp $$valuetosearch; \
					./Build/RB_Tree_OMP_CUDA3L1.exe $$node $(seed) $$opt $$omp $$valuetosearch; \
				done; \
			done; \
		done; \
	done;
	
cudastatistics:
	gcc -c -o ./Build/Callingpy.o ./Source/Callingpy.c 
	gcc -o ./Build/Callingpy.exe ./Build/Callingpy.o
	./Build/Callingpy.exe 1

cudaplots:
	gcc -c -o ./Build/Callingpy.o ./Source/Callingpy.c 
	gcc -o ./Build/Callingpy.exe ./Build/Callingpy.o
	./Build/Callingpy.exe 3

cudatables:
	Rscript ./Source/maketablesCUDA.R