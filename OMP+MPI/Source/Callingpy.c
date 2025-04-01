/*
 * Course: High Performance Computing 2023/2024
 *
 * Lecturer: Francesco Moscato	fmoscato@unisa.it
 *
 * Author:
 * Conato Christian		0622702273		c.conato@studenti.unisa.it
 * 
 * Provide a parallell version of the RB Tree search's algorithm.
 * 
 * Final Projects are assigned to each student. Student shall provide 
 * a parallel version of an algorithm with both "OpenMP + MPI" and 
 * "OpenMP + Cuda" approaches, comparing results with a known solution 
 * on single-processing node. 
 * Results and differences shall be discussed for different inputs 
 * (type and size). 
 * The parallel algorithm used in "OpenMP + MPI" solution 
 * could not be the same of the "OpenMP + CUDA" approach.
 * 
 * //This file is a utility files that calls python files
 * 
 * Copyright (C) 2024 - All Rights Reserved
 *
 * This file is part of ProjectHPC.
 *
 * This program is free software: you can redistribute it and/or modify it under the terms of
 * the GNU General Public License as published by the Free Software Foundation, either version
 * 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with ContestOMP.
 * If not, see <https://www.gnu.org/licenses/gpl-3.0.html>.
 */

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[]) {
    int mode = atoi(argv[1]);
    if(mode == 0) //If the modality is zero analize.py will be called
        system("python3 analize.py");
    if(mode == 1) //If the modality is zero analizeCuda.py will be called
        system("python3 analizeCuda.py");
    if(mode == 2) //If the modality is zero makeplotsMPI.py will be called
        system("python3 makeplotsMPI.py");
    if(mode == 3) //If the modality is zero makeplotsCuda.py will be called
        system("python3 makeplotsCuda.py");
    return 0;
}
