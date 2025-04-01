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
 * //This file is the parallel version of the RB Tree Search to find a random node in the tree
 * using an OMP+MPI approach
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

#include "../Headers/RB_Tree_OMP_MPI.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>

// Struct used for simplify the structure of a RB Tree node
struct SimpleNode{
    int key;
    int color;
};

/**
 * @brief Count the number of nodes in the Red-Black Tree.
 *
 * @param node The current node being checked.
 * @param sentinel The sentinel node representing NULL.
 * @return The number of nodes in the subtree rooted at the given node.
 */
int countNodes(struct Node* node, struct Node*sentinel) {
    if (node == NULL || node == sentinel) {
        return 0;
    }
    return countNodes(node->left,sentinel) + 1 + countNodes(node->right,sentinel);
}

/**
 * @brief Fill an array with the information of the Red-Black Tree in in-order traversal.
 *
 * @param node The current node being processed.
 * @param array The array to be filled.
 * @param index The current index in the array.
 * @param sentinel The sentinel node representing NULL.
 */
void fillArrayInOrder(struct Node* node, struct Node** array, int* index, struct Node* sentinel) {
    if (node != NULL && node != sentinel) {
        fillArrayInOrder(node->left, array, index, sentinel);
        array[*index] = node;
        (*index)++;
        fillArrayInOrder(node->right, array, index, sentinel);
    }
}

/**
 * @brief Create and fill an array with the information of the Red-Black Tree in sorted order.
 *
 * @param tree The Red-Black Tree.
 * @param arraySize Pointer to store the size of the created array.
 * @return The sorted array of nodes.
 */
struct Node** createSortedArray(struct RBTree* tree, int* arraySize) {
    int numNodes = countNodes(tree->root,tree->nil);
    struct Node** sortedArray = (struct Node**)malloc((numNodes+1)* sizeof(struct Node*));

    int index = 0;
    fillArrayInOrder(tree->root, sortedArray, &index, tree->nil);

    *arraySize = numNodes;
    return sortedArray;
}

/**
 * @brief Perform a binary search on a sorted array in parallel using OpenMP.
 *
 * @param localArray The local portion of the sorted array.
 * @param extralength The size of the local portion.
 * @param key The key to search for.
 * @return The result of the binary search.
 */
struct SimpleNode binarySearch(struct SimpleNode* localArray, int extralength, int key) {
    int end = extralength-1;
    struct SimpleNode result;
    result.key = -1;
    #pragma omp parallel shared(result)
    {
        int tid = omp_get_thread_num();
        int local_start = tid*(extralength/omp_get_num_threads());
        int local_end = (tid == omp_get_num_threads() -1) ? end : local_start + (extralength/omp_get_num_threads()) -1;
        while (local_start <= local_end) {
            int mid = local_start + (local_end - local_start) / 2;
            if (localArray[mid].key == key) {
                result = localArray[mid];
            }

            if (localArray[mid].key < key) {
                local_start = mid + 1;
            } else {
                local_end = mid - 1;
            }
        }
    }
    return result; // Chiave non trovata
}

/**
 * @brief MPI reduction operation to find the maximum key among SimpleNode elements.
 *
 * @param in Input data for a process.
 * @param inout Combined data for a process.
 * @param len Length of the data.
 * @param datatype MPI datatype.
 */
void maxKeyReducer(void *in, void *inout, int *len, MPI_Datatype *datatype){
    struct SimpleNode *inNode = (struct SimpleNode*)in;
    struct SimpleNode *inoutNode = (struct SimpleNode*)inout;

    if(inNode->key > inoutNode->key){
        *inoutNode = *inNode;
    }
}

/**
 * @brief Print performance information to a CSV file for parallel execution.
 *
 * @param n_nodes The number of nodes in the Red-Black Tree.
 * @param opt An optimization parameter.
 * @param numThreads The number of OpenMP threads used.
 * @param size The number of MPI processes.
 * @param creation_time Time taken for tree creation.
 * @param communicationtime Time taken for communication (MPI).
 * @param rbsearch_time Time taken for parallel binary search.
 * @param execution_time Total execution time.
 */
void printParToCSV(int n_nodes, int opt, int numThreads, int size, double creation_time, double communicationtime, double rbsearch_time, double execution_time){
    FILE *fp2; 
    char path2[200];
    sprintf(path2, "Informations/OMP_MPI/opt%d/%d.csv", opt, n_nodes); 
    char *filename2 = path2;
    fp2 = fopen(filename2, "a+");
    if (fp2 == NULL) {
        perror("Errore durante l'apertura del file");
        fprintf(stderr, "Impossibile aprire il file: %s\n", filename2);
    }
    fprintf(fp2, "OMP+MPI;%d;%d;%06f;%06f;%06f;%06f;\n", numThreads, size, creation_time, communicationtime, rbsearch_time, execution_time);
    fclose(fp2);
}

/**
 * @brief Function useful to print the result of the execution to a TXT file.
 * @param n_nodes The number of nodes in the Red-Black Tree.
 * @param result The result of the RB Search
 * @param valueToSearch key to search in the tree
*/
void printFoundToTXT(int n_nodes, struct SimpleNode globalResult, int valueToSearch){
    FILE *fp2; 
    char path2[200];
    sprintf(path2, "RB_Search_Report/nodes%d/%d.txt",n_nodes,n_nodes); 
    char *filename2 = path2;
    fp2 = fopen(filename2, "a+");
    if (fp2 == NULL) {
        perror("Errore durante l'apertura del file");
        fprintf(stderr, "Impossibile aprire il file: %s\n", filename2);
    }
    if(globalResult.key!=valueToSearch)
        fprintf(fp2, "OMP+MPI: Key %d NOT found;\n",valueToSearch);
    else
        fprintf(fp2, "OMP+MPI: Key %d found;\n",valueToSearch);
    fclose(fp2);
}

/**
 * @brief Generate a random number within a specified range.
 *
 * @param min The minimum value of the range.
 * @param max The maximum value of the range.
 * @return The randomly generated number.
 */
int getRandomNum(int min, int max) {
    return min + rand() % (max - min + 1);
}

/**
 * @brief Main function of the parallel program.
 *
 * @param argc Number of command-line arguments.
 * @param argv Array of command-line argument strings.
 * @return 0 on successful execution.
 */
int main(int argc, char* argv[]){
    // MPI Initialization
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int n_nodes = atoi(argv[1]); // Nodes of the tree
    int seed = atoi(argv[2]); //Seed for a random number generation
    int opt = atoi(argv[3]); //Optimization taken only for the purpose 
                            //of printing the csv file in the right folder
    int numThreads = atoi(argv[4]); //Threads OMP
    int valueToSearch = atoi(argv[5]); //key to search in the tree
    omp_set_num_threads(numThreads); //Setting the number of threads

    if(argc != 7){  
		printf("Error: rbsearch argc %d <process mpi> <number of nodes> <seed> <optimization> <Number of thread OMP used> <value to search>\n", argc);
		exit(0);
	}

    // Declaration of the variables for the times
    struct timeval execution_start, execution_stop, tree_creation_start, tree_creation_stop, rbsearch_start, rbsearch_stop;
    unsigned long long createstart, createstop, programstart, programstop, searchstart, searchstop;
    double executiontime, createtime, searchtime, communicationtime;

    // Start of the measuring of the total time of the program execution
    gettimeofday(&execution_start,NULL);

    // Start of the measuring of the creation time of the tree
    gettimeofday(&tree_creation_start,NULL);
    struct RBTree* tree = generateTree(n_nodes,seed,opt);
    gettimeofday(&tree_creation_stop,NULL); // End of the measuring of the creation time of the tree

    // Calculation of the creation time of the RB tree in seconds
    createstart = (unsigned long long)tree_creation_start.tv_sec * 1000000 + tree_creation_start.tv_usec;
	createstop = (unsigned long long)tree_creation_stop.tv_sec * 1000000 + tree_creation_stop.tv_usec;
    createtime = (double)(createstop - createstart)/(double)1000000;

    // Creation of the sorted array for containing all the nodes of the 
    // RB tree ordered in a ascending order
    int arraySize;
    struct Node** sortedArray = createSortedArray(tree, &arraySize);

    // Calculate local distribution of array elements for all the MPI process
    int local_n, start, end;
    if(size<=arraySize) {
        local_n = arraySize / size;
        start = rank * local_n;   
        end = start + local_n - 1;
        if (rank == size - 1) {
            end += arraySize % size; // Add the extra element to the final process
        }
    }
    else {
        fprintf(stderr,"ERRORE: IL NUMERO DI RANKS E'MAGGIORE DELLA DIMENSIONE DELL'ARRAY, SUDDIVISIONE DEL CARICO NON POSSIBILE!");
        return EXIT_FAILURE;
    }

    // Arrays useful for the Scatterv function
    int *displs = malloc((size+1)*sizeof(int)); // Calculate the displacement array for MPI_Scatterv
    int* counts = malloc((size+1)*sizeof(int)); // Calculate the counts array for MPI_Scatterv
    for (int i = 0; i < size; i++) {
        counts[i] = local_n;  
        displs[i] = i*local_n;
        if(i == size - 1)
            counts[i] += arraySize%size;
    }

    //Creation of a tempArray for containing the SimpleNodes of the RB tree
    //useful for the Scatterv function
    struct SimpleNode* tempArray = malloc((arraySize+1) * sizeof(struct SimpleNode));

    //Only the rank 0 will fill tempArray
    if (rank == 0){
        for (int i = 0; i < arraySize; i++) {
            tempArray[i].key = sortedArray[i]->key;
            tempArray[i].color = sortedArray[i]->color;
        }
    }

    /**
     * @brief Define a new MPI datatype for struct SimpleNode.
     *
     * Defines a new MPI datatype MPI_SimpleNode to facilitate communication
     * of the custom structure SimpleNode in MPI operations. It specifies the structure of the
     * datatype using the MPI_Type_create_struct function, which allows the description of complex
     * datatypes by specifying the lengths, offsets, and types of individual elements in the structure.
     * In this case, it describes a datatype composed of two integer fields: key and color.
     * 
     * @param MPI_SimpleNode The MPI datatype representing struct SimpleNode.
     */

    MPI_Datatype MPI_SimpleNode;
    int blocklengths[2] = {1, 1};
    MPI_Aint offsets[2] = {offsetof(struct SimpleNode, key), offsetof(struct SimpleNode, color)};
    MPI_Datatype types[2] = {MPI_INT, MPI_INT};
    MPI_Type_create_struct(2, blocklengths, offsets, types, &MPI_SimpleNode);
    MPI_Type_commit(&MPI_SimpleNode);
    
    //Calculation of the extrasize for the localArray for all the MPI process
    //since the tempArray could contain more elements if the division of the
    //whole structure is not perfect (the last MPI process will obtain more elements)
    long long int extrasize = local_n*((arraySize*2)%size*3);
    long long int extraelement = extrasize - counts[rank];
    if(extraelement>0)
        extrasize -= extraelement;
    else {
        extraelement = extraelement * -1;
        extrasize += extraelement;
    }

    struct SimpleNode* localArray = malloc((extrasize+1)*sizeof(struct SimpleNode)); //Creation of the array for every MPI process
    //Calculation the communication time
    double temptime;
    temptime = MPI_Wtime();
    int mpi_resultscatter = MPI_Scatterv(tempArray,counts,displs,MPI_SimpleNode,localArray,extrasize,MPI_SimpleNode,0,MPI_COMM_WORLD);
    communicationtime += MPI_Wtime() - temptime;

    if (mpi_resultscatter != MPI_SUCCESS) {
        fprintf(stderr, "Errore in MPI_Scatterv. Codice errore MPI: %d\n", mpi_resultscatter);
        MPI_Abort(MPI_COMM_WORLD, mpi_resultscatter);
    }

    // Start of the measuring of the time of the search function
    gettimeofday(&rbsearch_start,NULL); 

    // Search of the random node in the array with the binarySearch
    struct SimpleNode localResult = binarySearch(localArray, extrasize, valueToSearch);

    // MPI Reduction operation to find the maximum key
    MPI_Op maxKeyOp;
    MPI_Op_create(maxKeyReducer,1,&maxKeyOp);

    struct SimpleNode globalResult; //Variable that will contain the globalresult
    //Calculation the communication time
    temptime = MPI_Wtime();
    int mpi_resultreduce = MPI_Reduce(&localResult, &globalResult, 1, MPI_SimpleNode, maxKeyOp, 0, MPI_COMM_WORLD);

    if (mpi_resultreduce != MPI_SUCCESS) {
        fprintf(stderr, "Errore in MPI_Scatterv. Codice errore MPI: %d\n", mpi_resultreduce);
        MPI_Abort(MPI_COMM_WORLD, mpi_resultreduce);
    }

    //Final calculation of the communication time
    communicationtime += MPI_Wtime() - temptime;

    gettimeofday(&rbsearch_stop,NULL); // End of the measuring of the search time on the tree

    // Calculation of the search time on the RB tree in seconds
    searchstart = (unsigned long long)rbsearch_start.tv_sec * 1000000 + rbsearch_start.tv_usec;
    searchstop = (unsigned long long)rbsearch_stop.tv_sec * 1000000 + rbsearch_stop.tv_usec;
    searchtime = (double)(searchstop - searchstart)/(double)1000000;

    // Freeing memory allocated for the result and the tree
    free(sortedArray);
    free(tempArray);
    free(displs);
    free(counts);
    free(localArray);
    MPI_Type_free(&MPI_SimpleNode);

    gettimeofday(&execution_stop,NULL); // End of the measuring of the total execution time of the program
    
    //Calculation of the total time of the program execution
    programstart = (unsigned long long)execution_start.tv_sec * 1000000 + execution_start.tv_usec;
    programstop = (unsigned long long)execution_stop.tv_sec * 1000000 + execution_stop.tv_usec;
    executiontime = (double)(programstop - programstart)/(double)1000000;

    // Call from the rank 0 to the function for printing the results
    if(rank==0){
        printParToCSV(n_nodes,opt,numThreads,size,createtime,communicationtime,searchtime,executiontime);
        printFoundToTXT(n_nodes,globalResult,valueToSearch);
    }
        
    MPI_Finalize();
    return 0;
}