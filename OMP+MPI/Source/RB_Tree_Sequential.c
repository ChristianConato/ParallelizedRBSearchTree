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
 * //This file is the sequential version of the RB Tree Search to find a random node in the tree
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

#include "../Headers/RB_Tree_Sequential.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>

/**
 * @brief Search for a key in a Red-Black Tree recursively.
 *
 * @param root The root of the Red-Black Tree.
 * @param key The key to be searched.
 * @return A pointer to the node containing the key if found, otherwise NULL.
 */
struct Node* RBSearch(struct Node* root, int key){
    if (root == NULL || root->key == key) {
        return root;
    }
    if (key < root->key) {
        return RBSearch(root->left, key);
    }
    return RBSearch(root->right, key);
}

/**
 * @brief Print performance information to a CSV file.
 *
 * @param n_nodes The number of nodes in the Red-Black Tree.
 * @param opt An optimization parameter.
 * @param creation_time Time taken for tree creation.
 * @param communicationtime Time taken for communication (not used in this context).
 * @param rbsearch_time Time taken for Red-Black Tree search.
 * @param execution_time Total execution time.
 */
void printSeqToCSV(int n_nodes, int opt, double creation_time, double communicationtime, double rbsearch_time, double execution_time){
    FILE *fp1;
    FILE *fp2; 
    char path1[200];
    char path2[200];
    sprintf(path1, "Informations/OMP_MPI/opt%d/%d.csv", opt, n_nodes);
    sprintf(path2, "Informations/OMP_CUDA/opt%d/%d.csv", opt, n_nodes);
    char *filename1 = path1;
    char *filename2 = path2;
    fp1 = fopen(filename1, "a+");
    if (fp1 == NULL) {
        perror("Errore durante l'apertura del file");
        fprintf(stderr, "Impossibile aprire il file: %s\n", filename2);
    }
    fprintf(fp1, "Sequential;0;0;%06f;%06f;%06f;%06f;\n",creation_time,0.0,rbsearch_time,execution_time);
    fp2 = fopen(filename2, "a+");
    if (fp2 == NULL) {
        perror("Errore durante l'apertura del file");
        fprintf(stderr, "Impossibile aprire il file: %s\n", filename2);
    }
    fprintf(fp2, "Sequential;0;0;%06f;%06f;%06f;\n",creation_time,rbsearch_time,execution_time);
    fclose(fp1);
    fclose(fp2);
}
/**
 * @brief Function useful to print the result of the execution to a TXT file.
 * @param n_nodes The number of nodes in the Red-Black Tree.
 * @param result The result of the RB Search
 * @param valueToSearch key to search in the tree
*/
void printFoundToTXT(int n_nodes, struct Node* result, int valueToSearch){
    FILE *fp2; 
    char path2[200];
    sprintf(path2, "RB_Search_Report/nodes%d/%d.txt",n_nodes,n_nodes); 
    char *filename2 = path2;
    fp2 = fopen(filename2, "a+");
    if (fp2 == NULL) {
        perror("Errore durante l'apertura del file");
        fprintf(stderr, "Impossibile aprire il file: %s\n", filename2);
    }
    if(result==NULL)
        fprintf(fp2, "Sequential: Key %d NOT found;\n",valueToSearch);
    else
        fprintf(fp2, "Sequential: Key %d found;\n",valueToSearch);
    fclose(fp2);
}

/**
 * @brief Main function of the program.
 *
 * @param argc Number of command-line arguments.
 * @param argv Array of command-line argument strings.
 * @return 0 on successful execution.
 */
int main(int argc, char* argv[]){
    // Input parameters taken from the command line
    int n_nodes = atoi(argv[1]); // Nodes of the tree
    int seed = atoi(argv[2]); //Seed for a random number generation
    int opt = atoi(argv[3]); //Optimization taken only for the purpose 
                            //of printing the csv file in the right folder
    int valueToSearch = atoi(argv[4]); //key to search in the tree
    
    if(argc != 5){  
		printf("Error: rbsearch argc %d <number of nodes> <seed> <optimization> <value to search>\n", argc);
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

    // Start of the measuring of the time of the search function
    gettimeofday(&rbsearch_start,NULL);

    // Search of the random node in the tree
    struct Node* result = RBSearch(tree->root, valueToSearch);

    gettimeofday(&rbsearch_stop,NULL); // End of the measuring of the search time on the tree

    // Calculation of the search time on the RB tree in seconds
    searchstart = (unsigned long long)rbsearch_start.tv_sec * 1000000 + rbsearch_start.tv_usec;
	searchstop = (unsigned long long)rbsearch_stop.tv_sec * 1000000 + rbsearch_stop.tv_usec;
    searchtime = (double)(searchstop - searchstart)/(double)1000000;

    printFoundToTXT(n_nodes,result,valueToSearch); // Call to the function for printing the search results
    // Freeing memory allocated for the result and the tree
    free(result);
    free(tree);
    gettimeofday(&execution_stop,NULL); // End of the measuring of the total execution time of the program
    //Calculation of the total time of the program execution
    programstart = (unsigned long long)execution_start.tv_sec * 1000000 + execution_start.tv_usec;
	programstop = (unsigned long long)execution_stop.tv_sec * 1000000 + execution_stop.tv_usec;
    executiontime = (double)(programstop - programstart)/(double)1000000;

    printSeqToCSV(n_nodes,opt,createtime,communicationtime,searchtime,executiontime); // Call to the functions for printing the time results

    return 0;
}