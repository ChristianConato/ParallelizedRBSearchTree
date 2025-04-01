/*
 * Course: High Performance Computing 2023/2024
 *
 * Lecturer: Francesco Moscato	fmoscato@unisa.it
 *
 * Author:
 * Conato Christian		0622702273		c.conato@studenti.unisa.it
 *
 * Copyright (C) 2024 - All Rights Reserved
 *
 * //This file is the parallel version with OMP and CUDA of the RB Tree Search to find a random node in the tree
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

#include "../Headers/RB_Tree_Generator.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <omp.h>
#include <time.h>

#define CUDA_CHECK(X) {\
 cudaError_t _m_cudaStat = X;\
 if(cudaSuccess != _m_cudaStat) {\
    fprintf(stderr,"\nCUDA_ERROR: %s in file %s line %d\n",\
    cudaGetErrorString(_m_cudaStat), __FILE__, __LINE__);\
    exit(1);\
 } }

//===========================================================RB TREE GENERATION================================================================//

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

int getRandomNum(int min, int max) {
    return min + rand() % (max - min + 1);
}

/**
 * @brief Create a new node with the given key and color.
 *
 * @param key The key value of the new node.
 * @param color The color of the new node (RED or BLACK).
 * @return A pointer to the newly created node.
 */
struct Node* createNode(int key, int color) {
    struct Node* newNode = (struct Node*)malloc(sizeof(struct Node));
    newNode->key = key;
    newNode->color = color;
    newNode->parent = NULL;
    newNode->left = NULL;
    newNode->right = NULL;
    return newNode;
}

/**
 * @brief Initialize a new Red-Black Tree.
 *
 * @return A pointer to the newly initialized Red-Black Tree.
 */
struct RBTree* initializeRBTree() {
    struct RBTree* tree = (struct RBTree*)malloc(sizeof(struct RBTree));
    tree->nil = createNode(0, BLACK);
    tree->root = tree->nil;
    return tree;
}

/**
 * @brief Perform a left rotation on the Red-Black Tree.
 *
 * @param tree The Red-Black Tree.
 * @param x The node to be rotated.
 */
void leftRotate(struct RBTree* tree, struct Node* x) {
    struct Node* y = x->right;
    x->right = y->left;

    if (y->left != tree->nil) {
        y->left->parent = x;
    }

    y->parent = x->parent;

    if (x->parent == tree->nil) {
        tree->root = y;
    } else if (x == x->parent->left) {
        x->parent->left = y;
    } else {
        x->parent->right = y;
    }

    y->left = x;
    x->parent = y;
}

/**
 * @brief Perform a right rotation on the Red-Black Tree.
 *
 * @param tree The Red-Black Tree.
 * @param y The node to be rotated.
 */
void rightRotate(struct RBTree* tree, struct Node* y) {
    struct Node* x = y->left;
    y->left = x->right;

    if (x->right != tree->nil) {
        x->right->parent = y;
    }

    x->parent = y->parent;

    if (y->parent == tree->nil) {
        tree->root = x;
    } else if (y == y->parent->left) {
        y->parent->left = x;
    } else {
        y->parent->right = x;
    }

    x->right = y;
    y->parent = x;
}

/**
 * @brief Fix the Red-Black Tree properties after node insertion.
 *
 * @param tree The Red-Black Tree.
 * @param z The newly inserted node.
 */
void RBInsertFixup(struct RBTree* tree, struct Node* z) {
        while (z->parent->color == RED) {
        if (z->parent == z->parent->parent->left) {
            struct Node* y = z->parent->parent->right;

            if (y->color == RED) {
                z->parent->color = BLACK;
                y->color = BLACK;
                z->parent->parent->color = RED;
                z = z->parent->parent;
            } else {
                if (z == z->parent->right) {
                    z = z->parent;
                    leftRotate(tree, z);
                }

                z->parent->color = BLACK;
                z->parent->parent->color = RED;
                rightRotate(tree, z->parent->parent);
            }
        } else {
            struct Node* y = z->parent->parent->left;

            if (y->color == RED) {
                z->parent->color = BLACK;
                y->color = BLACK;
                z->parent->parent->color = RED;
                z = z->parent->parent;
            } else {
                if (z == z->parent->left) {
                    z = z->parent;
                    rightRotate(tree, z);
                }

                z->parent->color = BLACK;
                z->parent->parent->color = RED;
                leftRotate(tree, z->parent->parent);
            }
        }
    }

    tree->root->color = BLACK;
}

/**
 * @brief Insert a new key into the Red-Black Tree.
 *
 * @param tree The Red-Black Tree.
 * @param key The key value to be inserted.
 */
void RBInsert(struct RBTree* tree, int key) {
    struct Node* z = createNode(key, RED);
    struct Node* y = tree->nil;
    struct Node* x = tree->root;

    while (x != tree->nil) {
        y = x;
        if (z->key < x->key) {
            x = x->left;
        } else {
            x = x->right;
        }
    }

    z->parent = y;

    if (y == tree->nil) {
        tree->root = z;
    } else if (z->key < y->key) {
        y->left = z;
    } else {
        y->right = z;
    }

    z->left = tree->nil;
    z->right = tree->nil;
    z->color = RED;

    RBInsertFixup(tree, z);
}

/**
 * @brief Generate a random number within the specified range [min, max].
 *
 * @param min The minimum value of the range.
 * @param max The maximum value of the range.
 * @return The generated random number.
 */
int getRandomNumber(int min, int max) {
    return min + rand() % (max - min + 1);
}

/**
 * @brief Generate and insert random elements into the Red-Black Tree.
 * The numbers generated are extracted randomly from the range (1,numElements)
 *
 * @param numElements The number of random elements to generate.
 * @param tree The Red-Black Tree to insert elements into.
 */
void randomElementGenerator(int numElements, struct RBTree* tree){
    for (int i = 0; i < numElements; ++i) {
        int randomNumber = getRandomNumber(1, numElements);
        RBInsert(tree, randomNumber);
    }
}

/**
 * @brief Generate a Red-Black Tree with random elements taken from the
 * srand() function.
 *
 * @param numElements The number of random elements to generate.
 * @param seed The seed for the random number generator.
 * @param opt Unused parameter.
 * @return A pointer to the generated Red-Black Tree.
 */
struct RBTree* generateTree(int numElements, int seed, int opt){
    srand(seed);
    struct RBTree* tree = initializeRBTree();
    randomElementGenerator(numElements,tree);
    //inOrderTraversal(tree->root);
    return tree;
}

//==============================================================OMP+CUDA========================================================================//

/**
 * @brief Prints information about the execution of the program to a CSV file.
 *
 * @param n_nodes Number of nodes in the Red-Black Tree.
 * @param opt Optimization option.
 * @param numThreads Number of OpenMP threads.
 * @param RB_creation_time Time taken for Red-Black Tree creation.
 * @param kernel_execution_time Time taken for the CUDA kernel execution.
 * @param execution_time Total execution time of the program.
 */
void printCUDAToCSV(int n_nodes, int opt, int numThreads, double RB_creation_time, float kernel_execution_time, double execution_time){
    FILE *fp2; 
    char path2[200];
    sprintf(path2, "Informations/OMP_CUDA/opt%d/%d.csv", opt, n_nodes); 
    char *filename2 = path2;
    fp2 = fopen(filename2, "a+");
    if (fp2 == NULL) {
        perror("Errore durante l'apertura del file");
        fprintf(stderr, "Impossibile aprire il file: %s\n", filename2);
    }
    #ifdef L1_CACHE
        fprintf(fp2, "OMP+CUDA_L1;%d;1024;%06f;%06f;%06f;\n", numThreads, RB_creation_time, kernel_execution_time, execution_time);
    #else
        fprintf(fp2, "OMP+CUDA;%d;1024;%06f;%06f;%06f;\n", numThreads, RB_creation_time, kernel_execution_time, execution_time);
    #endif
    fclose(fp2);
}

/**
 * @brief Function useful to print the result of the execution to a TXT file.
 * @param n_nodes The number of nodes in the Red-Black Tree.
 * @param result The result of the RB Search
 * @param valueToSearch key to search in the tree
*/
void printFoundToTXT(int n_nodes, int finalResult, int valueToSearch){
    FILE *fp2; 
    char path2[200];
    sprintf(path2, "RB_Search_Report/nodes%d/%d.txt",n_nodes,n_nodes); 
    char *filename2 = path2;
    fp2 = fopen(filename2, "a+");
    if (fp2 == NULL) {
        perror("Errore durante l'apertura del file");
        fprintf(stderr, "Impossibile aprire il file: %s\n", filename2);
    }
    if(finalResult==0)
        fprintf(fp2, "OMP+CUDA: Key %d NOT found;\n",valueToSearch);
    else
        fprintf(fp2, "OMP+CUDA: Key %d found;\n",valueToSearch);
    fclose(fp2);
}

/**
 * @brief Performs a binary search on a sorted array.
 *
 * @param dsortedArray Pointer to the sorted array in device memory.
 * @param start Starting index for the search.
 * @param end Ending index for the search.
 * @param valueToSearch Value to be searched in the array.
 * @return 1 if the value is found, 0 otherwise.
 */
__host__ __device__ int binarySearch(struct SimpleNode* dsortedArray, int start, int end, int valueToSearch){
    while (start <= end) {
        int mid = start + (end - start) / 2;
        if (dsortedArray[mid].key == valueToSearch) return 1;
        else if (dsortedArray[mid].key < valueToSearch) start = mid + 1;
        else end = mid - 1;
    }
    return 0;
}

/**
 * @brief function to find the minimum among two numbers.
 *
 * @param a first number.
 * @param b second number.
 */
__device__ int cuda_fmin(int a,int b){
    return (a<b) ? a : b;
}

/**
 * @brief CUDA kernel function to perform parallel binary search on the GPU.
 *
 * @param dsortedArray Pointer to the sorted array in device memory.
 * @param arraySizeGPU Size of the array in GPU memory.
 * @param valueToSearch Value to be searched in the array.
 * @param dfound Pointer to the variable storing the result.
 */
__global__ void RBSearchKernel(struct SimpleNode* dsortedArray, int arraySizeGPU,int valueToSearch, int* dfound) {
    //Section for dividing the workload among the threads.
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int elemPerThread = arraySizeGPU/(blockDim.x*gridDim.x);
    int remaining = arraySizeGPU % (blockDim.x*gridDim.x);
    int start = i *elemPerThread + cuda_fmin(i,remaining); //the first i elements will have 1 extra element
    elemPerThread += (i < remaining) ? 1 : 0;
    int end = start + elemPerThread - 1;

    if(start<arraySizeGPU){
        if(binarySearch(dsortedArray,start,end,valueToSearch)){
            atomicAdd(dfound,1); //function that add 1 to dfound if the element has been found
        }
    }
}

/**
 * @brief Searches for a value on the GPU using CUDA and on the CPU using OpenMP, and measures the time taken.
 *
 * @param sortedArray Array of nodes sorted by the Red-Black Tree.
 * @param arraysize Size of the sortedArray.
 * @param numThreadsOMP Number of OpenMP threads to use for CPU search.
 * @param valueToSearch Value to be searched in the array.
 * @param numThreadCUDA Number of CUDA threads per block for GPU search.
 * @return Time taken for the GPU and CPU combined search.
 */
float searchOnDevice(struct Node** sortedArray, int arraysize, int numThreadsOMP, int valueToSearch){
    // Division of the sortedArray into two arrays, the array of the host and the array of the GPU
    int arraySizeGPU=arraysize/2, arraySizeHost=arraysize/2+arraysize%2; //Calculating the sizes of the two arrays
    struct SimpleNode* dsortedArray ; //Array of the GPU
    struct SimpleNode* tempArrayGPU; //TempArray for the GPU copying the 
    struct SimpleNode* hsortedArray; //Array of the host
    tempArrayGPU = (struct SimpleNode*)malloc(arraySizeGPU* sizeof(struct SimpleNode*));
    hsortedArray = (struct SimpleNode*)malloc(arraySizeHost*sizeof(struct SimpleNode*));

    // Filling the host sortedArray
    for (int i = 0; i < arraySizeHost; i++) {
        hsortedArray[i].key = sortedArray[i]->key;
        hsortedArray[i].color = sortedArray[i]->color;
    }

    // Filling the gpu sortedArray
    for (int i = 0; i < arraySizeGPU; i++) {
        tempArrayGPU[i].key = sortedArray[arraySizeHost + i]->key;
        tempArrayGPU[i].color = sortedArray[arraySizeHost + i]->color;
    }

    // Allocation via the cudaMalloc of dsortedArray of the portion of the array for the GPU
    CUDA_CHECK(cudaMalloc((void**)&dsortedArray,arraySizeGPU*sizeof(struct SimpleNode)));
    // Copy of the elements in tempArrayGPU into dsortedArray
    CUDA_CHECK(cudaMemcpy(dsortedArray, tempArrayGPU, arraySizeGPU*sizeof(struct SimpleNode*),cudaMemcpyHostToDevice));
    
    int* dfound;
    CUDA_CHECK(cudaMalloc((void**)&dfound,sizeof(int))); //Allocation via cudaMalloc of dfound that will
                                                        //contain the result of the thread that will find the result 

    /*Useful section for determining the maximum number of blocksize

    int threadsperblock = 1024;
    int blockSizeLimit = 0;
    int dynamicMemSize = 0;
    int blocksize;
    int minGridSize;
    cudaOccupancyMaxPotentialBlockSize(&minGridSize,&blocksize,RBSearchKernel,dynamicMemSize,blockSizeLimit);
    printf("Block size:%d",blocksize);
    int blockspergrid = ((elemPerThread-1)/blocksize+1);*/

    int workload = arraySizeGPU/10; //Useful for dividing the workload among the threads and for 
                                    //lowering the gridSize not to exceed the limit exploiting all threads in a single block
    dim3 blockSize = 1024;
    dim3 gridSize = ((workload-1)/blockSize.x+1);

    cudaEvent_t start, stop;
    CUDA_CHECK(cudaEventCreate(&start));
    CUDA_CHECK(cudaEventCreate(&stop));
    
    #ifdef L1_CACHE
        cudaFuncSetCacheConfig(RBSearchKernel,cudaFuncCachePreferL1); //Function to prioritize the L1 cache to 48KB
    #endif

    CUDA_CHECK(cudaEventRecord(start,0));
    RBSearchKernel<<<gridSize,blockSize>>>(dsortedArray,arraySizeGPU,valueToSearch,dfound);
    CUDA_CHECK(cudaEventRecord(stop,0));

    int numElemPerThreadOMP = arraySizeHost/numThreadsOMP;
    int foundCPU = 0;

    if(arraySizeHost>=numThreadsOMP){
        #pragma omp parallel for reduction(+:foundCPU)
        for (int i = 0; i < numThreadsOMP; i++) {
            int start = i * numElemPerThreadOMP;
            int end = start + numElemPerThreadOMP - 1;
            if (i == numThreadsOMP - 1) 
                end = arraySizeHost - 1;
            if (binarySearch(hsortedArray, start, end, valueToSearch)) 
                foundCPU++;
        }
    } else
        printf("Insufficienti elementi, il valore %d non potrà essere cercato!",valueToSearch);

    CUDA_CHECK(cudaDeviceSynchronize());

    int foundGPU;
    CUDA_CHECK(cudaMemcpy(&foundGPU,dfound,sizeof(int),cudaMemcpyDeviceToHost));
    int finalResult = foundGPU + foundCPU;
   
    float elapsed;
    CUDA_CHECK(cudaEventElapsedTime(&elapsed,start,stop));
    elapsed=elapsed/1000.f;
    
    CUDA_CHECK(cudaEventDestroy(start));
    CUDA_CHECK(cudaEventDestroy(stop));
    CUDA_CHECK(cudaFree(dfound));
    CUDA_CHECK(cudaFree(dsortedArray));
    free(tempArrayGPU);
    free(hsortedArray);
    printFoundToTXT(arraysize, finalResult, valueToSearch);
    return elapsed;
}

/**
 * @brief Main function to execute the Red-Black Tree search program.
 *
 * @param argc Number of command-line arguments.
 * @param argv Array of command-line arguments.
 * @return 0 on successful execution.
 */
int main(int argc, char* argv[]){

    int numElements = atoi(argv[1]); // Nodes of the tree
    int seed = atoi(argv[2]); //Seed for a random number generation
    int opt = atoi(argv[3]); //Optimization taken only for the purpose 
    int numThreadsOMP = atoi(argv[4]); //Threads OMP
    int valueToSearch = atoi(argv[5]); //key to search in the tree
    omp_set_num_threads(numThreadsOMP); //Setting the number of threads

    if(argc != 6){       
		printf("Error: rbsearch argc %d <number of nodes> <seed> <optimization> <Number of thread OMP used> <value to search>\n", argc);
		exit(0);
	}
    
    // Declaration of the variables for the times
    double program_execution, create_time;
    float kernel_execution_time;

    clock_t start_execution = clock(); // Taking the start time of the program execution
    clock_t start_creation = clock(); // Taking the start time of the creation of the tree
    struct RBTree* tree = generateTree(numElements,seed,opt); //Generation of the tree
    clock_t end_creation = clock(); // Taking the end time of the creation of the tree

    // Creation of the sorted array for containing all the nodes of the 
    // RB tree ordered in a ascending order
    int arraysize;
    struct Node** sortedArray = createSortedArray(tree, &arraysize);

    //Calling the searchOnDevice function that returns the kernel execution time
    kernel_execution_time = searchOnDevice(sortedArray,arraysize,numThreadsOMP,valueToSearch);
    free(tree);
    free(sortedArray);
    clock_t end_execution = clock(); // Taking the end time of the program execution

    //Calculation of the program execution time and the creation time
    program_execution = (double)(end_execution-start_execution)/CLOCKS_PER_SEC;
    create_time = (double)(end_creation-start_creation)/CLOCKS_PER_SEC;

    // Call the function for printing the time results
    printCUDAToCSV(numElements, opt, numThreadsOMP, create_time, kernel_execution_time, program_execution);

    return 0;
}