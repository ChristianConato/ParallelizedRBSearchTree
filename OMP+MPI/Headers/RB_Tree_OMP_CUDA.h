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
 * This file is the header file of RB_Tree_OMP_CUDA.c.
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

#include <omp.h>
#include "RB_Tree_Generator.h"

#ifndef RBTREEOMPCUDA_H
#define RBTREEOMPCUDA_H

struct Node* createNode(int key, int color);
struct RBTree* initializeRBTree();
void leftRotate(struct RBTree* tree, struct Node* x);
void rightRotate(struct RBTree* tree, struct Node* y);
void RBInsertFixup(struct RBTree* tree, struct Node* z);
void RBInsert(struct RBTree* tree, int key);
int getRandomNumber(int min, int max);
struct RBTree* generateTree(int numElements,int seed,int opt);
void randomElementGenerator(int numElements, struct RBTree* tree);

int countNodes(struct Node* node, struct Node*sentinel);
void fillArrayInOrder(struct Node* node, struct Node** array, int* index, struct Node* sentinel);
struct Node** createSortedArray(struct RBTree* tree, int* arraySize);
struct SimpleNode binarySearch(struct SimpleNode* localArray, int giga, int key);

void printCUDAToCSV(int n_nodes, int opt, int numThreads, double RB_creation_time, float kernel_execution_time, double execution_time);
void printFoundToTXT(int n_nodes, int finalResult, int valueToSearch);
float searchOnDevice(struct Node** sortedArray, int arraysize, int numThreadsOMP, int valueToSearch);

#endif