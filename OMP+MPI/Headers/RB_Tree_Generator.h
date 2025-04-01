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
 * This file is the header file of RB_Tree_Generator.c.
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

#ifndef RBTREEGENERATOR_H
#define RBTREEGENERATOR_H

// Color of the nodes
#define RED   0
#define BLACK 1

// Structure of a Node in the RB_Tree
struct Node {
    int key;
    int color;
    struct Node* parent;
    struct Node* left;
    struct Node* right;
};

// Structure of the RB_Tree
struct RBTree {
    struct Node* root;
    struct Node* nil;  //Pointer to a leaf node
};

struct Node* createNode(int key, int color);
struct RBTree* initializeRBTree();
void leftRotate(struct RBTree* tree, struct Node* x);
void rightRotate(struct RBTree* tree, struct Node* y);
void RBInsertFixup(struct RBTree* tree, struct Node* z);
void RBInsert(struct RBTree* tree, int key);
void inOrderTraversal(struct Node* node);
int getRandomNumber(int min, int max);
struct RBTree* generateTree(int numElements,int seed,int opt);
void randomElementGenerator(int numElements, struct RBTree* tree);

#endif