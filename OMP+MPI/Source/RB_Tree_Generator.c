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
 * This file contains the RB Tree Generation 
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

#include "../Headers/RB_Tree_Generator.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>

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
 * @brief Perform an in-order traversal of the Red-Black Tree and print the nodes.
 *
 * @param node The starting node for the traversal.
 */
void inOrderTraversal(struct Node* node) {
    if (node != NULL) {
        if (node->left != NULL) {
            inOrderTraversal(node->left);
        }

        if (node->color != BLACK) {
            printf("(%d, RED) ", node->key);
        } else {
            printf("(%d, BLACK) ", node->key);
        }

        if (node->right != NULL) {
            inOrderTraversal(node->right);
        }
    }
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