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
 * This file is the header file of RB_Tree_Sequential.c.
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

#include "RB_Tree_Generator.h"

#ifndef RBTREESEQUENTIAL_H
#define RBTREESEQUENTIAL_H

void printSeqToCSV(int n_nodes, int opt, double creation_time, double communicationtime, double rbsearch_time, double execution_time);
void printFoundToTXT(int n_nodes, struct Node* result, int valueToSearch);
struct Node* RBSearch(struct Node* root, int key);

#endif