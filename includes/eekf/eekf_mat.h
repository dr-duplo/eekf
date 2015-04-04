/***********************************************************************************
 * The MIT License (MIT)                                                           *
 *                                                                                 *
 * Copyright (c) 2015 Christian Mei�ner                                            *
 *                                                                                 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy    *
 * of this software and associated documentation files (the "Software"), to deal   *
 * in the Software without restriction, including without limitation the rights    *
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell       *
 * copies of the Software, and to permit persons to whom the Software is           *
 * furnished to do so, subject to the following conditions:                        *
 *                                                                                 *
 * The above copyright notice and this permission notice shall be included in all  *
 * copies or substantial portions of the Software.                                 *
 *                                                                                 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR      *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,        *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE     *
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER          *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,   *
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE   *
 * SOFTWARE.                                                                       *
 ***********************************************************************************/

/**
 * Basic matrix handling and computation functions.
 *
 * @copyright	The MIT Licence
 * @file		eekf_mat.h
 * @author 		Christian Mei�ner
 */

#ifndef EEKF_MAT_H
#define EEKF_MAT_H

#include <stdint.h>

#define EEKF_MAT_LOG  log
#define EEKF_MAT_SQRT sqrt
#define EEKF_MAT_RAND rand
#define EEKF_MAT_RAND_MAX RAND_MAX

/// matrix value type
typedef double eekf_value;

/// base matrix structure type
typedef struct {
	eekf_value *elements;	//!< pointer to matrix elements (column major order)
	uint8_t rows;			//!< number of rows
	uint8_t cols;			//!< number of columns
} eekf_mat;

/// assign data to a matrix
#define EEKF_ASSIGN_MATRIX(matrix, data, rows, cols)\
	do {\
		(matrix).elements = data;\
		(matrix).rows = rows;\
		(matrix).cols = cols;\
	} while(0)

/// declare a matrix with a non constant size (e.g. function parameter dependent)
#define EEKF_DECL_MAT_DYN(name, rows, cols)\
	eekf_value name##_elements[(rows)*(cols)];\
	eekf_mat name = {name##_elements, (rows), (cols)};

/// declare a matrix with given size and value initialization
#define EEKF_DECL_MAT_INIT(name, rows, cols, ...)\
	eekf_value name##_elements[(rows)*(cols)] = {__VA_ARGS__};\
	eekf_mat name = {name##_elements, (rows), (cols)};

/// declare a matrix with given size and initialize elements to zero
#define EEKF_DECL_MAT(name, rows, cols) EEKF_DECL_MAT_INIT(name, rows, cols, 0)

/// get pointer to a given row of a matrix
#define EEKF_MAT_ROW(mat, i) ((mat).elements + (i))

/// get pointer to a given col of a matrix
#define EEKF_MAT_COL(mat, j) ((mat).elements + (mat).rows * (j))

/// get pointer to a given element of a matrix
#define EEKF_MAT_EL(mat, r, c) ((mat).elements + (c) * (mat).rows + (r))

/**
 * Multiply two matrices such that C = A * B.
 *
 * @param [out] C 	pointer to matrix to hold the result
 * @param [in]  A 	pointer to left matrix of multiplication
 * @param [in]  B  	pointer to right matrix of multiplication
 * @return returns the pointer to result matrix on success, NULL otherwise
 */
eekf_mat* eekf_mat_mul(eekf_mat *C, eekf_mat const *A, eekf_mat const *B);

/**
 * Adds two matrices  such that C = A + B.
 *
 * @param [out] C  	pointer to matrix to hold the result
 * @param [in]  A  	pointer to left matrix of addition
 * @param [in]  B  	pointer to right matrix of addition
 * @return returns the pointer to result matrix on success, NULL otherwise
 */
eekf_mat* eekf_mat_add(eekf_mat *C, eekf_mat const *A, eekf_mat const *B);

/**
 * Subtracts two matrices such that C = A - B.
 *
 * @param [out] C	pointer to matrix to hold the result
 * @param [in]  A   pointer to left matrix of subtraction
 * @param [in]  B	pointer to right matrix of subtraction
 * @return returns the pointer to result matrix on success, NULL otherwise
 */
eekf_mat* eekf_mat_sub(eekf_mat *C, eekf_mat const *A, eekf_mat const *B);

/**
 * Transposes a matrix  such that  At = A'.
 *
 * @param [out] At  pointer to matrix to hold the result
 * @param [in]  A  	pointer to the matrix to be transposed
 * @return returns the pointer to result matrix on success, NULL otherwise
 */
eekf_mat* eekf_mat_trs(eekf_mat *At, eekf_mat const *A);

/**
 * Computes the Cholesky Factorization for a matrices.
 *
 * Computes Cholesky Factorization such that A = L * L' whereat L is a lower
 * triangular and A a symmetric positive-definite matrix;
 *
 * @param [out] L pointer to matrix to hold the result
 * @param [in]  A pointer to matrix to be factorized
 * @return returns the pointer to result matrix on success, NULL otherwise
 */
eekf_mat* eekf_mat_chol(eekf_mat *L, eekf_mat const *A);

/**
 * Computes the Forward Substitution of a linear equation system.
 *
 * Computes the Forward Substitution of a linear equation system L * X = B such that
 * L is a lower triangular matrix and x and b are matrices with equal number of columns.
 *
 * @param [out] X pointer to matrix to hold the result
 * @param [in]  L pointer to lower triangular matrix
 * @param [in]  B pointer to right equation side matrix
 * @return returns the pointer to result matrix on success, NULL otherwise
 */
eekf_mat* eekf_mat_fw_sub(eekf_mat *X, eekf_mat const *L, eekf_mat const *B);

#endif /* EEKF_MAT_H */
