/***********************************************************************************
 * The MIT License (MIT)                                                           *
 *                                                                                 *
 * Copyright (c) 2015 Christian Meißner                                            *
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
 * @copyright   The MIT Licence
 * @file        eekf_mat.c
 * @author      Christian Meißner
 */

#include <eekf/eekf_mat.h>

#include <stddef.h>
#include <string.h>
#include <math.h>

eekf_mat* eekf_mat_mul(eekf_mat *C, eekf_mat const *A, eekf_mat const *B)
{
    if ( NULL == C || NULL == A || NULL == B || A->cols != B->rows
            || C->rows * C->cols != A->rows * B->cols)
    {
        return NULL;
    }

    uint8_t c;
    uint8_t r;
    uint8_t i;
    eekf_value *res;
    eekf_value *value1;
    eekf_value *value2;

    memset(C->elements, 0, sizeof(eekf_value) * C->rows * C->cols);

    C->rows = A->rows;
    C->cols = B->cols;

    for (c = 0; c < C->cols; c++)
    {
        for (r = 0; r < C->rows; r++)
        {
            res = C->elements + c * C->rows + r;
            value1 = A->elements + r;
            value2 = B->elements + c * B->rows;
            for (i = 0; i < A->cols; i++, value1 += A->rows, value2++)
            {
                *res += *value1 * *value2;
            }
        }
    }

    return C;
}

eekf_mat* eekf_mat_add(eekf_mat *C, eekf_mat const *A, eekf_mat const *B)
{
    if (NULL == C || NULL == A || NULL == B || A->rows != B->rows
            || A->cols != B->cols || C->rows != A->rows || C->cols != A->cols)
    {
        return NULL;
    }

    uint16_t N = A->rows * A->cols;
    eekf_value *value1 = A->elements;
    eekf_value *value2 = B->elements;
    eekf_value *res = C->elements;
    uint16_t i;

    for (i = 0; i < N; i++, value1++, value2++, res++)
    {
        *res = *value1 + *value2;
    }

    return C;
}

eekf_mat* eekf_mat_sub(eekf_mat *C, eekf_mat const *A, eekf_mat const *B)
{
    if (NULL == C || NULL == A || NULL == B || A->rows != B->rows
            || A->cols != B->cols || C->rows != A->rows || C->cols != A->cols)
    {
        return NULL;
    }

    uint16_t N = A->rows * A->cols;
    eekf_value *value1 = A->elements;
    eekf_value *value2 = B->elements;
    eekf_value *res = C->elements;
    uint16_t i;

    for (i = 0; i < N; i++, value1++, value2++, res++)
    {
        *res = *value1 - *value2;
    }

    return C;
}

eekf_mat* eekf_mat_trs(eekf_mat *At, eekf_mat const *A)
{
    if (NULL == At || NULL == A || A->rows * A->cols != At->cols * At->rows)
    {
        return NULL;
    }

    uint8_t r, c;
    eekf_value *res;
    eekf_value *value;

    At->rows = A->cols;
    At->cols = A->rows;

    for (c = 0; c < At->cols; c++)
    {
        res = At->elements + c * At->rows;
        value = A->elements + c;
        for (r = 0; r < At->rows; r++, res++, value += A->rows)
        {
            *res = *value;
        }
    }

    return At;
}

eekf_mat* eekf_mat_chol(eekf_mat *L, eekf_mat const *A)
{
    if (NULL == L || NULL == A || A->rows * A->cols != L->rows * L->cols
            || A->rows != A->cols)
    {
        return NULL;
    }

    uint16_t n, N, r, c;
    uint16_t offset;
    eekf_value *de, *value;

    N = A->cols;

    L->rows = N;
    L->cols = N;

    memset(L->elements, 0, sizeof(eekf_value) * N * N);
    for (n = 0; n < N; n++)
    {
        offset = n * N + n;
        memcpy(L->elements + offset, A->elements + offset,
                sizeof(eekf_value) * (N - n));
    }

    // @see http://www.seas.ucla.edu/~vandenbe/103/lectures/chol.pdf
    for (n = 0; n < N; n++)
    {
        // get diagonal element
        de = L->elements + n * N + n;
        // check element is positive definite
        if (*de <= 0)
        {
            return NULL;
        }
        // square root of diagonal element in place
        *de = EEKF_MAT_SQRT(*de);
        // divide lower column elements by diagonal element
        for (r = 1; r < N - n; r++)
        {
            *(de + r) /= *de;
        }
        // compose right submatrix
        for (c = n + 1; c < N; c++)
        {
            value = L->elements + c * N + c;
            for (r = c; r < N; r++, value++)
            {
                *value -= *(L->elements + n * N + r)
                        * *(L->elements + n * N + c);
            }
        }
    }

    return L;
}

eekf_mat* eekf_mat_fw_sub(eekf_mat *X, eekf_mat const *L, eekf_mat const *B)
{
    if (NULL == X || NULL == L || NULL == B || L->rows != B->rows
            || X->rows * X->cols != B->cols * L->cols)
    {
        return NULL;
    }

    // set result dimensions
    X->rows = L->cols;
    X->cols = B->cols;

    // loop vars
    uint8_t i, j, k;
    eekf_value *b_i, *x_i, *x_j, *diag, *row, *a_ij;

    // loop over cols of x and b
    for (k = 0; k < X->cols; k++)
    {
        // loop over x rows
        for (i = 0, x_i = EEKF_MAT_COL(*X, k), b_i = EEKF_MAT_COL(*B, k), diag =
                L->elements, row = L->elements; i < X->rows;
                i++, diag += L->rows + 1, x_i++, b_i++, row++)
        {
            *x_i = *b_i;
            // substitute up to (excluding) current x
            for (j = 0, x_j = EEKF_MAT_COL(*X, k), a_ij = row; j < i;
                    j++, x_j++, a_ij += L->rows)
            {
                *x_i -= *a_ij * *x_j;
            }
            // divide by diagonal element of L
            *x_i /= *diag;
        }
    }
    // return result
    return X;
}
