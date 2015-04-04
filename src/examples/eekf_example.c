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
 * Example program that uses the eekf.
 *
 * @copyright   The MIT Licence
 * @file        eekf_example.c
 * @author      Christian Meißner
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <eekf/eekf.h>

/// the state prediction function: linear case for simplicity
eekf_return transition(eekf_mat* xp, eekf_mat* Jf, eekf_mat const *x,
        eekf_mat const *u, void* userData)
{
    // the Jacobian of transition() at x
    *EEKF_MAT_EL(*Jf, 0, 0) = 1;
    *EEKF_MAT_EL(*Jf, 1, 0) = 0;
    *EEKF_MAT_EL(*Jf, 0, 1) = 1; // dt
    *EEKF_MAT_EL(*Jf, 1, 1) = 1;

    // predict state from current state
    if (NULL == eekf_mat_mul(xp, Jf, x))
    {
        return eEekfReturnComputationFailed;
    }

    return eEekfReturnOk;
}

/// the measurement prediction function
eekf_return measurement(eekf_mat* zp, eekf_mat* Jh, eekf_mat const *x,
        void* userData)
{
    // the Jacobian of measurement() at x
    *EEKF_MAT_EL(*Jh, 0, 0) = 0;
    *EEKF_MAT_EL(*Jh, 0, 1) = 1;

    // compute the measurement from state x
    *EEKF_MAT_EL(*zp, 0, 0) = *EEKF_MAT_EL(*x, 1, 0);

    return eEekfReturnOk;
}

int main(int argc, char **argv)
{

    // filter context
    eekf_context ctx;
    // state of the filter
    EEKF_DECL_MAT_INIT(x, 2, 1, 0);
    EEKF_DECL_MAT_INIT(P, 2, 2, 1, 0, 0, 1);
    // input and process noise variables
    EEKF_DECL_MAT_INIT(u, 1, 1, 0);
    EEKF_DECL_MAT_INIT(Q, 2, 2, 0.33 / 3., 0.33 / 2., 0.33 / 2., 0.33);
    // measurement and measurement noise variables
    EEKF_DECL_MAT_INIT(z, 1, 1, 0);
    EEKF_DECL_MAT_INIT(R, 1, 1, 0.1);

    // initialize the filter context
    eekf_init(&ctx, &x, &P, transition, measurement, NULL);

    /* now we are ready to use the filter */

    // initialize random number generator
    srand(0);

    // print out header
    printf("k x dx P11 P12 P21 P22 rx rdx z\n");
    // loop over time and present some measurements
    int k;
    eekf_value dz = 0, cz = 0, dd = 1;
    for (k = 0; k < 1000; k++)
    {
        // compute measurement
        dd = (k % 20) ? dd : -dd;
        dz = dd * ((k % 20) - 9.5) / 10.0;
        cz += dz;
        *EEKF_MAT_EL(z, 0, 0) = dz + eekf_randn() * 0.1;

        // correct the current filter state
        eekf_correct(&ctx, &z, &R);

        // print out
        printf("%d %f %f %f %f %f %f %f %f %f\n", k, *EEKF_MAT_EL(*ctx.x, 0, 0),
                *EEKF_MAT_EL(*ctx.x, 1, 0), *EEKF_MAT_EL(*ctx.P, 0, 0),
                *EEKF_MAT_EL(*ctx.P, 0, 1), *EEKF_MAT_EL(*ctx.P, 1, 0),
                *EEKF_MAT_EL(*ctx.P, 1, 1), cz, dz, *EEKF_MAT_EL(z, 0, 0));

        // predict the next filter state
        eekf_predict(&ctx, &u, &Q);
    }

    return 0;
}
