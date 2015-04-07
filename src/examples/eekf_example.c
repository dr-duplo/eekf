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

// constant acceleration
eekf_value a = 0.1;
// time step duration
eekf_value dT = 0.1;
// process noise standard deviation
eekf_value s_w = 0.2;
// measurement noise standard deviation
eekf_value s_z = 10;
	
/// the state prediction function: linear case for simplicity
eekf_return transition(eekf_mat* xp, eekf_mat* Jf, eekf_mat const *x,
        eekf_mat const *u, void* userData)
{
	EEKF_DECL_MAT_INIT(xu, 2, 1, 0);
	EEKF_DECL_MAT_INIT(B, 2, 1, dT * dT / 2, dT);
	
    // the Jacobian of transition() at x
    *EEKF_MAT_EL(*Jf, 0, 0) = 1;
    *EEKF_MAT_EL(*Jf, 1, 0) = 0;
    *EEKF_MAT_EL(*Jf, 0, 1) = dT;
    *EEKF_MAT_EL(*Jf, 1, 1) = 1;

	*EEKF_MAT_EL(B, 0, 0) = dT * dT / 2;
	*EEKF_MAT_EL(B, 1, 0) = dT;
	
    // predict state from current state
    if (NULL == eekf_mat_add(xp, eekf_mat_mul(xp, Jf, x), eekf_mat_mul(&xu, &B, u)))
    {
		printf("arg\n");
        return eEekfReturnComputationFailed;
    }
		
    return eEekfReturnOk;
}

/// the measurement prediction function
eekf_return measurement(eekf_mat* zp, eekf_mat* Jh, eekf_mat const *x,
        void* userData)
{
    // the Jacobian of measurement() at x
    *EEKF_MAT_EL(*Jh, 0, 0) = 1;
    *EEKF_MAT_EL(*Jh, 0, 1) = 0;

    // compute the measurement from state x
    *EEKF_MAT_EL(*zp, 0, 0) = *EEKF_MAT_EL(*x, 0, 0);

    return eEekfReturnOk;
}

int main(int argc, char **argv)
{
    // filter context
    eekf_context ctx;
    // state of the filter
    EEKF_DECL_MAT_INIT(x, 2, 1, 0);
    EEKF_DECL_MAT_INIT(P, 2, 2, 
		pow(s_w, 2) * pow(dT, 4) / 4, pow(s_w, 2) * pow(dT, 3) / 2,
		pow(s_w, 2) * pow(dT, 3) / 2, pow(s_w, 2) * pow(dT, 2));
    // input and process noise variables
    EEKF_DECL_MAT_INIT(u, 1, 1, a);
    EEKF_DECL_MAT_INIT(Q, 2, 2, 
		pow(s_w, 2) * pow(dT, 4) / 4, pow(s_w, 2) * pow(dT, 3) / 2,
		pow(s_w, 2) * pow(dT, 3) / 2, pow(s_w, 2) * pow(dT, 2));
    // measurement and measurement noise variables
    EEKF_DECL_MAT_INIT(z, 1, 1, 0);
    EEKF_DECL_MAT_INIT(R, 1, 1, s_z * s_z);

    // initialize the filter context
    eekf_init(&ctx, &x, &P, transition, measurement, NULL);

    /* now we are ready to use the filter */

    // initialize random number generator
    srand(0);

    // print out header
    printf("k x dx P11 P12 P21 P22 rx rdx z\n");
    // loop over time and present some measurements
    int k;
    eekf_value v = 0, p = 0;
    for (k = 0; k < 1000; k++)
    {
        // compute virtual measurement
        p = *EEKF_MAT_EL(u, 0, 0) / 2.0 * pow(k * dT,2.0);
        v = *EEKF_MAT_EL(u, 0, 0) * k * dT;
		*EEKF_MAT_EL(z, 0, 0) = p + eekf_randn() * s_z;

        // correct the current filter state
        eekf_correct(&ctx, &z, &R);

        // print out
        printf("%d %f %f %f %f %f %f %f %f %f\n", k, *EEKF_MAT_EL(*ctx.x, 0, 0),
                *EEKF_MAT_EL(*ctx.x, 1, 0), *EEKF_MAT_EL(*ctx.P, 0, 0),
                *EEKF_MAT_EL(*ctx.P, 0, 1), *EEKF_MAT_EL(*ctx.P, 1, 0),
                *EEKF_MAT_EL(*ctx.P, 1, 1), p, v, *EEKF_MAT_EL(z, 0, 0));

		*EEKF_MAT_EL(u, 0, 0) = a;				
        // predict the next filter state
        eekf_predict(&ctx, &u, &Q);
    }

    return 0;
}
