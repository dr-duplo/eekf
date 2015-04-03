/***********************************************************************************
 * The MIT License (MIT)                                                           *
 *                                                                                 *
 * Copyright (c) 2015 Christian Meiﬂner                                            *
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
 * Functions to compute filter states of an extended kalman filter.
 *
 * @copyright	The MIT Licence
 * @file		eekf.h
 * @author 		Christian Meiﬂner
 */

#ifndef EEKF_H
#define EEKF_H

#include <eekf/eekf_mat.h>

/// eekf function return values
typedef enum
{
	eEekfReturnOk = 0,				//!< function call succeded
	eEekfReturnCallbackFailed,		//!< a callback function failed
	eEekfReturnComputationFailed,	//!< a computation failed
	eEekfReturnParameterError,		//!< function parameters are invalid
} eekf_return;

/**
 * Function type to compute the predicted state of the filter and the linearization from current filter
 * state.
 *
 * The function f takes the current filter state and input variables to compute the predicted state
 * xp of the filter and the linearization Jf of f at state x. Jf is the Jacobian of f differentiated
 * with respect to x.
 *
 * @param [out] xp			pointer to the matrix that will hold the state prediction
 * @param [out] Jf			pointer to the matrix that will hold Jacobian of f
 * @param [in]	x			pointer to the matrix hold the current state of the filter
 * @param [in]	u			pointer to the matrix hold the current input variables
 * @param [in]	userData	pointer to the optional user data
 * @return should return eEekfReturnOk if computation succeeded
 */
typedef eekf_return (*ekkf_fun_f)(eekf_mat* xp, eekf_mat* Jf, eekf_mat const *x,
		eekf_mat const *u, void* userData);

/**
 * Function type to compute a measurement prediction and measurement linearization from current
 * filter state.
 *
 * The function h takes the current filter state x to compute the predicted measurement zp and
 * linearization Jh of h at state x. Jh is the Jacobian of h differentiated with respect to x.
 *
 * @param [out]	zp			pointer to the matrix that will hold the measurement prediction
 * @param [out]	Jh			pointer to the matrix that will hold Jacobian of h
 * @param [in]	x			pointer to the matrix hold the current state of the filter
 * @param [in]  userData	pointer to the optional user data
 * @return	should return eEekfReturnOk if computation succeeded
 */
typedef eekf_return (*ekkf_fun_h)(eekf_mat* zp, eekf_mat* Jh, eekf_mat const *x, void* userData);

/// the filter context
typedef struct
{
	eekf_mat *x; 	//!< predicted/corrected state
	eekf_mat *P;	//!< predicted/corrected covariance
	ekkf_fun_f f;	//!< state transition function
	ekkf_fun_h h;	//!< measurement prediction function
	void *userData; //!< pointer to user defined data
} eekf_context;

/**
 * Initialize Extended Kalman Filter context.
 *
 * The function initializes the eekf context for use with predict and correct functions.
 * It holds pointers to the state and covariance matrices, which the user provides to the function.
 * It also includes function pointers for the nonlinear state transition function f and the nonlinear
 * measurement prediction function h.
 *
 * The dimension of x and P must match. DIM(x) = N x 1, DIM(P) = N x N whereas N is the number of
 * states.
 *
 * @param [in/out] ctx 		pointer to the context to initialize
 * @param [in]	   x 		pointer to the matrix holding the current state
 * @param [in]	   P		pointer to the matrix holding the current covariance
 * @param [in]	   f		function pointer to the state transition function
 * @param [in]	   h		function pointer to the measurement prediction function
 * @param [in]	   userData	optional pointer to user data
 * @return returns eEekfReturnOk on success
 */
eekf_return eekf_init(eekf_context *ctx, eekf_mat *x, eekf_mat *P, ekkf_fun_f f,
		ekkf_fun_h h, void *userData);

/**
 * Predict the next filter state.
 *
 * This function predicts the next filter state using the current state, the given input u and
 * the given process covariance. It will call the function f internally to infer the process state.
 * The dimension of noise covariance Q should match the dimensions of context P.
 *
 * @param [in/out] ctx	pointer to the filter context
 * @param [in] 	   u	pointer to the matrix holding input values
 * @param [in]	   Q	pointer to the matrix holding the process covariance
 * @return returns eEekfReturnOk on success
 */
eekf_return eekf_predict(eekf_context *ctx, eekf_mat const *u,
		eekf_mat const *Q);

/**
 * Correct the current filter state.
 *
 * This function corrects the current filter state by using the given real measurement z and
 * the given measurement process noise covariance R.
 * The dimensions of z and R must match: DIM(z) = N x 1, DIM(R) = N x N whereas N is the number of
 * measurement variables.
 *
 * @param [in/out] ctx	pointer to the filter context
 * @param [in]	   z	pointer to the matrix holding the measurement values
 * @param [in]	   R	pointer to the matrix holding the measurement covariance
 * @return
 */
eekf_return eekf_correct(eekf_context *ctx, eekf_mat const *z,
		eekf_mat const *R);

/**
 * Compute a random number of a normal distribution with standard deviation of 1.
 *
 * Remember to initialize the random generator with srand(...) or an equivalent function up front!
 *
 * @return returns the random number
 */
eekf_value eekf_randn();

#endif /* EEKF_H */
