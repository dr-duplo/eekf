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
 * Functions to compute filter states of an extended kalman filter.
 *
 * @copyright   The MIT Licence
 * @file        eekf.c
 * @author      Christian Meißner
 */

#include <eekf/eekf.h>

#include <stdlib.h>
#include <string.h>
#include <math.h>

eekf_return eekf_init(eekf_context *ctx, eekf_mat *x, eekf_mat *P, ekkf_fun_f f,
        ekkf_fun_h h, void *userData)
{
    if (NULL == ctx || NULL == x || NULL == P || NULL == f || NULL == h
            || x->rows != P->rows || x->rows != P->cols)
    {
        return eEekfReturnParameterError;
    }

    // state
    ctx->x = x;
    ctx->P = P;

    // callbacks
    ctx->f = f;
    ctx->h = h;

    // user defined data
    ctx->userData = userData;

    return eEekfReturnOk;
}

eekf_return eekf_predict(eekf_context *ctx, eekf_mat const *u,
        eekf_mat const *Q)
{
    if (NULL == Q || NULL == u || NULL == ctx)
    {
        return eEekfReturnParameterError;
    }

    EEKF_DECL_MAT_DYN(Jf, ctx->x->rows, ctx->x->rows);
    EEKF_DECL_MAT_DYN(Jft, Jf.cols, Jf.rows);
    EEKF_DECL_MAT_DYN(JfP, Jf.rows, Jf.cols);
    EEKF_DECL_MAT_DYN(xp, ctx->x->rows, ctx->x->cols);

    // predict state and linearize system: x1 = f(x,u), Jf = df(x,u)/dx
    if (NULL != ctx->f
            && eEekfReturnOk != ctx->f(&xp, &Jf, ctx->x, u, ctx->userData))
    {
        return eEekfReturnCallbackFailed;
    }
    // copy prediction to state
    memcpy(ctx->x->elements, xp.elements,
            sizeof(eekf_value) * ctx->x->rows * ctx->x->cols);

    // predict covariance Pp = A*P*A' + Q
    if (NULL
            == eekf_mat_add(ctx->P,
                    eekf_mat_mul(ctx->P, eekf_mat_mul(&JfP, &Jf, ctx->P),
                            eekf_mat_trs(&Jft, &Jf)), Q))
    {
        return eEekfReturnComputationFailed;
    }

    return eEekfReturnOk;
}

eekf_return eekf_correct(eekf_context *ctx, eekf_mat const *z,
        eekf_mat const *R)
{
    if (NULL == R || NULL == z || NULL == ctx || z->rows != R->rows
            || z->rows != R->cols)
    {
        return eEekfReturnParameterError;
    }

    // predicted measurement
    EEKF_DECL_MAT_DYN(zp, z->rows, z->cols);
    // measurement linearization
    EEKF_DECL_MAT_DYN(Jh, z->rows, ctx->x->rows);
    // helper matrices
    EEKF_DECL_MAT_DYN(PJht, ctx->x->rows, z->rows);
    EEKF_DECL_MAT_DYN(L, z->rows, z->rows);
    EEKF_DECL_MAT_DYN(U, z->rows, ctx->x->rows);

    // predict measurement and linearize measurement: zp = h(x), Jh = dh(x)/dx
    if (NULL != ctx->h
            && eEekfReturnOk != ctx->h(&zp, &Jh, ctx->x, ctx->userData))
    {
        return eEekfReturnCallbackFailed;
    }

    {
        EEKF_DECL_MAT_DYN(Ct, Jh.cols, Jh.rows);
        // cross covariance
        if (NULL == eekf_mat_mul(&PJht, ctx->P, eekf_mat_trs(&Ct, &Jh)))
        {
            return eEekfReturnComputationFailed;
        }
    }

    // compute cholesky factorization L of innovation covariance S = (Jh*P*Jh' + R) = L*L'
    // for efficient inversion - assumes S is symmetric positive-definite.
    {
        EEKF_DECL_MAT_DYN(S, R->rows, R->cols);
        // cholesky factorization
        if (NULL == eekf_mat_chol(&L, eekf_mat_add( // innovation covariance
                &S, eekf_mat_mul(&S, &Jh, &PJht), R)))
        {
            return eEekfReturnComputationFailed;
        }
    }

    // compute intermediate matrix for computational efficiency
    // K = U / L -> U = (L \ PJh')'
    {
        EEKF_DECL_MAT_DYN(PCtt, PJht.cols, PJht.rows);
        EEKF_DECL_MAT_DYN(LPCtt, z->rows, ctx->x->rows);
        if (NULL
                == eekf_mat_trs(&U,
                        eekf_mat_fw_sub(&LPCtt, &L,
                                eekf_mat_trs(&PCtt, &PJht))))
        {
            return eEekfReturnComputationFailed;
        }
    }

    // correct state
    // x = xp + U * L \ (z - zp)
    {
        EEKF_DECL_MAT_DYN(dz, z->rows, z->cols);
        EEKF_DECL_MAT_DYN(Ldz, z->rows, ctx->x->cols);
        EEKF_DECL_MAT_DYN(cx, ctx->x->rows, ctx->x->cols);

        if (NULL
                == eekf_mat_add(ctx->x, ctx->x,
                        eekf_mat_mul(&cx, &U,
                                eekf_mat_fw_sub(&Ldz, &L,
                                        eekf_mat_sub(&dz, z, &zp)))))
        {
            return eEekfReturnComputationFailed;
        }
    }

    // correct covariance
    // P = Pp - U * U'
    {
        EEKF_DECL_MAT_DYN(Ut, U.cols, U.rows);
        EEKF_DECL_MAT_DYN(UUt, ctx->P->rows, ctx->P->cols);

        if (NULL
                == eekf_mat_sub(ctx->P, ctx->P,
                        eekf_mat_mul(&UUt, &U, eekf_mat_trs(&Ut, &U))))
        {
            return eEekfReturnComputationFailed;
        }
    }

    return eEekfReturnOk;
}

eekf_value eekf_randn()
{
    eekf_value x1, x2, w;
    do
    {
        x1 = 2.0 * (eekf_value) EEKF_MAT_RAND() / EEKF_MAT_RAND_MAX - 1.0;
        x2 = 2.0 * (eekf_value) EEKF_MAT_RAND() / EEKF_MAT_RAND_MAX - 1.0;
        w = x1 * x1 + x2 * x2;
    } while (w >= 1.0);

    return x1 * EEKF_MAT_SQRT((-2.0 * EEKF_MAT_LOG(w)) / w);
}
