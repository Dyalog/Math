/*****************************************************************************
  Copyright (c) 2014, Intel Corp.
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of Intel Corporation nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
  THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************
* Contents: Native middle-level C interface to LAPACK function dgesvdx
* Author: Intel Corporation
* Generated November, 2011
*****************************************************************************/

#include "lapacke_utils.h"

lapack_int LAPACKE_dgesvdx_work( int matrix_layout, char jobu, char jobvt, char range,
                           		lapack_int m, lapack_int n, double* a,
                          		lapack_int lda, lapack_int vl, lapack_int vu,
                           		lapack_int il, lapack_int iu, lapack_int ns,
                           		double* s, double* u, lapack_int ldu,
                           		double* vt, lapack_int ldvt,	
                                double* work, lapack_int lwork, lapack_int* iwork )
{
    lapack_int info = 0;
    if( matrix_layout == LAPACK_COL_MAJOR ) {
        /* Call LAPACK function and adjust info */
        LAPACK_dgesvdx( &jobu, &jobvt,  &range, &m, &n, a, &lda, &vl, &vu,
            			&il, &iu, &ns, s, u, &ldu, vt, &ldvt,
                        work, &lwork, iwork, &info );
        if( info < 0 ) {
            info = info - 1;
        }
    } else if( matrix_layout == LAPACK_ROW_MAJOR ) {
        lapack_int nrows_u = ( LAPACKE_lsame( jobu, 'a' ) ||
                             LAPACKE_lsame( jobu, 's' ) ) ? m : 1;
        lapack_int ncols_u = LAPACKE_lsame( jobu, 'a' ) ? m :
                             ( LAPACKE_lsame( jobu, 's' ) ? MIN(m,n) : 1);
        lapack_int nrows_vt = LAPACKE_lsame( jobvt, 'a' ) ? n :
                              ( LAPACKE_lsame( jobvt, 's' ) ? MIN(m,n) : 1);
        lapack_int lda_t = MAX(1,m);
        lapack_int ldu_t = MAX(1,nrows_u);
        lapack_int ldvt_t = MAX(1,nrows_vt);
        double* a_t = NULL;
        double* u_t = NULL;
        double* vt_t = NULL;
        /* Check leading dimension(s) */
        if( lda < n ) {
            info = -8;
            LAPACKE_xerbla( "LAPACKE_dgesvdx_work", info );
            return info;
        }
        if( ldu < ncols_u ) {
            info = -16;
            LAPACKE_xerbla( "LAPACKE_dgesvdx_work", info );
            return info;
        }
        if( ldvt < n ) {
            info = -18;
            LAPACKE_xerbla( "LAPACKE_dgesvdx_work", info );
            return info;
        }
        /* Query optimal working array(s) size if requested */
        if( lwork == -1 ) {
            LAPACK_dgesvdx( &jobu, &jobvt, &range, &m, &n, a, &lda_t, &vl, &vu,
            				&il, &iu, &ns, s, u, &ldu_t, vt,
                            &ldvt_t, work, &lwork, iwork, &info );
            return (info < 0) ? (info - 1) : info;
        }
        /* Allocate memory for temporary array(s) */
        a_t = (double*)LAPACKE_malloc( sizeof(double) * lda_t * MAX(1,n) );
        if( a_t == NULL ) {
            info = LAPACK_TRANSPOSE_MEMORY_ERROR;
            goto exit_level_0;
        }
        if( LAPACKE_lsame( jobu, 'a' ) || LAPACKE_lsame( jobu, 's' ) ) {
            u_t = (double*)
                LAPACKE_malloc( sizeof(double) * ldu_t * MAX(1,ncols_u) );
            if( u_t == NULL ) {
                info = LAPACK_TRANSPOSE_MEMORY_ERROR;
                goto exit_level_1;
            }
        }
        if( LAPACKE_lsame( jobvt, 'a' ) || LAPACKE_lsame( jobvt, 's' ) ) {
            vt_t = (double*)
                LAPACKE_malloc( sizeof(double) * ldvt_t * MAX(1,n) );
            if( vt_t == NULL ) {
                info = LAPACK_TRANSPOSE_MEMORY_ERROR;
                goto exit_level_2;
            }
        }
        /* Transpose input matrices */
        LAPACKE_dge_trans( matrix_layout, m, n, a, lda, a_t, lda_t );
        /* Call LAPACK function and adjust info */
        LAPACK_dgesvdx( &jobu, &jobvt, &range, &m, &n, a, &lda_t, &vl, &vu,
            				&il, &iu, &ns, s, u, &ldu_t, vt,
                            &ldvt_t, work, &lwork, iwork, &info );
        if( info < 0 ) {
            info = info - 1;
        }
        /* Transpose output matrices */
        LAPACKE_dge_trans( LAPACK_COL_MAJOR, m, n, a_t, lda_t, a, lda );
        if( LAPACKE_lsame( jobu, 'a' ) || LAPACKE_lsame( jobu, 's' ) ) {
            LAPACKE_dge_trans( LAPACK_COL_MAJOR, nrows_u, ncols_u, u_t, ldu_t,
                               u, ldu );
        }
        if( LAPACKE_lsame( jobvt, 'a' ) || LAPACKE_lsame( jobvt, 's' ) ) {
            LAPACKE_dge_trans( LAPACK_COL_MAJOR, nrows_vt, n, vt_t, ldvt_t, vt,
                               ldvt );
        }
        /* Release memory and exit */
        if( LAPACKE_lsame( jobvt, 'a' ) || LAPACKE_lsame( jobvt, 's' ) ) {
            LAPACKE_free( vt_t );
        }
exit_level_2:
        if( LAPACKE_lsame( jobu, 'a' ) || LAPACKE_lsame( jobu, 's' ) ) {
            LAPACKE_free( u_t );
        }
exit_level_1:
        LAPACKE_free( a_t );
exit_level_0:
        if( info == LAPACK_TRANSPOSE_MEMORY_ERROR ) {
            LAPACKE_xerbla( "LAPACKE_dgesvdx_work", info );
        }
    } else {
        info = -1;
        LAPACKE_xerbla( "LAPACKE_dgesvdx_work", info );
    }
    return info;
}
