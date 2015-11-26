#include "api.h"

static void dft_common(int *rank, const int *n, fftw_complex *inout, int sign)
{
	const unsigned flags = FFTW_ESTIMATE;
	X(plan) plan = X(plan_dft)(*rank, n, inout, inout, sign, flags);
	X(execute)(plan);
	X(destroy_plan)(plan);
}

FFTW_EXTERN void dft(int *rank, const int *n, fftw_complex *inout)
{
	dft_common(rank, n, inout, FFTW_FORWARD);
}

FFTW_EXTERN void idft(int *rank, const int *n, fftw_complex *inout)
{
	dft_common(rank, n, inout, FFTW_BACKWARD);
}
