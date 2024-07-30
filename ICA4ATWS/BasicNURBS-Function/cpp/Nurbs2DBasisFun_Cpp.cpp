#include ".\\headfiles\\Nurbs_Basis.hpp"
#include "mex.h"

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    /**
     * Return the 2D NURBS basis functions to matlab
     * All non-zero basis functions at point [xi,eta] are computed.
     */
    
    /* Create input variable */
    double* elemWgts = mxGetPr(prhs[0]);
    double* uknots = mxGetPr(prhs[1]);
    double* vknots = mxGetPr(prhs[2]);
    int p = mxGetScalar(prhs[3]);
    int q = mxGetScalar(prhs[4]);
    double* u = mxGetPr(prhs[5]);
    
    int mu = mxGetN(prhs[1]);
    int mv = mxGetN(prhs[2]);
    int Nocp = (p + 1) * (q + 1);
    
    /* Create output variable */
    plhs[0] = mxCreateDoubleMatrix(1, Nocp, mxREAL);
    
    double* R = mxGetPr(plhs[0]);
    
    /* Declare variable */
    int uspan = FindSpan(uknots, p, mu, u[0]);
    int vspan = FindSpan(vknots, q, mv, u[1]);
    
    Nurbs2DBasisFun(elemWgts, uknots, vknots, p, q, uspan, vspan, u, R);
}