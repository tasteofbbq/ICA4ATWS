#include ".\\headfiles\\Nurbs_Basis.hpp"
#include "mex.h"

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    /**
     * Return the 2D NURBS basis functions and first/second derivatives to matlab
     * All non-zero basis functions and derivatives at point [xi,eta] are computed.
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
    plhs[0] = mxCreateDoubleMatrix(6, Nocp, mxREAL);
    
    double* Rders = mxGetPr(plhs[0]);
    
    /* Declare variable */
    int i, j;
    int counter = 0;
    int uspan = FindSpan(uknots, p, mu, u[0]);
    int vspan = FindSpan(vknots, q, mv, u[1]);
    
    double** ddNw = New2DArray(6, Nocp);
    
    Nurbs2DBasisFun2ndDer(elemWgts, uknots, vknots, p, q, uspan, vspan, u, ddNw);
    
    for (i = 0;i < Nocp;i++)
    {
        for (j = 0;j < 6;j++)
        {
            Rders[counter] = ddNw[j][i];
            counter += 1;
        }
    }
    
    /* Free memory */
    Del2DArray(ddNw, 6);
}