// ludcmp = LU decomposition (with back substitution for inverse of a matrix)

// ECMAScript module

import * as hlao from 'matrix-computations';
//import * as hlao from '../matrix-computations/hlao.mjs';

//REF: Numerical recipes in C, page 46
function ludcmp(a,n){
    // in:
    // a
    // n
    
    // out:
    // a
    // indx
    // d 
    
    // Given a matrix a[1..n][1..n], this routine replaces it by the LU decomposition of a rowwise
    // permutation of itself. 'a' and 'n' are input. 'a' is output, arranged as in equation (2.3.14) above;
    // indx[1..n] is an output vector that records the row permutation effected by the partial
    // pivoting; d is output as ±1 depending on whether the number of row interchanges was even
    // or odd, respectively. This routine is used in combination with lubksb to solve linear equations
    // or invert a matrix.
    
    var indx = [];
    var vv=hlao.zeros_vector((n+1),'row'); // vv stores the implicit scaling of each row.
    var d=1.0; // No row interchanges yet.
    for(var i=1;i<=n;i++){ // Loop over rows to get the implicit scaling information.
        var big=0.0; 
        for(var j=1;j<=n;j++){
            var temp=Math.abs(a[i][j]);
            if(temp > big) big=temp;
        }
        //if(big == 0.0) alert("Singular matrix in routine ludcmp");
        if(big == 0.0) console.log("Singular matrix in routine ludcmp");
        // No nonzero largest element.
        vv[i]=1.0/big; // Save the scaling.
    }
    for(var j=1;j<=n;j++){ // This is the loop over columns of Crout’s method.
        for(var i=1;i<j;i++){ // This is equation (2.3.12) except for i = j.
            var sum=a[i][j];
            for(var k=1;k<i;k++) sum -= a[i][k]*a[k][j];
            a[i][j]=sum;
        }
        big=0.0; // Initialize for the search for largest pivot element.
        for(var i=j;i<=n;i++){ // This is i = j of equation (2.3.12) and i = j+1. . .N
            sum=a[i][j]; // of equation (2.3.13).
            for(var k=1;k<j;k++)
                sum -= a[i][k]*a[k][j];
            a[i][j]=sum;
            var dum=vv[i]*Math.abs(sum);
            if(dum >= big) {
                // Is the figure of merit for the pivot better than the best so far?
                big=dum;
                var imax=i;
            }
        }
        if(j != imax){ // Do we need to interchange rows?
            for(var k=1;k<=n;k++){ // Yes, do so...
                dum=a[imax][k];
                a[imax][k]=a[j][k];
                a[j][k]=dum;
            }
            d = -(d); // ...and change the parity of d.
            vv[imax]=vv[j]; // Also interchange the scale factor.
        }
        indx[j]=imax;
        if(a[j][j] == 0.0) a[j][j]=TINY;
        // If the pivot element is zero the matrix is singular (at least to the precision of the
        // algorithm). For some applications on singular matrices, it is desirable to substitute
        // TINY for zero.
        if (j != n){ // Now, finally, divide by the pivot element.
            dum=1.0/(a[j][j]);
            for(var i=j+1;i<=n;i++) a[i][j] *= dum;
        }
    } // Go back for the next column in the reduction.
    //free_vector(vv,1,n);
    
    return [a,indx,d];
}

//REF: Numerical recipes in C, page 47
function lubksb(a,n,indx,b){
    // in:
    // a
    // n
    // indx
    // b
    
    // out:
    // b
    
    // Solves the set of n linear equations A·X = B. Here a[1..n][1..n] is input, not as the matrix
    // A but rather as its LU decomposition, determined by the routine ludcmp. indx[1..n] is input
    // as the permutation vector returned by ludcmp. b[1..n] is input as the right-hand side vector
    // B, and returns with the solution vector X. a, n, and indx are not modified by this routine
    // and can be left in place for successive calls with different right-hand sides b. This routine takes
    // into account the possibility that b will begin with many zero elements, so it is efficient for use
    // in matrix inversion.

    // When 'ii' is set to a positive value, it will become the
    // index of the first nonvanishing element of 'b'. We now
    // do the forward substitution, equation (2.3.6). The
    // only new wrinkle is to unscramble the permutation
    // as we go.
    
    var ii=0;
    for(var i=1;i<=n;i++){
        var ip=indx[i];
        var sum=b[ip];
        b[ip]=b[i];
        if(ii)
            for(var j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
        else if(sum) ii=i; // A nonzero element was encountered, so from now on we
        b[i]=sum;          // will have to do the sums in the loop above.
    }
    for(var i=n;i>=1;i--){ // Now we do the backsubstitution, equation (2.3.7).
        sum=b[i];
        for(j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
        b[i]=sum/a[i][i]; // Store a component of the solution vector X.
    } //All done!
    
    return b;
}

/*
requires:
minimization/globals.js
decomposition/lubksb.js
decomposition/ludcmp.js
decomposition/matrixInverseLU.js
*/

//inverse of a matrix
// - a[1][1]... a[N][N] (hence 'a' is a square matrix)

//Example 1
//var a = [];
//a.push([,,,]);
//a.push([,1,0,2]);
//a.push([,-1,5,0]);
//a.push([,0,3,-9]);
//var N = 3;
//matrixInverseLU(a,N);

//Example 2
//var a = [];
//a.push([,,,]);
//a.push([,0,1,-1]);
//a.push([,1,-1,0]);
//a.push([,1,1,2]);
//var N = 3;

//MATLAB example 1
//X = [1 0 2; -1 5 0; 0 3 -9]
//Y = inv(X)
//Y =
//    0.8824   -0.1176    0.1961
//    0.1765    0.1765    0.0392
//    0.0588    0.0588   -0.0980

//MATLAB example 2
//X = [1 0 2; -1 5 0; 0 3 -9]
//Y = inv(X)
//Y =
//    0.5000    0.7500    0.2500
//    0.5000   -0.2500    0.2500
//   -0.5000   -0.2500    0.2500

//REF: Numerical recipes in C, page 48
//Inverse of a Matrix
function matrixInverseLU(a,N){
    // in:
    // N
    // a - square matrix
    
    // out:
    // y
    
    var y = hlao.zeros_matrix((N+1),(N+1));
    var col = hlao.zeros_vector((N+1),'row');
    var LUdecomp = ludcmp(a,N); // ludcmp returns '[a,indx,d]' - Decompose the matrix just once.
    a = LUdecomp[0];
    var indx = LUdecomp[1];
    var d = LUdecomp[2];   // Set by ludcmp() either +/-1
    for(var j=1;j<=N;j++){ // Find inverse by columns.
        for(var i=1;i<=N;i++) col[i]=0.0;
        col[j]=1.0;
        col = lubksb(a,N,indx,col);
        for(var i=1;i<=N;i++) y[i][j]=col[i];
    }
    
    return y;
}

//MUST CHANGE --> NDIM to suit problem
var NDIM = 19; //dimension
var NTAB = NDIM; //Sets maximum size of tableau.

//monitors
var nfunc = 0;
var ndfunc = 0;

var EPS = 3.0e-8; // Machine precision.
var STPMX = (40000000*EPS); // Scaled maximum step length allowed in line searches.
var ITMAX = 100; // Maximum allowed number of iterations.
var TOLX = (4*EPS); // Convergence criterion on x values.
var ALF = 1.0e-4;     // Ensures sufficient decrease in function value.

var BIG = 1.0e30;
var CON = 1.4;        // Stepsize is decreased by CON at each iteration.
var CON2 = (CON*CON);
var SAFE = 2.0;       //Return when error is SAFE worse than the best so far.
var GTOL = 1.0e-4;

var SWAP = 1; //see function funcFIXED()

var TINY = 1.0e-20; //A small number.

/*
#define FREEALL free_vector(xi,1,n);free_vector(pnew,1,n); \
    free_matrix(hessin,1,n,1,n);free_vector(hdg,1,n);free_vector(g,1,n); \
    free_vector(dg,1,n);
*/

//#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))
function FMAX(a,b){
    var maxarg1 = a;
    var maxarg2 = b;
    if(maxarg1 > maxarg2) return maxarg1;
    else return maxarg2;
}

//#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
function SQR(a){
    var sqrarg = a;
    if(sqrarg == 0.0) return 0.0;
    else return (sqrarg * sqrarg);
}

//#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
function SIGN(a,b){
    if(b >= 0.0) return Math.abs(a);
    else return -1.0*Math.abs(a);
}

//#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ? (iminarg1) : (iminarg2))
function IMIN(a,b){
    var iminarg1 = a;
    var iminarg2 = b;
    if(iminarg1 < iminarg2) return iminarg1;
    else return iminarg2;
}

export {
    ludcmp,
    lubksb,
    matrixInverseLU
};
