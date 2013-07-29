#include <R.h>
#include <Rinternals.h>

#include <string.h>
#include <stdio.h>

SEXP checkColsBinary(SEXP Xobj);
SEXP makeNameOrder2(SEXP Xnames, int* nameLengths, int j, int k);
SEXP makeNameOrder3(SEXP Xnames, int* nameLengths, int j, int k, int m);
SEXP expandToOrder2(SEXP Xobj, SEXP colsBinary, SEXP numBinary,
		    SEXP numNotBinary);
SEXP expandToOrder3(SEXP Xobj, SEXP X2obj, SEXP colsBinary, SEXP numBinary,
		    SEXP numNotBinary);
/* int factorial(int x); */

/* compute x factorial */

/* int factorial(int x) { */

/*   if(x < 0) error("error in factorial function -- negative input"); */

/*   if(x == 0) return 1; */

/*   return x * factorial(x - 1); */
/* } */


/* an integer vector is returned with length = ncol(Xobj)
   each element specifies how many unique values are in that
   column: 1 (1), 2 (2), or 3 or more (3) */

SEXP checkColsBinary(SEXP Xobj) {

  /* extract dimensions */

  SEXP Xdim = getAttrib(Xobj, R_DimSymbol);

  int n = INTEGER(Xdim)[0];
  int p = INTEGER(Xdim)[1];

  /* get a direct pointer to the data */

  double *X;

  X = REAL(Xobj);

  SEXP uniqueVals;

  PROTECT(uniqueVals = allocVector(INTSXP, p));
  
  int *uniVals = INTEGER(uniqueVals);

  double val1, val2;

  for(int j = 0; j < p; j++) {

    val1 = X[j * n];
    int numVals = 1;
    int i = 1;

    while(numVals == 1 && i < n) {

      if(X[i + j * n] != val1) {

	val2 = X[i + j * n];
	numVals = 2;
      }

      i++;
    }
    
    while(numVals == 2 && i < n) {

      double nextVal = X[i + j * n];

      if(nextVal != val1 && nextVal != val2) {

	numVals = 3;
      }

      i++;
    }

    uniVals[j] = numVals;
  }

  UNPROTECT(1);
  return uniqueVals;

}

SEXP makeNameOrder2(SEXP Xnames, int* nameLengths, int j, int k) {

  if(j == k) {

    int length = nameLengths[j] + 6; /* +1 for null char */
    
    char name[length];
    sprintf(name, "I(%s^2)", CHAR(STRING_ELT(Xnames, j)));

    return mkChar(name);
  }

  int length = nameLengths[j] + nameLengths[k] + 2;
  char name[length];
  sprintf(name, "%s:%s", CHAR(STRING_ELT(Xnames, j)),
	  CHAR(STRING_ELT(Xnames, k)));

  return mkChar(name);

}


SEXP makeNameOrder3(SEXP Xnames, int* nameLengths, int j, int k, int m) {

  if(j != k) {

    if(k != m) { /* (1:1:1) */

      int length = nameLengths[j] + nameLengths[k] + nameLengths[m] + 3;
      char name[length];
      sprintf(name, "%s:%s:%s", CHAR(STRING_ELT(Xnames, j)),
	      CHAR(STRING_ELT(Xnames, k)), CHAR(STRING_ELT(Xnames, m)));

      return mkChar(name);
    }

    /* k == m (1:2) */
  
    int length = nameLengths[j] + nameLengths[k] + 7;
    char name[length];
    sprintf(name, "%s:I(%s^2)", CHAR(STRING_ELT(Xnames, j)),
	    CHAR(STRING_ELT(Xnames, k)));

    return mkChar(name);
  }

  /* j == k */

  if(k != m) { /* (2:1) */

    int length = nameLengths[j] + nameLengths[m] + 7;
    char name[length];
    sprintf(name, "I(%s^2):%s", CHAR(STRING_ELT(Xnames, j)),
	    CHAR(STRING_ELT(Xnames, m)));

    return mkChar(name);

  }

  /* final case: all three the same */

  int length = nameLengths[j] + 6;
    
  char name[length];
  sprintf(name, "I(%s^3)", CHAR(STRING_ELT(Xnames, j)));

  return mkChar(name);
}

/* First argument should be a real matrix,
   last three arguments should be integer vectors
   
   colsBinary is partial output of checkColsBinary function */


SEXP expandToOrder2(SEXP Xobj, SEXP colsBinary, SEXP numBinary,
		    SEXP numNotBinary) {

  /* Do error checking in R */

  /* extract dimensions */

  SEXP Xdim = getAttrib(Xobj, R_DimSymbol);

  int n = INTEGER(Xdim)[0];
  int p = INTEGER(Xdim)[1];

  /* get number of binary and non binary columns */
  
  int p2 = INTEGER(numBinary)[0];
  int p3 = INTEGER(numNotBinary)[0];

  if( (p2 + p3) != p )
    error("problem with checking whether columns are binary -- p2 + p3 != p");

  /* calculate number of columns of X2 */

  int X2ncol = p2 * (p2 - 1) / 2 + p3 * (p3 + 1) / 2 + p2 * p3;

  /* allocate it */

  SEXP X2obj;

  PROTECT(X2obj = allocMatrix(REALSXP, n, X2ncol));

  /* Allocate the Column names for the result */

  SEXP X2colnames;

  PROTECT(X2colnames = allocVector(STRSXP, X2ncol));

  /* get pointer to X's colnames which is a STRSXP */

  SEXP Xcolnames = VECTOR_ELT(getAttrib(Xobj, R_DimNamesSymbol), 1);

  /* find out how long each name is */

  int colnameLengths[p];

  for(int i = 0; i < p; i++) {

    colnameLengths[i] = strlen(CHAR(STRING_ELT(Xcolnames, i)));
  }
  
  /* get direct pointers to the data */

  double *X = REAL(Xobj);

  double *X2 = REAL(X2obj);

  int *colsBin = INTEGER(colsBinary);

  /* int which points at a column of X2 */

  int X2colindex = 0;

  /* start looping */

  for(int j = 0; j < p; j++) {

    for(int k = j; k < p; k++) {

      /* don't add a column or advance when j == k and that column is binary */

      if(!(j == k && colsBin[j] == 2)) {

	SET_STRING_ELT(X2colnames, X2colindex,
		       makeNameOrder2(Xcolnames, colnameLengths, j, k));

	for(int i = 0; i < n; i++) {

	  X2[i + X2colindex * n] = X[i + j * n] * X[i + k * n];

	}

	X2colindex++;
      }
    }
  }

  /* add the colnames to X2obj */

  SEXP dimnames;

  PROTECT(dimnames = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(dimnames, 1, X2colnames);
  setAttrib(X2obj, R_DimNamesSymbol, dimnames);

  UNPROTECT(3);
  return X2obj;
}


/* First two arguments should be a real matrices,
   last three arguments should be integer vectors
   
   colsBinary is partial output of checkColsBinary function */

SEXP expandToOrder3(SEXP Xobj, SEXP X2obj,  SEXP colsBinary, SEXP numBinary,
		    SEXP numNotBinary) {

  /* Do error checking in R */

  /* extract dimensions */

  SEXP Xdim = getAttrib(Xobj, R_DimSymbol);

  int n = INTEGER(Xdim)[0];
  int p = INTEGER(Xdim)[1];

  /* get number of binary and non binary columns */
  
  int p2 = INTEGER(numBinary)[0];
  int p3 = INTEGER(numNotBinary)[0];

  /* extract ncol of X2obj */

  SEXP X2dim = getAttrib(X2obj, R_DimSymbol);

  int X2ncol = INTEGER(X2dim)[1];

  /* calculate number of columns of X3 */

  int justBin = p2 * (p2 - 1) * (p2 - 2) / 6;
  int justNonBin = p3 * (p3 + 1) * (p3 + 2) / 6;
  int twoBinOneNonBin = p2 * (p2 - 1) * p3 / 2;
  int oneBinTwoNonBin = p3 * (p3 + 1) * p2 / 2;

  int X3ncol = justBin + justNonBin + twoBinOneNonBin + oneBinTwoNonBin;

  /* allocate it */

  SEXP X3obj;

  PROTECT(X3obj = allocMatrix(REALSXP, n, X3ncol));

  /* Allocate the Column names for the result */

  SEXP X3colnames;

  PROTECT(X3colnames = allocVector(STRSXP, X3ncol));

  /* get pointer to X's colnames which is a STRSXP */

  SEXP Xcolnames = VECTOR_ELT(getAttrib(Xobj, R_DimNamesSymbol), 1);

  /* find out how long each name is */

  int colnameLengths[p];

  for(int i = 0; i < p; i++) {

    colnameLengths[i] = strlen(CHAR(STRING_ELT(Xcolnames, i)));
    /* Rprintf("%d\n", colnameLengths[i]); */
  }

  /* get direct pointers to the data */

  double *X = REAL(Xobj);
  double *X2 = REAL(X2obj);
  double *X3 = REAL(X3obj);

  int *colsBin = INTEGER(colsBinary);

  /* ints which point at the columns of X2 and X3 */

  int X2colindex = 0;
  int X3colindex = 0;

  /* start looping */

  for(int j = 0; j < p; j++) {

    for(int k = j; k < p; k++) {

      if(!(j == k && colsBin[j] == 2)) {

	for(int m = k; m < p; m++) {

	  if(!(k == m && colsBin[k] == 2)
	     && !(j == m && colsBin[j] == 2)) {

	    SET_STRING_ELT(X3colnames, X3colindex,
			   makeNameOrder3(Xcolnames, colnameLengths, j, k, m));

	    for(int i = 0; i < n; i++) {

	      X3[i + X3colindex * n] = X2[i + X2colindex * n] *  X[i + m * n];
	    }

	    X3colindex++;
	  }

	}

	X2colindex++;
      }
    }
  }

  /* add the colnames to X2obj */

  SEXP dimnames;

  PROTECT(dimnames = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(dimnames, 1, X3colnames);
  setAttrib(X3obj, R_DimNamesSymbol, dimnames);

  UNPROTECT(3);
  return X3obj;
}
