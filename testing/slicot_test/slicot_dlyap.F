#include "fintrf.h"
C======================================================================
#if 0
C
C     timestwo.F
C     .F file needs to be preprocessed to generate .for equivalent
C
C     mex -R2018a slicot_dlyap.F -lslicot
C
#endif
C
C     timestwo.f
C
C     Computational function that takes a scalar and doubles it.

C     This is a MEX-file for MATLAB.
C     Copyright 1984-2018 The MathWorks, Inc.
C
C======================================================================

C     Gateway routine
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)

C     Declarations
C      implicit none

      external  SB03MD

C     mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs

C     Function declarations:

      mwPointer mxGetDoubles
      mwPointer mxCreateDoubleMatrix
      mwPointer mxGetM, mxGetN

      character*1  DICO, JOB, FACT, TRANA

      real*8  SCAL, SEP, FER

      integer IWORK(1), INFO

      integer  N, LDWORK
      real*8  A, C, U, WR, WI, DWORK, X
C-----------------------------------------------------------------------

      DICO  = 'D'
      JOB   = 'X'
      FACT  = 'N'
      TRANA = 'N'

      SCAL = 1.0
      SEP  = 0.0
      FER  = 0.0

      N = mxGetM(prhs(1))

      A = mxGetDoubles(prhs(1))
      C = mxGetDoubles(prhs(2))

      U = mxCreateDoubleMatrix(N,N,0)

      WR = mxCreateDoubleMatrix(N,N,0)
      WI = mxCreateDoubleMatrix(N,N,0)

      IWORK(1) = 0

      LDWORK = 2*N*N + 3*N
      DWORK = mxCreateDoubleMatrix(LDWORK,1,0)

      call SB03MD(DICO,JOB,FACT,TRANA,N,A,N,U,N,C,N,SCAL,SEP,
     +            FER,WR,WI,IWORK,DWORK,LDWORK,INFO)

      plhs(1) = mxCreateDoubleMatrix(N,N,0)

      X = mxGetDoubles(plhs(1))

      call mxCopyReal8ToPtr(C,X,N*N)

      call mxDestroyArray(DWORK)
      call mxDestroyArray(WI)
      call mxDestroyArray(WR)
      call mxDestroyArray(U)

      return
      end

C-----------------------------------------------------------------------
