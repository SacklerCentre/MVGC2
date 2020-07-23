C     Computational routine
C
C     mex tooti.c
C     ar r libtooti.a tooti.o

      subroutine tooti(y_output, x_input)
      real*8 x_input, y_output

      y_output = 2.0 * x_input

      return
      end
