      FUNCTION CLOSAP( P1, P2, V1, V2, DMIN, VERTEX)
C
C  FUNCTION TO FIND THE CLOSEST APPROACH OF TWO LINES
C  (CALLS FUNCTION DOT)
C
C  INPUT: P1,P2     A POINT ANYWHERE ON THE LINE
C         V1,V2     A VECTOR ALONG THE LINE (NOT NECESSARILY
C                                                   NORMALIZED)
C  OUTPUT: DMIN     DISTANCE BETWEEN LINES AT CLOSEST APPROACH
C          VERTEX   A POINT HALFWAY BETWEEN THE LINES AT THE
C                   CLOSEST APPROACH
C
C  FUNCTION VALUE:  NORMALLY TRUE BUT IF LINES ARE PARALLEL
C                   CLOSAP IS FALSE AND DMIN IS -1
C
      IMPLICIT NONE
      INTEGER I
      REAL*4 VERTEX,DMIN
      REAL*4 P1, P2, V1, V2
      REAL*4 R12
      REAL*4 A,B,C,D,E,DET,T1,T2,Q1,Q2
      REAL*4 DOT
      LOGICAL CLOSAP
      DIMENSION VERTEX(3) , P1(3) , P2(3) , V1(3) , V2(3) , R12(3)

c      write(*,*) 'In CLOSAP ',p1(1),p1(2),p1(3),p2(1),p2(2),p2(3),
c     +   v1(1),v1(2),v1(3),v2(1),v2(2),v2(3)

      DO 10 I = 1,3
        R12(I) = P1(I) - P2(I)
   10 CONTINUE
C
      A = DOT( V1, V1)
      B = DOT( V2, V2)
      C = DOT( V1, V2)
      D = DOT( R12, V1)
      E = DOT( R12, V2)
C
      DET = C*C - A*B
      CLOSAP = DET .NE. 0.
      DMIN = -1.
C
      IF (CLOSAP) THEN
C
        T1 = (B*D - C*E) / DET
        T2 = (C*D - A*E) / DET
C
        DO 20 I = 1,3
          Q1 = P1(I) + T1*V1(I)
          Q2 = P2(I) + T2*V2(I)
          VERTEX(I) = (Q1 + Q2) / 2.
          R12(I) = Q1 - Q2
   20   CONTINUE
        DMIN = SQRT( DOT(R12,R12) )
C
      ENDIF
C
c      write(*,*) 'Out CLOSAP ',DMIN,
c     +   VERTEX(1),VERTEX(2),VERTEX(3)
C
      RETURN
      END

      FUNCTION DOT( X, Y)
C
      IMPLICIT NONE
      REAL*4 X , Y , DOT
      DIMENSION X(3) , Y(3)
      DOT = X(1)*Y(1) + X(2)*Y(2) + X(3)*Y(3)
C
      RETURN
      END
