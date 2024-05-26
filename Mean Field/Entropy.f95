PROGRAM ENTROPY

   IMPLICIT NONE

   INTEGER(4), PARAMETER :: N = 3, dim = 2
   REAL(8) :: S(N, 2**N), R(N, dim**N)
   REAL(8) :: h, J, T, h_0
   REAL(8) :: m(N), Entropia(dim**N)
   CHARACTER(:), ALLOCATABLE :: filename

   J =  -1.00D0
   h =   0.00D0
   T =   0.10D0
   h_0 = 0.50D0

   CALL Spin(S, N)
   CALL RandomField(R, N, dim, h_0)

   IF(dim .EQ. 2) THEN
      filename = '/home/mateus/Code/Fortran/Mean Field/Data/Entropy/Bi.dat'
      OPEN(UNIT=2, FILE=filename)
   END IF

   IF(dim .EQ. 3) THEN
      filename = '/home/mateus/Code/Fortran/Mean Field/Data/Entropy/Tri.dat'
      OPEN(UNIT=3, FILE=filename)
   END IF

   DO WHILE (h <= 10)

      PRINT '(F6.1, A)', H/10*100, '% '

      CALL SelfConsistency(S, R, N, dim, h, J, T, m)
      CALL getEntropy(S, R, N, dim, h, J, T, m, Entropia)

      IF (dim .EQ. 2) WRITE(2, '(F6.4, 3F8.4)') h, m(1)/N, m(2)/N, m(3)/N

      IF (dim .EQ. 3) WRITE(3, '(F6.4, 3F8.4)') h, m(1)/N, m(2)/N, m(3)/N

      h = h + 0.1

   END DO

   IF (dim .EQ. 2) CLOSE(2)

   IF (dim .EQ. 3) CLOSE(3)

END PROGRAM


SUBROUTINE SelfConsistency(S, R, N, dim, h, J, T, m)

   INTEGER(4) :: N, dim
   REAL(8) :: S(N, 2**N), R(N, dim**N)
   REAL(8) :: h, J, T
   REAL(8) :: prec, error
   REAL(8) :: iter
   REAL(8) :: m(3), m0(3)

   m(1) =  0.00D0
   m(2) = -1.00D0
   m(3) =  1.00D0

   error = 1.00D0
   prec = 1.0e-8
   iter = 0.00D0

   DO WHILE (error >= prec .AND. iter <= 10000)

      CALL getMagnetization(N, dim, S, R, m, H, J, T, m0)

      error = MAXVAL(ABS(m0 - m))

      m = 0.1*m0 + 0.9*m

      !m = m0

      iter = iter + 1.00D0

   END DO

END SUBROUTINE

SUBROUTINE getMagnetization(N, dim, S, R, m, H, J, T, m0)

   INTEGER(4) :: ii, i, N, dim
   REAL(8) :: S(N,2**N), R(N, dim**N)
   REAL(8) :: E(dim**N, 2**N), Z(dim**N)
   REAL(8) :: h, T, J
   REAL(8) :: m(3), X(dim**N, N), m0(3)

   DO ii = 1, dim**N
      DO i = 1, 2**N
         E(ii, i) = - J*(S(1,i)*S(2,i) + S(2,i)*S(3,i) + S(3,i)*S(1,i)) &
            - h*SUM(S(:,i)) - SUM(R(:,ii)*S(:,i)) &
            - J*((2*(m(2)+m(3))*S(1, i)) + (2*(m(1)+m(3))*S(2, i)) + (2*(m(1)+m(2))*S(3, i)))
      END DO
   END DO

   DO ii = 1, dim**N
      Z(ii) = SUM(exp(-E(ii,:)/T))
   END DO

   X = 0.0D0

   DO ii = 1, dim**N
      DO k = 1, N
         X(ii, k) = X(ii, k) + (1D0/Z(ii))*SUM((S(k,:)*exp(-(E(ii,:))/T)))
      END DO
   END DO

   m0(1) = SUM(X(:,1))/(dim**N)
   m0(2) = SUM(X(:,2))/(dim**N)
   m0(3) = SUM(X(:,3))/(dim**N)

END

SUBROUTINE getEntropy(S, R, N, dim, h, J, T, m, Entropia)

   INTEGER(4) :: ii, i, N, dim
   REAL(8) :: S(N, 2**N), R(N, dim**N)
   REAL(8) :: E(dim**N, 2**N), Zh(dim**N), Zl(dim**N)
   REAL(8) :: Entropia(dim**N), Fh(dim**N), Fl(dim**N)
   REAL(8) :: h, T, Tl, Th, J
   REAL(8) :: m(3)

   DO ii = 1, dim**N
      DO i = 1, 2**N
         E(ii, i) = - J*(S(1,i)*S(2,i) + S(2,i)*S(3,i) + S(3,i)*S(1,i)) &
            - h*SUM(S(:,i)) - SUM(R(:,ii)*S(:,i)) &
            - J*((2*(m(2)+m(3))*S(1, i)) + (2*(m(1)+m(3))*S(2, i)) + (2*(m(1)+m(2))*S(3, i)))
      END DO
   END DO

   ! MÉTODO DAS DIFERENÇAS FINITAS CENTRADAS PARA O CÁLCULO DA ENTROPIA

   h = 1.0e-4

   ! CALCULANDO PARA T + h (T superior)
   Th = T + h

   DO ii = 1, dim**N

      Zh(ii) = SUM(exp(-(E(ii, :) - MINVAL(E(ii, :)))/Th))

      Fh(ii) = - Th*LOG(Zh(ii))

   END DO

   ! CALCULANDO PARA T - h (T inferior)
   Tl = T - h

   DO ii = 1, dim**N

      Zl(ii) = SUM(exp(-(E(ii, :) - MINVAL(E(ii, :)))/Tl))

      Fl(ii) = - Tl*LOG(Zl(ii))

   END DO

   DO ii = 1, dim**N

      Entropia(ii) = - (Fh(ii) - Fl(ii)) / (2.0*h)

   END DO

END


SUBROUTINE Spin(S, N)

   IMPLICIT NONE

   REAL(8) :: S(N, 2**N)
   INTEGER(4) :: i, j, N

   S = 1.D0

   DO j = 1, 2**N

      DO i = 1, N

         IF(btest(j-1,i-1)) S(i,j) = -1

      END DO

   END DO

END

SUBROUTINE RandomField(R, N, dim, h_0)

   IMPLICIT NONE

   REAL(8) :: R(N, dim**N), H_0
   INTEGER(4) :: N, dim

   IF(dim .EQ. 2) CALL Bimodal(R, N, dim, H_0)

   IF(dim .EQ. 3) CALL Trimodal(R, N, dim, H_0)

END

SUBROUTINE Bimodal(R, N, dim, h_0)

   IMPLICIT NONE

   REAL(8) :: R(N, dim**N), H_0
   INTEGER(4) :: i, j, N, dim

   IF(dim .EQ. 2) THEN

      R = H_0

      DO j = 1, dim**N
         DO i = 1, N
            IF (btest(j-1,i-1)) R(i,j) = - H_0
         END DO
      END DO
   END IF

END

SUBROUTINE Trimodal(R, N, dim, h_0)

   IMPLICIT NONE

   REAL(8) :: R(N, dim**N), H_0
   INTEGER(4) :: j, N, dim
   INTEGER(4) :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14, i15

   IF(N .EQ. 3) then
      j=1
      DO i1= -1, 1
         DO i2= -1, 1
            DO i3= -1, 1
               R(3,j) = i3*H_0
               R(2,j) = i2*H_0
               R(1,j) = i1*H_0
               j= j+1
            END DO
         END DO
      END DO
   END IF

   IF(N .EQ. 6) then
      j=1
      DO i1= -1, 1
         DO i2= -1, 1
            DO i3= -1, 1
               DO i4 = -1, 1
                  DO i5 = -1, 1
                     DO i6 = -1, 1
                        R(6,j) = i6*H_0
                        R(5,j) = i5*H_0
                        R(4,j) = i4*H_0
                        R(3,j) = i3*H_0
                        R(2,j) = i2*H_0
                        R(1,j) = i1*H_0
                        j= j+1
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
   END IF

   IF(N .EQ. 9) then
      j=1
      DO i1= -1, 1
         DO i2= -1, 1
            DO i3= -1, 1
               DO i4 = -1, 1
                  DO i5 = -1, 1
                     DO i6 = -1, 1
                        DO i7 = -1, 1
                           DO i8 = -1, 1
                              DO i9 = -1, 1
                                 R(9, j) = i9*H_0
                                 R(8, j) = i8*H_0
                                 R(7,j) = i7*H_0
                                 R(6,j) = i6*H_0
                                 R(5,j) = i5*H_0
                                 R(4,j) = i4*H_0
                                 R(3,j) = i3*H_0
                                 R(2,j) = i2*H_0
                                 R(1,j) = i1*H_0
                                 j= j+1
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
   END IF

   IF(N .EQ. 15) then
      j=1
      DO i1= -1, 1
         DO i2= -1, 1
            DO i3= -1, 1
               DO i4 = -1, 1
                  DO i5 = -1, 1
                     DO i6 = -1, 1
                        DO i7 = -1, 1
                           DO i8 = -1, 1
                              DO i9 = -1, 1
                                 R(15, j) = i15*H_0
                                 R(14, j) = i14*H_0
                                 R(13, j) = i13*H_0
                                 R(12, j) = i12*H_0
                                 R(11, j) = i11*H_0
                                 R(10, j) = i10*H_0
                                 R(9, j) = i9*H_0
                                 R(8, j) = i8*H_0
                                 R(7,j) = i7*H_0
                                 R(6,j) = i6*H_0
                                 R(5,j) = i5*H_0
                                 R(4,j) = i4*H_0
                                 R(3,j) = i3*H_0
                                 R(2,j) = i2*H_0
                                 R(1,j) = i1*H_0
                                 j= j+1
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
   END IF

END

