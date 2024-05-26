PROGRAM MEANFIELD

   IMPLICIT NONE

   INTEGER(4), PARAMETER :: N = 3
   REAL(8) :: S(N, 2**N), R(N, 2**N)
   REAL(8) :: T, h, J, h_0
   REAL(8) :: m(2**N, N)

   J =  -1.0D0
   h =   0.00D0
   T =   0.1D0
   h_0 = 0.0D0

   ! Essa é uma subrotina
   CALL SPIN(S, N)
   CALL CA(R, N, h_0)

   DO WHILE (h <= 10.0D0)

      PRINT '(F6.1, A)', H/10*100, '% '

      CALL Self_Consistency(S, R, N, h, J, T, m)

      WRITE(1, '(F6.2, F8.2)') h, SUM(m)/(N*(2**N))

      h = h + 0.1D0

   END DO

END PROGRAM


SUBROUTINE Self_Consistency(S, R, N, h, J, T, m)

   INTEGER(4) :: i, ii, k
   REAL(8) :: S(N, 2**N), R(N, 2**N)
   REAL(8) :: E(2**N, 2**N), Z(2**N)
   REAL(8) :: h, J, T
   REAL(8) :: prec, error
   REAL(8) :: iter
   REAL(8) :: m(2**N, 3), mred(2**N, 3)

   m(:, 1) =   1.0D0
   m(:, 2) =  -1.0D0
   m(:, 3) =   1.0D0

   error = 1.0D0
   prec = 1.0e-8
   iter = 0.0D0

   DO WHILE (error >= prec .AND. iter <= 10000)

      E = 0.0D0

      DO ii = 1, 2**N
         DO i = 1, 2**N
            E(ii, i) = - J*(S(1,i)*S(2,i) + S(2,i)*S(3,i) + S(3,i)*S(1,i)) &
               - h*SUM(S(:,i)) - SUM(R(:,ii)*S(:,i)) &
               - J*((2*(m(ii, 2)+m(ii, 3))*S(1, i)) + (2*(m(ii, 1)+m(ii, 3))*S(2, i)) + (2*(m(ii, 1)+m(ii, 2))*S(3, i)))
         END DO
      END DO

      Z = 0.0D0

      DO ii = 1, 2**N
         Z(ii) = SUM(exp(-E(ii,:)/T))
      END DO

      !print*, H, invZ
      !read(*,*)

      mred = 0.0D0

      DO ii = 1, 2**N
         DO k = 1, N
            mred(ii, k) = mred(ii, k) + (1D0/Z(ii))*SUM(S(k,:)*exp(-(E(ii,:))/T))
         END DO
      END DO

      !print*, mred
      !read(*,*)

      error = MAXVAL(ABS(mred - m))

      !m = 0.5*mred + 0.5*m

      m = mred

      iter = iter + 1.0D0

   END DO

END SUBROUTINE

SUBROUTINE SPIN(S, N)

   IMPLICIT NONE

   REAL(8) :: S(N, 2**N)
   INTEGER(4) :: i, j, N

   S = 1.D0

   DO j = 1, 2**N

      DO i = 1, N

         IF(btest(j-1,i-1)) S(i,j) = -1

      END DO

   END DO

END SUBROUTINE

SUBROUTINE CA(R, N, h_0)

   IMPLICIT NONE

   REAL(8) :: R(N, 2**N), h_0
   INTEGER(4) :: i, j, N

   R = h_0

   DO j = 1, 2**N

      DO i = 1, N

         if (btest(j-1,i-1)) R(i,j) = - h_0

      END DO

   END DO

END SUBROUTINE




