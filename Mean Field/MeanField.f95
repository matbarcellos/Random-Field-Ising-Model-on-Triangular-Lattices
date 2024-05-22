!****************************************************************************
!
!  PROGRAMA: Aplicação da técnica do melhoramento de campo médio com clusters p/ o modelo de Ising
!  TAMANHO: REDE DE TRÊS SÍTIOS
!  OBJETIVO: Reproduzir resultados TCC RAMOS (202)
!
!****************************************************************************

PROGRAM MeanField

!***********************************************************************    
! DECLARAÇÃO DE VARIÁVEIS
!***********************************************************************
INTEGER(4), PARAMETER :: N = 3
REAL(8) :: S(N, N**2)
REAL(8) :: T, h, J
REAL(8) :: m(3)

!***********************************************************************


J = -1.0D0
h =  0.0D0
T =  0.1D0

m(1) =    0D0
m(2) =   -1D0
m(3) =    1D0

CALL BaseSpins(S, N)

DO WHILE (h <= 7)

    CALL Self_Consistency(S, N, h, J, T, m)

    WRITE(1,'(F10.2, F10.2)') h, sum(m)/N

    h = h + 0.01D0

END DO

END PROGRAM MeanField


SUBROUTINE Self_Consistency(S, N, h, J, T, m)

!***********************************************************************    
! DECLARAÇÃO DE VARIÁVEIS DA SUBROTINA
!***********************************************************************
    INTEGER(4) :: i, k
    REAL(8) :: S(N, N**2)
    REAL(8) :: E(2**N), Z
    REAL(8) :: h, J, T
    REAL(8) :: prec, error, minE, invZ
    REAL(8) :: iter
    REAL(8) :: m(3), mred(3)
!***********************************************************************
    error = 1.0D0
    prec = 1.0e-6
    iter = 0

    DO WHILE (error >= prec .AND. iter <= 10000)

      DO i = 1, 2**N
         E(i) = - J*(S(1,i)*S(2,i) + S(2,i)*S(3,i) + S(3,i)*S(1,i)) &
         - h*SUM(S(:,i)) - J*((2*(m(2)+m(3))*S(1, i)) + (2*(m(1)+m(3))*S(2, i)) + (2*(m(1)+m(2))*S(3, i)))
      END DO

      minE = MINVAL(E)
      
      Z = SUM(exp(-(E(:))/T))

      invZ = 1D0/Z
        
      !print*, Z, invZ
      !read(*,*)

      mred = 0.0D0

      DO k = 1, N
         mred(k) = mred(k) + invZ*SUM(S(k,:)*exp(-(E(:))/T))
      END DO

        !print*, mred
        !read(*,*)

      error = MAXVAL(ABS(mred - m))

      m = mred
        
      iter = iter + 1

    END DO

END SUBROUTINE


SUBROUTINE BaseSpins(S, N)

   IMPLICIT NONE

   REAL(8) :: S(N, 2**N)
   INTEGER(4) :: i, j, N

   S = 1.D0

   DO j = 1, 2**N

      DO i = 1, N

         if (btest(j-1,i-1)) S(i,j) = -1

      END DO

   END DO

END SUBROUTINE BaseSpins

