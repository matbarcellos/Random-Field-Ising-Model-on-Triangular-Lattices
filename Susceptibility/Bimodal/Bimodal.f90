PROGRAM CLUSTER_BIMODAL

   USE OMP_LIB
   IMPLICIT NONE

   ! DECLARAÇÃO DE VARIÁVEIS
   INTEGER, PARAMETER :: N = 3
   REAL(8) :: S(N, 2**N), R(N, 2**N)
   INTEGER(4) :: j, ii, ij, il, ik, it
   REAL(8) :: H, H_0
   REAL(8) :: SpinInteraction(2**N), X(2**N)
   CHARACTER(LEN=10) :: filename
   REAL(8), DIMENSION(5) :: T_values = [0.000001D0, 0.1D0, 0.2D0, 0.3D0, 0.5D0]

   ! EXECUÇÃO DO CÓDIGO
   H = 0.D0
   H_0 = 0.5D0

   call BaseSpins(S, N)             ! BASE DE SPINS
   call BaseCampos(R, N, H_0)       ! BASE DE CAMPOS ALEATÓRIOS

   !***********************************************************************
   DO j = 1, 2**N
      DO ii = 1, N

         call rede(N, ii, ij, il, ik) ! ij, il e ik são vizinhos de ii

         !print*, ii, ij, il, ik
         !read(*,*)

         IF(ij>0)SpinInteraction(j) = SpinInteraction(j) + S(ii, j)*S(ij, j)
         IF(il>0)SpinInteraction(j) = SpinInteraction(j) + S(ii, j)*S(il, j)
         IF(ik>0)SpinInteraction(j) = SpinInteraction(j) + S(ii, j)*S(ik, j)

      END DO
   END DO
   !***********************************************************************

   DO it = 1, 5

      WRITE(filename, "('T',I1,'.dat')") it
      OPEN(UNIT=1, FILE=filename)

      DO WHILE (H <= 10)

         PRINT '(F6.1, A)', H/10*100, '% '
         PRINT '(A)', '------------------------------------------------------'

         call Susceptibilidade(S, R, N, H, SpinInteraction, X, T_values(it))

         WRITE(1, *)  H, SUM(X)/(N*(2**N))

         H = H + 0.01D0

      END DO

      ! Reset H para o valor inicial
      H = 0.D0

      ! Fechar o arquivo antes de começar a próxima iteração
      CLOSE(UNIT=1)
   END DO


END PROGRAM CLUSTER_BIMODAL

SUBROUTINE Susceptibilidade(S, R, N, H, SpinInteraction, X, T)

   REAL(8) :: S(N, 2**N), R(N, 2**N), Energia(2**N, 2**N), Z(2**N), M(2**N), M_d(2**N), M_l(2**N), T
   REAL(8) :: X(2**N), SpinInteraction(2**N)
   REAL(8) :: H, H_d, H_l, dh

   dh = 1.0e-4

   H_d = H + dh

   call Energy(S, R, N, H_d, SpinInteraction, Energia)
   call PartitionFunction(N, T, Z, Energia)
   call Magnetization(N, S, T, Z, Energia, M)

   M_d(:) = M(:)

   H_l = H - dh

   call Energy(S, R, N, H_l, SpinInteraction, Energia)
   call PartitionFunction(N, T, Z, Energia)
   call Magnetization(N, S, T, Z, Energia, M)

   M_l(:) = M(:)

   X(:) = (M_d(:) - M_l(:))/(2D0*dh)

END SUBROUTINE Susceptibilidade

SUBROUTINE Energy(S, R, N, H, SpinInteraction, Energia)

   IMPLICIT NONE

   INTEGER(4) :: N
   INTEGER :: i, j
   REAL(8) :: Energia(2**N, 2**N), S(N, 2**N), R(N, 2**N), SpinInteraction(2**N), H

   !$OMP PARALLEL DO
   DO i = 1, 2**N ! LOOP DAS CONFIGURAÇÕES DE CAMPO ALEATÓRIO

      Energia(i, :) = SpinInteraction(:)

      DO j = 1, 2**N  ! LOOP DAS CONFIGURAÇÕES DE SPIN

         Energia(i, j) =  Energia(i, j) - H*SUM(S(:,j)) - SUM(R(:,i)*S(:,j))

      END DO
   END DO
   !$OMP END PARALLEL DO

   PRINT '(A, F6.1)', 'H =', H
   PRINT '(A)', 'ENERGIA'

END SUBROUTINE

SUBROUTINE PartitionFunction(N, T, Z, Energia)

   IMPLICIT NONE

   REAL(8) :: Energia(2**N, 2**N), Z(2**N), T
   INTEGER(4) :: N, i
   REAL(8) :: minEnergy

   !$OMP PARALLEL DO
   DO i = 1, 2**N

      minEnergy = MINVAL(Energia(i, :))
      Z(i) = SUM(exp(-(Energia(i, :) - minEnergy) / T))

   END DO
   !$OMP END PARALLEL DO

   PRINT '(A)', 'FUNCAO DE PARTICAO'

END SUBROUTINE PartitionFunction

SUBROUTINE Magnetization(N, S, T, Z, Energia, M)

   IMPLICIT NONE

   INTEGER(4) :: N, i
   REAL(8) :: Energia(2**N, 2**N), S(N, 2**N), M(2**N), Z(2**N), T
   REAL(8) :: inverse_Z, minEnergy

   !$OMP PARALLEL DO
   DO i = 1, 2**N ! LOOP DAS CONFIGURAÇÕES DE CAMPO ALEATÓRIO

      inverse_Z = 1.D0 / Z(i)
      minEnergy = MINVAL(Energia(i, :))

      M(i) = inverse_Z * DOT_PRODUCT(exp(-(Energia(i, :) - minEnergy) / T), SUM(S, DIM=1))

   END DO
   !$OMP END PARALLEL DO

   PRINT '(A, F6.1)', 'MAGNETIZACAO ', SUM(M)/(N*(2**N))

END SUBROUTINE Magnetization

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

SUBROUTINE BaseCampos(R, N, H_0)

   IMPLICIT NONE

   REAL(8) :: R(N, 2**N), H_0
   INTEGER(4) :: i, j, N

   R = H_0

   DO j = 1, 2**N

      DO i = 1, N

         if (btest(j-1,i-1)) R(i,j) = - H_0

      END DO

   END DO

END SUBROUTINE BaseCampos

subroutine rede(ns, i, j, l, k)
   !entrada: i sitio ; ns=numero de sitios
   !saída: j, k, l são vizinhos de i

   implicit none
   integer ( kind = 4 ) :: i, j, k,l,ns

   j=0
   k=0
   l=0

   !       Cluster
   !             1
   !            / \
   !           3 - 2

   if (ns.eq.3) then
      j=i+1
      if (i.eq.3) j=1

   end if

   !Cluster
   !             1
   !            / \
   !           6 - 2
   !          / \ / \
   !         5 - 4 - 3


   if (ns.eq.6) then
      j=i+1
      if (i.eq.6) j=1
      if (i.eq.2) k=4
      if (i.eq.4) k=6
      if (i.eq.6) k=2
   end if

   !Cluster ns=9
   !             1
   !            / \
   !           8 - 2
   !          / \ / \
   !         7 - 9 - 3
   !          \ / \ /
   !           6 - 4
   !            \ /
   !             5

   if (ns.eq.9) then
      if (i.eq.1) then
         j=2
         k=8
      end if

      if (i.eq.2) then
         j=3
         k=9
      end if

      if (i.eq.3) then
         j=4
         k=9
      end if

      if (i.eq.4) then
         j=6
         k=9
      end if

      if (i.eq.5) then
         j=4
         k=6
      end if
      if (i.eq.6) then
         j=7
         k=9
      end if
      if (i.eq.7) then
         j=8
         k=9
      end if
      if (i.eq.8) then
         j=2
         k=9
      end if
      !           print*, i, j, k
      !           read(*,*)

   end if

   ! Cluster com ns=15
   !               1
   !              / \
   !             9 - 2
   !            / \ / \
   !           8 - 10 -3
   !          / \ / \ / \
   !         7 - 6 - 5 - 4
   !        / \ / \ / \ / \
   !       11- 12- 13 -14 -15

   if (ns.eq.15) then

      if (i.eq.1) then
         j=9
         k=2
      end if

      if (i.eq.2) then
         j=10
         k=3
         l=9
      end if

      if (i.eq.3) then
         j=5
         k=4
         l=10
      end if

      if (i.eq.4) then
         j=14
         k=5
      end if

      if (i.eq.5) then
         j=13
         k=14
         l=6
      end if

      if (i.eq.6) then
         j=13
         k=12
         l=7
      end if

      if (i.eq.7) then
         j=12
         k=11
      end if

      if (i.eq.8) then
         j=7
         k=6
         l=10
      end if

      if (i.eq.9) then
         j=8
         k=10
      end if

      if (i.eq.10) then
         j=6
         k=5
      end if

      if (i.eq.11) then
         j=12
      end if

      if (i.eq.12) then
         j=13
      end if

      if (i.eq.13) then
         j=14
      end if

      if (i.eq.14) then
         j=15
      end if

      if (i.eq.15) then
         j=4
      end if
   end if

END
