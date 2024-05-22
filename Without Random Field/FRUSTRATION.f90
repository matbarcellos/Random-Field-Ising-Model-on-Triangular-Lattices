!  N15_Bimodal.f90
!
!  FUNCTIONS:
!  N15_Bimodal
!

!****************************************************************************
!
!  PROGRAM: N15_Bimodal
!
!****************************************************************************

PROGRAM CLUSTER_FRUSTRATION

   USE OMP_LIB
   IMPLICIT NONE

   ! DECLARAÇÃO DE VARIÁVEIS
   INTEGER, PARAMETER :: N = 15
   INTEGER(4) :: j, ii, ij, il, ik
   REAL(8), ALLOCATABLE :: S(:,:), Energia(:), SpinInteraction(:)
   REAL(8) :: H, T
   REAL(8) :: M, Z

   ALLOCATE(S(N, 2**N), Energia(2**N), SpinInteraction(2**N))

   ! EXECUÇÃO DO CÓDIGO
   H = 0.D0
   T = 0.5D0

   call BaseSpins(S, N)             ! BASE DE SPINS

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

   DO WHILE (H <= 10)

      PRINT '(F6.1, A)', H/10*100, '% '
      PRINT '(A)', '------------------------------------------------------'

      call Energy(S, N, H, SpinInteraction, Energia)                  ! AUTOENERGIAS
      call PartitionFunction(N, T, Z, Energia)                           ! FUNÇÃO DE PARTIÇÃO
      call Magnetization(N, S, T, Z, Energia, M)                         ! MAGNETIZAÇÃO

      WRITE(1, *)  H, M/N

      H = H + 0.01D0

   END DO

   DEALLOCATE(S, Energia, SpinInteraction)

END PROGRAM CLUSTER_FRUSTRATION

SUBROUTINE Energy(S, N, H, SpinInteraction, Energia)

   IMPLICIT NONE

   INTEGER(4) :: N
   INTEGER :: j
   REAL(8) :: Energia(2**N), S(N, 2**N), SpinInteraction(2**N), H

   Energia(:) = SpinInteraction(:)

   !$OMP PARALLEL DO
   DO j = 1, 2**N  ! LOOP DAS CONFIGURAÇÕES DE SPIN

      Energia(j) =  Energia(j) - H*SUM(S(:,j))

   END DO
   !$OMP END PARALLEL DO

   PRINT '(A, F6.1)', 'H =', H
   PRINT '(A)', 'ENERGIA'

END SUBROUTINE

SUBROUTINE PartitionFunction(N, T, Z, Energia)

   IMPLICIT NONE

   REAL(8) :: Energia(2**N), Z, T
   INTEGER(4) :: N
   REAL(8) :: minEnergy

   minEnergy = MINVAL(Energia(:))
   Z = SUM(exp(-(Energia(:) - minEnergy) / T))

   PRINT '(A)', 'FUNCAO DE PARTICAO'

END SUBROUTINE PartitionFunction

SUBROUTINE Magnetization(N, S, T, Z, Energia, M)

   IMPLICIT NONE

   INTEGER(4) :: N
   REAL(8) :: Energia(2**N), S(N, 2**N), M, Z, T
   REAL(8) :: inverse_Z, minEnergy

   inverse_Z = 1.D0 / Z
   minEnergy = MINVAL(Energia(:))

   M = inverse_Z * DOT_PRODUCT(exp(-(Energia(:) - minEnergy) / T), SUM(S, DIM=1))

   PRINT '(A, F6.1)', 'MAGNETIZACAO ', M/N

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

   !Cluster
   !               1
   !              / \
   !             3 - 2
   !            / \ / \
   !           4 - 5 - 6
   !          / \ / \ / \
   !        10 - 9 - 8 - 7
   !        / \ / \ / \ / \
   !       11-12- 13 -14 -15
   !      / \ / \ / \ / \ / \
   !     21- 20-19- 18 -17- 16

   if (ns.eq.21) then
      j=i+1
      if (j.eq.22) j=11
      if (i.eq.1) k=3
      if (i.eq.2) k=6
      if (i.eq.3) k=5
      if (i.eq.4) k=9
      if (i.eq.5) k=2
      if (i.eq.6) k=8
      if (i.eq.7) k=14
      if (i.eq.8) then
         k=5
         l=13
      endif
      if (i.eq.9) then
         k=5
         l=12
      endif
      if (i.eq.10)k=4
      if (i.eq.11)k=20
      if (i.eq.12)then
         k=10
         l=19
      endif
      if (i.eq.13) then
         k=9
         l=18
      endif
      if (i.eq.14) then
         k=8
         l=17
      endif
      if (i.eq.15)k=7

      if (i.eq.17)k=15
      if (i.eq.18)k=14
      if (i.eq.19)k=13
      if (i.eq.20)k=12


   end if

END


