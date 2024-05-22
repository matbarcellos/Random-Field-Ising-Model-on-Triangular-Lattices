PROGRAM CLUSTER_BIMODAL

   USE OMP_LIB
   IMPLICIT NONE

   ! DECLARAÇÃO DE VARIÁVEIS
   INTEGER, PARAMETER :: N = 9
   REAL(8) :: S(N, 2**N), R(N, 3**N)
   INTEGER(4) :: j, ii, ij, il, ik, it
   REAL(8) :: H
   REAL(8) :: SpinInteraction(2**N), Entropia(3**N), Energia(2**N, 3**N)
   CHARACTER(LEN=10) :: filename
   REAL(8), DIMENSION(5) :: T_values = [0.5D0, 0.5D0, 0.5D0, 0.5D0, 0.5D0]

   ! EXECUÇÃO DO CÓDIGO
   H = 0.D0

   call BaseSpins(S, N)             ! BASE DE SPINS
   call BaseCampos(R, N)       ! BASE DE CAMPOS ALEATÓRIOS

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

         call Energy(S, R, N, H, SpinInteraction, Energia)
         call Entropy(N, T_values(it), Energia, Entropia)

         WRITE(1, *)  H, SUM(Entropia)/(N*(3**N))

         H = H + 0.01D0

      END DO

      ! Reset H para o valor inicial
      H = 0.D0

      ! Fechar o arquivo antes de começar a próxima iteração
      CLOSE(UNIT=1)
   END DO


END PROGRAM CLUSTER_BIMODAL

SUBROUTINE Energy(S, R, N, H, SpinInteraction, Energia)

   IMPLICIT NONE

   INTEGER(4) :: N
   INTEGER :: i, j
   REAL(8) :: Energia(3**N, 2**N), S(N, 2**N), R(N, 3**N), SpinInteraction(2**N), H

   !$OMP PARALLEL DO
   DO i = 1, 3**N ! LOOP DAS CONFIGURAÇÕES DE CAMPO ALEATÓRIO

      Energia(i, :) = SpinInteraction(:)

      DO j = 1, 2**N  ! LOOP DAS CONFIGURAÇÕES DE SPIN

         Energia(i, j) =  Energia(i, j) - H*SUM(S(:,j)) - 0.5D0*SUM(R(:,i)*S(:,j))

      END DO
   END DO
   !$OMP END PARALLEL DO

   PRINT '(A, F6.1)', 'H =', H
   PRINT '(A)', 'ENERGIA'

END SUBROUTINE

SUBROUTINE PartitionFunction(N, T, Z, Energia)

   IMPLICIT NONE

   REAL(8) :: Energia(3**N, 2**N), Z(3**N), T
   INTEGER(4) :: N, i
   REAL(8) :: minEnergy

   !$OMP PARALLEL DO
   DO i = 1, 3**N

      minEnergy = MINVAL(Energia(i, :))
      Z(i) = SUM(exp(-(Energia(i, :) - minEnergy) / T))

   END DO
   !$OMP END PARALLEL DO

   PRINT '(A)', 'FUNCAO DE PARTICAO'

END SUBROUTINE PartitionFunction

SUBROUTINE Entropy(N, T, Energia, Entropia)

   IMPLICIT NONE

   INTEGER(4) :: N
   REAL(8) :: Energia(3**N, 2**N), T, T_h, T_l, h, Entropia(3**N), Ztot_h(3**N), F_h(3**N), Ztot_l(3**N), F_l(3**N)
   INTEGER(4) :: i, j

   ! MÉTODO DAS DIFERENÇAS FINITAS CENTRADAS PARA O CÁLCULO DA ENTROPIA

   h = 1.0e-4

   ! CALCULANDO PARA T + h (T superior)
   T_h = T + h

   !$OMP PARALLEL DO
   DO i = 1, 3**N

      Ztot_h(i) = 0.D0 ! EVITA SOMAR VALORES RESIDUAIS ARMAZENADOS NA MEMÓRIA

      DO j = 1, 2**N

         Ztot_h(i) = Ztot_h(i) + exp(-(Energia(i, j) - MINVAL(Energia(i, :)))/T_h)

      END DO

      F_h(i) = - T_h*LOG(Ztot_h(i))

   END DO
   !$OMP END PARALLEL DO

   ! CALCULANDO PARA T - h (T inferior)
   T_l = T - h

   !$OMP PARALLEL DO
   DO i = 1, 3**N

      Ztot_l(i) = 0.D0 ! EVITA SOMAR VALORES RESIDUAIS ARMAZENADOS NA MEMÓRIA

      DO j = 1, 2**N

         Ztot_l(i) = Ztot_l(i) + exp(-(Energia(i, j) - MINVAL(Energia(i, :)))/T_l)

      END DO

      F_l(i) = - T_l*LOG(Ztot_l(i))

   END DO
   !$OMP END PARALLEL DO

   ! PRIMEIRA DERIVADA DE F COM RELAÇÃO A T USANDO O MÉTODO DAS DIFERENÇAS FINITAS CENTRADAS
   !$OMP PARALLEL DO
   DO i = 1, 3**N

      Entropia(i) = - (F_h(i) - F_l(i)) / (2.0*h)

   END DO
   !$OMP END PARALLEL DO

END SUBROUTINE Entropy

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

SUBROUTINE BaseCampos(R, N)

   IMPLICIT NONE
   INTEGER(4) :: N
   INTEGER(4) :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14, i15, j
   REAL(8) :: R(N, 3**N)

   if(N.eq.3) then
      j=1 ! contador de configuração
      !base nunca ira mudar.
      DO i1= -1, 1
         DO i2= -1, 1
            DO i3= -1, 1
               R(3,j) = i3
               R(2,j) = i2
               R(1,j) = i1
               j= j+1
            END DO
         END DO
      END DO
   end if

   if(N.eq.6) then
      j=1 ! contador de configuração
      !base nunca ira mudar.
      DO i1= -1, 1
         DO i2= -1, 1
            DO i3= -1, 1
               DO i4 = -1, 1
                  DO i5 = -1, 1
                     DO i6 = -1, 1
                        R(6,j) = i6
                        R(5,j) = i5
                        R(4,j) = i4
                        R(3,j) = i3
                        R(2,j) = i2
                        R(1,j) = i1
                        j= j+1
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
   end if

   if(N.eq.9) then
      j=1 ! contador de configuração
      !base nunca ira mudar.
      DO i1= -1, 1
         DO i2= -1, 1
            DO i3= -1, 1
               DO i4 = -1, 1
                  DO i5 = -1, 1
                     DO i6 = -1, 1
                        DO i7 = -1, 1
                           DO i8 = -1, 1
                              DO i9 = -1, 1
                                 R(9, j) = i9
                                 R(8, j) = i8
                                 R(7,j) = i7
                                 R(6,j) = i6
                                 R(5,j) = i5
                                 R(4,j) = i4
                                 R(3,j) = i3
                                 R(2,j) = i2
                                 R(1,j) = i1
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
   end if

   if(N.eq.15) then
      j=1 ! contador de configuração
      !base nunca ira mudar.
      DO i1= -1, 1
         DO i2= -1, 1
            DO i3= -1, 1
               DO i4 = -1, 1
                  DO i5 = -1, 1
                     DO i6 = -1, 1
                        DO i7 = -1, 1
                           DO i8 = -1, 1
                              DO i9 = -1, 1
                                 R(15, j) = i15
                                 R(14, j) = i14
                                 R(13, j) = i13
                                 R(12, j) = i12
                                 R(11, j) = i11
                                 R(10, j) = i10
                                 R(9, j) = i9
                                 R(8, j) = i8
                                 R(7,j) = i7
                                 R(6,j) = i6
                                 R(5,j) = i5
                                 R(4,j) = i4
                                 R(3,j) = i3
                                 R(2,j) = i2
                                 R(1,j) = i1
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
   end if

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


