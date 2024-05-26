PROGRAM RandomFieldVariation

   USE OMP_LIB

   IMPLICIT NONE

   INTEGER, PARAMETER :: N = 3
   INTEGER, PARAMETER :: dim = 3
   INTEGER(4) :: j, ii, ij, il, ik, it
   REAL(8), ALLOCATABLE :: S(:,:), R(:,:), M(:), Z(:), E(:,:), SS(:)
   REAL(8) :: T, H
   CHARACTER(LEN=10) :: filename
   REAL(8), DIMENSION(6) :: H0_values = [0.00D0, 0.25D0, 0.50D0, 0.75D0, 1.00D0, 1.25D0]

   ALLOCATE(S(N, dim**N), R(N, dim**N), M(dim**N), Z(dim**N), E(dim**N, 2**N), SS(2**N))

   T = 0.00001D0
   H = 0.0D0

   CALL Spin(S, N)  ! Create an array with all spin configurations

   ! Calculates the interaction energy between the lattice spins
   DO j = 1, 2**N
      DO ii = 1, N

         ! Find the neighbors ij, il, ik of site ii
         CALL Lattice(N, ii, ij, il, ik)

         !********************************************
         !  print*, ii, ij, il, ik
         !  read(*,*)
         !********************************************

         IF(ij>0)SS(j) = SS(j) + S(ii, j)*S(ij, j)
         IF(il>0)SS(j) = SS(j) + S(ii, j)*S(il, j)
         IF(ik>0)SS(j) = SS(j) + S(ii, j)*S(ik, j)

      END DO
   END DO

   DO IT = 1, 6

      IF(dim .EQ. 2) THEN
         WRITE(filename, "('Bi',I1,'.dat')") IT
         OPEN(UNIT=it, FILE=filename)
      END IF

      IF(dim .EQ. 3) THEN
         WRITE(filename, "('Tri',I1,'.dat')") IT
         OPEN(UNIT=it, FILE=filename)
      END IF

      CALL RandomField(R, N, dim, H0_values(IT))

      H = 0.0D0

      DO WHILE (H <= 10)

         CALL Energy(S, R, N, dim, H, SS, E)
         CALL PartitionFunction(N, dim, T, Z, E)
         CALL Magnetization(N, dim, S, T, Z, E, M)

         WRITE(IT, *)  H, SUM(M)/(N*(dim**N))

         H = H + 0.01D0

      END DO
   END DO

   DEALLOCATE(S, R, M, Z, E, SS)

END PROGRAM

SUBROUTINE Energy(S, R, N, dim, H, SS, E)

   IMPLICIT NONE

   INTEGER(4) :: N, dim, i, j
   REAL(8) :: E(dim**N, 2**N), S(N, 2**N), R(N, dim**N), SS(2**N), H

   !$OMP PARALLEL DO
   DO i = 1, dim**N

      E(i, :) = SS(:)

      DO j = 1, 2**N
         E(i, j) = E(i, j) - H*SUM(S(:,j)) - SUM(R(:,i)*S(:,j))
      END DO

   END DO
   !$OMP END PARALLEL DO

END SUBROUTINE Energy

SUBROUTINE PartitionFunction(N, dim, T, Z, E)

   IMPLICIT NONE

   REAL(8) :: E(dim**N, 2**N), Z(dim**N), T
   INTEGER(4) :: N, dim, i
   REAL(8) :: minEnergy

   !$OMP PARALLEL DO
   DO i = 1, dim**N

      minEnergy = MINVAL(E(i, :))
      Z(i) = SUM(exp(-(E(i, :) - minEnergy) / T))

   END DO
   !$OMP END PARALLEL DO

END SUBROUTINE PartitionFunction

SUBROUTINE Magnetization(N, dim, S, T, Z, E, M)
   IMPLICIT NONE

   INTEGER(4) :: N, dim, i
   REAL(8) :: E(dim**N, 2**N), S(N, 2**N), M(dim**N), Z(dim**N), T
   REAL(8) :: inverse_Z, minEnergy

   DO i = 1, dim**N

      inverse_Z = 1.D0 / Z(i)
      minEnergy = MINVAL(E(i, :))

      M(i) = inverse_Z * DOT_PRODUCT(exp(-(E(i, :) - minEnergy) / T), SUM(S, DIM=1))

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

SUBROUTINE RandomField(R, N, dim, H_0)

   IMPLICIT NONE

   REAL(8) :: R(N, dim**N), H_0
   INTEGER(4) :: N, dim

   IF(dim .EQ. 2) CALL Bimodal(R, N, dim, H_0)

   IF(dim .EQ. 3) CALL Trimodal(R, N, dim, H_0)

END

SUBROUTINE Bimodal(R, N, dim, H_0)

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

SUBROUTINE Trimodal(R, N, dim, H_0)

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

SUBROUTINE Lattice(ns, i, j, l, k)
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

   IF (ns .EQ. 3) then
      j=i+1
      IF (i .EQ. 3) j=1

   END IF

   !Cluster
   !             1
   !            / \
   !           6 - 2
   !          / \ / \
   !         5 - 4 - 3


   IF (ns .EQ. 6) then
      j=i+1
      IF (i .EQ. 6) j=1
      IF (i .EQ. 2) k=4
      IF (i .EQ. 4) k=6
      IF (i .EQ. 6) k=2
   END IF

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

   IF (ns .EQ. 9) then
      IF (i .EQ. 1) then
         j=2
         k=8
      END IF

      IF (i .EQ. 2) then
         j=3
         k=9
      END IF

      IF (i .EQ. 3) then
         j=4
         k=9
      END IF

      IF (i .EQ. 4) then
         j=6
         k=9
      END IF

      IF (i .EQ. 5) then
         j=4
         k=6
      END IF
      IF (i .EQ. 6) then
         j=7
         k=9
      END IF
      IF (i .EQ. 7) then
         j=8
         k=9
      END IF
      IF (i .EQ. 8) then
         j=2
         k=9
      END IF
      !           print*, i, j, k
      !           read(*,*)

   END IF

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

   IF (ns .EQ. 15) then

      IF (i .EQ. 1) then
         j=9
         k=2
      END IF

      IF (i .EQ. 2) then
         j=10
         k=3
         l=9
      END IF

      IF (i .EQ. 3) then
         j=5
         k=4
         l=10
      END IF

      IF (i .EQ. 4) then
         j=14
         k=5
      END IF

      IF (i .EQ. 5) then
         j=13
         k=14
         l=6
      END IF

      IF (i .EQ. 6) then
         j=13
         k=12
         l=7
      END IF

      IF (i .EQ. 7) then
         j=12
         k=11
      END IF

      IF (i .EQ. 8) then
         j=7
         k=6
         l=10
      END IF

      IF (i .EQ. 9) then
         j=8
         k=10
      END IF

      IF (i .EQ. 10) then
         j=6
         k=5
      END IF

      IF (i .EQ. 11) then
         j=12
      END IF

      IF (i .EQ. 12) then
         j=13
      END IF

      IF (i .EQ. 13) then
         j=14
      END IF

      IF (i .EQ. 14) then
         j=15
      END IF

      IF (i .EQ. 15) then
         j=4
      END IF
   END IF

END

