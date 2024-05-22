PROGRAM M_H_RF

INTEGER(4), PARAMETER :: N = 6
INTEGER(4) :: j, ii, ij, il, ik
REAL(8) :: S(N, N**2), R(N, N**2)
REAL(8) :: E(N**2, N**2), Z(2**N), M(2**N)
REAL(8) :: SpinInteraction(2**N)
REAL(8) :: H_0, h, T

H_0 = 0.00D0
T = 0.000001D0

CALL BaseSpins(S, N)

!***********************************************************************
DO j = 1, 2**N
    DO ii = 1, N

       call Rede(N, ii, ij, il, ik) ! ij, il e ik são vizinhos de ii

       !********************************************
       !  TESTE DE FUNCIONAMENTO DA SUBROTINA REDE
       !  print*, ii, ij, il, ik
       !  read(*,*)
       !********************************************

       IF(ij>0)SpinInteraction(j) = SpinInteraction(j) + S(ii, j)*S(ij, j)
       IF(il>0)SpinInteraction(j) = SpinInteraction(j) + S(ii, j)*S(il, j)
       IF(ik>0)SpinInteraction(j) = SpinInteraction(j) + S(ii, j)*S(ik, j)

    END DO
END DO
!***********************************************************************

DO WHILE (H_0 <= 10)

    CALL BaseCampos(R, N, H_0)

    h = 0.00D0

    DO WHILE(h <= 10)

        CALL Energy(S, R, N, H, SpinInteraction, E)
        CALL PartitionFunction(N, T, Z, E)
        CALL Magnetization(N, S, T, Z, E, M)

        h = h + 0.1

        WRITE(15,'(F10.2, F10.2, F10.2)') H_0, H, SUM(M)/(N*(2**N))

    END DO

    H_0 = H_0 + 0.1

END DO


END PROGRAM M_H_RF


SUBROUTINE Energy(S, R, N, H, SpinInteraction, E)

    IMPLICIT NONE
 
    INTEGER(4) :: N, i, j
    REAL(8) :: E(2**N, 2**N), S(N, 2**N), R(N, 2**N), SpinInteraction(2**N), H
 
    DO i = 1, 2**N
 
       E(i, :) = SpinInteraction(:)
 
       DO j = 1, 2**N
          E(i, j) = E(i, j) - H*SUM(S(:,j)) - SUM(R(:,i)*S(:,j))
       END DO
 
    END DO
 
END SUBROUTINE Energy

SUBROUTINE PartitionFunction(N, T, Z, E)

    IMPLICIT NONE
 
    REAL(8) :: E(2**N, 2**N), Z(2**N), T
    INTEGER(4) :: N, i
    REAL(8) :: minEnergy
 
    DO i = 1, 2**N
 
       minEnergy = MINVAL(E(i, :))
       Z(i) = SUM(exp(-(E(i, :) - minEnergy) / T))
 
    END DO
 
END SUBROUTINE PartitionFunction

SUBROUTINE Magnetization(N, S, T, Z, E, M)
    IMPLICIT NONE
 
    INTEGER(4) :: N, i
    REAL(8) :: E(2**N, 2**N), S(N, 2**N), M(2**N), Z(2**N), T
    REAL(8) :: inverse_Z, minEnergy
 
    DO i = 1, 2**N
 
       inverse_Z = 1.D0 / Z(i)
       minEnergy = MINVAL(E(i, :))
 
       M(i) = inverse_Z * DOT_PRODUCT(exp(-(E(i, :) - minEnergy) / T), SUM(S, DIM=1))
 
    END DO
 
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
 
SUBROUTINE Rede(ns, i, j, l, k)
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
 
END SUBROUTINE
 