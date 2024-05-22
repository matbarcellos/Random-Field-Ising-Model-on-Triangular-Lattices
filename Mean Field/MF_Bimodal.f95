PROGRAM MEANFIELD
    
    IMPLICIT NONE

    INTEGER(4), PARAMETER :: N = 3
    REAL(8) :: S(N, N**2), R(N, 2**N)
    REAL(8) :: T, h, J, h_0
    REAL(8) :: m(2**N, 3)
    
    J = -1.0D0
    h =  0D0
    T =  0.1D0
    h_0 = 0D0
    
    m(:, 1) =    0.01D0
    m(:, 2) =   -1D0
    m(:, 3) =    1D0

   ! Essa é uma subrotina 
    CALL SPIN(S, N) 
    CALL CA(R, N, h_0)
    
    DO WHILE (h <= 7.0D0)
        
        PRINT '(F6.1, A)', H/10*100, '% '

        CALL Self_Consistency(S, R, N, h, J, T, m)
    
        WRITE(1, *) h, SUM(m)/(N*(2**N))
    
        h = h + 0.01D0
    
    END DO
    
END PROGRAM 
    
    
SUBROUTINE Self_Consistency(S, R, N, h, J, T, m)
    
    !***********************************************************************    
    ! DECLARAÇÃO DE VARIÁVEIS DA SUBROTINA
    !***********************************************************************
        INTEGER(4) :: i, ii, k
        REAL(8) :: S(N, N**2), R(N, 2**N)
        REAL(8) :: E(2**N, 2**N), Z(2**N), invZ(2**N)
        REAL(8) :: h, J, T
        REAL(8) :: prec, error, minE(2**N)
        REAL(8) :: iter
        REAL(8) :: m(2**N, 3), mred(2**N, 3)
    !***********************************************************************
        error = 1.0D0
        prec = 1.0e-8
        iter = 0
    
        !DO WHILE (error >= prec .AND. iter <= 10000)
           
        DO WHILE (error >= prec)

            DO ii = 1, 2**N
                DO i = 1, 2**N
                    E(ii, i) = - J*(S(1,i)*S(2,i) + S(2,i)*S(3,i) + S(3,i)*S(1,i)) &
                    - h*SUM(S(:,i)) - SUM(R(:,ii)*S(:,i)) &
                    - J*((2*(m(ii, 2)+m(ii, 3))*S(1, i)) + (2*(m(ii, 1)+m(ii, 3))*S(2, i)) + (2*(m(ii, 1)+m(ii, 2))*S(3, i)))
                END DO
            END DO
            
            DO ii = 1, 2**N

                minE(ii) = MINVAL(E(ii,:))
                
                Z(ii) = SUM(exp(-(E(ii,:))/T))
            
                invZ(ii) = 1D0/Z(ii)

            END DO
            
          !print*, H, invZ
          !read(*,*)
    
            mred = 0.0D0

            DO ii = 1, 2**N
                DO k = 1, N
                    mred(ii, k) = mred(ii, k) + invZ(ii)*SUM(S(k,:)*exp(-(E(ii,:))/T))
                END DO
            END DO
    
            !print*, mred
            !read(*,*)
            
          error = MAXVAL(ABS(mred - m))
    
          m = mred
            
          iter = iter + 1
    
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
 
    REAL(8) :: R(N, 2**N), H_0
    INTEGER(4) :: i, j, N
 
    R = H_0
 
    DO j = 1, 2**N
 
       DO i = 1, N
 
          if (btest(j-1,i-1)) R(i,j) = - H_0
 
       END DO
 
    END DO
 
END SUBROUTINE 




