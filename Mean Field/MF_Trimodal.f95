PROGRAM MEANFIELD
    
    IMPLICIT NONE

    INTEGER(4), PARAMETER :: N = 3
    REAL(8) :: S(N, 2**N), R(N, 3**N)
    REAL(8) :: T, h, J, h_0
    REAL(8) :: m(3**N, N)
    
    J =  -1.0D0
    h =   0.00D0
    T =   0.1D0
    h_0 = 0.5D0

   ! Essa é uma subrotina 
    CALL SPIN(S, N) 
    CALL CA(R, N, h_0)
    
    DO WHILE (h <= 10.0D0)
        
        PRINT '(F6.1, A)', H/10*100, '% '

        CALL Self_Consistency(S, R, N, h, J, T, m)
    
        WRITE(1,'(F6.2, F8.2)') h, SUM(m)/(N*(3**N))
    
        h = h + 0.01D0
    
    END DO
    
END PROGRAM 
    
    
SUBROUTINE Self_Consistency(S, R, N, h, J, T, m)
    
        INTEGER(4) :: i, ii, k
        REAL(8) :: S(N, 2**N), R(N, 3**N)
        REAL(8) :: E(3**N, 2**N), Z(3**N)
        REAL(8) :: h, J, T
        REAL(8) :: prec, error
        REAL(8) :: iter
        REAL(8) :: m(3**N, 3), mred(3**N, 3)

        m(:, 1) =   1.0D0
        m(:, 2) =  -1.0D0
        m(:, 3) =   1.0D0

        error = 1.0D0
        prec = 1.0e-6
        iter = 0.0D0
           
        DO WHILE (error >= prec .AND. iter <= 10000)

            E = 0.0D0

            DO ii = 1, 3**N
                DO i = 1, 2**N
                    E(ii, i) = - J*(S(1,i)*S(2,i) + S(2,i)*S(3,i) + S(3,i)*S(1,i)) &
                    - h*SUM(S(:,i)) - SUM(R(:,ii)*S(:,i)) &
                    - J*((2*(m(ii, 2)+m(ii, 3))*S(1, i)) + (2*(m(ii, 1)+m(ii, 3))*S(2, i)) + (2*(m(ii, 1)+m(ii, 2))*S(3, i)))
                END DO
            END DO

            Z = 0.0D0
            
            DO ii = 1, 3**N   
                Z(ii) = SUM(exp(-E(ii,:)/T))
            END DO
            
          !print*, H, invZ
          !read(*,*)
    
            mred = 0.0D0

            DO ii = 1, 3**N
                DO k = 1, N
                    mred(ii, k) = mred(ii, k) + (1D0/Z(ii))*SUM(S(k,:)*exp(-(E(ii,:))/T))
                END DO
            END DO
    
            !print*, mred
            !read(*,*)
            
            error = MAXVAL(ABS(mred - m))
    
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
    
SUBROUTINE CA(R, N, H_0)

    IMPLICIT NONE
    INTEGER(4) :: N
    INTEGER(4) :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14, i15, j
    REAL(8) :: H_0, R(N, 3**N)
    

    if(N.eq.3) then
       j=1 ! contador de configuração
       !base nunca ira mudar.
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
    end if
 
END SUBROUTINE 



