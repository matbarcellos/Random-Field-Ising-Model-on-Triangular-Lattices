PROGRAM BimodalGaussian

   IMPLICIT NONE

   ! DECLARAÇÃO DE PARAMETROS
   INTEGER, PARAMETER :: Ns = 3
   INTEGER :: dim_num, seed, eval_num

   ! DECLARAÇÃO DE VARIÁVEIS
   INTEGER :: i, j, ii, ij, il, ik
   REAL(8) :: S(Ns, 2**Ns), H, T, Energia(2**Ns), Z, M
   REAL(8) :: x(Ns), x1(Ns), x2(Ns), x3(Ns), SpinInteraction(2**Ns), result, sigma1, mu1, sigma2, mu2, sigma3, mu3, value

   ! EXECUÇÃO DO CÓDIGO
   dim_num = Ns
   seed = 123456
   eval_num = 10000
   H = 0.D+00
   T = 0.01D+00

   sigma1 = 5D0
   mu1 = 0.5D0

   sigma2 = 5D0
   mu2 = -0.5D0

   sigma3 = 1D0
   mu3 = 0D0

   call BaseSpins(S, Ns) ! BASE DE SPINS


   !***********************************************************************
   DO j = 1, 2**Ns
      DO ii = 1, Ns

         call rede(Ns, ii, ij, il, ik) ! ij, il e ik são vizinhos de ii

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



   DO WHILE (H <= 10)

      print*, H/10*100, '%'

      result = 0.0D+00

      DO j = 1, eval_num  ! LOOP DAS CONFIGURAÇÕES DE SPIN

         call r8vec_normal_01(dim_num, seed, x1)                    ! X: vetor de número aleatórios
         call r8vec_normal_02(dim_num, seed, x2)
         call r8vec_normal_03(dim_num, seed, x3)

         !x1 = sigma1*x1 + mu1
         !x2 = sigma2*x2 + mu2

         ! GERAR UM NUMÉRIO ALEATÓRIO ENTRE 0 E 1
         call RANDOM_SEED()
         call RANDOM_NUMBER(value)

         ! SELECIONA UM DOS VALORES
         IF (value <= 0.5D0) THEN
            DO i = 1, Ns
               x(i) = 0.5D0*(x1(i)*sigma1 + mu1)
            END DO
         ELSE IF (value <= 1D0) THEN
            DO i = 1, Ns
               x(i) = 0.5D0*(x2(i)*sigma2 + mu2)
            END DO
         ELSE
            DO i = 1, Ns
            x(i) = x3(i)*sigma3 + mu3
            END DO
         END IF

         call Energy(S, X, Ns, H, SpinInteraction, Energia)                 ! AUTOENERGIAS
         call PartitionFunction(Ns, T, Z, Energia)                   ! FUNÇÃO DE PARTIÇÃO
         call Magnetization(Ns, S, T, Z, Energia, M)                 ! MAGNETIZAÇÃO

         result = result + M

      END DO

      result = result / REAL(eval_num*Ns, kind = 8 )

      WRITE(1, *)  H, result                                          ! MAGNETIZAÇÃO MÉDIA POR SPIN E POR DESORDEM VS CAMPO MAGNÉTICO EXTERNO

      H = H + 0.1D0

   END DO

END PROGRAM BimodalGaussian

SUBROUTINE Energy(S, X, Ns, H, SpinInteraction, Energia)

   IMPLICIT NONE

   INTEGER :: Ns, j
   REAL(8) :: X(Ns), Energia(2**Ns), S(Ns, 2**Ns), SpinInteraction(2**Ns), H

   DO j = 1, 2**Ns  ! LOOP DAS CONFIGURAÇÕES DE SPIN


      Energia(j) = SpinInteraction(j) - H * SUM(S(:, j)) &
      - SUM(X(:)*S(:,j))

   END DO

END SUBROUTINE

SUBROUTINE PartitionFunction(Ns, T, Z, Energia)

   IMPLICIT NONE

   INTEGER :: Ns, j
   REAL(8) :: Energia(2**Ns), Z, T

   Z = 0.D0 ! EVITA SOMAR VALORES RESIDUAIS ARMAZENADOS NA MEMÓRIA

   DO j = 1, 2**Ns

      Z = Z + exp(-(Energia(j) - MINVAL(Energia))/T)

   END DO

END SUBROUTINE PartitionFunction

SUBROUTINE Magnetization(Ns, S, T, Z, Energia, M)

   IMPLICIT NONE

   INTEGER :: Ns, j
   REAL(8) :: Energia(2**Ns), S(Ns, 2**Ns), M, Z, T

   M = 0.0D+00

   DO j = 1, 2**Ns ! LOOP DAS CONFIGURAÇÕES DE SPIN

      M = M + (1/Z)*(SUM(S(:, j))*exp(-(Energia(j) - MINVAL(Energia))/T))

   END DO

END SUBROUTINE Magnetization

SUBROUTINE BaseSpins(S, Ns)

   IMPLICIT NONE

   INTEGER :: i, j, Ns
   REAL(8) :: S(Ns, 2**Ns)

   S = 1.D0

   DO j = 1, 2**Ns

      DO i = 1, Ns

         if (btest(j-1,i-1)) S(i,j) = -1.D0

      END DO

   END DO

END SUBROUTINE BaseSpins

SUBROUTINE r8vec_normal_01(n, seed, x1)

   implicit none

   integer ( kind = 4 ) n
   integer ( kind = 4 ) m
   integer ( kind = 4 ), save :: made = 0
   real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
   real ( kind = 8 ) r(n+1)
   real ( kind = 8 ) r8_uniform_01
   integer ( kind = 4 ), save :: saved = 0
   integer ( kind = 4 ) seed
   real ( kind = 8 ) x1(n)
   integer ( kind = 4 ) x1_hi_index1
   integer ( kind = 4 ) x1_lo_index1
   real ( kind = 8 ), save :: y = 0.0D+00

   if ( n < 0 ) then
      n = made
      made = 0
      saved = 0
      y = 0.0D+00
      return
   else if ( n == 0 ) then
      return
   end if
   !
   !  Record the range of x1 we need to fill in.
   !
   x1_lo_index1 = 1
   x1_hi_index1 = n
   !
   !  Use up the old value, if we have it.
   !
   if ( saved == 1 ) then
      x1(1) = y
      saved = 0
      x1_lo_index1 = 2
   end if
   !
   !  Maybe we don't need any more values.
   !
   if ( x1_hi_index1 - x1_lo_index1 + 1 == 0 ) then
      !
      !  If we need just one new value, do that here to avoid null arrays.
      !
   else if ( x1_hi_index1 - x1_lo_index1 + 1 == 1 ) then

      r(1) = r8_uniform_01 ( seed )

      if ( r(1) == 0.0D+00 ) then
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'R8VEC_NORMAL_01 - Fatal error!'
         write ( *, '(a)' ) '  R8_UNIFORM_01 returned a value of 0.'
         stop
      end if

      r(2) = r8_uniform_01 ( seed )

      x1(x1_hi_index1) = &
         sqrt ( -2.0D+00 * log ( r(1) ) ) * cos ( 2.0D+00 * pi * r(2) )
      y =      sqrt ( -2.0D+00 * log ( r(1) ) ) * sin ( 2.0D+00 * pi * r(2) )

      saved = 1

      made = made + 2
      !
      !  If we require an even number of values, that's easy.
      !
   else if ( mod ( x1_hi_index1 - x1_lo_index1 + 1, 2 ) == 0 ) then

      m = ( x1_hi_index1 - x1_lo_index1 + 1 ) / 2

      call r8vec_uniform_01 ( 2 * m, seed, r )

      x1(x1_lo_index1:x1_hi_index1-1:2) = &
         sqrt ( -2.0D+00 * log ( r(1:2*m-1:2) ) ) &
         * cos ( 2.0D+00 * pi * r(2:2*m:2) )

      x1(x1_lo_index1+1:x1_hi_index1:2) = &
         sqrt ( -2.0D+00 * log ( r(1:2*m-1:2) ) ) &
         * sin ( 2.0D+00 * pi * r(2:2*m:2) )

      made = made + x1_hi_index1 - x1_lo_index1 + 1

   else

      x1_hi_index1 = x1_hi_index1 - 1

      m = ( x1_hi_index1 - x1_lo_index1 + 1 ) / 2 + 1

      call r8vec_uniform_01 ( 2 * m, seed, r )

      x1(x1_lo_index1:x1_hi_index1-1:2) = &
         sqrt ( -2.0D+00 * log ( r(1:2*m-3:2) ) ) &
         * cos ( 2.0D+00 * pi * r(2:2*m-2:2) )

      x1(x1_lo_index1+1:x1_hi_index1:2) = &
         sqrt ( -2.0D+00 * log ( r(1:2*m-3:2) ) ) &
         * sin ( 2.0D+00 * pi * r(2:2*m-2:2) )

      x1(n) = sqrt ( -2.0D+00 * log ( r(2*m-1) ) ) &
         * cos ( 2.0D+00 * pi * r(2*m) )

      y = sqrt ( -2.0D+00 * log ( r(2*m-1) ) ) &
         * sin ( 2.0D+00 * pi * r(2*m) )

      saved = 1

      made = made + x1_hi_index1 - x1_lo_index1 + 2

   end if

   return
END SUBROUTINE r8vec_normal_01

FUNCTION r8_uniform_01(seed)


   implicit none

   integer ( kind = 4 ) k
   real ( kind = 8 ) r8_uniform_01
   integer ( kind = 4 ) seed

   if ( seed == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
      write ( *, '(a)' ) '  Input value of SEED = 0.'
      stop
   end if

   k = seed / 127773

   seed = 16807 * ( seed - k * 127773 ) - k * 2836

   if ( seed < 0 ) then
      seed = seed + 2147483647
   end if
   !
   !  Although SEED can be represented ex1actly as a 32 bit integer,
   !  it generally cannot be represented ex1actly as a 32 bit real number!
   !
   r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

   return
END FUNCTION r8_uniform_01

SUBROUTINE r8vec_uniform_01 (n, seed, r)

   implicit none

   integer ( kind = 4 ) n

   integer ( kind = 4 ) i
   integer ( kind = 4 ) k
   integer ( kind = 4 ) seed
   real ( kind = 8 ) r(n)

   if ( seed == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
      write ( *, '(a)' ) '  Input value of SEED = 0.'
      stop
   end if

   do i = 1, n

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
         seed = seed + huge ( seed )
      end if

      r(i) = real ( seed, kind = 8 ) * 4.656612875D-10

   end do

   return
END SUBROUTINE r8vec_uniform_01

SUBROUTINE r8vec_normal_02(n, seed, x2)

   implicit none

   integer ( kind = 4 ) n
   integer ( kind = 4 ) m
   integer ( kind = 4 ), save :: made = 0
   real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
   real ( kind = 8 ) r(n+1)
   real ( kind = 8 ) r8_uniform_02
   integer ( kind = 4 ), save :: saved = 0
   integer ( kind = 4 ) seed
   real ( kind = 8 ) x2(n)
   integer ( kind = 4 ) x2_hi_index2
   integer ( kind = 4 ) x2_lo_index2
   real ( kind = 8 ), save :: y = 0.0D+00

   if ( n < 0 ) then
      n = made
      made = 0
      saved = 0
      y = 0.0D+00
      return
   else if ( n == 0 ) then
      return
   end if
   !
   !  Record the range of x2 we need to fill in.
   !
   x2_lo_index2 = 1
   x2_hi_index2 = n
   !
   !  Use up the old value, if we have it.
   !
   if ( saved == 1 ) then
      x2(1) = y
      saved = 0
      x2_lo_index2 = 2
   end if
   !
   !  Maybe we don't need any more values.
   !
   if ( x2_hi_index2 - x2_lo_index2 + 1 == 0 ) then
      !
      !  If we need just one new value, do that here to avoid null arrays.
      !
   else if ( x2_hi_index2 - x2_lo_index2 + 1 == 1 ) then

      r(1) = r8_uniform_02 ( seed )

      if ( r(1) == 0.0D+00 ) then
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'R8VEC_NORMAL_01 - Fatal error!'
         write ( *, '(a)' ) '  R8_UNIFORM_02 returned a value of 0.'
         stop
      end if

      r(2) = r8_uniform_02 ( seed )

      x2(x2_hi_index2) = &
         sqrt ( -2.0D+00 * log ( r(1) ) ) * cos ( 2.0D+00 * pi * r(2) )
      y =      sqrt ( -2.0D+00 * log ( r(1) ) ) * sin ( 2.0D+00 * pi * r(2) )

      saved = 1

      made = made + 2
      !
      !  If we require an even number of values, that's easy.
      !
   else if ( mod ( x2_hi_index2 - x2_lo_index2 + 1, 2 ) == 0 ) then

      m = ( x2_hi_index2 - x2_lo_index2 + 1 ) / 2

      call r8vec_uniform_02 ( 2 * m, seed, r )

      x2(x2_lo_index2:x2_hi_index2-1:2) = &
         sqrt ( -2.0D+00 * log ( r(1:2*m-1:2) ) ) &
         * cos ( 2.0D+00 * pi * r(2:2*m:2) )

      x2(x2_lo_index2+1:x2_hi_index2:2) = &
         sqrt ( -2.0D+00 * log ( r(1:2*m-1:2) ) ) &
         * sin ( 2.0D+00 * pi * r(2:2*m:2) )

      made = made + x2_hi_index2 - x2_lo_index2 + 1

   else

      x2_hi_index2 = x2_hi_index2 - 1

      m = ( x2_hi_index2 - x2_lo_index2 + 1 ) / 2 + 1

      call r8vec_uniform_02 ( 2 * m, seed, r )

      x2(x2_lo_index2:x2_hi_index2-1:2) = &
         sqrt ( -2.0D+00 * log ( r(1:2*m-3:2) ) ) &
         * cos ( 2.0D+00 * pi * r(2:2*m-2:2) )

      x2(x2_lo_index2+1:x2_hi_index2:2) = &
         sqrt ( -2.0D+00 * log ( r(1:2*m-3:2) ) ) &
         * sin ( 2.0D+00 * pi * r(2:2*m-2:2) )

      x2(n) = sqrt ( -2.0D+00 * log ( r(2*m-1) ) ) &
         * cos ( 2.0D+00 * pi * r(2*m) )

      y = sqrt ( -2.0D+00 * log ( r(2*m-1) ) ) &
         * sin ( 2.0D+00 * pi * r(2*m) )

      saved = 1

      made = made + x2_hi_index2 - x2_lo_index2 + 2

   end if

   return
END SUBROUTINE r8vec_normal_02

FUNCTION r8_uniform_02(seed)


   implicit none

   integer ( kind = 4 ) k
   real ( kind = 8 ) r8_uniform_02
   integer ( kind = 4 ) seed

   if ( seed == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_UNIFORM_02 - Fatal error!'
      write ( *, '(a)' ) '  Input value of SEED = 0.'
      stop
   end if

   k = seed / 127773

   seed = 16807 * ( seed - k * 127773 ) - k * 2836

   if ( seed < 0 ) then
      seed = seed + 2147483647
   end if
   !
   !  Although SEED can be represented ex2actly as a 32 bit integer,
   !  it generally cannot be represented ex2actly as a 32 bit real number!
   !
   r8_uniform_02 = real ( seed, kind = 8 ) * 4.656612875D-10

   return
END FUNCTION r8_uniform_02

SUBROUTINE r8vec_uniform_02 (n, seed, r)

   implicit none

   integer ( kind = 4 ) n

   integer ( kind = 4 ) i
   integer ( kind = 4 ) k
   integer ( kind = 4 ) seed
   real ( kind = 8 ) r(n)

   if ( seed == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8VEC_UNIFORM_02 - Fatal error!'
      write ( *, '(a)' ) '  Input value of SEED = 0.'
      stop
   end if

   do i = 1, n

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
         seed = seed + huge ( seed )
      end if

      r(i) = real ( seed, kind = 8 ) * 4.656612875D-10

   end do

   return
END SUBROUTINE r8vec_uniform_02

SUBROUTINE r8vec_normal_03(n, seed, x3)

   implicit none

   integer ( kind = 4 ) n
   integer ( kind = 4 ) m
   integer ( kind = 4 ), save :: made = 0
   real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
   real ( kind = 8 ) r(n+1)
   real ( kind = 8 ) r8_uniform_03
   integer ( kind = 4 ), save :: saved = 0
   integer ( kind = 4 ) seed
   real ( kind = 8 ) x3(n)
   integer ( kind = 4 ) x3_hi_index3
   integer ( kind = 4 ) x3_lo_index3
   real ( kind = 8 ), save :: y = 0.0D+00

   if ( n < 0 ) then
      n = made
      made = 0
      saved = 0
      y = 0.0D+00
      return
   else if ( n == 0 ) then
      return
   end if
   !
   !  Record the range of x3 we need to fill in.
   !
   x3_lo_index3 = 1
   x3_hi_index3 = n
   !
   !  Use up the old value, if we have it.
   !
   if ( saved == 1 ) then
      x3(1) = y
      saved = 0
      x3_lo_index3 = 2
   end if
   !
   !  Maybe we don't need any more values.
   !
   if ( x3_hi_index3 - x3_lo_index3 + 1 == 0 ) then
      !
      !  If we need just one new value, do that here to avoid null arrays.
      !
   else if ( x3_hi_index3 - x3_lo_index3 + 1 == 1 ) then

      r(1) = r8_uniform_03 ( seed )

      if ( r(1) == 0.0D+00 ) then
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'R8VEC_NORMAL_03 - Fatal error!'
         write ( *, '(a)' ) '  R8_UNIFORM_03 returned a value of 0.'
         stop
      end if

      r(2) = r8_uniform_03 ( seed )

      x3(x3_hi_index3) = &
         sqrt ( -2.0D+00 * log ( r(1) ) ) * cos ( 2.0D+00 * pi * r(2) )
      y =      sqrt ( -2.0D+00 * log ( r(1) ) ) * sin ( 2.0D+00 * pi * r(2) )

      saved = 1

      made = made + 2
      !
      !  If we require an even number of values, that's easy.
      !
   else if ( mod ( x3_hi_index3 - x3_lo_index3 + 1, 2 ) == 0 ) then

      m = ( x3_hi_index3 - x3_lo_index3 + 1 ) / 2

      call r8vec_uniform_03 ( 2 * m, seed, r )

      x3(x3_lo_index3:x3_hi_index3-1:2) = &
         sqrt ( -2.0D+00 * log ( r(1:2*m-1:2) ) ) &
         * cos ( 2.0D+00 * pi * r(2:2*m:2) )

      x3(x3_lo_index3+1:x3_hi_index3:2) = &
         sqrt ( -2.0D+00 * log ( r(1:2*m-1:2) ) ) &
         * sin ( 2.0D+00 * pi * r(2:2*m:2) )

      made = made + x3_hi_index3 - x3_lo_index3 + 1

   else

      x3_hi_index3 = x3_hi_index3 - 1

      m = ( x3_hi_index3 - x3_lo_index3 + 1 ) / 2 + 1

      call r8vec_uniform_03 ( 2 * m, seed, r )

      x3(x3_lo_index3:x3_hi_index3-1:2) = &
         sqrt ( -2.0D+00 * log ( r(1:2*m-3:2) ) ) &
         * cos ( 2.0D+00 * pi * r(2:2*m-2:2) )

      x3(x3_lo_index3+1:x3_hi_index3:2) = &
         sqrt ( -2.0D+00 * log ( r(1:2*m-3:2) ) ) &
         * sin ( 2.0D+00 * pi * r(2:2*m-2:2) )

      x3(n) = sqrt ( -2.0D+00 * log ( r(2*m-1) ) ) &
         * cos ( 2.0D+00 * pi * r(2*m) )

      y = sqrt ( -2.0D+00 * log ( r(2*m-1) ) ) &
         * sin ( 2.0D+00 * pi * r(2*m) )

      saved = 1

      made = made + x3_hi_index3 - x3_lo_index3 + 2

   end if

   return
END SUBROUTINE r8vec_normal_03

FUNCTION r8_uniform_03(seed)


   implicit none

   integer ( kind = 4 ) k
   real ( kind = 8 ) r8_uniform_03
   integer ( kind = 4 ) seed

   if ( seed == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_UNIFORM_03 - Fatal error!'
      write ( *, '(a)' ) '  Input value of SEED = 0.'
      stop
   end if

   k = seed / 127773

   seed = 16807 * ( seed - k * 127773 ) - k * 2836

   if ( seed < 0 ) then
      seed = seed + 2147483647
   end if
   !
   !  Although SEED can be represented ex3actly as a 32 bit integer,
   !  it generally cannot be represented ex3actly as a 32 bit real number!
   !
   r8_uniform_03 = real ( seed, kind = 8 ) * 4.656612875D-10

   return
END FUNCTION r8_uniform_03

SUBROUTINE r8vec_uniform_03 (n, seed, r)

   implicit none

   integer ( kind = 4 ) n

   integer ( kind = 4 ) i
   integer ( kind = 4 ) k
   integer ( kind = 4 ) seed
   real ( kind = 8 ) r(n)

   if ( seed == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8VEC_UNIFORM_03 - Fatal error!'
      write ( *, '(a)' ) '  Input value of SEED = 0.'
      stop
   end if

   do i = 1, n

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
         seed = seed + huge ( seed )
      end if

      r(i) = real ( seed, kind = 8 ) * 4.656612875D-10

   end do

   return
END SUBROUTINE r8vec_uniform_03

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
