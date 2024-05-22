PROGRAM SpinRandomField_N3

   IMPLICIT NONE

   ! DECLARAÇÃO DE PARAMETROS
   INTEGER, PARAMETER :: Ns = 6
   INTEGER :: dim_num, seed, eval_num

   ! DECLARAÇÃO DE VARIÁVEIS
   INTEGER :: j
   REAL(8) :: S(Ns, 2**Ns), H, T, Energia(2**Ns), Z, M
   REAL(8) :: x(6), result, delta

   ! EXECUÇÃO DO CÓDIGO
   dim_num = 6
   seed = 123456
   eval_num = 10000
   H = 0.D+00
   T = 0.1D+00

   call BaseSpins(S, Ns) ! BASE DE SPINS

   DO WHILE (H <= 10)

      print*, H/10*100, '%'

      result = 0.0D+00

      DO j = 1, eval_num  ! LOOP DAS CONFIGURAÇÕES DE SPIN

         call r8vec_normal_01 (dim_num, seed, x)                    ! X: vetor de número aleatórios

         delta = 0.5D0

         x = x*delta

         call Energy(S, x, dim_num, Ns, H, Energia)                  ! AUTOENERGIAS
         call PartitionFunction(Ns, T, Z, Energia)                   ! FUNÇÃO DE PARTIÇÃO
         call Magnetization(Ns, S, T, Z, Energia, M)                 ! MAGNETIZAÇÃO

         result = result + M

      END DO

      result = result / REAL(eval_num*6.D0, kind = 8 )

      WRITE(25, *)  H, result                                          ! MAGNETIZAÇÃO MÉDIA POR SPIN E POR DESORDEM VS CAMPO MAGNÉTICO EXTERNO

      H = H + 0.1

   END DO

END PROGRAM SpinRandomField_N3

SUBROUTINE Energy(S, x, dim_num, Ns, H, Energia)

   IMPLICIT NONE

   INTEGER :: Ns, dim_num, j
   REAL(8) :: x(dim_num), Energia(2**Ns), S(Ns, 2**Ns), H

   DO j = 1, 2**Ns  ! LOOP DAS CONFIGURAÇÕES DE SPIN

      Energia(j) = SUM(S(1:Ns-1, j)*S(2:Ns, j)) &
         + S(2, j)*S(4, j) + S(2, j)*S(6, j) + S(4, j)*S(6, j) + S(6, j)*S(1, j) &
         - H * SUM(S(:, j)) - x(1)*S(1, j) - x(2)*S(2, j) - x(3)*S(3, j) - x(4)*S(4, j) - x(5)*S(5, j) - x(6)*S(6, j)

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

SUBROUTINE r8vec_normal_01 (n, seed, x)

   implicit none

   integer ( kind = 4 ) n
   integer ( kind = 4 ) m
   integer ( kind = 4 ), save :: made = 0
   real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
   real ( kind = 8 ) r(n+1)
   real ( kind = 8 ) r8_uniform_01
   integer ( kind = 4 ), save :: saved = 0
   integer ( kind = 4 ) seed
   real ( kind = 8 ) x(n)
   integer ( kind = 4 ) x_hi_index
   integer ( kind = 4 ) x_lo_index
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
   !  Record the range of X we need to fill in.
   !
   x_lo_index = 1
   x_hi_index = n
   !
   !  Use up the old value, if we have it.
   !
   if ( saved == 1 ) then
      x(1) = y
      saved = 0
      x_lo_index = 2
   end if
   !
   !  Maybe we don't need any more values.
   !
   if ( x_hi_index - x_lo_index + 1 == 0 ) then
      !
      !  If we need just one new value, do that here to avoid null arrays.
      !
   else if ( x_hi_index - x_lo_index + 1 == 1 ) then

      r(1) = r8_uniform_01 ( seed )

      if ( r(1) == 0.0D+00 ) then
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'R8VEC_NORMAL_01 - Fatal error!'
         write ( *, '(a)' ) '  R8_UNIFORM_01 returned a value of 0.'
         stop
      end if

      r(2) = r8_uniform_01 ( seed )

      x(x_hi_index) = &
         sqrt ( -2.0D+00 * log ( r(1) ) ) * cos ( 2.0D+00 * pi * r(2) )
      y =      sqrt ( -2.0D+00 * log ( r(1) ) ) * sin ( 2.0D+00 * pi * r(2) )

      saved = 1

      made = made + 2
      !
      !  If we require an even number of values, that's easy.
      !
   else if ( mod ( x_hi_index - x_lo_index + 1, 2 ) == 0 ) then

      m = ( x_hi_index - x_lo_index + 1 ) / 2

      call r8vec_uniform_01 ( 2 * m, seed, r )

      x(x_lo_index:x_hi_index-1:2) = &
         sqrt ( -2.0D+00 * log ( r(1:2*m-1:2) ) ) &
         * cos ( 2.0D+00 * pi * r(2:2*m:2) )

      x(x_lo_index+1:x_hi_index:2) = &
         sqrt ( -2.0D+00 * log ( r(1:2*m-1:2) ) ) &
         * sin ( 2.0D+00 * pi * r(2:2*m:2) )

      made = made + x_hi_index - x_lo_index + 1

   else

      x_hi_index = x_hi_index - 1

      m = ( x_hi_index - x_lo_index + 1 ) / 2 + 1

      call r8vec_uniform_01 ( 2 * m, seed, r )

      x(x_lo_index:x_hi_index-1:2) = &
         sqrt ( -2.0D+00 * log ( r(1:2*m-3:2) ) ) &
         * cos ( 2.0D+00 * pi * r(2:2*m-2:2) )

      x(x_lo_index+1:x_hi_index:2) = &
         sqrt ( -2.0D+00 * log ( r(1:2*m-3:2) ) ) &
         * sin ( 2.0D+00 * pi * r(2:2*m-2:2) )

      x(n) = sqrt ( -2.0D+00 * log ( r(2*m-1) ) ) &
         * cos ( 2.0D+00 * pi * r(2*m) )

      y = sqrt ( -2.0D+00 * log ( r(2*m-1) ) ) &
         * sin ( 2.0D+00 * pi * r(2*m) )

      saved = 1

      made = made + x_hi_index - x_lo_index + 2

   end if

   return
END SUBROUTINE r8vec_normal_01

FUNCTION r8_uniform_01 (seed)


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
   !  Although SEED can be represented exactly as a 32 bit integer,
   !  it generally cannot be represented exactly as a 32 bit real number!
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

