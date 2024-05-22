! integracao via Monte Carlos

PROGRAM integral

INTEGER :: nfunc, dim_num, eval_num, seed
REAL(8) :: result(1)

nfunc = 1
dim_num = 1
eval_num = 1000000
seed = 123456

call monte_carlo_multifuctions( nfunc, dim_num, eval_num, seed, result )

print*, result

END 

REAL(8) FUNCTION func(dim_num, x) 
INTEGER :: dim_num
REAL(8) :: x(dim_num)

func = 0D0
DO i = 1, dim_num
func = func + 3D0*dexp(-2D0*(3D0*x(i))**2)
END DO

END

  subroutine monte_carlo_multifuctions ( nfunc, dim_num, eval_num, seed, result )

      implicit none
    
      integer ( kind = 4 ) dim_num,nfunc,ii
      integer ( kind = 4 ) eval_num
      real ( kind = 8 ), external :: func
      integer ( kind = 4 ) i
      real ( kind = 8 ) result (nfunc)
      integer ( kind = 4 ) seed
      real ( kind = 8 ) x(dim_num)
    
      result = 0.0D+00
    
      do i = 1, eval_num
    
       call r8vec_normal_01 ( dim_num, seed, x ) ! -> x contendo números aleatórios com dimensão dim_num
       
       do ii= 1, nfunc
        result(ii) = result(ii) + func ( dim_num, x )
       enddo

      end do
    
      !volume = product ( b(1:dim_num) - a(1:dim_num) )
      !result = result * volume / real ( eval_num, kind = 8 )

      do i=1, nfunc
       result(i) = result(i)  / real ( eval_num, kind = 8 )
      enddo
    
      return
    end      
               
    
      ! integracao via Monte Carlos
      subroutine monte_carlo_nd ( func, dim_num, eval_num, seed, result )
    
      implicit none
    
      integer ( kind = 4 ) dim_num
    
      integer ( kind = 4 ) eval_num
      real ( kind = 8 ), external :: func
      integer ( kind = 4 ) i
      real ( kind = 8 ) result
      integer ( kind = 4 ) seed
      real ( kind = 8 ) x(dim_num)
    
      result = 0.0D+00
    
      do i = 1, eval_num
    
       call r8vec_normal_01 ( dim_num, seed, x )

       result = result + func ( dim_num, x )

      end do

      result = result  / real ( eval_num, kind = 8 )
    
      return
    end
    
    subroutine r8vec_normal_01 ( n, seed, x )
    
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
    end
    
    function r8_uniform_01 ( seed )
    
    
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
    end
         
    subroutine r8vec_uniform_01 ( n, seed, r )
    
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
    end