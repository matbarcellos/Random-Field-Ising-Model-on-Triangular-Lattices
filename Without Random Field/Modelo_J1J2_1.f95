program energias_cluster
!#####################################################################################################################################
!######################################################### VARIÁVEIS GLOBAIS #########################################################
!#####################################################################################################################################
 
    implicit none

    integer, parameter :: N = 4 ! Número de sítios -- rede quadrada
    integer, dimension(2**N,N) :: site
    real(8), dimension(N) :: m
    real(8), dimension(2**N) :: H_total, Z
    real(8) :: Temp, J1, J2, Z_soma, m01, m02, m03, m04

    J1 =  -1.0d0
    J2 =  0.3d0
    Temp = 0.1d0

    ! Chute inicial Magnetizacao

    m01 =  1.0d0
    m02 =  1.0d0
    m03 =  1.0d0
    m04 =  1.0d0

!#####################################################################################################################################
!######################################################## PROGRAMA PRINCIPAL #########################################################
!#####################################################################################################################################

! Primeira Parte do Diagrama 0.0 até 0.5

    !open(unit=4, file='dados_J2_0.0.dat', status='unknown')
    !open(unit=5, file='dados_J2_0.1.dat', status='unknown')
    !open(unit=6, file='dados_J2_0.2.dat', status='unknown')
    !open(unit=7, file='dados_J2_0.3.dat', status='unknown')
    !open(unit=8, file='dados_J2_0.4.dat', status='unknown')
    !open(unit=9, file='dados_J2_0.5.dat', status='unknown')

    !open(unit=90, file='dados_J2_0.5.25.dat', status='unknown')
    !open(unit=91, file='dados_J2_0.5.50.dat', status='unknown')
    !open(unit=92, file='dados_J2_0.5.75.dat', status='unknown')

    !open(unit=93, file='dados_J2_0.6.25.dat', status='unknown')
    !open(unit=94, file='dados_J2_0.6.50.dat', status='unknown')
    !open(unit=95, file='dados_J2_0.6.75.dat', status='unknown')

! Segunda Parte do Diagrama 0.6 até 1.0

    !open(unit=10, file='dados_J2_0.6.dat', status='unknown')
    !open(unit=11, file='dados_J2_0.6_m0.dat', status='unknown')
    !open(unit=12, file='dados_J2_0.7.dat', status='unknown')
    !open(unit=13, file='dados_J2_0.7_m0.dat', status='unknown')
    !open(unit=14, file='dados_J2_0.8.dat', status='unknown')
    !open(unit=15, file='dados_J2_0.8_m0.dat', status='unknown')
    !open(unit=16, file='dados_J2_0.9.dat', status='unknown')
    !open(unit=17, file='dados_J2_0.9_m0.dat', status='unknown')
    !open(unit=18, file='dados_J2_1.0.dat', status='unknown')
    !open(unit=19, file='dados_J2_1.0_m0.dat', status='unknown')

    !open(unit=20, file='diagrama_0.0-0.5.dat', status='unknown') ! Arquivo do Diagrama 0.0 - 0.5

! Transições Primeira Ordem: do Diagrama 0.5 até 0.7

    !open(unit=22, file='dados_T0.1_J2_0.5-0.7.dat', status='unknown')
    open(unit=23, file='dados_T0.1_J2_0.5-0.7_m1.dat', status='unknown')
    !open(unit=24, file='dados_T0.2_J2_0.5-0.7.dat', status='unknown')
    !open(unit=25, file='dados_T0.2_J2_0.5-0.7_m1.dat', status='unknown')
    !open(unit=26, file='dados_T0.3_J2_0.5-0.7.dat', status='unknown')
    !open(unit=27, file='dados_T0.3_J2_0.5-0.7_m1.dat', status='unknown')
    !open(unit=28, file='dados_T0.4_J2_0.5-0.7.dat', status='unknown')
    !open(unit=29, file='dados_T0.4_J2_0.5-0.7_m1.dat', status='unknown')
    !open(unit=30, file='dados_T0.5_J2_0.5-0.7.dat', status='unknown')
    !open(unit=31, file='dados_T0.5_J2_0.5-0.7_m1.dat', status='unknown')
    !open(unit=32, file='dados_T0.6_J2_0.5-0.7.dat', status='unknown')
    !open(unit=33, file='dados_T0.6_J2_0.5-0.7_m1.dat', status='unknown')
    !open(unit=34, file='dados_T0.7_J2_0.5-0.7.dat', status='unknown')
    !open(unit=35, file='dados_T0.7_J2_0.5-0.7_m1.dat', status='unknown')
    !open(unit=36, file='dados_T0.8_J2_0.5-0.7.dat', status='unknown')
    !open(unit=37, file='dados_T0.8_J2_0.5-0.7_m1.dat', status='unknown')
    !open(unit=38, file='dados_T0.9_J2_0.5-0.7.dat', status='unknown')
    !open(unit=39, file='dados_T0.9_J2_0.5-0.7_m1.dat', status='unknown')
    !open(unit=40, file='dados_T1.0_J2_0.5-0.7.dat', status='unknown')
    !open(unit=41, file='dados_T1.0_J2_0.5-0.7_m1.dat', status='unknown')
    !open(unit=42, file='dados_T1.1_J2_0.5-0.7.dat', status='unknown')
    !open(unit=43, file='dados_T1.1_J2_0.5-0.7_m1.dat', status='unknown')
    !open(unit=44, file='dados_T1.2_J2_0.5-0.7.dat', status='unknown')
    !open(unit=45, file='dados_T1.2_J2_0.5-0.7_m1.dat', status='unknown')

    call sub_base_spins(N, site) ! montando a base de spins

    do while (abs(J2) <= abs(0.7d0))
        
        !Temp = 0.1d0
        !print *, abs(J2)
        !do while (Temp <= 4.0d0)
            
            call sub_hamiltoniano(N, J1, J2, site, m, H_total)
            call sub_Z(N, H_total, Temp, Z, Z_soma)
            call sub_mag_MPF(N, site, J1, J2, Temp, m, m01, m02, m03, m04)

            !if (abs(m(1)) <= 0.01d0) then
                !print *, J2, Temp
                !write (20,*) J2, Temp ! Diagrama no intervalo: 0.5 - 1.0
                !exit
            !end if

            write (23,*) abs(J2), m(1), -Temp*log(Z_soma) 
            !write (12,*) Temp, m(1), -Temp*log(Z_soma)

            !Temp = Temp + 0.01d0 ! Passo da Temperatura

        !end do

        J2 = J2 + 0.01d0 ! Passo de J2
        
    end do

    !call system('gnuplot -p Mag_vs_Temp.plt')

    close(23)   

end program energias_cluster

!#####################################################################################################################################
!################################### DECLARAÇÃO DAS SUBROTINAS E FUNÇÕES PARA O PROGRAMA PRINCIPAL ###################################
!#####################################################################################################################################

    ! MONTANDO A BASE DE SPINS - LISTA DE TUPLAS (SPINS)

    subroutine sub_base_spins(N, site)

        implicit none

        integer :: N, i, spin, pos
        integer, dimension(N) :: s
        integer, dimension(2**N,N) :: site

        do i = 0,(2**N)-1 ! Linha
            do pos = 0,N-1 ! Coluna
                spin = -1
                ! btest(i, pos) O resultado é verdadeiro se o bit pos de i tiver o valor 1. O resultado é falso se pos tiver o valor zero.
                
                if (btest(i, pos)) spin = 1
                
                s(pos+1) = spin
            end do
            site(i+1,1:N) = s ! site(1:16, 1:4) --> site(linhas, colunas)
        end do

    end subroutine sub_base_spins

    ! CÁLCULO EXATO INTRA-CLUSTER + APROXIMAÇÃO POR CAMPO MÉDIO INTER-CLUSTER = H_total
    
    subroutine sub_hamiltoniano(N, J1, J2, site, m, H_total)

        implicit none

        integer :: N, i
        integer, dimension(2**N,N) :: site
        real(8), dimension(2**N) :: H_total
        real(8), dimension(N) :: m
        real(8) :: J1, J2, H_intra, H_inter
        
        do i = 1,2**N

            H_intra = +J1 * ((site(i,1)*site(i,2)) + (site(i,2)*site(i,3)) + (site(i,3)*site(i,4)) + (site(i,4)*site(i,1))) & ! INTRA-CLUSTER (EXATA)
                      +J2 * ((site(i,1)*site(i,3)) + (site(i,2)*site(i,4)))

            H_inter = +J1 * ((site(i,1) * m(2)) - (m(1) * (m(2) / 2.0)) + (site(i,1) * m(4)) - (m(1) * (m(4) / 2.0)))  &  ! INTER-CLUSTER
                      +J1 * ((site(i,2) * m(3)) - (m(2) * (m(3) / 2.0)) + (site(i,2) * m(1)) - (m(2) * (m(1) / 2.0)))  &
                      +J1 * ((site(i,3) * m(4)) - (m(3) * (m(4) / 2.0)) + (site(i,3) * m(2)) - (m(3) * (m(2) / 2.0)))  &
                      +J1 * ((site(i,4) * m(1)) - (m(4) * (m(1) / 2.0)) + (site(i,4) * m(3)) - (m(4) * (m(3) / 2.0)))  &
                      +J2 * 3 * (site(i,1) * m(3)) - (m(1) * (m(3) / 2.0)) &   ! -> J2
                      +J2 * 3 * (site(i,2) * m(4)) - (m(2) * (m(4) / 2.0)) &
                      +J2 * 3 * (site(i,3) * m(1)) - (m(3) * (m(1) / 2.0)) &
                      +J2 * 3 * (site(i,4) * m(2)) - (m(4) * (m(2) / 2.0))

            H_total(i) = H_intra + H_inter

        end do
        
    end subroutine sub_hamiltoniano

    ! FUNÇÃO DE PARTIÇÃO

    subroutine sub_Z(N, H_total, Temp, Z, Z_soma)

        implicit none

        integer:: N, i
        real(8), dimension(2**N) :: H_total, Z
        real(8) :: Temp, Z_soma

    
        Z_soma = 0.0d0

        do i = 1,2**N

            Z(i) = exp((-H_total(i))/Temp) ! esta lista será usada para o cálculo da magnetização
            Z_soma = Z_soma + Z(i)

        end do
    
    end subroutine sub_Z

    ! MAGNETIZAÇÃO LOCAL m1, m2, m3, m4 -- COM MÉTODO DO PONTO FIXO

    subroutine sub_mag_MPF(N, site, J1, J2, Temp, m, m01, m02, m03, m04)

        implicit none

        integer, dimension(2**N,N) :: site
        integer :: N, i, j, iter
        real(8), dimension(2**N) :: H_total, Z
        real(8), dimension(N) :: m, mag_teste
        real(8) :: tol, erro, Temp, J1, J2, Z_soma, m01, m02, m03, m04

        m = [m01, m02, m03, m04]
        erro = 1.0 ! Atribuí esse valor para que se inicie o "do while" a seguir
        tol = 1.0d-8 
        iter = 0

        do while(erro >= tol .and. iter<= 1000) ! 
            call sub_hamiltoniano(N, J1, J2, site, m, H_total)
            call sub_Z(N, H_total, Temp, Z, Z_soma)
            mag_teste = .0d0

            do i = 1,N
                
                do j = 1,2**N
                    
                    mag_teste(i) = mag_teste(i) + (site(j,i) * (Z(j)/Z_soma))
                    
                end do

            end do

            erro = maxval(abs(mag_teste - m))
            m = mag_teste
            iter = iter + 1
            !print *, m
        end do
    end subroutine sub_mag_MPF