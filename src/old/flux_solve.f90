!flux 2d flow solve module 
!max wood
!version : 0.0.1
!updated : 28-03-25

!module 
module flux_solve
use omp_lib
use edge_flux
use io_utilities
use flux_data_methods
contains 


!flow initialise ===============
subroutine flux_flow_initialise(mesh,options)
implicit none 

!variables - inout
type(flux_mesh) :: mesh 
type(flux_options) :: options 

!variables - local 
integer(in64) :: ii 

!set angle of atack in radians 
options%aoarad = options%aoadeg*(pi/180.0d0)

!set scaled density and pressure from the input values 
options%rhoinf = options%rhoinf/1.225d0
options%tinf = (options%tinf/273.15)*(1.0d0/(options%gamma*options%R))

!set initial flow variables at normalised scale
options%cinf = sqrt(options%gamma*options%R*options%tinf)
options%pinf = (options%cinf*options%cinf*options%rhoinf)/options%gamma
options%machinf = options%machinf 
options%velinf = options%machinf*options%cinf

!set freestram velocity components 
options%uinf = options%velinf*cos(options%aoarad)
options%vinf = options%velinf*sin(options%aoarad)

!set the freestream total conditions 
options%p0inf = total_pressure(options%pinf,options%machinf,options%gamma)
options%rho0inf = total_density(options%rhoinf,options%machinf,options%gamma)
options%t0inf = total_temperature(options%tinf,options%machinf,options%gamma)

!set the primative conditions in each cell 
do ii=1,mesh%ncell
    mesh%cells(ii)%rho = options%rhoinf
    mesh%cells(ii)%u = options%uinf
    mesh%cells(ii)%v = options%vinf
    mesh%cells(ii)%p = options%pinf
    mesh%cells(ii)%mach = options%machinf
    mesh%cells(ii)%e = energy(mesh%cells(ii)%p,mesh%cells(ii)%rho,options%velinf,options%gamma)
    mesh%cells(ii)%cp = pressure_coefficient(mesh%cells(ii)%p,options)
end do 

!set the conservative variables in each cell 
do ii=1,mesh%ncell
    call mesh%cells(ii)%prim2con(options%gamma)
end do 

!display the flow properties 
if (options%cdisplay) then
    ! write(*,'(A)') '    Actual freestream flow properties: '
    ! write(*,'(A,A,A)') '    {pressure (Pa): ',real2F0_Xstring(pinf_AC,6_in),'}'
    ! write(*,'(A,A,A)') '    {density (Kg/m^3): ',real2F0_Xstring(options%rhoinf,6_in),'}' 
    ! write(*,'(A,A,A)') '    {speed of sound (m/s): ',real2F0_Xstring(sosinf_AC,6_in),'}' 
    ! write(*,'(A,A,A)') '    {velocity (m/s): ',real2F0_Xstring(options%machinf*sosinf_AC,6_in),'}' 
    write(*,'(A)') '    scaled freestream flow properties: '
    write(*,'(A,A,A)') '    {pressure : ',real2F0_Xstring(options%pinf,6_in64),'}' 
    write(*,'(A,A,A)') '    {density : ',real2F0_Xstring(options%rhoinf,6_in64),'}'
    write(*,'(A,A,A)') '    {speed of sound : ',real2F0_Xstring(options%cinf,6_in64),'}'
    write(*,'(A,A,A)') '    {velocity : ',real2F0_Xstring(options%velinf,6_in64),'}'
end if 
return 
end subroutine flux_flow_initialise

!flow solve subroutine ===============
subroutine flux_flow_solve(mesh,options)
implicit none 

!variables - inout
type(flux_mesh) :: mesh 
type(flux_options) :: options 

!variables - local 
logical :: resconv,nanflag
integer(in64) :: fsitt,ii 
real(dp) :: rhores

!set the number of threads 
call omp_set_num_threads(options%num_threads)

!initialise flags
resconv = .false.
nanflag = .false.

!initialise the parallel region 
!$OMP parallel

!solve 
do fsitt=1,options%niter_max

    !convergence cycle condition
    if (resconv) then
        cycle
    end if

    !nan value cycle condition 
    if (nanflag) then 
        cycle 
    end if


    !evaluate cell spectral radii and timesteps 
    !$OMP do schedule(static)
    do ii=1,mesh%ncell
        call get_cell_spectral_radius(mesh%cells(ii),mesh,options)
        call get_cell_timestep(mesh%cells(ii),options)
    end do 
    !$OMP end do 

    !rk itterate for this timestep 
    call rk_itterate(mesh,options,nanflag)

    !process this step 
    !$OMP single
    rhores = log10(sqrt(sum((mesh%cells(:)%r1 + mesh%cells(:)%d1)**2)))
    !if (mod(fsitt,100) == 0) then
        print *, 'iter = ',fsitt,' rho_res = ', rhores
    !end if 
    if (rhores .LT. options%residual_convtol) then 
        resconv = .true.
    end if 

    !$OMP end single
end do 

!end the parallel region
!$OMP end parallel 
return 
end subroutine flux_flow_solve

!rk itterate ===============
subroutine rk_itterate(mesh,options,nanflag)
implicit none 

!variables - inout
logical :: nanflag
type(flux_mesh) :: mesh 
type(flux_options) :: options 

!variables - local 
integer(in64) :: ii,rr 
real(dp) :: cell_step


!$OMP single

!store the initial state of each cell 
mesh%cells(:)%w10 = mesh%cells(:)%w1
mesh%cells(:)%w20 = mesh%cells(:)%w2
mesh%cells(:)%w30 = mesh%cells(:)%w3
mesh%cells(:)%w40 = mesh%cells(:)%w4

! !initialise the dissipation 
! mesh%cells(:)%d1 = 0.0d0 
! mesh%cells(:)%d2 = 0.0d0 
! mesh%cells(:)%d3 = 0.0d0 
! mesh%cells(:)%d4 = 0.0d0 


!$OMP end single

!rk itterate
do rr=1,options%rk_niter

    !evaluate residual in each cell 
    !$OMP do schedule(static)
    do ii=1,mesh%ncell
        call cell_flux(mesh%cells(ii)%f1,mesh%cells(ii)%f2,mesh%cells(ii)%f3,mesh%cells(ii)%f4,&
                       mesh%cells(ii)%g1,mesh%cells(ii)%g2,mesh%cells(ii)%g3,mesh%cells(ii)%g4,mesh%cells(ii)) 
    end do 
    !$OMP end do 
    !$OMP do schedule(static)
    do ii=1,mesh%ncell
        call get_cell_residual(mesh%cells(ii),mesh,options)
    end do 
    !$OMP end do 

    !evaluate dissipation in each cell if required 
    if (options%rk_dissipation(rr)) then 
        !$OMP do schedule(static)
        do ii=1,mesh%ncell
            call get_cell_pressure_sensor(mesh%cells(ii),mesh)
            call get_cell_laplacian(mesh%cells(ii),mesh)
        end do 
        !$OMP end do 
        !$OMP do schedule(static)
        do ii=1,mesh%ncell
            call get_cell_dissipation(mesh%cells(ii),mesh,options)
        end do 
        !$OMP end do 
    end if 

    !step each cell and update the primitive variables in each cell 
    !$OMP do schedule(static) private(cell_step)
    do ii=1,mesh%ncell
        cell_step = rk4_alpha(rr)*(mesh%cells(ii)%timestep/mesh%cells(ii)%volume)
        mesh%cells(ii)%w1 = mesh%cells(ii)%w10 - cell_step*(mesh%cells(ii)%r1 + mesh%cells(ii)%d1)
        mesh%cells(ii)%w2 = mesh%cells(ii)%w20 - cell_step*(mesh%cells(ii)%r2 + mesh%cells(ii)%d2)
        mesh%cells(ii)%w3 = mesh%cells(ii)%w30 - cell_step*(mesh%cells(ii)%r3 + mesh%cells(ii)%d3)
        mesh%cells(ii)%w4 = mesh%cells(ii)%w40 - cell_step*(mesh%cells(ii)%r4 + mesh%cells(ii)%d4)
        call mesh%cells(ii)%con2prim(options%gamma,options)
        if (isnan(mesh%cells(ii)%w1)) then 
            print *, 'cell nan: ',ii
            nanflag = .true.
        end if 
    end do
    !$OMP end do 
end do 
return 
end subroutine rk_itterate

!get cell residual ===============
subroutine get_cell_residual(cell,mesh,options)
implicit none 

!variables - inout
type(flux_cell) :: cell
type(flux_mesh) :: mesh 
type(flux_options) :: options 
    
!variables - local 
integer(in64) :: ii 
real(dp) :: vnorm,normmach
real(dp) :: f1c,f2c,f3c,f4c,g1c,g2c,g3c,g4c
real(dp) :: f1a,f2a,f3a,f4a,g1a,g2a,g3a,g4a
real(dp) :: flux1,flux2,flux3,flux4

!accumulate across edges 
cell%r1 = 0.0d0 
cell%r2 = 0.0d0 
cell%r3 = 0.0d0 
cell%r4 = 0.0d0 
do ii=1,cell%nedge

    !initialise
    f1c = cell%f1
    f2c = cell%f2
    f3c = cell%f3
    f4c = cell%f4
    g1c = cell%g1
    g2c = cell%g2
    g3c = cell%g3
    g4c = cell%g4

    !evaluate edge flux
    if (cell%edges(ii,3) .GT. 0) then !internal cell 
        f1a = mesh%cells(cell%edges(ii,3))%f1
        f2a = mesh%cells(cell%edges(ii,3))%f2
        f3a = mesh%cells(cell%edges(ii,3))%f3
        f4a = mesh%cells(cell%edges(ii,3))%f4
        g1a = mesh%cells(cell%edges(ii,3))%g1
        g2a = mesh%cells(cell%edges(ii,3))%g2
        g3a = mesh%cells(cell%edges(ii,3))%g3
        g4a = mesh%cells(cell%edges(ii,3))%g4
    else !boundary condition 
        if (cell%edges(ii,3) == -1) then !wall 
            call wall_flux(f1a,f2a,f3a,f4a,g1a,g2a,g3a,g4a,cell)
            f1c = f1a
            f2c = f2a
            f3c = f3a
            f4c = f4a
            g1c = g1a
            g2c = g2a
            g3c = g3a
            g4c = g4a
        elseif (cell%edges(ii,3) == -2) then !far field
            vnorm = (cell%u*cell%edges_nx(ii) + cell%v*cell%edges_ny(ii))
            normmach = abs(vnorm)/speed_of_sound(cell%p,cell%rho,options%gamma)
            if (normmach .GE. 1.0d0) then !supersonic 
                if (vnorm .LE. 0.0d0) then !inflow
                    call supersonic_test(f1a,f2a,f3a,f4a,g1a,g2a,g3a,g4a,cell,options,1_in64)
                else !outflow
                    call supersonic_test(f1a,f2a,f3a,f4a,g1a,g2a,g3a,g4a,cell,options,-1_in64)
                end if 
            else !subsonic  
                if (vnorm .LE. 0.0d0) then !inflow
                    call far_field_subsonic_flux_characteristic(f1a,f2a,f3a,f4a,g1a,g2a,g3a,g4a,cell,ii,options,1_in64)
                else !outflow
                    call far_field_subsonic_flux_characteristic(f1a,f2a,f3a,f4a,g1a,g2a,g3a,g4a,cell,ii,options,-1_in64)
                end if 
            end if 


            ! if (vnorm .LE. 0.0d0) then !inflow
            !     ! call far_field_subsonic_flux_characteristic(f1a,f2a,f3a,f4a,g1a,g2a,g3a,g4a,cell,ii,options,1_in64)

            !     call supersonic_test(f1a,f2a,f3a,f4a,g1a,g2a,g3a,g4a,cell,options,1_in64)



            ! elseif (vnorm .GT. 0.0d0) then !outflow
            !     ! call far_field_subsonic_flux_characteristic(f1a,f2a,f3a,f4a,g1a,g2a,g3a,g4a,cell,ii,options,-1_in64)

            !     call supersonic_test(f1a,f2a,f3a,f4a,g1a,g2a,g3a,g4a,cell,options,-1_in64)

            ! end if


        elseif (cell%edges(ii,3) == -3) then !inflow - stagnation conditions 

        elseif (cell%edges(ii,3) == -4) then !inflow - freestream 

        elseif (cell%edges(ii,3) == -5) then !outflow - backpressure

        end if 
    end if 
    flux1 = 0.5d0*((f1c + f1a)*cell%edges_dy(ii) - (g1c + g1a)*cell%edges_dx(ii))
    flux2 = 0.5d0*((f2c + f2a)*cell%edges_dy(ii) - (g2c + g2a)*cell%edges_dx(ii))
    flux3 = 0.5d0*((f3c + f3a)*cell%edges_dy(ii) - (g3c + g3a)*cell%edges_dx(ii))
    flux4 = 0.5d0*((f4c + f4a)*cell%edges_dy(ii) - (g4c + g4a)*cell%edges_dx(ii))
    
    !accumulate to cell 
    cell%r1 = cell%r1 + flux1
    cell%r2 = cell%r2 + flux2
    cell%r3 = cell%r3 + flux3
    cell%r4 = cell%r4 + flux4
end do 



return 
end subroutine get_cell_residual

!get cell dissipation ===============
subroutine get_cell_dissipation(cell,mesh,options)
implicit none 

!variables - inout
type(flux_cell) :: cell
type(flux_mesh) :: mesh 
type(flux_options) :: options 

!variables - local 
integer(in64) :: ii 
real(dp) :: psi0,psi1,psi01,s2,s4,dk2,dk4,edge_specrad

!evaluate 
cell%d1 = 0.0d0 
cell%d2 = 0.0d0 
cell%d3 = 0.0d0 
cell%d4 = 0.0d0 

do ii=1,cell%nedge
    if (cell%edges(ii,3) .GT. 0) then 

    !cell edge quantity scaling values *******************
    s2 = 1.0 
    s4 = 0.25

    edge_specrad = cell%edges_specrad(ii)
    ! edge_specrad = max(cell%specrad,mesh%cells(cell%edges(ii,3))%specrad)


    !4th order psi coefficient 
    psi0 = sqrt(cell%specrad/(4.0d0*edge_specrad))
    psi1 = sqrt(mesh%cells(cell%edges(ii,3))%specrad/(4.0d0*edge_specrad))
    psi01 = (4.0d0*(psi0*psi1))/(psi0 + psi1)

    !coefficients
    ! dk2 = 0.5d0*(cell%psens + mesh%cells(cell%edges(ii,3))%psens)*s2*options%k2
    ! dk4 = max(0.0d0,(options%k4 - dK2))*S4    


    

    ! dk2 = max(cell%psens,mesh%cells(cell%edges(ii,3))%psens)*s2*options%k2*edge_specrad
    dk2 = 0.5d0*(cell%psens + mesh%cells(cell%edges(ii,3))%psens)*s2*options%k2*edge_specrad
    dk4 = max(0.0d0,(options%k4*edge_specrad - 2.0d0*dK2))*S4
    ! dpsilam = psi01


    ! dk2 = max(cell%psens,mesh%cells(cell%edges(ii,3))%psens)*s2*options%k2
    ! dk4 = max(0.0d0,(options%k4*cell%edges_specrad(ii) - 2.0d0*dK2))*S4
    ! dpsilam = cell%edges_specrad(ii)*psi01

    !accumulate dissipation for the edge 
    cell%d1 = cell%d1 + (dk2*(cell%w1 - mesh%cells(cell%edges(ii,3))%w1) + &
    dk4*(cell%lw1 - mesh%cells(cell%edges(ii,3))%lw1))*psi01
    cell%d2 = cell%d2 + (dk2*(cell%w2 - mesh%cells(cell%edges(ii,3))%w2) + &
    dk4*(cell%lw2 - mesh%cells(cell%edges(ii,3))%lw2))*psi01
    cell%d3 = cell%d3 + (dk2*(cell%w3 - mesh%cells(cell%edges(ii,3))%w3) + &
    dk4*(cell%lw3 - mesh%cells(cell%edges(ii,3))%lw3))*psi01
    cell%d4 = cell%d4 + (dk2*(cell%w4 - mesh%cells(cell%edges(ii,3))%w4) + &
    dk4*(cell%lw4 - mesh%cells(cell%edges(ii,3))%lw4))*psi01
    end if 
end do 

! do ii=1,cell%nedge
!     if (cell%edges(ii,3) .GT. 0) then 

!     !cell edge quantity scaling values *******************
!     s2 = 1.0 
!     s4 = 0.25

!     !4th order psi coefficient 
!     psi0 = sqrt(cell%specrad/(4.0d0*cell%edges_specrad(ii)))
!     psi1 = sqrt(mesh%cells(cell%edges(ii,3))%specrad/(4.0d0*cell%edges_specrad(ii)))
!     psi01 = (4.0d0*(psi0*psi1))/(psi0 + psi1)

!     !coefficients
!     ! dk2 = 0.5d0*(cell%psens + mesh%cells(cell%edges(ii,3))%psens)*s2*options%k2
!     ! dk4 = max(0.0d0,(options%k4 - dK2))*S4    


!     edge_specrad = max(cell%specrad,mesh%cells(cell%edges(ii,3))%specrad)

!     dk2 = max(cell%psens,mesh%cells(cell%edges(ii,3))%psens)*s2*options%k2
!     dk4 = max(0.0d0,(options%k4*edge_specrad - 2.0d0*dK2))*S4
!     dpsilam = edge_specrad*psi01


!     ! dk2 = max(cell%psens,mesh%cells(cell%edges(ii,3))%psens)*s2*options%k2
!     ! dk4 = max(0.0d0,(options%k4*cell%edges_specrad(ii) - 2.0d0*dK2))*S4
!     ! dpsilam = cell%edges_specrad(ii)*psi01

!     !accumulate dissipation for the edge 
!     cell%d1 = cell%d1 + (dk2*(cell%w1 - mesh%cells(cell%edges(ii,3))%w1) + &
!     dk4*(cell%lw1 - mesh%cells(cell%edges(ii,3))%lw1))*dpsilam
!     cell%d2 = cell%d2 + (dk2*(cell%w2 - mesh%cells(cell%edges(ii,3))%w2) + &
!     dk4*(cell%lw2 - mesh%cells(cell%edges(ii,3))%lw2))*dpsilam
!     cell%d3 = cell%d3 + (dk2*(cell%w3 - mesh%cells(cell%edges(ii,3))%w3) + &
!     dk4*(cell%lw3 - mesh%cells(cell%edges(ii,3))%lw3))*dpsilam
!     cell%d4 = cell%d4 + (dk2*(cell%w4 - mesh%cells(cell%edges(ii,3))%w4) + &
!     dk4*(cell%lw4 - mesh%cells(cell%edges(ii,3))%lw4))*dpsilam
!     end if 
! end do 
return 
end subroutine get_cell_dissipation

!get cell laplacian ===============
subroutine get_cell_laplacian(cell,mesh)
implicit none 

!variables - inout
type(flux_cell) :: cell
type(flux_mesh) :: mesh 

!variables - local 
integer(in64) :: ii 

!evaluate 
cell%lw1 = 0.0d0 
cell%lw2 = 0.0d0 
cell%lw3 = 0.0d0 
cell%lw4 = 0.0d0 
do ii=1,cell%nedge
    if (cell%edges(ii,3) .GT. 0) then 
        cell%lw1 = cell%lw1 + cell%w1 - mesh%cells(cell%edges(ii,3))%w1
        cell%lw2 = cell%lw2 + cell%w2 - mesh%cells(cell%edges(ii,3))%w2
        cell%lw3 = cell%lw3 + cell%w3 - mesh%cells(cell%edges(ii,3))%w3
        cell%lw4 = cell%lw4 + cell%w4 - mesh%cells(cell%edges(ii,3))%w4
    end if 
end do 
return 
end subroutine get_cell_laplacian

!get cell pressure sensor ===============
subroutine get_cell_pressure_sensor(cell,mesh)
implicit none 

!variables - inout
type(flux_cell) :: cell
type(flux_mesh) :: mesh 

!variables - local 
integer(in64) :: ii 
real(dp) :: pn,pdn 

!evaluate 
pn = 0.0d0 
pdn = 0.0d0 
do ii=1,cell%nedge
    if (cell%edges(ii,3) .GT. 0) then 
        pn = pn + (cell%p - mesh%cells(cell%edges(ii,3))%p)
        pdn = pdn + (cell%p + mesh%cells(cell%edges(ii,3))%p)
    end if 
end do 
cell%psens = abs(pn/pdn)
return 
end subroutine get_cell_pressure_sensor

!evaluate cell timestep ===============
subroutine get_cell_timestep(cell,options)
implicit none 

!variables - inout
type(flux_cell) :: cell
type(flux_options) :: options 

!evaluate
cell%timestep = options%cfl*(cell%volume/cell%specrad)
return 
end subroutine get_cell_timestep

!evaluate cell spectral radius ===============
subroutine get_cell_spectral_radius(cell,mesh,options)
implicit none 

!variables - inout
type(flux_cell) :: cell
type(flux_mesh) :: mesh 
type(flux_options) :: options 

!variables - local 
integer(in64) :: ii 

!evaluate 
cell%specrad = 0.0d0 
do ii=1,cell%nedge
    cell%edges_specrad(ii) = edge_spectral_radius(cell,cell%edges(ii,3),ii,mesh,options)
end do 
cell%specrad = sum(cell%edges_specrad)
return 
end subroutine get_cell_spectral_radius

!edge spectral radius ===============
function edge_spectral_radius(cell,cell_adjacent_index,edgeidx,mesh,options) result(lam_edge)
implicit none 

!variables - inout
integer(in64) :: cell_adjacent_index,edgeidx
real(dp) :: lam_edge
type(flux_cell) :: cell
type(flux_mesh) :: mesh 
type(flux_options) :: options 

!variables - local 
real(dp) :: ue,ve,sose

!evaluate 
if (cell_adjacent_index .GT. 0) then !cell 
    ue = 0.5d0*(cell%u + mesh%cells(cell_adjacent_index)%u)
    ve = 0.5d0*(cell%v + mesh%cells(cell_adjacent_index)%v)
    sose = 0.5d0*(speed_of_sound(cell%p,cell%rho,options%gamma) + &
    speed_of_sound(mesh%cells(cell_adjacent_index)%p,mesh%cells(cell_adjacent_index)%rho,options%gamma))
else !boundary condition 
    ue = cell%u
    ve = cell%v
    sose = speed_of_sound(cell%p,cell%rho,options%gamma)
end if 

!calculate edge flow spectral radius
lam_edge = (abs(cell%edges_nx(edgeidx)*ue + cell%edges_ny(edgeidx)*ve) + sose)*cell%edges_len(edgeidx)
return 
end function edge_spectral_radius

end module flux_solve
