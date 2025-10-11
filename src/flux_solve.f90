!flux 2d flow solve module 
!max wood
!version : 0.0.2
!updated : 11-10-25

!module 
module flux_solve
! use omp_lib
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
integer(in32) :: ii 

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
options%p0inf = stagnation_pressure(options%pinf,options%machinf,options%gamma)
options%rho0inf = stagnation_density(options%rhoinf,options%machinf,options%gamma)
options%t0inf = stagnation_temperature(options%tinf,options%machinf,options%gamma)

!set the primative conditions in each cell 
do ii=1,mesh%ncell
    mesh%rho(ii) = options%rhoinf
    mesh%u(ii) = options%uinf
    mesh%v(ii) = options%vinf
    mesh%p(ii) = options%pinf
    mesh%mach(ii) = options%machinf
    mesh%e(ii) = energy(mesh%p(ii),mesh%rho(ii),options%velinf,options%gamma)
    mesh%cp(ii) = pressure_coefficient(mesh%p(ii),options)
end do 

!set the conservative variables in each cell 
do ii=1,mesh%ncell
    call prim2con(mesh%rho(ii),mesh%u(ii),mesh%v(ii),mesh%p(ii),options%gamma,mesh%w1(ii),mesh%w2(ii),mesh%w3(ii),mesh%w4(ii))
end do 

!display the flow properties 
if (options%cdisplay) then
    ! write(*,'(A)') '    actual freestream flow properties: '
    ! write(*,'(A,A,A)') '    {pressure (Pa): ',real2F0_Xstring(pinf_AC,6_in64),'}'
    ! write(*,'(A,A,A)') '    {density (Kg/m^3): ',real2F0_Xstring(options%rhoinf,6_in64),'}' 
    ! write(*,'(A,A,A)') '    {speed of sound (m/s): ',real2F0_Xstring(sosinf_AC,6_in64),'}' 
    ! write(*,'(A,A,A)') '    {velocity (m/s): ',real2F0_Xstring(options%machinf*sosinf_AC,6_in64),'}' 
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
integer(in32) :: ee,cc
integer(in32) :: iteration
logical :: resconv,nanflag
real(dp) :: rhores

!initialise flags
resconv = .false.
nanflag = .false.

!solve 
do iteration=1,options%niter_max

    !convergence cycle condition
    if (resconv) then
        cycle
    end if

    !nan value cycle condition 
    if (nanflag) then 
        cycle 
    end if

    !evaluate edge spectral radii
    do ee=1,mesh%nedge
        call edge_spectral_radius(mesh,options,ee)
    end do 

    !evaluate cell spectral radii and timesteps 
    do cc=1,mesh%ncell
        mesh%cells_specrad(cc) = 0.0d0 
        do ee=1,mesh%cells(cc)%nedge
            mesh%cells_specrad(cc) = mesh%cells_specrad(cc) + mesh%edges_specrad(mesh%cells(cc)%edge(ee))
        end do
        mesh%cells_dt(cc) = options%cfl*(mesh%cells_volume(cc)/mesh%cells_specrad(cc))
    end do 

    !rk iterate this timestep 
    call rk_iterate(mesh,options,nanflag)



    !display 
    rhores = log10(sqrt(sum((mesh%residual)**2)))

    print *, 'iter = ',iteration,' rho_res = ', rhores


end do 

return 
end subroutine flux_flow_solve

!rk iterate ===============
subroutine rk_iterate(mesh,options,nanflag)
implicit none 

!variables - inout
type(flux_mesh) :: mesh 
type(flux_options) :: options 
logical :: nanflag

!variables - local 
integer(in32) :: rr,ee,cc
real(dp) :: r1,r2,r3,r4,cell_stepsize

!store the initial state of each cell 
mesh%w10 = mesh%w1
mesh%w20 = mesh%w2
mesh%w30 = mesh%w3
mesh%w40 = mesh%w4

!set the edge dissipation initial values 


!iterate
do rr=1,options%rk_niter

    !evaluate edge fluxes
    call get_edge_fluxes(mesh,options)

    !evaluate edge dissipations 
    if (options%rk_dissipation(rr)) then 
        do ee=1,mesh%nedge
            call edge_pressure_sensor(mesh,ee) 
            call edge_laplacian(mesh,ee)
        end do
        do cc=1,mesh%ncell
            call cell_pressure_sensor(mesh,cc) 
            call cell_laplacian(mesh,cc) 
        end do 
        call get_edge_dissipations(mesh,options)
    end if  

    !timestep each cell 
    do cc=1,mesh%ncell
        cell_stepsize = rk4_alpha(rr)*(mesh%cells_dt(cc)/mesh%cells_volume(cc))
        r1 = 0.0d0 
        r2 = 0.0d0 
        r3 = 0.0d0 
        r4 = 0.0d0 
        do ee=1,mesh%cells(cc)%nedge
            r1 = r1 + (mesh%edges_r1(mesh%cells(cc)%edge(ee)) + mesh%edges_d1(mesh%cells(cc)%edge(ee)))*mesh%cells(cc)%edge_sign(ee)
            r2 = r2 + (mesh%edges_r2(mesh%cells(cc)%edge(ee)) + mesh%edges_d2(mesh%cells(cc)%edge(ee)))*mesh%cells(cc)%edge_sign(ee)
            r3 = r3 + (mesh%edges_r3(mesh%cells(cc)%edge(ee)) + mesh%edges_d3(mesh%cells(cc)%edge(ee)))*mesh%cells(cc)%edge_sign(ee)
            r4 = r4 + (mesh%edges_r4(mesh%cells(cc)%edge(ee)) + mesh%edges_d4(mesh%cells(cc)%edge(ee)))*mesh%cells(cc)%edge_sign(ee)
        end do 
        mesh%w1(cc) = mesh%w10(cc) - cell_stepsize*r1
        mesh%w2(cc) = mesh%w20(cc) - cell_stepsize*r2
        mesh%w3(cc) = mesh%w30(cc) - cell_stepsize*r3
        mesh%w4(cc) = mesh%w40(cc) - cell_stepsize*r4
        mesh%residual(cc) = r1
        call con2prim(mesh%rho(cc),mesh%u(cc),mesh%v(cc),mesh%p(cc),mesh%e(cc),options%gamma,mesh%w1(cc),mesh%w2(cc),mesh%w3(cc),mesh%w4(cc))
        if (isnan(mesh%w1(cc))) then 
            print *, 'cell nan: ',cc
            nanflag = .true.
            exit 
        end if 
    end do 
    if (nanflag) then 
        exit 
    end if 
end do 
return 
end subroutine rk_iterate

!get edge fluxes ===============
subroutine get_edge_fluxes(mesh,options)
implicit none 

!variables - inout
type(flux_mesh) :: mesh 
type(flux_options) :: options 

!variables - local 
integer(in32) :: ee 
real(dp) :: f1c,f2c,f3c,f4c,g1c,g2c,g3c,g4c
real(dp) :: f1a,f2a,f3a,f4a,g1a,g2a,g3a,g4a
real(dp) :: velnorm,machnorm

!evaluate each edge flux 
do ee=1,mesh%nedge
    if (mesh%edges(ee)%c2 .GT. 0) then !internal cell
        call cell_flux(f1c,f2c,f3c,f4c,g1c,g2c,g3c,g4c,mesh,mesh%edges(ee)%c1)
        call cell_flux(f1a,f2a,f3a,f4a,g1a,g2a,g3a,g4a,mesh,mesh%edges(ee)%c2)
    elseif (mesh%edges(ee)%c2 == -1) then !wall       
        call wall_flux(f1c,f2c,f3c,f4c,g1c,g2c,g3c,g4c,mesh,mesh%edges(ee)%c1)
        f1a = f1c
        f2a = f2c
        f3a = f3c
        f4a = f4c
        g1a = g1c
        g2a = g2c
        g3a = g3c
        g4a = g4c
    elseif (mesh%edges(ee)%c2 == -2) then !farfield
        velnorm = mesh%u(mesh%edges(ee)%c1)*mesh%edges(ee)%nx + mesh%v(mesh%edges(ee)%c1)*mesh%edges(ee)%ny
        machnorm = abs(velnorm)/speed_of_sound(mesh%p(mesh%edges(ee)%c1),mesh%rho(mesh%edges(ee)%c1),options%gamma)
        call cell_flux(f1c,f2c,f3c,f4c,g1c,g2c,g3c,g4c,mesh,mesh%edges(ee)%c1)
        if (velnorm .LT. 0.0d0) then !inflow
            call supersonic_farfield_flux(f1a,f2a,f3a,f4a,g1a,g2a,g3a,g4a,mesh,mesh%edges(ee)%c1,options,1)
        else !outflow
            call supersonic_farfield_flux(f1a,f2a,f3a,f4a,g1a,g2a,g3a,g4a,mesh,mesh%edges(ee)%c1,options,-1)
        end if 
    else 

    !todo: add other boundary conditions here

    end if 
    mesh%edges_r1(ee) = 0.5d0*((f1c + f1a)*mesh%edges(ee)%dy - (g1c + g1a)*mesh%edges(ee)%dx)
    mesh%edges_r2(ee) = 0.5d0*((f2c + f2a)*mesh%edges(ee)%dy - (g2c + g2a)*mesh%edges(ee)%dx)
    mesh%edges_r3(ee) = 0.5d0*((f3c + f3a)*mesh%edges(ee)%dy - (g3c + g3a)*mesh%edges(ee)%dx)
    mesh%edges_r4(ee) = 0.5d0*((f4c + f4a)*mesh%edges(ee)%dy - (g4c + g4a)*mesh%edges(ee)%dx)
end do 
return 
end subroutine get_edge_fluxes

!get edge dissipations ===============
subroutine get_edge_dissipations(mesh,options)
implicit none 

!variables - inout
type(flux_mesh) :: mesh 
type(flux_options) :: options 

!variables - local 
integer(in32) :: ee 
integer(in32) :: c1,c2
real(dp) :: s2,s4,psi0,psi1,psi01,dk2,dk4

!evaluate each edge dissipation 
do ee=1,mesh%nedge
    if (mesh%edges(ee)%c2 .GT. 0) then

        !cells on this edge 
        c1 = mesh%edges(ee)%c1
        c2 = mesh%edges(ee)%c2

        !cell edge quantity scaling values ******************* update these
        s2 = 1.0d0
        s4 = 0.25d0 

        !4th order psi coefficient 
        psi0 = sqrt(mesh%cells_specrad(c1)/(4.0d0*mesh%edges_specrad(ee)))
        psi1 = sqrt(mesh%cells_specrad(c2)/(4.0d0*mesh%edges_specrad(ee)))
        psi01 = (4.0d0*(psi0*psi1))/(psi0 + psi1)

        !edge dissipation coefficients 
        dk2 = 0.5d0*(mesh%cells_psensor(c1) + mesh%cells_psensor(c2))*s2*options%k2*mesh%edges_specrad(ee)
        dk4 = max(0.0d0,(options%k4*mesh%edges_specrad(ee) - 2.0d0*dK2))*S4

        !evaluate edge dissipations 
        mesh%edges_d1(ee) = (dk2*mesh%edges_l1(ee) + dk4*(mesh%l1(mesh%edges(ee)%c1) - mesh%l1(mesh%edges(ee)%c2)))*psi01
        mesh%edges_d2(ee) = (dk2*mesh%edges_l2(ee) + dk4*(mesh%l2(mesh%edges(ee)%c1) - mesh%l2(mesh%edges(ee)%c2)))*psi01
        mesh%edges_d3(ee) = (dk2*mesh%edges_l3(ee) + dk4*(mesh%l3(mesh%edges(ee)%c1) - mesh%l3(mesh%edges(ee)%c2)))*psi01
        mesh%edges_d4(ee) = (dk2*mesh%edges_l4(ee) + dk4*(mesh%l4(mesh%edges(ee)%c1) - mesh%l4(mesh%edges(ee)%c2)))*psi01
    else
        mesh%edges_d1(ee) = 0.0d0 
        mesh%edges_d2(ee) = 0.0d0 
        mesh%edges_d3(ee) = 0.0d0 
        mesh%edges_d4(ee) = 0.0d0 
    end if 
end do 
return 
end subroutine get_edge_dissipations

!cell pressure sensor ===============
subroutine cell_pressure_sensor(mesh,c) 
implicit none

!variables - inout
integer(in32) :: c
type(flux_mesh) :: mesh 

!variables - local
integer(in32) :: ee
real(dp) :: pn,pd

!evaluate
pn = 0.0d0 
pd = 0.0d0 
do ee=1,mesh%cells(c)%nedge
    pn = pn + mesh%edges_pn(mesh%cells(c)%edge(ee))*mesh%cells(c)%edge_sign(ee)
    pd = pd + mesh%edges_pd(mesh%cells(c)%edge(ee))
end do 
mesh%cells_psensor(c) = abs(pn/pd)
return 
end subroutine cell_pressure_sensor

!edge pressure sensor ===============
subroutine edge_pressure_sensor(mesh,e) 
implicit none

!variables - inout
integer(in32) :: e
type(flux_mesh) :: mesh 

!evaluate 
if (mesh%edges(e)%c2 .GT. 0) then 
    mesh%edges_pn(e) = mesh%p(mesh%edges(e)%c1) - mesh%p(mesh%edges(e)%c2)
    mesh%edges_pd(e) = mesh%p(mesh%edges(e)%c1) + mesh%p(mesh%edges(e)%c2)
else
    mesh%edges_pn(e) = 0.0d0
    mesh%edges_pd(e) = 0.0d0
end if 
return 
end subroutine edge_pressure_sensor

!cell laplacian ===============
subroutine cell_laplacian(mesh,c) 
implicit none

!variables - inout
integer(in32) :: c
type(flux_mesh) :: mesh 

!variables - local
integer(in32) :: ee

!evaluate 
mesh%l1(c) = 0.0d0 
mesh%l2(c) = 0.0d0 
mesh%l3(c) = 0.0d0 
mesh%l4(c) = 0.0d0 
do ee=1,mesh%cells(c)%nedge
    mesh%l1(c) = mesh%l1(c) + mesh%edges_l1(mesh%cells(c)%edge(ee))*mesh%cells(c)%edge_sign(ee)
    mesh%l2(c) = mesh%l2(c) + mesh%edges_l2(mesh%cells(c)%edge(ee))*mesh%cells(c)%edge_sign(ee)
    mesh%l3(c) = mesh%l3(c) + mesh%edges_l3(mesh%cells(c)%edge(ee))*mesh%cells(c)%edge_sign(ee)
    mesh%l4(c) = mesh%l4(c) + mesh%edges_l4(mesh%cells(c)%edge(ee))*mesh%cells(c)%edge_sign(ee)
end do 
return 
end subroutine cell_laplacian

!edge laplacian ===============
subroutine edge_laplacian(mesh,e) 
implicit none

!variables - inout
integer(in32) :: e
type(flux_mesh) :: mesh 

!evaluate 
if (mesh%edges(e)%c2 .GT. 0) then 
    mesh%edges_l1(e) = mesh%w1(mesh%edges(e)%c1) - mesh%w1(mesh%edges(e)%c2)
    mesh%edges_l2(e) = mesh%w2(mesh%edges(e)%c1) - mesh%w2(mesh%edges(e)%c2)
    mesh%edges_l3(e) = mesh%w3(mesh%edges(e)%c1) - mesh%w3(mesh%edges(e)%c2)
    mesh%edges_l4(e) = mesh%w4(mesh%edges(e)%c1) - mesh%w4(mesh%edges(e)%c2)
else
    mesh%edges_l1(e) = 0.0d0 
    mesh%edges_l2(e) = 0.0d0 
    mesh%edges_l3(e) = 0.0d0 
    mesh%edges_l4(e) = 0.0d0 
end if 
return 
end subroutine edge_laplacian

!edge spectral radius ===============
subroutine edge_spectral_radius(mesh,options,e)
implicit none 

!variables - inout
integer(in32) :: e
type(flux_mesh) :: mesh 
type(flux_options) :: options 

!variables - local 
real(dp) :: ue,ve,sose

!evaluate 
if (mesh%edges(e)%c2 .GT. 0) then !cell 
    ue = 0.5d0*(mesh%u(mesh%edges(e)%c1) + mesh%u(mesh%edges(e)%c2))
    ve = 0.5d0*(mesh%v(mesh%edges(e)%c1) + mesh%v(mesh%edges(e)%c2))
    sose = 0.5d0*(speed_of_sound(mesh%p(mesh%edges(e)%c1),mesh%rho(mesh%edges(e)%c1),options%gamma) + speed_of_sound(mesh%p(mesh%edges(e)%c2),mesh%rho(mesh%edges(e)%c2),options%gamma))
else !boundary condition 
    ue = mesh%u(mesh%edges(e)%c1)
    ve = mesh%v(mesh%edges(e)%c1)
    sose = speed_of_sound(mesh%p(mesh%edges(e)%c1),mesh%rho(mesh%edges(e)%c1),options%gamma)
end if 

!calculate edge flow spectral radius
mesh%edges_specrad(e) = (abs(mesh%edges(e)%nx*ue + mesh%edges(e)%ny*ve + sose) + sose)*mesh%edges(e)%length
return 
end subroutine edge_spectral_radius

end module flux_solve
