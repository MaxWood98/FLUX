!flux 2d flow solve module 
!max wood
!version : 0.0.4
!updated : 02-03-26

!module 
module flux_solve
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
    write(*,'(A,A,A)') '    {pressure: ',real2F0_Xstring(options%pinf,6_in64),'}' 
    write(*,'(A,A,A)') '    {density: ',real2F0_Xstring(options%rhoinf,6_in64),'}'
    write(*,'(A,A,A)') '    {speed of sound: ',real2F0_Xstring(options%cinf,6_in64),'}'
    write(*,'(A,A,A)') '    {velocity: ',real2F0_Xstring(options%velinf,6_in64),'}'
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
real(dp) :: rhores,rhores0

!initialise flags
resconv = .false.
nanflag = .false.

!set the number of threads 
call omp_set_num_threads(options%num_threads)

!initialise the dissipation arrays to zero (ensures this is ignored for schemes with no dissipation)
mesh%edges_d1(:) = 0.0d0 
mesh%edges_d2(:) = 0.0d0 
mesh%edges_d3(:) = 0.0d0 
mesh%edges_d4(:) = 0.0d0 

!initial display
if (options%cdisplay) then     
    write(*,'(A)') '    ittn    |   rhores   |     cl     |     cd     '
    write(*,'(A)') '+-----------+------------+------------+------------+'
end if 

!initialise the parallel region 
!$OMP parallel

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
    !$OMP do schedule(guided,50)
    do ee=1,mesh%nedge
        call edge_spectral_radius(mesh,options,ee)
    end do 
    !$OMP end do 

    !evaluate cell spectral radii and timesteps 
    !$OMP do schedule(guided,50)
    do cc=1,mesh%ncell
        mesh%cells_specrad(cc) = 0.0d0 
        do ee=1,mesh%cells(cc)%nedge
            mesh%cells_specrad(cc) = mesh%cells_specrad(cc) + mesh%edges_specrad(mesh%cells(cc)%edge(ee))
        end do
        mesh%cells_dt(cc) = options%cfl*(mesh%cells_volume(cc)/mesh%cells_specrad(cc))
    end do 
    !$OMP end do 

    !rk iterate this timestep 
    call rk_iterate(mesh,options,nanflag)

    !initiate single threaded region 
    !$OMP single
    
    !evaluate forces
    call get_forces(mesh,options)

    !check convergence
    if (iteration == 1) then 
        rhores0 = sqrt(sum((mesh%residual)**2))
    end if 
    rhores = log10(sqrt(sum((mesh%residual)**2))/rhores0)
    mesh%rhores = rhores
    if (rhores .LE. options%residual_convtol) then 
        resconv = .true.
    end if 

    !display                 
    if (options%cdisplay) then     
        if (mod(iteration,100) == 0) then
            write(*,'(A)') '    ittn    |   rhores   |     cl     |     cd     '
            write(*,'(A)') '+-----------+------------+------------+------------+'
        end if 
            write(*,'(A,I8,A,F11.8,A,F11.8,A,F11.8)') '    ',iteration,'  ',rhores,'  ',mesh%cl,'  ',mesh%cd
    end if
    !$OMP end single
end do 

!end the parallel region
!$OMP end parallel 

!set derived flow properties
do cc=1,mesh%ncell
    mesh%mach(cc) = sqrt(mesh%u(cc)**2 + mesh%v(cc)**2)/speed_of_sound(mesh%p(cc),mesh%rho(cc),options%gamma)
    mesh%cp(cc) = pressure_coefficient(mesh%p(cc),options)
end do 

!evaluate forces
call get_forces(mesh,options)
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

!$OMP single

!store the initial state of each cell 
mesh%w10 = mesh%w1
mesh%w20 = mesh%w2
mesh%w30 = mesh%w3
mesh%w40 = mesh%w4

!set the edge dissipation initial values 

!$OMP end single


!iterate
do rr=1,options%rk_niter

    !evaluate edge fluxes
    call get_edge_fluxes(mesh,options,options%rk_dissipation(rr))

    !evaluate edge dissipations 
    if (options%rk_dissipation(rr)) then 
        !$OMP do schedule(guided,50)
        do cc=1,mesh%ncell
            call cell_pressure_sensor(mesh,cc) 
            call cell_laplacian(mesh,cc) 
        end do 
        !$OMP end do 
        call get_edge_dissipations(mesh,options)
    end if  

    !timestep each cell 
    !$OMP do schedule(guided,50)
    do cc=1,mesh%ncell
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
        cell_stepsize = rk4_alpha(rr)*(mesh%cells_dt(cc)/mesh%cells_volume(cc))
        mesh%w1(cc) = mesh%w10(cc) - cell_stepsize*r1
        mesh%w2(cc) = mesh%w20(cc) - cell_stepsize*r2
        mesh%w3(cc) = mesh%w30(cc) - cell_stepsize*r3
        mesh%w4(cc) = mesh%w40(cc) - cell_stepsize*r4
        mesh%residual(cc) = r1
        call con2prim(mesh%rho(cc),mesh%u(cc),mesh%v(cc),mesh%p(cc),mesh%e(cc),options%gamma,mesh%w1(cc),mesh%w2(cc),mesh%w3(cc),mesh%w4(cc))
        if (isnan(mesh%w1(cc))) then 
            print *, 'cell nan: ',cc
            nanflag = .true.
            !exit 
        end if 
    end do 
    !$OMP end do 
    if (nanflag) then 
        exit 
    end if 
end do 
return 
end subroutine rk_iterate

!get edge fluxes ===============
subroutine get_edge_fluxes(mesh,options,dissipation)
implicit none 

!variables - inout
logical :: dissipation
type(flux_mesh) :: mesh 
type(flux_options) :: options 

!variables - local 
integer(in32) :: ee,c1,c2
real(dp) :: velnorm,machnorm
real(dp) :: rho1,u1,v1,p1,e1,rho2,u2,v2,p2,e2,fx1,fx2,fx3,fx4

!evaluate each edge flux 
!$OMP do schedule(guided,50) 
do ee=1,mesh%nedge
    c1 = mesh%edges(ee)%c1
    c2 = mesh%edges(ee)%c2
    if (dissipation) then !dissipation components
        call edge_pressure_sensor(mesh,ee) 
        call edge_laplacian(mesh,ee)
    end if 
    if (c2 .GT. 0) then !internal cell

        !extract properties
        rho1 = mesh%rho(c1)
        u1 = mesh%u(c1)
        v1 = mesh%v(c1)
        p1 = mesh%p(c1)
        e1 = mesh%e(c1)
        rho2 = mesh%rho(c2) 
        u2 = mesh%u(c2)
        v2 = mesh%v(c2)
        p2 = mesh%p(c2)
        e2 = mesh%e(c2)

        !evaluate the flux
        if (options%flux_method == 'vanleer') then 
            call vanleer_flux(rho1,u1,v1,e1,p1,rho2,u2,v2,e2,p2,options%gamma,mesh%edges(ee)%nx,mesh%edges(ee)%ny,mesh%edges(ee)%length,fx1,fx2,fx3,fx4)
        elseif (options%flux_method == 'ausm') then 
            call ausm_flux(rho1,u1,v1,e1,p1,rho2,u2,v2,e2,p2,options%gamma,mesh%edges(ee)%nx,mesh%edges(ee)%ny,mesh%edges(ee)%length,fx1,fx2,fx3,fx4)
        elseif (options%flux_method == 'roe') then 
            call roe_flux(rho1,u1,v1,e1,p1,rho2,u2,v2,e2,p2,options%gamma,mesh%edges(ee)%nx,mesh%edges(ee)%ny,mesh%edges(ee)%length,fx1,fx2,fx3,fx4)
        elseif (options%flux_method == 'jameson') then 
            call jameson_flux(rho1,u1,v1,e1,p1,rho2,u2,v2,e2,p2,mesh%edges(ee)%nx,mesh%edges(ee)%ny,mesh%edges(ee)%length,fx1,fx2,fx3,fx4)
        end if 
    elseif (c2 == -1) then !wall   
        rho1 = 0.0
        u1 = 0.0
        v1 = 0.0
        p1 = mesh%p(c1)
        e1 = 0.0
        rho2 = 0.0
        u2 = 0.0
        v2 = 0.0
        p2 = mesh%p(c1)
        e2 = 0.0
        call jameson_flux(rho1,u1,v1,e1,p1,rho2,u2,v2,e2,p2,mesh%edges(ee)%nx,mesh%edges(ee)%ny,mesh%edges(ee)%length,fx1,fx2,fx3,fx4)
    elseif (c2 == -2) then !farfield
        rho1 = mesh%rho(c1)
        u1 = mesh%u(c1)
        v1 = mesh%v(c1)
        p1 = mesh%p(c1)
        e1 = mesh%e(c1)
        velnorm = mesh%u(c1)*mesh%edges(ee)%nx + mesh%v(c1)*mesh%edges(ee)%ny
        machnorm = abs(velnorm)/speed_of_sound(mesh%p(c1),mesh%rho(c1),options%gamma)
        if (machnorm .GE. 1.0) then !supersonic 
            if (velnorm .LT. 0.0d0) then !inflow
                call farfield_supersonic_bc_prescribed(rho2,u2,v2,p2,e2,mesh,c1,options,1)
            else !outflow
                call farfield_supersonic_bc_prescribed(rho2,u2,v2,p2,e2,mesh,c1,options,-1)
            end if 
        else !subsonic 
            if (velnorm .LT. 0.0d0) then !inflow
                call farfield_subsonic_bc_characteristic(rho2,u2,v2,p2,e2,mesh,c1,ee,options,1)
            else !outflow
                call farfield_subsonic_bc_characteristic(rho2,u2,v2,p2,e2,mesh,c1,ee,options,-1)
            end if 
        end if 
        call jameson_flux(rho1,u1,v1,e1,p1,rho2,u2,v2,e2,p2,mesh%edges(ee)%nx,mesh%edges(ee)%ny,mesh%edges(ee)%length,fx1,fx2,fx3,fx4)
    elseif (c2 == -3) then !stagnation inflow
        rho1 = mesh%rho(c1)
        u1 = mesh%u(c1)
        v1 = mesh%v(c1)
        p1 = mesh%p(c1)
        e1 = mesh%e(c1)
        call subsonic_stagnation_inflow_bc(rho2,u2,v2,p2,e2,mesh,c1,ee,options)
        call jameson_flux(rho1,u1,v1,e1,p1,rho2,u2,v2,e2,p2,mesh%edges(ee)%nx,mesh%edges(ee)%ny,mesh%edges(ee)%length,fx1,fx2,fx3,fx4) 
    elseif (c2 == -4) then !pressure outflow 
        rho1 = mesh%rho(c1)
        u1 = mesh%u(c1)
        v1 = mesh%v(c1)
        p1 = mesh%p(c1)
        e1 = mesh%e(c1)
        velnorm = mesh%u(c1)*mesh%edges(ee)%nx + mesh%v(c1)*mesh%edges(ee)%ny
        machnorm = abs(velnorm)/speed_of_sound(mesh%p(c1),mesh%rho(c1),options%gamma)
        if (machnorm .GE. 1.0) then !supersonic 
            call farfield_supersonic_bc_prescribed(rho2,u2,v2,p2,e2,mesh,c1,options,-1)
        else !subsonic 
            call subsonic_pressure_outflow_bc_characteristic(rho2,u2,v2,p2,e2,mesh,c1,ee,options,options%outflow_pratio)
        end if 
        call jameson_flux(rho1,u1,v1,e1,p1,rho2,u2,v2,e2,p2,mesh%edges(ee)%nx,mesh%edges(ee)%ny,mesh%edges(ee)%length,fx1,fx2,fx3,fx4)
    else 

        !todo: add other boundary conditions here

    end if 

    !evaluate the residuals
    mesh%edges_r1(ee) = fx1
    mesh%edges_r2(ee) = fx2
    mesh%edges_r3(ee) = fx3
    mesh%edges_r4(ee) = fx4
end do 
!$OMP end do 
return 
end subroutine get_edge_fluxes

!get edge fluxes (complex) ===============
subroutine get_edge_fluxes_cpx(mesh,options,dissipation)
implicit none 

!variables - inout
logical :: dissipation
type(flux_mesh_cpx) :: mesh 
type(flux_options) :: options 

!variables - local 
integer(in32) :: ee,c1,c2
complex(dp) :: velnorm,machnorm,gamma_cpx
complex(dp) :: rho1,u1,v1,p1,e1,rho2,u2,v2,p2,e2,fx1,fx2,fx3,fx4

!get complex gamma
gamma_cpx = complex(options%gamma,0.0d0)

!evaluate each edge flux 
!$OMP do schedule(guided,50) 
do ee=1,mesh%nedge
    c1 = mesh%edges(ee)%c1
    c2 = mesh%edges(ee)%c2
    if (dissipation) then !dissipation components
    !     call edge_pressure_sensor(mesh,ee) 
    !     call edge_laplacian(mesh,ee)
    end if 
    if (c2 .GT. 0) then !internal cell

        !extract properties
        rho1 = mesh%rho(c1)
        u1 = mesh%u(c1)
        v1 = mesh%v(c1)
        p1 = mesh%p(c1)
        e1 = mesh%e(c1)
        rho2 = mesh%rho(c2) 
        u2 = mesh%u(c2)
        v2 = mesh%v(c2)
        p2 = mesh%p(c2)
        e2 = mesh%e(c2)

        !evaluate the flux
        if (options%flux_method == 'vanleer') then 
            call vanleer_flux_cpx(rho1,u1,v1,e1,p1,rho2,u2,v2,e2,p2,gamma_cpx,mesh%edges(ee)%nx,mesh%edges(ee)%ny,mesh%edges(ee)%length,fx1,fx2,fx3,fx4)
        elseif (options%flux_method == 'ausm') then 
            ! call ausm_flux_cpx(rho1,u1,v1,e1,p1,rho2,u2,v2,e2,p2,gamma_cpx,mesh%edges(ee)%nx,mesh%edges(ee)%ny,mesh%edges(ee)%length,fx1,fx2,fx3,fx4)
        elseif (options%flux_method == 'roe') then 
            ! call roe_flux_cpx(rho1,u1,v1,e1,p1,rho2,u2,v2,e2,p2,gamma_cpx,mesh%edges(ee)%nx,mesh%edges(ee)%ny,mesh%edges(ee)%length,fx1,fx2,fx3,fx4)
        elseif (options%flux_method == 'jameson') then 
            call jameson_flux_cpx(rho1,u1,v1,e1,p1,rho2,u2,v2,e2,p2,mesh%edges(ee)%nx,mesh%edges(ee)%ny,mesh%edges(ee)%length,fx1,fx2,fx3,fx4)
        end if 
    elseif (c2 == -1) then !wall   
        rho1 = 0.0
        u1 = 0.0
        v1 = 0.0
        p1 = mesh%p(c1)
        e1 = 0.0
        rho2 = 0.0
        u2 = 0.0
        v2 = 0.0
        p2 = mesh%p(c1)
        e2 = 0.0
        call jameson_flux_cpx(rho1,u1,v1,e1,p1,rho2,u2,v2,e2,p2,mesh%edges(ee)%nx,mesh%edges(ee)%ny,mesh%edges(ee)%length,fx1,fx2,fx3,fx4)
    elseif (c2 == -2) then !farfield
        rho1 = mesh%rho(c1)
        u1 = mesh%u(c1)
        v1 = mesh%v(c1)
        p1 = mesh%p(c1)
        e1 = mesh%e(c1)
        velnorm = mesh%u(c1)*mesh%edges(ee)%nx + mesh%v(c1)*mesh%edges(ee)%ny
        machnorm = sqrt(velnorm**2)/speed_of_sound_cpx(mesh%p(c1),mesh%rho(c1),gamma_cpx)
        if (real(machnorm,dp) .GE. 1.0) then !supersonic 
            if (real(velnorm,dp) .LT. 0.0d0) then !inflow
                call farfield_supersonic_bc_prescribed_cpx(rho2,u2,v2,p2,e2,mesh,c1,options,1)
            else !outflow
                call farfield_supersonic_bc_prescribed_cpx(rho2,u2,v2,p2,e2,mesh,c1,options,-1)
            end if 
        else !subsonic 
            if (real(velnorm,dp) .LT. 0.0d0) then !inflow
                call farfield_subsonic_bc_characteristic_cpx(rho2,u2,v2,p2,e2,mesh,c1,ee,options,1)
            else !outflow
                call farfield_subsonic_bc_characteristic_cpx(rho2,u2,v2,p2,e2,mesh,c1,ee,options,-1)
            end if 
        end if 
        call jameson_flux_cpx(rho1,u1,v1,e1,p1,rho2,u2,v2,e2,p2,mesh%edges(ee)%nx,mesh%edges(ee)%ny,mesh%edges(ee)%length,fx1,fx2,fx3,fx4)
    elseif (c2 == -3) then !stagnation inflow
        rho1 = mesh%rho(c1)
        u1 = mesh%u(c1)
        v1 = mesh%v(c1)
        p1 = mesh%p(c1)
        e1 = mesh%e(c1)
        call subsonic_stagnation_inflow_bc_cpx(rho2,u2,v2,p2,e2,mesh,c1,ee,options)
        call jameson_flux_cpx(rho1,u1,v1,e1,p1,rho2,u2,v2,e2,p2,mesh%edges(ee)%nx,mesh%edges(ee)%ny,mesh%edges(ee)%length,fx1,fx2,fx3,fx4) 
    elseif (c2 == -4) then !pressure outflow 
        rho1 = mesh%rho(c1)
        u1 = mesh%u(c1)
        v1 = mesh%v(c1)
        p1 = mesh%p(c1)
        e1 = mesh%e(c1)
        velnorm = mesh%u(c1)*mesh%edges(ee)%nx + mesh%v(c1)*mesh%edges(ee)%ny
        machnorm = sqrt(velnorm**2)/speed_of_sound_cpx(mesh%p(c1),mesh%rho(c1),gamma_cpx)
        if (real(machnorm,dp) .GE. 1.0) then !supersonic 
            call farfield_supersonic_bc_prescribed_cpx(rho2,u2,v2,p2,e2,mesh,c1,options,-1)
        else !subsonic 
            call subsonic_pressure_outflow_bc_characteristic_cpx(rho2,u2,v2,p2,e2,mesh,c1,ee,options,complex(options%outflow_pratio,0.0d0))
        end if 
        call jameson_flux_cpx(rho1,u1,v1,e1,p1,rho2,u2,v2,e2,p2,mesh%edges(ee)%nx,mesh%edges(ee)%ny,mesh%edges(ee)%length,fx1,fx2,fx3,fx4)
    else 

        !todo: add other boundary conditions here

    end if 

    !evaluate the residuals
    mesh%edges_r1(ee) = fx1
    mesh%edges_r2(ee) = fx2
    mesh%edges_r3(ee) = fx3
    mesh%edges_r4(ee) = fx4
end do 
!$OMP end do 
return 
end subroutine get_edge_fluxes_cpx

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
!$OMP do schedule(guided,50)
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
        ! dk2 = abs(mesh%edges_pn(ee)/mesh%edges_pd(ee))*s2*options%k2*mesh%edges_specrad(ee)
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
!$OMP end do 
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

!get forces ===============
subroutine get_forces(mesh,options)
implicit none 

!variables - inout
type(flux_mesh) :: mesh 
type(flux_options) :: options 

!variables - local 
integer(in32) :: ee  

!initialise 
mesh%cx = 0.0d0 
mesh%cy = 0.0d0 

!accumulate 
do ee=1,mesh%nedge
    if (mesh%edges(ee)%c2 == -1) then 
        mesh%cx = mesh%cx + mesh%edges(ee)%dy*pressure_coefficient(mesh%p(mesh%edges(ee)%c1),options)
        mesh%cy = mesh%cy - mesh%edges(ee)%dx*pressure_coefficient(mesh%p(mesh%edges(ee)%c1),options)
    end if 
end do 

!evaluate
mesh%cd = cos(options%aoarad)*mesh%cx + sin(options%aoarad)*mesh%cy
mesh%cl = cos(options%aoarad)*mesh%cy - sin(options%aoarad)*mesh%cx
return 
end subroutine get_forces

end module flux_solve
