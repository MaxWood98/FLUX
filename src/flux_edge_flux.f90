!flux 2d edge flux evaluation module 
!max wood
!version : 0.0.2
!updated : 12-10-25

!module 
module edge_flux
use flux_data_methods
contains 

!cell flux ===============
subroutine cell_flux(f1,f2,f3,f4,g1,g2,g3,g4,mesh,c)
implicit none 

!variables - inout 
integer(in32) :: c
real(dp) :: f1,f2,f3,f4,g1,g2,g3,g4
type(flux_mesh) :: mesh 

!evaluate f
f1 = mesh%rho(c)*mesh%u(c)
f2 = mesh%rho(c)*mesh%u(c)*mesh%u(c) + mesh%p(c)
f3 = mesh%rho(c)*mesh%u(c)*mesh%v(c)
f4 = mesh%rho(c)*mesh%u(c)*mesh%e(c) + mesh%u(c)*mesh%p(c)

!evaluate g 
g1 = mesh%rho(c)*mesh%v(c)
g2 = mesh%rho(c)*mesh%v(c)*mesh%u(c)
g3 = mesh%rho(c)*mesh%v(c)*mesh%v(c) + mesh%p(c)
g4 = mesh%rho(c)*mesh%v(c)*mesh%e(c) + mesh%v(c)*mesh%p(c)
return 
end subroutine cell_flux

!wall flux ===============
subroutine wall_flux(f1,f2,f3,f4,g1,g2,g3,g4,mesh,c)
implicit none 

!variables - inout 
integer(in32) :: c
real(dp) :: f1,f2,f3,f4,g1,g2,g3,g4
type(flux_mesh) :: mesh 

!evaluate f
f1 = 0.0d0 
f2 = mesh%p(c) 
f3 = 0.0d0 
f4 = 0.0d0 

!evaluate g
g1 = 0.0d0
g2 = 0.0d0
g3 = mesh%p(c) 
g4 = 0.0d0
return 
end subroutine wall_flux

!supersonic farfield ===============
subroutine supersonic_farfield_flux(f1,f2,f3,f4,g1,g2,g3,g4,mesh,c,options,inflow_outflow)
implicit none 

!variables - inout 
integer(in32) :: c,inflow_outflow
real(dp) :: f1,f2,f3,f4,g1,g2,g3,g4
type(flux_mesh) :: mesh 
type(flux_options) :: options 

!variables - local 
real(dp) :: ub,vb,rhob,pb,eb

!get boundary state
if (inflow_outflow == 1) then !inflow
    ub = options%uinf
    vb = options%vinf
    rhob = options%rhoinf
    pb = options%pinf
elseif (inflow_outflow == -1) then !outflow
    ub = mesh%u(c)
    vb = mesh%v(c)
    rhob = mesh%rho(c)
    pb = mesh%p(c)
end if 
eb = energy(pb,rhob,sqrt(ub*ub + vb*vb),options%gamma)

!evaluate f
f1 = rhob*ub
f2 = rhob*ub*ub + pb 
f3 = rhob*ub*vb
f4 = rhob*ub*eb + ub*pb

!evaluate g 
g1 = rhob*vb
g2 = rhob*vb*ub
g3 = rhob*vb*vb + pb 
g4 = rhob*vb*eb + vb*pb
return 
end subroutine supersonic_farfield_flux

!far field characteristic boundary condition (subsonic) ===============
subroutine far_field_subsonic_flux_characteristic(f1,f2,f3,f4,g1,g2,g3,g4,mesh,c,e,options,inflow_outflow)
implicit none 

!variables - inout 
integer(in32) :: e,c,inflow_outflow
real(dp) :: f1,f2,f3,f4,g1,g2,g3,g4
type(flux_mesh) :: mesh 
type(flux_options) :: options 

!variables - local 
real(dp) :: a,du,dv,drho,dpr,ub,vb,rhob,pb,eb,c1,c2,c3,c4
real(dp) :: vel_in(2),vel_inf(2),vel_in_n(2),vel_inf_n(2),vel_b(2),vel_b_n(2)
real(dp) :: basis_bx(2),basis_by(2)
real(dp) :: Mb2a(2,2),Ma2b(2,2) 

!get the local speed of sound 
a = speed_of_sound(mesh%p(c),mesh%rho(c),options%gamma)

!get the local basis coordinate system wrt to the edge normal 
if (inflow_outflow == 1) then !inflow
    basis_bx(1) = -mesh%edges(e)%nx !nx
    basis_bx(2) = -mesh%edges(e)%ny !ny
elseif (inflow_outflow == -1) then !outflow
    basis_bx(1) = mesh%edges(e)%nx !nx
    basis_bx(2) = mesh%edges(e)%ny !ny
end if 
basis_by(1) = mesh%edges(e)%dx/mesh%edges(e)%length
basis_by(2) = mesh%edges(e)%dy/mesh%edges(e)%length
call get_basis_change_2d(Ma2b,Mb2a,basis_bx,basis_by)

!translate the velocity to this basis 
vel_in(1) = mesh%u(c)
vel_in(2) = mesh%v(c)
vel_inf(1) = options%uinf
vel_inf(2) = options%vinf 
vel_in_n = change_basis(Ma2b,vel_in)
vel_inf_n = change_basis(Ma2b,vel_inf)

!get the characteristic delta variables 
du = vel_in_n(1) - vel_inf_n(1)
dv = vel_in_n(2) - vel_inf_n(2)
drho = mesh%rho(c) - options%rhoinf 
dpr = mesh%p(c) - options%pinf 
c1 = -a*a*drho + dpr 
c2 = mesh%rho(c)*a*dv 
c3 = mesh%rho(c)*a*du + dpr 
c4 = -mesh%rho(c)*a*du + dpr  

!apply the boundary condition 
if (inflow_outflow == 1) then !inflow 
    c1 = 0.0d0 
    c2 = 0.0d0 
    c3 = 0.0d0 
elseif (inflow_outflow == -1) then !outflow 
    c4 = 0.0d0 
end if 

!get the boundary deltas 
drho = (-1.0d0/(a*a))*c1 + (1.0d0/(2.0d0*a*a))*c3 + (1.0d0/(2.0d0*a*a))*c4 
du = (1.0d0/(2.0d0*mesh%rho(c)*a))*(c3 - c4)
dv = (1.0d0/(mesh%rho(c)*a))*c2 
dpr = 0.5d0*(c3 + c4)

!evaluate the boundary velocity in the base coordinate system 
vel_b_n(1) = vel_inf_n(1) + du 
vel_b_n(2) = vel_inf_n(2) + dv
vel_b = change_basis(Mb2a,vel_b_n)

!evaluate the boundary state
ub = vel_b(1)
vb = vel_b(2)
rhob = drho + options%rhoinf
pb = dpr + options%pinf
eb = energy(pb,rhob,sqrt(ub*ub + vb*vb),options%gamma)

!evaluate f
f1 = rhob*ub
f2 = rhob*ub*ub + pb 
f3 = rhob*ub*vb
f4 = rhob*ub*eb + ub*pb

!evaluate g 
g1 = rhob*vb
g2 = rhob*vb*ub
g3 = rhob*vb*vb + pb 
g4 = rhob*vb*eb + vb*pb
return 
end subroutine far_field_subsonic_flux_characteristic

!stagnation inflow ===============
subroutine subsonic_stagnation_inflow_flux(f1,f2,f3,f4,g1,g2,g3,g4,mesh,c,e,options)
implicit none 

!variables - inout 
integer(in32) :: e,c
real(dp) :: f1,f2,f3,f4,g1,g2,g3,g4
type(flux_mesh) :: mesh 
type(flux_options) :: options 

!variables - local 
real(dp) :: sosi,veli2,Hi,jm,alpha,aq,bq,cq
real(dp) :: sosb,sosb1,sosb2,velb,machb,pb,tb,rhob,ub,vb,eb

!get the internal enthalpy upstream running riemann invarient 
sosi = speed_of_sound(mesh%p(c),mesh%rho(c),options%gamma)
veli2 = mesh%u(c)**2 + mesh%v(c)**2
Hi = ((sosi*sosi)/(options%gamma - 1.0)) + 0.5d0*veli2 
jm = -sqrt(veli2) + ((2.0d0*sosi)/(options%gamma - 1.0))

!solve for the boundary speed of sound 
alpha = options%gamma - 1.0d0 
aq = (1.0d0/alpha) + (2.0d0/(alpha*alpha))
bq = -2.0d0*jm/alpha
cq = 0.5d0*jm*jm - Hi 
sosb1 = (-bq + sqrt(bq*bq - 4.0*aq*cq))/(2.0d0*aq)
sosb2 = (-bq - sqrt(bq*bq - 4.0*aq*cq))/(2.0d0*aq)
if (sosb1 .GE. sosb2) then 
    sosb = sosb1
else
    sosb = sosb2
end if 

!get the boundary flow state
velb = ((2.0d0*sosb)/(options%gamma - 1.0)) - jm 
machb = velb/sosb
pb = options%p0inf/((1.0d0 + 0.5d0*(options%gamma - 1.0d0)*machb*machb)**(options%gamma/(options%gamma - 1.0d0)))
tb = options%t0inf/(1.0d0 + 0.5d0*(options%gamma - 1.0d0)*machb*machb)
rhob = pb/(options%R*tb)
ub = -mesh%edges(e)%nx*velb
vb = -mesh%edges(e)%ny*velb
eb = energy(pb,rhob,sqrt(ub*ub + vb*vb),options%gamma)

!evaluate boundary flux
if (velb .GE. 0.0d0) then 
    if (machb .GE. 1.0) then 

        print *, 'supersonic stagnation inflow needs implementing '

        stop 

    else

        !evaluate f
        f1 = rhob*ub
        f2 = rhob*ub*ub + pb 
        f3 = rhob*ub*vb
        f4 = rhob*ub*eb + ub*pb

        !evaluate g 
        g1 = rhob*vb
        g2 = rhob*vb*ub
        g3 = rhob*vb*vb + pb 
        g4 = rhob*vb*eb + vb*pb
    end if 
else !revert to wall 
    call wall_flux(f1,f2,f3,f4,g1,g2,g3,g4,mesh,c)
end if 
return 
end subroutine subsonic_stagnation_inflow_flux

!pressure outflow boundary condition (subsonic) ===============
subroutine subsonic_pressure_outflow_flux_characteristic(f1,f2,f3,f4,g1,g2,g3,g4,mesh,c,e,options,pratio)
implicit none 

!variables - inout 
integer(in32) :: e,c
real(dp) :: f1,f2,f3,f4,g1,g2,g3,g4,pratio
type(flux_mesh) :: mesh 
type(flux_options) :: options 

!variables - local 
real(dp) :: a,du,dv,drho,dpr,ub,vb,rhob,pb,eb,c1,c2,c3,c4
real(dp) :: vel_in(2),vel_in_n(2),vel_b(2),vel_b_n(2)
real(dp) :: basis_bx(2),basis_by(2)
real(dp) :: Mb2a(2,2),Ma2b(2,2) 

!get the local speed of sound 
a = speed_of_sound(mesh%p(c),mesh%rho(c),options%gamma)

!get the local basis coordinate system wrt to the edge normal 
basis_bx(1) = mesh%edges(e)%nx !nx
basis_bx(2) = mesh%edges(e)%ny !ny
basis_by(1) = mesh%edges(e)%dx/mesh%edges(e)%length
basis_by(2) = mesh%edges(e)%dy/mesh%edges(e)%length
call get_basis_change_2d(Ma2b,Mb2a,basis_bx,basis_by)

!translate the velocity to this basis 
vel_in(1) = mesh%u(c)
vel_in(2) = mesh%v(c)
vel_in_n = change_basis(Ma2b,vel_in)

!get the characteristic delta variables 
du = 0.0d0 
dv = 0.0d0 
drho = 0.0d0 
dpr = mesh%p(c) - options%pinf*pratio
c1 = -a*a*drho + dpr 
c2 = mesh%rho(c)*a*dv 
c3 = mesh%rho(c)*a*du + dpr 
c4 = -mesh%rho(c)*a*du + dpr  

!apply the boundary condition 
c4 = 0.0d0 

!get the boundary deltas 
drho = (-1.0d0/(a*a))*c1 + (1.0d0/(2.0d0*a*a))*c3 + (1.0d0/(2.0d0*a*a))*c4 
du = (1.0d0/(2.0d0*mesh%rho(c)*a))*(c3 - c4)
dv = (1.0d0/(mesh%rho(c)*a))*c2 
dpr = 0.5d0*(c3 + c4)

!evaluate the boundary velocity in the local coordinate system 
vel_b_n(1) = vel_in_n(1) + du 
vel_b_n(2) = vel_in_n(2) + dv

!evaluate boundary flux 
if (vel_b_n(1) .LT. 0.0d0) then !if inflow then revert to a wall condition
    call wall_flux(f1,f2,f3,f4,g1,g2,g3,g4,mesh,c)
else !normal

    !evaluate the boundary velocity in the base coordinate system 
    vel_b = change_basis(Mb2a,vel_b_n)

    !evaluate the boundary state
    ub = vel_b(1)
    vb = vel_b(2)
    rhob = drho + mesh%rho(c) 
    pb = dpr + options%pinf*pratio
    eb = energy(pb,rhob,sqrt(ub*ub + vb*vb),options%gamma)

    !evaluate f
    f1 = rhob*ub
    f2 = rhob*ub*ub + pb 
    f3 = rhob*ub*vb
    f4 = rhob*ub*eb + ub*pb

    !evaluate g 
    g1 = rhob*vb
    g2 = rhob*vb*ub
    g3 = rhob*vb*vb + pb 
    g4 = rhob*vb*eb + vb*pb
end if 
return 
end subroutine subsonic_pressure_outflow_flux_characteristic

end module edge_flux
