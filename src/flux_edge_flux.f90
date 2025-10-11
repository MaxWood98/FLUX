!flux 2d edge flux evaluation module 
!max wood
!version : 0.0.2
!updated : 11-10-25

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




end module edge_flux