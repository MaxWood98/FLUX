!flux 2d edge flux evaluation module 
!max wood
!version : 0.0.1
!updated : 28-03-25

!module 
module edge_flux
use flux_data_methods
contains 

!cell flux ===============
subroutine cell_flux(f1,f2,f3,f4,g1,g2,g3,g4,cell)
implicit none 

!variables - inout 
real(dp) :: f1,f2,f3,f4,g1,g2,g3,g4
type(flux_cell) :: cell 

!evaluate f
f1 = cell%rho*cell%u
f2 = cell%rho*cell%u*cell%u + cell%p 
f3 = cell%rho*cell%u*cell%v
f4 = cell%rho*cell%u*cell%e + cell%u*cell%p

!evaluate g 
g1 = cell%rho*cell%v
g2 = cell%rho*cell%v*cell%u
g3 = cell%rho*cell%v*cell%v + cell%p 
g4 = cell%rho*cell%v*cell%e + cell%v*cell%p

! !evaluate f
! f1 = cell%rho*cell%u
! f2 = cell%rho*cell%u*cell%u + cell%p 
! f3 = cell%rho*cell%u*cell%v
! f4 = cell%rho*cell%u*(cell%e + cell%p)

! !evaluate g 
! g1 = cell%rho*cell%v
! g2 = cell%rho*cell%v*cell%u
! g3 = cell%rho*cell%v*cell%v + cell%p 
! g4 = cell%rho*cell%v*(cell%e + cell%p)
return 
end subroutine cell_flux

!wall flux ===============
subroutine wall_flux(f1,f2,f3,f4,g1,g2,g3,g4,cell)
implicit none 

!variables - inout 
real(dp) :: f1,f2,f3,f4,g1,g2,g3,g4
type(flux_cell) :: cell 

!evaluate f
f1 = 0.0d0 
f2 = cell%p 
f3 = 0.0d0 
f4 = 0.0d0 

!evaluate g
g1 = 0.0d0
g2 = 0.0d0
g3 = cell%p
g4 = 0.0d0
return 
end subroutine wall_flux



!supersonic test ===============
subroutine supersonic_test(f1,f2,f3,f4,g1,g2,g3,g4,cell,options,inflow_outflow)
implicit none 

!variables - inout 
integer(in64) :: inflow_outflow
real(dp) :: f1,f2,f3,f4,g1,g2,g3,g4
type(flux_cell) :: cell 
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
    ub = cell%u
    vb = cell%v
    rhob = cell%rho
    pb = cell%p
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
end subroutine supersonic_test



!far field characteristic boundary condition (subsonic) ===============
subroutine far_field_subsonic_flux_characteristic(f1,f2,f3,f4,g1,g2,g3,g4,cell,eidx,options,inflow_outflow)
implicit none 

!variables - inout 
integer(in64) :: eidx,inflow_outflow
real(dp) :: f1,f2,f3,f4,g1,g2,g3,g4
type(flux_cell) :: cell 
type(flux_options) :: options 

!variables - local 
real(dp) :: c,du,dv,drho,dpr,ub,vb,rhob,pb,eb,c1,c2,c3,c4
real(dp) :: vel_in(2),vel_inf(2),vel_in_n(2),vel_inf_n(2),vel_b(2),vel_b_n(2)
real(dp) :: basis_bx(2),basis_by(2)
real(dp) :: Mb2a(2,2),Ma2b(2,2) 


!add tangential velocity component passthrough ? *******


!get the local speed of sound 
c = speed_of_sound(cell%p,cell%rho,options%gamma)

!get the local basis coordinate system wrt to the edge normal 
if (inflow_outflow == 1) then !inflow
    basis_bx(1) = -cell%edges_nx(eidx) !nx
    basis_bx(2) = -cell%edges_ny(eidx) !ny
elseif (inflow_outflow == -1) then !outflow
    basis_bx(1) = cell%edges_nx(eidx) !nx
    basis_bx(2) = cell%edges_ny(eidx) !ny
end if 
basis_by(1) = cell%edges_dx(eidx)/cell%edges_len(eidx) 
basis_by(2) = cell%edges_dy(eidx)/cell%edges_len(eidx) 
call get_basis_change_2d(Ma2b,Mb2a,basis_bx,basis_by)

!translate the velocity to this basis 
vel_in(1) = cell%u
vel_in(2) = cell%v
vel_inf(1) = options%uinf
vel_inf(2) = options%vinf 
vel_in_n = change_basis(Ma2b,vel_in)
vel_inf_n = change_basis(Ma2b,vel_inf)

!get the characteristic delta variables 
du = vel_in_n(1) - vel_inf_n(1)
dv = vel_in_n(2) - vel_inf_n(2)
drho = cell%rho - options%rhoinf 
dpr = cell%p - options%pinf 
c1 = -c*c*drho + dpr 
c2 = cell%rho*c*dv 
c3 = cell%rho*c*du + dpr 
c4 = -cell%rho*c*du + dpr  

!apply the boundary condition 
if (inflow_outflow == 1) then !inflow 
    c1 = 0.0d0 
    c2 = 0.0d0 
    c3 = 0.0d0 
elseif (inflow_outflow == -1) then !outflow 
    c4 = 0.0d0 
end if 

!get the boundary deltas 
drho = (-1.0d0/(c*c))*c1 + (1.0d0/(2.0d0*c*c))*c3 + (1.0d0/(2.0d0*c*c))*c4 
du = (1.0d0/(2.0d0*cell%rho*c))*(c3 - c4)
dv = (1.0d0/(cell%rho*c))*c2 
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










! !change basis =========================
! function change_basis(M,vec_b0) result(vec_b1)
! implicit none

! !Variables - Import
! real(dp) :: vec_b0(2),vec_b1(2)
! real(dp) :: M(2,2)

! !Evaluate
! vec_b1(1) = M(1,1)*vec_b0(1) + M(1,2)*vec_b0(2)
! vec_b1(2) = M(2,1)*vec_b0(1) + M(2,2)*vec_b0(2)
! return 
! end function change_basis

! !get basis change ===============
! subroutine get_basis_change(Mb2a,Ma2b,basis_bx,basis_by,basis_ax,basis_ay)
! implicit none 

! !variables - inout
! real(dp) :: basis_bx(2),basis_by(2),basis_ax(2),basis_ay(2)
! real(dp) :: Mb2a(2,2),Ma2b(2,2)

! !variables - local 
! real(dp) :: det

! !evaluate basis change from b -> a
! Mb2a(1,1) = dot_product(basis_bx,basis_ax)
! Mb2a(2,1) = dot_product(basis_bx,basis_ay)
! Mb2a(1,2) = dot_product(basis_by,basis_ax)
! Mb2a(2,2) = dot_product(basis_by,basis_ay)

! !evaluate basis change from a -> b
! det = Mb2a(1,1)*Mb2a(2,2) - Mb2a(1,2)*Mb2a(2,1)
! det = 1.0d0/det
! Ma2b(1,1) = Mb2a(2,2)*det
! Ma2b(2,1) = -Mb2a(2,1)*det
! Ma2b(1,2) = -Mb2a(1,2)*det
! Ma2b(2,2) = Mb2a(1,1)*det 
! return 
! end subroutine get_basis_change


end module edge_flux