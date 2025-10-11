!flux 2d data and methods module 
!max wood
!version : 0.0.1
!updated : 28-03-25

!module 
module flux_data_methods

!integer data types 
use ISO_FORTRAN_ENV, only: in16=>int16
use ISO_FORTRAN_ENV, only: in32=>int32
use ISO_FORTRAN_ENV, only: in64=>int64

!real data types 
use ISO_FORTRAN_ENV, only: sp=>real32 
use ISO_FORTRAN_ENV, only: dp=>real64

!define constant pi 
real(dp) :: pi = 4.0d0*atan(1.0d0)

!define rk 
real(dp) :: rk4_alpha(4) = (/1.0d0/4.0d0,1.0d0/3.0d0,1.0d0/2.0d0,1.0d0/) !rk4 

!options type
type flux_options 
    logical :: cdisplay
    logical, dimension(:), allocatable :: rk_dissipation
    integer(in32) :: niter_max,rk_niter,num_threads
    real(dp) :: cfl,aoadeg,aoarad,gamma,R,k2,k4,residual_convtol
    real(dp) :: machinf,tinf,rhoinf,cinf,pinf,velinf,uinf,vinf,p0inf,rho0inf,t0inf
    character(len=:), allocatable :: meshpath,meshname
end type flux_options 

!vertex type
type flux_vertex
    integer(in32) :: index
    real(dp) :: coordinate(2)
end type flux_vertex

!edge type
type flux_edge 
    integer(in32) :: index
    integer(in32) :: v1,v2,c1,c2
    real(dp) :: dx,dy,nx,ny,length
end type flux_edge

!cell type
type flux_cell 
    integer(in32) :: index
    integer(in32) :: nedge
    integer(in32), dimension(:), allocatable :: edgev1,edgev2,edgec,edge
    real(dp), dimension(:), allocatable :: edge_sign
end type flux_cell

!mesh type 
type flux_mesh
    integer(in32) :: nvertex,nedge,ncell 
    real(dp) :: cl,cd,cm,mflux_in,mflux_out
    real(dp), dimension(:), allocatable :: cells_specrad,cells_volume,cells_dt,cells_psensor
    real(dp), dimension(:), allocatable :: edges_specrad
    real(dp), dimension(:), allocatable :: rho,u,v,p,mach,e,cp
    real(dp), dimension(:), allocatable :: w1,w2,w3,w4,w10,w20,w30,w40,residual
    real(dp), dimension(:), allocatable :: l1,l2,l3,l4
    real(dp), dimension(:), allocatable :: edges_r1,edges_r2,edges_r3,edges_r4
    real(dp), dimension(:), allocatable :: edges_l1,edges_l2,edges_l3,edges_l4
    real(dp), dimension(:), allocatable :: edges_d1,edges_d2,edges_d3,edges_d4
    real(dp), dimension(:), allocatable :: edges_pn,edges_pd
    type(flux_vertex), dimension(:), allocatable :: vertices
    type(flux_edge), dimension(:), allocatable :: edges
    type(flux_cell), dimension(:), allocatable :: cells 
    contains 
        procedure :: get_edges_geometry
        procedure :: get_cells_volume
end type flux_mesh

!routines 
contains 

!======================================================
!class methods ========================================
!======================================================

!evaluate edges geometry ===============
subroutine get_edges_geometry(self) 
implicit none 

!variables - inout
class(flux_mesh) :: self 

!variables - local
integer(in32) :: ii 
real(dp) :: dx,dy 

!evaluate 
do ii=1,self%nedge
    dx = self%vertices(self%edges(ii)%v2)%coordinate(1) - self%vertices(self%edges(ii)%v1)%coordinate(1)
    dy = self%vertices(self%edges(ii)%v2)%coordinate(2) - self%vertices(self%edges(ii)%v1)%coordinate(2)
    self%edges(ii)%length = sqrt(dx*dx + dy*dy)
    self%edges(ii)%dx = dx 
    self%edges(ii)%dy = dy 
    self%edges(ii)%nx = dy/self%edges(ii)%length
    self%edges(ii)%ny = -dx/self%edges(ii)%length
end do 
return 
end subroutine get_edges_geometry

!evaluate cell volumes ===============
subroutine get_cells_volume(self)
implicit none 

!variables - inout
class(flux_mesh) :: self

!variables - local
integer(in32) :: ii,jj

!evaluate cell volumes 
do ii=1,self%ncell
    self%cells_volume(ii) = 0.0d0
    do jj=1,self%cells(ii)%nedge
        self%cells_volume(ii) = self%cells_volume(ii) + edge_volume(self%vertices(self%cells(ii)%edgev1(jj))%coordinate,self%vertices(self%cells(ii)%edgev2(jj))%coordinate) 
    end do 
end do 


! self%cells_volume(:) = 0.0d0 
! do ii=1,self%ncell
!     self%cells_volume(ii) = 0.0d0
!     do jj=1,self%cells(ii)%nedge
        
!         self%cells_volume(ii) = self%cells_volume(ii) + self%cells(ii)%edge_sign(jj)* &
!         edge_volume( self%vertices(self%edges(self%cells(ii)%edge(jj))%v1 )%coordinate,self%vertices( self%edges(self%cells(ii)%edge(jj))%v2 )%coordinate) 


!     end do 
! end do 
return 
end subroutine get_cells_volume







! !evaluate cell edge geometries ===============
! subroutine get_cell_edge_geometries(self)
! implicit none 

! !variables - inout
! class(flux_mesh) :: self

! !variables - local
! integer(in64) :: ii 

! !evaluate cell volumes 
! do ii=1,self%ncell
!     call self%cells(ii)%get_edge_geometry(self)
! end do 
! return 
! end subroutine get_cell_edge_geometries

! !evaluate cell edge geometry ===============
! subroutine get_edge_geometry(self,mesh) 
! implicit none 

! !variables - inout
! type(flux_mesh) :: mesh 
! class(flux_cell) :: self

! !variables - local
! integer(in64) :: ii 
! real(dp) :: dx,dy 

! !evaluate 
! do ii=1,self%nedge
!     dx = mesh%vertices(self%edges(ii,2),1) - mesh%vertices(self%edges(ii,1),1)
!     dy = mesh%vertices(self%edges(ii,2),2) - mesh%vertices(self%edges(ii,1),2)
!     self%edges_len(ii) = sqrt(dx*dx + dy*dy)
!     self%edges_dx(ii) = dx 
!     self%edges_dy(ii) = dy 
!     self%edges_nx(ii) = dy/self%edges_len(ii)
!     self%edges_ny(ii) = -dx/self%edges_len(ii)
! end do 
! return 
! end subroutine get_edge_geometry

! !primative to conservative variables ===============
! subroutine prim2con(self,gamma)
! implicit none 

! !variables - inout 
! real(dp) :: gamma
! class(flux_cell) :: self 

! !evaluate 
! self%w1 = self%rho
! self%w2 = self%rho*self%u 
! self%w3 = self%rho*self%v  
! self%w4 = self%rho*energy(self%p,self%rho,sqrt(self%u*self%u + self%v*self%v),gamma)
! return 
! end subroutine prim2con

! !conservative to primative variables ===============
! subroutine con2prim(self,gamma,options)
! implicit none 

! !variables - inout 
! real(dp) :: gamma
! class(flux_cell) :: self 
! type(flux_options) :: options 

! !variables - local 
! real(dp) :: vel2 

! !evaluate 
! self%rho = self%w1
! self%u = self%w2/self%w1
! self%v = self%w3/self%w1 
! self%e = self%w4/self%w1
! vel2 = self%u*self%u + self%v*self%v
! self%p = (self%e - 0.5d0*vel2)*(gamma - 1.0)*self%w1
! self%mach = sqrt(vel2)/sqrt(gamma*(self%p/self%rho))
! self%cp = pressure_coefficient(self%p,options)
! return 
! end subroutine con2prim

!======================================================
!general methods ======================================
!======================================================

!edge volume ===============
function edge_volume(v1,v2) result(aseg) !+ve area for CCW oriented shapes
implicit none 

!variables - import
real(dp) :: v1(2),v2(2)

!result
real(dp) :: aseg

!segment area
aseg = 0.5d0*(v1(1)*v2(2) - v2(1)*v1(2))
return 
end function edge_volume
   
!speed of sound ===============
function speed_of_sound(p,rho,gamma) result(c)
implicit none 

!variables - inout 
real(dp) :: p,rho,gamma,c 

!evaluate 
c = sqrt(gamma*(p/rho))
return 
end function speed_of_sound

!pressure coefficient ===============
function pressure_coefficient(p,options) result(cp)
implicit none 

!variables - inout 
real(dp) :: p,cp
type(flux_options) :: options 

!evaluate
cp =  ((p/options%pinf) - 1.0d0)/(0.5d0*options%gamma*options%machinf*options%machinf)
return 
end function pressure_coefficient 

!energy ===============
function energy(p,rho,vel,gamma) result(e)
implicit none 

!variables - inout 
real(dp) :: p,rho,vel,gamma,e 

!evaluate 
e = (p/((gamma - 1.0d0)*rho)) + 0.5d0*vel*vel
return 
end function energy

!stagnation pressure ===============
function stagnation_pressure(p,mach,gamma) result(p0)
implicit none 

!variables - inout 
real(dp) :: p,mach,gamma,p0 

!evaluate
p0 = p*(1.0d0 + 0.5d0*(gamma - 1.0d0)*mach*mach)**(gamma/(gamma - 1.0d0))
return 
end function stagnation_pressure 

!stagnation density ===============
function stagnation_density(rho,mach,gamma) result(rho0)
implicit none 

!variables - inout 
real(dp) :: rho,mach,gamma,rho0 

!evaluate
rho0 = rho*(1.0d0 + 0.5d0*(gamma - 1.0d0)*mach*mach)**(gamma/(gamma - 1.0d0))
return 
end function stagnation_density

!stagnation temperature ===============
function stagnation_temperature(temp,mach,gamma) result(temp0)
implicit none 

!variables - inout 
real(dp) :: temp,mach,gamma,temp0 

!evaluate
temp0 = temp*(1.0d0 + 0.5d0*(gamma - 1.0d0)*mach*mach)
return 
end function stagnation_temperature

!change basis =========================
function change_basis(Ma2b,va) result(vb)
implicit none

!variables - inout
real(dp) :: va(2),vb(2)
real(dp) :: Ma2b(2,2)

!evaluate 
vb(1) = Ma2b(1,1)*va(1) + Ma2b(1,2)*va(2)
vb(2) = Ma2b(2,1)*va(1) + Ma2b(2,2)*va(2)
return 
end function change_basis

!get change of basis matrices =========================
subroutine get_basis_change_2d(Ma2b,Mb2a,bx,by)
implicit none 

!variables - inout
real(dp) :: bx(2),by(2)
real(dp) :: Ma2b(2,2),Mb2a(2,2)

!variables - local 
real(dp) :: det 

!get transform from b2a (assuming bx and by are represented as vectors in basis a)
Mb2a(1,1) = bx(1)
Mb2a(1,2) = by(1)
Mb2a(2,1) = bx(2)
Mb2a(2,2) = by(2)

!get transform from a2b
det = 1.0d0/(Mb2a(1,1)*Mb2a(2,2) - Mb2a(1,2)*Mb2a(2,1))
Ma2b(1,1) = Mb2a(2,2)*det
Ma2b(1,2) = -Mb2a(1,2)*det
Ma2b(2,1) = -Mb2a(2,1)*det
Ma2b(2,2) = Mb2a(1,1)*det
return 
end subroutine get_basis_change_2d

!primative to conservative variables ===============
subroutine prim2con(rho,u,v,p,gamma,w1,w2,w3,w4)
implicit none 

!variables - inout
real(dp) :: rho,u,v,p,gamma,w1,w2,w3,w4

!evaluate
w1 = rho
w2 = rho*u 
w3 = rho*v  
w4 = rho*energy(p,rho,sqrt(u*u + v*v),gamma)
return 
end subroutine prim2con

!conservative to primative variables ===============
subroutine con2prim(rho,u,v,p,e,gamma,w1,w2,w3,w4)
implicit none 

!variables - inout
real(dp) :: rho,u,v,p,e,gamma,w1,w2,w3,w4

!variables - local 
real(dp) :: vel2

!evaluate
rho = w1
u = w2/w1
v = w3/w1
e = w4/w1
vel2 = u*u + v*v 
p = (e - 0.5d0*vel2)*(gamma - 1.0)*w1
return 
end subroutine con2prim



end module flux_data_methods