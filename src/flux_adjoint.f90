!flux 2d adjoint module 
!max wood
!version : 0.0.4
!updated : 02-03-26

module flux_adjoint
use flux_io
use flux_solve
use flux_data_methods
use ieee_arithmetic
contains 

!dMdotdW (mass flux objective) ===============
subroutine mflux_objective(mesh,options,dJdW,dJdX)
implicit none 

!variables - inout
type(flux_mesh) :: mesh 
type(flux_options) :: options 
real(dp), dimension(:), allocatable :: dJdW,dJdX

!variables - local 
integer(in32) :: cc 
real(dp) :: dmfluxw1,dmfluxw2,dmfluxw3,dmfluxw4

!allocate the jacobians
allocate(dJdW(4*mesh%ncell))
allocate(dJdX(2*mesh%nvertex))
dJdW(:) = 0.0d0 
dJdX(:) = 0.0d0 !this is zero on surfaces as mass flux is not a directly surface driven property (technically non-zero on the measurment surface but we will ignore this)

!get the conservative variables
do cc=1,mesh%ncell
    call prim2con(mesh%rho(cc),mesh%u(cc),mesh%v(cc),mesh%p(cc),options%gamma,mesh%w1(cc),mesh%w2(cc),mesh%w3(cc),mesh%w4(cc))
    mesh%cp(cc) = pressure_coefficient(mesh%p(cc),options)
end do 

!evaluate the jacobian
do cc=1,mesh%ncell

    !evaluate the primitive variable gradients
    call get_drhodw_dvelndw(mesh,options,cc,-3,dmfluxw1,dmfluxw2,dmfluxw3,dmfluxw4)

    !populate dJdW
    dJdW(cc) = dmfluxw1
    dJdW(cc+mesh%ncell) = dmfluxw2
    dJdW(cc+2*mesh%ncell) = dmfluxw3
    dJdW(cc+3*mesh%ncell) = dmfluxw4
end do 
return 
end subroutine mflux_objective

!get drhodw and dvelndw in cell c for edges with the specified boundary condition ===============
subroutine get_drhodw_dvelndw(mesh,options,c,bc,dmfluxw1,dmfluxw2,dmfluxw3,dmfluxw4)
implicit none 

!variables - inout
integer(in32) :: c,bc
! real(dp) :: grhow1,grhow2,grhow3,grhow4,gvelw1,gvelw2,gvelw3,gvelw4

real(dp) :: dmfluxw1,dmfluxw2,dmfluxw3,dmfluxw4

type(flux_mesh) :: mesh 
type(flux_options) :: options 

!variables - local 
integer(in32) :: ee 
integer(in32) :: etgt
real(dp) :: w1,w2,w3,w4,h
! real(dp) :: grho1,grho2,grho3,grho4
! real(dp) :: gu1,gu2,gu3,gu4
! real(dp) :: gv1,gv2,gv3,gv4
complex(dp) :: gmf1,gmf2,gmf3,gmf4
complex(dp) :: rhoc,uc,vc,pc,ec,gammac,velc

!set the complex step stepsize
h = 1e-40_dp  

!get complex gamma
gammac = complex(options%gamma,0.0d0)

!extract the conservative variables 
w1 = mesh%w1(c)
w2 = mesh%w2(c)
w3 = mesh%w3(c)
w4 = mesh%w4(c)

!evaluate the gradient components
gmf1 = 0.0d0 
call con2prim_cpx(rhoc,uc,vc,pc,ec,gammac,complex(w1,h),complex(w2,0.0d0),complex(w3,0.0d0),complex(w4,0.0d0))
do ee=1,mesh%cells(c)%nedge
    etgt = mesh%cells(c)%edge(ee)
    if (mesh%edges(etgt)%c2 == bc) then 
        velc = uc*complex(mesh%edges(etgt)%nx,0.0d0) + vc*complex(mesh%edges(etgt)%ny,0.0d0)
        gmf1 = gmf1 + rhoc*velc*complex(mesh%edges(etgt)%length,0.0d0)
    end if 
end do 
dmfluxw1 = aimag(gmf1)/h

gmf2 = 0.0d0 
call con2prim_cpx(rhoc,uc,vc,pc,ec,gammac,complex(w1,0.0d0),complex(w2,h),complex(w3,0.0d0),complex(w4,0.0d0))
do ee=1,mesh%cells(c)%nedge
    etgt = mesh%cells(c)%edge(ee)
    if (mesh%edges(etgt)%c2 == bc) then 
        velc = uc*complex(mesh%edges(etgt)%nx,0.0d0) + vc*complex(mesh%edges(etgt)%ny,0.0d0)
        gmf2 = gmf2 + rhoc*velc*complex(mesh%edges(etgt)%length,0.0d0)
    end if 
end do 
dmfluxw2 = aimag(gmf2)/h

gmf3 = 0.0d0 
call con2prim_cpx(rhoc,uc,vc,pc,ec,gammac,complex(w1,0.0d0),complex(w2,0.0d0),complex(w3,h),complex(w4,0.0d0))
do ee=1,mesh%cells(c)%nedge
    etgt = mesh%cells(c)%edge(ee)
    if (mesh%edges(etgt)%c2 == bc) then 
        velc = uc*complex(mesh%edges(etgt)%nx,0.0d0) + vc*complex(mesh%edges(etgt)%ny,0.0d0)
        gmf3 = gmf3 + rhoc*velc*complex(mesh%edges(etgt)%length,0.0d0)
    end if 
end do 
dmfluxw3 = aimag(gmf3)/h

gmf4 = 0.0d0 
call con2prim_cpx(rhoc,uc,vc,pc,ec,gammac,complex(w1,0.0d0),complex(w2,0.0d0),complex(w3,0.0d0),complex(w4,h))
do ee=1,mesh%cells(c)%nedge
    etgt = mesh%cells(c)%edge(ee)
    if (mesh%edges(etgt)%c2 == bc) then 
        velc = uc*complex(mesh%edges(etgt)%nx,0.0d0) + vc*complex(mesh%edges(etgt)%ny,0.0d0)
        gmf4 = gmf4 + rhoc*velc*complex(mesh%edges(etgt)%length,0.0d0)
    end if 
end do 
dmfluxw4 = aimag(gmf4)/h



! grho1 = aimag(rhoc)/h
! gu1 = aimag(uc)/h
! gv1 = aimag(vc)/h
! call con2prim_cpx(rhoc,uc,vc,pc,ec,gammac,complex(w1,0.0d0),complex(w2,h),complex(w3,0.0d0),complex(w4,0.0d0))
! ! grho2 = aimag(rhoc)/h
! ! gu2 = aimag(uc)/h
! ! gv2 = aimag(vc)/h
! call con2prim_cpx(rhoc,uc,vc,pc,ec,gammac,complex(w1,0.0d0),complex(w2,0.0d0),complex(w3,h),complex(w4,0.0d0))
! ! grho3 = aimag(rhoc)/h
! ! gu3 = aimag(uc)/h
! ! gv3 = aimag(vc)/h
! call con2prim_cpx(rhoc,uc,vc,pc,ec,gammac,complex(w1,0.0d0),complex(w2,0.0d0),complex(w3,0.0d0),complex(w4,h))
! grho4 = aimag(rhoc)/h
! gu4 = aimag(uc)/h
! gv4 = aimag(vc)/h

! !evaluate the sum of the gradients over all target edges in the cell 
! gvelw1 = 0.0d0 
! gvelw2 = 0.0d0 
! gvelw3 = 0.0d0 
! gvelw4 = 0.0d0 
! grhow1 = 0.0d0 
! grhow2 = 0.0d0 
! grhow3 = 0.0d0 
! grhow4 = 0.0d0 
! do ee=1,mesh%cells(c)%nedge
!     etgt = mesh%cells(c)%edge(ee)
!     if (mesh%edges(etgt)%c2 == bc) then 
!         gvelw1 = gvelw1 + (gu1*mesh%edges(etgt)%nx + gv1*mesh%edges(etgt)%ny)*mesh%edges(etgt)%length*mesh%rho(c)*mesh%cells(c)%edge_sign(ee)
!         gvelw2 = gvelw2 + (gu2*mesh%edges(etgt)%nx + gv2*mesh%edges(etgt)%ny)*mesh%edges(etgt)%length*mesh%rho(c)*mesh%cells(c)%edge_sign(ee)
!         gvelw3 = gvelw3 + (gu3*mesh%edges(etgt)%nx + gv3*mesh%edges(etgt)%ny)*mesh%edges(etgt)%length*mesh%rho(c)*mesh%cells(c)%edge_sign(ee)
!         gvelw4 = gvelw4 + (gu4*mesh%edges(etgt)%nx + gv4*mesh%edges(etgt)%ny)*mesh%edges(etgt)%length*mesh%rho(c)*mesh%cells(c)%edge_sign(ee)


!         grhow1 = grhow1 + grho1*(mesh%u(c)*mesh%edges(etgt)%nx + mesh%v(c)*mesh%edges(etgt)%ny)*mesh%edges(etgt)%length*mesh%cells(c)%edge_sign(ee)
!         grhow2 = grhow2 + grho2*(mesh%u(c)*mesh%edges(etgt)%nx + mesh%v(c)*mesh%edges(etgt)%ny)*mesh%edges(etgt)%length*mesh%cells(c)%edge_sign(ee)
!         grhow3 = grhow3 + grho3*(mesh%u(c)*mesh%edges(etgt)%nx + mesh%v(c)*mesh%edges(etgt)%ny)*mesh%edges(etgt)%length*mesh%cells(c)%edge_sign(ee)
!         grhow4 = grhow4 + grho4*(mesh%u(c)*mesh%edges(etgt)%nx + mesh%v(c)*mesh%edges(etgt)%ny)*mesh%edges(etgt)%length*mesh%cells(c)%edge_sign(ee)
!     end if 
! end do 

!debug fd check --
!
!debug fd check --
return 
end subroutine get_drhodw_dvelndw

!dCddW (drag objective) ===============
subroutine cd_objective(mesh,options,dJdW,dJdX)
implicit none 

!variables - inout
type(flux_mesh) :: mesh 
type(flux_options) :: options 
real(dp), dimension(:), allocatable :: dJdW,dJdX

!variables - local 
integer(in32) :: cc,ee 
real(dp) :: gradc,gw1,gw2,gw3,gw4,wcos,wsin
real(dp) :: dcddp(mesh%ncell)

!allocate the jacobians
allocate(dJdW(4*mesh%ncell))
allocate(dJdX(2*mesh%nvertex))
dJdW(:) = 0.0d0 
dJdX(:) = 0.0d0 

!get the pressure coefficient 
do cc=1,mesh%ncell
    call prim2con(mesh%rho(cc),mesh%u(cc),mesh%v(cc),mesh%p(cc),options%gamma,mesh%w1(cc),mesh%w2(cc),mesh%w3(cc),mesh%w4(cc))
    mesh%cp(cc) = pressure_coefficient(mesh%p(cc),options)
end do 

!set incidence weightings
wcos = cos(options%aoarad)
wsin = sin(options%aoarad)

!contribution from each surface edge (dcx/dcp)*(dcp/dp) = dcddp
mesh%cd = 0.0d0 
dcddp(:) = 0.0d0 
do ee=1,mesh%nedge
    if (mesh%edges(ee)%c2 == -1) then 

        !accumulate cd
        mesh%cd = mesh%cd + (wcos*mesh%edges(ee)%dy - wsin*mesh%edges(ee)%dx)*mesh%cp(mesh%edges(ee)%c1)

        !accumulate dcddp
        gradc = (wcos*mesh%edges(ee)%dy - wsin*mesh%edges(ee)%dx)*(2.0d0/(options%gamma*options%pinf*options%machinf*options%machinf))
        dcddp(mesh%edges(ee)%c1) = dcddp(mesh%edges(ee)%c1) + gradc

        !accumulate dJdX x (dcydx -> dcydy = 0)
        dJdX(mesh%edges(ee)%v1) = dJdX(mesh%edges(ee)%v1) + wsin*mesh%cp(mesh%edges(ee)%c1)
        dJdX(mesh%edges(ee)%v2) = dJdX(mesh%edges(ee)%v2) - wsin*mesh%cp(mesh%edges(ee)%c1) 

        !accumulate dJdX y (dcxdy -> dcxdx = 0)
        dJdX(mesh%edges(ee)%v1+mesh%nvertex) = dJdX(mesh%edges(ee)%v1+mesh%nvertex) - wcos*mesh%cp(mesh%edges(ee)%c1)
        dJdX(mesh%edges(ee)%v2+mesh%nvertex) = dJdX(mesh%edges(ee)%v2+mesh%nvertex) + wcos*mesh%cp(mesh%edges(ee)%c1)
    end if 
end do 

!evaluate the jacobian
do cc=1,mesh%ncell

    !evaluate the gradient of pressure wrt to each conservative variable
    call get_dpdw(mesh,options,cc,gw1,gw2,gw3,gw4)

    !populate dJdW
    dJdW(cc) = dcddp(cc)*gw1
    dJdW(cc+mesh%ncell) = dcddp(cc)*gw2
    dJdW(cc+2*mesh%ncell) = dcddp(cc)*gw3
    dJdW(cc+3*mesh%ncell) = dcddp(cc)*gw4
end do 

!display 
if (options%cdisplay) then
    write(*,'(A,f16.14,A)') '    {cd = ',mesh%cd,'}' 
end if 
return 
end subroutine cd_objective

!get dpdw in cell c ===============
subroutine get_dpdw(mesh,options,c,gw1,gw2,gw3,gw4)
implicit none 

!variables - inout
integer(in32) :: c
real(dp) :: gw1,gw2,gw3,gw4
type(flux_mesh) :: mesh 
type(flux_options) :: options 

!variables - local 
real(dp) :: w1,w2,w3,w4,h
complex(dp) :: rhoc,uc,vc,pc,ec,gammac

!set the complex step stepsize
h = 1e-40_dp  

!get complex gamma
gammac = complex(options%gamma,0.0d0)

!extract the conservative variables 
w1 = mesh%w1(c)
w2 = mesh%w2(c)
w3 = mesh%w3(c)
w4 = mesh%w4(c)

!evaluate the gradients 
call con2prim_cpx(rhoc,uc,vc,pc,ec,gammac,complex(w1,h),complex(w2,0.0d0),complex(w3,0.0d0),complex(w4,0.0d0))
gw1 = aimag(pc)/h
call con2prim_cpx(rhoc,uc,vc,pc,ec,gammac,complex(w1,0.0d0),complex(w2,h),complex(w3,0.0d0),complex(w4,0.0d0))
gw2 = aimag(pc)/h
call con2prim_cpx(rhoc,uc,vc,pc,ec,gammac,complex(w1,0.0d0),complex(w2,0.0d0),complex(w3,h),complex(w4,0.0d0))
gw3 = aimag(pc)/h
call con2prim_cpx(rhoc,uc,vc,pc,ec,gammac,complex(w1,0.0d0),complex(w2,0.0d0),complex(w3,0.0d0),complex(w4,h))
gw4 = aimag(pc)/h

!debug fd check --
! call con2prim(rhop,up,vp,pp,ep,options%gamma,w1+1e-8,w2,w3,w4)
! print *, (pp - mesh%p(c))/1e-8
!debug fd check --
return 
end subroutine get_dpdw

!colour cells ===============
subroutine colour_cells(mesh,options)
implicit none 

!variables - inout
type(flux_mesh) :: mesh 
type(flux_options) :: options 

!variables - local 
logical :: shared_colour
integer(in32) :: ii,cc,ee,ee2
integer(in32) :: colour,cadj,cadj2 

!colour cells 
mesh%cells_colour(:) = -1
do cc=1,mesh%ncell

    !initialise colour
    colour = 1

    !check against neighbours and update
    do ii=1,10*mesh%ncell
        shared_colour = .false. 
        do ee=1,mesh%cells(cc)%nedge
            if (mesh%edges(mesh%cells(cc)%edge(ee))%c1 == cc) then 
                cadj = mesh%edges(mesh%cells(cc)%edge(ee))%c2
            else
                cadj = mesh%edges(mesh%cells(cc)%edge(ee))%c1
            end if 
            if (cadj .GT. 0) then 
                if (mesh%cells_colour(cadj) == colour) then 
                    shared_colour = .true.
                    exit  
                end if 
                do ee2=1,mesh%cells(cadj)%nedge
                    if (mesh%edges(mesh%cells(cadj)%edge(ee2))%c1 == cadj) then 
                        cadj2 = mesh%edges(mesh%cells(cadj)%edge(ee2))%c2
                    else
                        cadj2 = mesh%edges(mesh%cells(cadj)%edge(ee2))%c1
                    end if 
                    if (cadj2 .GT. 0) then 
                        if (cadj2 .NE. cc) then 
                            if (mesh%cells_colour(cadj2) == colour) then 
                                shared_colour = .true.
                                exit  
                            end if 
                        end if 
                    end if 
                end do 
                if (shared_colour) then 
                    exit
                end if 
            end if
        end do 
        if (shared_colour) then 
            colour = colour + 1 
        else
            exit
        end if 
    end do 

    !set colour 
    mesh%cells_colour(cc) = colour
end do 

!display
if (options%cdisplay) then
    write(*,'(A,I0,A,I0,A)') '    {cell colour min/max = ',minval(mesh%cells_colour),'/',maxval(mesh%cells_colour),'}' 
end if 
return 
end subroutine colour_cells

!colour vertices ===============
subroutine colour_vertices(mesh,options)
implicit none 

!variables - inout
type(flux_mesh) :: mesh 
type(flux_options) :: options 

!variables - local 
logical :: shared_colour
integer(in32) :: ii,vv,ee,ee2
integer(in32) :: colour,cadj,vidx  

!colour vertices 
mesh%vertices_colour(:) = -1
do vv=1,mesh%nvertex

    !initialise colour
    colour = 1

    !check against neighbours and update
    do ii=1,10*mesh%nvertex
        shared_colour = .false. 
        do ee=1,mesh%vertices(vv)%ncell
            cadj = mesh%vertices(vv)%cells(ee)
            do ee2=1,mesh%cells(cadj)%nedge
                vidx = mesh%cells(cadj)%edgev1(ee2)
                if (vidx .NE. vv) then 
                    if (mesh%vertices_colour(vidx) == colour) then 
                        shared_colour = .true.
                        exit  
                    end if 
                end if 
            end do 
            if (shared_colour) then 
                exit
            end if 
        end do 
        if (shared_colour) then 
            colour = colour + 1 
        else
            exit
        end if 
    end do 

    !set colour 
    mesh%vertices_colour(vv) = colour
end do 

!display
if (options%cdisplay) then
    write(*,'(A,I0,A,I0,A)') '    {vertex colour min/max = ',minval(mesh%vertices_colour),'/',maxval(mesh%vertices_colour),'}' 
end if 
return 
end subroutine colour_vertices

!full flow jacobian ===============
subroutine build_flow_jacobian_full(mesh,options,dRdW)
implicit none

!variables - inout
type(flux_mesh) :: mesh 
type(flux_options) :: options 

real(dp), dimension(:,:), allocatable :: dRdW

!variables - local 
integer(in32) :: cc,vv,ee,cc2
integer(in32) :: cidx

real(dp) :: r1,r2,r3,r4,h
real(dp) :: w10(mesh%ncell),w20(mesh%ncell),w30(mesh%ncell),w40(mesh%ncell)
real(dp) :: r10(mesh%ncell),r20(mesh%ncell),r30(mesh%ncell),r40(mesh%ncell)
real(dp) :: pr1(mesh%ncell),pr2(mesh%ncell),pr3(mesh%ncell),pr4(mesh%ncell)

h = 1e-8_dp 

mesh%edges_d1(:) = 0.0d0 
mesh%edges_d2(:) = 0.0d0 
mesh%edges_d3(:) = 0.0d0 
mesh%edges_d4(:) = 0.0d0 
do cc=1,mesh%ncell
    call prim2con(mesh%rho(cc),mesh%u(cc),mesh%v(cc),mesh%p(cc),options%gamma,mesh%w1(cc),mesh%w2(cc),mesh%w3(cc),mesh%w4(cc))
end do 
w10 = mesh%w1 
w20 = mesh%w2 
w30 = mesh%w3 
w40 = mesh%w4 


call get_edge_fluxes(mesh,options,.False.)
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
    r10(cc) = r1
    r20(cc) = r2
    r30(cc) = r3
    r40(cc) = r4
end do 


!allocate the flow jacobian 
allocate(dRdW(4*mesh%ncell,4*mesh%ncell))
dRdW = 0.0d0 

!construct the jacobian 
do cc=1,mesh%ncell !perturb each cell
    write(*,'(A,I0,A,I0)') 'cell: ',cc,'/',mesh%ncell
    do vv=1,4 !perturb each conservative variable in this cell

        !update to complex step ==========================
        !perturb the conservative variable w_vv in cell cc
        mesh%w1 = w10 
        mesh%w2 = w20 
        mesh%w3 = w30 
        mesh%w4 = w40 
        if (vv == 1) then 
            mesh%w1(cc) = mesh%w1(cc) + h
        elseif (vv == 2) then 
            mesh%w2(cc) = mesh%w2(cc) + h
        elseif (vv == 3) then 
            mesh%w3(cc) = mesh%w3(cc) + h
        elseif (vv == 4) then 
            mesh%w4(cc) = mesh%w4(cc) + h
        end if 
        
        do cc2=1,mesh%ncell
            call con2prim(mesh%rho(cc2),mesh%u(cc2),mesh%v(cc2),mesh%p(cc2),mesh%e(cc2),options%gamma,mesh%w1(cc2),mesh%w2(cc2),mesh%w3(cc2),mesh%w4(cc2))
        end do 

        !evaluate the flow residual
        call get_edge_fluxes(mesh,options,.False.)
        do cc2=1,mesh%ncell
            pr1(cc2) = 0.0d0 
            pr2(cc2) = 0.0d0 
            pr3(cc2) = 0.0d0 
            pr4(cc2) = 0.0d0 
            do ee=1,mesh%cells(cc2)%nedge
                pr1(cc2) = pr1(cc2) + (mesh%edges_r1(mesh%cells(cc2)%edge(ee)) + mesh%edges_d1(mesh%cells(cc2)%edge(ee)))*mesh%cells(cc2)%edge_sign(ee)
                pr2(cc2) = pr2(cc2) + (mesh%edges_r2(mesh%cells(cc2)%edge(ee)) + mesh%edges_d2(mesh%cells(cc2)%edge(ee)))*mesh%cells(cc2)%edge_sign(ee)
                pr3(cc2) = pr3(cc2) + (mesh%edges_r3(mesh%cells(cc2)%edge(ee)) + mesh%edges_d3(mesh%cells(cc2)%edge(ee)))*mesh%cells(cc2)%edge_sign(ee)
                pr4(cc2) = pr4(cc2) + (mesh%edges_r4(mesh%cells(cc2)%edge(ee)) + mesh%edges_d4(mesh%cells(cc2)%edge(ee)))*mesh%cells(cc2)%edge_sign(ee)
            end do 
        end do 
        !update to complex step ==========================


        
        !set the column index
        cidx = cc + (vv - 1)*mesh%ncell

        !unpack residuals 
        dRdW(1:mesh%ncell,cidx) = (pr1(:) - r10(:))/h
        dRdW(mesh%ncell+1:2*mesh%ncell,cidx) = (pr2(:) - r20(:))/h
        dRdW(2*mesh%ncell+1:3*mesh%ncell,cidx) = (pr3(:) - r30(:))/h
        dRdW(3*mesh%ncell+1:4*mesh%ncell,cidx) = (pr4(:) - r40(:))/h
    end do
end do 
return
end subroutine build_flow_jacobian_full

!sparse flow jacobian ===============
subroutine build_flow_jacobian_sparse(mesh_cpx,options,dRdW)
implicit none

!variables - inout
type(flux_mesh_cpx) :: mesh_cpx 
type(flux_options) :: options 
type(csc_matrix) :: dRdW

!variables - local 
integer(in32) :: clr,cc,vv,ee,rr,aa
integer(in32) :: ncolour,r0,col_offset,row,col,cadj,nblock
integer(in32) :: column_offset_index
real(dp) :: h
complex(dp) :: w10(mesh_cpx%ncell),w20(mesh_cpx%ncell),w30(mesh_cpx%ncell),w40(mesh_cpx%ncell)
complex(dp) :: pr1c(mesh_cpx%ncell),pr2c(mesh_cpx%ncell),pr3c(mesh_cpx%ncell),pr4c(mesh_cpx%ncell)

!define the complex stepsize
h = 1e-40_dp 

!initialise dissipation
mesh_cpx%edges_d1(:) = complex(0.0d0,0.0d0) 
mesh_cpx%edges_d2(:) = complex(0.0d0,0.0d0) 
mesh_cpx%edges_d3(:) = complex(0.0d0,0.0d0) 
mesh_cpx%edges_d4(:) = complex(0.0d0,0.0d0) 

!initialise conservative variables
do cc=1,mesh_cpx%ncell
    call prim2con_cpx(mesh_cpx%rho(cc),mesh_cpx%u(cc),mesh_cpx%v(cc),mesh_cpx%p(cc),complex(options%gamma,0.0d0),mesh_cpx%w1(cc),mesh_cpx%w2(cc),mesh_cpx%w3(cc),mesh_cpx%w4(cc))
end do 
w10 = mesh_cpx%w1 
w20 = mesh_cpx%w2 
w30 = mesh_cpx%w3 
w40 = mesh_cpx%w4 

!allocate the sparse flow jacobian 
dRdW%nnz = (sum(mesh_cpx%cells_nadj) + mesh_cpx%ncell)*16 
if (options%cdisplay) then
    write(*,'(A,I0,A)') '    {flow jacobian nnz = ',dRdW%nnz,'}'
    write(*,'(A,F8.6,A)') '    {flow jacobian sparsity = ',(real(dRdW%nnz,dp)/real(16_in64*mesh_cpx%ncell*mesh_cpx%ncell,dp))*100.0d0,'%}'
end if 
dRdW%nrow = 4*mesh_cpx%ncell
dRdW%ncol = 4*mesh_cpx%ncell
allocate(dRdW%row(dRdW%nnz))
allocate(dRdW%column(dRdW%nnz))
allocate(dRdW%value(dRdW%nnz))
allocate(dRdW%col_pointer(4*mesh_cpx%ncell + 1))
dRdW%row(:) = 0
dRdW%column(:) = 0
dRdW%value(:) = 0.0d0 
dRdW%col_pointer(:) = 0 

!evaluate the jacobian 
dRdW%col_pointer(1) = 1
nblock = (sum(mesh_cpx%cells_nadj) + mesh_cpx%ncell)*4 
ncolour = maxval(mesh_cpx%cells_colour)
do clr=1,ncolour
    if (options%cdisplay) then 
        write(*,'(A,I0,A,I0)') '    colour: ',clr,'/',ncolour
    end if 
    do vv=1,4 !perturb each conservative variable in all cells of this colour

        !perturb variables 
        !$OMP parallel do 
        do cc=1,mesh_cpx%ncell
            if (mesh_cpx%cells_colour(cc) == clr) then 
                if (vv == 1) then 
                    mesh_cpx%w1(cc) = mesh_cpx%w1(cc) + complex(0.0d0,h)
                elseif (vv == 2) then 
                    mesh_cpx%w2(cc) = mesh_cpx%w2(cc) + complex(0.0d0,h)
                elseif (vv == 3) then 
                    mesh_cpx%w3(cc) = mesh_cpx%w3(cc) + complex(0.0d0,h)
                elseif (vv == 4) then 
                    mesh_cpx%w4(cc) = mesh_cpx%w4(cc) + complex(0.0d0,h)
                end if 
            end if  
        end do 
        !$OMP end parallel do 
        !$OMP parallel do 
        do cc=1,mesh_cpx%ncell
            call con2prim_cpx(mesh_cpx%rho(cc),mesh_cpx%u(cc),mesh_cpx%v(cc),mesh_cpx%p(cc),mesh_cpx%e(cc),complex(options%gamma,0.0d0),mesh_cpx%w1(cc),mesh_cpx%w2(cc),mesh_cpx%w3(cc),mesh_cpx%w4(cc))
        end do 
        !$OMP end parallel do 

        !evaluate the residual 
        !$OMP parallel
        call get_edge_fluxes_cpx(mesh_cpx,options,.False.)
        !$OMP end parallel 
        !$OMP parallel do private(ee)
        do cc=1,mesh_cpx%ncell
            pr1c(cc) = complex(0.0d0,0.0d0) 
            pr2c(cc) = complex(0.0d0,0.0d0) 
            pr3c(cc) = complex(0.0d0,0.0d0) 
            pr4c(cc) = complex(0.0d0,0.0d0) 
            do ee=1,mesh_cpx%cells(cc)%nedge
                pr1c(cc) = pr1c(cc) + (mesh_cpx%edges_r1(mesh_cpx%cells(cc)%edge(ee)) + mesh_cpx%edges_d1(mesh_cpx%cells(cc)%edge(ee)))*mesh_cpx%cells(cc)%edge_sign(ee)
                pr2c(cc) = pr2c(cc) + (mesh_cpx%edges_r2(mesh_cpx%cells(cc)%edge(ee)) + mesh_cpx%edges_d2(mesh_cpx%cells(cc)%edge(ee)))*mesh_cpx%cells(cc)%edge_sign(ee)
                pr3c(cc) = pr3c(cc) + (mesh_cpx%edges_r3(mesh_cpx%cells(cc)%edge(ee)) + mesh_cpx%edges_d3(mesh_cpx%cells(cc)%edge(ee)))*mesh_cpx%cells(cc)%edge_sign(ee)
                pr4c(cc) = pr4c(cc) + (mesh_cpx%edges_r4(mesh_cpx%cells(cc)%edge(ee)) + mesh_cpx%edges_d4(mesh_cpx%cells(cc)%edge(ee)))*mesh_cpx%cells(cc)%edge_sign(ee)
            end do 
        end do
        !$OMP end parallel do 

        !extract non-zero values and populate the flow jacobian 
        do cc=1,mesh_cpx%ncell
            if (mesh_cpx%cells_colour(cc) == clr) then 
                
                !reset the column index offset
                column_offset_index = 0 

                !get the location of the start of this column in the sparse structure
                r0 = (sum(mesh_cpx%cells_nadj(1:cc-1)) + cc - 1)*4 + nblock*(vv - 1) + 1

                !get the column corrsponding to this cell cc and this conservative variable vv
                col = cc + (vv - 1)*mesh_cpx%ncell

                !extract each residual for this cell
                do rr=1,4

                    !get the row index of this entry 
                    row = cc + (rr - 1)*mesh_cpx%ncell
                    
                    !get the col_offset of the entry for this cell in this column
                    col_offset = column_offset_index
                    column_offset_index = column_offset_index + 1

                    !insert the values
                    if (rr == 1) then 
                        dRdW%value(r0 + col_offset) = aimag(pr1c(cc))/h
                    elseif (rr == 2) then 
                        dRdW%value(r0 + col_offset) = aimag(pr2c(cc))/h
                    elseif (rr == 3) then 
                        dRdW%value(r0 + col_offset) = aimag(pr3c(cc))/h
                    elseif (rr == 4) then 
                        dRdW%value(r0 + col_offset) = aimag(pr4c(cc))/h
                    end if 
                    dRdW%column(r0 + col_offset) = col
                    dRdW%row(r0 + col_offset) = row
                end do 

                !extract each residual for this cells adjacent cells
                do aa=1,mesh_cpx%cells(cc)%nedge

                    !get the adjacent cell 
                    if (mesh_cpx%edges(mesh_cpx%cells(cc)%edge(aa))%c1 .NE. cc) then 
                        cadj = mesh_cpx%edges(mesh_cpx%cells(cc)%edge(aa))%c1 
                    else
                        cadj = mesh_cpx%edges(mesh_cpx%cells(cc)%edge(aa))%c2 
                    end if 
                    if (cadj .LT. 0) then !skip if boundary condition 
                        cycle 
                    end if 

                    !extract each residual for this cell
                    do rr=1,4
                        
                        !get the row index of this entry 
                        row = cadj + (rr - 1)*mesh_cpx%ncell

                        !get the col_offset of the entry for this cell in this column
                        col_offset = column_offset_index
                        column_offset_index = column_offset_index + 1

                        !insert the values
                        if (rr == 1) then 
                            dRdW%value(r0 + col_offset) = aimag(pr1c(cadj))/h
                        elseif (rr == 2) then 
                            dRdW%value(r0 + col_offset) = aimag(pr2c(cadj))/h
                        elseif (rr == 3) then 
                            dRdW%value(r0 + col_offset) = aimag(pr3c(cadj))/h
                        elseif (rr == 4) then 
                            dRdW%value(r0 + col_offset) = aimag(pr4c(cadj))/h
                        end if 
                        dRdW%column(r0 + col_offset) = col
                        dRdW%row(r0 + col_offset) = row
                    end do 
                end do 
                
                !set the column end pointer 
                dRdW%col_pointer(col+1) = r0 + col_offset + 1
            end if 
        end do 

        !reset variables 
        !$OMP parallel do 
        do cc=1,mesh_cpx%ncell
            if (mesh_cpx%cells_colour(cc) == clr) then 
                if (vv == 1) then 
                    mesh_cpx%w1(cc) = w10(cc)
                elseif (vv == 2) then 
                    mesh_cpx%w2(cc) = w20(cc)
                elseif (vv == 3) then 
                    mesh_cpx%w3(cc) = w30(cc)
                elseif (vv == 4) then 
                    mesh_cpx%w4(cc) = w40(cc)
                end if 
            end if  
        end do 
        !$OMP end parallel do 
    end do 
end do 

!reset the primative variables 
mesh_cpx%w1 = w10 
mesh_cpx%w2 = w20 
mesh_cpx%w3 = w30 
mesh_cpx%w4 = w40 
!$OMP parallel do 
do cc=1,mesh_cpx%ncell
    call con2prim_cpx(mesh_cpx%rho(cc),mesh_cpx%u(cc),mesh_cpx%v(cc),mesh_cpx%p(cc),mesh_cpx%e(cc),complex(options%gamma,0.0d0),mesh_cpx%w1(cc),mesh_cpx%w2(cc),mesh_cpx%w3(cc),mesh_cpx%w4(cc))
end do 
!$OMP end parallel do 
return 
end subroutine build_flow_jacobian_sparse

!sparse grid jacobian ===============
subroutine build_grid_jacobian_sparse(mesh_cpx,options,dRdX)
implicit none

!variables - inout
type(flux_mesh_cpx) :: mesh_cpx 
type(flux_options) :: options 
type(csc_matrix) :: dRdX

!variables - local 
integer(in32) :: clr,cc,vv,ee,rr,aa
integer(in32) :: ncolour,r0,col_offset,row,col,cadj,nblock
integer(in32) :: column_offset_index
real(dp) :: h
complex(dp) :: pr1c(mesh_cpx%ncell),pr2c(mesh_cpx%ncell),pr3c(mesh_cpx%ncell),pr4c(mesh_cpx%ncell)
complex(dp) :: coordinatetes0_x(mesh_cpx%nvertex),coordinatetes0_y(mesh_cpx%nvertex)

!define the complex stepsize
h = 1e-40_dp 

!initialise dissipation
mesh_cpx%edges_d1(:) = complex(0.0d0,0.0d0) 
mesh_cpx%edges_d2(:) = complex(0.0d0,0.0d0) 
mesh_cpx%edges_d3(:) = complex(0.0d0,0.0d0) 
mesh_cpx%edges_d4(:) = complex(0.0d0,0.0d0) 

!store the initial coordinates
do cc=1,mesh_cpx%nvertex
    coordinatetes0_x(cc) = mesh_cpx%vertices(cc)%coordinate(1)
    coordinatetes0_y(cc) = mesh_cpx%vertices(cc)%coordinate(2)
end do 

!allocate the sparse flow jacobian 
dRdX%nnz = sum(mesh_cpx%vertices%ncell)*8
if (options%cdisplay) then
    write(*,'(A,I0,A)') '    {grid jacobian nnz = ',dRdX%nnz,'}'
    write(*,'(A,F8.6,A)') '    {grid jacobian sparsity = ',(real(dRdX%nnz,dp)/real(8_in64*mesh_cpx%ncell*mesh_cpx%nvertex,dp))*100.0d0,'%}'
end if 
dRdX%nrow = 4*mesh_cpx%ncell
dRdX%ncol = 2*mesh_cpx%nvertex
allocate(dRdX%row(dRdX%nnz))
allocate(dRdX%column(dRdX%nnz))
allocate(dRdX%value(dRdX%nnz))
allocate(dRdX%col_pointer(2*mesh_cpx%nvertex + 1))
dRdX%row(:) = 0
dRdX%column(:) = 0
dRdX%value(:) = 0.0d0 
dRdX%col_pointer(:) = 0 
 
!evaluate the jacobian 
dRdX%col_pointer(1) = 1
nblock = sum(mesh_cpx%vertices%ncell)*4 
ncolour = maxval(mesh_cpx%vertices_colour)
do clr=1,ncolour
    if (options%cdisplay) then 
        write(*,'(A,I0,A,I0)') '    colour: ',clr,'/',ncolour
    end if 
    do vv=1,2 !perturb each coordinate at vertices of this colour

        !perturb variables 
        !$OMP parallel do 
        do cc=1,mesh_cpx%nvertex
            if (mesh_cpx%vertices_colour(cc) == clr) then 
                if (vv == 1) then 
                    mesh_cpx%vertices(cc)%coordinate(1) = mesh_cpx%vertices(cc)%coordinate(1) + complex(0.0d0,h)
                elseif (vv == 2) then 
                    mesh_cpx%vertices(cc)%coordinate(2) = mesh_cpx%vertices(cc)%coordinate(2) + complex(0.0d0,h)
                end if 
            end if  
        end do 
        !$OMP end parallel do 

        !update edge geometries 
        call mesh_cpx%get_edges_geometry()

        !evaluate the residual 
        !$OMP parallel
        call get_edge_fluxes_cpx(mesh_cpx,options,.False.)
        !$OMP end parallel
        !$OMP parallel do private(ee)
        do cc=1,mesh_cpx%ncell
            pr1c(cc) = complex(0.0d0,0.0d0) 
            pr2c(cc) = complex(0.0d0,0.0d0) 
            pr3c(cc) = complex(0.0d0,0.0d0) 
            pr4c(cc) = complex(0.0d0,0.0d0) 
            do ee=1,mesh_cpx%cells(cc)%nedge
                pr1c(cc) = pr1c(cc) + (mesh_cpx%edges_r1(mesh_cpx%cells(cc)%edge(ee)) + mesh_cpx%edges_d1(mesh_cpx%cells(cc)%edge(ee)))*mesh_cpx%cells(cc)%edge_sign(ee)
                pr2c(cc) = pr2c(cc) + (mesh_cpx%edges_r2(mesh_cpx%cells(cc)%edge(ee)) + mesh_cpx%edges_d2(mesh_cpx%cells(cc)%edge(ee)))*mesh_cpx%cells(cc)%edge_sign(ee)
                pr3c(cc) = pr3c(cc) + (mesh_cpx%edges_r3(mesh_cpx%cells(cc)%edge(ee)) + mesh_cpx%edges_d3(mesh_cpx%cells(cc)%edge(ee)))*mesh_cpx%cells(cc)%edge_sign(ee)
                pr4c(cc) = pr4c(cc) + (mesh_cpx%edges_r4(mesh_cpx%cells(cc)%edge(ee)) + mesh_cpx%edges_d4(mesh_cpx%cells(cc)%edge(ee)))*mesh_cpx%cells(cc)%edge_sign(ee)
            end do 
        end do
        !$OMP end parallel do 

        !extract non-zero values and populate the grid jacobian 
        do cc=1,mesh_cpx%nvertex
            if (mesh_cpx%vertices_colour(cc) == clr) then 

                !reset the column index offset
                column_offset_index = 0 

                !get the location of the start of this column in the sparse structure
                r0 = sum(mesh_cpx%vertices(1:cc-1)%ncell)*4 + nblock*(vv - 1) + 1

                !get the column corrsponding to this vertex cc and this coordinate vv
                col = cc + (vv - 1)*mesh_cpx%nvertex

                !extract each residual for this vertices adjacent cells
                do aa=1,mesh_cpx%vertices(cc)%ncell

                    !get the adjacent cell 
                    cadj = mesh_cpx%vertices(cc)%cells(aa)

                    !extract each residual for this cell
                    do rr=1,4
                        
                        !get the row index of this entry 
                        row = cadj + (rr - 1)*mesh_cpx%ncell

                        !get the col_offset of the entry for this cell in this column
                        col_offset = column_offset_index
                        column_offset_index = column_offset_index + 1

                        !insert the values
                        if (rr == 1) then 
                            dRdX%value(r0 + col_offset) = aimag(pr1c(cadj))/h
                        elseif (rr == 2) then 
                            dRdX%value(r0 + col_offset) = aimag(pr2c(cadj))/h
                        elseif (rr == 3) then 
                            dRdX%value(r0 + col_offset) = aimag(pr3c(cadj))/h
                        elseif (rr == 4) then 
                            dRdX%value(r0 + col_offset) = aimag(pr4c(cadj))/h
                        end if 
                        dRdX%column(r0 + col_offset) = col
                        dRdX%row(r0 + col_offset) = row
                    end do 
                end do 
                
                !set the column end pointer 
                dRdX%col_pointer(col+1) = r0 + col_offset + 1
            end if 
        end do 

        !reset variables 
        !$OMP parallel do 
        do cc=1,mesh_cpx%nvertex
            if (mesh_cpx%vertices_colour(cc) == clr) then 
                if (vv == 1) then 
                    mesh_cpx%vertices(cc)%coordinate(1) = coordinatetes0_x(cc)
                elseif (vv == 2) then 
                    mesh_cpx%vertices(cc)%coordinate(2) = coordinatetes0_y(cc)
                end if 
            end if  
        end do 
        !$OMP end parallel do 
    end do 
end do 
return 
end subroutine build_grid_jacobian_sparse

!rk iterate ===============
subroutine rk_iterate_adj(dRdW_sp,dJdW,mesh,options,cell_timestep,nanflag)
implicit none 

!variables - inout
logical :: nanflag
type(flux_mesh) :: mesh 
type(flux_options) :: options 
type(csc_matrix) :: dRdW_sp
real(dp) :: dJdW(4*mesh%ncell),cell_timestep(4*mesh%ncell)

!variables - local
integer(in32) :: rr,cc

!store the initial psi
!$OMP single
mesh%psi0 = mesh%psi 
!$OMP end single

!iterate
do rr=1,options%rk_niter

    !evaluate the residual 
    call csc_vector_product(dRdW_sp,mesh%psi,mesh%psi_dRdw_prd)

    !evaluate dissipation if enabled



    !step each cell 
    !$OMP do schedule(guided,50) 
    do cc=1,4*mesh%ncell
        mesh%psi(cc) = mesh%psi0(cc) - rk4_alpha(rr)*cell_timestep(cc)*(mesh%psi_dRdw_prd(cc) - dJdW(cc))
        if (isnan(mesh%psi(cc))) then 
            print *, 'cell nan: ',cc
            nanflag = .true.
        end if 
    end do 
    !$OMP end do 
end do 

!store the mesh residual
!$OMP single
mesh%residual = mesh%psi_dRdw_prd(1:mesh%ncell) - dJdW(1:mesh%ncell)
!$OMP end single
return 
end subroutine rk_iterate_adj

!csc vector product ===============
subroutine csc_vector_product(matrix,vector,result)
implicit none 

!variables inout
type(csc_matrix) :: matrix
real(dp) :: vector(matrix%nrow),result(matrix%ncol)

!variables local 
integer(in64) :: ii,jj

!evaluate
!$OMP do schedule(guided,50) private(jj)
do ii=1,matrix%ncol
    result(ii) = 0.0d0 
    do jj=matrix%col_pointer(ii),matrix%col_pointer(ii+1)-1
        result(ii) = result(ii) + matrix%value(jj)*vector(matrix%row(jj))
    end do 
end do 
!$OMP end do 
return 
end subroutine csc_vector_product

!construct complex mesh
subroutine construct_complex_mesh(mesh_cpx,mesh,options)
implicit none 

!variables - inout
type(flux_mesh) :: mesh 
type(flux_mesh_cpx) :: mesh_cpx 
type(flux_options) :: options 

!variables - local 
integer(in32) :: ii,jj
integer(in32) :: cindex,nedge

!allocate the complex mesh 
mesh_cpx%nvertex = mesh%nvertex
mesh_cpx%nedge = mesh%nedge
mesh_cpx%ncell = mesh%ncell

!allocate cell based properties
allocate(mesh_cpx%cells(mesh%ncell))
allocate(mesh_cpx%cells_specrad(mesh%ncell))
allocate(mesh_cpx%cells_volume(mesh%ncell))
allocate(mesh_cpx%rho(mesh%ncell))
allocate(mesh_cpx%u(mesh%ncell))
allocate(mesh_cpx%v(mesh%ncell))
allocate(mesh_cpx%p(mesh%ncell))
allocate(mesh_cpx%mach(mesh%ncell))
allocate(mesh_cpx%e(mesh%ncell))
allocate(mesh_cpx%cp(mesh%ncell))
allocate(mesh_cpx%w1(mesh%ncell))
allocate(mesh_cpx%w2(mesh%ncell))
allocate(mesh_cpx%w3(mesh%ncell))
allocate(mesh_cpx%w4(mesh%ncell))
allocate(mesh_cpx%w10(mesh%ncell))
allocate(mesh_cpx%w20(mesh%ncell))
allocate(mesh_cpx%w30(mesh%ncell))
allocate(mesh_cpx%w40(mesh%ncell))
allocate(mesh_cpx%cells_dt(mesh%ncell))
allocate(mesh_cpx%l1(mesh%ncell))
allocate(mesh_cpx%l2(mesh%ncell))
allocate(mesh_cpx%l3(mesh%ncell))
allocate(mesh_cpx%l4(mesh%ncell))
allocate(mesh_cpx%cells_psensor(mesh%ncell))
allocate(mesh_cpx%residual(mesh%ncell))
allocate(mesh_cpx%cells_colour(mesh%ncell))
allocate(mesh_cpx%cells_nadj(mesh%ncell))
do ii=1,mesh%ncell
    cindex = mesh%cells(ii)%index
    nedge = mesh%cells(cindex)%nedge
    mesh_cpx%cells(cindex)%nedge = nedge
    mesh_cpx%cells(cindex)%index = cindex
    allocate(mesh_cpx%cells(cindex)%edgev1(nedge))
    allocate(mesh_cpx%cells(cindex)%edgev2(nedge))
    allocate(mesh_cpx%cells(cindex)%edgec(nedge))
    allocate(mesh_cpx%cells(cindex)%edge(nedge))    
    allocate(mesh_cpx%cells(cindex)%edge_sign(nedge))
    do jj=1,nedge
        mesh_cpx%cells(cindex)%edgev1(jj) = mesh%cells(cindex)%edgev1(jj)
        mesh_cpx%cells(cindex)%edgev2(jj) = mesh%cells(cindex)%edgev2(jj)
        mesh_cpx%cells(cindex)%edgec(jj) = mesh%cells(cindex)%edgec(jj)
        mesh_cpx%cells(cindex)%edge(jj) = mesh%cells(cindex)%edge(jj)
    end do 
end do 

!allocate edge based properties
allocate(mesh_cpx%edges(mesh%nedge))
allocate(mesh_cpx%edges_specrad(mesh%nedge))
allocate(mesh_cpx%edges_r1(mesh%nedge))
allocate(mesh_cpx%edges_r2(mesh%nedge))
allocate(mesh_cpx%edges_r3(mesh%nedge))
allocate(mesh_cpx%edges_r4(mesh%nedge))
allocate(mesh_cpx%edges_l1(mesh%nedge))
allocate(mesh_cpx%edges_l2(mesh%nedge))
allocate(mesh_cpx%edges_l3(mesh%nedge))
allocate(mesh_cpx%edges_l4(mesh%nedge))
allocate(mesh_cpx%edges_d1(mesh%nedge))
allocate(mesh_cpx%edges_d2(mesh%nedge))
allocate(mesh_cpx%edges_d3(mesh%nedge))
allocate(mesh_cpx%edges_d4(mesh%nedge))
allocate(mesh_cpx%edges_pn(mesh%nedge))
allocate(mesh_cpx%edges_pd(mesh%nedge))
do ii=1,mesh%nedge
    mesh_cpx%edges(ii)%index = ii 
    mesh_cpx%edges(ii)%v1 = mesh%edges(ii)%v1
    mesh_cpx%edges(ii)%v2 = mesh%edges(ii)%v2
    mesh_cpx%edges(ii)%c1 = mesh%edges(ii)%c1
    mesh_cpx%edges(ii)%c2 = mesh%edges(ii)%c2
end do 

!set the edge signs for each cell
do ii=1,mesh%ncell
    do jj=1,mesh%cells(ii)%nedge
        mesh_cpx%cells(ii)%edge_sign(jj) = complex(mesh%cells(ii)%edge_sign(jj),0.0d0)
    end do 
end do 

!allocate vertices
allocate(mesh_cpx%vertices(mesh%nvertex))
do ii=1,mesh%nvertex
    mesh_cpx%vertices(ii)%coordinate(1) = complex(mesh%vertices(ii)%coordinate(1),0.0d0)
    mesh_cpx%vertices(ii)%coordinate(2) = complex(mesh%vertices(ii)%coordinate(2),0.0d0)
end do

!assign colours and adjoint properties
mesh_cpx%cells_colour = mesh%cells_colour
mesh_cpx%vertices_colour = mesh%vertices_colour
mesh_cpx%cells_nadj = mesh%cells_nadj
do ii=1,mesh%nvertex
    mesh_cpx%vertices(ii)%ncell = mesh%vertices(ii)%ncell
    allocate(mesh_cpx%vertices(ii)%cells(mesh%vertices(ii)%ncell))
    mesh_cpx%vertices(ii)%cells = mesh%vertices(ii)%cells
end do 

!set the primative conditions in each cell 
do ii=1,mesh%ncell
    mesh_cpx%rho(ii) = complex(mesh%rho(ii),0.0d0)
    mesh_cpx%u(ii) = complex(mesh%u(ii),0.0d0)
    mesh_cpx%v(ii) = complex(mesh%v(ii),0.0d0)
    mesh_cpx%p(ii) = complex(mesh%p(ii),0.0d0)
    mesh_cpx%mach(ii) = complex(mesh%mach(ii),0.0d0)
    mesh_cpx%e(ii) = complex(mesh%e(ii),0.0d0)
    mesh_cpx%cp(ii) = complex(mesh%cp(ii),0.0d0)
end do 

!set the conservative variables in each cell 
do ii=1,mesh_cpx%ncell
    call prim2con_cpx(mesh_cpx%rho(ii),mesh_cpx%u(ii),mesh_cpx%v(ii),mesh_cpx%p(ii),complex(options%gamma,0.0d0),mesh_cpx%w1(ii),mesh_cpx%w2(ii),mesh_cpx%w3(ii),mesh_cpx%w4(ii))
end do 

!get the mesh properties 
call mesh_cpx%get_edges_geometry()
return 
end subroutine construct_complex_mesh

!flux adjoint solve ===============
subroutine flux_adjoint_solve(mesh,options)
implicit none 

!variables - inout
type(flux_mesh) :: mesh 
type(flux_options) :: options 

!variables - local 
logical :: nanflag,resconv
integer(in32) :: iteration,cc,ee
real(dp) :: psirhores,psirhores0
real(dp) :: cell_timestep(4*mesh%ncell)
real(dp), dimension(:), allocatable :: dJdW,dJdX,dtotal
type(csc_matrix) :: dRdW_sp,dRdX_sp
type(flux_mesh_cpx) :: mesh_cpx 

!initialise flags
resconv = .false.
nanflag = .false.

!set the number of threads 
call omp_set_num_threads(options%num_threads)

!get the number of adjacent cells for each mesh cell
call mesh%get_cells_nadj()

!get the vertex connectivity
call mesh%get_vertex_cells()

!colour the cells
if (options%cdisplay) then
    write(*,'(A)') '--> colouring cells'
end if 
call colour_cells(mesh,options)

!colour the vertices
if (options%cdisplay) then
    write(*,'(A)') '--> colouring vertices'
end if 
call colour_vertices(mesh,options)

!construct the complex mesh 
if (options%cdisplay) then
    write(*,'(A)') '--> casting to complex variables '
end if 
call construct_complex_mesh(mesh_cpx,mesh,options)

!evaluate the flow jacobian 
if (options%cdisplay) then
    write(*,'(A)') '--> evaluating the flow jacobian'
end if 
call build_flow_jacobian_sparse(mesh_cpx,options,dRdW_sp)

!evaluate the flow objective jacobian 
if (options%cdisplay) then
    write(*,'(A,A)') '--> evaluating the objective gradient: ',options%adjoint_objective
end if 
if (options%adjoint_objective == 'cd') then 
    call cd_objective(mesh,options,dJdW,dJdX)
elseif (options%adjoint_objective == 'mflux') then 
    call mflux_objective(mesh,options,dJdW,dJdX)
else
    write(*,'(A,A)') '** invalid adjoint objective specified: ',options%adjoint_objective
    stop
end if 

!evalaute the adjoint timesteps for each cell
if (options%cdisplay) then
    write(*,'(A)') '--> evaluating cell adjoint timesteps'
end if 
do ee=1,mesh%nedge
    call edge_spectral_radius(mesh,options,ee)
end do 
do cc=1,mesh%ncell
    mesh%cells_specrad(cc) = 0.0d0 
    do ee=1,mesh%cells(cc)%nedge
        mesh%cells_specrad(cc) = mesh%cells_specrad(cc) + mesh%edges_specrad(mesh%cells(cc)%edge(ee))
    end do
    mesh%cells_dt(cc) = options%cfl*(mesh%cells_volume(cc)/mesh%cells_specrad(cc))
end do 
do cc=1,mesh%ncell
    cell_timestep(cc) = (mesh%cells_dt(cc)/mesh%cells_volume(cc))
    cell_timestep(cc + mesh%ncell) = cell_timestep(cc)
    cell_timestep(cc + 2*mesh%ncell) = cell_timestep(cc)
    cell_timestep(cc + 3*mesh%ncell) = cell_timestep(cc)
end do 

!display
if (options%cdisplay) then
    write(*,'(A)') '--> solving'
end if 

!initialise psi
allocate(mesh%psi0(4*mesh%ncell))
allocate(mesh%psi(4*mesh%ncell))
allocate(mesh%psi_dRdw_prd(4*mesh%ncell))
mesh%psi0 = 0.0d0 
mesh%psi = 0.0d0 
mesh%psi_dRdw_prd = 0.0d0 

!initial display
if (options%cdisplay) then     
    write(*,'(A)') '    ittn    |  psirhores   '
    write(*,'(A)') '+-----------+-------------+'
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

    !rk itterate
    call rk_iterate_adj(dRdW_sp,dJdW,mesh,options,cell_timestep,nanflag)

    !initiate single thread section
    !$OMP single    

    !check convergence
    if (iteration == 1) then 
        psirhores0 = sqrt(sum((mesh%residual)**2))
    end if 
    psirhores = log10(sqrt(sum((mesh%residual)**2))/psirhores0)
    mesh%psires = psirhores
    if (psirhores .LE. options%residual_convtol) then 
        resconv = .true.
    end if 

    !display 
    if (options%cdisplay) then     
        if (mod(iteration,100) == 0) then
            write(*,'(A)') '    ittn    |  psirhores   '
            write(*,'(A)') '+-----------+-------------+'
        end if 
            write(*,'(A,I8,A,F11.8)') '    ',iteration,'  ',psirhores
    end if

    !end single thread section
    !$OMP end single
end do 

!end the parallel region
!$OMP end parallel 

!extract the results to the mesh variables
allocate(mesh%psi_rho(mesh%ncell))
allocate(mesh%psi_u(mesh%ncell))
allocate(mesh%psi_v(mesh%ncell))
allocate(mesh%psi_e(mesh%ncell))
if (nanflag) then 
    mesh%psi = ieee_value(1.0d0,ieee_quiet_nan)
end if 
mesh%psi_rho = mesh%psi(1:mesh%ncell)
mesh%psi_u = mesh%psi(mesh%ncell+1:2*mesh%ncell)
mesh%psi_v = mesh%psi(2*mesh%ncell+1:3*mesh%ncell)
mesh%psi_e = mesh%psi(3*mesh%ncell+1:4*mesh%ncell)

!evaluate the grid jacobian 
call build_grid_jacobian_sparse(mesh_cpx,options,dRdX_sp)

!evaluate the total derivative
allocate(dtotal(2*mesh%nvertex))
!$OMP parallel
call csc_vector_product(dRdX_sp,mesh%psi,dtotal)
!$OMP end parallel 
dtotal = dtotal - dJdX
dtotal = -dtotal

!extract the total derivative
allocate(mesh%vertex_derivative_x(mesh%nvertex))
allocate(mesh%vertex_derivative_y(mesh%nvertex))
mesh%vertex_derivative_x = dtotal(1:mesh%nvertex)
mesh%vertex_derivative_y = dtotal(mesh%nvertex+1:2*mesh%nvertex)
return 
end subroutine flux_adjoint_solve

end module flux_adjoint
