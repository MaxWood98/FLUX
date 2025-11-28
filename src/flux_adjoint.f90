!flux 2d adjoint module 
!max wood
!version : 0.0.1
!updated : 27-11-25

module flux_adjoint
use flux_io
use flux_solve
use flux_data_methods
contains 

!dCddW ===============
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
! call con2prim(rhop,up,vp,pp,ep,options%gamma,w10+1e-8,w20,w30,w40)
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
subroutine build_flow_jacobian_sparse(mesh,options,dRdW)
implicit none

!variables - inout
type(flux_mesh) :: mesh 
type(flux_options) :: options 
type(csc_matrix) :: dRdW

!variables - local 
integer(in32) :: clr,cc,vv,ee,rr,aa
integer(in32) :: ncolour,r0,col_offset,row,col,cadj,nblock
integer(in32) :: column_offset_index

real(dp) :: r1,r2,r3,r4,h
real(dp) :: w10(mesh%ncell),w20(mesh%ncell),w30(mesh%ncell),w40(mesh%ncell)
real(dp) :: r10(mesh%ncell),r20(mesh%ncell),r30(mesh%ncell),r40(mesh%ncell)
real(dp) :: pr1(mesh%ncell),pr2(mesh%ncell),pr3(mesh%ncell),pr4(mesh%ncell)



h = 1e-6_dp 



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








!allocate the sparse flow jacobian 
dRdW%nnz = (sum(mesh%cells_nadj) + mesh%ncell)*16 
if (options%cdisplay) then
    write(*,'(A,I0,A)') '    {flow jacobian nnz = ',dRdW%nnz,'}'
    write(*,'(A,F8.6,A)') '    {flow jacobian sparsity = ',(real(dRdW%nnz,dp)/real(16*mesh%ncell*mesh%ncell,dp))*100.0d0,'%}'
end if 
dRdW%nrow = 4*mesh%ncell
dRdW%ncol = 4*mesh%ncell
allocate(dRdW%row(dRdW%nnz))
allocate(dRdW%column(dRdW%nnz))
allocate(dRdW%value(dRdW%nnz))
allocate(dRdW%col_pointer(4*mesh%ncell + 1))
dRdW%row(:) = 0
dRdW%column(:) = 0
dRdW%value(:) = 0.0d0 
dRdW%col_pointer(:) = 0 

!evaluate the jacobian 
dRdW%col_pointer(1) = 1
nblock = (sum(mesh%cells_nadj) + mesh%ncell)*4 
ncolour = maxval(mesh%cells_colour)
do clr=1,ncolour
    write(*,'(A,I0,A,I0)') '    colour: ',clr,'/',ncolour
    do vv=1,4 !perturb each conservative variable in all cells of this colour



        !update to complex step ==========================

        !perturb variables 
        mesh%w1 = w10 
        mesh%w2 = w20 
        mesh%w3 = w30 
        mesh%w4 = w40 
        do cc=1,mesh%ncell
            if (mesh%cells_colour(cc) == clr) then 
                if (vv == 1) then 
                    mesh%w1(cc) = mesh%w1(cc) + h
                elseif (vv == 2) then 
                    mesh%w2(cc) = mesh%w2(cc) + h
                elseif (vv == 3) then 
                    mesh%w3(cc) = mesh%w3(cc) + h
                elseif (vv == 4) then 
                    mesh%w4(cc) = mesh%w4(cc) + h
                end if 
            end if  
        end do 
        do cc=1,mesh%ncell
            call con2prim(mesh%rho(cc),mesh%u(cc),mesh%v(cc),mesh%p(cc),mesh%e(cc),options%gamma,mesh%w1(cc),mesh%w2(cc),mesh%w3(cc),mesh%w4(cc))
        end do 

        !evaluate the flow residual
        call get_edge_fluxes(mesh,options,.False.)
        do cc=1,mesh%ncell
            pr1(cc) = 0.0d0 
            pr2(cc) = 0.0d0 
            pr3(cc) = 0.0d0 
            pr4(cc) = 0.0d0 
            do ee=1,mesh%cells(cc)%nedge
                pr1(cc) = pr1(cc) + (mesh%edges_r1(mesh%cells(cc)%edge(ee)) + mesh%edges_d1(mesh%cells(cc)%edge(ee)))*mesh%cells(cc)%edge_sign(ee)
                pr2(cc) = pr2(cc) + (mesh%edges_r2(mesh%cells(cc)%edge(ee)) + mesh%edges_d2(mesh%cells(cc)%edge(ee)))*mesh%cells(cc)%edge_sign(ee)
                pr3(cc) = pr3(cc) + (mesh%edges_r3(mesh%cells(cc)%edge(ee)) + mesh%edges_d3(mesh%cells(cc)%edge(ee)))*mesh%cells(cc)%edge_sign(ee)
                pr4(cc) = pr4(cc) + (mesh%edges_r4(mesh%cells(cc)%edge(ee)) + mesh%edges_d4(mesh%cells(cc)%edge(ee)))*mesh%cells(cc)%edge_sign(ee)
            end do 
        end do

        !update to complex step ==========================




        !extract non-zero values and populate the flow jacobian 
        do cc=1,mesh%ncell
            if (mesh%cells_colour(cc) == clr) then 
                
                !reset the column index offset
                column_offset_index = 0 

                !get the location of the start of this column in the sparse structure
                r0 = (sum(mesh%cells_nadj(1:cc-1)) + cc - 1)*4 + nblock*(vv - 1) + 1

                !get the column corrsponding to this cell cc and this conservative variable vv
                col = cc + (vv - 1)*mesh%ncell

                !extract each residual for this cell
                do rr=1,4

                    !get the row index of this entry 
                    row = cc + (rr - 1)*mesh%ncell
                    
                    !get the col_offset of the entry for this cell in this column
                    col_offset = column_offset_index
                    column_offset_index = column_offset_index + 1

                    !insert the values
                    if (rr == 1) then 
                        dRdW%value(r0 + col_offset) = (pr1(cc) - r10(cc))/h
                    elseif (rr == 2) then 
                        dRdW%value(r0 + col_offset) = (pr2(cc) - r20(cc))/h
                    elseif (rr == 3) then 
                        dRdW%value(r0 + col_offset) = (pr3(cc) - r30(cc))/h
                    elseif (rr == 4) then 
                        dRdW%value(r0 + col_offset) = (pr4(cc) - r40(cc))/h
                    end if 
                    dRdW%column(r0 + col_offset) = col
                    dRdW%row(r0 + col_offset) = row
                end do 

                !extract each residual for this cells adjacent cells
                do aa=1,mesh%cells(cc)%nedge

                    !get the adjacent cell 
                    if (mesh%edges(mesh%cells(cc)%edge(aa))%c1 .NE. cc) then 
                        cadj = mesh%edges(mesh%cells(cc)%edge(aa))%c1 
                    else
                        cadj = mesh%edges(mesh%cells(cc)%edge(aa))%c2 
                    end if 
                    if (cadj .LT. 0) then !skip if boundary condition 
                        cycle 
                    end if 

                    !extract each residual for this cell
                    do rr=1,4
                        
                        !get the row index of this entry 
                        row = cadj + (rr - 1)*mesh%ncell

                        !get the col_offset of the entry for this cell in this column
                        col_offset = column_offset_index
                        column_offset_index = column_offset_index + 1

                        !insert the values
                        if (rr == 1) then 
                            dRdW%value(r0 + col_offset) = (pr1(cadj) - r10(cadj))/h
                        elseif (rr == 2) then 
                            dRdW%value(r0 + col_offset) = (pr2(cadj) - r20(cadj))/h
                        elseif (rr == 3) then 
                            dRdW%value(r0 + col_offset) = (pr3(cadj) - r30(cadj))/h
                        elseif (rr == 4) then 
                            dRdW%value(r0 + col_offset) = (pr4(cadj) - r40(cadj))/h
                        end if 
                        dRdW%column(r0 + col_offset) = col
                        dRdW%row(r0 + col_offset) = row
                    end do 
                end do 
                
                !set the column end pointer 
                dRdW%col_pointer(col+1) = r0 + col_offset + 1
            end if 
        end do 
    end do 
end do 
return 
end subroutine build_flow_jacobian_sparse

!sparse grid jacobian ===============
subroutine build_grid_jacobian_sparse(mesh,options,dRdX)
implicit none

!variables - inout
type(flux_mesh) :: mesh 
type(flux_options) :: options 
type(csc_matrix) :: dRdX

!variables - local 
integer(in32) :: clr,cc,vv,ee,rr,aa
integer(in32) :: ncolour,r0,col_offset,row,col,cadj,nblock
integer(in32) :: column_offset_index

real(dp) :: r1,r2,r3,r4,h
real(dp) :: w10(mesh%ncell),w20(mesh%ncell),w30(mesh%ncell),w40(mesh%ncell)
real(dp) :: r10(mesh%ncell),r20(mesh%ncell),r30(mesh%ncell),r40(mesh%ncell)
real(dp) :: pr1(mesh%ncell),pr2(mesh%ncell),pr3(mesh%ncell),pr4(mesh%ncell)
real(dp) :: coordinates0(mesh%nvertex,2)


h = 1e-6_dp 



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

do vv=1,mesh%nvertex
    coordinates0(vv,:) = mesh%vertices(vv)%coordinate(:)
end do 



!allocate the sparse flow jacobian 
dRdX%nnz = sum(mesh%vertices%ncell)*8
if (options%cdisplay) then
    write(*,'(A,I0,A)') '    {grid jacobian nnz = ',dRdX%nnz,'}'
    write(*,'(A,F8.6,A)') '    {grid jacobian sparsity = ',(real(dRdX%nnz,dp)/real(8_in64*mesh%ncell*mesh%nvertex))*100.0d0,'%}'
end if 
dRdX%nrow = 4*mesh%ncell
dRdX%ncol = 2*mesh%nvertex
allocate(dRdX%row(dRdX%nnz))
allocate(dRdX%column(dRdX%nnz))
allocate(dRdX%value(dRdX%nnz))
allocate(dRdX%col_pointer(2*mesh%nvertex + 1))
dRdX%row(:) = 0
dRdX%column(:) = 0
dRdX%value(:) = 0.0d0 
dRdX%col_pointer(:) = 0 
 
!evaluate the jacobian 
dRdX%col_pointer(1) = 1
nblock = sum(mesh%vertices%ncell)*4 
ncolour = maxval(mesh%vertices_colour)
do clr=1,ncolour
    write(*,'(A,I0,A,I0)') '    colour: ',clr,'/',ncolour
    do vv=1,2 !perturb each coordinare at vertices of this colour



        !update to complex step ==========================

        !perturb coordinates
        do cc=1,mesh%nvertex
            mesh%vertices(cc)%coordinate = coordinates0(cc,:)
            if (mesh%vertices_colour(cc) == clr) then 
                if (vv == 1) then 
                    mesh%vertices(cc)%coordinate(1) = mesh%vertices(cc)%coordinate(1) + h
                elseif (vv == 2) then 
                    mesh%vertices(cc)%coordinate(2) = mesh%vertices(cc)%coordinate(2) + h
                end if 
            end if  
        end do 


        !update edge geometries 
        call mesh%get_edges_geometry()

        !set variables 
        mesh%w1 = w10 
        mesh%w2 = w20 
        mesh%w3 = w30 
        mesh%w4 = w40 
        do cc=1,mesh%ncell
            call con2prim(mesh%rho(cc),mesh%u(cc),mesh%v(cc),mesh%p(cc),mesh%e(cc),options%gamma,mesh%w1(cc),mesh%w2(cc),mesh%w3(cc),mesh%w4(cc))
        end do 

        !evaluate the flow residual
        call get_edge_fluxes(mesh,options,.False.)
        do cc=1,mesh%ncell
            pr1(cc) = 0.0d0 
            pr2(cc) = 0.0d0 
            pr3(cc) = 0.0d0 
            pr4(cc) = 0.0d0 
            do ee=1,mesh%cells(cc)%nedge
                pr1(cc) = pr1(cc) + (mesh%edges_r1(mesh%cells(cc)%edge(ee)) + mesh%edges_d1(mesh%cells(cc)%edge(ee)))*mesh%cells(cc)%edge_sign(ee)
                pr2(cc) = pr2(cc) + (mesh%edges_r2(mesh%cells(cc)%edge(ee)) + mesh%edges_d2(mesh%cells(cc)%edge(ee)))*mesh%cells(cc)%edge_sign(ee)
                pr3(cc) = pr3(cc) + (mesh%edges_r3(mesh%cells(cc)%edge(ee)) + mesh%edges_d3(mesh%cells(cc)%edge(ee)))*mesh%cells(cc)%edge_sign(ee)
                pr4(cc) = pr4(cc) + (mesh%edges_r4(mesh%cells(cc)%edge(ee)) + mesh%edges_d4(mesh%cells(cc)%edge(ee)))*mesh%cells(cc)%edge_sign(ee)
            end do 
        end do

        !update to complex step ==========================

        !extract non-zero values and populate the grid jacobian 
        do cc=1,mesh%nvertex
            if (mesh%vertices_colour(cc) == clr) then 


                ! r0 = (sum(mesh%cells_nadj(1:cc-1)) + cc - 1)*4 + nblock*(vv - 1) + 1
                
                !reset the column index offset
                column_offset_index = 0 

                !get the location of the start of this column in the sparse structure
                r0 = sum(mesh%vertices(1:cc-1)%ncell)*4 + nblock*(vv - 1) + 1

                !get the column corrsponding to this vertex cc and this coordinate vv
                col = cc + (vv - 1)*mesh%nvertex

                !extract each residual for this vertices adjacent cells
                do aa=1,mesh%vertices(cc)%ncell

                    !get the adjacent cell 
                    cadj = mesh%vertices(cc)%cells(aa)

                    !extract each residual for this cell
                    do rr=1,4
                        
                        !get the row index of this entry 
                        row = cadj + (rr - 1)*mesh%ncell

                        !get the col_offset of the entry for this cell in this column
                        col_offset = column_offset_index
                        column_offset_index = column_offset_index + 1

                        !insert the values
                        if (rr == 1) then 
                            dRdX%value(r0 + col_offset) = (pr1(cadj) - r10(cadj))/h
                        elseif (rr == 2) then 
                            dRdX%value(r0 + col_offset) = (pr2(cadj) - r20(cadj))/h
                        elseif (rr == 3) then 
                            dRdX%value(r0 + col_offset) = (pr3(cadj) - r30(cadj))/h
                        elseif (rr == 4) then 
                            dRdX%value(r0 + col_offset) = (pr4(cadj) - r40(cadj))/h
                        end if 
                        dRdX%column(r0 + col_offset) = col
                        dRdX%row(r0 + col_offset) = row
                    end do 
                end do 
                
                !set the column end pointer 
                dRdX%col_pointer(col+1) = r0 + col_offset + 1
            end if 
        end do 
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
    !!$OMP single
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
real(dp) :: vector(matrix%nrow),result(matrix%nrow)

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

!flux adjoint solve ===============
subroutine flux_adjoint_solve(mesh,options)
implicit none 

!variables - inout
type(flux_mesh) :: mesh 
type(flux_options) :: options 

!variables - local 
logical :: nanflag,resconv
integer(in32) :: iteration,cc,ee
real(dp) :: psirhores
real(dp) :: cell_timestep(4*mesh%ncell)
real(dp), dimension(:), allocatable :: dJdW,dJdX,dtotal
type(csc_matrix) :: dRdW_sp,dRdX_sp

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

!evaluate the flow jacobian 
if (options%cdisplay) then
    write(*,'(A)') '--> evaluating the flow jacobian'
end if 
call build_flow_jacobian_sparse(mesh,options,dRdW_sp)


!validate sparse jacobian
! print *, 'data ',dRdW_sp%col_pointer(1:4)
! do cc=1,65
!     print *,cc,dRdW_sp%row(cc),dRdW_sp%column(cc)
! end do 
! stop 


! do rr=1,dRdW_sp%nrow
!     do ii=dRdW_sp%row_pointer(rr),dRdW_sp%row_pointer(rr+1)

!         ! print *, 'row = ',rr,' col = ',dRdW_sp%column(ii)
!         dRdW(rr,dRdW_sp%column(ii)) = 0.0d0
!         ! dRdW(rr,dRdW_sp%column(ii)) = dRdW(rr,dRdW_sp%column(ii)) - dRdW_sp%value(ii)
!         ! dRdW(rr,:) = 0.0d0 
!         ! print *, dRdW(rr,dRdW_sp%column(ii))
!     end do
! end do 

! do ii=1,dRdW_sp%nnz

!     ! print *,dRdW_sp%row(ii),dRdW_sp%column(ii)
!     ! print *, dRdW(dRdW_sp%row(ii),dRdW_sp%column(ii)) - dRdW_sp%value(ii)
!     dRdW(dRdW_sp%row(ii),dRdW_sp%column(ii)) = 0.0d0 !dRdW(dRdW_sp%row(ii),dRdW_sp%column(ii)) - dRdW_sp%value(ii)
! end do

! print *, '========================='
! do ii=1,mesh%ncell*4

!     print *, dRdW(ii,1)

! end do 

! print *, 'shape = ',shape(dRdW)

! print *, 'dRdW error = ',sum(dRdW)


!evaluate the flow objective jacobian 
if (options%cdisplay) then
    write(*,'(A)') '--> evaluating the objective gradient'
end if 

call read_flow_field('flowfield',mesh)

call cd_objective(mesh,options,dJdW,dJdX)

print *, 'cd = ',mesh%cd




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
    psirhores = log10(sqrt(sum((mesh%residual)**2)))
    if (psirhores .LE. options%residual_convtol) then 
        resconv = .true.
    end if 

    !display 
    print *, 'iter = ',iteration,' psirho_res = ', psirhores

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
mesh%psi_rho = mesh%psi(1:mesh%ncell)
mesh%psi_u = mesh%psi(mesh%ncell+1:2*mesh%ncell)
mesh%psi_v = mesh%psi(2*mesh%ncell+1:3*mesh%ncell)
mesh%psi_e = mesh%psi(3*mesh%ncell+1:4*mesh%ncell)

!evaluate the grid jacobian 
call build_grid_jacobian_sparse(mesh,options,dRdX_sp)


! !validate sparse jacobian
! print *, 'data ',dRdX_sp%col_pointer(1:4)
! do cc=1,dRdX_sp%nnz
!     print *,cc,dRdX_sp%row(cc),dRdX_sp%column(cc),dRdX_sp%value(cc)
! end do 
! stop 

!evaluate the total derivative
allocate(dtotal(2*mesh%nvertex))
!$OMP parallel
call csc_vector_product(dRdX_sp,mesh%psi,dtotal)
!$OMP end parallel 
dtotal = dtotal + dJdX




return 
end subroutine flux_adjoint_solve



end module flux_adjoint