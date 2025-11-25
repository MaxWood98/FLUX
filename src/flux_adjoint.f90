!flux 2d adjoint module 
!max wood
!version : 0.0.1
!updated : 21-11-25


module flux_adjoint


use flux_io

use flux_solve

use flux_data_methods
contains 

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
    write(*,'(A,I0,A,I0,A)') '    {colour min/max = ',minval(mesh%cells_colour),' / ',maxval(mesh%cells_colour),'}' 
end if 
return 
end subroutine colour_cells

!full flow jacobian ===============
subroutine build_flow_jacobian_full(mesh,options,dRdW)
implicit none

!variables - inout
type(flux_mesh) :: mesh 
type(flux_options) :: options 

real(dp), dimension(:,:), allocatable :: dRdW

!variables - local 
integer(in32) :: ii,cc,vv,ee,cc2,aa,rr
integer(in32) :: cidx,ridx,cadj,row,col

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
    write(*,'(A,I0,A,I0)') 'cell: ',cc,' / ',mesh%ncell
    do vv=1,4 !perturb each conservative variable in this cell

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



        !-> results in dR_vv by dw_vv in cell cc for all other cells
        !   this is a full column of dRdW
        !   column index = (cc - 1)*4 + vv
        
        !set the column index
        cidx = (vv - 1)*mesh%ncell + cc

        !unpack residuals 
        dRdW(1:mesh%ncell,cidx) = (pr1(:) - r10(:))/h

        dRdW(mesh%ncell+1:2*mesh%ncell,cidx) = (pr2(:) - r20(:))/h

        dRdW(2*mesh%ncell+1:3*mesh%ncell,cidx) = (pr3(:) - r30(:))/h

        dRdW(3*mesh%ncell+1:4*mesh%ncell,cidx) = (pr4(:) - r40(:))/h

        ! ridx = 1 
        ! do ii=1,mesh%ncell
        !     dRdW(ridx,cidx) = abs(pr1(ii) - r10(ii))/h
        !     dRdW(ridx+1,cidx) = abs(pr2(ii) - r20(ii))/h
        !     dRdW(ridx+2,cidx) = abs(pr3(ii) - r30(ii))/h
        !     dRdW(ridx+3,cidx) = abs(pr4(ii) - r40(ii))/h
        !     ridx = ridx + 4
        ! end do 

    end do
end do 

! ridx = 1
! do rr=1,4
!     do vv=1,4
!         do cc=1,mesh%ncell

!             row = cc + (rr - 1)*mesh%ncell
!             col = cc + (vv - 1)*mesh%ncell
!             dRdW(row,col) = 1

!             print *, dRdW(row,col) 

!             do aa=1,mesh%cells(cc)%nedge

!                 !get the adjacent cell 
!                 if (mesh%edges(mesh%cells(cc)%edge(aa))%c1 .NE. cc) then 
!                     cadj = mesh%edges(mesh%cells(cc)%edge(aa))%c1 
!                 else
!                     cadj = mesh%edges(mesh%cells(cc)%edge(aa))%c2 
!                 end if 
!                 if (cadj .LT. 0) then !skip if boundary condition 
!                     cycle 
!                 end if 

!                 row = cadj + (rr - 1)*mesh%ncell
!                 col = cc + (vv - 1)*mesh%ncell
!                 dRdW(row,col) = 1
                
!                 !  print *, dRdW(row,col) 


!             end do
!         end do 
!     end do 
! end do 

! print *, dRdW

return
end subroutine build_flow_jacobian_full

!sparse flow jacobian ===============
subroutine build_flow_jacobian_sparse(mesh,options,dRdW)
implicit none

!variables - inout
type(flux_mesh) :: mesh 
type(flux_options) :: options 
type(csr_matrix) :: dRdW

!variables - local 
integer(in32) :: clr,cc,vv,ee,rr,aa
integer(in32) :: ncolour,b0,r0,row_offset,row,col,cadj,nblock
integer(in32) :: row_offset_index(4*mesh%ncell)

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
allocate(dRdW%row_pointer(4*mesh%ncell + 1))


b0 = 0 

!evaluate the jacobian 
nblock = (sum(mesh%cells_nadj) + mesh%ncell)*4
row_offset_index(:) = 0
ncolour = maxval(mesh%cells_colour)
do clr=1,ncolour
    write(*,'(A,I0,A,I0)') '    colour: ',clr,' / ',ncolour
    do vv=1,4 !perturb each conservative variable in all cells of this colour

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

        !extract non zero values and populate the flow jacobian 
        do cc=1,mesh%ncell
            if (mesh%cells_colour(cc) == clr) then 


                !extract each residual for this cell
                do rr=1,4

                    !get the start index of this row in the sparse jacobian
                    r0 = (sum(mesh%cells_nadj(1:cc-1)) + cc - 1)*4 + nblock*(rr - 1) + 1

                    !get the row and column index of this entry 
                    row = cc + (rr - 1)*mesh%ncell
                    col = cc + (vv - 1)*mesh%ncell

                    

                    !get the row_offset of the entry for this cell in this row 
                    row_offset = row_offset_index(row) 
                    row_offset_index(row) = row_offset_index(row) + 1

                    ! print *, r0 + row_offset,48480

                    !insert the values
                    if (rr == 1) then 
                        dRdW%value(r0 + row_offset) = (pr1(cc) - r10(cc))/h
                    elseif (rr == 2) then 
                        dRdW%value(r0 + row_offset) = (pr2(cc) - r20(cc))/h
                    elseif (rr == 3) then 
                        dRdW%value(r0 + row_offset) = (pr3(cc) - r30(cc))/h
                    elseif (rr == 4) then 
                        dRdW%value(r0 + row_offset) = (pr4(cc) - r40(cc))/h
                    end if 
                    dRdW%column(r0 + row_offset) = col
                    dRdW%row(r0 + row_offset) = row
                    b0 = b0 + 1

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

                        !get the start index of this row in the sparse jacobian
                        r0 = (sum(mesh%cells_nadj(1:cadj-1)) + cadj - 1)*4 + nblock*(rr - 1) + 1

                        !get the row and column index of this entry 
                        row = cadj + (rr - 1)*mesh%ncell
                        col = cc + (vv - 1)*mesh%ncell

                        !get the row_offset of the entry for this cell in this row 
                        row_offset = row_offset_index(row) 
                        row_offset_index(row) = row_offset_index(row) + 1

                        ! print *, r0,row_offset

                        !insert the values
                        if (rr == 1) then 
                            dRdW%value(r0 + row_offset) = (pr1(cadj) - r10(cadj))/h
                        elseif (rr == 2) then 
                            dRdW%value(r0 + row_offset) = (pr2(cadj) - r20(cadj))/h
                        elseif (rr == 3) then 
                            dRdW%value(r0 + row_offset) = (pr3(cadj) - r30(cadj))/h
                        elseif (rr == 4) then 
                            dRdW%value(r0 + row_offset) = (pr4(cadj) - r40(cadj))/h
                        end if 
                        dRdW%column(r0 + row_offset) = col
                        dRdW%row(r0 + row_offset) = row
                        b0 = b0 + 1

                    end do 
                end do 
            end if 
        end do 

        ! !extract non zero values and populate the flow jacobian 
        ! do cc=1,mesh%ncell
        !     if (mesh%cells_colour(cc) == clr) then 

        !         !get the block start index in the sparse jacobian 
        !         b0 = (sum(mesh%cells_nadj(1:cc-1)) + cc - 1)*16 
                
        !         !extract each residual for this cell
        !         do rr=1,4
                
        !             !get the index of the start of the row for this residual in the sparse jacobian 
        !             r0 = b0 + (mesh%cells_nadj(cc) + 1)*4*(rr - 1)
                    
        !             !get the row and column index of this entry 
        !             row = (cc - 1)*4 + rr
        !             col = (cc - 1)*4 + vv 

        !             !get the row_offset of the entry for this cell in this row 
        !             row_offset = row_offset_index(row) 
        !             row_offset_index(row) = row_offset_index(row) + 1
                    
        !             !insert the values
        !             !goes in position r0 + row_offset
        !             if (rr == 1) then 
        !                 dRdW%value(r0 + row_offset) = abs(pr1(cc) - r10(cc))/h
        !             elseif (rr == 2) then 
        !                 dRdW%value(r0 + row_offset) = abs(pr2(cc) - r20(cc))/h
        !             elseif (rr == 3) then 
        !                 dRdW%value(r0 + row_offset) = abs(pr3(cc) - r30(cc))/h
        !             elseif (rr == 4) then 
        !                 dRdW%value(r0 + row_offset) = abs(pr4(cc) - r40(cc))/h
        !             end if 
        !             dRdW%column(r0 + row_offset) = col
        !         end do 
          

        !         !extract each residual for this cells adjacent cells
        !         do aa=1,mesh%cells(cc)%nedge

        !             !get the adjacent cell 
        !             if (mesh%edges(mesh%cells(cc)%edge(aa))%c1 .NE. cc) then 
        !                 cadj = mesh%edges(mesh%cells(cc)%edge(aa))%c1 
        !             else
        !                 cadj = mesh%edges(mesh%cells(cc)%edge(aa))%c2 
        !             end if 
        !             if (cadj .LT. 0) then !skip if boundary condition 
        !                 cycle 
        !             end if 



        !             !get the block start index in the sparse jacobian 
        !             b0 = (sum(mesh%cells_nadj(1:cadj-1)) + cadj - 1)*16 
                    
        !             !extract each residual for this cell
        !             do rr=1,4
                    
        !                 !get the index of the start of the row for this residual in the sparse jacobian 
        !                 r0 = b0 + (mesh%cells_nadj(cadj) + 1)*4*(rr - 1)
                        
        !                 !get the row and column index of this entry 
        !                 row = (cadj - 1)*4 + rr
        !                 col = (cc - 1)*4 + vv

        !                 !get the row_offset of the entry for this cell in this row 
        !                 row_offset = row_offset_index(row) 
        !                 row_offset_index(row) = row_offset_index(row) + 1
                        
        !                 !insert the values
        !                 !goes in position r0 + row_offset
        !                 if (rr == 1) then 
        !                     dRdW%value(r0 + row_offset) = abs(pr1(cadj) - r10(cadj))/h
        !                 elseif (rr == 2) then 
        !                     dRdW%value(r0 + row_offset) = abs(pr2(cadj) - r20(cadj))/h
        !                 elseif (rr == 3) then 
        !                     dRdW%value(r0 + row_offset) = abs(pr3(cadj) - r30(cadj))/h
        !                 elseif (rr == 4) then 
        !                     dRdW%value(r0 + row_offset) = abs(pr4(cadj) - r40(cadj))/h
        !                 end if 
        !                 dRdW%column(r0 + row_offset) = col
        !             end do 
                    



        !         end do 



        !     end if
        ! end do 


    end do 
end do 


print *,'visited = ',b0





! !set the row end pointers
! row = 2
! dRdW%row_pointer(1) = 1
! do rr=1,4
!     do cc=1,mesh%ncell
!         ! r0 = (sum(mesh%cells_nadj(1:cc-1)) + cc - 1)*16*(rr - 1)  (sum(mesh%cells_nadj) + mesh%ncell)*16

!         print *, 'rindex = ',(sum(mesh%cells_nadj(1:cc-1)) + cc - 1)*4 + nblock*(rr - 1)      + (mesh%cells_nadj(cc) + 1)*4

!         ! print *, 'rindex = ',(sum(mesh%cells_nadj(1:cc-1)) + cc - 1)*16*(rr - 1) 

!     end do 
! end do 

print *, 'here'


! !set the row end pointers
! row = 2
! dRdW%row_pointer(1) = 1
! do cc=1,mesh%ncell

!     b0 = (sum(mesh%cells_nadj(1:cc-1)) + cc - 1)*16 

!     do rr=1,4

!         r0 = b0 + (mesh%cells_nadj(cc) + 1)*4*(rr - 1)

!         dRdW%row_pointer(row) = r0 + (mesh%cells_nadj(cc) + 1)*4
!         row = row + 1
!     end do 
! end do 

! do cc=1,4*mesh%ncell+1
!     print *, dRdW%row_pointer(cc)
! end do


return 
end subroutine build_flow_jacobian_sparse



!flux adjoint solve ===============
subroutine flux_adjoint_solve(mesh,options)
implicit none 

!variables - inout
type(flux_mesh) :: mesh 
type(flux_options) :: options 

!variables - local 
integer(in32) :: rr,ii 
real(dp), dimension(:,:), allocatable :: dRdW
type(csr_matrix) :: dRdW_sp






!get the number of adjacent cells for each mesh cell
call mesh%get_cells_nadj()

!colour the cells
if (options%cdisplay) then
    write(*,'(A)') '--> colouring mesh'
end if 
call colour_cells(mesh,options)

!evaluate the flow jacobian 
if (options%cdisplay) then
    write(*,'(A)') '--> evaluating the flow jacobian'
end if 

call read_flow_field('flowfield',mesh)
call build_flow_jacobian_full(mesh,options,dRdW)

call read_flow_field('flowfield',mesh)
call build_flow_jacobian_sparse(mesh,options,dRdW_sp)


!validate sparse jacobian






! do rr=1,dRdW_sp%nrow
!     do ii=dRdW_sp%row_pointer(rr),dRdW_sp%row_pointer(rr+1)

!         ! print *, 'row = ',rr,' col = ',dRdW_sp%column(ii)
!         dRdW(rr,dRdW_sp%column(ii)) = 0.0d0
!         ! dRdW(rr,dRdW_sp%column(ii)) = dRdW(rr,dRdW_sp%column(ii)) - dRdW_sp%value(ii)
!         ! dRdW(rr,:) = 0.0d0 
!         ! print *, dRdW(rr,dRdW_sp%column(ii))
!     end do
! end do 

do ii=1,dRdW_sp%nnz
    ! print *,dRdW_sp%row(ii),dRdW_sp%column(ii)
    ! print *, dRdW(dRdW_sp%row(ii),dRdW_sp%column(ii)) - dRdW_sp%value(ii)
    dRdW(dRdW_sp%row(ii),dRdW_sp%column(ii)) = 0.0d0 !dRdW(dRdW_sp%row(ii),dRdW_sp%column(ii)) - dRdW_sp%value(ii)
end do

do ii=1,mesh%ncell*4

    print *, dRdW(2,ii)

end do 

print *, 'shape = ',shape(dRdW)

print *, 'dRdW error = ',sum(dRdW)




return 
end subroutine flux_adjoint_solve



end module flux_adjoint