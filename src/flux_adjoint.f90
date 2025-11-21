!flux 2d adjoint module 
!max wood
!version : 0.0.1
!updated : 21-11-25


module flux_adjoint
use edge_flux
use io_utilities
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
subroutine build_flow_jacobian_full(mesh,options)
implicit none

!variables - inout
type(flux_mesh) :: mesh 
type(flux_options) :: options 

real(dp), dimension(:,:), allocatable :: dRdW

!variables - local 
integer(in32) :: ii,cc,vv
integer(in32) :: cidx,ridx

!allocate the flow jacobian 
allocate(dRdW(4*mesh%ncell,4*mesh%ncell))
dRdW = 0.0d0 


!construct the jacobian 
do cc=1,mesh%ncell !perturb each cell
    write(*,'(A,I0)') 'cell: ',cc 
    do vv=1,4 !perturb each conservative variable in this cell

        !perturb the conservative variable w_vv in cell cc

        !evaluate the flow residual
        
        !-> results in dR_vv by dw_vv in cell cc for all other cells
        !   this is a full column of dRdW
        !   column index = (cc - 1)*4 + vv
        
        !set the column index
        cidx = (cc - 1)*4 + vv

        !unpack residuals 
        ridx = 1 
        do ii=1,mesh%ncell
            dRdW(ridx,cc) = 1.0!mesh%r1(cc)
            dRdW(ridx+1,cc) = 1.0!mesh%r2(cc)
            dRdW(ridx+2,cc) = 1.0!mesh%r3(cc)
            dRdW(ridx+3,cc) = 1.0!mesh%r4(cc)
            ridx = ridx + 4
        end do 


    end do
end do 


return
end subroutine build_flow_jacobian_full

!sparse flow jacobian ===============





end module flux_adjoint