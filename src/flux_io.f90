!flux 2d io module 
!max wood
!version : 0.0.1
!updated : 28-03-25

!module 
module flux_io
use io_utilities
use flux_data_methods
contains 

!import mesh ===============
subroutine import_mesh(mesh,filename)
implicit none 

!variables - inout
character(*) :: filename
type(flux_mesh) :: mesh 

!variables - local
integer(in32) :: ii,jj
integer(in32) :: iostatus,cindex,nedge
character(len=100) :: rtemp 

!check if file exists 
if (.NOT. file_exists(filename)) then 
    write(*,'(A)') '** cannot locate mesh file: '//trim(filename)
    stop 
end if 

!open mesh file
open(11,file=filename) 

!read cells
iostatus = 0 
do while (iostatus == 0)
    read(11,'(A)',iostat=iostatus) rtemp
    if (rtemp(1:5) == 'ncell') then 
        read(rtemp(9:len_trim(rtemp)),*) mesh%ncell
        allocate(mesh%cells(mesh%ncell))
        allocate(mesh%cells_specrad(mesh%ncell))
        allocate(mesh%cells_volume(mesh%ncell))
        allocate(mesh%rho(mesh%ncell))
        allocate(mesh%u(mesh%ncell))
        allocate(mesh%v(mesh%ncell))
        allocate(mesh%p(mesh%ncell))
        allocate(mesh%mach(mesh%ncell))
        allocate(mesh%e(mesh%ncell))
        allocate(mesh%cp(mesh%ncell))
        allocate(mesh%w1(mesh%ncell))
        allocate(mesh%w2(mesh%ncell))
        allocate(mesh%w3(mesh%ncell))
        allocate(mesh%w4(mesh%ncell))
        allocate(mesh%w10(mesh%ncell))
        allocate(mesh%w20(mesh%ncell))
        allocate(mesh%w30(mesh%ncell))
        allocate(mesh%w40(mesh%ncell))
        allocate(mesh%cells_dt(mesh%ncell))
        allocate(mesh%l1(mesh%ncell))
        allocate(mesh%l2(mesh%ncell))
        allocate(mesh%l3(mesh%ncell))
        allocate(mesh%l4(mesh%ncell))
        allocate(mesh%cells_psensor(mesh%ncell))
        allocate(mesh%residual(mesh%ncell))
        

        do ii=1,mesh%ncell
            read(11,*) cindex,nedge 
            mesh%cells(cindex)%nedge = nedge
            mesh%cells(cindex)%index = cindex
            allocate(mesh%cells(cindex)%edgev1(nedge))
            allocate(mesh%cells(cindex)%edgev2(nedge))
            allocate(mesh%cells(cindex)%edgec(nedge))
            allocate(mesh%cells(cindex)%edge(nedge))    
            allocate(mesh%cells(cindex)%edge_sign(nedge))




            do jj=1,nedge
                read(11,*) mesh%cells(cindex)%edgev1(jj),mesh%cells(cindex)%edgev2(jj),mesh%cells(cindex)%edgec(jj),mesh%cells(cindex)%edge(jj)
            end do 
        end do 
        exit 
    end if 
end do
rewind(11)

!read edges 
iostatus = 0 
do while (iostatus == 0)
    read(11,'(A)',iostat=iostatus) rtemp
    if (rtemp(1:5) == 'nedge') then 
        read(rtemp(8:len_trim(rtemp)),*) mesh%nedge
        allocate(mesh%edges(mesh%nedge))
        allocate(mesh%edges_specrad(mesh%nedge))
        allocate(mesh%edges_r1(mesh%nedge))
        allocate(mesh%edges_r2(mesh%nedge))
        allocate(mesh%edges_r3(mesh%nedge))
        allocate(mesh%edges_r4(mesh%nedge))
        allocate(mesh%edges_l1(mesh%nedge))
        allocate(mesh%edges_l2(mesh%nedge))
        allocate(mesh%edges_l3(mesh%nedge))
        allocate(mesh%edges_l4(mesh%nedge))
        allocate(mesh%edges_d1(mesh%nedge))
        allocate(mesh%edges_d2(mesh%nedge))
        allocate(mesh%edges_d3(mesh%nedge))
        allocate(mesh%edges_d4(mesh%nedge))
        allocate(mesh%edges_pn(mesh%nedge))
        allocate(mesh%edges_pd(mesh%nedge))


        do ii=1,mesh%nedge
            mesh%edges(ii)%index = ii 
            read(11,*) mesh%edges(ii)%v1,mesh%edges(ii)%v2,mesh%edges(ii)%c1,mesh%edges(ii)%c2
        end do 
        exit
    end if 
end do 
rewind(11)

!read vertices 
iostatus = 0 
do while (iostatus == 0)
    read(11,'(A)',iostat=iostatus) rtemp
    if (rtemp(1:7) == 'nvertex') then 
        read(rtemp(11:len_trim(rtemp)),*) mesh%nvertex
        allocate(mesh%vertices(mesh%nvertex))
        do ii=1,mesh%nvertex
            mesh%vertices(ii)%index = ii 
            read(11,*) mesh%vertices(ii)%coordinate
        end do 
        exit 
    end if 
end do 

!close mesh file 
close(11)

!set the sign for each edge corresponding to each cell 
do ii=1,mesh%ncell
    do jj=1,mesh%cells(ii)%nedge
        if (ii == mesh%edges(mesh%cells(ii)%edge(jj))%c1) then 
            mesh%cells(ii)%edge_sign(jj) = 1.0
        elseif (ii == mesh%edges(mesh%cells(ii)%edge(jj))%c2) then 
            mesh%cells(ii)%edge_sign(jj) = -1.0
        else
            write(*,'(A)') '** invalid edge cell link '
        end if 
    end do 
end do 

!evaluate the edge geometries
call mesh%get_edges_geometry()

!evaluate the cell volumes 
call mesh%get_cells_volume()
return 
end subroutine import_mesh

!write vtk ===============
subroutine write_vtk(mesh,options,filename)
implicit none 

!variables - inout
character(*) :: filename
type(flux_mesh) :: mesh 
type(flux_options) :: options 

!variables - local
integer(in64) :: ii,jj
integer(in64) :: celldlen

!open vtk file
open(11,file=filename) 

!write header
write(11,'(A)') '# vtk DataFile Version 2.0'
write(11,'(A)') 'flowfield data'
write(11,'(A)') 'ASCII'
write(11,'(A)') 'DATASET UNSTRUCTURED_GRID'
write(11,'(A)') ''

!write vertices
write(11,'(A,I0,A)') 'POINTS ',mesh%nvertex,' double'
do ii=1,mesh%nvertex
    write(11,'(E17.10,A,E17.10,A,E17.10)') mesh%vertices(ii)%coordinate(1),' ',mesh%vertices(ii)%coordinate(2),' ',0.0d0
end do 
write(11,'(A)') ''

!write cells 
celldlen = 0
do ii=1,mesh%ncell
    celldlen = celldlen + mesh%cells(ii)%nedge + 1
end do 
write(11,'(A,I0,A,I0)') 'CELLS ',mesh%ncell,' ',celldlen
do ii=1,mesh%ncell
    write(11,'(I0)',advance='no') mesh%cells(ii)%nedge
    do jj=1,mesh%cells(ii)%nedge
        write(11,'(A,I0)',advance='no') ' ',mesh%cells(ii)%edgev1(jj) - 1
    end do 
    write(11,'(A)') '' 
end do 
write(11,'(A)') ''

!write cell types (7 polygon)
write(11,'(A,I0)') 'CELL_TYPES ',mesh%ncell
do ii=1,mesh%ncell
    write(11,'(I0)') 7 
end do 

!write cell based data
write(11,'(A,I0)') 'CELL_DATA ',mesh%ncell
write(11,'(A)') 'SCALARS Cp double' !cp
write(11,'(A)') 'LOOKUP_TABLE default'
do ii=1,mesh%ncell
    write(11,'(E17.10)') pressure_coefficient(mesh%p(ii),options)
end do 

write(11,'(A)') 'SCALARS p double' !p
write(11,'(A)') 'LOOKUP_TABLE default'
do ii=1,mesh%ncell
    write(11,'(E17.10)') mesh%p(ii)
end do 

write(11,'(A)') 'SCALARS Mach double' !mach
write(11,'(A)') 'LOOKUP_TABLE default'
do ii=1,mesh%ncell
    write(11,'(E17.10)') mesh%mach(ii)
end do 

write(11,'(A)') 'SCALARS u double' !u
write(11,'(A)') 'LOOKUP_TABLE default'
do ii=1,mesh%ncell
    write(11,'(E17.10)') mesh%u(ii)
end do 

write(11,'(A)') 'SCALARS v double' !v
write(11,'(A)') 'LOOKUP_TABLE default'
do ii=1,mesh%ncell
    write(11,'(E17.10)') mesh%v(ii)
end do 

! write(11,'(A)') 'SCALARS density_residual double' !rhores
! write(11,'(A)') 'LOOKUP_TABLE default'
! do ii=1,mesh%ncell
!     write(11,'(E17.10)') mesh%cells(ii)%r1 + mesh%cells(ii)%d1
! end do 

!close vtk file 
close(11)
return 
end subroutine write_vtk


! !read restart file =========================
! subroutine read_restart_file(filename,mesh,options)
! implicit none 

! !variables - inout
! character(*), intent(in) :: filename
! type(flux_options) :: options 
! type(flux_mesh) :: mesh 

! !variables - Local 
! integer(in64) :: cc 

! !read file 
! open(11,file=filename)
! do cc=1,mesh%ncell
!     read(11,*) mesh%cells(cc)%u,mesh%cells(cc)%v,&
!     mesh%cells(cc)%mach,mesh%cells(cc)%p,mesh%cells(cc)%rho,mesh%cells(cc)%cp,mesh%cells(cc)%e
! end do 
! close(11)

! !set the conservative variables in each cell 
! do cc=1,mesh%ncell
!     call mesh%cells(cc)%prim2con(options%gamma)
! end do 
! return 
! end subroutine read_restart_file


end module flux_io
