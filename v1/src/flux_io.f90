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
integer(in64) :: ii,jj
integer(in64) :: iostatus,cindex,nedge
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
        do ii=1,mesh%ncell
            read(11,*) cindex,nedge 
            mesh%cells(cindex)%nedge = nedge
            mesh%cells(cindex)%index = cindex
            allocate(mesh%cells(cindex)%edges(nedge,3))
            allocate(mesh%cells(cindex)%edges_dx(nedge))
            allocate(mesh%cells(cindex)%edges_dy(nedge))
            allocate(mesh%cells(cindex)%edges_nx(nedge))
            allocate(mesh%cells(cindex)%edges_ny(nedge))
            allocate(mesh%cells(cindex)%edges_len(nedge))
            allocate(mesh%cells(cindex)%edges_specrad(nedge))
            do jj=1,nedge
                read(11,*) mesh%cells(cindex)%edges(jj,:)
            end do 
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
        allocate(mesh%vertices(mesh%nvertex,2))
        do ii=1,mesh%nvertex
            read(11,*) mesh%vertices(ii,:)
        end do 
        exit 
    end if 
end do 

!close mesh file 
close(11)

!evaluate cell volumes 
call mesh%get_cell_volumes()

!evaluate cell edge geometries 
call mesh%get_cell_edge_geometries()
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
    write(11,'(E17.10,A,E17.10,A,E17.10)') mesh%vertices(ii,1),' ',mesh%vertices(ii,2),' ',0.0d0
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
        write(11,'(A,I0)',advance='no') ' ',mesh%cells(ii)%edges(jj,1) - 1
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
    write(11,'(E17.10)') pressure_coefficient(mesh%cells(ii)%p,options)
end do 

write(11,'(A)') 'SCALARS p double' !p
write(11,'(A)') 'LOOKUP_TABLE default'
do ii=1,mesh%ncell
    write(11,'(E17.10)') mesh%cells(ii)%p
end do 

write(11,'(A)') 'SCALARS Mach double' !mach
write(11,'(A)') 'LOOKUP_TABLE default'
do ii=1,mesh%ncell
    write(11,'(E17.10)') mesh%cells(ii)%mach
end do 

write(11,'(A)') 'SCALARS u double' !u
write(11,'(A)') 'LOOKUP_TABLE default'
do ii=1,mesh%ncell
    write(11,'(E17.10)') mesh%cells(ii)%u
end do 

write(11,'(A)') 'SCALARS v double' !v
write(11,'(A)') 'LOOKUP_TABLE default'
do ii=1,mesh%ncell
    write(11,'(E17.10)') mesh%cells(ii)%v
end do 

write(11,'(A)') 'SCALARS density_residual double' !rhores
write(11,'(A)') 'LOOKUP_TABLE default'
do ii=1,mesh%ncell
    write(11,'(E17.10)') mesh%cells(ii)%r1 + mesh%cells(ii)%d1
end do 

!close vtk file 
close(11)
return 
end subroutine write_vtk


!read restart file =========================
subroutine read_restart_file(filename,mesh,options)
implicit none 

!variables - inout
character(*), intent(in) :: filename
type(flux_options) :: options 
type(flux_mesh) :: mesh 

!variables - Local 
integer(in64) :: cc 

!read file 
open(11,file=filename)
do cc=1,mesh%ncell
    read(11,*) mesh%cells(cc)%u,mesh%cells(cc)%v,&
    mesh%cells(cc)%mach,mesh%cells(cc)%p,mesh%cells(cc)%rho,mesh%cells(cc)%cp,mesh%cells(cc)%e
end do 
close(11)

!set the conservative variables in each cell 
do cc=1,mesh%ncell
    call mesh%cells(cc)%prim2con(options%gamma)
end do 
return 
end subroutine read_restart_file


end module flux_io
