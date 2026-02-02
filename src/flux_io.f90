!flux 2d io module 
!max wood
!version : 0.0.4
!updated : 02-02-26

!module 
module flux_io
use io_utilities
use flux_data_methods
contains 

!read command arguments subroutine ===============
subroutine get_command_arguments(options)
implicit none

!variables - import
type(flux_options) :: options 

!variables - local 
integer(in32) :: nargs

!check and process supplied command arguments 
nargs = command_argument_count()

!set mode
if (nargs == 0) then !default to solve mode
    options%mode = 'solve'
else !set mode from first argument 
    options%mode = get_command_argument_n_str(1)
end if 
return 
end subroutine get_command_arguments

!read options ===============
subroutine read_options(options)
implicit none

!variables - inout
type(flux_options) :: options 

!variables - Local
integer(in32) :: ii,fh
integer(in64) :: itemtemp
integer(in64), dimension(:), allocatable :: dissflags

!check if file exists 
if (.NOT. file_exists('flux_options')) then 
    write(*,'(A)') '** cannot locate options file: '//trim('flux_options')
    stop 
end if 

!open file 
open(newunit=fh,file='flux_options')

!set options 
call set_log_opt(options%cdisplay,fh,'console_display')
call set_str_opt(options%meshpath,fh,'meshpath')
call set_real_opt(options%aoadeg,fh,'aoadeg')
call set_real_opt(options%machinf,fh,'machinf')
call set_real_opt(options%gamma,fh,'gamma')
call set_real_opt(options%R,fh,'R')
call set_real_opt(options%tinf,fh,'tinf')
call set_real_opt(options%rhoinf,fh,'rhoinf')
call set_real_opt(options%outflow_pratio,fh,'outflow_pratio')
itemtemp = -1 
call set_int_opt(itemtemp,fh,'niter_max')
if (itemtemp .GE. 0) then 
    options%niter_max = int(itemtemp,in32)
end if 
call set_real_opt(options%cfl,fh,'cfl')
call set_real_opt(options%k2,fh,'k2')
call set_real_opt(options%k4,fh,'k4')
call set_str_opt(options%flux_method,fh,'flux_method')
itemtemp = -1 
call set_int_opt(itemtemp,fh,'num_threads')
if (itemtemp .GE. 0) then 
    options%num_threads = int(itemtemp,in32)
end if 
call set_real_opt(options%residual_convtol,fh,'residual_convtol')
call set_int_opt_arr1(dissflags,fh,'RKstagediss_eval')
if (allocated(dissflags)) then 
    if (allocated(options%rk_dissipation)) then 
        deallocate(options%rk_dissipation)
    end if 
    allocate(options%rk_dissipation(size(dissflags)))
    do ii=1,size(dissflags)
        if (dissflags(ii) == 1) then 
            options%rk_dissipation(ii) = .true.
        else
            options%rk_dissipation(ii) = .false.
        end if 
    end do 
    deallocate(dissflags)
end if 

!close file 
close(fh)
return 
end subroutine read_options

!import mesh ===============
subroutine import_mesh(mesh,options,filename)
implicit none 

!variables - inout
character(*) :: filename
type(flux_mesh) :: mesh 
type(flux_options) :: options 

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
        if (options%mode == 'adjoint') then 
            allocate(mesh%cells_colour(mesh%ncell))
            allocate(mesh%cells_nadj(mesh%ncell))
        end if   
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
        if (options%mode == 'adjoint') then 
            allocate(mesh%vertices_colour(mesh%nvertex))
        end if 
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

!write vtu ===============
subroutine write_vtu(mesh,options,filename)
implicit none 

!variables - inout
character(*) :: filename
type(flux_mesh) :: mesh 
type(flux_options) :: options 

!variables - local
integer(in64) :: ii,jj,offset

!open vtu file
open(11,file=filename) 

!write headers
write(11,'(A)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
write(11,'(A)') ' <UnstructuredGrid>'

!write mesh data header
write(11,'(A,I0,A,I0,A)') '  <Piece NumberOfPoints="',mesh%nvertex,'" NumberOfCells="',mesh%ncell,'">'

!write vertices
write(11,'(A)') '   <Points>'
write(11,'(A)') '    <DataArray type="Float64" NumberOfComponents="3" Format="ascii">'
do ii=1,mesh%nvertex
    write(11,'(A,E17.10,A,E17.10,A,E17.10)') '     ',mesh%vertices(ii)%coordinate(1),' ',mesh%vertices(ii)%coordinate(2),' ',0.0d0
end do 
write(11,'(A)') '    </DataArray>'
write(11,'(A)') '   </Points>'

!write cells header
write(11,'(A)') '   <Cells>'

!write cells
write(11,'(A)') '    <DataArray type="Int32" Name="connectivity" Format="ascii">'
do ii=1,mesh%ncell
    do jj=1,mesh%cells(ii)%nedge
        write(11,'(A,I0)') '     ',mesh%cells(ii)%edgev1(jj) - 1
    end do 
end do 
write(11,'(A)') '    </DataArray>'

!write offsets
write(11,'(A)') '    <DataArray type="Int32" Name="offsets" Format="ascii">'
offset = 0
do ii=1,mesh%ncell
    offset = offset + mesh%cells(ii)%nedge
    write(11,'(A,I0)') '     ',offset
end do 
write(11,'(A)') '    </DataArray>'

!write cell types 
write(11,'(A)') '    <DataArray type="UInt8" Name="types" Format="ascii">'
do ii=1,mesh%ncell
    if (mesh%cells(ii)%nedge == 3) then 
        write(11,'(A)') '     5'
    elseif (mesh%cells(ii)%nedge == 4) then 
        write(11,'(A)') '     9'
    else
        write(11,'(A)') '     7'
    end if 
end do 
write(11,'(A)') '    </DataArray>'

!write cells footer
write(11,'(A)') '   </Cells>'

!write cell data  ====================
!write cell data header
write(11,'(A)') '   <CellData>'

!mach
write(11,'(A)') '    <DataArray type="Float64" Name="Mach" format="ascii" NumberOfComponents="1">'
do ii=1,mesh%ncell
    mesh%mach(ii) = sqrt(mesh%u(ii)**2 + mesh%v(ii)**2)/speed_of_sound(mesh%p(ii),mesh%rho(ii),options%gamma)
    write(11,'(A,E17.10)')'     ',mesh%mach(ii)
end do 
write(11,'(A)') '    </DataArray>'

!cp
write(11,'(A)') '    <DataArray type="Float64" Name="Cp" format="ascii" NumberOfComponents="1">'
do ii=1,mesh%ncell
    mesh%cp(ii) = pressure_coefficient(mesh%p(ii),options)
    write(11,'(A,E17.10)')'     ',mesh%cp(ii)
end do 
write(11,'(A)') '    </DataArray>'

!p
write(11,'(A)') '    <DataArray type="Float64" Name="Pressure" format="ascii" NumberOfComponents="1">'
do ii=1,mesh%ncell
    write(11,'(A,E17.10)')'     ',mesh%p(ii)
end do 
write(11,'(A)') '    </DataArray>'

!rho 
write(11,'(A)') '    <DataArray type="Float64" Name="Density" format="ascii" NumberOfComponents="1">'
do ii=1,mesh%ncell
    write(11,'(A,E17.10)')'     ',mesh%rho(ii)
end do 
write(11,'(A)') '    </DataArray>'

!velocity 
write(11,'(A)') '    <DataArray type="Float64" Name="Velocity" format="ascii" NumberOfComponents="3">'
do ii=1,mesh%ncell
    write(11,'(A,E17.10,A,E17.10,A,E17.10)')'     ',mesh%u(ii),' ',mesh%v(ii),' ',0.0d0 
end do 
write(11,'(A)') '    </DataArray>'

!rhores
if (options%mode == 'solve') then 
    write(11,'(A)') '    <DataArray type="Float64" Name="Density Resdidual" format="ascii" NumberOfComponents="1">'
    do ii=1,mesh%ncell
        write(11,'(A,E17.10)')'     ',mesh%residual(ii)
    end do 
    write(11,'(A)') '    </DataArray>'
end if 

!write adjoint variables
if (options%mode == 'adjoint') then 

    !psirho 
    write(11,'(A)') '    <DataArray type="Float64" Name="psi_rho" format="ascii" NumberOfComponents="1">'
    do ii=1,mesh%ncell
        write(11,'(A,E17.10)')'     ',mesh%psi_rho(ii)
    end do 
    write(11,'(A)') '    </DataArray>'

    !psiuv 
    write(11,'(A)') '    <DataArray type="Float64" Name="psi_uv" format="ascii" NumberOfComponents="3">'
    do ii=1,mesh%ncell
        write(11,'(A,E17.10,A,E17.10,A,E17.10)')'     ',mesh%psi_u(ii),' ',mesh%psi_v(ii),' ',0.0d0 
    end do 
    write(11,'(A)') '    </DataArray>'

    !psie 
    write(11,'(A)') '    <DataArray type="Float64" Name="psi_e" format="ascii" NumberOfComponents="1">'
    do ii=1,mesh%ncell
        write(11,'(A,E17.10)')'     ',mesh%psi_e(ii)
    end do 
    write(11,'(A)') '    </DataArray>'
end if 

!write cell data footer
write(11,'(A)') '   </CellData>'
!write cell data  ====================

!write point data  ====================
!write point data header
write(11,'(A)') '   <PointData>'

!write adjoint variables
if (options%mode == 'adjoint') then 

    !vertex derivative 
    write(11,'(A)') '    <DataArray type="Float64" Name="vertex_derivative" format="ascii" NumberOfComponents="3">'
    do ii=1,mesh%nvertex
        write(11,'(A,E17.10,A,E17.10,A,E17.10)')'     ',mesh%vertex_derivative_x(ii),' ',mesh%vertex_derivative_y(ii),' ',0.0d0 
    end do 
    write(11,'(A)') '    </DataArray>'
end if 

!write point data footer
write(11,'(A)') '   </PointData>'
!write point data  ====================

!write footers
write(11,'(A)') '  </Piece>'
write(11,'(A)') ' </UnstructuredGrid>'
write(11,'(A)') '</VTKFile>'

!close vtu file 
close(11)
return 
end subroutine write_vtu


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

write(11,'(A,I0)') 'CELL_DATA ',mesh%ncell

write(11,'(A)') 'SCALARS cp double' !cp
write(11,'(A)') 'LOOKUP_TABLE default'
do ii=1,mesh%ncell
    mesh%cp(ii) = pressure_coefficient(mesh%p(ii),options)
    write(11,'(E17.10)') mesh%cp(ii)
end do 

write(11,'(A)') 'SCALARS p double' !p
write(11,'(A)') 'LOOKUP_TABLE default'
do ii=1,mesh%ncell
    write(11,'(E17.10)') mesh%p(ii)
end do 

write(11,'(A)') 'SCALARS rho double' !rho
write(11,'(A)') 'LOOKUP_TABLE default'
do ii=1,mesh%ncell
    write(11,'(E17.10)') mesh%rho(ii)
end do 

write(11,'(A)') 'SCALARS mach double' !mach
write(11,'(A)') 'LOOKUP_TABLE default'
do ii=1,mesh%ncell
    mesh%mach(ii) = sqrt(mesh%u(ii)**2 + mesh%v(ii)**2)/speed_of_sound(mesh%p(ii),mesh%rho(ii),options%gamma)
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

if (options%mode == 'solve') then 
    write(11,'(A)') 'SCALARS density_residual double' !rhores
    write(11,'(A)') 'LOOKUP_TABLE default'
    do ii=1,mesh%ncell
        write(11,'(E17.10)') mesh%residual(ii)
    end do 
end if 

!write adjoint variables
if (options%mode == 'adjoint') then 

    write(11,'(A)') 'SCALARS cell_colour integer' !cell colour
    write(11,'(A)') 'LOOKUP_TABLE default'
    do ii=1,mesh%ncell
        write(11,'(I0)') mesh%cells_colour(ii)
    end do 

    write(11,'(A)') 'SCALARS psi_rho double' !psi_rho
    write(11,'(A)') 'LOOKUP_TABLE default'
    do ii=1,mesh%ncell
        write(11,'(E17.10)') mesh%psi_rho(ii)
    end do 

    write(11,'(A)') 'SCALARS psi_u double' !psi_u
    write(11,'(A)') 'LOOKUP_TABLE default'
    do ii=1,mesh%ncell
        write(11,'(E17.10)') mesh%psi_u(ii)
    end do 

    write(11,'(A)') 'SCALARS psi_v double' !psi_v
    write(11,'(A)') 'LOOKUP_TABLE default'
    do ii=1,mesh%ncell
        write(11,'(E17.10)') mesh%psi_v(ii)
    end do 

    write(11,'(A)') 'SCALARS psi_e double' !psi_e
    write(11,'(A)') 'LOOKUP_TABLE default'
    do ii=1,mesh%ncell
        write(11,'(E17.10)') mesh%psi_e(ii)
    end do 

    write(11,'(A,I0)') 'POINT_DATA ',mesh%nvertex

    write(11,'(A)') 'SCALARS vertex_colour integer' !vertex colour
    write(11,'(A)') 'LOOKUP_TABLE default'
    do ii=1,mesh%nvertex
        write(11,'(I0)') mesh%vertices_colour(ii)
    end do 

    write(11,'(A)') 'SCALARS vertex_derivative_x double' !vertex derivative x
    write(11,'(A)') 'LOOKUP_TABLE default'
    do ii=1,mesh%nvertex
        write(11,'(E17.10)') mesh%vertex_derivative_x(ii)
    end do 

    write(11,'(A)') 'SCALARS vertex_derivative_y double' !vertex derivative y
    write(11,'(A)') 'LOOKUP_TABLE default'
    do ii=1,mesh%nvertex
        write(11,'(E17.10)') mesh%vertex_derivative_y(ii)
    end do 

    write(11,'(A)') 'VECTORS vertex_derivative double' !vertex derivative 
    do ii=1,mesh%nvertex
        write(11,'(E17.10,A,E17.10,A,E17.10)') mesh%vertex_derivative_x(ii),' ',mesh%vertex_derivative_y(ii),' ',0.0d0
    end do 
end if 

!close vtk file 
close(11)
return 
end subroutine write_vtk

!write flow field =========================
subroutine write_flow_field(filename,mesh)
implicit none 

!variables - inout
character(*), intent(in) :: filename
type(flux_mesh) :: mesh 

!variables local
integer(in32) :: cc 

!write flow data 
open(11,file=filename) 
do cc=1,mesh%ncell
    write(11,'(E17.10,A,E17.10,A,E17.10,A,E17.10,A,E17.10)') mesh%rho(cc),' ',mesh%u(cc),' ',mesh%v(cc),' ',mesh%p(cc),' ',mesh%e(cc)
end do 
close(11)
return 
end subroutine write_flow_field

!read flow field =========================
subroutine read_flow_field(filename,mesh)
implicit none 

!variables - inout
character(*), intent(in) :: filename
type(flux_mesh) :: mesh 

!variables local
integer(in32) :: cc 

!write flow data 
open(11,file=filename) 
do cc=1,mesh%ncell
    read(11,*) mesh%rho(cc),mesh%u(cc),mesh%v(cc),mesh%p(cc),mesh%e(cc)
end do 
close(11)
return 
end subroutine read_flow_field

!write gradient =========================
subroutine write_gradient(filename,mesh)
implicit none 

!variables - inout
character(*), intent(in) :: filename
type(flux_mesh) :: mesh 

!variables local
integer(in32) :: vv  

!write
open(11,file=filename) 
do vv=1,mesh%nvertex
    write(11,'(E17.10,A,E17.10)') mesh%vertex_derivative_x(vv),' ',mesh%vertex_derivative_y(vv)
end do 
close(11)
return 
end subroutine write_gradient

!write forces =========================
subroutine write_forces(filename,mesh)
implicit none 

!variables - inout
character(*), intent(in) :: filename
type(flux_mesh) :: mesh 

!write
open(11,file=filename) 
write(11,'(A,E17.10)') 'cl = ',mesh%cl
write(11,'(A,E17.10)') 'cd = ',mesh%cd
close(11)
return 
end subroutine write_forces

end module flux_io
