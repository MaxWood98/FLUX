!flux 2d - a cell based euler solver
!max wood
!version : 0.0.7
!updated : 05-03-25

!program 
program flux2d
use flux_io
use flux_solve
use flux_adjoint
implicit none 

!variables 
type(flux_mesh) :: mesh 
type(flux_options) :: options 

!set default options 
call set_default_options(options)

!read command arguments 
call get_command_arguments(options)

!read options 
call read_options(options)

!display splash 
if (options%cdisplay) then
    write(*,'(A)') ' '
    write(*,'(A)')'+--------------------------------------------+'
    write(*,'(A)')'|                   flux 2d                  |'
    write(*,'(A)')'|      2d unstructured euler flow solver     |'
    write(*,'(A)')'|       Version 0.0.7 || 05/03/2026          |'
    write(*,'(A)')'|                 Max Wood                   |'
    write(*,'(A)')'|           University of Bristol            |'
    write(*,'(A)')'|    Department of Aerospace Engineering     |'
    write(*,'(A)')'+--------------------------------------------+'
    if (options%mode == 'solve') then !primal solve
        write(*,'(A)') '                {primal solve}                 '
    elseif (options%mode == 'adjoint') then !adjoint solve
        write(*,'(A)') '                {adjoint solve}                 '
    end if 
    write(*,'(A)') ' '
end if

!import mesh 
if (options%cdisplay) then
    write(*,'(A)') '--> importing mesh: '//trim(options%meshpath)//trim(options%meshname)
end if 
call import_mesh(mesh,options,options%meshpath)
if (options%cdisplay) then
    write(*,'(A,I0,A)') '    {ncell = ',mesh%ncell,'}' 
    write(*,'(A,I0,A)') '    {nedge = ',mesh%nedge,'}' 
    write(*,'(A,I0,A)') '    {nvertex = ',mesh%nvertex,'}' 
    write(*,'(A,E12.6,A,E12.6,A)') '    {cell volume (max/min) = ',maxval(mesh%cells_volume),' / ',&
    minval(mesh%cells_volume),'}' 
end if 

!initialise the flowfield
if (options%cdisplay) then
    write(*,'(A)') '--> initialising'
end if 
call flux_flow_initialise(mesh,options)

!mode switch 
if (options%mode == 'solve') then !primal solve

    !solve 
    if (options%cdisplay) then
        write(*,'(A)') '--> solving'
    end if 
    call flux_flow_solve(mesh,options)

    !post-process 
    ! call write_vtk(mesh,options,'flow.vtk')
    call write_vtu(mesh,options,'flow.vtu')
    call write_flow_field('flowfield',mesh)
    call write_forces('forces',mesh)
elseif (options%mode == 'adjoint') then !adjoint solve

    !UPDATE INTEGER TO IN64 FOR ALL METHODS INDEXING INTO SPARSE JACOBIANS TO AVOID OVERFLOW ERRORS
    !ADD ADJOINT DISSIPATION 

    !adjoint dev =====================

    !import solution
    call read_flow_field('flowfield',mesh)

    !adjoint solve 
    call flux_adjoint_solve(mesh,options)

    !post-process 
    ! call write_vtk(mesh,options,'flow.vtk')
    call write_vtu(mesh,options,'flow.vtu')
    call write_gradient('gradient',mesh)
    call write_forces('forces',mesh)
end if 
if (options%cdisplay) then
    print *, 'COMPLETE'
end if 
stop 
end program flux2d 
