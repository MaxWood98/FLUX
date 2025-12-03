!flux 2d - a cell based euler solver
!max wood
!version : 0.0.1
!updated : 28-03-25

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
! call set_default_options(options)

!read command arguments 
call get_command_arguments(options)





!import options 
options%cdisplay = .true.
! options%meshpath = 'test_meshes/naca0012'
! options%meshpath = 'test_meshes/kh0p1'
! options%meshpath = 'test_meshes/cv8x8op'
! options%meshpath = 'test_meshes/cv8x8rnd'
options%meshpath = 'grid_cell_naca'



options%aoadeg = 0.0d0 
options%machinf = 0.95d0 
options%gamma = 1.4d0 
options%R = 287.058d0
options%tinf = 288.0d0
options%rhoinf = 1.225d0 

options%outflow_pratio = 0.8d0


options%niter_max = 10000
options%cfl = 4.0

options%rk_niter = 4
options%k2 = 0.005d0 
options%k4 = 0.001d0

options%num_threads = 8
options%residual_convtol = -8.0 

options%flux_method = 'vanleer'




allocate(options%rk_dissipation(options%rk_niter))
options%rk_dissipation(1) = .false.
options%rk_dissipation(2) = .false.
options%rk_dissipation(3) = .false.
options%rk_dissipation(4) = .false.




!display splash 
if (options%cdisplay) then
    write(*,'(A)') ' '
    write(*,'(A)')'+--------------------------------------------+'
    write(*,'(A)')'|                   flux 2d                  |'
    write(*,'(A)')'|      2d unstructured euler flow solver     |'
    write(*,'(A)')'|       Version 0.0.2 || 27/11/2025          |'
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
    print *, 'COMPLETE'

    !post-process 
    call write_vtk(mesh,options,'flow.vtk')
    call write_flow_field('flowfield',mesh)
elseif (options%mode == 'adjoint') then !adjoint solve

    !UPDATE INTEGER TO IN64 FOR ALL METHODS INDEXING INTO SPARSE JACOBIANS TO AVOID OVERFLOW ERRORS
    !ADD ADJOINT DISSIPATION 
    !UPDATE JACOBIAN ROUTINES TO USE AND PARALELLISM

    !adjoint dev =====================

    !import solution
    call read_flow_field('flowfield',mesh)

    !adjoint solve 
    call flux_adjoint_solve(mesh,options)
    print *, 'COMPLETE'

    !post-process 
    call write_vtk(mesh,options,'flow.vtk')
    call write_gradient('gradient',mesh)
end if 
stop 
end program flux2d 
