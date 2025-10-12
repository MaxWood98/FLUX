!flux 2d - a cell based euler solver
!max wood
!version : 0.0.1
!updated : 28-03-25

!program 
program flux2d
use flux_io
use flux_solve
implicit none 

!variables 
type(flux_mesh) :: mesh 
type(flux_options) :: options 



!import options 
options%cdisplay = .true.
! options%meshpath = 'test_meshes/naca0012'
! options%meshpath = 'test_meshes/kh0p1'
! options%meshpath = 'test_meshes/cv8x8op'
! options%meshpath = 'test_meshes/cv8x8rnd'
options%meshpath = 'grid_cell'



options%aoadeg = 45.0d0 
options%machinf = 0.5d0 
options%gamma = 1.4d0 
options%R = 287.058d0
options%tinf = 288.0d0
options%rhoinf = 1.225d0 

options%niter_max = 10000
options%cfl = 1.5d0 

options%rk_niter = 4
options%k2 = 1.0d0 
options%k4 = 0.1d0

options%num_threads = 16
options%residual_convtol = -12.0 






allocate(options%rk_dissipation(options%rk_niter))
options%rk_dissipation(1) = .true.
options%rk_dissipation(2) = .false.
options%rk_dissipation(3) = .false.
options%rk_dissipation(4) = .false.




!display splash 
if (options%cdisplay) then
    write(*,'(A)') ' '
    write(*,'(A)')'+--------------------------------------------+'
    write(*,'(A)')'|                   flux 2d                  |'
    write(*,'(A)')'|      2d unstructured euler flow solver     |'
    write(*,'(A)')'|       Version 0.0.2 || 11/10/2025          |'
    write(*,'(A)')'|                 Max Wood                   |'
    write(*,'(A)')'|           University of Bristol            |'
    write(*,'(A)')'|    Department of Aerospace Engineering     |'
    write(*,'(A)')'+--------------------------------------------+'
    write(*,'(A)') ' '
end if

!import mesh 
if (options%cdisplay) then
    write(*,'(A)') '--> importing mesh: '//trim(options%meshpath)//trim(options%meshname)
end if 
call import_mesh(mesh,options%meshpath)
if (options%cdisplay) then
    write(*,'(A,I0,A)') '    {ncell = ',mesh%ncell,'}' 
    write(*,'(A,I0,A)') '    {nedge = ',mesh%nedge,'}' 
    write(*,'(A,I0,A)') '    {nvertex = ',mesh%nvertex,'}' 
    write(*,'(A,E12.6,A,E12.6,A)') '    {cell volume (max/min) = ',maxval(mesh%cells_volume),' / ',&
    minval(mesh%cells_volume),'}' 
end if 

!initialise the flow 
if (options%cdisplay) then
    write(*,'(A)') '--> initialising'
end if 
call flux_flow_initialise(mesh,options)


! ! call read_restart_file('flowfield',mesh,options)


!solve 
if (options%cdisplay) then
    write(*,'(A)') '--> solving'
end if 
call flux_flow_solve(mesh,options)
print *, 'COMPLETE'

!post-process 
call write_vtk(mesh,options,'flow.vtk')

! !export data 

stop 
end program flux2d 
