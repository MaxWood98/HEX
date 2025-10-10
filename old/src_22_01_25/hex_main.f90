!hex mesh generator main program
!max wood
!version : 0.0.1
!updated : 20-12-24

!program 
program hex 
use hex2d
implicit none 

!variables 
type(hex_options) :: options 
type(halfedge) :: geometry
type(hex_mesh) :: mesh2d

!set default options 

options%cdisplay = .true.
options%geomname = 'test_geometries/shopt_0p1.fv'

options%ncell_max = 1000000
options%farfield_r = 20.0d0 
options%nrefine = 8
options%nflood = 6


!display
if (options%cdisplay) then
    write(*,'(A)') ' '
    write(*,'(A)')'+--------------------------------------------+'
    write(*,'(A)')'|                    hex                     |'
    write(*,'(A)')'|  2d/3d unstructured volume mesh generator  |'
    write(*,'(A)')'|        Version 0.0.1 || 20/12/2024         |'
    write(*,'(A)')'|                 Max Wood                   |'
    write(*,'(A)')'|           University of Bristol            |'
    write(*,'(A)')'|    Department of Aerospace Engineering     |'
    write(*,'(A)')'+--------------------------------------------+'
    write(*,'(A)') ' '
end if

!import geometry 
if (options%cdisplay) then 
    write(*,'(A)') '--> importing geometry '
end if 
call geometry%build_from_fv(options%geomname)
if (options%cdisplay) then 
    write(*,'(A,I0,A)') '    {geometry is ',geometry%ndim,' dimensional}'
    write(*,'(A,I0,A)') '    {imported = ',geometry%nvertex,' vertices}'
    write(*,'(A,I0,A)') '    {imported = ',geometry%nedge,' halfedges}'
    write(*,'(A,I0,A)') '    {imported = ',geometry%nface,' faces}'
end if 

!dimension switch 
if (geometry%ndim == 2) then 
    mesh2d = hex2d_mesh(geometry,options)
elseif (geometry%ndim == 3) then 

end if 
stop 
end program hex 