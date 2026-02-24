!hex mesh generator main program
!max wood
!version : 0.0.4
!updated : 24-02-26

!TODO =====================
! add initial fv geometry normals calculation to set halfedge normals such that it can deal with nested geometry 
!TODO =====================

!program 
program hex 
use hex2d
use hex_io
use hex_gradient
use hex_utilities
implicit none 

!variables 
logical :: is_selfintersecting
type(halfedge) :: geometry
type(facevertex) :: geometry_fv
type(hex_options) :: options 
type(hex_mesh), target :: mesh2d
real(dp), dimension(:,:), allocatable :: surface_gradient

!set default options 
call set_default_options(options)

!read command arguments 
call get_command_arguments(options)

!read options file
call read_hex_options(options,trim(options%optionspath)//trim(options%optionsname))

!display
if (options%cdisplay) then
    write(*,'(A)') ' '
    write(*,'(A)')'+--------------------------------------------+'
    write(*,'(A)')'|                    hex                     |'
    write(*,'(A)')'|  2d/3d unstructured volume mesh generator  |'
    write(*,'(A)')'|       Version 0.0.4 || 25/02/2026         |'
    write(*,'(A)')'|                 Max Wood                   |'
    write(*,'(A)')'|           University of Bristol            |'
    write(*,'(A)')'|    Department of Aerospace Engineering     |'
    write(*,'(A)')'+--------------------------------------------+'
    write(*,'(A)') ' '
end if

!options display
if (options%cdisplay) then 
    write(*,'(A)') '--> options summary: '
    write(*,'(A,A,A)') '    {meshing geometry ',options%mesh_in_out,'}'
    write(*,'(A,A,A)') '    {refinement tree type: ',options%mesh_treetype,'}'
    write(*,'(A,A,A)') '    {mesh-tree relation: ',options%tree_mesh_relation,'}'
    write(*,'(A,I0,A)') '    {refinement levels: ',options%nrefine,'}'
    write(*,'(A,I0,A)') '    {maximum number of cells: ',options%ncell_max,'}'
    write(*,'(A,I0,A)') '    {imported ',options%nrefzone,' refinement zones}'
    write(*,'(A,I0,A)') '    {imported ',options%nclipplane,' clipping planes}'
    write(*,'(A,I0,A)') '    {imported ',options%nbczone,' boundary condition tag zones}'
    write(*,'(A,I0,A)') '    {imported ',options%nbcremzone,' boundary condition remove connected zones}'
end if 

!import geometry 
if (options%cdisplay) then 
    write(*,'(A)') '--> importing geometry: '//trim(options%geompath)//trim(options%geomname)
end if 
call geometry_fv%read_fv_file(trim(options%geompath)//trim(options%geomname))
if (geometry_fv%nedge == 0) then !build edges if not supplied
    call geometry_fv%build_edges_from_faces()
end if 
if (geometry_fv%ndim == 2) then !initialise pseudo-faces if mesh is 2d
    call geometry_fv%build_edge_loop_faces()
end if 
call geometry%build_from_fv(geometry_fv)
if (geometry_fv%ndim == 2) then !if 2d then set the edge normals from the fv geometry 
    call set_2d_halfedge_normals_from_facevertex(geometry,geometry_fv)
end if 
if (options%cdisplay) then 
    write(*,'(A,I0,A)') '    {geometry is ',geometry%ndim,' dimensional}'
    write(*,'(A,I0,A)') '    {imported = ',geometry%nvertex,' vertices}'
    write(*,'(A,I0,A)') '    {imported = ',geometry%nedge,' halfedges}'
    write(*,'(A,I0,A)') '    {imported = ',geometry%nface,' faces}'
end if 

!process
if (options%mode == 'check') then 
    if (options%cdisplay) then
        write(*,'(A)') '--> checking for self intersections'
    end if
    if (geometry%ndim == 2) then 
        is_selfintersecting = is_self_intersecting(geometry,options)
        if (options%cdisplay) then
            write(*,'(A,L,A)') '    {is self intersecting: ',is_selfintersecting,'}'
        end if
    elseif (geometry%ndim == 3) then 
    
    end if 
    call write_geometry_check(is_selfintersecting,options)
elseif (options%mode == 'mesh') then 
    if (geometry%ndim == 2) then 
        mesh2d = hex2d_mesh(geometry,options)
        if (options%cdisplay) then 
            write(*,'(A)') '--> exporting mesh'
        end if 
        call write_hex_cell_mesh_2d(mesh2d,options%meshpath//options%meshname)
        if (options%export_vtu) then 
            if (options%cdisplay) then 
                write(*,'(A)') '--> exporting mesh vtu'
            end if 
            call write_mesh_vtu_2d(mesh2d,options%meshpath//options%meshname//'.vtu')
        end if 
    elseif (geometry%ndim == 3) then 
    
    end if 
elseif (options%mode == 'project') then 

    !import data
    if (options%cdisplay) then 
        write(*,'(A)') '--> importing mesh: '//trim(options%meshpath)//trim(options%meshname)
    end if 
    call read_hex_cell_mesh_2d(mesh2d,options%meshpath//options%meshname)
    if (options%cdisplay) then
        write(*,'(A,I0,A)') '    {cells: ',mesh2d%ncell,'}'
        write(*,'(A,I0,A)') '    {edges: ',mesh2d%nedge,'}'
        write(*,'(A,I0,A)') '    {vertices: ',mesh2d%nvertex,'}'
        write(*,'(A,I0,A)') '    {surface intersection vertices: ',mesh2d%nvertex_surfint,'}'
    end if 
    if (options%cdisplay) then 
        write(*,'(A)') '--> importing gradient: '//trim(options%gradientpath)//trim(options%gradientname)
    end if 
    call read_gradient(mesh2d,options%gradientpath//options%gradientname)

    !project gradient
    if (options%cdisplay) then 
        write(*,'(A)') '--> projecting gradient'
    end if 
    call project_gradient_mesh2geometry_chain_2d(options,geometry,mesh2d,surface_gradient)

    !write surface gradients 
    call write_gradient_2d(surface_gradient,trim(options%gradientpath)//'surface_'//trim(options%gradientname))
end if 





!testing write cell based mesh
! if (options%cdisplay) then 
!     write(*,'(A)') '--> exporting mesh'
! end if 

!DEBUG-TESTING write mesh 
! call write_hex_mesh_2d(mesh2d,'grid')
! call write_hex_cell_mesh_2d(mesh2d,mesh2d,options%meshpath//options%meshname)



!display complete
if (options%cdisplay) then 
    write(*,'(A)') '    {complete}'
end if 
stop 
end program hex 