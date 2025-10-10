!hex gradient module 
!max wood
!version : 0.0.1
!updated : 20-02-25

!module 
module hex_gradient
use hex_geometry
use halfedge_mesh
use hex_data_methods

!routines
contains 

!gradient projection (chain) 2d ==============
subroutine project_gradient_mesh2geometry_chain_2d(options,geometry,mesh,surface_gradient)
implicit none 

!variables - inout
type(halfedge) :: geometry
type(hex_options) :: options 
type(hex_mesh), target :: mesh
real(dp), dimension(:,:), allocatable :: surface_gradient

!variables - local 
integer(in64) :: ii 
integer(in64) :: vs1,vs2
real(dp) :: ef

!allocate gradient 
allocate(surface_gradient(geometry%nvertex,2))
surface_gradient(:,:) = 0.0d0 

!evaluate mesh derivative 
if (options%gradient_method == 'mesh_derivative') then 
    call surface_mesh_derivative_2d(geometry,mesh)
end if 

!project gradient
do ii=1,mesh%nvertex
    if (mesh%vertex(ii)%flag) then 

        !edge fraction
        ef = mesh%vertex(ii)%rdata

        !target surface vertices 
        vs1 = mesh%vertex(ii)%ivdata(1)
        vs2 = mesh%vertex(ii)%ivdata(2)

        !project
        if (options%gradient_method == 'interpolate_linear') then 
            surface_gradient(vs1,1) = surface_gradient(vs1,1) + (1.0d0 - ef)*mesh%vertex(ii)%gradient(1)
            surface_gradient(vs1,2) = surface_gradient(vs1,2) + (1.0d0 - ef)*mesh%vertex(ii)%gradient(2)
            surface_gradient(vs2,1) = surface_gradient(vs2,1) + ef*mesh%vertex(ii)%gradient(1)
            surface_gradient(vs2,2) = surface_gradient(vs2,2) + ef*mesh%vertex(ii)%gradient(2)
        elseif (options%gradient_method == 'interpolate_rbf') then 

        elseif (options%gradient_method == 'mesh_derivative') then 
            surface_gradient(vs1,1) = surface_gradient(vs1,1) + mesh%vertex(ii)%mesh_gradient(1,1)*mesh%vertex(ii)%gradient(1)
            surface_gradient(vs1,1) = surface_gradient(vs1,1) + mesh%vertex(ii)%mesh_gradient(1,2)*mesh%vertex(ii)%gradient(2)
            surface_gradient(vs1,2) = surface_gradient(vs1,2) + mesh%vertex(ii)%mesh_gradient(1,3)*mesh%vertex(ii)%gradient(1)
            surface_gradient(vs1,2) = surface_gradient(vs1,2) + mesh%vertex(ii)%mesh_gradient(1,4)*mesh%vertex(ii)%gradient(2)
            surface_gradient(vs2,1) = surface_gradient(vs2,1) + mesh%vertex(ii)%mesh_gradient(2,1)*mesh%vertex(ii)%gradient(1)
            surface_gradient(vs2,1) = surface_gradient(vs2,1) + mesh%vertex(ii)%mesh_gradient(2,2)*mesh%vertex(ii)%gradient(2)
            surface_gradient(vs2,2) = surface_gradient(vs2,2) + mesh%vertex(ii)%mesh_gradient(2,3)*mesh%vertex(ii)%gradient(1)
            surface_gradient(vs2,2) = surface_gradient(vs2,2) + mesh%vertex(ii)%mesh_gradient(2,4)*mesh%vertex(ii)%gradient(2)
        end if 
    end if  
end do 
return
end subroutine project_gradient_mesh2geometry_chain_2d


!build surface mesh derivative 2d ==============
subroutine surface_mesh_derivative_2d(geometry,mesh)
implicit none 

!variables - inout
type(halfedge) :: geometry
type(hex_mesh), target :: mesh

!variables - local 
integer(in64) :: ii,jj
integer(in64) :: etgt
real(dp) :: nvedge
real(dp) :: v1m(2),v2m(2),v1g(2),v2g(2),grad_v1(4),grad_v2(4)

!build v2e
call mesh%index_edges()
call mesh%get_v2e()

!evaluate 
do ii=1,mesh%nvertex
    if (mesh%vertex(ii)%flag) then 

        !initialise 
        nvedge = 0.0d0 
        allocate(mesh%vertex(ii)%mesh_gradient(2,4)) !row1 = v1 | row2 = v2
        mesh%vertex(ii)%mesh_gradient(:,:) = 0.0d0 

        !loop surface attached edges
        do jj=1,mesh%max_valence
            if (mesh%v2e(ii,jj) .GT. 0) then
                etgt = mesh%v2e(ii,jj)
                if ((mesh%edge(etgt)%cell1 .NE. -1) .AND. (mesh%edge(etgt)%cell2 .NE. -1)) then 

                    !increment volume edge count
                    nvedge = nvedge + 1.00 

                    !geometry vertices
                    v1g = geometry%vertex(mesh%vertex(ii)%ivdata(1))%coordinate(1:2)
                    v2g = geometry%vertex(mesh%vertex(ii)%ivdata(2))%coordinate(1:2)

                    !mesh vertices
                    v1m = mesh%vertex(ii)%coordinate(1:2)
                    if (associated(mesh%edge(etgt)%vertex1,mesh%vertex(ii))) then 
                        v2m = mesh%edge(etgt)%vertex2%coordinate(1:2)
                    else
                        v2m = mesh%edge(etgt)%vertex1%coordinate(1:2)
                    end if 

                    !accumulate derivative
                    call line_line_intersect_gradient(v1g,v2g,v1m,v2m,grad_v1,grad_v2)
                    mesh%vertex(ii)%mesh_gradient(1,:) = mesh%vertex(ii)%mesh_gradient(1,:) + grad_v1
                    mesh%vertex(ii)%mesh_gradient(2,:) = mesh%vertex(ii)%mesh_gradient(2,:) + grad_v2
                end if 
            end if 
        end do 

        !normalise gradient 
        if (nvedge .NE. 0.0d0) then 
            mesh%vertex(ii)%mesh_gradient(:,:) = mesh%vertex(ii)%mesh_gradient(:,:)/nvedge
        end if 
    end if 
end do 
return 
end subroutine surface_mesh_derivative_2d

end module hex_gradient