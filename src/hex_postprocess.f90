!hex postprocess module
!max wood
!version : 0.0.6
!updated : 22-06-25

!module 
module hex_postprocess
use hex_geometry
use hex_data_methods

!routines 
contains 


!remove cells on target boundary condition =========================
subroutine remove_cells_on_boundary_condition(mesh,boundary_condition)
implicit none 

!variables - inout 
integer(in64) :: boundary_condition
type(hex_mesh), target :: mesh

!variables - local
integer(in64) :: ii,ee,cc,ff
integer(in64) :: nfront,nfrontn,etgt,vindex1,vindex2,cindex2
integer(in64), dimension(:), allocatable :: cfront,cfrontn

!build the mesh cells 
call mesh%get_cells()

!flag cells on the target boundary condition and build the front 
nfront = 0 
allocate(cfront(mesh%ncell))
allocate(cfrontn(mesh%ncell))
do ee=1,mesh%nedge 
    if ((mesh%edge(ee)%cell1 == boundary_condition) .OR. (mesh%edge(ee)%cell2 == boundary_condition)) then 
        if (mesh%edge(ee)%cell1 .GT. 0) then 
            if (.NOT. mesh%cell(mesh%edge(ee)%cell1)%flag) then 
                mesh%cell(mesh%edge(ee)%cell1)%flag = .true.
                nfront = nfront + 1
                cfront(nfront) = mesh%edge(ee)%cell1
            end if 
        end if 
        if (mesh%edge(ee)%cell2 .GT. 0) then 
            if (.NOT. mesh%cell(mesh%edge(ee)%cell2)%flag) then 
                mesh%cell(mesh%edge(ee)%cell2)%flag = .true.
                nfront = nfront + 1
                cfront(nfront) = mesh%edge(ee)%cell2
            end if 
        end if
    end if 
end do

!flood the flag
do ff=1,2*mesh%ncell
    nfrontn = 0
    do cc=1,nfront
        do ee=1,mesh%cell(cfront(cc))%nedge
            etgt = mesh%cell(cfront(cc))%edges(ee)
            if (mesh%edge(etgt)%cell1 .GT. 0) then 
                if (.NOT. mesh%cell(mesh%edge(etgt)%cell1)%flag) then 
                    nfrontn = nfrontn + 1
                    cfrontn(nfrontn) = mesh%edge(etgt)%cell1
                    mesh%cell(mesh%edge(etgt)%cell1)%flag = .true.
                end if 
            end if 
            if (mesh%edge(etgt)%cell2 .GT. 0) then 
                if (.NOT. mesh%cell(mesh%edge(etgt)%cell2)%flag) then 
                    nfrontn = nfrontn + 1
                    cfrontn(nfrontn) = mesh%edge(etgt)%cell2
                    mesh%cell(mesh%edge(etgt)%cell2)%flag = .true.
                end if 
            end if 
        end do 
    end do 
    if (nfrontn == 0) then 
        exit
    end if 
    cfront(1:nfrontn) = cfrontn(1:nfrontn)
    nfront = nfrontn
end do 

!tag vertices and edges to remove 
do ii=1,mesh%nvertex
    mesh%vertex(ii)%flag = .false.
end do 
do ii=1,mesh%nedge
    mesh%edge(ii)%flag = .false.
end do 
do cc=1,mesh%ncell
    if (mesh%cell(cc)%flag) then 
        do ee=1,mesh%cell(cc)%nedge
            etgt = mesh%cell(cc)%edges(ee)
            mesh%edge(etgt)%flag = .true.
            mesh%edge(etgt)%vertex1%flag = .true.
            mesh%edge(etgt)%vertex2%flag = .true.
        end do 
    end if 
end do 

!unflag items on any retained cells 
do cc=1,mesh%ncell
    if (.NOT. mesh%cell(cc)%flag) then 
        do ee=1,mesh%cell(cc)%nedge
            etgt = mesh%cell(cc)%edges(ee)
            mesh%edge(etgt)%flag = .false.
            mesh%edge(etgt)%vertex1%flag = .false.
            mesh%edge(etgt)%vertex2%flag = .false.
        end do 
    end if 
end do 

!update the boundary condition on new boundary edges to farfield and re-orient them if neccesary
do ee=1,mesh%nedge
    if ((mesh%edge(ee)%cell1 .GT. 0) .AND. (mesh%edge(ee)%cell2 .GT. 0)) then 
        if ((mesh%cell(mesh%edge(ee)%cell1)%flag) .AND. (.NOT. mesh%cell(mesh%edge(ee)%cell2)%flag)) then !reasign boundary condition and flip 
            cindex2 = mesh%edge(ee)%cell2
            vindex1 = mesh%edge(ee)%vertex1%index
            vindex2 = mesh%edge(ee)%vertex2%index
            mesh%edge(ee)%cell1 = cindex2
            mesh%edge(ee)%cell2 = -2
            mesh%edge(ee)%vertex1 => mesh%vertex(vindex2)
            mesh%edge(ee)%vertex2 => mesh%vertex(vindex1)
        elseif((mesh%cell(mesh%edge(ee)%cell2)%flag) .AND. (.NOT. mesh%cell(mesh%edge(ee)%cell1)%flag)) then !reasign boundary condition
            mesh%edge(ee)%cell2 = -2
        end if 
    end if
end do 

!remove edges
call remove_flagged_edges(mesh)

!remove vertices
call remove_flagged_vertices(mesh)

!reindex mesh 
call mesh%index_vertices()
call mesh%index_edges()
call mesh%index_cells()
return 
end subroutine remove_cells_on_boundary_condition


!set boundary condition zone =========================
subroutine set_boundary_condition_inzone(mesh,zone)
implicit none 

!variables - inout 
type(hex_bczone) :: zone
type(hex_mesh), target :: mesh

!variables - local
integer(in64) :: ee
real(dp) :: xmin,xmax,ymin,ymax

!set edges 
do ee=1,mesh%nedge
    if ((mesh%edge(ee)%cell1 == -2) .OR. (mesh%edge(ee)%cell2 == -2)) then 
        xmin = min(mesh%edge(ee)%vertex1%coordinate(1),mesh%edge(ee)%vertex2%coordinate(1))
        xmax = max(mesh%edge(ee)%vertex1%coordinate(1),mesh%edge(ee)%vertex2%coordinate(1))
        ymin = min(mesh%edge(ee)%vertex1%coordinate(2),mesh%edge(ee)%vertex2%coordinate(2))
        ymax = max(mesh%edge(ee)%vertex1%coordinate(2),mesh%edge(ee)%vertex2%coordinate(2))
        if ((xmin .LE. zone%xmax) .AND. (xmax .GE. zone%xmin)) then 
            if ((ymin .LE. zone%ymax) .AND. (ymax .GE. zone%ymin)) then 
                if (mesh%edge(ee)%cell1 == -2) then 
                    mesh%edge(ee)%cell1 = zone%bctag
                end if 
                if (mesh%edge(ee)%cell2 == -2) then 
                    mesh%edge(ee)%cell2 = zone%bctag
                end if
            end if 
        end if
    end if 
end do 
return 
end subroutine set_boundary_condition_inzone


!clip to plane =========================
subroutine clip_to_plane(mesh,vp1,vp2)
implicit none 

!variables - inout 
real(dp) :: vp1(2),vp2(2)
type(hex_mesh), target :: mesh

!variables - local 
logical :: is_clipped
integer(in64) :: ii,ee
integer(in64) :: vindex1,vindex2,cindex2
real(dp) :: pnorm(2),pref(2),vref(2),cmid(3)

!set the plane normal an reference point 
call get_point_normal_plane_from_two_points(vp1,vp2,pref,pnorm)

!get cells 
call mesh%index_vertices()
call mesh%index_edges()
call mesh%get_cells()
call mesh%get_cell_edges()

!flag cells to remove
do ii=1,mesh%ncell
    mesh%cell(ii)%flag = .false.
    cmid = mesh%cell(ii)%get_midpoint(mesh) 
    vref = cmid(1:2) - pref
    vref = vref/(norm2(vref) + 1e-12_dp)
    if (dot_product(vref,pnorm) .LT. 0.0d0) then 
        mesh%cell(ii)%flag = .true.
    end if 
end do 

!flag edges and vertices to remove 
do ii=1,mesh%nvertex
    mesh%edge(ii)%flag = .false.
end do 
do ii=1,mesh%nedge
    mesh%edge(ii)%flag = .false.
end do 
do ii=1,mesh%ncell
    if (mesh%cell(ii)%flag) then
        do ee=1,mesh%cell(ii)%nedge
            mesh%edge(mesh%cell(ii)%edges(ee))%flag = .true.
            mesh%edge(mesh%cell(ii)%edges(ee))%vertex1%flag = .true.
            mesh%edge(mesh%cell(ii)%edges(ee))%vertex2%flag = .true.
        end do 
    end if 
end do
do ii=1,mesh%ncell !retain items shared with retained cells
    if (.NOT. mesh%cell(ii)%flag) then
        do ee=1,mesh%cell(ii)%nedge
            mesh%edge(mesh%cell(ii)%edges(ee))%flag = .false.
            mesh%edge(mesh%cell(ii)%edges(ee))%vertex1%flag = .false.
            mesh%edge(mesh%cell(ii)%edges(ee))%vertex2%flag = .false.
        end do 
    end if
end do 

!update positions of the vertices on clipped edges 
call mesh%get_v2e()
do ii=1,mesh%nedge
    if ((mesh%edge(ii)%cell1 .GT. 0) .AND. (mesh%edge(ii)%cell2 .GT. 0)) then 
        is_clipped = .false.
        if ((mesh%cell(mesh%edge(ii)%cell1)%flag) .AND. (.NOT. mesh%cell(mesh%edge(ii)%cell2)%flag)) then
            is_clipped = .true.
        elseif ((mesh%cell(mesh%edge(ii)%cell2)%flag) .AND. (.NOT. mesh%cell(mesh%edge(ii)%cell1)%flag)) then
            is_clipped = .true.
        end if 
        if (is_clipped) then 
            call snap_vertex_to_plane(mesh,mesh%edge(ii)%vertex1,vp1,vp2)
            call snap_vertex_to_plane(mesh,mesh%edge(ii)%vertex2,vp1,vp2)
        end if 
    end if 
end do 

!update the boundary conditions of the clipped edges
do ii=1,mesh%nedge
    if ((mesh%edge(ii)%cell1 .GT. 0) .AND. (mesh%edge(ii)%cell2 .GT. 0)) then 
        if ((mesh%cell(mesh%edge(ii)%cell1)%flag) .AND. (.NOT. mesh%cell(mesh%edge(ii)%cell2)%flag)) then !reasign connectivity and flip 
            cindex2 = mesh%edge(ii)%cell2
            vindex1 = mesh%edge(ii)%vertex1%index
            vindex2 = mesh%edge(ii)%vertex2%index
            mesh%edge(ii)%cell1 = cindex2
            mesh%edge(ii)%cell2 = -2
            mesh%edge(ii)%vertex1 => mesh%vertex(vindex2)
            mesh%edge(ii)%vertex2 => mesh%vertex(vindex1)
        elseif ((mesh%cell(mesh%edge(ii)%cell2)%flag) .AND. (.NOT. mesh%cell(mesh%edge(ii)%cell1)%flag)) then !reasign connectivity
             mesh%edge(ii)%cell2 = -2
        end if 
    end if 
end do 

!remove edges
call remove_flagged_edges(mesh)

!remove vertices
call remove_flagged_vertices(mesh)

!reindex mesh 
call mesh%index_vertices()
call mesh%index_edges()
call mesh%index_cells()
return 
end subroutine clip_to_plane


!snap vertex to plane =========================
subroutine snap_vertex_to_plane(mesh,vertex,vp1,vp2)
implicit none 

!variables - inout 
real(dp) :: vp1(2),vp2(2)
type(hex_vertex) :: vertex
type(hex_mesh), target :: mesh

!variables - local 
integer(in64) :: ee 
integer(in64) :: surface_edge_target,vtx_target
integer(in64), dimension(:), allocatable :: surface_edges
real(dp) :: pnorm(2),pref(2),vref(2),vnew(2)

!get the plane reference point and normal 
call get_point_normal_plane_from_two_points(vp1,vp2,pref,pnorm)

!snap 
if (vertex%is_on_surface(mesh)) then 
    surface_edges = vertex%get_surface_edges(mesh)
    surface_edge_target = 0 
    do ee=1,size(surface_edges,1)
        if (mesh%edge(surface_edges(ee))%vertex1%index == vertex%index) then 
            vtx_target = mesh%edge(surface_edges(ee))%vertex2%index 
        else
            vtx_target = mesh%edge(surface_edges(ee))%vertex1%index 
        end if 
        vref = mesh%vertex(vtx_target)%coordinate(1:2) - pref
        vref = vref/(norm2(vref) + 1e-12_dp)
        if (dot_product(vref,pnorm) .GT. 0.0d0) then 
            surface_edge_target = surface_edges(ee)
            exit 
        end if 
    end do 
    if (surface_edge_target .NE. 0) then !located a retained surface edge -> intersect this edge with the plane and snap to this point
        vnew = line_line_intersect(vp1,vp2,mesh%edge(surface_edge_target)%vertex1%coordinate(1:2),&
        mesh%edge(surface_edge_target)%vertex2%coordinate(1:2))
        vertex%coordinate(1:2) = vnew
    else !no retained surface edge located -> snap to the closest point on the plane
        vnew = closest_point_on_line_to_point(vp1,vp2,vertex%coordinate(1:2))
        vertex%coordinate(1:2) = vnew
    end if 
else
    vnew = closest_point_on_line_to_point(vp1,vp2,vertex%coordinate(1:2))
    vertex%coordinate(1:2) = vnew
end if 
return 
end subroutine snap_vertex_to_plane


!remove double boundary condition edges =========================
subroutine remove_double_bc_edges(mesh,ndoublebc)
implicit none 

!variables - inout 
integer(in64) :: ndoublebc
type(hex_mesh), target :: mesh

!variables - local 
integer(in64) :: ii

!set edge flags
do ii=1,mesh%nedge
    mesh%edge(ii)%flag = .false.
end do 

!flag double bounday condition edges 
ndoublebc = 0 
do ii=1,mesh%nedge
    if ((mesh%edge(ii)%cell1 .LT. 0) .AND. (mesh%edge(ii)%cell2 .LT. 0)) then 
        mesh%edge(ii)%flag = .true.
        ndoublebc = ndoublebc + 1 
    end if 
end do 

!remove edges if neccesary
if (ndoublebc .NE. 0) then 

    !remove flagged edges 
    call remove_flagged_edges(mesh)

    !flag unused vertices 
    do ii=1,mesh%nvertex
        mesh%vertex(ii)%flag = .true. 
    end do
    do ii=1,mesh%nedge
        mesh%edge(ii)%vertex1%flag = .false.
        mesh%edge(ii)%vertex2%flag = .false.
    end do

    !remove flagged vertices 
    call remove_flagged_vertices(mesh)

    !index edges 
    call mesh%index_edges() 
end if 
return 
end subroutine remove_double_bc_edges


!remove unlinked surface edges =========================
subroutine remove_unlinked_surface_edges(mesh,nunassociated)
implicit none 

!variables - inout 
integer(in64) :: nunassociated
type(hex_mesh), target :: mesh

!variables - local 
integer(in64) :: ii

!set edge flags
do ii=1,mesh%nedge
    mesh%edge(ii)%flag = .false.
end do 

!flag unassociated surface edges 
nunassociated = 0 
do ii=1,mesh%nedge
    if ((mesh%edge(ii)%cell1 == -1) .OR. (mesh%edge(ii)%cell2 == -1)) then 
        if ((mesh%edge(ii)%cell1 == 0) .OR. (mesh%edge(ii)%cell2 == 0)) then 
            mesh%edge(ii)%flag = .true.
            nunassociated = nunassociated + 1
        end if 
    end if 
end do 

!remove edges if neccesary
if (nunassociated .NE. 0) then 

    !remove flagged edges 
    call remove_flagged_edges(mesh)

    !flag unused vertices 
    do ii=1,mesh%nvertex
        mesh%vertex(ii)%flag = .true. 
    end do
    do ii=1,mesh%nedge
        mesh%edge(ii)%vertex1%flag = .false.
        mesh%edge(ii)%vertex2%flag = .false.
    end do

    !remove flagged vertices 
    call remove_flagged_vertices(mesh)

    !index edges 
    call mesh%index_edges() 
end if 
return 
end subroutine remove_unlinked_surface_edges


!far field smoothing =========================
subroutine far_field_smoothing(mesh,options)
implicit none 

!variables - inout 
type(hex_options) :: options 
type(hex_mesh), target :: mesh
    
!variables - local 
integer(in64) :: ii,aa,ss
integer(in64) :: vlncc
real(dp) :: vlnccr
real(dp) :: coordinate_new(mesh%nvertex,2)

!set vertex indecies
call mesh%index_vertices()

!get mesh v2e
call mesh%get_v2e()

!set vertex tags and flags
do ii=1,mesh%nvertex
    mesh%vertex(ii)%tag = 0 
    mesh%vertex(ii)%flag = .false.
end do 

!tag far field vertices 
do ii=1,mesh%nedge
    if ((mesh%edge(ii)%cell1 == -2) .OR. (mesh%edge(ii)%cell2 == -2)) then  
        mesh%edge(ii)%vertex1%tag = 1
        mesh%edge(ii)%vertex2%tag = 1
    end if
end do 

!flag valence 2 far field vertices that are not on geometry surfaces
do ii=1,mesh%nvertex
    if (mesh%vertex(ii)%tag == 1) then 
        vlncc = 0 
        do aa=1,mesh%max_valence
            if (mesh%v2e(ii,aa) .NE. 0) then 
                vlncc = vlncc + 1
            else
                exit 
            end if 
        end do 
        if (vlncc == 2) then 
            if (.NOT. mesh%vertex(ii)%is_on_surface(mesh)) then 
                mesh%vertex(ii)%flag = .true.
            end if 
        end if 
    end if 
end do 

!smooth flagged far-field vertices 
coordinate_new(:,:) = 0.0d0 
do ss=1,options%nsmooth_farfield 
    do ii=1,mesh%nvertex
        if (mesh%vertex(ii)%flag) then 
            vlnccr = 0.0d0 
            coordinate_new(ii,1:2) = 0.0d0  
            do aa=1,mesh%max_valence
                if (mesh%v2e(ii,aa) .NE. 0) then 
                    if (mesh%edge(mesh%v2e(ii,aa))%vertex1%index == ii) then 
                        if (mesh%vertex(mesh%edge(mesh%v2e(ii,aa))%vertex2%index)%tag == 1) then 
                            coordinate_new(ii,1:2) = coordinate_new(ii,1:2) + mesh%edge(mesh%v2e(ii,aa))%vertex2%coordinate(1:2)
                            vlnccr = vlnccr + 1.0d0 
                        end if 
                    else
                        if (mesh%vertex(mesh%edge(mesh%v2e(ii,aa))%vertex1%index)%tag == 1) then 
                            coordinate_new(ii,1:2) = coordinate_new(ii,1:2) + mesh%edge(mesh%v2e(ii,aa))%vertex1%coordinate(1:2)
                            vlnccr = vlnccr + 1.0d0 
                        end if 
                    end if 
                end if 
            end do 
            coordinate_new(ii,1:2) = coordinate_new(ii,1:2)/vlnccr
        end if 
    end do 
    do ii=1,mesh%nvertex
        if (mesh%vertex(ii)%flag) then 
            mesh%vertex(ii)%coordinate(1:2) = coordinate_new(ii,1:2)
        end if
    end do 
end do 
return 
end subroutine far_field_smoothing


!inter-layer smoothing =========================
subroutine inter_layer_smoothing(mesh,options)
implicit none 

!variables - inout 
type(hex_options) :: options 
type(hex_mesh), target :: mesh
    
!variables - local 
integer(in64) :: ii,ee,ss
real(dp) :: dv,dv_tol,valencei
real(dp) :: coordinate_new(mesh%nvertex,2)

!set layer volume delta tollerance
dv_tol = 1e-12

!assign edges to current cells 
call mesh%get_cells()

!get cell volumes 
call mesh%get_cell_areas()

!get mesh v2e
call mesh%get_v2e()

!get mesh v2c
call mesh%get_v2c()

!get cell midpoints
do ii=1,mesh%ncell
    mesh%cell(ii)%midpoint = mesh%cell(ii)%get_midpoint(mesh)
end do 

!set cell flags and tags
do ii=1,mesh%ncell
    mesh%cell(ii)%tag = 0
    mesh%cell(ii)%flag = .false.
end do 

!flag boundary contition cells 
do ii=1,mesh%nedge
    if (mesh%edge(ii)%cell1 .LT. 0) then 
        mesh%cell(mesh%edge(ii)%cell2)%flag = .true.
    elseif (mesh%edge(ii)%cell2 .LT. 0) then 
        mesh%cell(mesh%edge(ii)%cell1)%flag = .true.
    end if 
end do 

!tag vertices in inter-layer cells 
do ii=1,mesh%nvertex
    mesh%vertex(ii)%tag = 0 
    mesh%vertex(ii)%flag = .false.
end do 
do ii=1,mesh%nedge
    if ((mesh%edge(ii)%cell1 .GT. 0) .AND. (mesh%edge(ii)%cell2 .GT. 0)) then  
        if ((mesh%cell(mesh%edge(ii)%cell1)%flag) .OR. (mesh%cell(mesh%edge(ii)%cell2)%flag)) then  
            cycle 
        end if 
        dv = abs(mesh%cell(mesh%edge(ii)%cell2)%volume - mesh%cell(mesh%edge(ii)%cell1)%volume)
        if (dv .GT. dv_tol) then 
            mesh%cell(mesh%edge(ii)%cell1)%tag = 1
            mesh%cell(mesh%edge(ii)%cell2)%tag = 1
            do ee=1,mesh%cell(mesh%edge(ii)%cell1)%nedge
                mesh%edge(mesh%cell(mesh%edge(ii)%cell1)%edges(ee))%vertex1%flag = .true.
                mesh%edge(mesh%cell(mesh%edge(ii)%cell1)%edges(ee))%vertex2%flag = .true.
            end do 
            do ee=1,mesh%cell(mesh%edge(ii)%cell2)%nedge
                mesh%edge(mesh%cell(mesh%edge(ii)%cell2)%edges(ee))%vertex1%flag = .true.
                mesh%edge(mesh%cell(mesh%edge(ii)%cell2)%edges(ee))%vertex2%flag = .true.
            end do 
        end if 
    end if 
end do 

!tag and un-flag surface and far field vertices
do ii=1,mesh%nedge
    if ((mesh%edge(ii)%cell1 .LT. 0) .OR. (mesh%edge(ii)%cell2 .LT. 0)) then  
        mesh%edge(ii)%vertex1%tag = 2
        mesh%edge(ii)%vertex2%tag = 2
        mesh%edge(ii)%vertex1%flag = .false.
        mesh%edge(ii)%vertex2%flag = .false.
    end if
end do 

! !flood vertex flags if requested
! do ss=1,1
!     do ii=1,mesh%nvertex
!         if (mesh%vertex(ii)%flag) then 
!             mesh%vertex(ii)%tag = 1 
!         end if  
!     end do 
!     do ii=1,mesh%nvertex
!         if ((mesh%vertex(ii)%flag) .AND. (mesh%vertex(ii)%tag == 1)) then 
!             do ee=1,mesh%max_valence
!                 if (mesh%v2e(ii,ee) .NE. 0) then 
!                     if (mesh%edge(mesh%v2e(ii,ee))%vertex1%index == ii) then 
!                         if (mesh%vertex(mesh%edge(mesh%v2e(ii,ee))%vertex2%index)%tag == 0) then 
!                             mesh%edge(mesh%v2e(ii,ee))%vertex2%flag = .true.
!                         end if 
!                     else
!                         if (mesh%vertex(mesh%edge(mesh%v2e(ii,ee))%vertex1%index)%tag == 0) then 
!                             mesh%edge(mesh%v2e(ii,ee))%vertex1%flag = .true.
!                         end if 
!                     end if 
!                 end if 
!             end do 
!         end if 
!     end do 
! end do 

!smooth mesh
coordinate_new(:,:) = 0.0d0 
do ss=1,options%nsmooth_interlayer
    do ii=1,mesh%nvertex
        if (mesh%vertex(ii)%flag) then 
            valencei = 0.0d0 
            coordinate_new(ii,:) = 0.0d0 
            do ee=1,mesh%max_cvalence
                if (mesh%v2c(ii,ee) .NE. 0) then 
                    valencei = valencei + mesh%cell(mesh%v2c(ii,ee))%volume
                    coordinate_new(ii,:) = coordinate_new(ii,:) + &
                    mesh%cell(mesh%v2c(ii,ee))%midpoint(1:2)*mesh%cell(mesh%v2c(ii,ee))%volume
                end if 
            end do 
            coordinate_new(ii,:) = coordinate_new(ii,:)/valencei
        end if 
    end do 
    do ii=1,mesh%nvertex
        if (mesh%vertex(ii)%flag) then 
            mesh%vertex(ii)%coordinate(1:2) = coordinate_new(ii,:)
        end if 
    end do 
    do ii=1,mesh%ncell
        if (mesh%cell(ii)%tag == 1) then 
            mesh%cell(ii)%midpoint = mesh%cell(ii)%get_midpoint(mesh)
        end if 
    end do 
    call mesh%get_cell_areas() !UPDATE TO ONLY DO ON TARGETED CELLS FOR EFFICIENCY
end do 
return 
end subroutine inter_layer_smoothing


!collapse small or degenerate cells =========================
subroutine collapse_small_or_degenerate_cells(mesh,options,ntagged,method)
implicit none 

!variables - inout 
character(*) :: method
integer(in64) :: ntagged
type(hex_options) :: options 
type(hex_mesh), target :: mesh

!variables - local 
logical :: cell_invalid
integer(in64) :: ii,ee
integer(in64) :: etgt,ecell,cnew
real(dp) :: etgt_len,ecell_len

!get cells 
call mesh%get_cells()

!get cell areas
call mesh%get_cell_areas()

!method switch 
if (method == 'small') then !check for and flag small cells 
    ntagged = 0 
    do ii=1,mesh%ncell
        mesh%cell(ii)%flag = .false.
        if (mesh%cell(ii)%volume .LE. options%cellvol_min) then 
            mesh%cell(ii)%flag = .true.
            ntagged = ntagged + 1
        end if  
    end do 
elseif (method == 'degenerate') then !check for and flag degenerate cells 
    ! call mesh%get_cell_edges()
    do ii=1,mesh%ncell
        mesh%cell(ii)%flag = .false.
        if (mesh%cell(ii)%nedge .LE. 2) then 
            mesh%cell(ii)%flag = .true.
            ntagged = ntagged + 1
        end if 
    end do 
end if 

!eliminate tagged cells
if (ntagged .NE. 0) then 

    !set edge flags
    do ii=1,mesh%nedge
        mesh%edge(ii)%flag = .false. 
    end do 

    ! !assign edges to current cells 
    ! call mesh%get_cell_edges()

    !tag edges to eliminate to merge cells 
    ntagged = 0 
    do ii=1,mesh%ncell
        if (mesh%cell(ii)%flag) then 

            !skip this cell if any edges within it are already tagged to remove 
            cell_invalid = .false.
            do ee=1,mesh%cell(ii)%nedge
                ecell = mesh%cell(ii)%edges(ee)
                if (mesh%edge(ecell)%flag) then 
                    cell_invalid = .true.
                    exit 
                end if
            end do 
            if (cell_invalid) then 
                cycle 
            end if 
            
            !find the longest edge on this cell that points to a valid cell to eliminate
            etgt = 0 
            cnew = 0
            etgt_len = 0.0d0 
            do ee=1,mesh%cell(ii)%nedge
                ecell = mesh%cell(ii)%edges(ee)
                if (.NOT. mesh%edge(ecell)%flag) then !not already removed
                    if ((mesh%edge(ecell)%cell1 .GT. 0) .AND. (mesh%edge(ecell)%cell2 .GT. 0)) then !internal edge 
                        if (mesh%edge(ecell)%cell1 == ii) then  
                            if (mesh%cell(mesh%edge(ecell)%cell2)%flag) then !skip as adjacent cell is also small
                                cycle
                            end if 
                        elseif (mesh%edge(ecell)%cell2 == ii) then 
                            if (mesh%cell(mesh%edge(ecell)%cell1)%flag) then !skip as adjacent cell is also small
                                cycle
                            end if 
                        end if 
                        ecell_len = norm2(mesh%edge(ecell)%vertex2%coordinate(1:2) - mesh%edge(ecell)%vertex1%coordinate(1:2))
                        if (ecell_len .GT. etgt_len) then 
                            etgt = ecell 
                            etgt_len = ecell_len
                            if (mesh%edge(ecell)%cell1 == ii) then 
                                cnew = mesh%edge(ecell)%cell2
                            else
                                cnew = mesh%edge(ecell)%cell1
                            end if 
                        end if 
                    end if 
                end if 
            end do 

            !tag the edge to remove and reassign cell adjacencies in the current cell if a valid edge has been found 
            if (etgt .NE. 0) then 

                !increment removed cell count 
                ntagged = ntagged + 1

                !remove target edge
                mesh%edge(etgt)%flag = .true.

                !reassign cell tags for edges in cell ii 
                do ee=1,mesh%cell(ii)%nedge
                    ecell = mesh%cell(ii)%edges(ee)
                    if (mesh%edge(ecell)%cell1 == ii) then 
                        mesh%edge(ecell)%cell1 = cnew
                    end if
                    if (mesh%edge(ecell)%cell2 == ii) then 
                        mesh%edge(ecell)%cell2 = cnew
                    end if
                end do 
            end if 
        end if 
    end do 

    !remove flagged edges 
    call remove_flagged_edges(mesh)
end if 
return 
end subroutine collapse_small_or_degenerate_cells


!collapse degenerate edges =========================
subroutine collapse_degenerate_edges(mesh,ndegenrate)
implicit none 

!variables - inout 
integer(in64) :: ndegenrate
type(hex_mesh), target :: mesh

!variables - local 
integer(in64) :: ii

!check for and flag degenerate edges 
ndegenrate = 0 
do ii=1,mesh%nedge
    mesh%edge(ii)%flag = .false. 
end do 
do ii=1,mesh%nedge
    if ((mesh%edge(ii)%vertex1%index == mesh%edge(ii)%vertex2%index)) then 
        mesh%edge(ii)%flag = .true.
        ndegenrate = ndegenrate + 1
    end if 
end do 

!collapse degenerate edges
if (ndegenrate .NE. 0) then 
    call remove_flagged_edges(mesh)
end if 
return 
end subroutine collapse_degenerate_edges


!eliminate short edges =========================
subroutine collapse_short_edges(mesh,options,nshort)
implicit none 

!variables - inout 
integer(in64) :: nshort
type(hex_options) :: options 
type(hex_mesh), target :: mesh

!variables - local 
integer(in64) :: ii,aa
integer(in64) :: etgt,ncollapse
real(dp) :: elen
type(hex_vertex), pointer :: vretain,vremove

!index edges
call mesh%index_edges()

!check for and flag short edges 
nshort = 0 
do ii=1,mesh%nedge
    mesh%edge(ii)%flag = .false. 
end do 
do ii=1,mesh%nedge
    elen = norm2(mesh%edge(ii)%vertex2%coordinate(1:2) - mesh%edge(ii)%vertex1%coordinate(1:2))
    if (elen .LE. options%edgelength_min) then 
        mesh%edge(ii)%flag = .true. 
        nshort = nshort + 1
    end if 
end do 

!collapse short edges
ncollapse = 0 
if (nshort .NE. 0) then 

    !get v2e
    call mesh%get_v2e()

    !tag surface vertices
    do ii=1,mesh%nvertex
        mesh%vertex(ii)%flag = .false. 
    end do
    do ii=1,mesh%nedge
        if ((mesh%edge(ii)%cell1 == -1) .OR. (mesh%edge(ii)%cell2 == -1)) then 
            mesh%edge(ii)%vertex1%flag = .true.
            mesh%edge(ii)%vertex2%flag = .true.
        end if 
    end do 

    !reset vertex tags
    do ii=1,mesh%nvertex
        mesh%vertex(ii)%tag = 0 
    end do

    !collapse short edges 
    do ii=1,mesh%nedge
        if (mesh%edge(ii)%flag) then 

            !skip if either vertex is tagged
            if ((mesh%edge(ii)%vertex1%tag == 1) .OR. (mesh%edge(ii)%vertex2%tag == 1)) then 
                cycle 
            end if 

            !tag both vertices on this edge 
            mesh%edge(ii)%vertex1%tag = 1 
            mesh%edge(ii)%vertex2%tag = 1 

            !select vertex to retain and vertex to remove 
            if (mesh%edge(ii)%vertex1%flag) then 
                vretain => mesh%edge(ii)%vertex1
                vremove => mesh%edge(ii)%vertex2
            elseif (mesh%edge(ii)%vertex2%flag) then 
                vretain => mesh%edge(ii)%vertex2
                vremove => mesh%edge(ii)%vertex1
            else
                vretain => mesh%edge(ii)%vertex1
                vremove => mesh%edge(ii)%vertex2
            end if 

            !update the end vertex of all edges on vremove
            do aa=1,mesh%max_valence
                if (mesh%v2e(vremove%index,aa) .NE. 0) then 
                    etgt = mesh%v2e(vremove%index,aa)
                    if (mesh%edge(etgt)%vertex1%index == vremove%index) then 
                        mesh%edge(etgt)%vertex1 => mesh%vertex(vretain%index)
                    elseif (mesh%edge(etgt)%vertex2%index == vremove%index) then 
                        mesh%edge(etgt)%vertex2 => mesh%vertex(vretain%index)
                    end if 
                end if 
            end do 

            !count collapsed edges 
            ncollapse = ncollapse + 1
        end if 
    end do 

    !check for any unflagged edges with the same vertex assigned to both ends and flag them to remove 
    do ii=1,mesh%nedge
        if ((mesh%edge(ii)%vertex1%index == mesh%edge(ii)%vertex2%index) .AND. (.NOT. mesh%edge(ii)%flag)) then 
            mesh%edge(ii)%flag = .true.
        end if 
    end do 

    !remove flagged edges 
    call remove_flagged_edges(mesh)

    !flag unused vertices 
    do ii=1,mesh%nvertex
        mesh%vertex(ii)%flag = .true. 
    end do
    do ii=1,mesh%nedge
        mesh%edge(ii)%vertex1%flag = .false.
        mesh%edge(ii)%vertex2%flag = .false.
    end do

    !remove flagged vertices 
    call remove_flagged_vertices(mesh)

    !set return collapsed count 
    nshort = ncollapse

    !index edges 
    call mesh%index_edges() 
end if 
return 
end subroutine collapse_short_edges


!partition bisected cells =========================
subroutine partition_bisected_cells(mesh,nbisected)
implicit none 

!variables - inout 
integer(in64) :: nbisected
type(hex_mesh), target :: mesh

!variables - local 
logical :: bisected
integer(in64) :: ii,ee,ff,aa,ncellN
integer(in64) :: eadj,nupdate

!get v2e
call mesh%index_edges() 
call mesh%get_v2e()

!get cells 
call mesh%get_cells()
ncellN = mesh%ncell

!set all edge tags and flags
do ii=1,mesh%nedge
    mesh%edge(ii)%tag = 0 
    mesh%edge(ii)%flag = .false.
end do 

!check for and partition any bisected cells 
nbisected = 0
do ii=1,mesh%ncell

    !initialise 
    bisected = .false.

    !set all edge flags and tags in this cell 
    do ee=1,mesh%cell(ii)%nedge
        mesh%edge(mesh%cell(ii)%edges(ee))%tag = 1
        mesh%edge(mesh%cell(ii)%edges(ee))%flag = .false.
    end do 

    !find the edge valence for each vertex in this cell 
    do ee=1,mesh%cell(ii)%nedge
        mesh%edge(mesh%cell(ii)%edges(ee))%vertex1%idata = 0 
        mesh%edge(mesh%cell(ii)%edges(ee))%vertex2%idata = 0 
    end do 
    do ee=1,mesh%cell(ii)%nedge
        mesh%edge(mesh%cell(ii)%edges(ee))%vertex1%idata = mesh%edge(mesh%cell(ii)%edges(ee))%vertex1%idata + 1
        mesh%edge(mesh%cell(ii)%edges(ee))%vertex2%idata = mesh%edge(mesh%cell(ii)%edges(ee))%vertex2%idata + 1
    end do 

    !flood all edges in this cell from the first
    mesh%edge(mesh%cell(ii)%edges(1))%flag = .true.
    do ff=1,4*mesh%cell(ii)%nedge

        !tag adjacent edges in this cell 
        nupdate = 0 
        do ee=1,mesh%cell(ii)%nedge
            if (mesh%edge(mesh%cell(ii)%edges(ee))%flag) then
                if (mesh%edge(mesh%cell(ii)%edges(ee))%vertex1%idata .LE. 2) then !exclude pinch points
                    do aa=1,mesh%max_valence !across v1
                        eadj = mesh%v2e(mesh%edge(mesh%cell(ii)%edges(ee))%vertex1%index,aa)
                        if (eadj .GT. 0) then 
                            if ((mesh%edge(eadj)%tag == 1) .AND. (.NOT. mesh%edge(eadj)%flag)) then 
                                mesh%edge(eadj)%flag = .true.
                                nupdate = nupdate + 1
                            end if 
                        end if 
                    end do  
                end if 
                if (mesh%edge(mesh%cell(ii)%edges(ee))%vertex2%idata .LE. 2) then !exclude pinch points
                    do aa=1,mesh%max_valence !across v2
                        eadj = mesh%v2e(mesh%edge(mesh%cell(ii)%edges(ee))%vertex2%index,aa)
                        if (eadj .GT. 0) then 
                            if ((mesh%edge(eadj)%tag == 1) .AND. (.NOT. mesh%edge(eadj)%flag)) then 
                                mesh%edge(eadj)%flag = .true.
                                nupdate = nupdate + 1
                            end if 
                        end if 
                    end do  
                end if 
            end if 
        end do 

        !exit if no edges tagged
        if (nupdate == 0) then 
            exit 
        end if 
    end do 

    !tag cell as bisected if not all edges have been tagged during the flood 
    do ee=1,mesh%cell(ii)%nedge
        if (.NOT. mesh%edge(mesh%cell(ii)%edges(ee))%flag) then 
            bisected = .true.
            exit 
        end if 
    end do 

    !partition if bisected 
    if (bisected) then

        !increment bisected count 
        nbisected = nbisected + 1

        !increment cell count 
        ncellN = ncellN + 1

        !update cell adjacency on all flag = .false. edges in this cell
        do ee=1,mesh%cell(ii)%nedge
            if (.NOT. mesh%edge(mesh%cell(ii)%edges(ee))%flag) then 
                if (mesh%edge(mesh%cell(ii)%edges(ee))%cell1 == ii) then 
                    mesh%edge(mesh%cell(ii)%edges(ee))%cell1 = ncellN
                end if 
                if (mesh%edge(mesh%cell(ii)%edges(ee))%cell2 == ii) then 
                    mesh%edge(mesh%cell(ii)%edges(ee))%cell2 = ncellN
                end if 
            end if 
        end do 
    end if 

    !reset all edge tags and flags in the original cell 
    do ee=1,mesh%cell(ii)%nedge
        mesh%edge(mesh%cell(ii)%edges(ee))%tag = 1
        mesh%edge(mesh%cell(ii)%edges(ee))%flag = .false.
    end do 
end do 

!update cell indexing 
call mesh%index_cells()
return 
end subroutine partition_bisected_cells


!remove flagged edges =========================
subroutine remove_flagged_edges(mesh)
implicit none 

!variables - inout 
type(hex_mesh), target :: mesh

!variables - local 
integer(in64) :: ii 
integer(in64) :: nedgeN
type(hex_edge) :: edge_temp(size(mesh%edge,dim=1)) 

!count retained edges
nedgeN = 0
do ii=1,mesh%nedge
    if (.NOT. mesh%edge(ii)%flag) then 
        nedgeN = nedgeN + 1
    end if 
end do 

!eliminate flagged edges from the mesh 
edge_temp = mesh%edge
deallocate(mesh%edge)
allocate(mesh%edge(nedgeN))
nedgeN = 0 
do ii=1,mesh%nedge
    if (.NOT. edge_temp(ii)%flag) then 
        nedgeN = nedgeN + 1
        mesh%edge(nedgeN) = edge_temp(ii)
    end if  
end do 
mesh%nedge = nedgeN

!reindex edges 
call mesh%index_edges()
return 
end subroutine remove_flagged_edges


!remove flagged vertices =========================
subroutine remove_flagged_vertices(mesh)
implicit none 

!variables - inout 
type(hex_mesh), target :: mesh

!variables - local 
integer(in64) :: ii 
integer(in64) :: nvertexN
integer(in64) :: edge_vertices(size(mesh%edge,dim=1),2)
type(hex_vertex) :: vertex_temp(size(mesh%vertex,dim=1)) 

!count retained vertices and reindex them 
nvertexN = 0
do ii=1,mesh%nvertex
    if (.NOT. mesh%vertex(ii)%flag) then 
        nvertexN = nvertexN + 1
        mesh%vertex(ii)%index = nvertexN
    else
        mesh%vertex(ii)%index = 0
    end if 
end do 

!store edge vertices 
edge_vertices(:,:) = 0 
do ii=1,mesh%nedge
    edge_vertices(ii,1) = mesh%edge(ii)%vertex1%index 
    edge_vertices(ii,2) = mesh%edge(ii)%vertex2%index
end do 

!reallocate mesh vertices
vertex_temp = mesh%vertex
deallocate(mesh%vertex)
allocate(mesh%vertex(nvertexN))
do ii=1,mesh%nvertex
    if (.NOT. vertex_temp(ii)%flag) then 
        mesh%vertex(vertex_temp(ii)%index)%external = vertex_temp(ii)%external 
        mesh%vertex(vertex_temp(ii)%index)%flag = vertex_temp(ii)%flag 
        mesh%vertex(vertex_temp(ii)%index)%tag = vertex_temp(ii)%tag 
        mesh%vertex(vertex_temp(ii)%index)%index = vertex_temp(ii)%index 
        mesh%vertex(vertex_temp(ii)%index)%coordinate = vertex_temp(ii)%coordinate 
    end if 
end do 

!re-assign edge pointers
do ii=1,mesh%nedge
    mesh%edge(ii)%vertex1 => mesh%vertex(edge_vertices(ii,1))
    mesh%edge(ii)%vertex2 => mesh%vertex(edge_vertices(ii,2))
end do 

!set new mesh vertex count 
mesh%nvertex = nvertexN
return 
end subroutine remove_flagged_vertices


end module hex_postprocess