!hex tritree module
!max wood
!version : 0.0.1
!updated : 09-02-25

!module 
module hex_tritree
use halfedge_mesh
use geometry_kdtree
use hex_data_methods
use hex_geometry

!tritree type
type tritree
integer(in64) :: nnode,nvertex
integer(in64), dimension(:), allocatable :: vertex_tag
real(dp), dimension(:,:), allocatable :: vertices 
type(trinode), dimension(:), allocatable :: node
contains 
    procedure :: refine_node => refine_node_tri
end type tritree

!tri node type
type trinode
logical :: flag = .false.,visited = .false.
integer(in64) :: level = 1 
integer(in64) :: index = 0
integer(in64) :: vertex1 = 0,vertex2 = 0,vertex3 = 0
integer(in64) :: evertex1 = 0,evertex2 = 0,evertex3 = 0
integer(in64), pointer :: aevertex1 => null(),aevertex2 => null(),aevertex3 => null()
type(trinode), pointer :: parent => null()
type(trinode), pointer :: child1 => null()
type(trinode), pointer :: child2 => null()
type(trinode), pointer :: child3 => null()
type(trinode), pointer :: child4 => null()
type(trinode), pointer :: adjacent1 => null()
type(trinode), pointer :: adjacent2 => null()
type(trinode), pointer :: adjacent3 => null()
contains 
    procedure :: map_adjacent_node => map_adjacent_node_tri
    procedure :: cascade_adjacency => cascade_adjacency_tri
end type trinode

!methods ==================================================
contains 


!build full mesh (primal) =========================
subroutine build_full_mesh_from_tritree_primal(mesh_full,tri_tree)
implicit none 

!variables - import
type(hex_mesh), target :: mesh_full
type(tritree), target :: tri_tree 

!variables - local 
logical :: edge_valid
integer(in64) :: ii,eins,vins 
integer(in64) :: v2v(tri_tree%nvertex,6)
type(hex_edge), dimension(:), allocatable :: edge_temp 


!build vertices 
mesh_full%nvertex = tri_tree%nvertex
allocate(mesh_full%vertex(mesh_full%nvertex))
do ii=1,tri_tree%nvertex
    mesh_full%vertex(ii)%index = ii 
    mesh_full%vertex(ii)%coordinate = tri_tree%vertices(ii,:)
    mesh_full%vertex(ii)%coordinate(3) = 0.0d0 
end do 

!index retained nodes 
vins = 0 
do ii=1,tri_tree%nnode 
    if (.NOT. associated(tri_tree%node(ii)%child1)) then 
        vins = vins + 1
        tri_tree%node(ii)%index = vins 
    end if 
end do 

!build and index edges 
v2v(:,:) = 0
allocate(edge_temp(tri_tree%nnode*6))
eins = 0 
do ii=1,tri_tree%nnode 
    if (.NOT. associated(tri_tree%node(ii)%child1)) then 
        edge_valid = .false.
        if (associated(tri_tree%node(ii)%adjacent1)) then 
            if (.NOT. associated(tri_tree%node(ii)%adjacent1%child1)) then 
                edge_valid = .true.
            end if 
        else
            edge_valid = .true.
        end if 
        if (is_connected_v2v(v2v(tri_tree%node(ii)%vertex1,:),tri_tree%node(ii)%vertex2)) then 
            edge_valid = .false.
        end if 
        if (edge_valid) then 

            eins = eins + 1
            edge_temp(eins)%index = eins

            edge_temp(eins)%vertex2 => mesh_full%vertex(tri_tree%node(ii)%vertex1)
            edge_temp(eins)%vertex1 => mesh_full%vertex(tri_tree%node(ii)%vertex2)
            edge_temp(eins)%cell1 = tri_tree%node(ii)%index
            if (associated(tri_tree%node(ii)%adjacent1)) then
                edge_temp(eins)%cell2 = tri_tree%node(ii)%adjacent1%index
            else
                edge_temp(eins)%cell2 = -2
            end if 
            call add_connection(v2v(tri_tree%node(ii)%vertex1,:),v2v(tri_tree%node(ii)%vertex2,:),&
                                tri_tree%node(ii)%vertex1,tri_tree%node(ii)%vertex2) 
        end if 
        edge_valid = .false.
        if (associated(tri_tree%node(ii)%adjacent2)) then 
            if (.NOT. associated(tri_tree%node(ii)%adjacent2%child1)) then 
                edge_valid = .true.
            end if 
        else
            edge_valid = .true.
        end if 
        if (is_connected_v2v(v2v(tri_tree%node(ii)%vertex2,:),tri_tree%node(ii)%vertex3)) then 
            edge_valid = .false.
        end if 
        if (edge_valid) then 

            eins = eins + 1
            edge_temp(eins)%index = eins

            edge_temp(eins)%vertex2 => mesh_full%vertex(tri_tree%node(ii)%vertex2)
            edge_temp(eins)%vertex1 => mesh_full%vertex(tri_tree%node(ii)%vertex3)
            edge_temp(eins)%cell1 = tri_tree%node(ii)%index
            if (associated(tri_tree%node(ii)%adjacent2)) then
                edge_temp(eins)%cell2 = tri_tree%node(ii)%adjacent2%index
            else
                edge_temp(eins)%cell2 = -2
            end if 
            call add_connection(v2v(tri_tree%node(ii)%vertex2,:),v2v(tri_tree%node(ii)%vertex3,:),&
                                tri_tree%node(ii)%vertex2,tri_tree%node(ii)%vertex3) 
        end if 
        edge_valid = .false.
        if (associated(tri_tree%node(ii)%adjacent3)) then 
            if (.NOT. associated(tri_tree%node(ii)%adjacent3%child1)) then 
                edge_valid = .true.
            end if 
        else
            edge_valid = .true.
        end if 
        if (is_connected_v2v(v2v(tri_tree%node(ii)%vertex3,:),tri_tree%node(ii)%vertex1)) then 
            edge_valid = .false.
        end if 
        if (edge_valid) then 

            eins = eins + 1
            edge_temp(eins)%index = eins

            edge_temp(eins)%vertex2 => mesh_full%vertex(tri_tree%node(ii)%vertex3)
            edge_temp(eins)%vertex1 => mesh_full%vertex(tri_tree%node(ii)%vertex1)
            edge_temp(eins)%cell1 = tri_tree%node(ii)%index
            if (associated(tri_tree%node(ii)%adjacent3)) then
                edge_temp(eins)%cell2 = tri_tree%node(ii)%adjacent3%index
            else
                edge_temp(eins)%cell2 = -2
            end if 
            call add_connection(v2v(tri_tree%node(ii)%vertex3,:),v2v(tri_tree%node(ii)%vertex1,:),&
                                tri_tree%node(ii)%vertex3,tri_tree%node(ii)%vertex1) 
        end if 
    end if 
end do 
mesh_full%nedge = eins

!trim edges 
allocate(mesh_full%edge(eins))
mesh_full%edge(:) = edge_temp(1:eins)

!set edge directions 
do ii=1,eins
    mesh_full%edge(ii)%direction = mesh_full%edge(ii)%vertex2%coordinate - mesh_full%edge(ii)%vertex1%coordinate
    mesh_full%edge(ii)%direction = mesh_full%edge(ii)%direction/norm2(mesh_full%edge(ii)%direction)
end do 
return 
end subroutine build_full_mesh_from_tritree_primal


!build full mesh (dual) =========================
subroutine build_full_mesh_from_tritree_dual(mesh_full,tri_tree)
implicit none 

!variables - import
type(hex_mesh), target :: mesh_full
type(tritree), target :: tri_tree 

!variables - local 
integer(in64) :: ii,eins,vins 
integer(in64) :: v2v(tri_tree%nnode,6)
integer(in64) :: edge_swap(4)
type(hex_edge), dimension(:), allocatable :: edge_temp 

!build vertices
vins = 0 
do ii=1,tri_tree%nnode 
    if (.NOT. associated(tri_tree%node(ii)%child1)) then 
        vins = vins + 1
        tri_tree%node(ii)%index = vins 
    end if 
end do 
mesh_full%nvertex = vins
vins = 0 
allocate(mesh_full%vertex(mesh_full%nvertex))
do ii=1,tri_tree%nnode 
    if (.NOT. associated(tri_tree%node(ii)%child1)) then 
        vins = vins + 1
        mesh_full%vertex(vins)%index = vins 
        mesh_full%vertex(vins)%coordinate = (tri_tree%vertices(tri_tree%node(ii)%vertex1,:) + &
                                             tri_tree%vertices(tri_tree%node(ii)%vertex2,:) + &
                                             tri_tree%vertices(tri_tree%node(ii)%vertex3,:))/3.0d0 
        mesh_full%vertex(vins)%coordinate(3) = 0.0d0 
    end if 
end do 

!tag far-field adjacent vertices in the tritree
do ii=1,tri_tree%nnode 
    if (.NOT. associated(tri_tree%node(ii)%child1)) then 
        if (.NOT. associated(tri_tree%node(ii)%adjacent1)) then 
            tri_tree%vertex_tag(tri_tree%node(ii)%vertex1) = -2
            tri_tree%vertex_tag(tri_tree%node(ii)%vertex2) = -2
        end if
        if (.NOT. associated(tri_tree%node(ii)%adjacent2)) then 
            tri_tree%vertex_tag(tri_tree%node(ii)%vertex2) = -2
            tri_tree%vertex_tag(tri_tree%node(ii)%vertex3) = -2
        end if
        if (.NOT. associated(tri_tree%node(ii)%adjacent3)) then 
            tri_tree%vertex_tag(tri_tree%node(ii)%vertex3) = -2
            tri_tree%vertex_tag(tri_tree%node(ii)%vertex1) = -2
        end if
    end if
end do 

!build and index edges 
v2v(:,:) = 0
allocate(edge_temp(tri_tree%nnode*6))
eins = 0 
do ii=1,tri_tree%nnode 
    if (.NOT. associated(tri_tree%node(ii)%child1)) then 
        if (associated(tri_tree%node(ii)%adjacent1)) then 
            if (.NOT. associated(tri_tree%node(ii)%adjacent1%child1)) then 
                if (.NOT. is_connected_v2v(v2v(tri_tree%node(ii)%index,:),tri_tree%node(ii)%adjacent1%index)) then 
                    eins = eins + 1
                    edge_temp(eins)%index = eins
                    edge_temp(eins)%vertex1 => mesh_full%vertex(tri_tree%node(ii)%index)
                    edge_temp(eins)%vertex2 => mesh_full%vertex(tri_tree%node(ii)%adjacent1%index)
                    
                    if (tri_tree%vertex_tag(tri_tree%node(ii)%vertex1) < 0) then 
                        edge_temp(eins)%cell1 = tri_tree%vertex_tag(tri_tree%node(ii)%vertex1)
                    else    
                        edge_temp(eins)%cell1 = tri_tree%node(ii)%vertex1
                    end if 
                    if (tri_tree%vertex_tag(tri_tree%node(ii)%vertex2) < 0) then 
                        edge_temp(eins)%cell2 = tri_tree%vertex_tag(tri_tree%node(ii)%vertex2)
                    else
                        edge_temp(eins)%cell2 = tri_tree%node(ii)%vertex2
                    end if 
                    call add_connection(v2v(tri_tree%node(ii)%index,:),v2v(tri_tree%node(ii)%adjacent1%index,:),&
                                        tri_tree%node(ii)%index,tri_tree%node(ii)%adjacent1%index) 
                    if (edge_temp(eins)%cell1 < 0) then !orient far-field edge
                        edge_swap(1) = edge_temp(eins)%vertex1%index
                        edge_swap(2) = edge_temp(eins)%vertex2%index
                        edge_swap(3) = edge_temp(eins)%cell1
                        edge_swap(4) = edge_temp(eins)%cell2
                        edge_temp(eins)%vertex1 => mesh_full%vertex(edge_swap(2))
                        edge_temp(eins)%vertex2 => mesh_full%vertex(edge_swap(1))
                        edge_temp(eins)%cell1 = edge_swap(4)
                        edge_temp(eins)%cell2 = edge_swap(3)
                    end if 
                end if 
            end if 
        end if 
        if (associated(tri_tree%node(ii)%adjacent2)) then 
            if (.NOT. associated(tri_tree%node(ii)%adjacent2%child1)) then 
                if (.NOT. is_connected_v2v(v2v(tri_tree%node(ii)%index,:),tri_tree%node(ii)%adjacent2%index)) then 
                    eins = eins + 1
                    edge_temp(eins)%index = eins
                    edge_temp(eins)%vertex1 => mesh_full%vertex(tri_tree%node(ii)%index)
                    edge_temp(eins)%vertex2 => mesh_full%vertex(tri_tree%node(ii)%adjacent2%index)

                    if (tri_tree%vertex_tag(tri_tree%node(ii)%vertex2) < 0) then 
                        edge_temp(eins)%cell1 = tri_tree%vertex_tag(tri_tree%node(ii)%vertex2)
                    else
                        edge_temp(eins)%cell1 = tri_tree%node(ii)%vertex2
                    end if 
                    if (tri_tree%vertex_tag(tri_tree%node(ii)%vertex3) < 0) then 
                        edge_temp(eins)%cell2 = tri_tree%vertex_tag(tri_tree%node(ii)%vertex3)
                    else
                        edge_temp(eins)%cell2 = tri_tree%node(ii)%vertex3
                    end if 
                    call add_connection(v2v(tri_tree%node(ii)%index,:),v2v(tri_tree%node(ii)%adjacent2%index,:),&
                                        tri_tree%node(ii)%index,tri_tree%node(ii)%adjacent2%index) 
                    if (edge_temp(eins)%cell1 < 0) then !orient far-field edge
                        edge_swap(1) = edge_temp(eins)%vertex1%index
                        edge_swap(2) = edge_temp(eins)%vertex2%index
                        edge_swap(3) = edge_temp(eins)%cell1
                        edge_swap(4) = edge_temp(eins)%cell2
                        edge_temp(eins)%vertex1 => mesh_full%vertex(edge_swap(2))
                        edge_temp(eins)%vertex2 => mesh_full%vertex(edge_swap(1))
                        edge_temp(eins)%cell1 = edge_swap(4)
                        edge_temp(eins)%cell2 = edge_swap(3)
                    end if
                end if
            end if 
        end if 
        if (associated(tri_tree%node(ii)%adjacent3)) then 
            if (.NOT. associated(tri_tree%node(ii)%adjacent3%child1)) then 
                if (.NOT. is_connected_v2v(v2v(tri_tree%node(ii)%index,:),tri_tree%node(ii)%adjacent3%index)) then 
                    eins = eins + 1
                    edge_temp(eins)%index = eins
                    edge_temp(eins)%vertex1 => mesh_full%vertex(tri_tree%node(ii)%index)
                    edge_temp(eins)%vertex2 => mesh_full%vertex(tri_tree%node(ii)%adjacent3%index)

                    if (tri_tree%vertex_tag(tri_tree%node(ii)%vertex3) < 0) then 
                        edge_temp(eins)%cell1 = tri_tree%vertex_tag(tri_tree%node(ii)%vertex3)
                    else
                        edge_temp(eins)%cell1 = tri_tree%node(ii)%vertex3
                    end if 
                    if (tri_tree%vertex_tag(tri_tree%node(ii)%vertex1) < 0) then 
                        edge_temp(eins)%cell2 = tri_tree%vertex_tag(tri_tree%node(ii)%vertex1)
                    else
                        edge_temp(eins)%cell2 = tri_tree%node(ii)%vertex1
                    end if 
                    call add_connection(v2v(tri_tree%node(ii)%index,:),v2v(tri_tree%node(ii)%adjacent3%index,:),&
                                        tri_tree%node(ii)%index,tri_tree%node(ii)%adjacent3%index) 
                    if (edge_temp(eins)%cell1 < 0) then !orient far-field edge
                        edge_swap(1) = edge_temp(eins)%vertex1%index
                        edge_swap(2) = edge_temp(eins)%vertex2%index
                        edge_swap(3) = edge_temp(eins)%cell1
                        edge_swap(4) = edge_temp(eins)%cell2
                        edge_temp(eins)%vertex1 => mesh_full%vertex(edge_swap(2))
                        edge_temp(eins)%vertex2 => mesh_full%vertex(edge_swap(1))
                        edge_temp(eins)%cell1 = edge_swap(4)
                        edge_temp(eins)%cell2 = edge_swap(3)
                    end if
                end if
            end if 
        end if 
    end if 
end do 
mesh_full%nedge = eins

!trim edges 
allocate(mesh_full%edge(eins))
mesh_full%edge(:) = edge_temp(1:eins)

!set edge directions 
do ii=1,eins
    mesh_full%edge(ii)%direction = mesh_full%edge(ii)%vertex2%coordinate - mesh_full%edge(ii)%vertex1%coordinate
    mesh_full%edge(ii)%direction = mesh_full%edge(ii)%direction/norm2(mesh_full%edge(ii)%direction)
end do 
return 
end subroutine build_full_mesh_from_tritree_dual


!add connection =========================
subroutine add_connection(v2v_v1,v2v_v2,v1,v2) 
implicit none 

!variables - inout
integer(in64) :: v1,v2
integer(in64), dimension(:) :: v2v_v1,v2v_v2

!variables - local 
integer(in64) :: ii

!add
do ii=1,size(v2v_v1)
    if (v2v_v1(ii) == 0) then
        v2v_v1(ii) = v2 
        exit 
    end if 
end do 
do ii=1,size(v2v_v2)
    if (v2v_v2(ii) == 0) then
        v2v_v2(ii) = v1 
        exit 
    end if 
end do 
return 
end subroutine add_connection


!is connected check =========================
function is_connected_v2v(v2v_v1,v2) result(connected)
implicit none 

!variables - inout
logical :: connected
integer(in64) :: v2
integer(in64), dimension(:) :: v2v_v1

!variables - local 
integer(in64) :: ii

!check
connected = .false.
do ii=1,size(v2v_v1)
    if (v2v_v1(ii) == v2) then
        connected = .true.
        exit 
    end if 
end do 
return 
end function is_connected_v2v


!build mesh tritree =========================
subroutine build_tritree(options,tri_tree,kdtree) 
implicit none 

!variables - inout 
type(tritree), target :: tri_tree 
type(gkdtree), target :: kdtree
type(hex_options) :: options 

!variables - local 
logical :: refzone_limited,inzone
integer(in64) :: ii,rr,ff,nn,zz
integer(in64) :: nnodeC,nadjf,ntgtnode,nupdate,nrefine,reffloodidx
integer(in64) :: nref_flood(max(options%nrefine,1))
real(dp) :: xmax,xmin,ymax,ymin,xmid,ymid,radius,rf !,xmaxg,xming,ymaxg,yming
real(dp) :: vid(3)
type(edge), pointer :: edgec
type(node) :: tgt_nodes(kdtree%nnode)

!initialise tree structure
call initialise_tritree_hexagon(tri_tree,options)

!initialise vertex tags
tri_tree%vertex_tag(:) = 0 

!interpolate number of refinement flood iterations at each refinement level 
nref_flood(:) = 0 
if (options%nrefine .GT. 0) then 
    do rr=1,options%nrefine 
        if (rr .LT. nint(real(options%nrefine,dp)/2.0)) then 
            rf = real(rr,dp)/(real(options%nrefine,dp)/2.0)
            nref_flood(rr) = nint(real(options%nflood_mid,dp)*rf + real(options%nflood_coarse,dp)*(1.0 - rf))
        elseif (rr == nint(real(options%nrefine,dp)/2.0)) then 
            nref_flood(rr) = options%nflood_mid
        elseif (rr .GT. nint(real(options%nrefine,dp)/2.0)) then 
            rf = (real(rr,dp) - (real(options%nrefine,dp)/2.0))/(real(options%nrefine,dp)/2.0)
            nref_flood(rr) = nint(real(options%nflood_fine,dp)*rf + real(options%nflood_mid,dp)*(1.0 - rf))
        end if 
    end do 
    nref_flood(1) = options%nflood_coarse
    nref_flood(options%nrefine) = options%nflood_fine
end if 

!refine tree 
do rr=1,10*options%nrefine 

    !reset flags
    do ii=1,tri_tree%nnode
        tri_tree%node(ii)%flag = .false.
        tri_tree%node(ii)%visited = .false.
    end do 

    !select nodes to refine
    do ii=1,tri_tree%nnode 
        if (tri_tree%node(ii)%level == rr) then 

            !cell bounding box and midpoint
            xmax = max(tri_tree%vertices(tri_tree%node(ii)%vertex1,1),&
            tri_tree%vertices(tri_tree%node(ii)%vertex2,1),&
            tri_tree%vertices(tri_tree%node(ii)%vertex3,1))
            xmin = min(tri_tree%vertices(tri_tree%node(ii)%vertex1,1),&
            tri_tree%vertices(tri_tree%node(ii)%vertex2,1),&
            tri_tree%vertices(tri_tree%node(ii)%vertex3,1))
            ymax = max(tri_tree%vertices(tri_tree%node(ii)%vertex1,2),&
            tri_tree%vertices(tri_tree%node(ii)%vertex2,2),&
            tri_tree%vertices(tri_tree%node(ii)%vertex3,2))
            ymin = min(tri_tree%vertices(tri_tree%node(ii)%vertex1,2),&
            tri_tree%vertices(tri_tree%node(ii)%vertex2,2),&
            tri_tree%vertices(tri_tree%node(ii)%vertex3,2))
            xmid = (1.0d0/3.0d0)*(tri_tree%vertices(tri_tree%node(ii)%vertex1,1) + &
            tri_tree%vertices(tri_tree%node(ii)%vertex2,1) + tri_tree%vertices(tri_tree%node(ii)%vertex3,1))
            ymid = (1.0d0/3.0d0)*(tri_tree%vertices(tri_tree%node(ii)%vertex1,2) + &
            tri_tree%vertices(tri_tree%node(ii)%vertex2,2) + tri_tree%vertices(tri_tree%node(ii)%vertex3,2))

            !flag nodes to refine ======
            !check if this node is within a refinement zone 
            refzone_limited = .false.
            if (options%nrefzone .NE. 0) then 
                do zz=1,options%nrefzone
                    inzone = .false.
                    if (options%refinement_zones(zz)%type == 'point') then 
                        radius = sqrt((options%refinement_zones(zz)%xmid - xmid)**2 + (options%refinement_zones(zz)%ymid - ymid)**2)
                        if (radius .LE. options%refinement_zones(zz)%radius) then 
                            inzone = .true.
                        else
                            if ((options%refinement_zones(zz)%xmid .GE. xmin) .AND. (options%refinement_zones(zz)%xmid .LE. xmax)) then
                                if ((options%refinement_zones(zz)%ymid .GE. ymin) .AND. (options%refinement_zones(zz)%ymid .LE. ymax)) then
                                    inzone = .true.
                                end if 
                            end if 
                        end if 
                    elseif (options%refinement_zones(zz)%type == 'line') then 
                        vid = min_dist_point_to_line(options%refinement_zones(zz)%v1(1:2),options%refinement_zones(zz)%v2(1:2),(/xmid,ymid/))
                        radius = vid(3)
                        if (radius .LE. options%refinement_zones(zz)%radius) then 
                            inzone = .true.
                        elseif (line_triangle_intersect_bool(options%refinement_zones(zz)%v1(1:2),options%refinement_zones(zz)%v2(1:2),&
                                tri_tree%vertices(tri_tree%node(ii)%vertex1,1:2),&
                                tri_tree%vertices(tri_tree%node(ii)%vertex2,1:2),&
                                tri_tree%vertices(tri_tree%node(ii)%vertex3,1:2))) then
                            inzone = .true.
                        end if 
                    elseif (options%refinement_zones(zz)%type == 'quad') then 
                        if ((options%refinement_zones(zz)%xmax .GE. xmin) .AND. (options%refinement_zones(zz)%xmin .LE. xmax)) then
                            if ((options%refinement_zones(zz)%ymax .GE. ymin) .AND. (options%refinement_zones(zz)%ymin .LE. ymax)) then
                                inzone = .true.
                            end if 
                        end if 
                    end if 
                    if (inzone) then 
                        if (tri_tree%node(ii)%level .LE. options%refinement_zones(zz)%rlevel) then 
                            tri_tree%node(ii)%flag = .true.
                            tri_tree%node(ii)%visited = .true.
                        elseif (tri_tree%node(ii)%level .GE. options%refinement_zones(zz)%rlevel) then 
                            refzone_limited = .true.
                        end if 
                    end if 
                end do 
            end if  

            !check if this node overlaps the geometry 
            if ((.NOT. refzone_limited) .AND. (tri_tree%node(ii)%level .LE. options%nrefine)) then 
                call kdtree%nodes_overlapping_region_2d(xmin,xmax,ymin,ymax,ntgtnode,tgt_nodes)
                if (ntgtnode .GT. 0) then 
                    do nn=1,ntgtnode
                        do ff=1,tgt_nodes(nn)%nitem
                            edgec => tgt_nodes(nn)%edge(ff)%edge
                            if (line_triangle_intersect_bool(edgec%origin%coordinate(1:2),edgec%opposite%origin%coordinate(1:2),&
                                tri_tree%vertices(tri_tree%node(ii)%vertex1,1:2),&
                                tri_tree%vertices(tri_tree%node(ii)%vertex2,1:2),&
                                tri_tree%vertices(tri_tree%node(ii)%vertex3,1:2))) then !intersect check
                                tri_tree%node(ii)%flag = .true.
                                tri_tree%node(ii)%visited = .true.
                                exit
                            end if 

                            ! edgec => tgt_nodes(nn)%edge(ff)%edge
                            ! xming = min(edgec%origin%coordinate(1),edgec%opposite%origin%coordinate(1))
                            ! xmaxg = max(edgec%origin%coordinate(1),edgec%opposite%origin%coordinate(1))
                            ! yming = min(edgec%origin%coordinate(2),edgec%opposite%origin%coordinate(2))
                            ! ymaxg = max(edgec%origin%coordinate(2),edgec%opposite%origin%coordinate(2))
                            ! if ((xmaxg .GE. xmin) .AND. (xming .LE. xmax)) then !overlap check 
                            !     if ((ymaxg .GE. ymin) .AND. (yming .LE. ymax)) then 
                            !         tri_tree%node(ii)%flag = .true.
                            !         tri_tree%node(ii)%visited = .true.
                            !         exit
                            !     end if 
                            ! end if 
                        end do 
                        if (tri_tree%node(ii)%flag) then 
                            exit 
                        end if 
                    end do 
                end if 
            end if 
            !flag nodes to refine ======
        end if 
    end do 

    !flood refinement
    reffloodidx = min(options%nrefine,rr)
    do ff=1,nref_flood(reffloodidx)
        do ii=1,tri_tree%nnode 
            if ((tri_tree%node(ii)%flag) .AND. (tri_tree%node(ii)%visited)) then 
                if (associated(tri_tree%node(ii)%adjacent1)) then 
                    if (.NOT. associated(tri_tree%node(ii)%adjacent1%child1)) then 
                        tri_tree%node(ii)%adjacent1%flag = .true.
                    end if 
                end if 
                if (associated(tri_tree%node(ii)%adjacent2)) then 
                    if (.NOT. associated(tri_tree%node(ii)%adjacent2%child1)) then 
                        tri_tree%node(ii)%adjacent2%flag = .true.
                    end if 
                end if 
                if (associated(tri_tree%node(ii)%adjacent3)) then 
                    if (.NOT. associated(tri_tree%node(ii)%adjacent3%child1)) then 
                        tri_tree%node(ii)%adjacent3%flag = .true.
                    end if 
                end if 
            end if 
        end do 
        do ii=1,tri_tree%nnode 
            if (tri_tree%node(ii)%flag) then 
                tri_tree%node(ii)%visited = .true.
            end if 
        end do 
    end do 
    
    !clean refinement 
    do ff=1,10
        nupdate = 0 
        do ii=1,tri_tree%nnode 
            if (.NOT. associated(tri_tree%node(ii)%child1)) then 
                if (.NOT. tri_tree%node(ii)%flag) then 
                    nadjf = 0
                    if (associated(tri_tree%node(ii)%adjacent1)) then
                        if (tri_tree%node(ii)%adjacent1%flag) then 
                            nadjf = nadjf + 1
                        end if 
                    end if 
                    if (associated(tri_tree%node(ii)%adjacent2)) then
                        if (tri_tree%node(ii)%adjacent2%flag) then 
                            nadjf = nadjf + 1
                        end if 
                    end if 
                    if (associated(tri_tree%node(ii)%adjacent3)) then
                        if (tri_tree%node(ii)%adjacent3%flag) then 
                            nadjf = nadjf + 1
                        end if 
                    end if 
                    if (nadjf .GE. 2) then 
                        tri_tree%node(ii)%flag = .true.
                        nupdate = nupdate + 1
                    end if 
                end if 
            end if 
        end do 
        if (nupdate == 0) then 
            exit
        end if
    end do 

    !refine nodes
    nrefine = 0 
    nnodeC = tri_tree%nnode 
    do ii=1,nnodeC
        if (tri_tree%node(ii)%flag) then 
            call tri_tree%refine_node(tri_tree%node(ii))
            nrefine = nrefine + 1
        end if 
    end do 

    !display
    if (options%cdisplay) then 
        write(*,'(A,I0,A,I0,A,I0,A,I0)') '    iteration: ',rr,' - ncell: ',tri_tree%nnode,' - nflood: ',nref_flood(reffloodidx),' - nrefine: ',nrefine
    end if

    !exit at zero refinements 
    if (nrefine == 0) then 
        exit
    end if 
end do 
return
end subroutine build_tritree


!refine trinode =========================
subroutine refine_node_tri(tri_tree,tnode)
implicit none 

!variables - import 
class(tritree), target :: tri_tree
class(trinode), target :: tnode

!variables - local 
type(trinode), pointer :: child1,child2,child3,child4 
type(trinode) :: node_map1,node_map2,node_map3

!build adjacent node maps
if (associated(tnode%adjacent1)) then 
    call tnode%map_adjacent_node(tnode%adjacent1,node_map1)
end if 
if (associated(tnode%adjacent2)) then 
    call tnode%map_adjacent_node(tnode%adjacent2,node_map2)
end if 
if (associated(tnode%adjacent3)) then 
    call tnode%map_adjacent_node(tnode%adjacent3,node_map3)
end if 

!build new edge vertices 
if (tnode%evertex1 == 0) then 
    tri_tree%nvertex = tri_tree%nvertex + 1
    tri_tree%vertices(tri_tree%nvertex,:) = 0.5d0*(tri_tree%vertices(tnode%vertex1,:) + tri_tree%vertices(tnode%vertex2,:))
    tnode%evertex1 = tri_tree%nvertex
    if (associated(node_map1%aevertex3)) then 
        node_map1%aevertex3 = tri_tree%nvertex
    end if 
end if 
if (tnode%evertex2 == 0) then 
    tri_tree%nvertex = tri_tree%nvertex + 1
    tri_tree%vertices(tri_tree%nvertex,:) = 0.5d0*(tri_tree%vertices(tnode%vertex2,:) + tri_tree%vertices(tnode%vertex3,:))
    tnode%evertex2 = tri_tree%nvertex
    if (associated(node_map2%aevertex3)) then 
        node_map2%aevertex3 = tri_tree%nvertex
    end if 
end if 
if (tnode%evertex3 == 0) then 
    tri_tree%nvertex = tri_tree%nvertex + 1
    tri_tree%vertices(tri_tree%nvertex,:) = 0.5d0*(tri_tree%vertices(tnode%vertex3,:) + tri_tree%vertices(tnode%vertex1,:))
    tnode%evertex3 = tri_tree%nvertex
    if (associated(node_map3%aevertex3)) then 
        node_map3%aevertex3 = tri_tree%nvertex
    end if 
end if 

!build child nodes
tri_tree%nnode = tri_tree%nnode + 1
child1 => tri_tree%node(tri_tree%nnode)
tri_tree%nnode = tri_tree%nnode + 1
child2 => tri_tree%node(tri_tree%nnode)
tri_tree%nnode = tri_tree%nnode + 1
child3 => tri_tree%node(tri_tree%nnode)
tri_tree%nnode = tri_tree%nnode + 1
child4 => tri_tree%node(tri_tree%nnode)

!set child 1
tnode%child1 => child1
child1%level = tnode%level + 1
child1%parent => tnode  
child1%vertex1 = tnode%vertex1
child1%vertex2 = tnode%evertex1
child1%vertex3 = tnode%evertex3
child1%adjacent2 => child4 
if (associated(node_map1%child1)) then 
    child1%adjacent1 => node_map1%child1
    call node_map1%child1%cascade_adjacency(child1%parent,child1)
else
    child1%adjacent1 => tnode%adjacent1
end if 
if (associated(node_map3%child3)) then 
    child1%adjacent3 => node_map3%child3
    call node_map3%child3%cascade_adjacency(child1%parent,child1)
else
    child1%adjacent3 => tnode%adjacent3
end if 

!set child 2
tnode%child2 => child2
child2%level = tnode%level + 1
child2%parent => tnode  
child2%vertex1 = tnode%evertex1
child2%vertex2 = tnode%vertex2
child2%vertex3 = tnode%evertex2
child2%adjacent3 => child4
if (associated(node_map1%child3)) then 
    child2%adjacent1 => node_map1%child3
    call node_map1%child3%cascade_adjacency(child2%parent,child2)
else
    child2%adjacent1 => tnode%adjacent1
end if 
if (associated(node_map2%child1)) then 
    child2%adjacent2 => node_map2%child1
    call node_map2%child1%cascade_adjacency(child2%parent,child2)
else
    child2%adjacent2 => tnode%adjacent2
end if 

!set child 3
tnode%child3 => child3
child3%level = tnode%level + 1
child3%parent => tnode  
child3%vertex1 = tnode%evertex3
child3%vertex2 = tnode%evertex2
child3%vertex3 = tnode%vertex3
child3%adjacent1 => child4
if (associated(node_map2%child3)) then 
    child3%adjacent2 => node_map2%child3
    call node_map2%child3%cascade_adjacency(child3%parent,child3)
else
    child3%adjacent2 => tnode%adjacent2
end if 
if (associated(node_map3%child1)) then 
    child3%adjacent3 => node_map3%child1
    call node_map3%child1%cascade_adjacency(child3%parent,child3)
else
    child3%adjacent3 => tnode%adjacent3
end if 

!set child 4
tnode%child4 => child4
child4%level = tnode%level + 1
child4%parent => tnode  
child4%vertex1 = tnode%evertex2
child4%vertex2 = tnode%evertex3
child4%vertex3 = tnode%evertex1
child4%adjacent1 => child3
child4%adjacent2 => child1
child4%adjacent3 => child2
return 
end subroutine refine_node_tri


!build oriented adjacent node map =========================
subroutine map_adjacent_node_tri(self,adjnode,node_map)
implicit none 

!variables - import 
class(trinode), target :: self
type(trinode), pointer :: adjnode
type(trinode) :: node_map

!find relative orientation to build map
if (associated(adjnode)) then 
    if (associated(adjnode%adjacent1,self)) then 
        node_map%child1 => adjnode%child2
        node_map%child2 => adjnode%child3
        node_map%child3 => adjnode%child1
        node_map%aevertex1 => adjnode%evertex2
        node_map%aevertex2 => adjnode%evertex3
        node_map%aevertex3 => adjnode%evertex1
    elseif (associated(adjnode%adjacent2,self)) then 
        node_map%child1 => adjnode%child3
        node_map%child2 => adjnode%child1
        node_map%child3 => adjnode%child2
        node_map%aevertex1 => adjnode%evertex3
        node_map%aevertex2 => adjnode%evertex1
        node_map%aevertex3 => adjnode%evertex2
    elseif (associated(adjnode%adjacent3,self)) then 
        node_map%child1 => adjnode%child1
        node_map%child2 => adjnode%child2
        node_map%child3 => adjnode%child3
        node_map%aevertex1 => adjnode%evertex1
        node_map%aevertex2 => adjnode%evertex2
        node_map%aevertex3 => adjnode%evertex3
    end if 
else
    node_map%child1 => null()
    node_map%child2 => null()
    node_map%child3 => null()
    node_map%aevertex1 => null()
    node_map%aevertex2 => null()
    node_map%aevertex3 => null()
end if 
return 
end subroutine map_adjacent_node_tri


!cascade adjacency update =========================
recursive subroutine cascade_adjacency_tri(self,adj_tgt,adj_new)
implicit none

!variables - import 
class(trinode) :: self
type(trinode), pointer :: adj_tgt,adj_new

!variables - local 
logical :: updated

!update if present 
updated = .false.
if (associated(self%adjacent1,adj_tgt)) then 
    self%adjacent1 => adj_new
end if 
if (associated(self%adjacent2,adj_tgt)) then 
    self%adjacent2 => adj_new
end if 
if (associated(self%adjacent3,adj_tgt)) then 
    self%adjacent3 => adj_new
end if 

!call on children if an update is made
if (updated) then 
    if (associated(self%child1)) then 
        call self%child1%cascade_adjacency(adj_tgt,adj_new)
    end if 
    if (associated(self%child2)) then 
        call self%child2%cascade_adjacency(adj_tgt,adj_new)
    end if 
    if (associated(self%child3)) then 
        call self%child3%cascade_adjacency(adj_tgt,adj_new)
    end if 
end if 
return 
end subroutine cascade_adjacency_tri   


!initialise tree tri =========================
subroutine initialise_tritree_tri(tri_tree,options)
implicit none 

!variables - inout 
type(tritree), target :: tri_tree 
type(hex_options) :: options 


!allocate tree
tri_tree%nnode = 1
tri_tree%nvertex = 3
allocate(tri_tree%node(options%ncell_max))
allocate(tri_tree%vertices(2*options%ncell_max,3))
allocate(tri_tree%vertex_tag(2*options%ncell_max))

!build initial domain
tri_tree%vertices(1,:) = (/-options%farfield_r,-options%farfield_r,0.0d0/)
tri_tree%vertices(2,:) = (/0.0d0,options%farfield_r,0.0d0/)
tri_tree%vertices(3,:) = (/options%farfield_r,-options%farfield_r,0.0d0/)
tri_tree%node(1)%index = 1
tri_tree%node(1)%adjacent1 => null()
tri_tree%node(1)%adjacent2 => null()
tri_tree%node(1)%adjacent3 => null()
tri_tree%node(1)%vertex1 = 1
tri_tree%node(1)%vertex2 = 2
tri_tree%node(1)%vertex3 = 3
return 
end subroutine initialise_tritree_tri


!initialise tree quad =========================
subroutine initialise_tritree_quad(tri_tree,options)
implicit none 

!variables - inout 
type(tritree), target :: tri_tree 
type(hex_options) :: options 


!allocate tree
tri_tree%nnode = 4
tri_tree%nvertex = 5
allocate(tri_tree%node(options%ncell_max))
allocate(tri_tree%vertices(2*options%ncell_max,3))
allocate(tri_tree%vertex_tag(2*options%ncell_max))

!build initial domain
tri_tree%vertices(1,:) = (/-options%farfield_r,-options%farfield_r,0.0d0/)
tri_tree%vertices(2,:) = (/-options%farfield_r,options%farfield_r,0.0d0/)
tri_tree%vertices(3,:) = (/options%farfield_r,options%farfield_r,0.0d0/)
tri_tree%vertices(4,:) = (/options%farfield_r,-options%farfield_r,0.0d0/)
tri_tree%vertices(5,:) = (/0.0d0,0.0d0,0.0d0/)
tri_tree%node(1)%index = 1
tri_tree%node(1)%adjacent1 => tri_tree%node(2)
tri_tree%node(1)%adjacent2 => tri_tree%node(4)
tri_tree%node(1)%adjacent3 => null()
tri_tree%node(1)%vertex1 = 2
tri_tree%node(1)%vertex2 = 5
tri_tree%node(1)%vertex3 = 1
tri_tree%node(2)%index = 2
tri_tree%node(2)%adjacent1 => tri_tree%node(3)
tri_tree%node(2)%adjacent2 => tri_tree%node(1)
tri_tree%node(2)%adjacent3 => null()
tri_tree%node(2)%vertex1 = 3
tri_tree%node(2)%vertex2 = 5
tri_tree%node(2)%vertex3 = 2
tri_tree%node(3)%index = 3
tri_tree%node(3)%adjacent1 => tri_tree%node(4)
tri_tree%node(3)%adjacent2 => tri_tree%node(2)
tri_tree%node(3)%adjacent3 => null()
tri_tree%node(3)%vertex1 = 4
tri_tree%node(3)%vertex2 = 5
tri_tree%node(3)%vertex3 = 3
tri_tree%node(4)%index = 4
tri_tree%node(4)%adjacent1 => tri_tree%node(1)
tri_tree%node(4)%adjacent2 => tri_tree%node(3)
tri_tree%node(4)%adjacent3 => null()
tri_tree%node(4)%vertex1 = 1
tri_tree%node(4)%vertex2 = 5
tri_tree%node(4)%vertex3 = 4
return 
end subroutine initialise_tritree_quad


!initialise tree hexagon =========================
subroutine initialise_tritree_hexagon(tri_tree,options)
implicit none 

!variables - inout 
type(tritree), target :: tri_tree 
type(hex_options) :: options 

!allocate tree
tri_tree%nnode = 6
tri_tree%nvertex = 7
allocate(tri_tree%node(options%ncell_max))
allocate(tri_tree%vertices(2*options%ncell_max,3))
allocate(tri_tree%vertex_tag(2*options%ncell_max))

!build initial domain
tri_tree%vertices(1,:) = (/-options%farfield_r,0.0d0,0.0d0/)
tri_tree%vertices(2,:) = (/-0.5d0*options%farfield_r,(sqrt(3.0)/2.0)*options%farfield_r,0.0d0/)
tri_tree%vertices(3,:) = (/0.5d0*options%farfield_r,(sqrt(3.0)/2.0)*options%farfield_r,0.0d0/)
tri_tree%vertices(4,:) = (/options%farfield_r,0.0d0,0.0d0/)
tri_tree%vertices(5,:) = (/0.5d0*options%farfield_r,-(sqrt(3.0)/2.0)*options%farfield_r,0.0d0/)
tri_tree%vertices(6,:) = (/-0.5d0*options%farfield_r,-(sqrt(3.0)/2.0)*options%farfield_r,0.0d0/)
tri_tree%vertices(7,:) = (/0.0d0,0.0d0,0.0d0/)
tri_tree%node(1)%index = 1
tri_tree%node(1)%adjacent1 => tri_tree%node(2)
tri_tree%node(1)%adjacent2 => tri_tree%node(6)
tri_tree%node(1)%adjacent3 => null()
tri_tree%node(1)%vertex1 = 2
tri_tree%node(1)%vertex2 = 7
tri_tree%node(1)%vertex3 = 1
tri_tree%node(2)%index = 2
tri_tree%node(2)%adjacent1 => tri_tree%node(3)
tri_tree%node(2)%adjacent2 => tri_tree%node(1)
tri_tree%node(2)%adjacent3 => null()
tri_tree%node(2)%vertex1 = 3
tri_tree%node(2)%vertex2 = 7
tri_tree%node(2)%vertex3 = 2
tri_tree%node(3)%index = 3
tri_tree%node(3)%adjacent1 => tri_tree%node(4)
tri_tree%node(3)%adjacent2 => tri_tree%node(2)
tri_tree%node(3)%adjacent3 => null()
tri_tree%node(3)%vertex1 = 4
tri_tree%node(3)%vertex2 = 7
tri_tree%node(3)%vertex3 = 3
tri_tree%node(4)%index = 4
tri_tree%node(4)%adjacent1 => tri_tree%node(5)
tri_tree%node(4)%adjacent2 => tri_tree%node(3)
tri_tree%node(4)%adjacent3 => null()
tri_tree%node(4)%vertex1 = 5
tri_tree%node(4)%vertex2 = 7
tri_tree%node(4)%vertex3 = 4
tri_tree%node(5)%index = 5
tri_tree%node(5)%adjacent1 => tri_tree%node(6)
tri_tree%node(5)%adjacent2 => tri_tree%node(4)
tri_tree%node(5)%adjacent3 => null()
tri_tree%node(5)%vertex1 = 6
tri_tree%node(5)%vertex2 = 7
tri_tree%node(5)%vertex3 = 5
tri_tree%node(6)%index = 6
tri_tree%node(6)%adjacent1 => tri_tree%node(1)
tri_tree%node(6)%adjacent2 => tri_tree%node(5)
tri_tree%node(6)%adjacent3 => null()
tri_tree%node(6)%vertex1 = 1
tri_tree%node(6)%vertex2 = 7
tri_tree%node(6)%vertex3 = 6
return 
end subroutine initialise_tritree_hexagon

end module hex_tritree