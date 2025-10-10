!hex mesh generator 2d main program
!max wood
!version : 0.0.1
!updated : 20-12-24

!module 
module hex2d
use hex_data
use hex_geometry
use halfedge_mesh
use geometry_kdtree

!routines 
contains 

!build mesh 2d =========================
function hex2d_mesh(geometry,options) result(mesh)
implicit none 

!variables - inout 
type(hex_mesh) mesh
type(halfedge) :: geometry
type(hex_options) :: options 

!variables - local 
integer(in64) :: ii
type(hex_mesh), target :: mesh_full
type(gkdtree), target :: kdtree
type(tritree), target :: tri_tree 

!set mesh valididy flag
mesh%isvalid = .true.

!extend the geometry arrays to accomodate new entries from splitting edges 
call geometry%extend()
call geometry%extend()
call geometry%extend()



!build kdtree on the geometry
call kdtree%build_tree(geometry,10_in64,4_in64,.true.)

!build refined tritree
call build_tritree(options,tri_tree,kdtree)

!construct the full mesh 
call build_full_mesh(mesh_full,tri_tree)

!extend the mesh arrays to accomodate new items from clipping to the geometry 
call mesh_full%extend()

!clip the full mesh to the geometry
call clip_to_geometry(mesh,mesh_full,geometry,kdtree)

!postprocess and clean the clipped mesh 








print *, '## debug'

print *, 'nvertex = ',tri_tree%nvertex 
open(11,file='vertices')
do ii=1,tri_tree%nvertex 
    write(11,'(F0.16,A,F0.16,A,F0.16)') tri_tree%vertices(ii,1),' ',tri_tree%vertices(ii,2),' ',tri_tree%vertices(ii,3)
end do 
close(11)

open(11,file='faces')
do ii=1,tri_tree%nnode 
    if (.NOT. associated(tri_tree%node(ii)%child1)) then 
        write(11,'(I0,A,I0,A,I0)') tri_tree%node(ii)%vertex1,' ',tri_tree%node(ii)%vertex2,' ',tri_tree%node(ii)%vertex3
    end if 
end do 
close(11)



open(11,file='verticesf')
do ii=1,mesh_full%nvertex
    write(11,'(F0.16,A,F0.16,A,F0.16)') mesh_full%vertex(ii)%coordinate(1),' ',&
    mesh_full%vertex(ii)%coordinate(2),' ',mesh_full%vertex(ii)%coordinate(3)
end do 
close(11)

open(11,file='facesf')
do ii=1,mesh_full%nedge
    write(11,'(I0,A,I0)') mesh_full%edge(ii)%vertex1%index,' ',mesh_full%edge(ii)%vertex2%index
end do 
close(11)


print *, 'complete'


return 
end function hex2d_mesh


!clip to geometry =========================
subroutine clip_to_geometry(mesh,mesh_full,geometry,kdtree)
implicit none 

!variables - import
type(hex_mesh) :: mesh
type(hex_mesh), target :: mesh_full
type(halfedge), target :: geometry
type(gkdtree), target :: kdtree

!variables - local 
integer(in64) :: ii,jj,nn,ee,nfront,nfrontN
integer(in64) :: nedge_mesh,ntgtnode,nedgecheck,n_split_vertex,nvidx,medge_tgt,vnew,etgt,vtgt
integer(in64) :: neidx(2),edge_split_vertex(mesh_full%nvertex)
integer(in64), dimension(:), allocatable :: efront,efrontN
integer(in64), dimension(:,:), allocatable :: mesh_v2e
real(dp) :: xmax,xmin,ymax,ymin
real(dp) :: vi(3),edge_split_vertex_distance(mesh_full%nvertex)
type(node) :: tgt_nodes(kdtree%nnode)
type(edge) :: edge_check(geometry%nedge)

!set geometry flags
do ii=1,geometry%nedge
    geometry%edge(ii)%flag = .false.
end do 

!set mesh vertex tags to zero (a non-zero tag indicates the vertex is shared with one on the clipped geometry)
do ii=1,mesh_full%nvertex
    mesh_full%vertex(ii)%tag = 0 
end do 

!clip mesh edges to the geometry 
ntgtnode = 0 
edge_split_vertex(:) = 0 
edge_split_vertex_distance(:) = 0.0d0 
nedge_mesh = mesh_full%nedge
do ii=1,nedge_mesh

    !check if this edge overlaps the geometry 
    xmin = min(mesh_full%edge(ii)%vertex1%coordinate(1),mesh_full%edge(ii)%vertex2%coordinate(1)) 
    xmax = max(mesh_full%edge(ii)%vertex1%coordinate(1),mesh_full%edge(ii)%vertex2%coordinate(1)) 
    ymin = min(mesh_full%edge(ii)%vertex1%coordinate(2),mesh_full%edge(ii)%vertex2%coordinate(2)) 
    ymax = max(mesh_full%edge(ii)%vertex1%coordinate(2),mesh_full%edge(ii)%vertex2%coordinate(2)) 
    call kdtree%nodes_overlapping_region_2d(xmin,xmax,ymin,ymax,ntgtnode,tgt_nodes)
    if (ntgtnode == 0) then 
        cycle
    end if 

    !build a list of geometry edges to check 
    nedgecheck = 0 
    do nn=1,ntgtnode
        do jj=1,tgt_nodes(nn)%nitem
            nedgecheck = nedgecheck + 1
            edge_check(nedgecheck) = tgt_nodes(nn)%edge(jj)%edge
            geometry%edge(edge_check(nedgecheck)%index)%flag = .false.
            geometry%edge(edge_check(nedgecheck)%index)%opposite%flag = .false.
        end do 
    end do 

    !split each geometry edge the mesh edge interects and store the new geometry vertices 
    !(update the kdtree as each edge is split)
    n_split_vertex = 0 
    do ee=1,nedgecheck
        if ((.NOT. edge_check(ee)%flag) .AND. (.NOT. edge_check(ee)%opposite%flag)) then 
            vi(1:2) = line_line_intersect(edge_check(ee)%origin%coordinate(1:2),edge_check(ee)%opposite%origin%coordinate(1:2),&
            mesh_full%edge(ii)%vertex1%coordinate(1:2),mesh_full%edge(ii)%vertex2%coordinate(1:2))
            if (is_in_line_segment(edge_check(ee)%origin%coordinate(1:2),edge_check(ee)%opposite%origin%coordinate(1:2),vi) .AND. &
                is_in_line_segment(mesh_full%edge(ii)%vertex1%coordinate(1:2),mesh_full%edge(ii)%vertex2%coordinate(1:2),vi)) then 

                !split the geometry edge 
                call geometry%edge(edge_check(ee)%index)%split(geometry,vi,nvidx,neidx)
                geometry%edge(neidx(1))%flag = .true.
                geometry%edge(neidx(2))%flag = .true.

                !add the new edges to the geometry kdtree
                geometry%edge(neidx(1))%tag = geometry%edge(edge_check(ee)%index)%tag
                geometry%edge(neidx(2))%tag = geometry%edge(edge_check(ee)%index)%tag
                call kdtree%node(geometry%edge(edge_check(ee)%index)%tag)%add_item_2d(geometry%edge(neidx(1)))
                call kdtree%node(geometry%edge(edge_check(ee)%index)%tag)%add_item_2d(geometry%edge(neidx(2)))

                !flag this edge and its opposite
                geometry%edge(edge_check(ee)%index)%flag = .true.
                geometry%edge(edge_check(ee)%index)%opposite%flag = .true.

                !store the new geometry vertex to use to split the mesh edge 
                n_split_vertex = n_split_vertex + 1
                edge_split_vertex(n_split_vertex) = nvidx
            end if 
        end if 
    end do 

    !cycle if no proper intersections 
    if (n_split_vertex == 0) then 
        cycle
    end if 

    !order edge_split_vertex by distance from the origin of the mesh edge 
    if (n_split_vertex .GT. 1) then 
        do ee=1,n_split_vertex
            edge_split_vertex_distance(ee) = norm2(geometry%vertex(edge_split_vertex(ee))%coordinate(1:2) - &
                                                   mesh_full%edge(ii)%vertex1%coordinate(1:2))
        end do 
        call order_items_int_by_real(edge_split_vertex(1:n_split_vertex),edge_split_vertex_distance(1:n_split_vertex))
    end if 
     
    !split the mesh edge into n segments using each unique intersection with the geometry
    medge_tgt = ii 
    do ee=1,n_split_vertex
        call mesh_full%edge(medge_tgt)%split(mesh_full,geometry%vertex(edge_split_vertex(ee))%coordinate,medge_tgt,vnew)
        mesh_full%vertex(vnew)%tag = edge_split_vertex(ee)
    end do 
end do 

!build mesh v2e
allocate(mesh_v2e(mesh_full%nvertex,6))
mesh_v2e(:,:) = 0
do ee=1,mesh_full%nedge
    do ii=1,6
        if (mesh_v2e(mesh_full%edge(ee)%vertex1%index,ii) == 0) then 
            mesh_v2e(mesh_full%edge(ee)%vertex1%index,ii) = mesh_full%edge(ee)%index
            exit 
        end if 
    end do 
    do ii=1,6
        if (mesh_v2e(mesh_full%edge(ee)%vertex2%index,ii) == 0) then 
            mesh_v2e(mesh_full%edge(ee)%vertex2%index,ii) = mesh_full%edge(ee)%index
            exit
        end if 
    end do 
end do 

!flood the geometry external portion of the mesh (not flooding through tagged vertices)
nfront = 0
nfrontN = 0
allocate(efront(mesh_full%nedge))
allocate(efrontN(mesh_full%nedge))
do ee=1,mesh_full%nedge
    mesh_full%edge(ee)%external = .false.
end do 
do ee=1,mesh_full%nvertex
    mesh_full%vertex(ee)%external = .false.
end do 
do ee=1,mesh_full%nedge
    if ((mesh_full%edge(ee)%cell1 == -2) .OR. (mesh_full%edge(ee)%cell2 == -2)) then 
        nfront = 1
        efront(nfront) = ee
    end if 
end do 
do nn=1,mesh_full%nedge
    nfrontN = 0 
    do ee=1,nfront

        !edge
        etgt = efront(ee)

        !tag
        mesh_full%edge(etgt)%external = .true.
        mesh_full%vertex(mesh_full%edge(etgt)%vertex1%index)%external = .true.
        mesh_full%vertex(mesh_full%edge(etgt)%vertex2%index)%external = .true.

        !flood across vertex 1
        vtgt = mesh_full%edge(etgt)%vertex1%index
        if (mesh_full%vertex(vtgt)%tag == 0) then 
            do ii=1,6
                if (mesh_v2e(vtgt,ii) .NE. 0) then 
                    if (.NOT. mesh_full%edge(mesh_v2e(vtgt,ii))%external) then 
                        nfrontN = nfrontN + 1
                        efrontN(nfrontN) = mesh_v2e(vtgt,ii)
                        mesh_full%edge(mesh_v2e(vtgt,ii))%external = .true.
                    end if 
                end if 
            end do 
        end if 

        !flood across vertex 2
        vtgt = mesh_full%edge(etgt)%vertex2%index
        if (mesh_full%vertex(vtgt)%tag == 0) then 
            do ii=1,6
                if (mesh_v2e(vtgt,ii) .NE. 0) then 
                    if (.NOT. mesh_full%edge(mesh_v2e(vtgt,ii))%external) then 
                        nfrontN = nfrontN + 1
                        efrontN(nfrontN) = mesh_v2e(vtgt,ii)
                        mesh_full%edge(mesh_v2e(vtgt,ii))%external = .true.
                    end if 
                end if 
            end do 
        end if 
    end do 

    ! print *, nfrontN
    ! print *, efront(1:nfront)

    if (nfrontN == 0) then 
        exit 
    end if 
    nfront = nfrontN
    efront(1:nfront) = efrontN(1:nfront) 
end do 



!extract the desired (internal or external) portion of the mesh 



!build the complete final mesh by unifying the clipped mesh with the geometry mesh



return 
end subroutine clip_to_geometry


!build full mesh =========================
subroutine build_full_mesh(mesh_full,tri_tree)
implicit none 

!variables - import
type(hex_mesh), target :: mesh_full
type(tritree), target :: tri_tree 

!variables - local 
integer(in64) :: ii,eins,vins 
integer(in64) :: v2v(tri_tree%nnode,6)
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

                    ! edge_temp(eins)%cell1 = tri_tree%node(ii)%vertex1
                    ! edge_temp(eins)%cell2 = tri_tree%node(ii)%vertex2


                    call add_connection(v2v(tri_tree%node(ii)%index,:),v2v(tri_tree%node(ii)%adjacent1%index,:),&
                                        tri_tree%node(ii)%index,tri_tree%node(ii)%adjacent1%index) 
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

                    ! edge_temp(eins)%cell1 = tri_tree%node(ii)%vertex2
                    ! edge_temp(eins)%cell2 = tri_tree%node(ii)%vertex3


                    call add_connection(v2v(tri_tree%node(ii)%index,:),v2v(tri_tree%node(ii)%adjacent2%index,:),&
                                        tri_tree%node(ii)%index,tri_tree%node(ii)%adjacent2%index) 
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

                    ! edge_temp(eins)%cell1 = tri_tree%node(ii)%vertex3
                    ! edge_temp(eins)%cell2 = tri_tree%node(ii)%vertex1


                    call add_connection(v2v(tri_tree%node(ii)%index,:),v2v(tri_tree%node(ii)%adjacent3%index,:),&
                                        tri_tree%node(ii)%index,tri_tree%node(ii)%adjacent3%index) 
                end if
            end if 
        end if 
    end if 
end do 
mesh_full%nedge = eins

!trim edges 
allocate(mesh_full%edge(eins))
mesh_full%edge(:) = edge_temp(1:eins)


! !build and index edges 
! v2v(:,:) = 0
! allocate(mesh_full%edge(tri_tree%nnode*6))
! eins = 0 
! do ii=1,tri_tree%nnode 
!     if (.NOT. associated(tri_tree%node(ii)%child1)) then 
!         if (associated(tri_tree%node(ii)%adjacent1)) then 
!             if (.NOT. associated(tri_tree%node(ii)%adjacent1%child1)) then 
!                 if (.NOT. is_connected_v2v(v2v(tri_tree%node(ii)%index,:),tri_tree%node(ii)%adjacent1%index)) then 
!                     eins = eins + 1
!                     mesh_full%edge(eins)%vertex1 => mesh_full%vertex(tri_tree%node(ii)%index)
!                     mesh_full%edge(eins)%vertex2 => mesh_full%vertex(tri_tree%node(ii)%adjacent1%index)
!                     mesh_full%edge(eins)%cell1 = tri_tree%node(ii)%vertex1
!                     mesh_full%edge(eins)%cell2 = tri_tree%node(ii)%vertex2
!                     ! mesh_full%edge(eins,1) = tri_tree%node(ii)%index
!                     ! mesh_full%edge(eins,2) = tri_tree%node(ii)%adjacent1%index
!                     ! mesh_full%edge(eins,3) = tri_tree%node(ii)%vertex1
!                     ! mesh_full%edge(eins,4) = tri_tree%node(ii)%vertex2
!                     call add_connection(v2v(tri_tree%node(ii)%index,:),v2v(tri_tree%node(ii)%adjacent1%index,:),&
!                                         tri_tree%node(ii)%index,tri_tree%node(ii)%adjacent1%index) 
!                 end if 
!             end if 
!         end if 
!         if (associated(tri_tree%node(ii)%adjacent2)) then 
!             if (.NOT. associated(tri_tree%node(ii)%adjacent2%child1)) then 
!                 if (.NOT. is_connected_v2v(v2v(tri_tree%node(ii)%index,:),tri_tree%node(ii)%adjacent2%index)) then 
!                     eins = eins + 1
!                     mesh_full%edge(eins)%vertex1 => mesh_full%vertex(tri_tree%node(ii)%index)
!                     mesh_full%edge(eins)%vertex2 => mesh_full%vertex(tri_tree%node(ii)%adjacent2%index)
!                     mesh_full%edge(eins)%cell1 = tri_tree%node(ii)%vertex2
!                     mesh_full%edge(eins)%cell2 = tri_tree%node(ii)%vertex3
!                     ! mesh_full%edge(eins,1) = tri_tree%node(ii)%index
!                     ! mesh_full%edge(eins,2) = tri_tree%node(ii)%adjacent2%index
!                     ! mesh_full%edge(eins,3) = tri_tree%node(ii)%vertex2
!                     ! mesh_full%edge(eins,4) = tri_tree%node(ii)%vertex3
!                     call add_connection(v2v(tri_tree%node(ii)%index,:),v2v(tri_tree%node(ii)%adjacent2%index,:),&
!                                         tri_tree%node(ii)%index,tri_tree%node(ii)%adjacent2%index) 
!                 end if
!             end if 
!         end if 
!         if (associated(tri_tree%node(ii)%adjacent3)) then 
!             if (.NOT. associated(tri_tree%node(ii)%adjacent3%child1)) then 
!                 if (.NOT. is_connected_v2v(v2v(tri_tree%node(ii)%index,:),tri_tree%node(ii)%adjacent3%index)) then 
!                     eins = eins + 1
!                     mesh_full%edge(eins)%vertex1 => mesh_full%vertex(tri_tree%node(ii)%index)
!                     mesh_full%edge(eins)%vertex2 => mesh_full%vertex(tri_tree%node(ii)%adjacent3%index)
!                     mesh_full%edge(eins)%cell1 = tri_tree%node(ii)%vertex3
!                     mesh_full%edge(eins)%cell2 = tri_tree%node(ii)%vertex1
!                     ! mesh_full%edge(eins,1) = tri_tree%node(ii)%index
!                     ! mesh_full%edge(eins,2) = tri_tree%node(ii)%adjacent3%index
!                     ! mesh_full%edge(eins,3) = tri_tree%node(ii)%vertex3
!                     ! mesh_full%edge(eins,4) = tri_tree%node(ii)%vertex1
!                     call add_connection(v2v(tri_tree%node(ii)%index,:),v2v(tri_tree%node(ii)%adjacent3%index,:),&
!                                         tri_tree%node(ii)%index,tri_tree%node(ii)%adjacent3%index) 
!                 end if
!             end if 
!         end if 
!     end if 
! end do 
! mesh_full%nedge = eins


return 
end subroutine build_full_mesh


!order items =========================
subroutine order_items_int_by_real(array,order)
implicit none 

!variables - inout
integer(in64), dimension(:) :: array
real(dp), dimension(:) :: order

!variables - local 
integer(in64) :: ii 
integer(in64) :: item_ins,item_tgt
integer(in64) :: array_ordered(size(array,dim=1))
logical :: mask(size(array,dim=1))

!assign order 
item_ins = 1 
mask(:) = .true.
do ii=1,size(array,dim=1)
    item_tgt = minloc(order,dim=1,mask=mask)
    array_ordered(item_ins) = array(item_tgt)
    mask(item_tgt) = .false.
    item_ins = item_ins + 1
end do 
array = array_ordered
return 
end subroutine order_items_int_by_real


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
integer(in64) :: ii,rr,ff,nn
integer(in64) :: nnodeC,nadjf,ntgtnode
real(dp) :: xmax,xmin,ymax,ymin,xmaxg,xming,ymaxg,yming
type(edge), pointer :: edge
type(node) :: tgt_nodes(kdtree%nnode)

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

!initialise vertex tags
tri_tree%vertex_tag(:) = 0 

!refine tree 
do rr=1,options%nrefine 

    !reset flags
    do ii=1,tri_tree%nnode
        tri_tree%node(ii)%flag = .false.
        tri_tree%node(ii)%visited = .false.
    end do 

    !select nodes to refine
    do ii=1,tri_tree%nnode 
        if (tri_tree%node(ii)%level == rr) then 

            !cell bounding box
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

            !check if this node overlaps the geometry
            call kdtree%nodes_overlapping_region_2d(xmin,xmax,ymin,ymax,ntgtnode,tgt_nodes)
            if (ntgtnode .GT. 0) then 
                do nn=1,ntgtnode
                    do ff=1,tgt_nodes(nn)%nitem
                        edge => tgt_nodes(nn)%edge(ff)%edge
                        xming = min(edge%origin%coordinate(1),edge%opposite%origin%coordinate(1))
                        xmaxg = max(edge%origin%coordinate(1),edge%opposite%origin%coordinate(1))
                        yming = min(edge%origin%coordinate(2),edge%opposite%origin%coordinate(2))
                        ymaxg = max(edge%origin%coordinate(2),edge%opposite%origin%coordinate(2))
                        if ((xmaxg .GE. xmin) .AND. (xming .LE. xmax)) then 
                            if ((ymaxg .GE. ymin) .AND. (yming .LE. ymax)) then 
                                tri_tree%node(ii)%flag = .true.
                                tri_tree%node(ii)%visited = .true.
                                exit
                            end if 
                        end if 
                    end do 
                    if (tri_tree%node(ii)%flag) then 
                        exit 
                    end if 
                end do 
            end if 
        end if 
    end do 

    !flood refinement
    do ff=1,options%nflood 
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
    do ii=1,tri_tree%nnode 
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
            if (nadjf == 2) then 
                tri_tree%node(ii)%flag = .true.
            end if 
        end if 
    end do 

    !refine nodes
    nnodeC = tri_tree%nnode 
    do ii=1,nnodeC
        if (tri_tree%node(ii)%flag) then 
            call tri_tree%refine_node(tri_tree%node(ii))
        end if 
    end do 
end do 
return
end subroutine build_tritree


end module hex2d