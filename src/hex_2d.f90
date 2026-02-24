!hex mesh generator 2d main program
!max wood
!version : 0.1.2
!updated : 24-02-26

!module 
module hex2d
use hex_tritree
use hex_geometry
use halfedge_mesh
use geometry_kdtree
use hex_postprocess
use hex_data_methods

!routines 
contains 

!build mesh 2d =========================
function hex2d_mesh(geometry,options) result(mesh)
implicit none 

!variables - inout 
type(hex_mesh), target :: mesh
type(halfedge) :: geometry
type(hex_options) :: options 

!variables - local 
integer(in64) :: ii
integer(in64) :: nbisected,nshort,nsmall,ndegenerate,ndedge,nunassociated,ndoublebc,nperturb,ntagged,nupdate
type(hex_mesh), target :: mesh_full
type(gkdtree), target :: kdtree
type(tritree), target :: tri_tree 

!set mesh validity flag
mesh%isvalid = .true.

!initialise properties
mesh_full%ncell = 0
mesh_full%nedge = 0
mesh_full%nvertex = 0
mesh_full%nvertex_surfint = 0
mesh%ncell = 0
mesh%nedge = 0
mesh%nvertex = 0
mesh%nvertex_surfint = 0

!extend the geometry arrays to accomodate new entries from splitting edges 
call geometry%extend(options%ncell_max)

!build kdtree on the geometry
if (options%cdisplay) then 
    write(*,'(A)') '--> building geometry kdtree'
end if 
call kdtree%build_tree(geometry,20_in64,4_in64,options%cdisplay)

!get the tree offset to centre the tree on the geometry centroid
options%tree_offset = geometry%get_centroid()

!build refined tree
if (options%cdisplay) then 
    write(*,'(A)') '--> refining mesh tree of type: '//options%mesh_treetype
end if 
if (options%mesh_treetype == 'quadrilateral') then 

elseif (options%mesh_treetype == 'triangle') then 
    call build_tritree(options,tri_tree,kdtree)
end if 

!construct the full mesh 
if (options%cdisplay) then 
    write(*,'(A)') '--> building mesh: '//options%tree_mesh_relation
end if 
if (options%mesh_treetype == 'quadrilateral') then 
    if (options%tree_mesh_relation == 'primal') then 
  
    elseif (options%tree_mesh_relation == 'dual') then 

    end if 
elseif (options%mesh_treetype == 'triangle') then 
    if (options%tree_mesh_relation == 'primal') then 
        call build_full_mesh_from_tritree_primal(mesh_full,tri_tree)
    elseif (options%tree_mesh_relation == 'dual') then 
        call build_full_mesh_from_tritree_dual(mesh_full,tri_tree)
    end if 
end if 

!extend the mesh arrays to accomodate new items from clipping to the geometry 
call mesh_full%extend(options%ncell_max)

!perturb mesh_full vertices that lie on the geometry 
call perturb_ongeom_vertices(mesh_full,kdtree,options%vtx_proximity_tolerance,nperturb)
if (options%cdisplay) then 
    write(*,'(A,I0,A)') '    {perturbed ',nperturb,' geometry co-incident vertices}'
end if

!clip the full mesh to the geometry
if (options%cdisplay) then 
    write(*,'(A)') '--> clipping the mesh to the geometry'
end if 
call clip_to_geometry(mesh,mesh_full,geometry,kdtree,options)

!index the mesh 
call mesh%index_edges()
call mesh%index_cells()

!postprocess
if (options%allow_postprocess) then 

    !postprocess display
    if (options%cdisplay) then 
        write(*,'(A)') '--> post-processing the mesh'
    end if

    !remove double boundary condition surface edges
    !(these occur when the geometry extends beyond the mesh far field)
    call remove_double_bc_edges(mesh,ndoublebc)
    if (options%cdisplay) then 
        write(*,'(A,I0,A)') '    {removed ',ndoublebc,' double boundary condition edges}'
    end if

    !remove geometry surface edges that were not accociated with any cell
    !(this happens if the object they belong to is fully contained within a single mesh cell and does not intersect any mesh edges)
    call remove_unlinked_surface_edges(mesh,nunassociated)
    if (options%cdisplay) then 
        write(*,'(A,I0,A)') '    {removed ',nunassociated,' unassociated geometry surface edges}'
    end if

    !divide partitioned cells 
    nbisected = 0 
    do ii=1,mesh%nedge
        call partition_bisected_cells(mesh,nupdate)
        nbisected = nbisected + nupdate
        if (nupdate == 0) then 
            exit 
        end if 
    end do 
    if (options%cdisplay) then 
        write(*,'(A,I0,A)') '    {partitioned ',nbisected,' geometry bisected cells}'
    end if

    !clip to any specified planes 
    if (options%nclipplane .GT. 0) then 
        do ii=1,options%nclipplane
            call clip_to_plane(mesh,options%clip_planes(ii)%v1(1:2),options%clip_planes(ii)%v2(1:2))
        end do 
        if (options%cdisplay) then 
        write(*,'(A,I0,A)') '    {clipped mesh to ',options%nclipplane,' planes}'
    end if
    end if 

    !check for and eliminate zero length edges 
    nshort = 0 
    do ii=1,mesh%nedge
        call collapse_short_edges(mesh,options,nupdate)
        nshort = nshort + nupdate
        if (nupdate == 0) then 
            exit 
        end if 
    end do 
    if (options%cdisplay) then 
        write(*,'(A,I0,A)') '    {collapsed ',nshort,' short edges}'
    end if

    !smooth mesh 
    if ((options%nsmooth_interlayer .GT. 0) .OR. (options%nsmooth_farfield .GT. 0)) then 
        if (options%cdisplay) then 
            write(*,'(A)') '    {smoothing}'
        end if
        if (options%nsmooth_interlayer .GT. 0) then 
            call inter_layer_smoothing(mesh,options)
        end if 
        if (options%nsmooth_farfield .GT. 0) then 
            call far_field_smoothing(mesh,options)
        end if 
    end if 

    !check for and eliminate small and degenerate cells along with multiple edges between the same two cells
    ndegenerate = 0 
    nsmall = 0 
    ndedge = 0
    do ii=1,mesh%nedge
        nupdate = 0
        call remove_edges_with_shared_cells(mesh,ntagged)
        ndedge = ndedge + ntagged
        nupdate = nupdate + ntagged
        call collapse_small_or_degenerate_cells(mesh,options,ntagged,'small')
        nsmall = nsmall + ntagged
        nupdate = nupdate + ntagged
        call collapse_small_or_degenerate_cells(mesh,options,ntagged,'degenerate')
        ndegenerate = ndegenerate + ntagged
        nupdate = nupdate + ntagged
        if (nupdate == 0) then 
            exit 
        end if 
    end do 
    if (options%cdisplay) then 
        write(*,'(A,I0,A)') '    {merged ',ndegenerate,' degenerate cells}'
    end if
    if (options%cdisplay) then 
        write(*,'(A,I0,A)') '    {merged ',nsmall,' small cells}'
    end if
    if (options%cdisplay) then 
        write(*,'(A,I0,A)') '    {eliminated ',ndedge,' cell-duplicate edges}'
    end if

    !set any boundary condition zones 
    if (options%nbczone .GT. 0) then 
        do ii=1,options%nbczone
            call set_boundary_condition_inzone(mesh,options%boundarycondition_zones(ii))
        end do 
    end if 

    !simplify surfaces
    if (options%simplify_surfaces) then 
        call simplify_surfaces(mesh)
    end if 

    !remove regions connected to specified boundary conditions 
    if (options%nbcremzone .GT. 0) then 
        do ii=1,options%nbcremzone
            call remove_cells_on_boundary_condition(mesh,options%boundarycondition_remove_zones(ii))
        end do 
    end if 

    !remove isolated mesh regions connected only to wall boundary conditions 



end if

!index mesh
call mesh%index_vertices()
call mesh%index_edges()
call mesh%get_v2e()
call mesh%get_cells()

!evaluate cell areas
call mesh%get_cell_areas()

!flag surface intersecting vertices 
call flag_surface_intersecting_vertices(mesh)

!display mesh properties
if (options%cdisplay) then 
    write(*,'(A)') '--> mesh construction completed:'
    write(*,'(A,I0,A)') '    {cells: ',mesh%ncell,'}'
    write(*,'(A,I0,A)') '    {edges: ',mesh%nedge,'}'
    write(*,'(A,I0,A)') '    {vertices: ',mesh%nvertex,'}'
    write(*,'(A,I0,A)') '    {surface intersection vertices: ',mesh%nvertex_surfint,'}'
    write(*,'(A,E11.5,A,E11.5,A,E11.5,A)') '    {cell volume (max/min/total) : ',maxval(mesh%cell(:)%volume),' / ',&
    minval(mesh%cell(:)%volume),' / ',sum(mesh%cell(:)%volume),'}'
end if


! !write tritree ======================
! print *, '## debug'

! print *, 'nvertex = ',tri_tree%nvertex 
! open(11,file='vertices')
! do ii=1,tri_tree%nvertex 
!     write(11,'(F0.16,A,F0.16,A,F0.16)') tri_tree%vertices(ii,1),' ',tri_tree%vertices(ii,2),' ',tri_tree%vertices(ii,3)
! end do 
! close(11)

! open(11,file='faces')
! do ii=1,tri_tree%nnode 
!     if (.NOT. associated(tri_tree%node(ii)%child1)) then 
!         write(11,'(I0,A,I0,A,I0)') tri_tree%node(ii)%vertex1,' ',tri_tree%node(ii)%vertex2,' ',tri_tree%node(ii)%vertex3
!     end if 
! end do 
! close(11)
! !write tritree ======================


! !write mesh_full =========================
! open(11,file='verticesf')
! do ii=1,mesh_full%nvertex
!     write(11,'(F0.16,A,F0.16,A,F0.16)') mesh_full%vertex(ii)%coordinate(1),' ',&
!     mesh_full%vertex(ii)%coordinate(2),' ',mesh_full%vertex(ii)%coordinate(3)
! end do 
! close(11)
! open(11,file='facesf')
! do ii=1,mesh_full%nedge
!     if (mesh_full%edge(ii)%flag) then 
!         write(11,'(I0,A,I0)') mesh_full%edge(ii)%vertex1%index,' ',mesh_full%edge(ii)%vertex2%index
!     end if 
! end do 
! close(11)
! print *, 'complete'
! !write mesh_full =========================


! !write mesh =========================
! print *, mesh%nvertex,' / ',mesh%nedge
! open(11,file='verticesf')
! do ii=1,mesh%nvertex
!     write(11,'(F0.16,A,F0.16,A,F0.16)') mesh%vertex(ii)%coordinate(1),' ',&
!     mesh%vertex(ii)%coordinate(2),' ',mesh%vertex(ii)%coordinate(3)
! end do 
! close(11)
! open(11,file='facesf')
! do ii=1,mesh%nedge
!     write(11,'(I0,A,I0,A,I0,A,I0)') mesh%edge(ii)%vertex1%index,' ',mesh%edge(ii)%vertex2%index,' ',&
!     mesh%edge(ii)%cell1,' ',mesh%edge(ii)%cell2
! end do 
! close(11)
! print *, 'complete'
! !write mesh =========================

return 
end function hex2d_mesh


!clip to geometry =========================
subroutine clip_to_geometry(mesh,mesh_full,geometry,kdtree,options)
implicit none 

!variables - import
type(hex_mesh), target :: mesh
type(hex_mesh), target :: mesh_full
type(halfedge), target :: geometry
type(gkdtree), target :: kdtree
type(hex_options) :: options 

!variables - local 
logical :: vtx_already_split
integer(in64) :: ii,jj,nn,ee,nfront,nfrontN,ndegenrate
integer(in64) :: nedge_mesh,ntgtnode,nedgecheck,n_split_vertex,nvidx,medge_tgt,vnew,etgt,vtgt,vins,eins,edge0
integer(in64) :: neidx(2),edge_split_vertex(mesh_full%nvertex)
integer(in64), dimension(:), allocatable :: efront,efrontN
integer(in64), dimension(:,:), allocatable :: geometry_ve0ends
real(dp) :: xmax,xmin,ymax,ymin,de_or,de_op,zelen_tol,edpv,splitdist1,splitdist2
real(dp) :: vi(3),vnormal(3),emidvec(3),edge_split_vertex_distance(mesh_full%nvertex)
type(node) :: tgt_nodes(kdtree%nnode)
type(edge) :: edgea
type(edge), dimension(:), allocatable :: edge_check

!set zero edge length tollerance 
zelen_tol = options%edgelength_min

!get the full mesh cell volumes 
call mesh_full%get_cells()
call mesh_full%get_cell_areas()

!set mesh vertex tags to zero (a non-zero tag indicates the vertex is shared with one on the clipped geometry)
do ii=1,size(mesh_full%vertex,dim=1)
    mesh_full%vertex(ii)%tag = 0 
end do 

!initialise external tags in the volume mesh
do ee=1,size(mesh_full%edge,dim=1)
    mesh_full%edge(ee)%external = .false.
end do 
do ee=1,size(mesh_full%vertex,dim=1)
    mesh_full%vertex(ee)%external = .false.
end do 

!set geometry flags
do ii=1,size(geometry%edge,dim=1)
    geometry%edge(ii)%flag = .false.
end do 

!initialise the intersection vertex
vi(:) = 0.0d0 

!clip mesh edges to the geometry 
ntgtnode = 0 
edge_split_vertex(:) = 0 
edge_split_vertex_distance(:) = 0.0d0 
nedge_mesh = mesh_full%nedge
allocate(edge_check(2*geometry%nedge))
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
    if (geometry%nedge .GT. size(edge_check,dim=1)) then 
        deallocate(edge_check)
        allocate(edge_check(2*geometry%nedge))
    end if 
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
    !(update the kdtree as each geometry edge is split)
    n_split_vertex = 0 
    do ee=1,nedgecheck
        if ((.NOT. edge_check(ee)%flag) .AND. (.NOT. edge_check(ee)%opposite%flag)) then 
            vi(1:2) = line_line_intersect(edge_check(ee)%origin%coordinate(1:2),edge_check(ee)%opposite%origin%coordinate(1:2),&
            mesh_full%edge(ii)%vertex1%coordinate(1:2),mesh_full%edge(ii)%vertex2%coordinate(1:2))
            if ((ieee_is_nan(vi(1))) .OR. (ieee_is_nan(vi(2)))) then !skip if invalid intersection
                cycle 
            end if  
            if (is_in_line_segment(edge_check(ee)%origin%coordinate(1:2),edge_check(ee)%opposite%origin%coordinate(1:2),vi) .AND. &
                is_in_line_segment(mesh_full%edge(ii)%vertex1%coordinate(1:2),mesh_full%edge(ii)%vertex2%coordinate(1:2),vi)) then 

                !check if the split vertex is within tollerance of either end of the geometry edge 
                de_or = norm2(vi(1:2) - geometry%edge(edge_check(ee)%index)%origin%coordinate(1:2))
                de_op = norm2(vi(1:2) - geometry%edge(edge_check(ee)%index)%opposite%origin%coordinate(1:2))

                !process the edge split
                if (de_or .LE. zelen_tol) then !intersects at edge%origin -> do not split the geometry edge
                    nvidx = geometry%edge(edge_check(ee)%index)%origin%index
                elseif (de_op .LE. zelen_tol) then !intersects at edge%opposite%origin -> do not split the geometry edge
                    nvidx = geometry%edge(edge_check(ee)%index)%opposite%origin%index
                else !intersects within the edge -> split the geometry edge

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
                end if 

                !store the new geometry vertex to use to split the mesh edge if it is not here already
                vtx_already_split = .false.
                do jj=1,n_split_vertex
                    if (edge_split_vertex(jj) == nvidx) then 
                        vtx_already_split = .true.
                        ! print *, '** warning: vertex already selected to split this edge'
                        exit
                    end if 
                end do 
                if (.NOT. vtx_already_split) then 
                    n_split_vertex = n_split_vertex + 1
                    edge_split_vertex(n_split_vertex) = nvidx
                end if 
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

        !check distance to each end of the mesh edge from the split vertex
        splitdist1 = norm2(geometry%vertex(edge_split_vertex(ee))%coordinate(1:2) - mesh_full%edge(ii)%vertex1%coordinate(1:2))
        splitdist2 = norm2(geometry%vertex(edge_split_vertex(ee))%coordinate(1:2) - mesh_full%edge(ii)%vertex2%coordinate(1:2))
        
        !split edge or link vertices depending on intersection location 
        if (splitdist1 .LE.  zelen_tol) then !do not split but set external state and tag vertex

            !tag this mesh vertex
            mesh_full%edge(medge_tgt)%vertex1%tag = edge_split_vertex(ee)

            !surface normal at this geometry vertex 
            vnormal = 0.5d0*(geometry%vertex(edge_split_vertex(ee))%edge%normal + &
                             geometry%vertex(edge_split_vertex(ee))%edge%previous%normal)

            !set the external status of the relevent section(s) of the edge 
            emidvec = 0.5d0*(mesh_full%edge(medge_tgt)%vertex1%coordinate + mesh_full%edge(medge_tgt)%vertex2%coordinate) - &
            geometry%vertex(edge_split_vertex(ee))%coordinate
            edpv = dot_product(emidvec,vnormal)
            if (edpv .GT. 0.0d0) then 
                mesh_full%edge(medge_tgt)%external = .true.
            end if 
        elseif (splitdist2 .LE.  zelen_tol) then !do not split but set external state and tag vertex
            
            !tag this mesh vertex
            mesh_full%edge(medge_tgt)%vertex2%tag = edge_split_vertex(ee)

            !surface normal at this geometry vertex 
            vnormal = 0.5d0*(geometry%vertex(edge_split_vertex(ee))%edge%normal + &
                             geometry%vertex(edge_split_vertex(ee))%edge%previous%normal)

            !set the external status of the relevent section(s) of the edge 
            emidvec = 0.5d0*(mesh_full%edge(medge_tgt)%vertex1%coordinate + mesh_full%edge(medge_tgt)%vertex2%coordinate) - &
            geometry%vertex(edge_split_vertex(ee))%coordinate
            edpv = dot_product(emidvec,vnormal)
            if (edpv .GT. 0.0d0) then 
                mesh_full%edge(medge_tgt)%external = .true.
            end if 
        else !split edge 

            !split the edge and tag the intersection vertex with its corresponding geometry vertex
            edge0 = medge_tgt
            call mesh_full%edge(medge_tgt)%split(mesh_full,geometry%vertex(edge_split_vertex(ee))%coordinate,medge_tgt,vnew)
            mesh_full%vertex(vnew)%tag = edge_split_vertex(ee)

            !surface normal at this geometry vertex 
            vnormal = 0.5d0*(geometry%vertex(edge_split_vertex(ee))%edge%normal + &
                             geometry%vertex(edge_split_vertex(ee))%edge%previous%normal)

            !set the external status of the relevent section(s) of the edge 
            emidvec = 0.5d0*(mesh_full%edge(edge0)%vertex1%coordinate + mesh_full%edge(edge0)%vertex2%coordinate) - &
            geometry%vertex(edge_split_vertex(ee))%coordinate
            edpv = dot_product(emidvec,vnormal)
            if (edpv .GT. 0.0d0) then 
                mesh_full%edge(edge0)%external = .true.
            elseif (edpv == 0.0d0) then 
                emidvec = 0.5d0*(mesh_full%edge(medge_tgt)%vertex1%coordinate + mesh_full%edge(medge_tgt)%vertex2%coordinate) - &
                geometry%vertex(edge_split_vertex(ee))%coordinate
                edpv = dot_product(emidvec,vnormal)
                if (edpv .LT. 0.0d0) then !opposite as checking against the next segment 
                    mesh_full%edge(edge0)%external = .true.
                elseif (edpv == 0.0d0) then 
                    print *, '** indeterminate volume-surface edge direction relation'
                    mesh%isvalid = .false.
                end if 
            end if 
            if (ee == n_split_vertex) then !set the final edge segment 
                emidvec = 0.5d0*(mesh_full%edge(medge_tgt)%vertex1%coordinate + mesh_full%edge(medge_tgt)%vertex2%coordinate) - &
                geometry%vertex(edge_split_vertex(ee))%coordinate
                edpv = dot_product(emidvec,vnormal)
                if (edpv .GT. 0.0d0) then 
                    mesh_full%edge(medge_tgt)%external = .true.
                elseif (edpv == 0.0d0) then 
                    emidvec = 0.5d0*(mesh_full%edge(edge0)%vertex1%coordinate + mesh_full%edge(edge0)%vertex2%coordinate) - &
                    geometry%vertex(edge_split_vertex(ee))%coordinate
                    edpv = dot_product(emidvec,vnormal)
                    if (edpv .LT. 0.0d0) then !opposite as checking against the previous segment 
                        mesh_full%edge(medge_tgt)%external = .true.
                    elseif (edpv == 0.0d0) then 
                        print *, '** indeterminate volume-surface edge direction relation (final segment)'
                        mesh%isvalid = .false.
                    end if 
                end if 
            end if
        end if 
    end do 
end do 

!tag vertices on already tagged external edges as external 
do ee=1,mesh_full%nedge
    if (mesh_full%edge(ee)%external) then 
        mesh_full%edge(ee)%vertex1%external = .true.
        mesh_full%edge(ee)%vertex2%external = .true.
    end if 
end do 

!for each vertex in the geometry find its fraction along its original geometry edge and store the indecies of the end two vertices 
allocate(geometry_ve0ends(geometry%nvertex,2))
do ii=1,geometry%nvertex
    edgea = geometry%vertex(ii)%edge
    geometry%vertex(ii)%rdata = fraction_along_line(edgea%origin0%coordinate(1:2),edgea%opposite%origin0%coordinate(1:2),&
                                                    geometry%vertex(ii)%coordinate(1:2))
    geometry_ve0ends(ii,1) = edgea%origin0%index
    geometry_ve0ends(ii,2) = edgea%opposite%origin0%index
end do 

!flag all the 'external' halfedges in the geometry, those with the correctly oriented normal vector to be part of the mesh 
call flag_external_halfedges(geometry)

!build mesh v2e
call mesh_full%index_edges()
call mesh_full%get_v2e()

!flood the geometry external portion of the mesh (not flooding through tagged vertices)
nfront = 0
nfrontN = 0
allocate(efront(mesh_full%nedge))
allocate(efrontN(mesh_full%nedge))
do ee=1,mesh_full%nedge
    if (mesh_full%edge(ee)%external) then 
        nfront = nfront + 1
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
            do ii=1,mesh_full%max_valence
                if (mesh_full%v2e(vtgt,ii) .NE. 0) then 
                    if (.NOT. mesh_full%edge(mesh_full%v2e(vtgt,ii))%external) then 
                        nfrontN = nfrontN + 1
                        efrontN(nfrontN) = mesh_full%v2e(vtgt,ii)
                        mesh_full%edge(mesh_full%v2e(vtgt,ii))%external = .true.
                    end if 
                end if 
            end do 
        end if 

        !flood across vertex 2
        vtgt = mesh_full%edge(etgt)%vertex2%index
        if (mesh_full%vertex(vtgt)%tag == 0) then 
            do ii=1,mesh_full%max_valence
                if (mesh_full%v2e(vtgt,ii) .NE. 0) then 
                    if (.NOT. mesh_full%edge(mesh_full%v2e(vtgt,ii))%external) then 
                        nfrontN = nfrontN + 1
                        efrontN(nfrontN) = mesh_full%v2e(vtgt,ii)
                        mesh_full%edge(mesh_full%v2e(vtgt,ii))%external = .true.
                    end if 
                end if 
            end do 
        end if 
    end do 
    if (nfrontN == 0) then 
        exit 
    end if 
    nfront = nfrontN
    efront(1:nfront) = efrontN(1:nfront) 
end do 

!flag the items to retain 
do ee=1,mesh_full%nedge
    mesh_full%edge(ee)%flag = .false.
end do 
do ee=1,mesh_full%nvertex
    mesh_full%vertex(ee)%flag = .false.
end do 
if (options%mesh_in_out == 'external') then 
    do ee=1,mesh_full%nedge
        if (mesh_full%edge(ee)%external) then 
            mesh_full%edge(ee)%flag = .true.
        end if 
    end do 
    do ee=1,mesh_full%nvertex
        if (mesh_full%vertex(ee)%external) then 
            mesh_full%vertex(ee)%flag = .true.
        end if 
    end do 
elseif (options%mesh_in_out == 'internal') then 
    do ee=1,mesh_full%nedge
        if (.NOT. mesh_full%edge(ee)%external) then 
            mesh_full%edge(ee)%flag = .true.
        end if 
    end do 
    do ee=1,mesh_full%nvertex
        if (.NOT. mesh_full%vertex(ee)%external) then 
            mesh_full%vertex(ee)%flag = .true.
        end if 
    end do
end if 

!untag geometry linked vertices so as to not double-copy them 
do ee=1,mesh_full%nvertex
    if (mesh_full%vertex(ee)%tag .NE. 0) then 
        mesh_full%vertex(ee)%flag = .false.
    end if
end do 

!count items to build the final mesh 
mesh%nvertex = 0
do ii=1,mesh_full%nvertex
    if (mesh_full%vertex(ii)%flag) then 
        mesh%nvertex = mesh%nvertex + 1 
    end if 
end do 
do ii=1,geometry%nvertex
    mesh%nvertex = mesh%nvertex + 1 
end do 
mesh%nedge = 0
do ii=1,mesh_full%nedge
    if (mesh_full%edge(ii)%flag) then 
        mesh%nedge = mesh%nedge + 1 
    end if 
end do 
do ii=1,geometry%nedge
    if (geometry%edge(ii)%flag) then 
        mesh%nedge = mesh%nedge + 1 
    end if 
end do 
mesh%ncell = 0

!build the complete final mesh by unifying mesh_full with with the oriented geometry mesh
allocate(mesh%vertex(mesh%nvertex))
allocate(mesh%edge(mesh%nedge))
vins = 0 
do ii=1,mesh_full%nvertex
    if (mesh_full%vertex(ii)%flag) then 
        vins = vins + 1
        mesh_full%vertex(ii)%index = vins
    end if 
end do 
do ii=1,geometry%nvertex
    vins = vins + 1
    geometry%vertex(ii)%index = vins
end do 
eins = 0
do ii=1,mesh_full%nedge
    if (mesh_full%edge(ii)%flag) then 
        eins = eins + 1
        mesh_full%edge(ii)%index = eins
    end if 
end do 
do ii=1,geometry%nedge
    if (geometry%edge(ii)%flag) then 
        eins = eins + 1
        geometry%edge(ii)%index = eins
    end if 
end do 
do ii=1,mesh_full%nvertex
    if (mesh_full%vertex(ii)%flag) then 
        vins = mesh_full%vertex(ii)%index
        mesh%vertex(vins)%flag = .false.
        mesh%vertex(vins)%tag = 0 
        mesh%vertex(vins)%index = vins
        mesh%vertex(vins)%coordinate = mesh_full%vertex(ii)%coordinate
        mesh%vertex(vins)%coordinate(3) = 0.0d0 
        mesh%vertex(vins)%rdata = -1.0d0 
    end if 
end do 
do ii=1,geometry%nvertex
    vins = geometry%vertex(ii)%index
    mesh%vertex(vins)%flag = .false.
    mesh%vertex(vins)%tag = 0 
    mesh%vertex(vins)%index = vins
    mesh%vertex(vins)%coordinate = geometry%vertex(ii)%coordinate
    mesh%vertex(vins)%coordinate(3) = 0.0d0 
    mesh%vertex(vins)%rdata = geometry%vertex(ii)%rdata 
    allocate(mesh%vertex(vins)%ivdata(2))
    mesh%vertex(vins)%ivdata = geometry_ve0ends(ii,:)
end do 
do ii=1,mesh_full%nedge
    if (mesh_full%edge(ii)%flag) then 
        eins = mesh_full%edge(ii)%index
        mesh%edge(eins)%flag = .false.
        mesh%edge(eins)%tag = 0 
        mesh%edge(eins)%index = eins
        mesh%edge(eins)%cell1 = mesh_full%edge(ii)%cell1
        mesh%edge(eins)%cell2 = mesh_full%edge(ii)%cell2
        mesh%edge(eins)%direction = mesh_full%edge(ii)%direction
        if (mesh_full%edge(ii)%vertex1%tag == 0) then 
            mesh%edge(eins)%vertex1 => mesh%vertex(mesh_full%edge(ii)%vertex1%index)
        else
            mesh%edge(eins)%vertex1 => mesh%vertex(geometry%vertex(mesh_full%edge(ii)%vertex1%tag)%index)
        end if 
        if (mesh_full%edge(ii)%vertex2%tag == 0) then 
            mesh%edge(eins)%vertex2 => mesh%vertex(mesh_full%edge(ii)%vertex2%index)
        else
            mesh%edge(eins)%vertex2 => mesh%vertex(geometry%vertex(mesh_full%edge(ii)%vertex2%tag)%index)
        end if 
    end if 
end do 
do ii=1,geometry%nedge
    if (geometry%edge(ii)%flag) then 
        eins = geometry%edge(ii)%index 
        mesh%edge(eins)%flag = .false.
        mesh%edge(eins)%tag = 0 
        mesh%edge(eins)%index = eins
        mesh%edge(eins)%cell1 = 0
        mesh%edge(eins)%cell2 = -1
        if (options%mesh_in_out == 'external') then 
            mesh%edge(eins)%vertex2 => mesh%vertex(geometry%edge(ii)%origin%index)
            mesh%edge(eins)%vertex1 => mesh%vertex(geometry%edge(ii)%opposite%origin%index)
            mesh%edge(eins)%direction = -geometry%edge(ii)%direction
        elseif (options%mesh_in_out == 'internal') then 
            mesh%edge(eins)%vertex1 => mesh%vertex(geometry%edge(ii)%origin%index)
            mesh%edge(eins)%vertex2 => mesh%vertex(geometry%edge(ii)%opposite%origin%index)
            mesh%edge(eins)%direction = geometry%edge(ii)%direction
        end if 
    end if 
end do 

!eliminate any degenerate edges joining the same vertex
call collapse_degenerate_edges(mesh,ndegenrate)
if (options%cdisplay) then 
    write(*,'(A,I0,A)') '    {eliminated ',ndegenrate,' degenerate edges}'
end if

!initialise the cell reference volumes 
mesh%ncell = mesh_full%ncell 
allocate(mesh%cell(mesh%ncell))
do ii=1,mesh%ncell
    mesh%cell(ii)%volume_ref = mesh_full%cell(ii)%volume 
end do 

!set the cell association of each geometry surface edge 
if (options%geom_cell_link_method == 'dot_dir') then 
    call set_geometry_edge_cell_associations_dpdir(mesh,geometry)
elseif (options%geom_cell_link_method == 'sflood') then 
    call set_geometry_edge_cell_associations_sflood(mesh,geometry)
elseif (options%geom_cell_link_method == 'auto') then 
    call set_geometry_edge_cell_associations_auto(mesh,geometry)
end if 
return 
end subroutine clip_to_geometry


!set geometry surface cell association (automatic method selection) =========================
subroutine set_geometry_edge_cell_associations_auto(mesh,geometry)
implicit none 

!variables - import 
type(hex_mesh), target :: mesh
type(halfedge), target :: geometry
    
!variables - local 
integer(in64) :: ii,ee
integer(in64) :: vedge_idx,avedge_idx,vtx_attached0,vtx_attached1
integer(in64) :: vertex0,edgeF,vertexF,cshared
real(dp) :: dp_dir
real(dp) :: enormal(2),edir(2)
type(hex_edge) :: vol_edge

!build mesh v2e
call mesh%index_edges()
call mesh%get_v2e()

!index mesh cell edges 
call mesh%get_cells()

!flag all vertices on volume edges 
do ii=1,mesh%nvertex
    mesh%vertex(ii)%flag = .false.
end do 
do ee=1,mesh%nedge
    if ((mesh%edge(ee)%cell1 .NE. -1) .AND. (mesh%edge(ee)%cell2 .NE. -1)) then 
        mesh%edge(ee)%vertex1%flag = .true.
        mesh%edge(ee)%vertex2%flag = .true.
    end if 
end do 

!set association of geometry surface edges that are attached to external edges 
vedge_idx = 0
vtx_attached0 = 0
vtx_attached1 = 0
do ee=1,mesh%nedge
    if (mesh%edge(ee)%cell2 == -1) then 

        !check if this edge is connected to a volume edge (and find the connected edge with the minimum angle)
        call get_min_angle_vedge(mesh,mesh%edge(ee),0_in64,vol_edge,vtx_attached0)

        !if a volume edge has been found 
        if (vol_edge%index .GT. 0) then 

            !index of the edge 
            vedge_idx = vol_edge%index

            !find orientation of this edge WRT the volume edge normal 
            enormal(1) = -mesh%edge(vedge_idx)%direction(2)
            enormal(2) = mesh%edge(vedge_idx)%direction(1)
            if (vtx_attached0 == 1) then 
                edir = mesh%edge(ee)%direction(1:2)
            elseif (vtx_attached0 == 2) then
                edir = -mesh%edge(ee)%direction(1:2)
            end if 
            dp_dir = dot_product(edir,enormal)

            !set the cell of the surface edge 
            if (dp_dir .GT. 0.0d0) then !set from direction
                mesh%edge(ee)%cell1 = mesh%edge(vedge_idx)%cell1
            elseif (dp_dir .LT. 0.0d0) then !set from direction
                mesh%edge(ee)%cell1 = mesh%edge(vedge_idx)%cell2
            else !set from cell flood if the direction is ambigous 

                !traverse the surface until another volume connected edge is located 
                if (vtx_attached0 == 1) then 
                    vertex0 = mesh%edge(ee)%vertex1%index
                else
                    vertex0 = mesh%edge(ee)%vertex2%index
                end if 
                call traverse_surface_to_flagged_vertex(mesh,ee,vertex0,edgeF,vertexF)
                if (vertexF == -1) then 
                    print *, '** invalid surface connectivity detected (geometry cell linking)'
                    mesh%isvalid = .false.
                    cycle 
                end if 
                
                !find the minimum angle volume edge attached to this vertex
                if (vertexF == mesh%edge(edgeF)%vertex1%index) then 
                    call get_min_angle_vedge(mesh,mesh%edge(edgeF),1_in64,vol_edge,vtx_attached1)
                else
                    call get_min_angle_vedge(mesh,mesh%edge(edgeF),2_in64,vol_edge,vtx_attached1)
                end if 

                !index of the edge 
                avedge_idx = vol_edge%index

                !switch method on edge-cell association 
                if ((mesh%edge(vedge_idx)%cell1 == mesh%edge(avedge_idx)%cell1) .AND. &
                    (mesh%edge(vedge_idx)%cell2 == mesh%edge(avedge_idx)%cell2)) then !both volume edges join the same two cells (this geometry section partitions one edge from the un-clipped mesh)

                    !check the relative connectivity orientation of the volume and surface edges to select a cell 
                    cshared = 0
                    if (vtx_attached0 == 1) then !initial vertex is the start of the surface edge 
                        if (mesh%edge(vedge_idx)%vertex1%index == vertex0) then !initial vertex is the start of the volume edge 
                            cshared = mesh%edge(vedge_idx)%cell2
                        else !initial vertex is the end of the volume edge 
                            cshared = mesh%edge(vedge_idx)%cell1
                        end if  
                    else !initial vertex is the end of the surface edge 
                        if (mesh%edge(vedge_idx)%vertex1%index == vertex0) then !initial vertex is the start of the volume edge 
                            cshared = mesh%edge(vedge_idx)%cell1
                        else !initial vertex is the end of the volume edge 
                            cshared = mesh%edge(vedge_idx)%cell2
                        end if 
                    end if 

                    !set this cell for the surface edge 
                    mesh%edge(ee)%cell1 = cshared
                else !each volume edge joins two pairs of different cells, thus the cell this geometry belongs to is the one cell shared between both edges

                    !find the cell shared between these two edges 
                    cshared = 0
                    if (mesh%edge(vedge_idx)%cell1 == mesh%edge(avedge_idx)%cell1) then 
                        cshared = mesh%edge(vedge_idx)%cell1
                    elseif (mesh%edge(vedge_idx)%cell2 == mesh%edge(avedge_idx)%cell2) then 
                        cshared = mesh%edge(vedge_idx)%cell2
                    elseif (mesh%edge(vedge_idx)%cell1 == mesh%edge(avedge_idx)%cell2) then 
                        cshared = mesh%edge(vedge_idx)%cell1
                    elseif (mesh%edge(vedge_idx)%cell2 == mesh%edge(avedge_idx)%cell1) then 
                        cshared = mesh%edge(vedge_idx)%cell2 
                    else
                        print *, '** pair of edges that do not share a cell selected (geometry cell linking)'
                        print '(A,I0)', '** unable to determine the cell association of edge: ',ee
                        ! print *, vedge_idx,avedge_idx
                        ! print *, ee,edgeF
                        mesh%isvalid = .false.
                        cycle 
                    end if 

                    !set this cell for the surface edge 
                    mesh%edge(ee)%cell1 = cshared
                end if 
            end if 
        end if 
    end if 
end do 

!flood geometry edge cell associations through the surface mesh 
call flood_geometry_edge_cells(mesh,geometry)
return 
end subroutine set_geometry_edge_cell_associations_auto


!set geometry surface cell association (from flood along surface) =========================
subroutine set_geometry_edge_cell_associations_sflood(mesh,geometry)
implicit none 

!variables - import 
type(hex_mesh), target :: mesh
type(halfedge), target :: geometry
    
!variables - local 
integer(in64) :: ii,ee
integer(in64) :: vedge_idx,avedge_idx,vtx_attached0,vtx_attached1
integer(in64) :: vertex0,edgeF,vertexF,cshared
type(hex_edge) :: vol_edge

!build mesh v2e
call mesh%index_edges()
call mesh%get_v2e()

!index mesh cell edges 
call mesh%get_cells()

!flag all vertices on volume edges 
do ii=1,mesh%nvertex
    mesh%vertex(ii)%flag = .false.
end do 
do ee=1,mesh%nedge
    if ((mesh%edge(ee)%cell1 .NE. -1) .AND. (mesh%edge(ee)%cell2 .NE. -1)) then 
        mesh%edge(ee)%vertex1%flag = .true.
        mesh%edge(ee)%vertex2%flag = .true.
    end if 
end do 

!set association of geometry surface edges that are attached to external edges 
vedge_idx = 0
vtx_attached0 = 0
vtx_attached1 = 0
do ee=1,mesh%nedge
    if (mesh%edge(ee)%cell2 == -1) then 

        !check if this edge is connected to a volume edge (and find the connected edge with the minimum angle)
        call get_min_angle_vedge(mesh,mesh%edge(ee),0_in64,vol_edge,vtx_attached0)

        !if a volume edge has been found 
        if (vol_edge%index .GT. 0) then 

            !index of the edge 
            vedge_idx = vol_edge%index

            !traverse the surface until another volume connected edge is located 
            if (vtx_attached0 == 1) then 
                vertex0 = mesh%edge(ee)%vertex1%index
            else
                vertex0 = mesh%edge(ee)%vertex2%index
            end if 
            call traverse_surface_to_flagged_vertex(mesh,ee,vertex0,edgeF,vertexF)
            if (vertexF == -1) then 
                print *, '** invalid surface connectivity detected (geometry cell linking)'
                mesh%isvalid = .false.
                cycle 
            end if 
            
            !find the minimum angle volume edge attached to this vertex
            if (vertexF == mesh%edge(edgeF)%vertex1%index) then 
                call get_min_angle_vedge(mesh,mesh%edge(edgeF),1_in64,vol_edge,vtx_attached1)
            else
                call get_min_angle_vedge(mesh,mesh%edge(edgeF),2_in64,vol_edge,vtx_attached1)
            end if 

            !index of the edge 
            avedge_idx = vol_edge%index

            !switch method on edge-cell association 
            if ((mesh%edge(vedge_idx)%cell1 == mesh%edge(avedge_idx)%cell1) .AND. &
                (mesh%edge(vedge_idx)%cell2 == mesh%edge(avedge_idx)%cell2)) then !both volume edges join the same two cells (this geometry section partitions one edge from the un-clipped mesh)

                !check the relative connectivity orientation of the volume and surface edges to select a cell 
                cshared = 0
                if (vtx_attached0 == 1) then !initial vertex is the start of the surface edge 
                    if (mesh%edge(vedge_idx)%vertex1%index == vertex0) then !initial vertex is the start of the volume edge 
                        cshared = mesh%edge(vedge_idx)%cell2
                    else !initial vertex is the end of the volume edge 
                        cshared = mesh%edge(vedge_idx)%cell1
                    end if  
                else !initial vertex is the end of the surface edge 
                    if (mesh%edge(vedge_idx)%vertex1%index == vertex0) then !initial vertex is the start of the volume edge 
                        cshared = mesh%edge(vedge_idx)%cell1
                    else !initial vertex is the end of the volume edge 
                        cshared = mesh%edge(vedge_idx)%cell2
                    end if 
                end if 

                !set this cell for the surface edge 
                mesh%edge(ee)%cell1 = cshared
            else !each volume edge joins two pairs of different cells, thus the cell this geometry belongs to is the one cell shared between both edges

                !find the cell shared between these two edges 
                cshared = 0
                if (mesh%edge(vedge_idx)%cell1 == mesh%edge(avedge_idx)%cell1) then 
                    cshared = mesh%edge(vedge_idx)%cell1
                elseif (mesh%edge(vedge_idx)%cell2 == mesh%edge(avedge_idx)%cell2) then 
                    cshared = mesh%edge(vedge_idx)%cell2
                elseif (mesh%edge(vedge_idx)%cell1 == mesh%edge(avedge_idx)%cell2) then 
                    cshared = mesh%edge(vedge_idx)%cell1
                elseif (mesh%edge(vedge_idx)%cell2 == mesh%edge(avedge_idx)%cell1) then 
                    cshared = mesh%edge(vedge_idx)%cell2 
                else
                    print *, '** pair of edges that do not share a cell selected (geometry cell linking)'
                    print '(A,I0)', '** unable to determine the cell association of edge: ',ee
                    print *, mesh%edge(vedge_idx)%cell1,mesh%edge(vedge_idx)%cell2
                    print *, mesh%edge(avedge_idx)%cell1,mesh%edge(avedge_idx)%cell2
                    print *, vedge_idx,avedge_idx
                    ! print *, ee,edgeF
                    mesh%isvalid = .false.
                    cycle 
                end if 

                !set this cell for the surface edge 
                mesh%edge(ee)%cell1 = cshared
            end if 
        end if 
    end if 
end do 

!flood geometry edge cell associations through the surface mesh 
call flood_geometry_edge_cells(mesh,geometry)
return 
end subroutine set_geometry_edge_cell_associations_sflood


!traverse surface until flagged vertex =========================
subroutine traverse_surface_to_flagged_vertex(mesh,edge0,vertex0,edgeF,vertexF)
implicit none 

!variables - import 
integer(in64) :: edge0,vertex0,edgeF,vertexF
type(hex_mesh), target :: mesh

!variables - local 
integer(in64) :: ii,ff
integer(in64) :: edgec,vertexc,edgen,vertexn,etgt

!initialise final edge and vertex
edgeF = -1 
vertexF = -1

!check other end of base edge 
if (vertex0 == mesh%edge(edge0)%vertex1%index) then 
    if (mesh%edge(edge0)%vertex2%flag) then 
        edgeF = edge0
        vertexF = mesh%edge(edge0)%vertex2%index
    end if 
else
    if (mesh%edge(edge0)%vertex1%flag) then 
        edgeF = edge0 
        vertexF = mesh%edge(edge0)%vertex1%index
    end if
end if 

!flood if neccesary
if (edgeF == -1) then 

    !set current edge and vertex
    edgec = edge0 
    vertexc = vertex0

    !flood 
    do ff=1,mesh%nedge*10

        !initialise 
        edgen = 0
        vertexn = 0 

        !find next vertex
        if (vertexc == mesh%edge(edgec)%vertex1%index) then 
            vertexn = mesh%edge(edgec)%vertex2%index
        else
            vertexn = mesh%edge(edgec)%vertex1%index
        end if 

        !check if vertexn is flagged then exit if so 
        if (mesh%vertex(vertexn)%flag) then 
            edgeF = edgec 
            vertexF = vertexn
            exit 
        end if 

        !find next edge 
        do ii=1,mesh%max_valence
            etgt = mesh%v2e(vertexn,ii)
            if (etgt .NE. 0) then 
                if (etgt .NE. edgec) then 
                    if ((mesh%edge(etgt)%cell1 == -1) .OR. (mesh%edge(etgt)%cell2 == -1)) then
                        edgen = etgt
                        exit 
                    end if 
                end if 
            end if 
        end do 

        !update
        edgec = edgen 
        vertexc = vertexn
    end do 
end if 
return 
end subroutine traverse_surface_to_flagged_vertex


!set geometry surface cell association (from geometric directions) =========================
subroutine set_geometry_edge_cell_associations_dpdir(mesh,geometry)
implicit none 

!variables - import 
type(hex_mesh), target :: mesh
type(halfedge), target :: geometry

!variables - local 
integer(in64) :: ee
integer(in64) :: vedge_idx,vtx_attached
real(dp) :: dp_dir
real(dp) :: enormal(2),edir(2)
type(hex_edge) :: vol_edge

!build mesh v2e
call mesh%index_edges()
call mesh%get_v2e()

!set association of geometry surface edges that are attached to external edges 
vedge_idx = 0
do ee=1,mesh%nedge
    if (mesh%edge(ee)%cell2 == -1) then 

        !check if this edge is connected to a volume edge (and find the connected edge with the minimum angle)
        call get_min_angle_vedge(mesh,mesh%edge(ee),0_in64,vol_edge,vtx_attached)

        !if a volume edge has been found 
        if (vol_edge%index .GT. 0) then 

            !index of the edge 
            vedge_idx = vol_edge%index

            !find orientation of this edge WRT the volume edge normal 
            enormal(1) = -mesh%edge(vedge_idx)%direction(2)
            enormal(2) = mesh%edge(vedge_idx)%direction(1)
            if (vtx_attached == 1) then 
                edir = mesh%edge(ee)%direction(1:2)
            elseif (vtx_attached == 2) then
                edir = -mesh%edge(ee)%direction(1:2)
            end if 
            dp_dir = dot_product(edir,enormal)

            !set the cell of the surface edge 
            if (dp_dir .GT. 0.0d0) then 
                mesh%edge(ee)%cell1 = mesh%edge(vedge_idx)%cell1
            elseif (dp_dir .LT. 0.0d0) then 
                mesh%edge(ee)%cell1 = mesh%edge(vedge_idx)%cell2
            else
                print *, '** indeterminate volume-surface edge direction relation (geometry cell linking)'
                mesh%isvalid = .false.
            end if 
        end if 
    end if 
end do 

!flood geometry edge cell associations through the surface mesh 
call flood_geometry_edge_cells(mesh,geometry)
return 
end subroutine set_geometry_edge_cell_associations_dpdir


!get minimum angle volume edge attached to a surface edge  =========================
subroutine get_min_angle_vedge(mesh,tgt_edge,search_end,vol_edge,vtx_attached)
implicit none

!variables - import 
integer(in64) :: search_end,vtx_attached
type(hex_mesh), target :: mesh
type(hex_edge), target :: tgt_edge,vol_edge

!find edge on target end switch 
vtx_attached = 0 
if (search_end == 0) then !search both ends 
    call get_minimal_angle_vedge_on_svtx(mesh,tgt_edge,tgt_edge%vertex1,vol_edge)
    vtx_attached = 1
    if (vol_edge%index == -1) then 
        call get_minimal_angle_vedge_on_svtx(mesh,tgt_edge,tgt_edge%vertex2,vol_edge)
        vtx_attached = 2
    end if 
elseif (search_end == 1) then !search end 1
    call get_minimal_angle_vedge_on_svtx(mesh,tgt_edge,tgt_edge%vertex1,vol_edge)
    vtx_attached = 1
elseif (search_end == 2) then !search end 2
    call get_minimal_angle_vedge_on_svtx(mesh,tgt_edge,tgt_edge%vertex2,vol_edge)
    vtx_attached = 2
end if 
return 
end subroutine get_min_angle_vedge


!get the minimal angle volume edge connected to vertex 'tgt_vertex' on surface edge 'edge' =========================
subroutine get_minimal_angle_vedge_on_svtx(mesh,tgt_edge,tgt_vertex,vol_edge)
implicit none 

!variables - import 
type(hex_edge), target :: tgt_edge,vol_edge
type(hex_vertex), target :: tgt_vertex
type(hex_mesh), target :: mesh

!variables - local
integer(in64) :: ii 
integer(in64) :: etgt,aedge_idx
real(dp) :: vedge_angle
real(dp) :: surfdir(3),voldir(3)

!initialise vol_edge  
vol_edge%flag = .false.
vol_edge%index = -1
vol_edge%vertex1 => null()
vol_edge%vertex2 => null()

!initialise directions 
voldir(:) = 0.0d0 
surfdir(:) = 0.0d0 

!set surface direction 
if (associated(tgt_edge%vertex1,tgt_vertex)) then 
    surfdir(1:2) = tgt_edge%direction(1:2)
else
    surfdir(1:2) = -tgt_edge%direction(1:2)
end if 

!search about vertex
aedge_idx = 0 
vedge_angle = 10000.0d0 
do ii=1,mesh%max_valence
    etgt = mesh%v2e(tgt_vertex%index,ii) 
    if (etgt > 0) then 
        ! if ((mesh%edge(etgt)%cell1 > 0) .AND. (mesh%edge(etgt)%cell2 > 0)) then !exclude far field tagged edges 
        if ((mesh%edge(etgt)%cell1 .NE. -1) .AND. (mesh%edge(etgt)%cell2 .NE. -1)) then !include far field tagged edges 
            if (associated(mesh%edge(etgt)%vertex1,tgt_vertex)) then 
                voldir(1:2) = mesh%edge(etgt)%direction(1:2) 
            else
                voldir(1:2) = -mesh%edge(etgt)%direction(1:2)
            end if 
            if (aedge_idx == 0) then 
                aedge_idx = etgt
                vol_edge = mesh%edge(etgt)
                vol_edge%flag = .true. !tag edge found
                vedge_angle = vec2vec_angle(surfdir,voldir)
            else
                if (vec2vec_angle(surfdir,voldir) .LT. vedge_angle) then 
                    aedge_idx = etgt
                    vol_edge = mesh%edge(etgt)
                    vol_edge%flag = .true. !tag edge found
                    vedge_angle = vec2vec_angle(surfdir,voldir)
                end if 
            end if 
        end if 
    end if 
end do 
return 
end subroutine get_minimal_angle_vedge_on_svtx


!flood geometry cell associations =========================
subroutine flood_geometry_edge_cells(mesh,geometry)
implicit none 

!variables - import 
type(hex_mesh), target :: mesh
type(halfedge), target :: geometry

!variables - local 
integer(in64) :: ii,ee 
integer(in64) :: nupdate,meidx,meidxa

!flood 
do ii=1,geometry%nedge
    nupdate = 0 
    do ee=1,geometry%nedge
        if (geometry%edge(ee)%flag) then 
            meidx = geometry%edge(ee)%index 
            if (mesh%edge(meidx)%cell1 .NE. 0) then 
                if (geometry%edge(ee)%previous%flag) then 
                    meidxa = geometry%edge(ee)%previous%index
                    if (mesh%edge(meidxa)%cell1 == 0) then 
                        mesh%edge(meidxa)%cell1 = mesh%edge(meidx)%cell1
                        nupdate = nupdate + 1
                    end if 
                end if 
                if (geometry%edge(ee)%next%flag) then 
                    meidxa = geometry%edge(ee)%next%index
                    if (mesh%edge(meidxa)%cell1 == 0) then 
                        mesh%edge(meidxa)%cell1 = mesh%edge(meidx)%cell1
                        nupdate = nupdate + 1
                    end if 
                end if 
            end if 
        end if 
    end do 
    if (nupdate == 0) then 
        exit 
    end if 
end do 
return 
end subroutine flood_geometry_edge_cells


!perturb on geometry mesh vertices =========================
subroutine perturb_ongeom_vertices(mesh,kdtree,dv,nperturb)
implicit none 

!variables - import
integer(in64) :: nperturb
real(dp) :: dv 
type(hex_mesh), target :: mesh
type(gkdtree), target :: kdtree

!variables - local 
integer(in64) :: ii,jj,nn
integer(in64) :: ntgtnode
real(dp) :: xmin,xmax,ymin,ymax,dmin,dx,dy,enorm
real(dp) :: vidmin(3)
type(edge) :: etgt,emin 
type(node) :: tgt_nodes(kdtree%nnode)

!check and perturb vertices 
nperturb = 0
vidmin(:) = 0.0d0 
do ii=1,mesh%nvertex

    !check if this vertex is near the geometry 
    xmin = mesh%vertex(ii)%coordinate(1) - 2.0d0*dv
    xmax = mesh%vertex(ii)%coordinate(1) + 2.0d0*dv
    ymin = mesh%vertex(ii)%coordinate(2) - 2.0d0*dv
    ymax = mesh%vertex(ii)%coordinate(2) + 2.0d0*dv
    call kdtree%nodes_overlapping_region_2d(xmin,xmax,ymin,ymax,ntgtnode,tgt_nodes)
    if (ntgtnode == 0) then 
        cycle
    end if 

    !check if within dv of any edge 
    dmin = huge(dp)
    do nn=1,ntgtnode
        do jj=1,tgt_nodes(nn)%nitem
            etgt = tgt_nodes(nn)%edge(jj)%edge
            vidmin = min_dist_point_to_line(etgt%origin%coordinate(1:2),etgt%opposite%origin%coordinate(1:2),&
                                            mesh%vertex(ii)%coordinate(1:2))
            if (vidmin(3) .LT. dmin) then 
                dmin = vidmin(3)
                emin = etgt
            end if 
        end do 
    end do 

    !if within proximity then perturb 
    if (dmin .LE. dv) then 
        if (norm2(mesh%vertex(ii)%coordinate(1:2) - etgt%origin%coordinate(1:2)) .LE. 2.0d0*dv) then 
            dx = etgt%normal(1) + etgt%previous%normal(1)
            dy = etgt%normal(2) + etgt%previous%normal(2)
            enorm = sqrt(dx*dx + dy*dy)
            dx = dx/enorm 
            dy = dy/enorm 
        elseif (norm2(mesh%vertex(ii)%coordinate(1:2) - etgt%opposite%origin%coordinate(1:2)) .LE. 2.0d0*dv) then 
            dx = etgt%normal(1) + etgt%next%normal(1)
            dy = etgt%normal(2) + etgt%next%normal(2)
            enorm = sqrt(dx*dx + dy*dy)
            dx = dx/enorm 
            dy = dy/enorm 
        else
            dx = etgt%normal(1)
            dy = etgt%normal(2)
        end if 
        mesh%vertex(ii)%coordinate(1) = mesh%vertex(ii)%coordinate(1) + (1.1d0*dx)*dv
        mesh%vertex(ii)%coordinate(2) = mesh%vertex(ii)%coordinate(2) + (1.1d0*dy)*dv
        nperturb = nperturb + 1
    end if 
end do 
return 
end subroutine perturb_ongeom_vertices


!build oriented geometry mesh =========================
subroutine flag_external_halfedges(geometry)
implicit none 

!variables - import
type(halfedge), target :: geometry

!variables - local 
integer(in64) :: ii
real(dp) :: dx,dy
real(dp) :: normal(3)

!set geometry flags
do ii=1,geometry%nedge
    geometry%edge(ii)%flag = .false.
end do 

!set flag on halfedges with normals that point in the correct direction to true 
normal(:) = 0.0d0 
do ii=1,geometry%nedge
    dx = geometry%edge(ii)%opposite%origin%coordinate(1) - geometry%edge(ii)%origin%coordinate(1)
    dy = geometry%edge(ii)%opposite%origin%coordinate(2) - geometry%edge(ii)%origin%coordinate(2)
    normal(1) = dy 
    normal(2) = -dx 
    if (dot_product(normal,geometry%edge(ii)%normal) > 0.0) then 
        geometry%edge(ii)%flag = .true.
    end if 
end do 
return 
end subroutine flag_external_halfedges


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


!flag surface intersecting vertices =========================
subroutine flag_surface_intersecting_vertices(mesh)
implicit none 

!variables - inout 
type(hex_mesh), target :: mesh

!variables - local 
integer(in64) :: ii,aa
integer(in64) :: etgt 

!build v2e
call mesh%get_v2e()

!tag vertices
mesh%nvertex_surfint = 0 
do ii=1,mesh%nvertex
    mesh%vertex(ii)%flag = .false.
    if (allocated(mesh%vertex(ii)%ivdata)) then 
        do aa=1,mesh%max_valence
            etgt = mesh%v2e(ii,aa) 
            if (etgt > 0) then 
                ! if ((mesh%edge(etgt)%cell1 > 0) .AND. (mesh%edge(etgt)%cell2 > 0)) then !exclude far field tagged edges 
                if ((mesh%edge(etgt)%cell1 .NE. -1) .AND. (mesh%edge(etgt)%cell2 .NE. -1)) then !include far field tagged edges 
                    mesh%vertex(ii)%flag = .true.
                    mesh%nvertex_surfint = mesh%nvertex_surfint + 1
                    exit
                end if 
            end if
        end do 
    end if 
end do 
return 
end subroutine flag_surface_intersecting_vertices


!set default options =========================
subroutine set_default_options(options)
implicit none 

!variables - inout 
type(hex_options) :: options 

!set default console display to true 
options%cdisplay = .true.

!set default vtu export
options%export_vtu = .false.

!set the default mesh path 
options%meshpath = ''
options%meshname = 'grid'

!set default options 
options%mesh_treetype = 'triangle'
options%geompath = ''
options%geomname = 'geometry.fv'
options%optionspath = ''
options%optionsname = 'hex_options'
options%mesh_in_out = 'external'
options%tree_mesh_relation = 'primal'
options%mode = 'mesh'
options%ncell_max = 1000000
options%farfield_r = 20.0d0 
options%nrefine = 8
options%nflood_coarse = 4
options%nflood_mid = 10
options%nflood_fine = 8
options%nsmooth_interlayer = 1
options%nsmooth_farfield = 2
options%simplify_surfaces = .false.
options%allow_postprocess = .true.
options%nclipplane = 0 
options%nrefzone = 0 
options%nbczone = 0 

!set default constants 
options%edgelength_min = 1e-12_dp
options%cellvol_min = 1e-12_dp
options%vtx_proximity_tolerance = 4.0*atan(1.0d0)*(1e-12_dp)

!gradient projection options 
options%gradientpath = ''
options%gradientname = 'gradient'

!set internal parameters
options%geom_cell_link_method = 'sflood'
! options%geom_cell_link_method = 'dot_dir' !testing
! options%geom_cell_link_method = 'auto' !testing
options%gradient_method = 'mesh_derivative'
! options%gradient_method = 'interpolate_linear' !testing

!set flags
options%options_from_commandline = .false.
options%geom_from_commandline = .false.
return 
end subroutine set_default_options


!check geometry for self intersections function ===========================
function is_self_intersecting(geometry,options) result(is_selfintersecting)
implicit none 

!variables - inout
logical :: is_selfintersecting
type(halfedge) :: geometry
type(hex_options) :: options 

!variables - local 
integer(in64) :: ii,jj,nn
integer(in64) :: ntgtnode,nedgecheck
real(dp) :: gpad,epad,xmin,xmax,ymin,ymax
real(dp) :: v1(2),v2(2),vt1(2),vt2(2)
type(gkdtree), target :: kdtree
type(edge) :: edge_check(geometry%nedge)
type(node), dimension(:), allocatable :: tgt_nodes

!set the global padding factor
gpad = 1e-8_dp

!initialise intersection state
is_selfintersecting = .false.

!build kdtree on the geometry
call kdtree%build_tree(geometry,20_in64,4_in64,options%cdisplay)
allocate(tgt_nodes(kdtree%nnode))

!flag the geometry external halfedges 
call flag_external_halfedges(geometry)

!check for self intersections with all edges apart from each edges next and previous 
do ii=1,geometry%nedge

    !skip internal edges
    if (.NOT.geometry%edge(ii)%flag) then 
        cycle
    end if 

    !edge ends
    v1(:) = geometry%edge(ii)%origin%coordinate(1:2)
    v2(:) = geometry%edge(ii)%opposite%origin%coordinate(1:2)

    !Set padding size
    epad = norm2(v2 - v1)

    !intersection bounding box
    xmin = min(v1(1),v2(1)) - epad*gpad
    xmax = max(v1(1),v2(1)) + epad*gpad
    ymin = min(v1(2),v2(2)) - epad*gpad
    ymax = max(v1(2),v2(2)) + epad*gpad
    call kdtree%nodes_overlapping_region_2d(xmin,xmax,ymin,ymax,ntgtnode,tgt_nodes)

    !build a list of geometry edges to check 
    nedgecheck = 0 
    do nn=1,ntgtnode
        do jj=1,tgt_nodes(nn)%nitem
            if (.NOT.tgt_nodes(nn)%edge(jj)%edge%flag) then 
                cycle
            end if 
            if (tgt_nodes(nn)%edge(jj)%edge%index == geometry%edge(ii)%next%index) then 
                cycle
            elseif (tgt_nodes(nn)%edge(jj)%edge%index == geometry%edge(ii)%previous%index) then 
                cycle
            elseif (tgt_nodes(nn)%edge(jj)%edge%index == geometry%edge(ii)%index) then 
                cycle
            end if 
            nedgecheck = nedgecheck + 1
            edge_check(nedgecheck) = tgt_nodes(nn)%edge(jj)%edge
        end do 
    end do 

    !check for intersections in the selected edges 
    do jj=1,nedgecheck

        !get the vertices on this edge 
        vt1(:) = edge_check(jj)%origin%coordinate(1:2)
        vt2(:) = edge_check(jj)%opposite%origin%coordinate(1:2)

        !check 
        if (line_line_intersect_internal_bool(v1,v2,vt1,vt2)) then 
            is_selfintersecting = .true.
            exit 
        end if 
    end do 

    !exit if self intersection found 
    if (is_selfintersecting) then 
        exit 
    end if 
end do 
return 
end function is_self_intersecting

end module hex2d