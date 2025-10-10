!hex io module
!max wood
!version : 0.0.3
!updated : 23-03-25

!module 
module hex_io
use io_utilities
use hex_data_methods

!routines 
contains 

!read command arguments subroutine ===========================
subroutine get_command_arguments(options)
implicit none

!variables - import
type(hex_options) :: options 

!variables - local 
integer(in32) :: ii
integer(in32) :: nargs
integer(in64) :: slashpos
character(len=:), allocatable :: argcurr

!check and process supplied command arguments 
nargs = command_argument_count()

!set mode
if (nargs == 0) then 
    options%mode = 'mesh'
    write(*,'(A)') '** defaulting to mesh generation mode'
    ! write(*,'(A)') '** at least one argument [mode] must be supplied,& 
    ! & optionaly followed by:'
    ! print *, '-o [options file with path]'
    ! print *, '-g [geometry file with path]'
    ! stop
else

    !get first argument 
    argcurr = get_command_argument_n_str(1)

    !process argument
    if (argcurr(1:1) .NE. '-') then !set mode
        
        !get operation mode 
        options%mode = get_command_argument_n_str(1)

        !error if incorrect mode specified 
        if (options%mode == 'check') then 
            !do nothing
        elseif (options%mode == 'mesh') then 
            !do nothing
        elseif (options%mode == 'project') then 
            !do nothing
        else
            write(*,'(A,A)') '** unknown mode option requested : ',options%mode
            stop
        end if 
    else !default mode
        options%mode = 'mesh'
        write(*,'(A)') '** defaulting to mesh generation mode'
    end if 
end if 

!scan for additional arguments 
do ii=1,nargs-1

    !current argument 
    argcurr = get_command_argument_n_str(ii)

    !If option tag
    if (argcurr == '-o') then !next is options filename with path 

        !read options filepath and extract optionspath
        argcurr = get_command_argument_n_str(ii+1)
        call str_splitfslash(options%optionspath,options%optionsname,slashpos,argcurr)
        if (slashpos == -1) then 
            options%optionspath = ''
            options%optionsname = trim(argcurr)
        end if 

        !set flag 
        options%options_from_commandline = .true.
    elseif (argcurr == '-g') then !next is geometry filename with path

        !read geometry filepath and extract geompath 
        argcurr = get_command_argument_n_str(ii+1)
        call str_splitfslash(options%geompath,options%geomname,slashpos,argcurr)
        if (slashpos == -1) then 
            options%geompath = ''
            options%geomname = trim(argcurr)
        end if 

        !set flag 
        options%geom_from_commandline = .true.
    end if 
end do 
return 
end subroutine get_command_arguments


!read input options =========================
subroutine read_hex_options(options,filename)
implicit none 

!variables - inout
type(hex_options) :: options 
character(*) :: filename

!variables - local 
integer(in32) :: iostatus
integer(in64) :: slashpos,index
character(len=1000) :: rtemp 
character(len=:), allocatable :: str_temp
type(character_array), dimension(:), allocatable :: rline 

!check if file exists 
if (.NOT. file_exists(filename)) then 
    write(*,'(A)') '** cannot locate options file: '//trim(filename)
    stop 
end if 

!open options file
open(11,file=filename) 

!read options -----------
!set console display
call set_log_opt(options%cdisplay,11,'console_display')

!set geometry filepaths
if (.NOT. options%geom_from_commandline) then 
    str_temp = scan_opt_str(11,'geompath')
    if (str_temp .NE. 'na') then 
        call set_str_opt(str_temp,11,'geompath')
        call str_splitfslash(options%geompath,options%geomname,slashpos,str_temp)
        if (slashpos == -1) then 
            options%geompath = ''
            options%geomname = str_temp
        end if 
    end if 
end if 

!set mesh filepaths 
str_temp = scan_opt_str(11,'meshpath')
if (str_temp .NE. 'na') then 
    call set_str_opt(str_temp,11,'meshpath')
    call str_splitfslash(options%meshpath,options%meshname,slashpos,str_temp)
    if (slashpos == -1) then 
        options%meshpath = ''
        options%meshname = str_temp
    end if 
end if 

!set gradient filepaths 
str_temp = scan_opt_str(11,'gradientpath')
if (str_temp .NE. 'na') then 
    call set_str_opt(str_temp,11,'gradientpath')
    call str_splitfslash(options%gradientpath,options%gradientname,slashpos,str_temp)
    if (slashpos == -1) then 
        options%gradientpath = ''
        options%gradientname = str_temp
    end if
end if 

!set mesh structure options
call set_str_opt(options%mesh_in_out,11,'mesh_in_out')
call set_str_opt(options%mesh_treetype,11,'tree_type')
call set_str_opt(options%tree_mesh_relation,11,'tree_mesh_relation')
call set_int_opt(options%nrefine,11,'tree_nrefine')
call set_int_opt(options%nflood_coarse,11,'tree_nflood_coarse')
call set_int_opt(options%nflood_mid,11,'tree_nflood_mid')
call set_int_opt(options%nflood_fine,11,'tree_nflood_fine')

!set mesh property options 
call set_real_opt(options%edgelength_min,11,'mesh_edgelength_min')
call set_real_opt(options%cellvol_min,11,'mesh_cellvol_min')
call set_int_opt(options%ncell_max,11,'mesh_ncell_max')
call set_real_opt(options%farfield_r,11,'mesh_farfield_r')

!set postprocessing options 
call set_log_opt(options%allow_postprocess,11,'postprocess_allow_postprocess')
call set_int_opt(options%nsmooth_interlayer,11,'postprocess_nsmooth_interlayer')
call set_int_opt(options%nsmooth_farfield,11,'postprocess_nsmooth_farfield')

!read refinement zones 
index = 0 
iostatus = 0 
options%nrefzone = 0 
rewind(11)
do while (iostatus == 0)
    read(11,'(A)',iostat=iostatus) rtemp
    if (iostatus == -1) then 
        exit
    end if 
    if (rtemp(1:5) == 'rzone') then 
        options%nrefzone = options%nrefzone + 1
    end if 
end do 
allocate(options%refinement_zones(options%nrefzone))
iostatus = 0 
rewind(11)
do while (iostatus == 0)
    read(11,'(A)',iostat=iostatus) rtemp
    if (iostatus == -1) then 
        exit
    end if 
    if (rtemp(1:5) == 'rzone') then 
        index = index + 1
        if (rtemp(7:11) == 'point') then !read point zone
            rline = split(rtemp,' ')
            options%refinement_zones(index)%type = rline(2)%string
            options%refinement_zones(index)%xmid = str2real(rline(3)%string)
            options%refinement_zones(index)%ymid = str2real(rline(4)%string)
            options%refinement_zones(index)%zmid = str2real(rline(5)%string)
            options%refinement_zones(index)%radius = str2real(rline(6)%string)
            options%refinement_zones(index)%rlevel = str2int(rline(7)%string)
        elseif (rtemp(7:10) == 'line') then !read line zone 
            rline = split(rtemp,' ')
            options%refinement_zones(index)%type = rline(2)%string
            options%refinement_zones(index)%xmin = str2real(rline(3)%string)
            options%refinement_zones(index)%ymin = str2real(rline(4)%string)
            options%refinement_zones(index)%zmin = str2real(rline(5)%string)
            options%refinement_zones(index)%xmax = str2real(rline(6)%string)
            options%refinement_zones(index)%ymax = str2real(rline(7)%string)
            options%refinement_zones(index)%zmax = str2real(rline(8)%string)
            options%refinement_zones(index)%radius = str2real(rline(9)%string)
            options%refinement_zones(index)%rlevel = str2int(rline(10)%string)
            options%refinement_zones(index)%v1(1) = options%refinement_zones(index)%xmin
            options%refinement_zones(index)%v1(2) = options%refinement_zones(index)%ymin
            options%refinement_zones(index)%v1(3) = options%refinement_zones(index)%zmin
            options%refinement_zones(index)%v2(1) = options%refinement_zones(index)%xmax
            options%refinement_zones(index)%v2(2) = options%refinement_zones(index)%ymax
            options%refinement_zones(index)%v2(3) = options%refinement_zones(index)%zmax
        elseif (rtemp(7:10) == 'quad') then !read quad zone 
            rline = split(rtemp,' ')
            options%refinement_zones(index)%type = rline(2)%string
            options%refinement_zones(index)%xmin = str2real(rline(3)%string)
            options%refinement_zones(index)%ymin = str2real(rline(4)%string)
            options%refinement_zones(index)%zmin = str2real(rline(5)%string)
            options%refinement_zones(index)%xmax = str2real(rline(6)%string)
            options%refinement_zones(index)%ymax = str2real(rline(7)%string)
            options%refinement_zones(index)%zmax = str2real(rline(8)%string)
            options%refinement_zones(index)%rlevel = str2int(rline(9)%string)
            options%refinement_zones(index)%v1(1) = options%refinement_zones(index)%xmin
            options%refinement_zones(index)%v1(2) = options%refinement_zones(index)%ymin
            options%refinement_zones(index)%v1(3) = options%refinement_zones(index)%zmin
            options%refinement_zones(index)%v2(1) = options%refinement_zones(index)%xmax
            options%refinement_zones(index)%v2(2) = options%refinement_zones(index)%ymax
            options%refinement_zones(index)%v2(3) = options%refinement_zones(index)%zmax
        end if 
    end if 
end do 

!read clipping planes 
index = 0 
iostatus = 0 
options%nclipplane = 0 
rewind(11)
do while (iostatus == 0)
    read(11,'(A)',iostat=iostatus) rtemp
    if (iostatus == -1) then 
        exit
    end if 
    if (rtemp(1:6) == 'cplane') then 
        options%nclipplane = options%nclipplane + 1
    end if 
end do 
allocate(options%clip_planes(options%nclipplane))
iostatus = 0 
rewind(11)
do while (iostatus == 0)
    read(11,'(A)',iostat=iostatus) rtemp
    if (iostatus == -1) then 
        exit
    end if 
    if (rtemp(1:6) == 'cplane') then 
        index = index + 1
        rline = split(rtemp,' ')
        options%clip_planes(index)%v1(1) = str2real(rline(2)%string)
        options%clip_planes(index)%v1(2) = str2real(rline(3)%string)
        options%clip_planes(index)%v1(3) = str2real(rline(4)%string)
        options%clip_planes(index)%v2(1) = str2real(rline(5)%string)
        options%clip_planes(index)%v2(2) = str2real(rline(6)%string)
        options%clip_planes(index)%v2(3) = str2real(rline(7)%string)
    end if 
end do 

!read boundary contition tag zones
index = 0 
iostatus = 0 
options%nbczone = 0 
rewind(11)
do while (iostatus == 0)
    read(11,'(A)',iostat=iostatus) rtemp
    if (iostatus == -1) then 
        exit
    end if 
    if (rtemp(1:6) == 'bczone') then 
        options%nbczone = options%nbczone + 1
    end if 
end do 
allocate(options%boundarycondition_zones(options%nbczone))
iostatus = 0 
rewind(11)
do while (iostatus == 0)
    read(11,'(A)',iostat=iostatus) rtemp
    if (iostatus == -1) then 
        exit
    end if 
    if (rtemp(1:6) == 'bczone') then 
        index = index + 1
        rline = split(rtemp,' ')
        options%boundarycondition_zones(index)%xmin = str2real(rline(2)%string)
        options%boundarycondition_zones(index)%ymin = str2real(rline(3)%string)
        options%boundarycondition_zones(index)%zmin = str2real(rline(4)%string)
        options%boundarycondition_zones(index)%xmax = str2real(rline(5)%string)
        options%boundarycondition_zones(index)%ymax = str2real(rline(6)%string)
        options%boundarycondition_zones(index)%zmax = str2real(rline(7)%string)
        options%boundarycondition_zones(index)%bctag = str2int(rline(8)%string)
    end if 
end do 

!read boundary condition tag removal zones 
index = 0 
iostatus = 0 
options%nbcremzone = 0 
rewind(11)
do while (iostatus == 0)
    read(11,'(A)',iostat=iostatus) rtemp
    if (iostatus == -1) then 
        exit
    end if 
    if (rtemp(1:8) == 'bcremove') then 
        options%nbcremzone = options%nbcremzone + 1
    end if 
end do 
allocate(options%boundarycondition_remove_zones(options%nbcremzone))
iostatus = 0 
rewind(11)
do while (iostatus == 0)
    read(11,'(A)',iostat=iostatus) rtemp
    if (iostatus == -1) then 
        exit
    end if 
    if (rtemp(1:8) == 'bcremove') then 
        index = index + 1
        rline = split(rtemp,' ')
        options%boundarycondition_remove_zones(index) = str2int(rline(2)%string)
    end if 
end do 
!read options -----------

!close file 
close(11)
return 
end subroutine read_hex_options


!hex_cell mesh format export subroutine 2d =========================
subroutine write_hex_cell_mesh_2d(mesh,filename)
implicit none 

!variables - import
character(*) :: filename
type(hex_mesh), target :: mesh

!variables - local 
logical :: invalid_cell
integer(in64) :: ii,jj
integer(in64) :: etgt
integer(in64) :: loop(mesh%nedge)

! !get mesh cells 
! call mesh%index_vertices()
! call mesh%index_edges()
! call mesh%index_cells()
! call mesh%get_cell_edges()

!write mesh to file  
invalid_cell = .false.
open(11,file=filename) !mesh file
write(11,'(A,I0)') 'ncell = ',mesh%ncell 
do ii=1,mesh%ncell

    !get ordered loop of edges for this cell 
    loop(1:mesh%cell(ii)%nedge) = mesh%cell(ii)%get_edge_loop(mesh)

    !debug 
    ! print *, 'cell -> ', mesh%cell(ii)%index,' || ',mesh%cell(ii)%nedge
    ! print *, mesh%cell(ii)%nedge ,' -> ',loop(1:mesh%cell(ii)%nedge)

    !write index and number of edges in this cell 
    write(11,'(I0,A,I0)') mesh%cell(ii)%index,' ',mesh%cell(ii)%nedge 

    !write cell edges (v1 v2 adjacent_cell)
    do jj=1,mesh%cell(ii)%nedge 
        etgt = loop(jj)
        if (etgt == 0) then 
            write(*,'(A,I0,A)') '    ** zero index in vertex loop of cell ',ii,', this cell is likely bisected'
            invalid_cell = .true.
            exit
        end if 
        if (mesh%edge(etgt)%cell1 == ii) then 
            write(11,'(I0,A,I0,A,I0)') mesh%edge(etgt)%vertex1%index,' ',mesh%edge(etgt)%vertex2%index,' ',mesh%edge(etgt)%cell2
        else
            write(11,'(I0,A,I0,A,I0)') mesh%edge(etgt)%vertex2%index,' ',mesh%edge(etgt)%vertex1%index,' ',mesh%edge(etgt)%cell1
        end if 
    end do 
    if (invalid_cell) then 
        exit 
    end if 
end do 
if (.NOT. invalid_cell) then 
    write(11,'(A,I0)') 'nvertex = ',mesh%nvertex 
    do ii=1,mesh%nvertex 
        write(11,'(E17.10,A,E17.10)') mesh%vertex(ii)%coordinate(1),' ',mesh%vertex(ii)%coordinate(2)
    end do 
end if 
close(11)
if (invalid_cell) then 
    write(*,'(A)') '    ** invalid mesh not written to file'
end if
return 
end subroutine write_hex_cell_mesh_2d


!hex mesh format export subroutine 2d =========================
subroutine write_hex_mesh_2d(mesh,filename)
implicit none 

!variables - import
character(*) :: filename
type(hex_mesh), target :: mesh

!variables - local 
integer(in64) :: ii

!write mesh to file  
open(11,file=filename) !mesh file
write(11,'(I0,A,I0,A,I0)') mesh%ncell,' ',mesh%nedge,' ',mesh%nvertex 
do ii=1,mesh%nedge !edges
    write(11,'(I0,A,I0,A,I0,A,I0)') mesh%edge(ii)%vertex2%index,' ',mesh%edge(ii)%vertex1%index,' ',&
                                    mesh%edge(ii)%cell2,' ',mesh%edge(ii)%cell1 
end do
do ii=1,mesh%nvertex !vertices
    write(11,'(I0,A,E17.10,A,E17.10)') ii,' ',mesh%vertex(ii)%coordinate(1),' ',mesh%vertex(ii)%coordinate(2)
end do
write(11,'(I0)') mesh%nvertex_surfint
do ii=1,mesh%nvertex !surface link data 
    if (mesh%vertex(ii)%flag) then 
        write(11,'(I0,A,E17.10,A,I0,A,I0)') ii,' ',mesh%vertex(ii)%rdata,' ',&
                                            mesh%vertex(ii)%ivdata(1),' ',mesh%vertex(ii)%ivdata(2) !mesh vertex index | fraction | geometry vertex 1 | geometry vertex 2
    end if 
end do
close(11)
return 
end subroutine write_hex_mesh_2d


!import mesh 2d =========================
subroutine read_hex_mesh_2d(mesh,filename)
implicit none 

!variables - import
character(*) :: filename
type(hex_mesh), target :: mesh

!variables - local 
integer(in64) :: ii
integer(in64) :: vidx,ev1,ev2,nvsurf
real(dp) :: ef

!check if file exists 
if (.NOT. file_exists(filename)) then 
    write(*,'(A)') '** cannot locate mesh file: '//trim(filename)
    stop 
end if 

!open file 
open(11,file=filename)

!read item quantities 
read(11,*) mesh%ncell,mesh%nedge,mesh%nvertex

!allocate
allocate(mesh%edge(mesh%nedge))
allocate(mesh%vertex(mesh%nvertex))

!read edges 
do ii=1,mesh%nedge
    read(11,*) ev2,ev1,mesh%edge(ii)%cell2,mesh%edge(ii)%cell1
    mesh%edge(ii)%vertex1 => mesh%vertex(ev1)
    mesh%edge(ii)%vertex2 => mesh%vertex(ev2)
end do 

!read vertices and initialise flags and indecies
do ii=1,mesh%nvertex
    read(11,*) vidx,mesh%vertex(ii)%coordinate(1:2)
    mesh%vertex(ii)%coordinate(3) = 0.0d0 
    mesh%vertex(ii)%index = ii
    mesh%vertex(ii)%flag = .false.
    mesh%vertex(ii)%rdata = 0.0d0 
end do 

!read surface links and set vertex flags
read(11,*) nvsurf
do ii=1,nvsurf
    read(11,*) vidx,ef,ev1,ev2
    allocate(mesh%vertex(vidx)%ivdata(2))
    mesh%vertex(vidx)%rdata = ef
    mesh%vertex(vidx)%ivdata(1) = ev1
    mesh%vertex(vidx)%ivdata(2) = ev2
    mesh%vertex(vidx)%flag = .true.
end do 

!close
close(11)
return 
end subroutine read_hex_mesh_2d


!import gradient =========================
subroutine read_gradient(mesh,filename)
implicit none 

!variables - import
character(*) :: filename
type(hex_mesh), target :: mesh

!variables - local 
integer(in64) :: ii

!check if file exists 
if (.NOT. file_exists(filename)) then 
    write(*,'(A)') '** cannot locate gradient file: '//trim(filename)
    stop 
end if 

!read file 
open(11,file=filename)
do ii=1,mesh%nvertex
    read(11,*) mesh%vertex(ii)%gradient(1:2)
    mesh%vertex(ii)%gradient(3) = 0.0d0 
end do 
close(11)
return 
end subroutine read_gradient


!export gradient 2d =========================
subroutine write_gradient_2d(surface_gradient,filename)
implicit none 

!variables - import
character(*) :: filename
real(dp), dimension(:,:), allocatable :: surface_gradient

!variables - local
integer(in64) :: ii

!write
open(11,file=filename)
do ii=1,size(surface_gradient,dim=1)
    write(11,'(E17.10,A,E17.10)') surface_gradient(ii,1),' ',surface_gradient(ii,2)
end do 
close(11)
return 
end subroutine write_gradient_2d

end module hex_io