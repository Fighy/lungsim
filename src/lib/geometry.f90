
module geometry
!*Brief Description:* This module handles all geometry read/write/generation.
!
!*LICENSE:*
!
!
!
!*Full Description:*
!
!This module handles all geometry read/write/generation.
  use arrays
  use diagnostics
  use indices
  use mesh_utilities
  use other_consts ! currently has pi
  use precision ! sets dp for precision

  implicit none

  !Module parameters

  !Module types

  !Module variables

  !Interfaces
  private
  public add_mesh
  public add_matching_mesh
  public append_units
  public define_1d_elements
  public define_elem_geometry_2d
  public define_mesh_geometry_test
  public define_node_geometry
  public define_node_geometry_2d
  public define_data_geometry
  public group_elem_parent_term
  public define_rad_from_file
  public define_rad_from_geom
  public element_connectivity_1d
  public element_connectivity_2d
  public inlist
  public evaluate_ordering
  public get_final_real
  public get_local_node_f
  public make_data_grid
  public reallocate_node_elem_arrays
  public set_initial_volume
  public list_mesh_statistics
  public triangles_from_surface
  public volume_of_mesh
  public get_final_integer
  public get_four_nodes

contains
!
!!!#####################################################################
!
!*allocate_node_arrays:* allocate memory for arrays associated with 1D trees
  subroutine allocate_node_arrays(num_nodes)

    integer,intent(in) :: num_nodes
    
    if(.not.allocated(nodes)) allocate (nodes(num_nodes))
    if(.not.allocated(node_xyz)) allocate (node_xyz(3,num_nodes))
    if(.not.allocated(node_field)) allocate (node_field(num_nj,num_nodes))
    if(.not.allocated(elems_at_node)) allocate(elems_at_node(num_nodes,0:3))
    nodes = 0 !initialise node index values
    node_xyz = 0.0_dp !initialise
    node_field = 0.0_dp !initialise

  end subroutine allocate_node_arrays

!###################################################################################
!
!*add_mesh:* Reads in an ipmesh file and adds this mesh to the terminal branches of an existing tree geometry
  subroutine add_mesh(AIRWAY_MESHFILE)
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_ADD_MESH" :: ADD_MESH

    character(len=MAX_FILENAME_LEN), intent(in) :: AIRWAY_MESHFILE
    ! Local parameters
    character(len=100) :: buffer
    integer, parameter :: fh = 15
    integer :: ios
    integer :: line
    integer :: i,ibeg,iend,i_ss_end,j,ne,ne0,ne_global,ne_parent,ne_start, &
         ngen_parent,np,np0,np_global,&
         num_elems_new,num_elems_to_add,num_nodes_new,nunit,nlabel
    integer,dimension(1000) :: element_temp,generation, &
         parent_element,symmetry_temp
    real(dp),dimension(1000) :: length,radius,a_A
    character(len=60) :: sub_name

    sub_name = 'add_mesh'
    call enter_exit(sub_name,1)

    ios = 0
    line = 0
    open(fh, file=AIRWAY_MESHFILE)

    ! ios is negative if an end of record condition is encountered or if
    ! an endfile condition was detected.  It is positive if an error was
    ! detected.  ios is zero otherwise.
    ios=0
    line=0
    i=0 ! count the number of elements read in

    do while (ios == 0)
       read(fh, '(A)', iostat=ios) buffer
       ! line contains: element, parent element, generation,
       !                symmetry, length, outer radius, a/A ratio
       ! note that a/A ratio is always 1 for the conducting airways
       if (ios == 0) then
          line = line + 1
          i=i+1
          i_ss_end = len(buffer)

          do nlabel = 1,7
             ibeg = index(buffer," ") + 1 !get location of first integer beyond ws in string
             buffer = adjustl(buffer(ibeg:i_ss_end)) ! get info beyond ws, remove leading ws
             iend = index(buffer," ") !get location of next ws in sub-string
             select case (nlabel)
             case (1)
                read (buffer(1:iend-1), *, iostat=ios) element_temp(i) ! not used??
             case(2)
                read (buffer(1:iend-1), *, iostat=ios) parent_element(i)
             case(3)
                read (buffer(1:iend-1), *, iostat=ios) generation(i)
             case(4)
                read (buffer(1:iend-1), *, iostat=ios) symmetry_temp(i)
             case(5)
                read (buffer(1:iend-1), *, iostat=ios) length(i)
             case(6)
                read (buffer(1:iend-1), *, iostat=ios) radius(i)
             case(7)
                read (buffer(1:iend-1), *, iostat=ios) a_A(i)
             end select
          enddo
       endif
    enddo
    close(fh)

    num_elems_to_add = i

!!! increase the size of node and element arrays to accommodate the additional elements
    ! the number of nodes after adding mesh will be:
    num_nodes_new = num_nodes + num_units*num_elems_to_add
    ! the number of elems after adding mesh will be:
    num_elems_new = num_elems + num_units*num_elems_to_add
    call reallocate_node_elem_arrays(num_elems_new,num_nodes_new)

    ne = num_elems ! the starting local element number
    ne_global = elems(ne) ! assumes this is the highest element number (!!!)
    np = num_nodes ! the starting local node number
    np_global = nodes(np) ! assumes this is the highest node number (!!!)

    do nunit = 1,num_units ! for all terminal branches, append the mesh

       ne_parent = units(nunit) ! local element number of terminal, to append to
       ngen_parent = elem_ordrs(1,ne_parent)
       ne_start = ne !starting element number for the unit

       do i=1,num_elems_to_add

          if(parent_element(i).eq.0)then
             ne_parent = units(nunit)
          else
             ne_parent = ne_start+parent_element(i)
          endif

          ne0 = ne_parent
          np0 = elem_nodes(2,ne0)

          ne_global = ne_global + 1 ! new global element number
          ne = ne + 1 ! new local element number
          np_global = np_global + 1 !new global node number
          np = np + 1 ! new local node number

          nodes(np) = np_global
          elems(ne) = ne_global

          elem_nodes(1,ne) = np0
          elem_nodes(2,ne) = np

          elem_ordrs(1,ne) = ngen_parent + generation(i)
          elem_symmetry(ne) = symmetry_temp(i)+1 ! uses 0/1 in file; 1/2 in code

          ! record the element connectivity
          elem_cnct(-1,0,ne) = 1 ! one parent branch
          elem_cnct(-1,1,ne) = ne0 ! store parent element
          elem_cnct(1,0,ne0) = elem_cnct(1,0,ne0) + 1
          elem_cnct(1,elem_cnct(1,0,ne0),ne0) = ne

          ! record the direction and location of the branch
          do j=1,3
             elem_direction(j,ne) = elem_direction(j,ne0)
             node_xyz(j,np) = node_xyz(j,np0) + &
                  elem_direction(j,ne)*length(i)
          enddo !j

          elem_field(ne_length,ne) = length(i)
          elem_field(ne_radius,ne) = radius(i)
          elem_field(ne_a_A,ne) = a_A(i)
          elem_field(ne_vol,ne) = PI*radius(i)**2*length(i)

       enddo !i
    enddo !nunit

    num_nodes = np
    num_elems = ne

    call element_connectivity_1d
    call evaluate_ordering ! calculate new ordering of tree

    call enter_exit(sub_name,2)

  end subroutine add_mesh

!
!###################################################################################
!
!*add_matching_mesh:*
  subroutine add_matching_mesh()
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_ADD_MATCHING_MESH" :: ADD_MATCHING_MESH

    !Parameters to become inputs
    real(dp) :: offset(3)
    logical :: REVERSE=.TRUE.
    character(len=60) :: mesh_type='terminal'
    !local variables
    integer :: num_nodes_new,num_elems_new,ne,ne_global,np,np_global,np0,nonode,np_m
    integer :: nj,ne_m,noelem,ne0,n,nindex,ne1,noelem0,nu,cap_conns,cap_term,np1,np2
    integer, allocatable :: np_map(:)
    character(len=60) :: sub_name

    sub_name = 'add_matching_mesh'
    call enter_exit(sub_name,1)
    !Ultimately offset should be an input argument
    offset(1)=0.0_dp
    offset(2)=1e-6_dp
    offset(3)=0.0_dp


    allocate(np_map(num_nodes))
!!! increase the size of node and element arrays to accommodate the additional elements
    ! the number of nodes after adding mesh will be:
    num_nodes_new = 2*num_nodes
    ! the number of elems after adding mesh will be:
    if(mesh_type.eq.'basic')then
      num_elems_new = 2*num_elems
    elseif(mesh_type.eq.'terminal')then
      num_elems_new = 2*num_elems + num_units
    elseif(mesh_type.eq.'ladder')then
        print *, "NOT YET IMPLEMENTED"
        call exit(0)
    endif
    call reallocate_node_elem_arrays(num_elems_new,num_nodes_new)
    noelem0=0
    ne0 = num_elems ! the starting local element number
    ne_global = elems(ne0) ! assumes this is the highest element number (!!!)
    np0 = num_nodes ! the starting local node number
    np_global = nodes(np0) ! assumes this is the highest node number (!!!)

    do nonode=1,num_nodes
       np=np_global+nonode
       np_m=nodes(nonode)
       np_map(np_m)=np !maps new to old node numbering
       nodes(np0+nonode)=np
       do nj=1,3
         node_xyz(nj,np)=node_xyz(nj,np_m)+offset(nj)
       enddo
       elems_at_node(np,0)=0 !initialise
     !Doesnt map versions, would be added here
    enddo

    do noelem=1,num_elems
        ne=ne_global+noelem
        elem_field(ne_group,ne)=2.0_dp!VEIN
        ne_m=elems(noelem)
        elem_field(ne_group,ne_m)=0.0_dp!ARTERY
        elems(ne0+noelem)=ne
        if(.NOT.REVERSE)then
          elem_nodes(1,ne)=np_map(elem_nodes(1,ne_m))
          elem_nodes(2,ne)=np_map(elem_nodes(2,ne_m))
          elem_cnct(1,0,ne)=elem_cnct(1,0,ne_m)!The numberdownstream are the number downstream
          elem_cnct(-1,0,ne)=elem_cnct(-1,0,ne_m)
          do n=1,elem_cnct(1,0,ne)
            elem_cnct(1,n,ne)=elem_cnct(1,n,ne_m)+ne0
          enddo
          do n=1,elem_cnct(-1,0,ne)
            elem_cnct(-1,n,ne)=elem_cnct(-1,n,ne_m)+ne0
          enddo
        else
          elem_nodes(1,ne)=np_map(elem_nodes(2,ne_m))
          elem_nodes(2,ne)=np_map(elem_nodes(1,ne_m))
          elem_cnct(-1,0,ne)=elem_cnct(1,0,ne_m) !The number upstream are the number downstream
          elem_cnct(1,0,ne)=elem_cnct(-1,0,ne_m)!The number downstream are the number upstream
          do n=1,elem_cnct(1,0,ne)
            elem_cnct(1,n,ne)=elem_cnct(-1,n,ne_m)+ne0
          enddo
          do n=1,elem_cnct(-1,0,ne)
            elem_cnct(-1,n,ne)=elem_cnct(1,n,ne_m)+ne0
          enddo
        endif
        !if worrying about regions and versions do it here
        elems_at_node(elem_nodes(1,ne),0)=elems_at_node(elem_nodes(1,ne),0)+1
        elems_at_node(elem_nodes(1,ne),elems_at_node(elem_nodes(1,ne),0))=ne
        elems_at_node(elem_nodes(2,ne),0)=elems_at_node(elem_nodes(2,ne),0)+1
        elems_at_node(elem_nodes(2,ne),elems_at_node(elem_nodes(2,ne),0))=ne
        nindex=no_gen
        elem_ordrs(nindex,ne)=elem_ordrs(nindex,ne_m)
        nindex=no_sord
        elem_ordrs(nindex,ne)=elem_ordrs(nindex,ne_m)
        nindex=no_hord
        elem_ordrs(nindex,ne)=elem_ordrs(nindex,ne_m)
      enddo

     !update current no of nodes and elements to determine connectivity
     np0=np !current highest node
     ne1=ne !current highest element
     noelem0=num_elems+noelem0
     if(mesh_type.eq.'ladder')then
        !To be implemented
     elseif(mesh_type.eq.'terminal')then
       cap_conns=0
       cap_term=0
       do nu=1,num_units
         ne=units(nu)
         cap_term=cap_term+1
         np1=elem_nodes(2,ne)
         np2=np_map(np1)
         noelem0=noelem0+1
         ne1=ne1+1
         elems(noelem0)=ne1
         elem_nodes(1,ne1)=np1
         elem_nodes(2,ne1)=np2
         elems_at_node(np1,0)=elems_at_node(np1,0)+1
         elems_at_node(np1,elems_at_node(np1,0))=ne1
         elems_at_node(np2,0)=elems_at_node(np2,0)+1
         elems_at_node(np2,elems_at_node(np2,0))=ne1
         elem_cnct(1,elem_cnct(1,0,ne)+1,ne)=ne1
         elem_cnct(1,0,ne)=elem_cnct(1,0,ne)+1
         elem_cnct(-1,elem_cnct(-1,0,ne+ne_global)+1,ne+ne_global)=ne1
         elem_cnct(-1,0,ne+ne_global)=elem_cnct(-1,0,ne+ne_global)+1
         elem_cnct(-1,0,ne1)=1
         elem_cnct(1,0,ne1)=1
         elem_cnct(-1,1,ne1)=ne
         elem_cnct(1,1,ne1)=ne+ne0
         nindex=no_gen
         elem_ordrs(nindex,ne1)=elem_ordrs(nindex,ne_m)
         nindex=no_sord
         elem_ordrs(nindex,ne1)=elem_ordrs(nindex,ne_m)
         nindex=no_hord
         elem_ordrs(nindex,ne1)=elem_ordrs(nindex,ne_m)
         elem_field(ne_group,ne1)=1.0_dp!connection between meshes
       enddo
        print *, 'Number of connections', cap_term
     endif
    num_nodes=num_nodes_new
    num_elems=num_elems_new
    deallocate(np_map)
    call enter_exit(sub_name,2)

  end subroutine add_matching_mesh

!
!###################################################################################
!
!*append_units:* Appends terminal units at the end of a tree structure
  subroutine append_units()
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_APPEND_UNITS" :: APPEND_UNITS

    integer :: ne,ne0,nu
    character(len=60) :: sub_name

    sub_name = 'append_units'
    call enter_exit(sub_name,1)

    num_units = 0
    DO ne=1,num_elems
       IF(elem_cnct(1,0,ne).eq.0)THEN
          num_units=num_units+1
       ENDIF
    ENDDO

    if(allocated(units))then !increasing the array size; just overwrite
       deallocate(units)
       deallocate(unit_field)
    endif
    allocate(units(num_units))
    allocate(unit_field(num_nu,num_units))

    unit_field=0.0_dp
    units=0
    elem_units_below(1:num_elems) = 0 !initialise the number of terminal units below a branch

    nu=0
    DO ne=1,num_elems
       IF(elem_cnct(1,0,ne).eq.0)THEN
          nu=nu+1
          units(nu)=ne     !Set up units array containing terminals
          elem_units_below(ne)=1
       ENDIF
    ENDDO

    ! count the effective number of elements below each branch
    do ne=num_elems,2,-1
       ne0=elem_cnct(-1,1,ne)
       elem_units_below(ne0) = elem_units_below(ne0) &
            + elem_units_below(ne)*elem_symmetry(ne)
    enddo !ne

    call enter_exit(sub_name,2)

  end subroutine append_units

!
!###################################################################################
!
  !*define_1d_elements:* Reads in an 1D element ipelem file to define a geometry
  subroutine define_1d_elements(ELEMFILE)
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_DEFINE_1D_ELEMENTS" :: DEFINE_1D_ELEMENTS

    character(len=MAX_FILENAME_LEN), intent(in) :: ELEMFILE
    !     Local Variables
    integer :: ibeg,iend,ierror,i_ss_end,j,ne,ne_global,&
         nn,np,np1,np2,np_global
    character(LEN=132) :: ctemp1
    character(LEN=40) :: sub_string
    character(len=60) :: sub_name

    sub_name = 'define_1d_elements'
    call enter_exit(sub_name,1)

    open(10, file=ELEMFILE, status='old')

    read_number_of_elements : do
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1
       if(index(ctemp1, "elements")> 0) then
          num_elems = get_final_integer(ctemp1) !return the final integer
          exit read_number_of_elements
       endif
    end do read_number_of_elements

!!! allocate memory for element arrays
    if(allocated(elems)) deallocate(elems)
    allocate(elems(num_elems))
    if(allocated(elem_cnct)) deallocate(elem_cnct)
    allocate(elem_cnct(-1:1,0:2,0:num_elems))
    if(allocated(elem_nodes)) deallocate(elem_nodes)
    allocate(elem_nodes(2,num_elems))
    if(allocated(elem_ordrs)) deallocate(elem_ordrs)
    allocate(elem_ordrs(num_ord,num_elems))
    if(allocated(elem_symmetry)) deallocate(elem_symmetry)
    allocate(elem_symmetry(num_elems))
    if(allocated(elem_units_below)) deallocate(elem_units_below)
    allocate(elem_units_below(num_elems))
    if(allocated(elems_at_node)) deallocate(elems_at_node)
    allocate(elems_at_node(num_nodes,0:3))
    if(allocated(elem_field)) deallocate(elem_field)
    allocate(elem_field(num_ne,num_elems))
    if(allocated(elem_direction)) deallocate(elem_direction)
    allocate(elem_direction(3,num_elems))
    if(model_type.eq.'gas_mix')then
      if(allocated(expansile)) deallocate(expansile)
      allocate(expansile(num_elems))
    endif

!!! initialise element arrays
    elems=0
    elem_nodes=0
    elem_symmetry = 1
    elem_field = 0.0_dp
    if(model_type.eq.'gas_mix')expansile = .false.

    ne=0

    read_an_element : do
       !.......read element number
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1
       if(index(ctemp1, "Element")> 0) then
          ne_global = get_final_integer(ctemp1) !return the final integer
          ne=ne+1
          elems(ne)=ne_global

          read_element_nodes : do
             read(unit=10, fmt="(a)", iostat=ierror) ctemp1
             if(index(ctemp1, "global")> 0) then !found the correct line
                iend=len(ctemp1)
                ibeg=index(ctemp1,":")+1 !get location of first integer in string
                sub_string = adjustl(ctemp1(ibeg:iend)) ! get the characters beyond : remove leading blanks
                i_ss_end=len(sub_string) !get the end location of the sub-string
                ibeg=1
                do nn=1,2
                   iend=index(sub_string," ") !get location of first blank in sub-string
                   read (sub_string(ibeg:iend-1), '(i7)' ) np_global
                   call get_local_node(np_global,np) ! get local node np for global node
                   elem_nodes(nn,ne)=np ! the local node number, not global
                   sub_string = adjustl(sub_string(iend:i_ss_end)) ! get chars beyond blank, remove leading blanks
                end do
                exit read_element_nodes
             endif !index
          end do read_element_nodes
          if(ne.ge.num_elems) exit read_an_element
       endif

    end do read_an_element

    close(10)

    ! calculate the element lengths and directions
    do ne=1,num_elems
       np1=elem_nodes(1,ne)
       np2=elem_nodes(2,ne)
       elem_field(ne_length,ne) = DSQRT((node_xyz(1,np2) - &
            node_xyz(1,np1))**2 + (node_xyz(2,np2) - &
            node_xyz(2,np1))**2 + (node_xyz(3,np2) - &
            node_xyz(3,np1))**2)
       do j=1,3
          elem_direction(j,ne) = (node_xyz(j,np2) - &
               node_xyz(j,np1))/elem_field(ne_length,ne)
       enddo !j
    enddo

    call element_connectivity_1d
    call evaluate_ordering

    call enter_exit(sub_name,2)

  END subroutine define_1d_elements


!!!##################################################

  subroutine define_elem_geometry_2d(ELEMFILE,sf_option)
      !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_DEFINE_ELEM_GEOMETRY_2D" :: DEFINE_ELEM_GEOMETRY_2D

    ! Reads in 2D ipelem file.

    character(len=*) :: ELEMFILE
    character(len=4) :: sf_option
    
    !     Local Variables
    integer :: ierror,ne,ne_global,nn,np,number_of_elements
    character(len=132) :: ctemp1,readfile
    character(len=60) :: sub_name = 'define_elem_geometry_2d'
    
    call enter_exit(sub_name,1)
    
    !readfile = trim(elemfile)//'.ipelem'
    open(10, file=elemfile, status='old')
    
    read_number_of_elements : do
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1
       if(index(ctemp1, "elements")> 0) then
          number_of_elements = get_final_integer(ctemp1) !return the final integer
          exit read_number_of_elements
       endif
    end do read_number_of_elements

    num_elems_2d=number_of_elements
    if(.not.allocated(elems_2d)) allocate(elems_2d(num_elems_2d))
    if(.not.allocated(elem_nodes_2d)) allocate(elem_nodes_2d(4,num_elems_2d))
    if(.not.allocated(elem_versn_2d)) allocate(elem_versn_2d(4,num_elems_2d))
    
    ne = 0
    
    read_an_element : do 
       !.......read element number
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1
       if(index(ctemp1, "Element")> 0) then
          ne_global = get_final_integer(ctemp1) !return the final integer
          ne = ne + 1
          elems_2d(ne) = ne_global
          
          read_element_nodes : do 
             read(unit=10, fmt="(a)", iostat=ierror) ctemp1
             if(index(ctemp1, "global")> 0) then !found the correct line
                call get_four_nodes(ne,ctemp1) !number of versions for node np
                ! note that only the ne'th data of elem_nodes_2d is passed to 'get_four_nodes'
                do nn=1,4
                   np=elem_nodes_2d(nn,ne)
                   if(node_versn_2d(np).gt.1)then
                      read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !contains version# for njj=1
                      read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !contains version# for njj=1
                      read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !contains version# for njj=1
                      elem_versn_2d(nn,ne) = get_final_integer(ctemp1) !return the final integer
                   else
                      elem_versn_2d(nn,ne)= 1
                   endif !nversions
                enddo !nn
                exit read_element_nodes
             endif !index
          end do read_element_nodes
          
          if(ne.ge.number_of_elements) exit read_an_element
       endif
       
    end do read_an_element
    
    close(10)
    
    call element_connectivity_2d
    call line_segments_for_2d_mesh(sf_option)

    call enter_exit(sub_name,2)
    
  end subroutine define_elem_geometry_2d

!!!###############################################################

!*define_mesh_geometry_test:*
  subroutine define_mesh_geometry_test()
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_DEFINE_MESH_GEOMETRY_TEST" :: DEFINE_MESH_GEOMETRY_TEST

    !     Local Variables
    integer :: j,ne,np,np1,np2

    num_elems = 400
    num_nodes = num_elems + 1

!!! allocate memory
    if(.not.allocated(nodes)) allocate (nodes(num_nodes))
    if(.not.allocated(node_xyz)) allocate (node_xyz(3,num_nodes))
    if(.not.allocated(node_field)) allocate (node_field(num_nj,num_nodes))
    if(.not.allocated(elems)) allocate(elems(num_elems))
    if(.not.allocated(elem_cnct)) allocate(elem_cnct(-1:1,0:2,0:num_elems))
    if(.not.allocated(elem_nodes)) allocate(elem_nodes(2,num_elems))
    if(.not.allocated(elem_ordrs)) allocate(elem_ordrs(num_ord,num_elems))
    if(.not.allocated(elem_symmetry)) allocate(elem_symmetry(num_elems))
    if(.not.allocated(elem_units_below)) allocate(elem_units_below(num_elems))
    if(.not.allocated(elems_at_node)) allocate(elems_at_node(num_nodes,0:3))
    if(.not.allocated(elem_field)) allocate(elem_field(num_ne,num_elems))
    if(.not.allocated(elem_direction)) allocate(elem_direction(3,num_elems))
    if(.not.allocated(expansile)) allocate(expansile(num_elems))

!!! initialise array values
    nodes = 0 !initialise node index values
    node_xyz = 0.0_dp !initialise
    node_field = 0.0_dp !initialise
    elems=0
    elem_nodes=0
    elem_symmetry = 1
    elem_field = 0.0_dp
    expansile = .false.

!!! set up node arrays
    nodes(1) = 1
    do np=2,101
       nodes(np) = np
       node_xyz(3,np) = node_xyz(3,np-1) - 1.0_dp
    enddo

    np=102
    node_xyz(1,np) = node_xyz(1,np-1) - 0.5_dp
    do np=102,151
       nodes(np) = np
       node_xyz(3,np) = node_xyz(3,np-1) - 0.5_dp
       node_xyz(1,np) = node_xyz(1,np-1) - 0.5_dp
    enddo

    np=152
    nodes(np) = np
    node_xyz(1,np) = node_xyz(1,101) + 0.5_dp
    node_xyz(3,np) = node_xyz(3,101) - 0.5_dp
    do np=153,201
       nodes(np) = np
       node_xyz(3,np) = node_xyz(3,np-1) - 0.5_dp
       node_xyz(1,np) = node_xyz(1,np-1) + 0.5_dp
    enddo

    np=202
    nodes(np) = np
    node_xyz(1,np) = node_xyz(1,151) - 0.5_dp
    node_xyz(3,np) = node_xyz(3,151) - 0.5_dp
    do np=203,251
       nodes(np) = np
       node_xyz(3,np) = node_xyz(3,np-1) - 0.5_dp
       node_xyz(1,np) = node_xyz(1,np-1) - 0.5_dp
    enddo

    np=252
    nodes(np) = np
    node_xyz(1,np) = node_xyz(1,151) + 0.5_dp
    node_xyz(3,np) = node_xyz(3,151) - 0.5_dp
    do np=253,301
       nodes(np) = np
       node_xyz(3,np) = node_xyz(3,np-1) - 0.5_dp
       node_xyz(1,np) = node_xyz(1,np-1) + 0.5_dp
    enddo

    np=302
    nodes(np) = np
    node_xyz(1,np) = node_xyz(1,201) - 0.5_dp
    node_xyz(3,np) = node_xyz(3,201) - 0.5_dp
    do np=303,351
       nodes(np) = np
       node_xyz(3,np) = node_xyz(3,np-1) - 0.5_dp
       node_xyz(1,np) = node_xyz(1,np-1) - 0.5_dp
    enddo

    np=352
    nodes(np) = np
    node_xyz(1,np) = node_xyz(1,201) + 0.5_dp
    node_xyz(3,np) = node_xyz(3,201) - 0.5_dp
    do np=353,401
       nodes(np) = np
       node_xyz(3,np) = node_xyz(3,np-1) - 0.5_dp
       node_xyz(1,np) = node_xyz(1,np-1) + 0.5_dp
    enddo

!!! set up element arrays
    do ne=1,num_elems
       elems(ne) = ne
       elem_nodes(1,ne) = ne
       elem_nodes(2,ne) = ne+1
    enddo

    elem_nodes(1,151) = 101
    elem_nodes(1,201) = 151
    elem_nodes(1,251) = 151
    elem_nodes(1,301) = 201
    elem_nodes(1,351) = 201

    elem_field(ne_radius,1:100) = 10.0_dp
    elem_field(ne_radius,101:200) = 5.0_dp
    elem_field(ne_radius,201:400) = dsqrt(elem_field(ne_radius,101)**2/2.0_dp)

    ! calculate the element lengths and directions
    do ne=1,num_elems
       np1=elem_nodes(1,ne)
       np2=elem_nodes(2,ne)
       elem_field(ne_length,ne) = DSQRT((node_xyz(1,np2) - &
            node_xyz(1,np1))**2 + (node_xyz(2,np2) - &
            node_xyz(2,np1))**2 + (node_xyz(3,np2) - &
            node_xyz(3,np1))**2)
       do j=1,3
          elem_direction(j,ne) = (node_xyz(j,np2) - &
               node_xyz(j,np1))/elem_field(ne_length,ne)
       enddo !j
       ! element volume
       elem_field(ne_vol,ne) = PI * elem_field(ne_radius,ne)**2 * &
            elem_field(ne_length,ne)
       elem_field(ne_a_A,ne) = 1.0_dp ! set default for ratio a/A
    enddo

    call element_connectivity_1d
    call evaluate_ordering

  end subroutine define_mesh_geometry_test
!
!###################################################################################
!
  subroutine define_node_geometry(NODEFILE)
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_DEFINE_NODE_GEOMETRY" :: DEFINE_NODE_GEOMETRY

  !*define_node_geometry:* Reads in an ipnode file to define a tree geometry

    character(len=MAX_FILENAME_LEN), intent(in) :: NODEFILE !Input nodefile
    !     Local Variables
    integer :: i,ierror,np,np_global,num_nodes_temp,num_versions,nv,NJT=0
    character(len=300) :: ctemp1,readfile
    real(dp) :: point
    character(len=60) :: sub_name
    logical :: overwrite = .false. ! initialised

    sub_name = 'define_node_geometry'
    call enter_exit(sub_name,1)

    readfile = trim(nodefile)//'.ipnode'
    open(10, file=NODEFILE, status='old')

    if(num_nodes.gt.0) overwrite = .true.
    
    !.....read in the total number of nodes. read each line until one is found
    !.....that has the correct keyword (nodes). then return the integer that is
    !.....at the end of the line
    read_number_of_nodes : do !define a do loop name
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !read a line into ctemp1
       if(index(ctemp1, "nodes")> 0) then !keyword "nodes" is found in ctemp1
          num_nodes_temp = get_final_integer(ctemp1) !return the final integer
          exit read_number_of_nodes !exit the named do loop
       endif
    end do read_number_of_nodes

    if(.not.overwrite) call allocate_node_arrays(num_nodes_temp) ! don't allocate if just overwriting

    !.....read in the number of coordinates
    read_number_of_coords : do !define a do loop name
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !read a line into ctemp1
       if(index(ctemp1, "coordinates")> 0) then !keyword "coordinates" is found
          NJT = get_final_integer(ctemp1) !return the final integer
          exit read_number_of_coords !exit the named do loop
       endif
    end do read_number_of_coords
    
    ! note that only the first version of coordinate is currently read in   
    
    !.....read the coordinate, derivative, and version information for each node. 
    np=0
    read_a_node : do !define a do loop name
       !.......read node number
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1
       if(index(ctemp1, "Node")> 0) then
          np_global = get_final_integer(ctemp1) !get node number

          np = np+1
          nodes(np) = np_global
          !.......read coordinates
          do i=1,3 ! for the x,y,z coordinates
             num_versions=1
             read(unit=10, fmt="(a)", iostat=ierror) ctemp1
             if(index(ctemp1, "versions")> 0) then
                num_versions = get_final_integer(ctemp1)
                if(num_versions > 1)then
                   read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !temporary line
                   read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                   point = get_final_real(ctemp1)
                   do nv=2,num_versions
                      read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !temporary line
                      read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                   enddo
                else
                   read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                   point = get_final_real(ctemp1)
                endif
             else ! no prompting for versions
                point = get_final_real(ctemp1)
             endif
             node_xyz(i,np)=point
          end do !i
       endif !index
       if(np.ge.num_nodes_temp) exit read_a_node
    end do read_a_node

    if(.not.overwrite) num_nodes = num_nodes_temp
    
    close(10)

    call enter_exit(sub_name,2)

  end subroutine define_node_geometry

!!!##################################################

  subroutine define_node_geometry_2d(NODEFILE)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_DEFINE_NODE_GEOMETRY_2D" :: DEFINE_NODE_GEOMETRY_2D
  
  !*define_node_geometry_2d:* Reads in an ipnode file to define surface nodes

    character(len=*),intent(in) :: NODEFILE
    !     Local Variables
    integer :: i,ierror,np,np_global,&
         num_versions,nv
    integer,parameter :: num_derivs = 3
    character(len=132) :: ctemp1,readfile
    character(len=60) :: sub_name = 'define_node_geometry_2d'
    
    call enter_exit(sub_name,1)

    open(10, file=nodefile, status='old')

    !.....read in the total number of nodes. read each line until one is found
    !.....that has the correct keyword (nodes). then return the integer that is
    !.....at the end of the line
    read_number_of_nodes : do !define a do loop name
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !read a line into ctemp1
        if(index(ctemp1, "nodes")> 0) then !keyword "nodes" is found in ctemp1
          num_nodes_2d = get_final_integer(ctemp1) !return the final integer
          exit read_number_of_nodes !exit the named do loop
        endif
    end do read_number_of_nodes

    !write(*,*) 'Number of nodes are:',num_nodes_2d

!!!allocate memory to arrays that require node number
    if(.not.allocated(nodes_2d)) allocate(nodes_2d(num_nodes_2d))
    if(.not.allocated(node_xyz_2d)) allocate(node_xyz_2d(4,10,16,num_nodes_2d))
    if(.not.allocated(node_versn_2d)) allocate(node_versn_2d(num_nodes_2d))

   
    !.....read the coordinate, derivative, and version information for each node. 
    np=0
    read_a_node : do !define a do loop name
       !.......read node number
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1
       if(index(ctemp1, "Node")> 0) then
          np_global = get_final_integer(ctemp1) !get node number

          np=np+1
          nodes_2d(np) = np_global
          
          !.......read coordinates
          do i=1,3 ! for the x,y,z coordinates
             num_versions = 0
             read(unit=10, fmt="(a)", iostat=ierror) ctemp1
             if(index(ctemp1, "versions")> 0) num_versions = get_final_integer(ctemp1)
             node_versn_2d(np) = max(1,num_versions) !number of versions for node np
             do nv=1,node_versn_2d(np)
                if(num_versions > 1)then
                   read(unit=10, fmt="(a)", iostat=ierror) ctemp1 ! "For version number..."
                endif
                !...........coordinate          
                if(num_versions > 0) read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                node_xyz_2d(1,nv,i,np) = get_final_real(ctemp1)
                if(num_derivs.ge.1)then
                   !..........derivative 1
                   read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                   node_xyz_2d(2,nv,i,np) = get_final_real(ctemp1)
                endif
                if(num_derivs.ge.2)then
                   !..........derivative 2
                   read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                   node_xyz_2d(3,nv,i,np) = get_final_real(ctemp1)
                endif
                if(num_derivs.ge.3)then
                   !...........derivative 1&2
                   read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                   node_xyz_2d(4,nv,i,np) = get_final_real(ctemp1)
                endif
                if(num_derivs.ge.4)then
                   write(*,'(''This code is only valid for a surface geometry'')')
                   read(*,*)
                endif
             end do !nv
          end do !i
       endif !index
       if(np.ge.num_nodes_2d) exit read_a_node
    end do read_a_node
    
    close(10)
    !write(*,*) 'Here is:',node_xyz_2d(1,1,1,2),nodes_2d(2) 
    call enter_exit(sub_name,2)

  end subroutine define_node_geometry_2d


!!!##################################################

  subroutine define_data_geometry(datafile)
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_DEFINE_DATA_GEOMETRY" :: DEFINE_DATA_GEOMETRY


!!! read data points from a file
    
!!! dummy arguments
    character(len=*) :: datafile
!!! local variables
    integer :: iend,ierror,length_string,ncount,nj,itemp
    character(len=132) :: buffer,readfile
    character(len=60) :: sub_name = 'define_data_geometry'

    call enter_exit(sub_name,1)
    

    !set the counted number of data points to zero
    ncount = 0
    
    !readfile = trim(datafile)//'.ipdata'
    open(10, file=datafile, status='old')
    read(unit=10, fmt="(a)", iostat=ierror) buffer

!!! first run through to count the number of data points
    read_line_to_count : do
       read(unit=10, fmt="(a)", iostat=ierror) buffer
       if(ierror<0) exit !ierror<0 means end of file
       ncount = ncount + 1
    end do read_line_to_count
    num_data = ncount
    close (10)
    write(*,'('' Read'',I7,'' data points from file'')') num_data

!!! allocate arrays now that we know the size required
    if(allocated(data_xyz)) deallocate(data_xyz)
    if(allocated(data_weight)) deallocate(data_weight)
    allocate(data_xyz(3,num_data))
    allocate(data_weight(3,num_data))

!!! read the data point information
    !readfile = trim(datafile)//'.ipdata'
    open(10, file=datafile, status='old')
    read(unit=10, fmt="(a)", iostat=ierror) buffer

    !set the counted number of data points to zero
    ncount = 0
    read_line_of_data : do
       
       ! read the data #; z; y; z; wd1; wd2; wd3 for each data point
       read(unit=10, fmt="(a)", iostat=ierror) buffer
       if(ierror<0) exit !ierror<0 means end of file
       length_string = len_trim(buffer) !length of buffer, and removed trailing blanks
       
       ! read data number
       buffer=adjustl(buffer) !remove leading blanks
       iend=index(buffer," ",.false.)-1 !index returns location of first blank
       if(length_string == 0) exit
       ncount=ncount+1
       read (buffer(1:iend), '(i6)') itemp
       
       do nj=1,3
          ! read x,y,z coordinates
          buffer = adjustl(buffer(iend+1:length_string)) !remove data number from string
          buffer = adjustl(buffer) !remove leading blanks
          length_string = len(buffer) !new length of buffer
          iend=index(buffer," ",.false.)-1 !index returns location of first blank
          read (buffer(1:iend), '(D25.17)') data_xyz(nj,ncount)
       enddo !nj
       
       do nj=1,3
          !        ! read weightings
          !        buffer = adjustl(buffer(iend+1:length_string)) !remove data number from string
          !        buffer = adjustl(buffer) !remove leading blanks
          !        length_string = len(buffer) !new length of buffer
          !        iend=index(buffer," ",.false.)-1 !index returns location of first blank
          !        read (buffer(1:iend), '(D25.17)') data_weight(nj,ncount)
          data_weight(nj,ncount)=1.0_dp
       enddo !nj
       
    enddo read_line_of_data
    
    close(10)
    call enter_exit(sub_name,2)

  end subroutine define_data_geometry
!
!###################################################################################

  subroutine list_mesh_statistics(filename)
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_LIST_MESH_STATISTICS" :: LIST_MESH_STATISTICS

!     nindex(i) stores the order type (i=1,2,3,4 for generation, Horsfield
!     order, Strahler order, diameter-defined Strahler order
!     respectively) for the current element.

    use arrays,only: dp,data_xyz,data_weight,num_data,elem_ordrs,elem_units_below,&
         node_xyz,num_units
    use arrays,only: elem_cnct,elem_field,elem_symmetry,num_elems
    use indices,only: ne_length,ne_radius
    use diagnostics,only: enter_exit
    use other_consts,only: PI

    implicit none

     ! Parameters
     character(len=*),intent(in) :: filename

     ! Local Variables

     real(dp) :: translate(3)
     real(dp),allocatable :: data_temp(:,:),vertex_xyz(:,:)
     real(dp) :: boxrange(3),cofm_surfaces(3),offset=-2.0_dp,point_xyz(3)
     integer :: np0,np1,np2

     integer :: index(4),ntotaln(20)
     integer, allocatable :: nbranches(:,:)
     real(dp),allocatable ::  stats(:,:), branches(:,:), std_dev(:,:,:)
     real(dp),allocatable ::  sum_mean(:,:,:)
!    real(dp),allocatable ::  n_terminal(:,:,:),ntally(:,:,:)
     real(dp),allocatable ::  ntally(:,:,:)
     integer,allocatable ::  n_terminal(:), NORD(:,:)

     real(dp) :: xp1(3), xp2(3), xp3(3), sum_term
     integer :: i,j,k,nj,nline,nn,num_data_estimate,num_vertices
     integer :: ne, ne_child
     logical :: internal
     character(len=1) :: char1
     character(len=100) :: writefile

     integer :: ne_flow, N, num_ddp, num_llp
     integer :: ne0,ne_next,nmax_order,n_order,n_segments,num_in_order(40)
     integer :: average_term_gen,n_branch,ne1,ne2,ne_major,ne_minor,noelem,NPNE
     integer :: nbins(5), bins(5), nindex(4)
     real(dp) :: flow_proportion,mean_diameter,mean_sum_deadspace,mean_sum_pathlength,&
         ratio_L_D(40),sum_crosssec(40),sum_deadspace,&
         sum_diameter(40),sum_length(40),sum_pathlength,sum_volume(40),total_length,&
         max_bin,min_bin
     real(dp) :: means(20),sdt(20),ratios(4,3),r_squared(4,3)
     real(dp) :: angle, slope

     real(dp),allocatable :: sizes_from_entry(:,:),diameters(:),lengths(:)
     real(dp),allocatable :: YREGRESS(:,:)
     real(dp),parameter :: max_path=310.0_dp,min_path=130.0_dp
     real(dp),parameter :: max_dead=16.0_dp,min_dead=4.0_dp
     real(dp),parameter :: undefined = 10000_dp

     integer :: nbin,count_bins(2,100)
     integer,parameter :: num_bins = 20
     integer,parameter :: NUM_SCHEMES = 3
     logical :: add, DIAM_STRAHLER
     character(len=200) :: opfile

     character(len=60) :: sub_name = 'list_mesh_statistics'
     call enter_exit(sub_name,1)

     if(.not.allocated(sizes_from_entry)) allocate(sizes_from_entry(2,num_elems))
     if(.not.allocated(stats)) allocate(stats(21,num_elems))
     if(.not.allocated(branches)) allocate(branches(10,num_elems))

     if(.not.allocated(sum_mean)) allocate(sum_mean(4,6,num_elems))
     if(.not.allocated(std_dev)) allocate(std_dev(4,6,num_elems))

! QUESTION FOR MHT - are these still needed ?
! What is NBJ and NPNE ???
! MHT, Do we need still have nbranches ?
! QUESTION: WHat is NMAX_GEN and NUM_SCHEMES ?

     if(.not.allocated(n_terminal)) allocate(n_terminal(num_elems))
     if(.not.allocated(nbranches)) allocate(nbranches(5,num_elems))

     if(.not.allocated(diameters)) allocate(diameters(num_elems))
     if(.not.allocated(lengths)) allocate(lengths(num_elems))

     if(.not.allocated(YREGRESS)) allocate(YREGRESS(num_elems,3))
     if(.not.allocated(ntally)) allocate(ntally(4,6,num_elems))

     sizes_from_entry = 0.0_dp

     num_in_order(:) = 0
     nmax_order = 0
     sum_crosssec = 0.0_dp
     sum_diameter = 0.0_dp
     sum_length = 0.0_dp
     sum_volume = 0.0_dp
     ratio_L_D = 0.0_dp

     n_branch = 0
     num_ddp = 0
     num_llp = 0

     do ne = 1,num_elems
       ne0 = elem_cnct(-1,1,ne)
       add = .false.
       if(ne0.eq.0)then !stem of tree
          add = .true.
       else if(ne0.ne.0.and.elem_ordrs(2,ne0).ne.elem_ordrs(2,ne))then  !only count as extra branch if at start
          add = .true.
       endif

       if(add)then
          n_branch = n_branch + 1
          n_order = elem_ordrs(2,ne)
          nmax_order = max(n_order,nmax_order)
          num_in_order(n_order) = num_in_order(n_order) + 1
          ! Add length of all segments along branch, calculate mean diameter
          n_segments = 1
          mean_diameter = elem_field(ne_radius,ne) * 2.0_dp
          total_length = elem_field(ne_length,ne)
          ne_next=ne
          do while(elem_cnct(1,0,ne_next).eq.1.and.elem_symmetry(ne_next).eq.1)
            ne_next = elem_cnct(1,1,ne_next) !next segment
            total_length = total_length + elem_field(ne_length,ne_next) !sum lengths
            mean_diameter = mean_diameter  & 
                     + elem_field(ne_radius,ne_next) * 2.0
            n_segments = n_segments + 1 !count number of segments in branch

            diameters(ne) = elem_field(ne_radius,ne_next)*2.0
            lengths(ne) =   elem_field(ne_length,ne_next)
          enddo
          mean_diameter = mean_diameter/dble(n_segments)

          ! mean branch diameter
          sum_diameter(n_order) = sum_diameter(n_order) + mean_diameter

          ! total branch length
          sum_length(n_order) = sum_length(n_order) + total_length

          ! total branch volume
          sum_volume(n_order) = sum_volume(n_order) +  & 
                      total_length*PI*(mean_diameter**2)/4.0

          ! total branch cross-sectional area
          sum_crosssec(n_order) = sum_crosssec(n_order) +  & 
                      PI*(mean_diameter**2)/4.0_dp

          ! branch L:D ratio
          ratio_L_D(n_order) = ratio_L_D(n_order) + total_length/mean_diameter

!! HARI MIGRATED
!         IF(INDEX(1).GT.1)THEN
!HBK       nb=NBJ(1,ne)
           np0=elem_cnct(1,0,ne) !start of parent
           np1=elem_cnct(1,1,ne) !start node
           np2=elem_cnct(1,2,ne) !end node
           DO nj=1,3
            xp1(nj)=node_xyz(nj,np0)
            xp2(nj)=node_xyz(nj,np1)
            xp3(nj)=node_xyz(nj,np2)
           ENDDO !nj
           call mesh_angle(angle,xp1,xp2,xp3)
           stats(2,ne)=angle !temporary storage of angle
           branches(3,N)=angle*180.d0/PI !store the branching angle to parent

           ntotal=ntotal+1
           IF(diameters(ne0).GT.0.d0.AND.diameters(ne).GT.0.d0)THEN
            IF(diameters(ne)/diameters(ne0).LE.1.d0)THEN
              num_ddp = num_ddp + 1
            ENDIF
           ENDIF

           IF(diameters(ne0).GE.4.d0) THEN
             nbins(1)=nbins(1)+1
             bins(1)=bins(1)+angle
           ELSE IF(diameters(ne0).GE.2.d0) THEN
             nbins(2)=nbins(2)+1
             bins(2)=bins(2)+angle
           ELSE IF(diameters(ne0).GE.1.d0) THEN
             nbins(3)=nbins(3)+1
             bins(3)=bins(3)+angle
           ELSE IF(diameters(ne0).GE.0.7d0) THEN
             nbins(4)=nbins(4)+1
             bins(4)=bins(4)+angle
           ENDIF
!          ELSE
!            branches(3,N)=undefined
!          ENDIF

!        ne_next=ne
!        do while(elem_cnct(1,0,ne_next).eq.1.and.elem_symmetry(ne_next).eq.1)

           stats(5,ne)=lengths(ne)
           ne_next=ne
           if(elem_cnct(1,0,ne_next).EQ.1) then
             ne_next=elem_cnct(1,1,ne_next)
             stats(5,ne)=stats(5,ne)+lengths(ne_next)
           endif

          if(elem_cnct(1,0,ne).EQ.0) then
            n_terminal(1)=n_terminal(1)+1
          endif

          !CC... Geometric properties of mesh
          branches(4,N)=undefined !initialise to no rotation angle
          IF(elem_cnct(-1,0,ne).GT.0.AND.elem_cnct(1,0,ne).GT.1) then
            ne0=elem_cnct(-1,1,ne)
            IF(elem_cnct(1,0,ne0).GT.1) then
              CALL mesh_plane_angle(NBJ,ne,NPNE,angle,XP)
              branches(4,N)=angle*180.d0/PI !rotation angle
            ENDIF
          ENDIF

        endif
     enddo !ne

!   DO ne=1,num_elems
!      ne0=elem_cnct(-1,1,ne) !parent
!      np1=elem_nodes(1,ne)
!      np2=elem_nodes(2,ne)
!      if(elem_cnct(1,0,ne).eq.0)then !terminal
!       elems(ne) = ne_global


     !=======================
     !=======================
     do noelem = 1,num_elems
!       ne=NELIST(noelem)
        ne0=elem_cnct(-1,1,noelem) !parent

        DO j=11,21 !initialise values for the summary statistics
          stats(j,ne)=undefined
        ENDDO !j
        ne0=elem_cnct(-1,1,ne)
        IF(ne0.NE.0)THEN
          IF(elem_ordrs(1,ne0).NE.elem_ordrs(1,ne))THEN
            IF(stats(5,ne)/stats(5,ne0).LE.1.d0)THEN
              num_llp=num_llp+1
            ENDIF
            stats(19,ne)=stats(5,ne)/stats(5,ne0) !L/Lparent
            IF(diameters(ne0).GT.0.d0.AND.diameters(ne).GT.0.d0)THEN
              stats(16,ne)=diameters(ne)/diameters(ne0) !D/Dparent
            ENDIF
          ENDIF
        ENDIF

        IF(elem_cnct(1,0,ne).GE.2) THEN !'bifurcations only
            ne1=elem_cnct(1,1,ne) !first child
            ne2=elem_cnct(1,2,ne) !second child

            !XX.. Summary statistics
            IF(stats(6,ne1).LT.undefined &
                       .AND.stats(6, ne2).LT.undefined)THEN
              IF(stats(6,ne1).GE.stats(6,ne2))THEN !diameter classification
                ne_major=ne1
                ne_minor=ne2
              ELSE
                ne_major=ne2
                ne_minor=ne1
              ENDIF
              IF(stats(2,ne_minor).LT.undefined.AND.stats(2,ne_major).LT.undefined)THEN
                stats(11,ne)=stats(2,ne_minor)*180.d0/PI
                stats(12,ne)=stats(2,ne_major)*180.d0/PI
              ENDIF

              IF(diameters(ne_minor).GT.0.d0.AND. &
                diameters(ne_major).GT.0.d0)THEN
                stats(13,ne)=stats(5,ne_minor)/diameters(ne_minor) !L/D minor
                stats(14,ne)=stats(5,ne_major)/diameters(ne_major) !L/D major
                stats(15,ne)=diameters(ne_minor)/diameters(ne_major) !minor D / major D
                stats(17,ne)=diameters(ne_minor)/diameters(ne) !minor D / D parent
                stats(18,ne)=diameters(ne_major)/diameters(ne) !major D / D parent
              ENDIF
              IF(stats(5,ne1).LE.stats(5,ne2))THEN !length classification
                ne_major=ne1
                ne_minor=ne2
              ELSE
                ne_major=ne2
                ne_minor=ne1
              ENDIF !length criteria
              stats(20,ne)=stats(5,ne_major)/stats(5,ne_minor)
            ENDIF

!!          IF(DOP)THEN
!!            WRITE(OP_STRING,'('' ne'',I5,''  angle'',F6.2,'
!!   &          //'''  L/Dmin'',F5.2,''  L/Dmaj'',F5.2)')
!!   &          ne,stats(2,ne)*180.d0/PI,stats(13,ne),stats(14,ne)
!!            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
!!          ENDIF !DOP

        ENDIF !elem_cnct
      ENDDO !noelem
      !=======================
      !=======================


      !C -------------------------------------------------------------
      !C .. Calculate mean branching statistics from values in branches
      DO N=1,n_branch

        !i : for generations, Horsfield orders, Strahler orders
        DO i=1,NUM_SCHEMES
        nindex(i)=nbranches(i,N)

         !.............................................
         !C... length and diameter
         DO j=1,2
           sum_mean(i,j,nindex(i))=sum_mean(i,j,nindex(i))+ &
             branches(j,N)
           IF(i.EQ.3.AND.j.EQ.1)THEN
             IF(nindex(i).NE.nbranches(5,N))THEN !not same as parent 
               ntally(i,j,nindex(i))=ntally(i,j,nindex(i))+1
             ENDIF
           ELSE
             ntally(i,j,nindex(i))=ntally(i,j,nindex(i))+1
           ENDIF
         ENDDO !j
         !C.. length and diameter
         !.............................................


         !...........................................
         !CCC. branching angle and rotation angle(j=3,4)
         DO j=3,4
           IF(branches(j,N).LT.undefined)THEN
             sum_mean(i,j,nindex(i))=sum_mean(i,j,nindex(i))+branches(j,N)
             ntally(i,j,nindex(i))=ntally(i,j,nindex(i))+1
           ENDIF
         ENDDO !j
         !CCC . branching angle and rotation angle 
         !...........................................


         !..........................................
         !. ratio of L:D (j=5)
         j=5
         IF(branches(j,N).LT.undefined)THEN
           sum_mean(i,j,nindex(i))=sum_mean(i,j,nindex(i))+branches(j,N)
           ntally(i,j,nindex(i))=ntally(i,j,nindex(i))+1
         ENDIF
         !..........................................

         ENDDO !i

         !C.. Summary statistics from branches
         DO j=3,5 !branching angle, rotation angle, L/D
           IF(branches(j,N).LT.undefined)THEN
             means(j-2)=means(j-2)+branches(j,N)
           ENDIF
         ENDDO !j

      ENDDO !N
     !C------------------------------------------------------------
     !C.Calculate mean branching statistics from values in branches

      DO N=1,GENM
        DO j=3,5
          ntotaln(j-2)=ntotaln(j-2)+ntally(1,j,N)
       ENDDO !j
      ENDDO !N

      DO N=1,GENM
        DO i=1,NUM_SCHEMES
          DO j=1,5
            IF(ntally(i,j,N).GT.0)THEN
              sum_mean(i,j,N)=sum_mean(i,j,N)/DBLE(ntally(i,j,N))
              NMAX_GEN(i)=N
             ELSE
              sum_mean(i,j,N)=0.d0
            ENDIF
          ENDDO !j
        ENDDO !i
      ENDDO !N

      DO N=1,5
        IF(nbins(N).NE.0)THEN
          bins(N)=bins(N)/DBLE(nbins(N))*180.d0/PI
        ENDIF
      ENDDO !N

       !CC. Summary statistics from branches
       DO j=3,5 !branching angle, rotation angle, L/D
         IF(ntotaln(j-2).NE.0)THEN
           means(j-2)=means(j-2)/DBLE(ntotaln(j-2))
         ELSE
           means(j-2)=0.d0
         ENDIF
       ENDDO !j

       i=2 !Horsfield orders
       j=6 !Nw/Nw-1
       DO N=1,GENM-1
         IF(ntally(i,1,N).GT.0.AND.ntally(i,1,N+1).GT.0)THEN
           sum_mean(i,j,N)=DBLE(ntally(i,1,N))/DBLE(ntally(i,1,N+1))
         ELSE
           sum_mean(i,j,N)=0.d0
         ENDIF
       ENDDO !N

       !CC. Summary statistics from CE
       DO noelem=1,NELIST(0)
         ne=NELIST(noelem)
         DO j=11,21
           IF(stats(j,ne).LT.undefined)THEN
             means(j-7)=means(j-7)+stats(j,ne)
             ntotaln(j-7)=ntotaln(j-7)+1
           ENDIF
         ENDDO !j
       ENDDO !noelem

       DO j=11,21
         IF(ntotaln(j-7).GT.0)THEN
           means(j-7)=means(j-7)/DBLE(ntotaln(j-7))
         ENDIF
       ENDDO !j
       !C... End of mean calculation
       !---------------------------------------

       !!--------------------------------------------------
       !C... Calculate the standard deviations
       !C... Sum of (value-mean)^2
       DO N=1,n_branch
          DO i=1,NUM_SCHEMES !for generations, Horsfield orders, Strahler orders
            nindex(i)=nbranches(i,N)
            DO j=1,5 !length, diameter, branching angle, rotation angle, L/D
              IF(branches(j,N).LT.undefined)THEN
                std_dev(i,j,nindex(i))=std_dev(i,j,nindex(i)) &
                  +(branches(j,N)-sum_mean(i,j,nindex(i)))**2
              ENDIF
            ENDDO !j
          ENDDO !i
          DO j=3,5 !branching angle, rotation angle, L/D
            IF(branches(j,N).LT.undefined)THEN
              sdt(j-2)=sdt(j-2) + &
                 (branches(j,N)-means(j-2))**2
            ENDIF
          ENDDO !j
       ENDDO !N

       do noelem = 1,num_elems
          ne=NELIST(noelem)
          DO j=11,21
            IF(stats(j,ne).LT.undefined)THEN
              sdt(j-7)=sdt(j-7)+(stats(j,ne)-means(j-7))**2
            ENDIF
          ENDDO !j
       ENDDO !noelem
       !!---------------------------------------------------


      !!-----------------------------------------------------
      !CC. SD = sqrt(1/(n-1)*sum)
      DO N=1,GENM
         DO i=1,NUM_SCHEMES !for generations, Horsfield orders, Strahler orders
           DO j=1,5 !length, diameter, branching angle, rotation angle, L/D
             IF(ntally(i,j,N).GT.1)THEN
               std_dev(i,j,N)=DSQRT(std_dev(i,j,N)/DBLE(ntally(i,j,N)-1))
             ELSE
               std_dev(i,j,N)=0.d0
             ENDIF
           ENDDO !j
         ENDDO !i
      ENDDO !N

      DO j=1,13
        IF(ntotaln(j).GT.1)THEN
          sdt(j)=DSQRT(sdt(j)/DBLE(ntotaln(j)-1))
        ENDIF
      ENDDO !j
      !CCC...... End of standard deviation calculation
      !----------------------------------------------------



     opfile = trim(filename)//'.opmesh'
     open(10, file=opfile, status='replace')

     write(*,'(/'' Horsfield   #branches     Length     Diameter    ratio L:D    Volume   CrossSecArea'')')
     write(*,'(''   order                    (mm)        (mm)                     (mL)      (cm^2)'')')
     write(10,'(/'' Horsfield   #branches     Length     Diameter    ratio L:D    Volume   CrossSecArea'')')
     write(10,'(''   order                    (mm)        (mm)                     (mL)      (cm^2)'')')

     do n_order = 1,nmax_order
       write(*,'(i6,i12,f14.2,f12.2,f12.2,f12.2,f12.2)') &
            n_order,&
            num_in_order(n_order),&
            sum_length(n_order)/dble(num_in_order(n_order)),&
            sum_diameter(n_order)/dble(num_in_order(n_order)),&
            ratio_L_D(n_order)/dble(num_in_order(n_order)), &
            sum_volume(n_order)/1.0e+3_dp,&
            sum_crosssec(n_order)/100.0_dp

       write(10,'(i6,i12,f14.2,f12.2,f12.2,f12.2,f12.2)') &
            n_order,&
            num_in_order(n_order),&
            sum_length(n_order)/dble(num_in_order(n_order)),&
            sum_diameter(n_order)/dble(num_in_order(n_order)),&
            ratio_L_D(n_order)/dble(num_in_order(n_order)), &
            sum_volume(n_order)/1.0e+3_dp,&
            sum_crosssec(n_order)/100.0_dp
     enddo


       sum_term = 0
       average_term_gen = 0
       !=======
       i=1
       !=======
       DO N=1,NMAX_GEN(i)
          WRITE(10,'(3(I10),5(F8.2,'' ('',F6.2,'')''))') &
            N,ntally(i,1,N),n_terminal(N),&
            sum_mean(i,1,N),std_dev(i,1,N), &
            sum_mean(i,2,N),std_dev(i,2,N),&
            sum_mean(i,3,N),std_dev(i,3,N),&
            sum_mean(i,4,N),std_dev(i,4,N),&
            sum_mean(i,5,N),std_dev(i,5,N)

!         CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          average_term_gen=average_term_gen + n_terminal(N)*N
          sum_term=sum_term + n_terminal(N)
        ENDDO
        IF(SUM_TERM.GT.0)THEN
          average_term_gen=average_term_gen/DBLE(sum_term)
        ELSE
          average_term_gen=0.d0
        ENDIF

        WRITE(10,'(/''Horsfield   #branches  Length  Diameter  Branching  Rotation  ratio L:D   Nw/Nw-1'')')
!       CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        WRITE(10,'('' order   (mm)   (mm)  angle(deg)   angle(deg)'')')
!       CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        WRITE(10,'(115(''-''))')
!       CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        !=======
        i=2
        !=======
        DO N=1,NMAX_GEN(i)
          WRITE(10,'(2(I10),5(F8.2,'' ('',F6.2,'')''),F8.2)') &
            N,ntally(2,1,N), &
            sum_mean(i,1,N),std_dev(i,1,N), &
            sum_mean(i,2,N),std_dev(i,2,N), &
            sum_mean(i,3,N),std_dev(i,3,N), &
            sum_mean(i,4,N),std_dev(i,4,N), &
            sum_mean(i,5,N),std_dev(i,5,N), &
            sum_mean(i,6,N)
        ENDDO

        WRITE(10,*) 'Strahler  #branches   Length   Diameter  Branching  Rotation  ratio L:D'
!       CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(10,*) ' order  (mm)   (mm)  angle(deg)  angle(deg)'
!       CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(10,*) '------------------------------------------'
!       CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        !=======
        i=3
        !=======
        DO N=1,NMAX_GEN(i)
          WRITE(*,'(2(I10),5(F8.2,'' ('',F6.2,'')''))') &
            N,ntally(3,1,N),&
            sum_mean(i,1,N),std_dev(i,1,N),&
            sum_mean(i,2,N),std_dev(i,2,N),&
            sum_mean(i,3,N),std_dev(i,3,N),&
            sum_mean(i,4,N),std_dev(i,4,N),&
            sum_mean(i,5,N),std_dev(i,5,N)
        ENDDO  
        

        !---------------------------------------------------
        IF(DIAM_STRAHLER) THEN
          WRITE(10,*)  &
            'Diam-Def  Strahler  #branches   Length  Diameter  Branching   Rotation  ratio  L:D'
!         CALL WRITES(IOFI,OP_STRING)
          WRITE(10,*) 'order    (mm)    (mm)   angle(deg)   angle(deg)'
!         CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

          WRITE(10,'(115(''-''))')
!         CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          
          i=4
          DO N=1,NMAX_GEN(i)
            WRITE(10,'(2(I10),5(F8.2,'' ('',F6.2,'')''))') &
              N,ntally(i,1,N), &
              sum_mean(i,1,N),std_dev(i,1,N), &
              sum_mean(i,2,N),std_dev(i,2,N), &
              sum_mean(i,3,N),std_dev(i,3,N), &
              sum_mean(i,4,N),std_dev(i,4,N), &
              sum_mean(i,5,N),std_dev(i,5,N)
!           CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO            
        ENDIF !DIAM_STRAHLER
        !---------------------------------------------------
        
        !-------------------------------------------------
        DO i=2,NUM_SCHEMES !Horsfield and Strahler orders
         DO N=1,NMAX_GEN(i)
          X(N)=N
          YREGRESS(N,1)=DLOG10(DBLE(ntally(i,1,N)))
          YREGRESS(N,2)=DLOG10(sum_mean(i,1,N))
          YREGRESS(N,3)=DLOG10(sum_mean(i,2,N))
         ENDDO !N
         DO j=1,3 !number of branches, length, diameter
!          CALL LINREGRESS(NMAX_GEN(i),r_squared(i,j),slope,X,YREGRESS(1,j))
           CALL linear_regression(NMAX_GEN(i),r_squared(i,j),slope,X,YREGRESS(1,j))
           ratios(i,j)=10.d0**DABS(slope)
         ENDDO !j
        ENDDO !i
        !----------------------------------------------------

       write(10,*) 'terminal generation  =', average_term_gen
       write(10,*) 'branching angle  =', means(1),sdt(1)
       write(10,*) 'rotation angle  =', means(2),sdt(2)
       write(10,*) 'minor angle  =', means(4),sdt(4)
       write(10,*) 'major angle  =', means(5),sdt(5)
       write(10,*) 'L/D ratio  =', means(3),sdt(3)
       write(10,*) 'L/D minor ratio  =', means(6),sdt(6)
       write(10,*) 'L/D major ratio  =', means(7),sdt(7)
       write(10,*) 'minorD/majorD=', means(8),sdt(8)
       write(10,*) 'D/Dparent=', means(9),sdt(9)
       write(10,*) '%D/Dparent<1 =', DBLE(num_ddp)/DBLE(ntotal)*100.d0
       write(10,*) 'Dmin/Dparent =', means(10),sdt(10)
       write(10,*) 'Dmax/Dparent =', means(11),sdt(11)
       write(10,*) 'L/Lp =', means(12),sdt(12)
       write(10,*) '%L/Lp<1 =', DBLE(num_llp)/DBLE(NELIST(0)-1)*100.d0
       write(10,*) 'L1/L2 (L1 < L2) =', means(13),sdt(13)

!          ' branching angle      = '',F7.3,'' ('',F6.3,'')'','
!          ' rotation angle       = '',F7.3,'' ('',F6.3,'')'','
!          ' minor angle          = '',F7.3,'' ('',F6.3,'')'','
!          ' major angle          = '',F7.3,'' ('',F6.3,'')'','
!          ' L/D                  = '',F7.3,'' ('',F6.3,'')'','
!          ' L/D minor child      = '',F7.3,'' ('',F6.3,'')'','
!          ' L/D major child      = '',F7.3,'' ('',F6.3,'')'','
!          ' minor D/major D      = '',F7.3,'' ('',F6.3,'')'','
!          ' D/Dparent            = '',F7.3,'' ('',F6.3,'')'','
!          ' %D/Dparent  < 1      = '',F7.3,'
!          ' Dmin/Dparent         = '',F7.3,'' ('',F6.3,'')'','
!          ' Dmaj/Dparent         = '',F7.3,'' ('',F6.3,'')'','
!          ' L/Lp                 = '',F7.3,'' ('',F6.3,'')'','
!          ' %L/Lp < 1            = '',F7.3,'
!          ' L1/L2 (L1 < L2)      = '',F7.3,'' ('',F6.3,'')'')')
!          average_term_gen,means(1),sdt(1),&
!          means(2),sdt(2),&
!          means(4),sdt(4),&
!          means(5),sdt(5),&
!          means(3),sdt(3),means(6),sdt(6),&
!          means(7),sdt(7),means(8),sdt(8),means(9),sdt(9),&
!          DBLE(num_ddp)/DBLE(ntotal)*100.d0,means(10),sdt(10),&
!          means(11),sdt(11),&
!          means(12),sdt(12),&
!          DBLE(num_llp)/DBLE(NELIST(0)-1)*100.d0,&
!          means(13),sdt(13)

!       CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        WRITE(10,*) ' Rb Horsfield    = ',ratios(2,1),r_squared(2,1)
        WRITE(10,*) ' Rl Horsfield    = ',ratios(2,2),r_squared(2,2)
        WRITE(10,*) ' Rd Horsfield    = ',ratios(2,3),r_squared(2,3)
        WRITE(10,*) ' Rb Strahler    = ',ratios(3,1),r_squared(3,1)
        WRITE(10,*) ' Rl Strahler    = ',ratios(3,2),r_squared(3,2)
        WRITE(10,*) ' Rd Strahler    = ',ratios(3,3),r_squared(3,3)

!       WRITE(OP_STRING,'('' Rb Strahler    = '',F7.3,'
!    '    //''' Rsq ='',F6.3,'//
!    '    '/'' Rl Strahler          = '',F7.3,'' Rsq = '',F6.3,'//
!    '    '/'' Rd Strahler          = '',F7.3,'' Rsq = '',F6.3,'//
!    '    '/'' Rb Horsfield         = '',F7.3,'' Rsq = '',F6.3,'//
!    '    '/'' Rl Horsfield         = '',F7.3,'' Rsq = '',F6.3,'//
!    '    '/'' Rd Horsfield         = '',F7.3,'' Rsq = '',F6.3)')
!    '    ratios(3,1),r_squared(3,1),ratios(3,2),r_squared(3,2),
!         ratios(3,3),r_squared(3,3),&
!         ratios(2,1),r_squared(2,1),&
!         ratios(2,2),r_squared(2,2),&
!         ratios(2,3),r_squared(2,3)

! CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        IF(DIAM_STRAHLER) THEN
         write(10,*) 'Rb diam-def Strahler=',ratios(4,1),'Rsq=',r_squared(4,1)
         write(10,*) 'Rl diam-def Strahler=',ratios(4,2),'Rsq=',r_squared(4,2)
         write(10,*) 'Rd diam-def Strahler=',ratios(4,3),'Rsq=',r_squared(4,3)
!         WRITE(OP_STRING,
!    &      '('' Rb diam-def Strahler  ='',F7.3,'' Rsq = '',F6.3,'//
!    &      '/'' Rl diam-def Strahler  ='',F7.3,'' Rsq = '',F6.3,'//
!    &      '/'' Rd diam-def Strahler  ='',F7.3,'' Rsq = '',F6.3)')
!    &      ratios(4,1),r_squared(4,1),&
!           ratios(4,2),r_squared(4,2),&
!           ratios(4,3),r_squared(4,3)
!         CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF

        WRITE(10,*) 'mean angle Dp 4.0+   = ',bins(1)
        WRITE(10,*) 'mean angle Dp 3.0+   = ',bins(2)
        WRITE(10,*) 'mean angle Dp 2.0+   = ',bins(3)
        WRITE(10,*) 'mean angle Dp 1.0+   = ',bins(4)
        WRITE(10,*) 'mean angle Dp 0.7+   = ',bins(5)

!       WRITE(OP_STRING,'('' mean angle Dp 4.0+   = '',F7.3,'//
!    &    '/''  mean angle Dp 3.0+   = '',F7.3,'//
!    &    '/''  mean angle Dp 2.0+   = '',F7.3,'//
!    &    '/''  mean angle Dp 1.0+   = '',F7.3,'//
!    &    '/''  mean angle Dp 0.7+   = '',F7.3)')
!    '    (bins(j),j=1,5)


     write(*,'('' Total anatomical dead space = '',f8.2, '' mL'')') &
            sum(sum_volume)/1.0e+3_dp*3.0_dp/2.0_dp

     mean_sum_pathlength = 0.0_dp
     mean_sum_deadspace = 0.0_dp

     count_bins = 0

     do ne = 1,num_elems
       if(elem_cnct(1,0,ne).eq.0)then !terminal
          sum_pathlength = elem_field(ne_length,ne)
!!        sum_deadspace = elem_field(ne_vol,ne)
!!        flow_proportion = 1.0_dp
          ne0 = elem_cnct(-1,1,ne)
          ne_child = ne
          do while (ne0.ne.0) !while not the stem
             sum_pathlength = sum_pathlength + elem_field(ne_length,ne0)
!!           flow_proportion = flow_proportion * elem_field(ne_flow,ne_child)/ &
!!                elem_field(ne_flow,ne0)
!!           sum_deadspace = sum_deadspace + elem_field(ne_vol,ne0) * flow_proportion
             ne_child = ne0
             ne0 = elem_cnct(-1,1,ne0)
          enddo !while not the stem

          mean_sum_pathlength = mean_sum_pathlength + sum_pathlength
          mean_sum_deadspace = mean_sum_deadspace + sum_deadspace
          do nbin = 1,num_bins
             min_bin = min_path+(nbin-1)*(max_path-min_path)/num_bins
             max_bin = min_path+nbin*(max_path-min_path)/num_bins
             if(sum_pathlength.gt.min_bin.and.sum_pathlength.le.max_bin) count_bins(1,nbin) = count_bins(1,nbin)+1
             min_bin = min_dead+(nbin-1)*(max_dead-min_dead)/num_bins
             max_bin = min_dead+nbin*(max_dead-min_dead)/num_bins
             if(sum_deadspace.gt.min_bin.and.sum_deadspace.le.max_bin) count_bins(2,nbin) = count_bins(2,nbin)+1
          enddo
       endif !terminal
     enddo !ne
     mean_sum_pathlength = mean_sum_pathlength/elem_units_below(1)
     mean_sum_deadspace = mean_sum_deadspace/elem_units_below(1)

     write(*,'('' Mean pathlength to acinus = '',f8.2, '' mm; ['',f8.2,''-'',f8.2,'']'')') &
         mean_sum_pathlength,min_path,max_path
     write(*,'('' Mean deadspace to acinus  = '',f8.2, '' mL; ['',f8.2,''-'',f8.2,'']'')') &
         mean_sum_deadspace/1.0+3_dp,min_dead/1.0+3_dp,max_dead/1.0+3_dp


     write(*,'(''........................... '')')
     write(*,'('' N_bin     N_PL_bin    N_Vd_bin'')')
     do nbin = 1,num_bins
       write(*,'(i4,i12,i12)') nbin,count_bins(1,nbin),count_bins(2,nbin)
     enddo

     close(10)
     deallocate(sizes_from_entry)
     deallocate(stats)
     deallocate(branches)

  end subroutine list_mesh_statistics


! 
  subroutine triangles_from_surface(num_triangles,num_vertices,surface_elems,triangle,vertex_xyz)
      
     !#This subroutine is called by make_data_grid subroutine.
     !#Application: to grow a grid in 2D surface.
     
     integer,intent(in) :: surface_elems(:)
     integer,allocatable :: triangle(:,:)
     real(dp),allocatable :: vertex_xyz(:,:)
     integer,parameter :: ndiv = 3
     integer :: i,j,ne,num_surfaces,num_triangles,num_tri_vert,num_vertices
     integer :: index1,index2,nmax_1,nmax_2,step_1,step_2,nvertex_row
     real(dp) :: X(3),xi(3)
     character(len=3) :: repeat
     logical :: four_nodes
     character(len=60) :: sub_name = 'triangles_from_surface'
     call enter_exit(sub_name,1)
     
     if(.not.allocated(triangle)) allocate(triangle(3,2*num_elems_2d*ndiv**2))
     if(.not.allocated(vertex_xyz)) allocate(vertex_xyz(3,num_elems_2d*(ndiv+1)**2))
     triangle = 0
     vertex_xyz = 0.0_dp
     num_surfaces = count(surface_elems.ne.0)
     num_triangles = 0
     num_vertices = 0
     num_tri_vert = 0 


     do ne=1,num_elems_2d
        four_nodes = .false.
        repeat = '0_0'
        if(elem_nodes_2d(1,ne).eq.elem_nodes_2d(2,ne)) repeat = '1_0'
        if(elem_nodes_2d(1,ne).eq.elem_nodes_2d(3,ne)) repeat = '2_0'
        if(elem_nodes_2d(2,ne).eq.elem_nodes_2d(4,ne)) repeat = '2_1'
        if(elem_nodes_2d(3,ne).eq.elem_nodes_2d(4,ne)) repeat = '1_1'

        select case(repeat)
        case ('0_0')
        
           nmax_1 = ndiv+1 ! ndiv+1 vertices in xi1 direction
           step_1 = 0      ! # of vertices in xi1 is constant
           nmax_2 = ndiv+1 ! ndiv+1 vertices in xi2 direction
           step_2 = 0      ! # of vertices in xi2 is constant
           index1 = 1
           index2 = 2
           four_nodes = .true.

        case ('1_0')
               
           nmax_1 = 1      ! start with 1 vertex in xi1 direction
           step_1 = 1      ! increase # of vertices in xi1 with each step in xi2
           nmax_2 = ndiv+1 ! ndiv+1 vertices in xi2 direction
           step_2 = 0      ! # of vertices in xi2 is constant
           index1 = 1
           index2 = 2

        case ('1_1')

           nmax_1 = ndiv+1 ! start with ndiv+1 vertices in xi1 direction
           step_1 = -1     ! decrease # of vertices in xi1 with each step in xi2
           nmax_2 = ndiv+1 ! ndiv+1 vertices in xi2 direction
           step_2 = 0      ! # of vertices in xi2 is constant
           index1 = 1
           index2 = 2

        case ('2_0')

           nmax_2 = ndiv+1 ! ndiv+1 vertices in xi1 direction
           step_2 = 0      ! # of vertices in xi1 is constant
           nmax_1 = 1      ! start with 1 vertex in xi2 direction
           step_1 = 1      ! increase # of vertices in xi2 with each step in xi1          
           index2 = 1
           index1 = 2

        case ('2_1')

           nmax_2 = ndiv+1 ! ndiv+1 vertices in xi1 direction
           step_2 = 0      ! # of vertices in xi1 is constant
           nmax_1 = ndiv+1 ! start with ndiv+1 vertices in xi2 direction
           step_1 = -1     ! decrease # of vertices in xi2 with each step in xi1
           index2 = 1
           index1 = 2
        end select

        xi(index2) = 0.0_dp
        do i = 1,nmax_2
           xi(index1) = 0.0_dp
           do j = 1,nmax_1
              num_vertices = num_vertices + 1
              X = coord_at_xi(ne,xi,'hermite')
              vertex_xyz(1:3,num_vertices) = X(1:3)
              if(nmax_1.gt.1) xi(index1) = xi(index1) + 1.0_dp/(nmax_1-1)
              if(i.gt.1.and.j.gt.1)then
                 num_triangles = num_triangles + 1
                 triangle(1,num_triangles) = num_vertices
                 triangle(2,num_triangles) = num_vertices-1
                 triangle(3,num_triangles) = nvertex_row+j-1
                 if(four_nodes.or.(.not.four_nodes.and.j.lt.nmax_1).or.step_1.eq.-1)then
                    num_triangles = num_triangles + 1
                    triangle(1,num_triangles) = num_vertices
                    triangle(2,num_triangles) = nvertex_row+j
                    triangle(3,num_triangles) = nvertex_row+j-1
                 endif
                 if(step_1.eq.-1.and.j.eq.nmax_1)then
                    num_triangles = num_triangles + 1
                    triangle(1,num_triangles) = num_vertices
                    triangle(2,num_triangles) = nvertex_row+j+1
                    triangle(3,num_triangles) = nvertex_row+j
                 endif
              else if(step_1.eq.-1.and.i.eq.nmax_2.and.j.eq.1)then
                      num_triangles = num_triangles + 1
                      triangle(1,num_triangles) = num_vertices
                      triangle(2,num_triangles) = num_vertices-1
                      triangle(3,num_triangles) = num_vertices-2
              endif
           enddo !j
           nvertex_row = num_vertices-nmax_1 !record first vertex # in previous row
           if(nmax_2.gt.1) xi(index2) = xi(index2) + 1.0_dp/(nmax_2-1)
           nmax_1 = nmax_1 + step_1
       enddo !i
     end do

     write(*,'('' Made'',I8,'' triangles to cover'',I6,'' surface elements'')') &
           num_triangles,num_elems_2d

     call enter_exit(sub_name,2)
  
  end subroutine triangles_from_surface

! 
!###################################################################################
! 
  subroutine make_data_grid(surface_elems,spacing,to_export,filename,groupname)
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_MAKE_DATA_GRID" :: MAKE_DATA_GRID

     ! Parameters
     integer,intent(in) :: surface_elems(:)
     real(dp),intent(in) :: spacing
     logical,intent(in) :: to_export
     character(len=*),intent(in) :: filename
     character(len=*),intent(in) :: groupname

     ! Local Variables
     real(dp) :: max_bound(3),min_bound(3),scale_mesh,translate(3)
     integer,allocatable :: triangle(:,:)
     real(dp),allocatable :: data_temp(:,:),vertex_xyz(:,:)
     real(dp) :: boxrange(3),cofm_surfaces(3),offset=-2.0_dp,point_xyz(3)
     integer :: i,j,k,ne,nj,nline,nn,num_data_estimate,num_triangles,num_vertices
     logical :: internal
     character(len=1) :: char1
     character(len=100) :: writefile
     character(len=60) :: sub_name = 'make_data_grid'

     call enter_exit(sub_name,1)

     call triangles_from_surface(num_triangles,num_vertices,surface_elems,triangle,vertex_xyz)

     if(offset.gt.0.0_dp)then
   !!! generate within a scaled mesh, then return to original size afterwards
   !!! this option is used when we don't want the seed points too close to the boundary
       scale_mesh = 1.0_dp-(offset/100.0_dp)
       cofm_surfaces = sum(vertex_xyz,dim=2)/size(vertex_xyz,dim=2)
       translate = cofm_surfaces * (scale_mesh - 1.0_dp)
       forall (i=1:num_vertices) vertex_xyz(1:3,i) = vertex_xyz(1:3,i)*scale_mesh + translate(1:3)
     endif

  !!! find the bounding coordinates for the surface mesh

     min_bound = minval(vertex_xyz,2)
     max_bound = maxval(vertex_xyz,2)
     boxrange = max_bound - min_bound
     num_data_estimate = int(dble((boxrange(1)/spacing+1.0_dp) * &
                         (boxrange(2)/spacing+1.0_dp) * (boxrange(3))/spacing+1.0_dp) * &
                         volume_internal_to_surface(triangle,vertex_xyz)/(boxrange(1)*boxrange(2)*boxrange(3)))


  !!! allocate arrays based on estimated number of data points

     if(allocated(data_xyz)) deallocate(data_xyz)
       allocate(data_xyz(3,num_data_estimate))
       i=0
       num_data = 0
       point_xyz = min_bound + 0.5_dp*spacing 
       do while(point_xyz(3).le.max_bound(3)) ! for z direction
          i=i+1
          j=0 
          do while(point_xyz(2).le.max_bound(2)) ! for y direction
             j=j+1
             k=0
             internal = .true.
             do while(point_xyz(1).le.max_bound(1)) ! for x direction
                k=k+1
                internal = point_internal_to_surface(num_vertices,triangle,point_xyz,vertex_xyz)
                if(internal)then
                  num_data = num_data+1
                  if(num_data.le.num_data_estimate)then
                     data_xyz(:,num_data) = point_xyz
                  else
                     num_data_estimate = num_data_estimate + 1000
                     allocate(data_temp(3,num_data-1))
                     data_temp = data_xyz ! copy to temporary array
                     deallocate(data_xyz) !deallocate initially allocated memory
                     allocate(data_xyz(3,num_data_estimate))
                     data_xyz(:,1:num_data-1) = data_temp(:,1:num_data-1)
                     deallocate(data_temp) !deallocate the temporary array

                     write(*,'('' WARNING: number of data is'',I6,''; increased array size'')') num_data
                     data_xyz(:,num_data) = point_xyz
                  endif
                endif
                ! increment the location in x direction by 'spacing'
                point_xyz(1) = point_xyz(1) + spacing
             enddo
             ! increment the location in y direction by 'spacing'
             point_xyz(1) = min_bound(1) + 0.5_dp*spacing
             point_xyz(2) = point_xyz(2) + spacing
         enddo
         ! increment the location in z direction by 'spacing'
         point_xyz(1:2) = min_bound(1:2) + 0.5_dp*spacing
         point_xyz(3) = point_xyz(3) + spacing
      enddo

      write(*,'('' Made'',I7,'' data points inside surface elements'')') num_data

      if(allocated(data_weight)) deallocate(data_weight)
         allocate(data_weight(3,num_data))
         data_weight(:,1:num_data) = 1.0_dp
      if(offset .gt. 0.0_dp)then
         !!! return the triangles to their original size and location
         forall (i=1:3) vertex_xyz(1:3,i) = (vertex_xyz(1:3,i) - translate(1:3))/scale_mesh
      endif

      if(to_export)then
         !!! export vertices as nodes
         writefile = trim(filename)//'.exnode'
         open(10, file = writefile, status='replace')
         !**    write the group name
         write(10,'( '' Group name: '',A)') trim(groupname)
         !*** Exporting Geometry
         !*** Write the field information
         write(10,'( '' #Fields=1'' )')
         write(10,'('' 1) coordinates, coordinate, rectangular cartesian, #Components=3'')')
         do nj=1,3
            if(nj.eq.1) write(10,'(2X,''x.  '')',advance="no")
            if(nj.eq.2) write(10,'(2X,''y.  '')',advance="no")
            if(nj.eq.3) write(10,'(2X,''z.  '')',advance="no")
            write(10,'(''Value index='',I2,'', #Derivatives='',I1)',advance="no") nj,0
            write(10,'()')
         enddo
         do i = 1,num_vertices
           !***    write the node
           write(10,'(1X,''Node: '',I12)') i
           write(10,'(2x,3(f12.6))') vertex_xyz(:,i)
         enddo
         close(10)

      !!! export the triangles as surface elements
      writefile = trim(filename)//'.exelem'
      open(10, file=writefile, status='replace')
      !**     write the group name
      write(10,'( '' Group name: '',a10)') groupname
      !**     write the lines
      write(10,'( '' Shape. Dimension=1'' )')
      nline=0
      do ne = 1,num_triangles
        write(10,'( '' Element: 0 0 '',I5)') nline+1
        write(10,'( '' Element: 0 0 '',I5)') nline+2
        write(10,'( '' Element: 0 0 '',I5)') nline+3
        nline=nline+3
      enddo !ne

      !**        write the elements
      write(10,'( '' Shape. Dimension=2, line*line'' )')
      write(10,'( '' #Scale factor sets=1'' )')
      write(10,'( '' l.Lagrange*l.Lagrange, #Scale factors=4'' )')
      write(10,'( '' #Nodes= '',I2 )') 4
      write(10,'( '' #Fields=1'' )')
      write(10,'( '' 1) coordinates, coordinate, rectangular cartesian, #Components=3'')')

      do nj=1,3
        if(nj==1) char1='x'; if(nj==2) char1='y'; if(nj==3) char1='z';
        write(10,'(''  '',A2,''. l.Lagrange*l.Lagrange, no modify, standard node based.'')') char1
        write(10,'( ''   #Nodes=4'')')
        do nn=1,4
          write(10,'(''   '',I1,''. #Values=1'')') nn
          write(10,'(''     Value indices: '',I4)') 1
          write(10,'(''     Scale factor indices: '',I4)') nn
        enddo !nn
      enddo !nj

      nline=0
      do ne=1,num_triangles
        !**         write the element
        write(10,'(1X,''Element: '',I12,'' 0 0'' )') ne
        !**          write the faces
        WRITE(10,'(3X,''Faces: '' )')
        WRITE(10,'(5X,''0 0'',I6)') nline+1
        WRITE(10,'(5X,''0 0'',I6)') nline+2
        WRITE(10,'(5X,''0 0'',I6)') nline+3
        WRITE(10,'(5X,''0 0'',I6)') 0
        nline=nline+3
        !**          write the nodes
        write(10,'(3X,''Nodes:'' )')
        write(10,'(4X,4(1X,I12))') triangle(:,ne),triangle(3,ne)
        !**          write the scale factors
        write(10,'(3X,''Scale factors:'' )')
        write(10,'(4X,4(1X,E12.5))') 1.0_dp,1.0_dp,1.0_dp,1.0_dp
      enddo
      close(10)
      endif

      deallocate(triangle)
      deallocate(vertex_xyz)

      call enter_exit(sub_name,2)

      end subroutine make_data_grid
! 
!###################################################################################
! 
  subroutine define_rad_from_file(FIELDFILE, radius_type_in)
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_DEFINE_RAD_FROM_FILE" :: DEFINE_RAD_FROM_FILE

  !*define_rad_from_file:* reads in a radius field associated with an aiway tree
! and assigns radius information to each element, also calculates volume of each
! element

    character(len=MAX_FILENAME_LEN), intent(in) :: FIELDFILE
    character(len=MAX_STRING_LEN), optional ::  radius_type_in
    !     Local Variables
    integer :: ierror,ne,np,np1,np2,np_global,surround
    character(len=MAX_STRING_LEN) ::  radius_type
    character(LEN=132) :: ctemp1
    LOGICAL :: versions
    real(dp) :: constrict,radius
    character(len=60) :: sub_name

    sub_name = 'define_rad_from_file'
    call enter_exit(sub_name,1)
    
    versions = .TRUE.
    if(present(radius_type_in))then
      radius_type = radius_type_in
    else
      radius_type = 'no_taper'
    endif

!!! note that 'constrict' should not be used here (so is set to 1.0).
!!! this should be specified and used as part of a simulation, not when
!!! reading in airway geometry
    constrict = 1.0_dp
    open(10, file=FIELDFILE, status='old')

    !.....check whether versions are prompted (>1)
    read_versions : do !define a do loop name
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !read a line into ctemp1
       if(index(ctemp1, "different")> 0) then !keyword "different" is found
          if(index(ctemp1, " N")> 0) then !keyword " N" is found
             versions=.false.
          endif
          exit read_versions !exit the named do loop
       endif
    end do read_versions

    np = 0
    !.....read the coordinate, derivative, and version information for each node.
    read_a_node : do !define a do loop name
       !.......read node number
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1
       if(index(ctemp1, "Node")> 0) then
          np_global = get_final_integer(ctemp1) !get global node number
          ! find the corresponding local node number
          call get_local_node(np_global,np) ! get local node np for global node
          surround=elems_at_node(np,0)         !get number of surrounding elems
          ne=elems_at_node(np,1)  !First element at this node
          if(surround==1)then !only one element at this node so either a terminal or inlet
             if(radius_type.eq.'taper')then !inlet radius needs to be defined
               if(elem_cnct(-1,0,ne).eq.0)then!Inlet as it has no parent need to set up radius into this vessel
                 read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                 read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                 if(index(ctemp1, "version number")>0) then
                    read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                 endif
                 if(index(ctemp1, "value")> 0) then
                    radius = get_final_real(ctemp1)
                    elem_field(ne_radius_in,ne)=constrict*radius
                 endif
               endif
             endif
             if(elem_cnct(-1,0,ne).eq.0)cycle      !No parent therefore inlet. Skip and go to the next
             read(unit=10, fmt="(a)", iostat=ierror) ctemp1
             read(unit=10, fmt="(a)", iostat=ierror) ctemp1
             if(index(ctemp1, "version number")>0) then
                read(unit=10, fmt="(a)", iostat=ierror) ctemp1
             endif
             if(index(ctemp1, "value")> 0) then
                radius = get_final_real(ctemp1)
                 if(radius_type.eq.'taper')then
                   elem_field(ne_radius_out,ne)=constrict*radius
                 else
                  elem_field(ne_radius,ne)=constrict*radius
                 endif
             endif
          elseif(surround.gt.1)then !Terminal airway - use first radius
             read(unit=10, fmt="(a)", iostat=ierror) ctemp1
             read(unit=10, fmt="(a)", iostat=ierror) ctemp1
             read(unit=10, fmt="(a)", iostat=ierror) ctemp1
             if(index(ctemp1, "value")> 0) then
                radius = get_final_real(ctemp1)
                  if(radius_type.eq.'taper')then
                    elem_field(ne_radius_out,ne)=constrict*radius
                  else
                    elem_field(ne_radius,ne)=constrict*radius
                  endif
             endif
          endif
       endif !index
       if(np.ge.num_nodes) exit read_a_node
    end do read_a_node

!If airway type calculate element volume
    ! calculate the element volumes
    do ne=1,num_elems
       if(radius_type.eq.'taper')then
         if(elem_cnct(-1,0,ne).ne.0)then !radius in is radius of upstream vessel
            elem_field(ne_radius_in,ne)=elem_field(ne_radius_out,elem_cnct(-1,1,ne))
         endif
         elem_field(ne_radius,ne)=(elem_field(ne_radius_in,ne)+elem_field(ne_radius_out,ne))/2
       endif
       !       ne_global=elems(noelem)
       np1=elem_nodes(1,ne)
       np2=elem_nodes(2,ne)
       ! element volume
       elem_field(ne_vol,ne) = PI * elem_field(ne_radius,ne)**2 * &
            elem_field(ne_length,ne)
       elem_field(ne_a_A,ne) = 1.0_dp ! set default for ratio a/A
    enddo
    call enter_exit(sub_name,2)

  END subroutine define_rad_from_file
!
!##################################################################################
!
!*define_rad_from_geom:* Defines vessel or airway radius based on their geometric structure
  subroutine define_rad_from_geom(ORDER_SYSTEM, CONTROL_PARAM, START_FROM, START_RAD, group_type_in, group_option_in)
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_DEFINE_RAD_FROM_GEOM" :: DEFINE_RAD_FROM_GEOM

    character(LEN=100), intent(in) :: ORDER_SYSTEM,START_FROM
   character(LEN=100), optional :: group_type_in, group_option_in
   real(dp), intent(in) :: CONTROL_PARAM,START_RAD
   !Input options ORDER_SYSTEM=STRAHLER (CONTROL_PARAM=RDS), HORSFIELD (CONTROL_PARAM=RDH)

   !Local variables
   character(LEN=100) :: group_type, group_options
   integer :: ne_min,ne_max,nindex,ne,n_max_ord,n,ne_start,&
      inlet_count
   real(dp) :: radius
   character(len=60) :: sub_name

   sub_name = 'define_rad_from_geom'
    call enter_exit(sub_name,1)
    !define list of elements you are going to operate on
    if(present(group_type_in))then
      group_type = group_type_in
    else!default to all
      group_type='all'
    endif
    if(group_type.eq.'all')then
       ne_min=1
       ne_max=num_elems
    elseif(group_type.eq.'efield')then
    elseif(group_type.eq.'list')then
      read (START_FROM,'(I10)') ne_min
      read (group_option_in,'(I10)') ne_max
    endif
    !Define start element
    if(START_FROM.eq.'inlet')then
      inlet_count=0
      do ne=ne_min,ne_max
         if(elem_cnct(-1,0,ne).eq.0)then
           inlet_count=inlet_count+1
           ne_start=ne
         endif
         if(inlet_count.gt.1)then
            WRITE(*,*) ' More than one inlet in this group, using last found, ne = ',ne
         endif
      enddo
    else!element number defined
       read (START_FROM,'(I10)') ne_start
    endif

    !Strahler and Horsfield ordering system
    if(ORDER_SYSTEM(1:5).EQ.'strah')THEN
      nindex = no_sord !for Strahler ordering
    else if(ORDER_SYSTEM(1:5).eq.'horsf')then
      nindex = no_hord !for Horsfield ordering
    endif

    ne=ne_start
    n_max_ord=elem_ordrs(nindex,ne)
    elem_field(ne_radius,ne)=START_RAD

    do ne=ne_min,ne_max
     radius=10.0_dp**(log10(CONTROL_PARAM)*dble(elem_ordrs(nindex,ne)-n_max_ord)&
        +log10(START_RAD))
     elem_field(ne_radius,ne)=radius
     if(ne_radius_in.gt.0)then
        elem_field(ne_radius_in,ne)=radius
        elem_field(ne_radius_out,ne)=radius
     endif
    enddo

    call enter_exit(sub_name,2)

  END subroutine define_rad_from_geom
!
!###########################################################################
!
!*element_connectivity_1d:*  Calculates element connectivity in 1D and stores in elelem_cnct
  subroutine element_connectivity_1d()
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_ELEMENT_CONNECTIVITY_1D" :: ELEMENT_CONNECTIVITY_1D

    !     Local Variables
    integer :: ne,ne2,nn,noelem,np,np2,np1
    integer,parameter :: NNT=2
    character(len=60) :: sub_name

    sub_name = 'element_connectivity_1d'
    call enter_exit(sub_name,1)

    elem_cnct = 0 !initialise

    ! calculate elems_at_node array: stores the elements that nodes are in
    elems_at_node(1:num_nodes,0) = 0 !initialise number of adjacent elements

    DO ne=1,num_elems
       DO nn=1,2
          np=elem_nodes(nn,ne)
          elems_at_node(np,0)=elems_at_node(np,0)+1
          elems_at_node(np,elems_at_node(np,0))=ne ! local element that np is in
       ENDDO !nn
    ENDDO !noelem

    ! calculate elem_cnct array: stores the connectivity of all elements

    elem_cnct=0 !initialise all elem_cnct

    DO ne=1,num_elems
       !     ne_global=elems(noelem)
       IF(NNT == 2) THEN !1d
          np1=elem_nodes(1,ne) !first local node
          np2=elem_nodes(2,ne) !second local node
          DO noelem=1,elems_at_node(np2,0)
             ne2=elems_at_node(np2,noelem)
             IF(ne2 /= ne)THEN
                elem_cnct(-1,0,ne2)=elem_cnct(-1,0,ne2)+1
                elem_cnct(-1,elem_cnct(-1,0,ne2),ne2)=ne !previous element
                elem_cnct(1,0,ne)=elem_cnct(1,0,ne)+1
                elem_cnct(1,elem_cnct(1,0,ne),ne)=ne2
             ENDIF !ne2
          ENDDO !noelem2

       ENDIF
    ENDDO

    call enter_exit(sub_name,2)

  END subroutine element_connectivity_1d


!
!###########################################################################
!

  subroutine element_connectivity_2d

!!! calculates the connectivity of elements in each xi direction

!!! local variables
    integer :: ne,ne2,nn,np,noelem2,np_list(4),np_list_2(4)
    integer,parameter :: num_elem_nodes = 4
    
    if(.not.allocated(elem_cnct_2d)) allocate(elem_cnct_2d(-2:2,0:10,num_elems_2d))
    if(.not.allocated(elems_at_node_2d)) allocate(elems_at_node_2d(num_nodes_2d,0:10))

    ! calculate elems_at_node_2d array: stores the elements that nodes are in
    
    elems_at_node_2d = 0 !initialise all
    
    do ne = 1,num_elems_2d
       do nn = 1,num_elem_nodes
          np = elem_nodes_2d(nn,ne)
          elems_at_node_2d(np,0) = elems_at_node_2d(np,0)+1
          elems_at_node_2d(np,elems_at_node_2d(np,0)) = ne !element that np is in
       enddo !nn
    enddo !noelem
    
    ! calculate elem_cnct_2d array: stores the connectivity of all elements
    
    elem_cnct_2d = 0 !initialise all elem_cnct_2d
    
    do ne = 1,num_elems_2d ! for each of the 2d elements
       np_list(1:4) = elem_nodes_2d(1:4,ne) ! the list of nodes in the element (including repeated)
!!! check the elements attached to the 1st node
       do noelem2 = 1,elems_at_node_2d(np_list(1),0) ! for each element attached to the 1st node
          ne2 = elems_at_node_2d(np_list(1),noelem2) ! attached element number
          if(ne2.ne.ne)then
             np_list_2(1:4) = elem_nodes_2d(1:4,ne2) !list of nodes in attached element
             if(np_list(2).ne.np_list(1))then ! only if first two nodes are not repeated
                if(inlist(np_list(2),np_list_2))then
                   elem_cnct_2d(-2,0,ne) = elem_cnct_2d(-2,0,ne)+1
                   elem_cnct_2d(-2,elem_cnct_2d(-2,0,ne),ne) = ne2 
                endif
             endif
             if(np_list(3).ne.np_list(1))then ! only if the two nodes are not repeated
                if(inlist(np_list(3),np_list_2))then
                   elem_cnct_2d(-1,0,ne) = elem_cnct_2d(-1,0,ne)+1
                   elem_cnct_2d(-1,elem_cnct_2d(-1,0,ne),ne) = ne2 
                endif
             endif
          endif
       enddo
!!! check the elements attached to the 4th node
       do noelem2 = 1,elems_at_node_2d(np_list(4),0) ! for each element attached to the 4th node
          ne2 = elems_at_node_2d(np_list(4),noelem2) ! attached element number
          if(ne2.ne.ne)then
             np_list_2(1:4) = elem_nodes_2d(1:4,ne2) !list of nodes in attached element
             if(np_list(2).ne.np_list(4))then ! only if two nodes are not repeated
                if(inlist(np_list(2),np_list_2))then
                   elem_cnct_2d(2,0,ne) = elem_cnct_2d(2,0,ne)+1
                   elem_cnct_2d(2,elem_cnct_2d(2,0,ne),ne) = ne2 
                endif
             endif
             if(np_list(3).ne.np_list(4))then ! only if the two nodes are not repeated
                if(inlist(np_list(3),np_list_2))then
                   elem_cnct_2d(1,0,ne) = elem_cnct_2d(1,0,ne)+1
                   elem_cnct_2d(1,elem_cnct_2d(1,0,ne),ne) = ne2 
                endif
             endif
          endif
       enddo !noelem2
    enddo

  end subroutine element_connectivity_2d


!!! ##########################################################################      

  subroutine line_segments_for_2d_mesh(sf_option)

!!! sets up the line segment arrays for a 2d mesh
    
    character(len=4),intent(in) :: sf_option
!!! local variables
    integer :: ne,ne_adjacent,ni1,nj,nl,nl_adj,npn(2)
    logical :: MAKE
    logical :: based_on_elems = .true.
    
    if(.not.allocated(elem_lines_2d)) allocate(elem_lines_2d(4,num_elems_2d))
    if(.not.allocated(scale_factors_2d)) allocate(scale_factors_2d(16,num_elems_2d))

    elem_lines_2d=0
    num_lines_2d = 0

    if(based_on_elems)then
!!! estimate number of lines, for allocating memory to arrays
!!! before setting up arrays, count the number of lines required
       do ne=1,num_elems_2d
          MAKE=.FALSE.
          if(elem_cnct_2d(-1,0,ne) == 0) MAKE=.TRUE. !exterior, make line
          ne_adjacent=elem_cnct_2d(-1,1,ne)
          if(ne_adjacent > 0)then
             if(elem_lines_2d(4,ne_adjacent) == 0) MAKE=.TRUE.
          endif
          if(MAKE) num_lines_2d = num_lines_2d+1
          MAKE=.FALSE.
          if(elem_cnct_2d(-2,0,ne) == 0) MAKE=.TRUE. !exterior, make line
          ne_adjacent=elem_cnct_2d(-2,1,ne)
          if(ne_adjacent > 0)then
             if(elem_lines_2d(2,ne_adjacent) == 0) MAKE=.TRUE.
          endif
          if(MAKE) num_lines_2d=num_lines_2d+1
          num_lines_2d = num_lines_2d+2
          elem_lines_2d(2,ne) = 1 ! at this stage just to tag it for conditional above
          elem_lines_2d(4,ne) = 1 ! at this stage just to tag it for conditional above
       enddo !ne
       
       elem_lines_2d = 0
       
       if(.not.allocated(lines_2d)) allocate(lines_2d(0:num_lines_2d))
       if(.not.allocated(line_versn_2d)) allocate(line_versn_2d(2,3,num_lines_2d))
       if(.not.allocated(lines_in_elem)) allocate(lines_in_elem(0:4,num_lines_2d))
       if(.not.allocated(nodes_in_line)) allocate(nodes_in_line(3,0:3,num_lines_2d))
       if(.not.allocated(arclength)) allocate(arclength(3,num_lines_2d)) 
       lines_in_elem=0
       lines_2d=0
       nodes_in_line=0
       line_versn_2d=0
       num_lines_2d = 0 ! reset to zero for loop below
       
!!! Now run through the same as above, and set up the arrays
       do ne=1,num_elems_2d
          !check whether to make a line
          MAKE=.FALSE.
          if(elem_cnct_2d(-1,0,ne) == 0) MAKE=.TRUE. !exterior, make line
          ne_adjacent=elem_cnct_2d(-1,1,ne)
          if(ne_adjacent.gt.0)then
             if(elem_lines_2d(4,ne_adjacent) == 0) MAKE=.TRUE.
          endif
          if(MAKE)then
             num_lines_2d = num_lines_2d+1
             lines_2d(num_lines_2d) = num_lines_2d !record a new line number
             lines_in_elem(0,num_lines_2d) = lines_in_elem(0,num_lines_2d)+1
             lines_in_elem(lines_in_elem(0,num_lines_2d),num_lines_2d) = ne !line num_lines_2d is in element ne
             elem_lines_2d(3,ne) = num_lines_2d !num_lines_2d is global line # corresponding to local line 3 of ne
             npn(1) = 1
             npn(2) = 3
             nodes_in_line(2,1,num_lines_2d)=elem_nodes_2d(1,ne) !records 1st node in line
             nodes_in_line(3,1,num_lines_2d)=elem_nodes_2d(3,ne) !records 2nd node in line
             nodes_in_line(1,0,num_lines_2d)=2 !Xi-direction of line segment num_lines_2d
             do nj=1,3
                nodes_in_line(1,nj,num_lines_2d)=4 !type of basis function (1 for linear,4 for cubicHermite)
                do ni1=1,2
                   line_versn_2d(ni1,nj,num_lines_2d)=elem_versn_2d(npn(ni1),ne)
                enddo !n
             enddo !nj
          else !get adjacent element line number
             !WARNING:: this only works if all Xi directions are consistent!!!!
             ne_adjacent=elem_cnct_2d(-1,1,ne)
             elem_lines_2d(3,ne)=elem_lines_2d(4,ne_adjacent)
          endif
          
          !check whether to make a line
          MAKE=.FALSE.
          if(elem_cnct_2d(-2,0,ne) == 0) MAKE=.TRUE. !exterior, make line
          ne_adjacent=elem_cnct_2d(-2,1,ne)
          if(ne_adjacent.gt.0)then
             if(elem_lines_2d(2,ne_adjacent) == 0) MAKE=.TRUE.
          endif

          if(MAKE)then
             num_lines_2d=num_lines_2d+1
             lines_2d(num_lines_2d)=num_lines_2d !record a new line number
             lines_in_elem(0,num_lines_2d)=lines_in_elem(0,num_lines_2d)+1
             lines_in_elem(lines_in_elem(0,num_lines_2d),num_lines_2d)=ne !line num_lines_2d is in element ne
             elem_lines_2d(1,ne)=num_lines_2d !num_lines_2d is global line # corresponding to local line 1 of ne
             npn(1)=1
             npn(2)=2
             nodes_in_line(2,1,num_lines_2d)=elem_nodes_2d(1,ne) !records 1st node in line
             nodes_in_line(3,1,num_lines_2d)=elem_nodes_2d(2,ne) !records 2nd node in line
             nodes_in_line(1,0,num_lines_2d)=1 !Xi-direction of line segment num_lines_2d
             do nj=1,3
                nodes_in_line(1,nj,num_lines_2d)=4 !type of basis function (1 for linear,4 for cubicHermite)
                do ni1=1,2
                   line_versn_2d(ni1,nj,num_lines_2d)=elem_versn_2d(npn(ni1),ne)
                enddo !n
             enddo !nj
          else !get adjacent element line number
             !WARNING:: this only works if all Xi directions are consistent!!!!
             ne_adjacent = elem_cnct_2d(-2,1,ne)
             do nl_adj = 1,4
                nl = elem_lines_2d(nl_adj,ne_adjacent)
                if(nl /= 0)then
                   if(nodes_in_line(2,1,nl) == elem_nodes_2d(1,ne) .and. &
                        nodes_in_line(3,1,nl) == elem_nodes_2d(2,ne))then
                      elem_lines_2d(1,ne) = nl
                   elseif(nodes_in_line(2,1,nl) == elem_nodes_2d(2,ne) .and. &
                        nodes_in_line(3,1,nl) == elem_nodes_2d(1,ne))then
                      elem_lines_2d(1,ne) = nl
                   endif
                endif
             enddo
!             elem_lines_2d(1,ne)=elem_lines_2d(2,ne_adjacent)
          endif
          
          !*! new:       
          MAKE=.TRUE.
          ne_adjacent=elem_cnct_2d(1,1,ne)
          if(ne_adjacent.gt.0)then
             if(elem_lines_2d(3,ne_adjacent) /= 0) MAKE=.FALSE.
          endif
          
          if(MAKE)then
             num_lines_2d=num_lines_2d+1
             lines_2d(num_lines_2d)=num_lines_2d !record a new line number
             lines_in_elem(0,num_lines_2d)=lines_in_elem(0,num_lines_2d)+1
             lines_in_elem(lines_in_elem(0,num_lines_2d),num_lines_2d)=ne !line num_lines_2d is in element ne
             elem_lines_2d(4,ne)=num_lines_2d !num_lines_2d is global line # corresponding to local line 4 of ne
             npn(1)=2
             npn(2)=4
             nodes_in_line(2,1,num_lines_2d)=elem_nodes_2d(2,ne) !records 1st node in line
             nodes_in_line(3,1,num_lines_2d)=elem_nodes_2d(4,ne) !records 2nd node in line
             nodes_in_line(1,0,num_lines_2d)=2 !Xi-direction of line segment num_lines_2d
             do nj=1,3
                nodes_in_line(1,nj,num_lines_2d)=4 !type of basis function (1 for linear,4 for cubicHermite)
                do ni1=1,2
                   line_versn_2d(ni1,nj,num_lines_2d)=elem_versn_2d(npn(ni1),ne)
                enddo !n
             enddo !nj
          else !get adjacent element line number
             !WARNING:: this only works if all Xi directions are consistent!!!!
             ne_adjacent=elem_cnct_2d(1,1,ne)
             elem_lines_2d(4,ne)=elem_lines_2d(3,ne_adjacent)
          endif
          
          MAKE=.TRUE.
          ne_adjacent=elem_cnct_2d(2,1,ne)
          if(ne_adjacent.gt.0)then
             if(elem_lines_2d(1,ne_adjacent) /= 0) MAKE=.FALSE.
          endif
          
          if(MAKE)then
             num_lines_2d = num_lines_2d+1
             lines_2d(num_lines_2d) = num_lines_2d !record a new line number
             lines_in_elem(0,num_lines_2d) = lines_in_elem(0,num_lines_2d)+1
             lines_in_elem(lines_in_elem(0,num_lines_2d),num_lines_2d) = ne !line num_lines_2d is in element ne
             elem_lines_2d(2,ne)=num_lines_2d !num_lines_2d is global line # corresponding to local line 2 of ne
             npn(1) = 3
             npn(2) = 4
             nodes_in_line(2,1,num_lines_2d)=elem_nodes_2d(3,ne) !records 1st node in line
             nodes_in_line(3,1,num_lines_2d)=elem_nodes_2d(4,ne) !records 2nd node in line
             nodes_in_line(1,0,num_lines_2d)=1 !Xi-direction of line segment num_lines_2d
             do nj=1,3
                nodes_in_line(1,nj,num_lines_2d)=4 !type of basis function (1 for linear,4 for cubicHermite)
                do ni1=1,2
                   line_versn_2d(ni1,nj,num_lines_2d)=elem_versn_2d(npn(ni1),ne)
                enddo !n
             enddo !nj
          else !get adjacent element line number
             !WARNING:: this only works if all Xi directions are consistent!!!!
             ne_adjacent=elem_cnct_2d(2,1,ne)
             elem_lines_2d(2,ne)=elem_lines_2d(1,ne_adjacent)
          endif
          
       enddo !ne
    endif

    call calc_scale_factors_2d(sf_option)
    
  end subroutine line_segments_for_2d_mesh


!
!###################################################################################
!
!*evaluate_ordering:* calculates generations, Horsfield orders, Strahler orders for a given tree
  subroutine evaluate_ordering()
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_EVALUATE_ORDERING" :: EVALUATE_ORDERING

    integer :: INLETS,ne,ne0,ne2,noelem2,np,np2,nn, &
         num_attach,n_children,n_generation, &
         n_horsfield,OUTLETS,STRAHLER,STRAHLER_ADD,temp1
    LOGICAL :: DISCONNECT,DUPLICATE
    character(len=60) :: sub_name
    sub_name = 'evaluate_ordering'
    call enter_exit(sub_name,1)

    !Calculate generations, Horsfield orders, Strahler orders
    !.....Calculate branch generations

    ! Initialise memory to zero.
    do nn=1,num_ord
      do ne=1,num_elems
        elem_ordrs(nn,ne) = 0
      enddo
    enddo

    maxgen=1
    DO ne=1,num_elems
       ne0=elem_cnct(-1,1,ne) !parent
       IF(ne0.NE.0)THEN
          n_generation=elem_ordrs(1,ne0) !parent generation
          IF(elem_cnct(1,0,ne0).EQ.1)THEN !single daughter
             elem_ordrs(1,ne)=n_generation + (elem_symmetry(ne)-1)
          ELSE IF(elem_cnct(1,0,ne0).GE.2)THEN
             elem_ordrs(1,ne)=n_generation+1
          ENDIF
       ELSE
          elem_ordrs(1,ne)=1 !generation 1
       ENDIF
       maxgen=max(maxgen,elem_ordrs(1,ne))
    ENDDO !noelem

    !.....Calculate the branch orders
    DO ne=num_elems,1,-1
       n_horsfield=MAX(elem_ordrs(2,ne),1)
       n_children=elem_cnct(1,0,ne) !number of child branches
       IF(n_children.EQ.1)THEN
          IF(elem_ordrs(1,elem_cnct(1,1,ne)).EQ.0)  n_children=0
       ENDIF
       STRAHLER=0
       STRAHLER_ADD=1
       IF(n_children.GE.2)THEN !branch has two or more daughters
          STRAHLER=elem_ordrs(3,elem_cnct(1,1,ne)) !first daughter
          DO noelem2=1,n_children !for all daughters
             ne2=elem_cnct(1,noelem2,ne) !global element # of daughter
             temp1=elem_ordrs(2,ne2) !Horsfield order of daughter
             IF(temp1.GT.n_horsfield)then
               n_horsfield=temp1
             endif
             IF(elem_ordrs(3,ne2).LT.STRAHLER)THEN
                STRAHLER_ADD=0
             ELSE IF(elem_ordrs(3,ne2).GT.STRAHLER)THEN
                STRAHLER_ADD=0
                STRAHLER=elem_ordrs(3,ne2) !highest daughter
             ENDIF
          ENDDO !noelem2 (ne2)
          n_horsfield=n_horsfield+1 !Horsfield ordering
       ELSE IF(n_children.EQ.1)THEN
          ne2=elem_cnct(1,1,ne) !local element # of daughter
          n_horsfield=elem_ordrs(2,ne2)+(elem_symmetry(ne)-1)

          STRAHLER_ADD=elem_ordrs(3,ne2)+(elem_symmetry(ne)-1)
       ENDIF !elem_cnct
       elem_ordrs(2,ne)=n_horsfield !store the Horsfield order
       elem_ordrs(3,ne)=STRAHLER+STRAHLER_ADD !Strahler order
    ENDDO !noelem

    !       Check for disconnected nodes and number of inlets and outlets
    DUPLICATE=.FALSE.
    DO ne=1,num_elems
       np=elem_nodes(1,ne)
       np2=elem_nodes(2,ne)
       IF(np.EQ.np2)THEN
          DUPLICATE=.TRUE.
       ENDIF
    ENDDO

    DISCONNECT=.FALSE.
    INLETS=0
    OUTLETS=0
    DO np=1,num_nodes
       num_attach=elems_at_node(np,0)
       IF(num_attach.EQ.0)THEN
          DISCONNECT=.TRUE.
       ELSEIF(num_attach.EQ.1)THEN
          ne=elems_at_node(np,1)
          IF(elem_cnct(1,0,ne).EQ.0) OUTLETS=OUTLETS+1
         IF(elem_cnct(-1,0,ne).EQ.0) INLETS=INLETS+1
       ELSEIF(num_attach.GT.3)THEN
          WRITE(*,*) ' Node ',np,' attached to',num_attach,' elements'
       ENDIF
    ENDDO

    call enter_exit(sub_name,2)

  end subroutine evaluate_ordering
!
!###################################################################################
!
!>*set_initial_volume:* assigns a volume to terminal units appended on a tree structure
!>based on an assumption of a linear gradient in the gravitational direction with max
!> min and COV values defined.
  subroutine set_initial_volume(Gdirn,COV,total_volume,Rmax,Rmin)
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_SET_INITIAL_VOLUME" :: SET_INITIAL_VOLUME

    !     Parameter List
    integer,intent(in) :: Gdirn
    real(dp),intent(in) :: COV,total_volume,Rmax,Rmin
    !     Local parameters
    integer :: ne,np2,nunit
    real(dp) ::  factor_adjust,max_z,min_z,random_number,range_z,&
         volume_estimate,volume_of_tree,Vmax,Vmin,Xi
    character(len=60) :: sub_name

    sub_name = 'set_initial_volume'
    call enter_exit(sub_name,1)
    
    volume_estimate = 1.0_dp
    volume_of_tree = 0.0_dp

    call volume_of_mesh(volume_estimate,volume_of_tree)

    random_number=-1.1_dp

    Vmax = Rmax * (total_volume-volume_estimate)/elem_units_below(1)
    Vmin = Rmin * (total_volume-volume_estimate)/elem_units_below(1)

!!! for each elastic unit find the maximum and minimum coordinates in the Gdirn direction
    max_z=-1.0e+6_dp
    min_z=1.0e+6_dp
    do nunit=1,num_units
       ne=units(nunit)
       np2=elem_nodes(2,ne)
       max_z=MAX(max_z,node_xyz(Gdirn,np2))
       min_z=MIN(min_z,node_xyz(Gdirn,np2))
    enddo !nunit

    range_z=DABS(max_z-min_z)
    if(DABS(range_z).le.1.0e-5_dp) range_z=1.0_dp

!!! for each elastic unit allocate a size based on a gradient in the Gdirn direction, and
!!! perturb by a user-defined COV. This should be calling a random number generator.
    do nunit=1,num_units
       ne=units(nunit)
       np2=elem_nodes(2,ne) !end node
       Xi=(node_xyz(Gdirn,np2)-min_z)/range_z
       random_number=random_number+0.1_dp
       IF(random_number.GT.1.0_dp) random_number=-1.1_dp
       unit_field(nu_vol,nunit)=(Vmax*Xi+Vmin*(1.0_dp-Xi))*(1.0_dp+COV*random_number)
       unit_field(nu_vt,nunit)=0.0_dp !initialise the tidal volume to a unit
    enddo !nunit

    ! correct unit volumes such that total volume is exactly as specified
    call volume_of_mesh(volume_estimate,volume_of_tree)
    factor_adjust = (total_volume-volume_of_tree)/(volume_estimate-volume_of_tree)
    do nunit=1,num_units
       unit_field(nu_vol,nunit) = unit_field(nu_vol,nunit)*factor_adjust
    enddo

    write(*,'('' Number of elements is '',I5)') num_elems
    write(*,'('' Initial volume is '',F6.2,'' L'')') total_volume/1.0e+6_dp
    write(*,'('' Deadspace volume is '',F6.1,'' mL'')') volume_of_tree/1.0e+3_dp

    call enter_exit(sub_name,2)

  end subroutine set_initial_volume

!
!###################################################################################
!
!*volume_of_mesh:* calculates the volume of an airway mesh including conducting and respiratory airways
  subroutine volume_of_mesh(volume_model,volume_tree)
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_VOLUME_OF_MESH" :: VOLUME_OF_MESH

    real(dp) :: volume_model,volume_tree
    !     Local Variables
    integer :: ne,ne0,nunit
    real(dp),allocatable :: vol_anat(:),vol_below(:)
    character(len=60) :: sub_name

    sub_name = 'volume_of_mesh'
    call enter_exit(sub_name,1)

    if(.not.allocated(vol_anat)) allocate(vol_anat(num_elems))
    if(.not.allocated(vol_below)) allocate(vol_below(num_elems))

    vol_anat = elem_field(ne_vol,1:num_elems) !initialise to branch volume
    vol_below = elem_field(ne_vol,1:num_elems) !initialise to branch volume

    do nunit=1,num_units
       ne=units(nunit)
       vol_below(ne) = vol_below(ne) + unit_field(nu_vol,nunit) !add elastic unit volume
    enddo !nunit

    do ne=num_elems,2,-1
       ne0=elem_cnct(-1,1,ne)
       vol_anat(ne0) = vol_anat(ne0) + DBLE(elem_symmetry(ne))*vol_anat(ne)
       vol_below(ne0) = vol_below(ne0) + DBLE(elem_symmetry(ne))*vol_below(ne)
    enddo !noelem

    volume_model = vol_below(1)
    volume_tree = vol_anat(1)

    deallocate(vol_anat)
    deallocate(vol_below)

    call enter_exit(sub_name,2)

  end subroutine volume_of_mesh

!
!###################################################################################
!
!>get_final_real
  
  function get_final_real(string)

    character :: string*(*)
    integer :: ibeg,iend
    real(dp) :: rsign,rtemp,get_final_real
    character :: sub_string*(40)
    
    iend=len(string) !get the length of the string
    ibeg=index(string,":")+1 !get location of real in string
    sub_string = adjustl(string(ibeg:iend)) ! get the characters beyond :
    iend=len(sub_string) !get the length of the sub-string
    if(sub_string(1:1).eq.'-')then !check whether negative
       rsign=-1.0_dp
       ibeg=2
    else
       rsign=1.0_dp
       ibeg=1
    endif
   
    read (sub_string(ibeg:iend), * ) rtemp !get real value
    rtemp=rtemp*rsign !apply sign to number
    
    get_final_real=rtemp !return the real value
    
  end function get_final_real

!
!###################################################################################
!
!*get_final_string*
  subroutine get_final_string(string,rtemp)

    character, intent(in) :: string*(132)
    integer :: ibeg,iend
    real(dp), intent(out) :: rtemp
    real(dp) :: rsign
    character :: sub_string*(40)

    iend=len(string) !get the length of the string
    ibeg=index(string,":")+1 !get location of real in string
    sub_string = adjustl(string(ibeg:iend)) ! get the characters beyond :
    iend=len(sub_string) !get the length of the sub-string
    if(sub_string(1:1).eq.'-')then !check whether negative
       rsign=-1.0_dp
       ibeg=2
    else
       rsign=1.0_dp
       ibeg=1
    endif
    read (sub_string(ibeg:iend), '(D25.17)' ) rtemp !get real value
    rtemp=rtemp*rsign !apply sign to number

  end subroutine get_final_string

!
!###################################################################################
!
!*get_local_node*
  subroutine get_local_node(np_global,np_local)

    integer,intent(in) :: np_global
    integer,intent(out) :: np_local

    integer :: np
    logical :: found

    np=1
    found=.false.
    do while (.not.found)
       if(nodes(np).eq.np_global)then
          found=.true.
       elseif(np.gt.num_nodes)then
          found = .true.
          write(*,'('' Global node '',I6,'' not in node list'')') np_global
          read(*,*)
       else
          np=np+1
       endif
    enddo

    np_local = np
    return

  end subroutine get_local_node


!
!###################################################################################
!
!*group_elem_parent_term*
   subroutine group_elem_parent_term(ne_parent)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_GROUP_ELEM_PARENT_TERM" :: GROUP_ELEM_PARENT_TERM

   ! Parameters
   integer,intent(in) :: ne_parent  ! Number of element that feed that specific area

   ! Local Variables
   integer :: ne,ne_count,noelem,num_parents
   integer,allocatable :: templist(:)
   character(len=60) :: sub_name='group_elem_parent_term'

   call enter_exit(sub_name,1)

   allocate(templist(num_elems))
   if(.not.allocated(parentlist)) allocate(parentlist(num_elems))

   !reset the list of parent elements to zero
   call group_elem_by_parent(ne_parent,templist)

   ne_count=count(templist.ne.0)
   num_parents=0
   parentlist=0

   do noelem=1,ne_count
     ne = templist(noelem)
     if(elem_cnct(1,0,ne).eq.0)then
        num_parents=num_parents+1
        parentlist(num_parents)=ne
     endif !elem_cnct
   enddo !noelem

   deallocate(templist)

   call enter_exit(sub_name,2)

   end subroutine group_elem_parent_term
!
!###################################################################################
!
!*group_elem_by_parent*
   subroutine group_elem_by_parent(ne_parent,elemlist)

   ! Parameters
   integer,intent(in) :: ne_parent  ! Number of element that feed that specific area
   integer :: elemlist(:)

   ! Local Variables
   integer :: nt_bns,ne_count,num_nodes,m,n,ne0
   integer,allocatable :: ne_old(:),ne_temp(:)
   character(len=60) :: sub_name='group_elem_by_parent'

   call enter_exit(sub_name,1)

   elemlist = 0
   allocate(ne_old(size(elemlist)))
   allocate(ne_temp(size(elemlist)))

   nt_bns=1
   ne_old(1) = ne_parent
   ne_count = 1
   elemlist(ne_count)=ne_parent

   do while(nt_bns.ne.0)
     num_nodes=nt_bns
     nt_bns=0
     do m=1,num_nodes
       ne0=ne_old(m) !parent global element number
       do n=1,elem_cnct(1,0,ne0) !for each daughter branch
         nt_bns=nt_bns+1
         ne_temp(nt_bns)=elem_cnct(1,n,ne0)
       enddo !n
     enddo !m
     do n=1,nt_bns
       ne_old(n)=ne_temp(n) !updates list of previous generation element numbers
       ne_count=ne_count+1
       elemlist(ne_count)=ne_temp(n)
     enddo !n
   enddo !while

   deallocate(ne_old)
   deallocate(ne_temp)

   call enter_exit(sub_name,2)

   end subroutine group_elem_by_parent

!
!###################################################################################
!
   
!
!#####################################################################
!
!*reallocate_node_elem_arrays:* Reallocates the size of geometric arrays when modifying geometries
  subroutine reallocate_node_elem_arrays(num_elems_new,num_nodes_new)

!!! Parameters
    integer,intent(in) :: num_elems_new,num_nodes_new

!!! Local variables
    integer,allocatable :: nodelem_temp(:),enodes_temp(:,:),enodes_temp2(:,:,:)
    real(dp),allocatable :: xyz_temp(:,:),rnodes_temp(:,:)
    logical,allocatable :: exp_temp(:)
    character(len=60) :: sub_name

    sub_name = 'reallocate_node_elem_arrays'
    call enter_exit(sub_name,1)

    allocate(nodelem_temp(num_nodes))
    nodelem_temp = nodes ! copy to temporary array
    deallocate(nodes) !deallocate initially allocated memory
    allocate(nodes(num_nodes_new))
    nodes(1:num_nodes)=nodelem_temp(1:num_nodes)
    deallocate(nodelem_temp) !deallocate the temporary array

    allocate(xyz_temp(3,num_nodes))
    xyz_temp=node_xyz
    deallocate(node_xyz)
    allocate(node_xyz(3,num_nodes_new))
    node_xyz(1:3,1:num_nodes)=xyz_temp(1:3,1:num_nodes)

    allocate(nodelem_temp(num_elems))
    nodelem_temp = elems ! copy to temporary array
    deallocate(elems) !deallocate initially allocated memory
    allocate(elems(num_elems_new))
    elems(1:num_elems)=nodelem_temp(1:num_elems)
    deallocate(nodelem_temp) !deallocate the temporary array

    allocate(enodes_temp(2,num_elems))
    enodes_temp=elem_nodes
    deallocate(elem_nodes)
    allocate(elem_nodes(2,num_elems_new))
    elem_nodes(1:2,1:num_elems)=enodes_temp(1:2,1:num_elems)
    deallocate(enodes_temp)

    if(allocated(elem_field).and.num_ne.gt.0)then
        allocate(rnodes_temp(num_ne,num_elems))
        rnodes_temp=elem_field
        deallocate(elem_field)
        allocate(elem_field(num_ne,num_elems_new))
        elem_field(1:num_ne,1:num_elems)=rnodes_temp(1:num_ne,1:num_elems)
        deallocate(rnodes_temp)
        elem_field(1:num_ne,num_elems+1:num_elems_new) = 0.0_dp
    endif

    allocate(rnodes_temp(3,num_elems))
    rnodes_temp=elem_direction
    deallocate(elem_direction)
    allocate(elem_direction(3,num_elems_new))
    elem_direction(1:3,1:num_elems)=rnodes_temp(1:3,1:num_elems)
    deallocate(rnodes_temp)
    elem_direction(1:3,num_elems+1:num_elems_new) = 0.0_dp

    if(allocated(node_field).and.num_nj.gt.0)then
      allocate(rnodes_temp(num_nj,num_nodes))
      rnodes_temp=node_field
      deallocate(node_field)
      allocate(node_field(num_nj,num_nodes_new))
      node_field(1:num_nj,1:num_nodes)=rnodes_temp(1:num_nj,1:num_nodes)
      deallocate(rnodes_temp)
      node_field(1:num_nj,num_nodes+1:num_nodes_new)=0.0_dp
    endif

    allocate(nodelem_temp(num_elems))
    nodelem_temp = elem_symmetry ! copy to temporary array
    deallocate(elem_symmetry) !deallocate initially allocated memory
    allocate(elem_symmetry(num_elems_new))
    elem_symmetry(1:num_elems)=nodelem_temp(1:num_elems)
    deallocate(nodelem_temp) !deallocate the temporary array
    elem_symmetry(num_elems+1:num_elems_new)=1

    allocate(enodes_temp2(-1:1,0:2,0:num_elems))
    enodes_temp2=elem_cnct
    deallocate(elem_cnct)
    allocate(elem_cnct(-1:1,0:2,0:num_elems_new))
    elem_cnct(-1:1,0:2,0:num_elems)=enodes_temp2(-1:1,0:2,0:num_elems)
    deallocate(enodes_temp2)
    elem_cnct(-1:1,0:2,num_elems+1:num_elems_new) = 0

    allocate(enodes_temp(num_ord,num_elems))
    enodes_temp=elem_ordrs
    deallocate(elem_ordrs)
    allocate(elem_ordrs(num_ord,num_elems_new))
    elem_ordrs(1:num_ord,1:num_elems)=enodes_temp(1:num_ord,1:num_elems)
    deallocate(enodes_temp)
    elem_ordrs(1:num_ord,num_elems+1:num_elems_new) = 0

    if(allocated(elem_units_below).and.num_nu.gt.0)then
      allocate(nodelem_temp(num_elems))
      nodelem_temp=elem_units_below
      deallocate(elem_units_below)
      allocate(elem_units_below(num_elems_new))
      elem_units_below(1:num_elems)=nodelem_temp(1:num_elems)
      deallocate(nodelem_temp)
      elem_units_below(num_elems+1:num_elems_new)=0
    endif

    allocate(enodes_temp(num_nodes,0:3))
    enodes_temp=elems_at_node
    deallocate(elems_at_node)
    allocate(elems_at_node(num_nodes_new,0:3))
    elems_at_node(1:num_nodes,0:3)=enodes_temp(1:num_nodes,0:3)
    deallocate(enodes_temp)
    elems_at_node(num_nodes+1:num_nodes_new,0:3)=0

    if(model_type.eq.'gas_mix')then
      allocate(exp_temp(num_elems))
      exp_temp = expansile
      deallocate(expansile)
      allocate(expansile(num_elems_new))
      expansile(1:num_elems)=exp_temp(1:num_elems)
      deallocate(exp_temp)
      expansile(num_elems+1:num_elems_new)=.false.
    endif

    call enter_exit(sub_name,2)

  end subroutine reallocate_node_elem_arrays

!!!###############################################################

  function get_local_node_f(ndimension,np_global) result(get_local_node)

!! dummy arguments
    integer,intent(in) :: ndimension,np_global
!!! local variables
    integer :: np
    integer :: get_local_node
    logical :: found

    found = .false.
    np = 0

    select case (ndimension)

    case(1)
       do while (.not.found)
          np=np+1
          if(np.gt.num_nodes) then
             found = .true.
             write(*,'('' Warning: local node not found for global node'',I6)') np_global
          endif
          if(np_global.eq.nodes(np))then
             get_local_node = np
             found = .true.
          endif
       enddo

    case(2)
       do while (.not.found)
          np=np+1
          if(np.gt.num_nodes_2d) then
             found = .true.
             write(*,'('' Warning: local node not found for global node'',I6)') np_global
             read(*,*)
          endif
           if(np_global.eq.nodes_2d(np))then
             get_local_node = np
             found = .true.
          endif
       enddo

    end select


  end function get_local_node_f


!###################################################################################
!
!*get_final_integer*

  function get_final_integer(string)
    
    character,intent(in) :: string*(*)
    integer :: ibeg,iend,ierror,nsign,ntemp
    character :: sub_string*(40)
    
    integer :: get_final_integer
    
    ! ###########################################################################
    
    iend=len(string) !get the length of the string
    ibeg=index(string,":")+1 !get location of integer in string, follows ":"
    sub_string = adjustl(string(ibeg:iend)) ! get the characters beyond ":"
    iend=len(sub_string) !length of the sub-string
    if(sub_string(1:1).eq.'-')then !check for negative sign
       nsign=-1
       ibeg=2
    else
       nsign=1
       ibeg=1
    endif
    
    read (sub_string(ibeg:iend), '(i10)', iostat=ierror ) ntemp !get integer values
    if(ierror.gt.0)then
       !... something wrong with data
       write(*,'(''Data read error'')')
       write(*,'(a)') sub_string(ibeg:iend)
    endif
    ntemp=ntemp*nsign !apply sign to number
    
    get_final_integer=ntemp !return the integer value
    
  end function get_final_integer
  
!!!##################################################

  subroutine get_four_nodes(ne,string)

    !     Parameter List
    integer, INTENT(IN) :: ne
    character(len=132), INTENT(IN) :: string
    integer :: ibeg,iend,i_ss_end,nn,np_global
    character(len=40) :: sub_string
    

    iend=len(string)
    ibeg=index(string,":")+1 !get location of first integer in string
    sub_string = adjustl(string(ibeg:iend)) ! get the characters beyond : and remove the leading blanks
    i_ss_end=len(sub_string) !get the end location of the sub-string
    
    ibeg=1
    do nn=1,4
       iend=index(sub_string," ") !get location of first blank in sub-string
       read (sub_string(ibeg:iend-1), '(i7)' ) np_global
       sub_string = adjustl(sub_string(iend:i_ss_end)) ! get chars beyond blank, remove leading blanks
       elem_nodes_2d(nn,ne)=get_local_node_f(2,np_global)
    enddo ! nn
    
  end subroutine get_four_nodes

! 
! ##########################################################################      
!
  subroutine mesh_angle(ANGLE,XP1,XP2,XP3)

!#### Subroutine: MESH_ANGLE
!###  Description:
!###    MESH_ANGLE calculates the angle between vectors XP1-XP2 and XP2-XP3
    use arrays,only: dp,elems,elem_cnct,elem_direction,elem_field
    use diagnostics,only: enter_exit

    IMPLICIT NONE
!   INCLUDE 'geom00.cmn'
!   INCLUDE 'tol00.cmn'

    real(dp),INTENT(IN) :: ANGLE
    real(dp),INTENT(IN) :: XP1(3),XP2(3),XP3(3)

    !Local variables
    INTEGER :: nj
    REAL(dp) :: SCALAR,U(3),V(3)
    character(len=60) :: sub_name

    sub_name = 'mesh_angle'
    call enter_exit(sub_name,1)


     DO nj=1,3
       U(nj)=XP2(nj)-XP1(nj) !end - start
       V(nj)=XP3(nj)-XP2(nj) !end - start
     ENDDO !nj

      CALL normalise(3,U)
      CALL normalise(3,V)


      ANGLE=scalar_f(3,U,V)
      ANGLE=MAX(-1.d0,ANGLE)
      ANGLE=MIN(1.d0,ANGLE)
      ANGLE=DACOS(ANGLE)

    call enter_exit(sub_name,2)
  end subroutine mesh_angle


! ##########################################################################
!
  subroutine mesh_plane_angle(NBJ,ne,NPNE,ANGLE,XP)

!#### Subroutine: MESH_PLANE_ANGLE
!###  Description:
!###  MESH_PLANE_ANGLE calculates the rotation angle between the
!###  branching planes of a parent (ne) and its child branches.

    use arrays,only: dp,elem_cnct

      IMPLICIT NONE

!!     INCLUDE 'b00.cmn'
!!     INCLUDE 'b01.cmn'
!!     INCLUDE 'cbdi02.cmn'
!!     INCLUDE 'geom00.cmn'
!!     INCLUDE 'tol00.cmn'

!!  Parameter list
      INTEGER :: NBJ(NJM,NEM),ne,NPNE(NNM,NBFM,NEM),&
        NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 ANGLE,XP(NKM,NVM,NJM,NPM)

!! Local variables
      INTEGER nb,ne0,ne1,nj,np1,np2,np3,np4,np5
      REAL*8 norm_1(4),norm_2(4),SCALAR,XPOINT(3,5),temp
      REAL*8 CALC_ANGLE

      nb=NBJ(1,ne)
      ne0=NXI(-1,1,ne) !parent element
      np1=NPNE(1,nb,ne) !start node
      np2=NPNE(2,nb,ne) !end node
      ne1=NXI(1,1,ne0)

!!     IF(DOP)THEN
!!       WRITE(OP_STRING,'(5(I6))') nb,ne0,np1,np2,ne1
!!       CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
!!     ENDIF

      IF(ne1.EQ.ne) ne1=NXI(1,2,ne0) !sibling
      np3=NPNE(2,nb,ne1) !end node of sibling
      np4=NPNE(2,nb,NXI(1,1,ne)) !end node of first child
      np5=NPNE(2,nb,NXI(1,2,ne)) !end node of second child

      DO nj=1,3
        XPOINT(nj,1)=XP(1,1,nj,np1)
        XPOINT(nj,2)=XP(1,1,nj,np2)
        XPOINT(nj,3)=XP(1,1,nj,np3)
        XPOINT(nj,4)=XP(1,1,nj,np4)
        XPOINT(nj,5)=XP(1,1,nj,np5)
      ENDDO !nj

      CALL plane_from_3_points(norm_1,2,XPOINT(1,1),XPOINT(1,2),XPOINT(1,3))
      CALL normalise2(3,norm_1,temp) !unit vector

      CALL plane_from_3_points(norm_2,2,XPOINT(1,2),XPOINT(1,4),XPOINT(1,5))
      CALL normalise2(3,norm_2,temp)

      ANGLE=CALC_ANGLE(norm_1,norm_2)

  end subroutine mesh_plane_angle

!



! ##########################################################################      
! 

  subroutine normalise(NUMCMPTS,VECTOR)

!#### Subroutine: NORMALISE
!###  Description:
!###  NORMALISE divides the components of VECTOR by it's length

    use arrays,only: dp
    use diagnostics,only: enter_exit
    character(len=60) :: sub_name

    implicit none
!!  INCLUDE 'cbdi02.cmn'
!!  INCLUDE 'tol00.cmn'

!!  Parameter List
    integer :: NUMCMPTS
    real(dp) :: VECTOR(*)

!!  Local Variables
    INTEGER :: ni
    REAL(dp) :: VECTOR_LENGTH

    sub_name = 'normalise'

    VECTOR_LENGTH=0.0d0
     DO ni=1,NUMCMPTS
       VECTOR_LENGTH=VECTOR_LENGTH+VECTOR(ni)*VECTOR(ni)
     ENDDO !ni
     VECTOR_LENGTH=DSQRT(VECTOR_LENGTH)

     IF(VECTOR_LENGTH.GT.ZERO_TOL) THEN
       DO ni=1,NUMCMPTS
         VECTOR(ni)=VECTOR(ni)/VECTOR_LENGTH
       ENDDO !ni
     ELSE
       WRITE(*,*) ' >>WARNING: Cannot normalise a zero length vector'
     ENDIF

    call enter_exit(sub_name,2)
  end subroutine normalise

!
! ##########################################################################      
!
  subroutine linear_regression(N,r_squared,slope,X,Y)

!#### Subroutine: LINREGRESS
!###  Description:
!###    Calculates linear regression equation and r-squared 
!###    correlation coefficient for a set of data.

!*** Created by Kelly Burrowes, March 2003.
    use arrays,only: dp
    use diagnostics,only: enter_exit

    implicit none
!!  INCLUDE 'b01.cmn'
!!  INCLUDE 'cbdi02.cmn'
!!  Parameter list
    integer, INTENT(IN) :: N
    real(dp), INTENT(IN) :: X(N),Y(N)
    real(dp), INTENT(OUT) :: r_squared,slope

!!  Local variables
    integer :: i
    REAL(dp) :: AX,AY,intercept,R,SXX,SXY,SYY
    REAL(dp) :: XSUM,XT,XXSUM,XYSUM,YSUM,YT

    sub_name = 'linear_regression'
    call enter_exit(sub_name,1)

      YSUM=0.d0
      XSUM=0.d0
      XXSUM=0.d0
      XYSUM=0.d0
      DO i=1,N
        YSUM=YSUM+Y(i)
        XSUM=XSUM+X(i)
        XYSUM=XYSUM+X(i)*Y(i)
        XXSUM=XXSUM+X(i)*X(i)
      ENDDO !N

!... calculate least squares estimate of straight line thru solution
      slope=(XYSUM-XSUM*YSUM/N)/(XXSUM-XSUM*XSUM/N)
      intercept=(YSUM/N)-(slope*XSUM/N)

!... calculate r-squared correlation coefficient
!... see Numerical Recipes, Fortran 77, 2nd edition, page 632.      
      AX=XSUM/N !mean of X
      AY=YSUM/N !mean of Y
      SXX=0.d0
      SYY=0.d0
      SXY=0.d0
      DO i=1,N
        XT=X(i)-AX
        YT=Y(i)-AY
        SXX=SXX+XT**2.d0
        SYY=SYY+YT**2.d0
        SXY=SXY+XT*YT
      ENDDO
      R=SXY/DSQRT(SXX*SYY)
      r_squared=R**2
!     IF(DOP)THEN
!       WRITE(OP_STRING,'('' Gradient: '',F8.3,'' Intercept: '','
!    '    //' F8.3,''R-squared correlation'',F8.3)')
!    '    slope,intercept,r_squared
!       CALL WRITES(IODI,OP_STRING,ERROR,*9999)
!     ENDIF

    call enter_exit(sub_name,2)
  end subroutine linear_regression


!
! ##########################################################################      
!

  subroutine normalise2(NUMCMPTS,VECTOR,VECTOR_LENGTH)

!#### Subroutine: NORMALISE2
!###  Description:
!###    NORMALISE divides the components of VECTOR by its length
!###    and returns the original length

      IMPLICIT NONE
!!     INCLUDE 'cbdi02.cmn'
!!     INCLUDE 'tol00.cmn'

!!!     Parameter List
      INTEGER NUMCMPTS
      REAL*8 VECTOR(*)

!!     Local Variables
      INTEGER ni
      REAL*8 VECTOR_LENGTH
    sub_name = 'normalise2'
    call enter_exit(sub_name,1)

      VECTOR_LENGTH=0.0d0
      DO ni=1,NUMCMPTS
        VECTOR_LENGTH=VECTOR_LENGTH+VECTOR(ni)*VECTOR(ni)
      ENDDO !ni
      VECTOR_LENGTH=DSQRT(VECTOR_LENGTH)

      IF(VECTOR_LENGTH.GT.ZERO_TOL) THEN
        DO ni=1,NUMCMPTS
          VECTOR(ni)=VECTOR(ni)/VECTOR_LENGTH
        ENDDO !ni
      ELSE
        WRITE(*,*) 'WARNING: Cannot normalise a zero length vector'
      ENDIF

  end subroutine normalise2


! ##########################################################################      
! 
  subroutine plane_from_3_points(NORML,NORMALTYPE,POINT1,POINT2,POINT3)

!#### Subroutine: PLANE_FROM_3_PTS
!###  Description:
!###    PLANE_FROM_3_PTS finds the equation of a plane in three
!###    dimensions and a vector normal to the plane from three
!###    non-collinear points.
!###    NORMALTYPE=1 for raw normal and plane equation
!###    NORMALTYPE=2 for unit normal and plane equation
!###    The coefficients represent aX + bY + cZ + d = 0
!###    NORML(1)=a,NORML(2)=b,NORML(3)=c,NORML(4)=d
!***  Created by Martin Buist, Jan 1997


!!    INCLUDE 'cbdi02.cmn'
!!    INCLUDE 'geom00.cmn'
!!    INCLUDE 'tol00.cmn'

    use arrays,only: dp
    use diagnostics,only: enter_exit

    IMPLICIT NONE

!!  Parameter list
    INTEGER :: NORMALTYPE
    REAL*8 POINT1(3),POINT2(3),POINT3(3),NORML(4)
    CHARACTER ERROR*(*)

!!  Local variables
    INTEGER :: nj
    REAL(dp) :: DIFF1(3),DIFF2(3),NORMSIZE
    LOGICAL COLLINEAR
    character(len=60) :: sub_name

    sub_name = 'plane_from_3_points'
    call enter_exit(sub_name,1)


! Check for collinearity
      COLLINEAR=.FALSE.
      call check_collinear(POINT1,POINT2,POINT3,collinear)
      IF(.NOT.collinear) THEN
        DO nj=1,NJT
          DIFF1(nj)=POINT2(nj)-POINT1(nj)
          DIFF2(nj)=POINT2(nj)-POINT3(nj)
        ENDDO !nj

        NORML(1)=(DIFF1(2)*DIFF2(3))-(DIFF1(3)*DIFF2(2))
        NORML(2)=(DIFF1(3)*DIFF2(1))-(DIFF1(1)*DIFF2(3))
        NORML(3)=(DIFF1(1)*DIFF2(2))-(DIFF1(2)*DIFF2(1))

        IF(NORMALTYPE.EQ.2) THEN
          NORMSIZE=0.0d0
          DO nj=1,NJT
            NORMSIZE=NORMSIZE+(NORML(nj)**2.0d0)
          ENDDO !nj
          NORMSIZE=DSQRT(NORMSIZE)
          DO nj=1,NJT
            NORML(nj)=NORML(nj)/NORMSIZE
          ENDDO !nj
        ENDIF

        NORML(4)=0.0d0
        DO nj=1,NJT
          NORML(4)=NORML(4)-(NORML(nj)*POINT1(nj))
        ENDDO !nj

      ELSE !Colinear

        DO nj=1,4
          NORML(nj)=0.0d0
        ENDDO !nj

      ENDIF

    call enter_exit(sub_name,2)
  end subroutine plane_from_3_points
! 
! ##########################################################################      


  function inlist(item,ilist)
!!! dummy arguments
    integer :: item,ilist(:)
! local variables
    integer :: n
    logical :: inlist

    inlist = .false.
    do n=1,size(ilist)
       if(item == ilist(n)) inlist = .true.
    enddo

  end function inlist
! 
!###########################################################################################
! 
  function coord_at_xi(ne,xi,basis)

    ! Parameters
    integer,intent(in) :: ne
    real(dp),intent(in) :: xi(:)
    character(len=*),intent(in) :: basis

    ! Local Variables
    integer :: nn,nv
    real(dp) :: phi(4),phi_10(2),phi_11(2),phi_20(2),phi_21(2),x(4,3,4)
    real(dp) :: coord_at_xi(3)

    select case (basis)
      case('linear')
         forall (nn=1:4) x(1,1:3,nn) = node_xyz_2d(1,1,1:3,elem_nodes_2d(nn,ne))
         phi(1) = (1.0_dp - xi(1))*(1.0_dp - xi(2))
         phi(2) = xi(1)*(1.0_dp - xi(2))
         phi(3) = (1.0_dp - xi(1))*xi(2)
         phi(4) = xi(1)*xi(2)

         coord_at_xi(1:3) = phi(1)*x(1,1:3,1)+phi(2)*x(1,1:3,2)+phi(3)*x(1,1:3,3)+phi(4)*x(1,1:3,4)

      case('hermite')
         do nn=1,4
           nv = elem_versn_2d(nn,ne)
           x(1:4,1:3,nn) = node_xyz_2d(1:4,nv,1:3,elem_nodes_2d(nn,ne))
         enddo
         phi_10(1) = (2.0_dp*xi(1)-3.0_dp)*xi(1)*xi(1)+1.0_dp  ! 2xi^3-3xi^2+1
         phi_10(2) = (2.0_dp*xi(2)-3.0_dp)*xi(2)*xi(2)+1.0_dp  ! 2xi^3-3xi^2+1
         phi_20(1) = xi(1)*xi(1)*(3.0_dp-2.0_dp*xi(1))         ! -2xi^3+3xi^2
         phi_20(2) = xi(2)*xi(2)*(3.0_dp-2.0_dp*xi(2))         ! -2xi^3+3xi^2
         phi_11(1) = ((xi(1)-2.0_dp)*xi(1)+1.0_dp)*xi(1)       ! xi^3-2xi^2+xi
         phi_11(2) = ((xi(2)-2.0_dp)*xi(2)+1.0_dp)*xi(2)       ! xi^3-2xi^2+xi
         phi_21(1) = xi(1)*xi(1)*(xi(1)-1.0_dp)                ! xi^3-xi^2
         phi_21(2) = xi(2)*xi(2)*(xi(2)-1.0_dp)                ! xi^3-xi^2
         coord_at_xi(1:3) = phi_10(1)*phi_10(2)*x(1,1:3,1) &
                            + phi_20(1)*phi_10(2)*x(1,1:3,2) &
                            + phi_10(1)*phi_20(2)*x(1,1:3,3) &
                            + phi_20(1)*phi_20(2)*x(1,1:3,4) &
                            + phi_11(1)*phi_10(2)*x(2,1:3,1) &
                            + phi_21(1)*phi_10(2)*x(2,1:3,2) &
                            + phi_11(1)*phi_20(2)*x(2,1:3,3) &
                            + phi_21(1)*phi_20(2)*x(2,1:3,4) &
                            + phi_10(1)*phi_11(2)*x(3,1:3,1) &
                            + phi_20(1)*phi_11(2)*x(3,1:3,2) &
                            + phi_10(1)*phi_21(2)*x(3,1:3,3) &
                            + phi_20(1)*phi_21(2)*x(3,1:3,4) &
                            + phi_11(1)*phi_11(2)*x(4,1:3,1) &
                            + phi_21(1)*phi_11(2)*x(4,1:3,2) &
                            + phi_11(1)*phi_21(2)*x(4,1:3,3) &
                            + phi_21(1)*phi_21(2)*x(4,1:3,4)
     end select

  end function coord_at_xi

! 
!###########################################################################################
! 
end module geometry
