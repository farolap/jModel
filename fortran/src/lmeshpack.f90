module lmeshpack

  use constants
  use datastruct

  implicit none

  contains

    subroutine loadpts(mesh, fileg)
      implicit none
      integer :: nq,ne,nt,i
      type(grid_structure),intent(inout) :: mesh
      character(len=100)                 :: filen
      character(len=100),intent(in)      :: fileg

      write(filen,'(2A)') trim(fileg),'/vertices.dat'             ! Vertices
      open(unit = 100, file = filen);

!      write(filen,'(2A)') trim(fileg),'/edges.dat';               ! Edges
!      open(unit = 101, file = filen);

      write(filen,'(2A)') trim(fileg),'/triangles.dat';           ! Vertices forming triangle
      open(unit = 102, file = filen);

      write(filen,'(2A)') trim(fileg),'/etri.dat';                ! Edges forming triangle
      open(unit = 103, file = filen);

      write(filen,'(2A)') trim(fileg),'/edgetrishare.dat';        ! Triangles sharing edge
      open(unit = 104, file = filen);

      write(filen,'(2A)') trim(fileg),'/edgevershare.dat';        ! Vertices sharing edge
      open(unit = 105, file = filen);

      write(filen,'(2A)') trim(fileg),'/vertexedgshare.dat';      ! Edges neighbour to vertex
      open(unit = 106, file = filen);

      write(filen,'(2A)') trim(fileg),'/vertextrishare.dat';      ! Triangle neighbour to vertex
      open(unit = 107, file = filen);

      write(filen,'(2A)') trim(fileg),'/vtrcount.dat';            ! Number of triangles sharing 
      open(unit = 108, file = filen);

      write(filen,'(2A)') trim(fileg),'/vedcount.dat';
      open(unit = 109, file = filen);

      read(100,*) mesh%nv,nq; ! vtx pts
!      read(101,*) mesh%ne,ne; ! triangles
      read(102,*) mesh%nt,nt; ! edge pts

      mesh%ne = mesh%nt + mesh%nv-2 

      do i=103,109
        read(i,*);
      enddo

      allocate(mesh%tr(mesh%nt), mesh%v(mesh%nv),mesh%ed(mesh%ne), &
          mesh%edhx(mesh%ne),mesh%hx(mesh%nv));

      do i = 1,mesh%nv
        read(100,*) mesh%v(i)%p           ! Vertices
        read(109,*) mesh%v(i)%nnb         !
        read(108,*) mesh%v(i)%nnbt        !

        allocate(mesh%v(i)%ed(mesh%v(i)%nnb),mesh%v(i)%tr(mesh%v(i)%nnbt), &
                mesh%hx(i)%trhx_areag(mesh%v(i)%nnbt),mesh%hx(i)%ctrhx_areag(mesh%v(i)%nnbt));
        allocate(mesh%hx(i)%cor(mesh%v(i)%nnb),mesh%hx(i)%cord(mesh%v(i)%nnb,2));
        allocate(mesh%hx(i)%trskwg(mesh%v(i)%nnb,mesh%v(i)%nnb));

        read(106,*) mesh%v(i)%ed           !
        read(107,*) mesh%v(i)%tr           !
      enddo
      do i = 1,mesh%ne
!        read(101,*) mesh%ed(i)%c%p         ! Edges
        read(104,*) mesh%ed(i)%sh          ! Triangles sharing edge
        read(105,*) mesh%ed(i)%v           ! vertices  sharing edge

      enddo
      
      do i = 1,mesh%nt
        read(102,*) mesh%tr(i)%v           ! Vertices of triangles
        read(103,*) mesh%tr(i)%ed          ! edges of triangles
      enddo


      do i =100,125
        close(i)
      enddo

    end subroutine loadpts

    subroutine loadall(fileg,mesh)
      implicit none
      integer :: nq,ne,nt,i
      type(grid_structure),intent(inout) :: mesh
      character(len=100)                 :: filen
      character(len=100),intent(in)      :: fileg

      call loadpts(mesh,fileg)

      write(filen,'(2A)') trim(fileg),'/edges.dat'
      open(unit = 100, file = filen);

      write(filen,'(2A)') trim(fileg),'/edges_hx.dat'
      open(unit = 101, file = filen);

      write(filen,'(2A)') trim(fileg),'/edge_nr.dat'
      open(unit = 102, file = filen);

      write(filen,'(2A)') trim(fileg),'/edge_tg.dat'
      open(unit = 103, file = filen);

      write(filen,'(2A)') trim(fileg),'/edge_hx_nr.dat'
      open(unit = 104, file = filen);

      write(filen,'(2A)') trim(fileg),'/edge_hx_tg.dat'
      open(unit = 105, file = filen);

      write(filen,'(2A)') trim(fileg),'/edge_leng.dat'
      open(unit = 106, file = filen);

      write(filen,'(2A)') trim(fileg),'/edge_hx_leng.dat'
      open(unit = 107, file = filen);

      write(filen,'(2A)') trim(fileg),'/tri_area.dat'
      open(unit = 108, file = filen);

      write(filen,'(2A)') trim(fileg),'/tri_edge_area.dat'
      open(unit = 109, file = filen);

      write(filen,'(2A)') trim(fileg),'/circumcenter.dat'
      open(unit = 110, file = filen);

      write(filen,'(2A)') trim(fileg),'/centroid.dat'
      open(unit = 111, file = filen);

      write(filen,'(2A)') trim(fileg),'/hx_area.dat'
      open(unit = 112, file = filen);

      do i =1,mesh%ne
        read(100,*) mesh%ed(i)%c%p
        read(101,*) mesh%edhx(i)%c%p
        read(102,*) mesh%ed(i)%nr
        read(103,*) mesh%ed(i)%tg
        read(104,*) mesh%edhx(i)%nr
        read(105,*) mesh%edhx(i)%tg
        read(106,*) mesh%ed(i)%leng
        read(107,*) mesh%edhx(i)%leng
      enddo

      do i=1,mesh%nt
        read(108,*) mesh%tr(i)%areag
        read(109,*) mesh%tr(i)%ctrve_areag
        read(110,*) mesh%tr(i)%c%p
        read(111,*) mesh%tr(i)%b%p
      enddo

      do i=1,mesh%nv
        read(112,*) mesh%hx(i)%areag
      enddo

      do i=100,112
          close(i)
      enddo

    endsubroutine

end module lmeshpack
