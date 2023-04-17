module aux_funcs
  !=============================================================================
  !  Global operators for shallow water model
  !
  ! Pedro da Silva Peixoto (pedrosp@ime.usp.br)
  ! Oct 2018
  !=============================================================================

  !Use main grid data structures
  use constants
  use datastruct, only: &
    grid_structure, &
    scalar_field, &
    vector_variable_uv, &
    vector_field_cart

  use smeshpack
  use legendre_gauss

  contains
    subroutine test(mesh)
      implicit none
      type(grid_structure), intent(in)        :: mesh
      real*8  :: dy,y(150),wts(150)
      integer :: i
      
      call p_quadrature_rule ( 150, y, wts )
      
      do i=1,150
        write(*,*) y(i),wts(i)
      enddo
        write(*,*) sum(wts)
      stop

    end subroutine test

    function velTC9(lat)
      real*8, intent(in) :: lat
      real*8,parameter      :: u0=80d0, &
                               lat0 = pi/7, &
                               lat1 = pi/2-lat0, &
                               en=dexp(-4/(lat1-lat0)**2)

      if (lat <= lat0) then
        velTC9 = 0d0;
      elseif (lat < lat1) then
        velTC9 = u0/en*dexp(1d0/((lat-lat0)*(lat-lat1)));
      else
        velTC9 = 0d0;
      endif

      return

    end function velTC9

    function hTC9(lat)
      implicit none
      real*8             :: hTC9
      real*8,intent(in)  :: lat
      real*8             :: dy,y,u,f,h,ni
      integer :: i,j,sgn

      ni = 10000;

      h = 10000;

      dy = abs(lat)/ni;
      sgn = 1;
      if (lat < 0) then
        sgn = -1;
      endif
      y = 0d0;

      do i =1,ni
        y = y + sgn*dy;
        u = velTC9(y);
        f = 2d0*omega*dsin(y)
        h = h - rad*u*(f + dtan(y)*u/rad)*dy/g;
      enddo

      hTC9 = h;
      return

    end function hTC9


    subroutine integrate_vector(h,ni,mesh)
      implicit none
      type(grid_structure), intent(in)        :: mesh
      integer, intent(in)                     :: ni
      real*8, intent(inout)                   :: h(ni)
      real*8  :: dy,lat(ni),y,u,f
      integer :: i,j,sgn

      h = 10000;

      dy = pi/(ni-1);
      do i =1,ni
        lat(i) = -pi/2+(i-1)*dy
      enddo

      do i = 1,ni
        dy = abs(lat(i))/ni;
        sgn = 1;
        if (lat(i) < 0) then
          sgn = -1;
        endif
        y = 0d0;
        do j =1,ni
          y = y + sgn*dy;
          u = velTC9(y);
          f = 2d0*omega*dsin(y)
          h(i) = h(i) - rad*u*(f + dtan(y)*u/rad)*dy/g;
        enddo
      enddo
!      h = h + 269.28178505923097d0

!      write(*,*) sum(h)/(ni)
!      stop

    end subroutine integrate_vector

!    subroutine integrate_vector(func,lat,wts,funcint,mesh)
!      implicit none
!      type(grid_structure), intent(in)        :: mesh
!      real*8,intent(in)                       :: func(:),lat(:),wts(:)
!      real*8,intent(inout)                    :: funcint(:)
!      real*8  :: dy
!      integer :: i,idx

!      dy = lat(2) - lat(1);

!      do i =2,size(func,1)
!!        funcint(i) = dy*func(i)*wts(i)
!        funcint(i) = sum(func(1:i)*wts(1:i))
!      enddo

!!      do i =2,size(func,1)
!!        funcint(i) = sum(funcint(1:i));
!!      enddo


!    end subroutine integrate_vector


end module aux_funcs
