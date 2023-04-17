module swma_operators
  use constants
  use datastruct, only: &
    grid_structure, &
    scalar_field, &
    vector_variable_uv, &
    vector_field_cart

  use basic_funcs
  
  implicit none
  contains

    subroutine scalar_ve2tr(fve,ftr,mesh)
      implicit none

      type(grid_structure), intent(in)  :: mesh
      type(scalar_field), intent(in)    :: fve
      type(scalar_field), intent(inout) :: ftr
      integer :: i

      !$omp parallel do &
      !$omp default(none) &
      !$omp shared(mesh, ftr, fve) &
      !$omp private(i) &
      !$omp schedule(static)
      do i =1,mesh%nt
          ftr%f(i) = (fve%f(mesh%tr(i)%v(1))*mesh%tr(i)%ctrve_areag(2) + &
              fve%f(mesh%tr(i)%v(2))*mesh%tr(i)%ctrve_areag(3) + &
              fve%f(mesh%tr(i)%v(3))*mesh%tr(i)%ctrve_areag(1))/mesh%tr(i)%areag;
      enddo
      !$omp end parallel do
    end subroutine scalar_ve2tr

    subroutine scalar_tr2ed(ftr,fed,mesh)
      implicit none

      type(grid_structure), intent(in)  :: mesh
      type(scalar_field), intent(in)    :: ftr
      type(scalar_field), intent(inout) :: fed
      integer :: i

      !$omp parallel do &
      !$omp default(none) &
      !$omp shared(mesh, ftr, fed) &
      !$omp private(i) &
      !$omp schedule(static)
      do i=1,mesh%ne
          fed%f(i) = (ftr%f(mesh%ed(i)%sh(1)) + ftr%f(mesh%ed(i)%sh(2)))*.5d0;
      enddo
      !$omp end parallel do

    end subroutine scalar_tr2ed

    subroutine vector_ve2tr(fve,ftr,mesh)
      implicit none

      type(grid_structure), intent(in)        :: mesh
      type(vector_field_cart), intent(in)    :: fve
      type(vector_field_cart), intent(inout) :: ftr

      real*8  :: vec(3)
      integer :: i,j
      integer :: idx(3)

      idx = (/2,3,1/);

      !$omp parallel do &
      !$omp default(none) &
      !$omp shared(mesh, fve, ftr,idx) &
      !$omp private(i,j,vec) &
      !$omp schedule(static)
      do i=1,mesh%nt
        vec = 0d0;
        do j=1,3
          vec = vec + fve%p(mesh%tr(i)%v(j))%v &
              *mesh%tr(i)%ctrve_areag(idx(j))/mesh%tr(i)%areag;
        enddo
        ftr%p(i)%v = vec;
      enddo
      !$omp end parallel do
      return

    end subroutine vector_ve2tr

    subroutine vector_tr2ed(ftr,fed,mesh)
      implicit none
      type(grid_structure), intent(in)        :: mesh
      type(vector_field_cart), intent(in)    :: ftr
      type(vector_field_cart), intent(inout) :: fed
      real*8  :: vec(3)
      integer :: i,j

      !$omp parallel do &
      !$omp default(none) &
      !$omp shared(mesh, ftr, fed) &
      !$omp private(i,j,vec) &
      !$omp schedule(static)
      do i=1,mesh%ne
        vec = 0d0;
        do j=1,2
          vec = vec + ftr%p(mesh%ed(i)%sh(j))%v/2d0;
        enddo
        fed%p(i)%v = proj_vec_sphere(vec, mesh%edhx(i)%c%p);
      enddo
      !$omp end parallel do
    end subroutine vector_tr2ed

    subroutine div_tr2ve(fu,div,mesh)
      implicit none
      type(grid_structure), intent(in)      :: mesh
      type(vector_field_cart), intent(in)  :: fu
      type(scalar_field), intent(inout)     :: div
      integer :: i,j,k,ide,sh
      real*8  :: uint(3),A,nec(3),vec

      div%f =0d0

      !$omp parallel do &
      !$omp default(none) &
      !$omp shared(mesh, fu, div) &
      !$omp private(i,j,ide,sh,uint,A,nec,vec) &
      !$omp schedule(static)
      do i =1,mesh%nv
        A = mesh%hx(i)%areag;
        do j =1,mesh%v(i)%nnb
          ide  = mesh%v(i)%ed(j);
          nec = mesh%edhx(ide)%nr;
          uint = 0d0;
          do k = 1,2
            sh = mesh%ed(ide)%sh(k);
            uint = uint + fu%p(sh)%v*.5d0;
          enddo
          
          vec = dot_product(nec,uint);
          ! vec = dot_product(nec,proj_vec_sphere(uint, mesh%edhx(ide)%c%p));
          div%f(i) = div%f(i) + &
                   vec*mesh%edhx(ide)%leng*mesh%hx(i)%cor(j);

        enddo
        div%f(i) = div%f(i)/(A*rad);
      enddo
      !$omp end parallel do
    end subroutine div_tr2ve

    subroutine vort_tr2ve(fu,zeta,mesh)
      implicit none
      type(grid_structure), intent(in)      :: mesh
      type(vector_field_cart), intent(in)  :: fu
      type(scalar_field), intent(inout)     :: zeta
      integer :: i,j,k,ide,sh
      real*8  :: uint(3),A,tgu(3),vec

      zeta%f =0d0

      !$omp parallel do &
      !$omp default(none) &
      !$omp shared(mesh, fu, zeta) &
      !$omp private(i,j,ide,sh,uint,A,tgu,vec) &
      !$omp schedule(static)
      do i =1,mesh%nv
        A = mesh%hx(i)%areag;
        do j =1,mesh%v(i)%nnb
          ide  = mesh%v(i)%ed(j);
          tgu = mesh%edhx(ide)%tg;

          uint = 0d0;
          do k = 1,2
            sh = mesh%ed(ide)%sh(k);
            uint = uint + fu%p(sh)%v*.5d0;
          enddo
          vec = dot_product(tgu,uint);
          zeta%f(i) = zeta%f(i) + vec*mesh%edhx(ide)%leng*mesh%hx(i)%cor(j);

        enddo
        zeta%f(i) = -zeta%f(i)/(A*rad);
      enddo
      !$omp end parallel do
    end subroutine vort_tr2ve

    subroutine ke_energy(u,ke)
      implicit none
      type(vector_field_cart), intent(in)    :: u
      type(scalar_field), intent(inout)       :: ke
      integer :: i,j
      
      ke%f = 0d0;
      do i =1,3
          !$omp parallel do &
          !$omp default(none) &
          !$omp shared(ke, u,i) &
          !$omp private(j) &
          !$omp schedule(static)
          do j=1,size(ke%f,1)
              ke%f(j) = ke%f(j) + (u%p(j)%v(i)**2)*.5d0;
          enddo
          !$omp end parallel do
      enddo
    
    end subroutine ke_energy

    subroutine grad_tr2ve(he,h,grad,mesh)
      implicit none
      type(grid_structure), intent(in)        :: mesh
      type(scalar_field), intent(in)          :: he,h
      type(vector_field_cart), intent(inout) :: grad
      integer :: i,j,sh(2),ide,cor
      real*8  :: A,nec(3),vec(3),fint,de

      do i= 1,3
        grad%p(:)%v(i) =0d0
      enddo

      !$omp parallel do &
      !$omp default(none) &
      !$omp shared(mesh, he,h,grad) &
      !$omp private(A,i,j,sh,nec,vec,fint,ide,de,cor) &
      !$omp schedule(static)
      do i =1,mesh%nv
        A   = mesh%hx(i)%areag*rad;
        vec = 0d0
        do j =1,mesh%v(i)%nnb
          ide  = mesh%v(i)%ed(j);
          sh = mesh%ed(ide)%sh;
          de = mesh%edhx(ide)%leng;
          fint = he%f(ide)-h%f(i);

          nec = mesh%edhx(ide)%nr;

          cor    = mesh%hx(i)%cor(j)

          vec = vec + fint*nec*cor*de;

        enddo

        grad%p(i)%v = proj_vec_sphere(vec,mesh%v(i)%p)/A;
      enddo
      !$omp end parallel do


    end subroutine grad_tr2ve

    subroutine visc2x(u,lap,mesh)
      implicit none
      type(grid_structure), intent(in)        :: mesh
      type(vector_field_cart), intent(in)    :: u
      type(vector_field_cart), intent(inout) :: lap
      real*8  :: A,pq(3),peu(3),vecr(3),vec(3), &
                 veci(3)
      integer :: i,j,vedi,idv(2),cor

      do i= 1,3
        lap%p(:)%v(i) =0d0
      enddo

      !$omp parallel do &
      !$omp default(none) &
      !$omp shared(mesh, u,lap) &
      !$omp private(i,j,A,pq) &
      !$omp private(peu,vedi,idv,vecr,vec,veci) &
      !$omp private(cor) &
      !$omp schedule(static)
      do i =1,mesh%nv
        A = mesh%hx(i)%areag
        pq = mesh%v(i)%p;
        vec = 0d0;
        do j =1,mesh%v(i)%nnb
          vedi = mesh%v(i)%ed(j);
          peu  = mesh%ed(vedi)%c%p
          idv  = mesh%ed(vedi)%v(:);

          vecr = peu - pq;
          veci = mesh%v(idv(2))%p - mesh%v(idv(1))%p;
          veci = veci/norm(veci);
          cor  = int(dsign(1d0,dot_product(veci,vecr)));
          vec  = vec + (u%p(idv(2))%v - u%p(idv(1))%v)*cor*mesh%edhx(vedi)%leng/mesh%ed(vedi)%leng;

        enddo
        lap%p(i)%v = proj_vec_sphere(vec,mesh%v(i)%p)/(A*rad**2);
      enddo
      !$omp end parallel do

    end subroutine visc2x

    subroutine grad_vector_ve2tr(u,uperp,gradu,gradv,mesh)
      implicit none
      type(grid_structure), intent(in)        :: mesh
      type(vector_field_cart), intent(in)    :: u,uperp
      type(vector_field_cart), intent(inout) :: gradu,gradv
      integer :: i,j,edv(2),ide
      real*8  :: us(3),up(3),uiu(3),uiv(3),de,nec(3)
      
      do i =1,3
        gradu%p(:)%v(i) = 0d0;
        gradv%p(:)%v(i) = gradu%p(:)%v(i);
      enddo
      
      !$omp parallel do &
      !$omp default(none) &
      !$omp shared(mesh, u,uperp,gradu,gradv) &
      !$omp private(i,j,de,edv,ide) &
      !$omp private(us,up) &
      !$omp private(uiu,uiv,nec) &
      !$omp schedule(static)
      do i = 1,mesh%nt
        uiu = 0d0;
        uiv = 0d0;
        do j = 1,3
          ide = mesh%tr(i)%ed(j);
          nec = mesh%ed(ide)%nr;
          edv = mesh%ed(ide)%v;

          de = mesh%ed(ide)%leng;

          us = (u%p(edv(1))%v + u%p(edv(2))%v)*.5d0*mesh%tr(i)%cor(j);
          up = (uperp%p(edv(1))%v + uperp%p(edv(2))%v)*.5d0*mesh%tr(i)%cor(j);

          uiu = uiu + us*de;
          uiv = uiv + up*de;


        enddo
        gradu%p(i)%v(:) = uiu/(mesh%tr(i)%areag*rad);
        gradv%p(i)%v(:) = uiv/(mesh%tr(i)%areag*rad);

        gradu%p(i)%v(:) = proj_vec_sphere(gradu%p(i)%v,mesh%tr(i)%c%p);
        gradv%p(i)%v(:) = proj_vec_sphere(gradv%p(i)%v,mesh%tr(i)%c%p);
      enddo
      !$omp end parallel do
    end subroutine grad_vector_ve2tr

    subroutine loadini(u,h,j,filei,mesh)
        type(grid_structure),intent(in) :: mesh
        type(scalar_field),intent(inout) :: h
        type(vector_field_cart),intent(inout) :: u
        integer,intent(in) :: j
        integer :: i
        character(len=100),intent(in) :: filei
        character(len=100)            :: filen

        print*, 'Loading',j

        call charsave(filei,filen,'h',j,4)
        open(unit = 110, file = filen);

        call charsave(filei,filen,'u',j,4)
        open(unit = 111, file = filen);

        do i=1,mesh%nv
            read(110,*) h%f(i);
            read(111,*) u%p(i)%v;
        enddo
        do i=110,111
            close(i);
        enddo
    endsubroutine loadini

    subroutine inicond(j,u,utr,h,htr,hed,bath,div,gradb,zeta,f,mesh,resn)
        implicit none
        type(grid_structure), intent(in) :: mesh
        type(scalar_field), intent(inout):: h,htr,div,zeta,f,bath,hed
        type(vector_field_cart), intent(inout) :: u,utr,gradb
        integer,intent(in) :: j
        integer            :: i,m,n
        real*8             :: u0,b0,cte,dfdx,dfdy,ulon, &
            dudx,dudy,dvdx,dvdy,ulat, &
            h0,R0,lonc,latc,r,lon,lat,alpha, &
            p(3),ufv(3),hf,uf,ff,lonvec(3), &
            nygg,dygg,l1,beta, &
            cc1,cc2,cc3,cc4, &
            den,lat0,lat1,umen, &
            utmp,vectmp(3), &
            u00,en,lat2,l2,clat,e1,e2
        integer :: jy
        real*8,allocatable :: hgg(:)
        real*8             :: vec(3)
        character(len=100) ::resn

        bath%f       = 0d0;
        zeta%fexact  = 0d0;
        div%fexact   = 0d0;
        h%fexact     = 0d0;
        f%f          = 0d0;
        u%p(:)%v(1)  = 0d0;
        u%p(:)%v(2)  = 0d0;
        u%p(:)%v(3)  = 0d0;

        do i = 1,3
            gradb%pexact(:)%v(i) = 0d0;
            gradb%pexact(:)%v(i) = 0d0;
            u%pexact(:)%v(i)     = 0d0;
            u%pexact(:)%v(i)     = 0d0;
        enddo

        if (j==2) then
            write(resn,*) '../resultA/TC2/'

            u0 = 2d1;
            h0 = 5960d0;
            b0 = 2000d0;
            R0 = pi/9d0;
            lonc = -pi*.5d0;
            latc = pi/6d0;
            cte = 1d0/g*(rad*omega*u0+u0**2*.5d0);

            do i=1,mesh%nv
                lon = mesh%v(i)%lon;
                lat = mesh%v(i)%lat;

                r   = dsqrt(min(R0**2,(lon-lonc)**2+(lat-latc)**2));
                bath%f(i) = b0*(1d0-r/R0);
                h%f(i)    = h0-cte*dsin(lat)**2 - bath%f(i);

                f%f(i) = 2*omega*dsin(lat);

                ulon = u0*dcos(lat);
                ulat = 0d0;

                u%p(i)%v = ulon*mesh%v(i)%c%lonvec + ulat*mesh%v(i)%c%latvec
            enddo

        elseif (j==3) then
            write(resn,*) '../resultA/TC3/'

            m=1; n=1;

            do i=1,mesh%nv
                lon = mesh%v(i)%lon;
                lat = mesh%v(i)%lat;

                h%f(i) = dcos(m*lon)*dcos(n*lat)**4;

                ulon = -m*dcos(n*lat)**4/dcos(lat)*dsin(lon)*dsin(m*lon);
                ulat = -4*n*dcos(n*lat)**3*dsin(n*lat)*dsin(lon)*dsin(m*lon);

                u%p(i)%v = ulon*mesh%v(i)%c%lonvec + ulat*mesh%v(i)%c%latvec

                dfdx = -m*dsin(m*lon)*dcos(n*lat)**4;
                dfdy = dcos(m*lon)*(-4*n*dcos(n*lat)**3*dsin(n*lat));

                dudx = -m*dcos(n*lat)**4/dcos(lat)*(dcos(lon)*dsin(m*lon) + m*dcos(m*lon)*dsin(lon));
                dudy = -m*(-4*n*dcos(n*lat)**3*dsin(n*lat)/dcos(lat) + &
                            dcos(n*lat)**4/dcos(lat)**2*dsin(lat))*dsin(lon)*dsin(m*lon);
                dvdx = -4*n*dcos(n*lat)**3*dsin(n*lat)*(dcos(lon)*dsin(m*lon) + m*dsin(lon)*dcos(m*lon));
                dvdy = -4*n*(-3*n*dcos(n*lat)**2*dsin(n*lat)**2 + &
                            n*dcos(n*lat)**4)*dsin(lon)*dsin(m*lon);

                div%fexact(i) = 1/(rad*dcos(lat))*((dvdy*dcos(lat)-ulat*dsin(lat))+dudx);

                vec = 1d0/(rad*dcos(lat))*dfdx*mesh%v(i)%c%lonvec + 1d0/rad*dfdy*mesh%v(i)%c%latvec;
                gradb%pexact(i)%v = g*vec;

                zeta%fexact(i) = -1d0/(rad*dcos(lat))*(dcos(lat)*dudy-dsin(lat)*ulon-dvdx);

            enddo
            do i=1,mesh%nt
                lon = mesh%tr(i)%c%lon;
                lat = mesh%tr(i)%c%lat;

                htr%fexact(i) = dcos(m*lon)*dcos(n*lat)**4;

                ulon = -m*dcos(n*lat)**4/dcos(lat)*dsin(lon)*dsin(m*lon);
                ulat = -4*n*dcos(n*lat)**3*dsin(n*lat)*dsin(lon)*dsin(m*lon);

                utr%pexact(i)%v = ulon*mesh%tr(i)%c%lonvec + ulat*mesh%tr(i)%c%latvec;
            enddo

        elseif (j==4) then
            write(resn,*) '../resultA/TC4/';
            ! saveenergy=.true.

            alpha = 0*.5d0*pi/2d0;
            u0 = 2d0*pi*rad/(86400d0*12);
            b0 = 2.94d4/g;

            cte = 1d0/g*(rad*omega*u0+u0**2*.5d0)

            do i =1,mesh%nv
                lat = datan2(mesh%v(i)%p(3),dsqrt(mesh%v(i)%p(1)**2+mesh%v(i)%p(2)**2));

                f%f(i)  = 2d0*omega*dsin(lat);

                ! h%f(i) = b0 - cte*(dsin(lat)**2*dcos(alpha));
                bath%f(i) = b0 - cte*(-dcos(lon)*dcos(lat)*dsin(alpha)+dsin(lat)*dcos(alpha))**2;

                ulon = u0*(dcos(lat)*dcos(alpha)+dcos(lon)*dsin(lat)*dsin(alpha));
                ulat = -u0*dsin(lon)*dsin(alpha);

                u%p(i)%v = ulon*mesh%v(i)%c%lonvec + ulat*mesh%v(i)%c%latvec;

                dfdx = 0d0;
                dfdy = -cte*(2d0*dsin(lat)*dcos(lat));! - u0**2/g*(2d0*dsin(lat)*dcos(lat))/2d0;

                dudx = 0d0;
                dudy = -u0*dsin(lat);
                dvdx = 0d0;
                dvdy = 0d0;

                div%fexact(i) = 1/(rad*dcos(lat))*((dvdy*dcos(lat)-ulat*dsin(lat))+dudx);

                vec = 1d0/(rad*dcos(lat))*dfdx*mesh%v(i)%c%lonvec + &
                                    1d0/rad*dfdy*mesh%v(i)%c%latvec;
                gradb%pexact(i)%v = g*vec;

                zeta%fexact(i) = -(dcos(lat)*dudy-dsin(lat)*ulon-dvdx)/(rad*dcos(lat));
            enddo
            h%f = Hm

            h%fexact = h%f;
            do i =1,3
                u%pexact(:)%v(i) = u%p(:)%v(i);
            enddo

            do i =1,mesh%ne
                lat = mesh%edhx(i)%c%lat;
                lon = mesh%edhx(i)%c%lon;
                hed%fexact(i) = b0 - cte*(-dcos(lon)*dcos(lat)*dsin(alpha)+dsin(lon)*dcos(alpha))**2;
            enddo

            do i =1,mesh%nt
                lat = mesh%tr(i)%c%lat;
                lon = mesh%tr(i)%c%lon;
                lat = datan2(mesh%tr(i)%c%p(3),dsqrt(mesh%tr(i)%c%p(1)**2 &
                    +mesh%tr(i)%c%p(2)**2));

                htr%fexact(i) = b0 - cte*(-dcos(lon)*dcos(lat)*dsin(alpha) &
                    +dsin(lat)*dcos(alpha))**2;
                htr%fexact(i) = b0 - cte*(dsin(lat)*dcos(alpha))**2;

                ulon = u0*(dcos(lat)*dcos(alpha)+dcos(lon)*dsin(lat)*dsin(alpha));
                ulat = -u0*dsin(lon)*dsin(alpha);

                utr%pexact(i)%v(:) = ulon*mesh%tr(i)%c%lonvec & 
                                    +ulat*mesh%tr(i)%c%latvec;
            enddo
        elseif (j==8) then
            write(resn,*) '../resultA/TC8/';
            ! saveenergy=.true.

            alpha = 0*.5d0*pi/2d0;
            u0 = 2d0*pi*rad/(86400d0*12);
            b0 = 2.94d4/g;

            cte = 1d0/g*(rad*omega*u0+u0**2*.5d0)

            do i =1,mesh%nv
                lat = datan2(mesh%v(i)%p(3),dsqrt(mesh%v(i)%p(1)**2+mesh%v(i)%p(2)**2));

                f%f(i)  = 2d0*omega*dsin(lat);

                ! h%f(i) = b0 - cte*(dsin(lat)**2*dcos(alpha));
                h%f(i) = b0 - cte*(-dcos(lon)*dcos(lat)*dsin(alpha)+dsin(lat)*dcos(alpha))**2;

                ulon = u0*(dcos(lat)*dcos(alpha)+dcos(lon)*dsin(lat)*dsin(alpha));
                ulat = -u0*dsin(lon)*dsin(alpha);

                u%p(i)%v = ulon*mesh%v(i)%c%lonvec + ulat*mesh%v(i)%c%latvec;

                dfdx = 0d0;
                dfdy = -cte*(2d0*dsin(lat)*dcos(lat));! - u0**2/g*(2d0*dsin(lat)*dcos(lat))/2d0;

                dudx = 0d0;
                dudy = -u0*dsin(lat);
                dvdx = 0d0;
                dvdy = 0d0;

                div%fexact(i) = 1/(rad*dcos(lat))*((dvdy*dcos(lat)-ulat*dsin(lat))+dudx);

                vec = 1d0/(rad*dcos(lat))*dfdx*mesh%v(i)%c%lonvec + &
                                    1d0/rad*dfdy*mesh%v(i)%c%latvec;
                gradb%pexact(i)%v = g*vec;

                zeta%fexact(i) = -(dcos(lat)*dudy-dsin(lat)*ulon-dvdx)/(rad*dcos(lat));
            enddo

            h%fexact = h%f;
            do i =1,3
                u%pexact(:)%v(i) = u%p(:)%v(i);
            enddo

            do i =1,mesh%ne
                lat = mesh%edhx(i)%c%lat;
                lon = mesh%edhx(i)%c%lon;
                hed%fexact(i) = b0 - cte*(-dcos(lon)*dcos(lat)*dsin(alpha)+dsin(lon)*dcos(alpha))**2;
            enddo

            do i =1,mesh%nt
                lat = mesh%tr(i)%c%lat;
                lon = mesh%tr(i)%c%lon;
                lat = datan2(mesh%tr(i)%c%p(3),dsqrt(mesh%tr(i)%c%p(1)**2 &
                    +mesh%tr(i)%c%p(2)**2));

                htr%fexact(i) = b0 - cte*(-dcos(lon)*dcos(lat)*dsin(alpha) &
                    +dsin(lat)*dcos(alpha))**2;
                htr%fexact(i) = b0 - cte*(dsin(lat)*dcos(alpha))**2;

                ulon = u0*(dcos(lat)*dcos(alpha)+dcos(lon)*dsin(lat)*dsin(alpha));
                ulat = -u0*dsin(lon)*dsin(alpha);

                utr%pexact(i)%v(:) = ulon*mesh%tr(i)%c%lonvec & 
                                    +ulat*mesh%tr(i)%c%latvec;
            enddo
        elseif (j==9) then
            write(resn,*) '../resultA/TC9/';
            ! saveenergy=.true.

            call TC9mfield(hgg,mesh)

            nygg = 4*FLOOR(SQRT(REAL(mesh%nv)))
            dygg = pi/nygg

            do i =1,mesh%nv
                p = mesh%v(i)%p

                lat = datan2(p(3),dsqrt(p(1)**2+p(2)**2))
                lon = datan2(p(2),p(1))

                f%f(i) = 2d0*omega*dsin(lat)

                l1 = lat + piby2
                jy = floor(l1/dygg) + 1
                beta = (l1 - (jy - 1)*dygg)/dygg
                if (jy == 1 .or. jy == nygg) then
                    cc2 = 1.0D0 - beta
                    cc3 = beta
                    h%f(i)=(cc2*hgg(jy) + cc3*hgg(jy+1))
                else
                    cc1 = -beta*(beta - 1.0D0) &
                        *(beta - 2.0D0)/6.0D0
                    cc2 = 0.5D0*(beta + 1.0D0) &
                        *(beta - 1.0D0)*(beta - 2.0D0)
                    cc3 = -0.5D0*(beta + 1.0D0) &
                        *beta*(beta - 2.0D0)
                    cc4 = (beta + 1.0D0) &
                        *beta*(beta - 1.0D0)/6.0D0
                    h%f(i) = (cc1*hgg(jy-1) + cc2*hgg(jy) &
                        + cc3*hgg(jy+1) + cc4*hgg(jy+2))
                endif

                lat0 = pi/7d0
                lat1 = pi/2d0-lat0
                en = exp(-4/(lat1 - lat0)**2)
                u00=80
                umen = u00/en
                utmp =0

                den = (lat - lat0)*(lat - lat1)
                if (den < 0.0D0) then
                    utmp = umen*exp(1.0D0/den)
                end if
                call convert_vec_sph2cart(utmp, 0d0, p, vectmp)
                u%p(i)%v = vectmp
    
            enddo

            h0 = 120.0D0
            alpha = 1.0D0/3.0D0
            beta = 1.0D0/15.0D0
            lat2 = 0.5D0*piby2
            do i=1,mesh%nv
                p = mesh%v(i)%p

                lat = datan2(p(3),dsqrt(p(1)**2+p(2)**2))
                lon = datan2(p(2),p(1))

                l1 = lon
                l2 = lat
                clat = dcos(l2)
                !IF (l1 > pi) l1 = l1 - 2.0d0*pi
                e1 = EXP(-(l1/alpha)**2)
                e2 = EXP(-((lat2 - l2)/beta)**2)
                h%f(i) = h%f(i)+h0*clat*e1*e2
        

            enddo


        elseif (j==10) then
            write(resn,*) '../resultA/TC10/';

            do i=1,mesh%nv
                p = mesh%v(i)%p
                p = rotation(p)

                lat = datan2(p(3),dsqrt(p(1)**2+p(2)**2))
                lon = datan2(p(2),p(1))

                call TC10(hf,uf,ff,dfdy,lat)
                ufv = uf*lonvec

                call convert_vec_sph2cart(uf, 0d0, p, ufv)
                ufv = Irotation(ufv)
                u%p(i)%v = ufv;

                h%f(i)   = hf;
                f%f(i)  = ff;
            enddo
        endif
    end subroutine inicond

    !------------------------------------
    subroutine TC9mfield(hgg,mesh)
        implicit none
        type(grid_structure),intent(in) :: mesh
        real*8,allocatable,intent(out) :: hgg(:)
        real*8 :: f1,f2
        real*8 :: l1,l2,den,u1,u2, &
            totvol,totarea,dygg
        real*8,parameter :: u00=80,lat0=pi/7d0
        real*8,parameter :: lat1=pi/2d0-lat0
        real*8,parameter :: en = dexp(-4/(lat1 - lat0)**2)
        real*8,parameter :: umen = u00/en
        integer :: j,nygg


        
        
        nygg = 4*FLOOR(SQRT(REAL(mesh%nv)))
        dygg = pi/nygg
        allocate(hgg(nygg+1))
        hgg(1)=0d0

        totvol = 0.00
        totarea = 0.0D0
        do j = 2, nygg
            l1 = (j-2)*dygg - piby2
            den = (l1 - lat0)*(l1 - lat1)
            if (den < 0.0D0) then
                u1 = umen*exp(1.0D0/den)
            else
                u1 = 0.0D0
            endif
          l2 = (j-1)*dygg - piby2
          den = (l2 - lat0)*(l2 - lat1)
          if (den < 0.0D0) then
              u2 = umen*exp(1.0D0/den)
          else
              u2 = 0.0D0
          endif
          f1 = 2d0*omega*dsin(l1)
          f2 = 2d0*omega*dsin(l2)
          u1 = u1*(f1 + dTAN(l1)*u1/rad)
          u2 = u2*(f2 + dTAN(l2)*u2/rad)
          hgg(j) = hgg(j-1) - rad*0.5d0*(u1 + u2)*dygg

          totarea = totarea + COS(l2)*dygg
          totvol = totvol + hgg(j)*COS(l2)*dygg
        enddo
        hgg(nygg+1) = hgg(nygg)
        totvol = totvol/(totarea*g)
        hgg = hgg + (1.0D4 - totvol)*g !potential phi2
        hgg=hgg/g !Height
    endsubroutine TC9mfield

    subroutine TC10(h,u,f,dhdy,lat)
        real*8,intent(out) :: h,u,f,dhdy
        real*8,intent(in) :: lat
        integer,parameter :: k=160
        integer,parameter :: ni=2d0*k+2d0
        real*8 :: F0,C
        real*8,parameter :: h0=1d5/g
        real*8,parameter :: u0=dsqrt(ni*h0*g)

        f = 2d0*omega

        if (lat>0) then
            h = h0*(2d0-dsin(lat)**ni)
            dhdy = -h0*ni*dsin(lat)**(ni-1)*dcos(lat)/rad*g
            F0 = dcos(lat)/dsin(lat)*f*rad
            C = -g*h0*ni*dsin(lat)**(ni-2)*dcos(lat)**2
            u=(-F0+dsqrt(F0**2-4*C))*.5d0
        else
            h = 2*h0
            u = 0d0
            dhdy = 0
        endif
    endsubroutine TC10

    function rotation(pq) result(pqr)
        real*8 :: pq(3)
        real*8 :: pqr(3)
        real*8 :: Rmat(3,3)
        real*8,parameter :: lat=1*pi/180,&
            lon=3*pi/180

        ! lon = datan2(pq(2),pq(1))
        ! lat = datan2(pq(3),dsqrt(pq(1)**2+pq(2)**2))

        Rmat(1:3,1:3)=0._r8
        Rmat(1,1)=dsin(lat)*dcos(lon)
        Rmat(1,2)=dsin(lat)*dsin(lon)
        Rmat(1,3)=-dcos(lat)
        Rmat(2,1)=-dsin(lon)
        Rmat(2,2)=dcos(lon)
        Rmat(3,1)=dcos(lat)*dcos(lon)
        Rmat(3,2)=dcos(lat)*dsin(lon)
        Rmat(3,3)=dsin(lat)

        ! Rmat = roty1(pi*.5d0)


        pqr = matmul(Rmat,pq)

    end function rotation

    function Irotation(pqr) result(pq)
        real*8 :: pq(3)
        real*8 :: pqr(3)
        real*8 :: Rmat(3,3),RmatT(3,3)
        real*8,parameter :: lat=1*pi/180,&
            lon=3*pi/180

            Rmat(1:3,1:3)=0._r8
            Rmat(1,1)=dsin(lat)*dcos(lon)
            Rmat(1,2)=dsin(lat)*dsin(lon)
            Rmat(1,3)=-dcos(lat)
            Rmat(2,1)=-dsin(lon)
            Rmat(2,2)=dcos(lon)
            Rmat(3,1)=dcos(lat)*dcos(lon)
            Rmat(3,2)=dcos(lat)*dsin(lon)
            Rmat(3,3)=dsin(lat)

        RmatT = transpose(Rmat)
        pq = matmul(RmatT,pqr)

    end function Irotation
    !------------------------------------

    subroutine calc_op(u,utr,h,bath,ke,zeta,div,gradb,uhq_perp,lap,bih,f,mesh)
        implicit none
        type(grid_structure), intent(in)      :: mesh
        type(vector_field_cart),intent(in)    :: u
        type(scalar_field),intent(in)         :: f,h,bath
        type(vector_field_cart),intent(inout) :: gradb,utr,uhq_perp,lap,bih
        type(scalar_field),intent(inout)      :: ke,zeta,div
        type(scalar_field)                    :: phitr,phied,phi
        type(vector_field_cart)               :: uhtr,uh
        integer                               :: i
        real*8                                :: pq(mesh%nv,3), &
            uv(mesh%nv,3),uperpv(mesh%nv,3)

        allocate(phitr%f(mesh%nt),phied%f(mesh%ne));
        allocate(phi%f(mesh%nv));
        allocate(uhtr%p(mesh%nt));
        ! allocate(ke%f(mesh%nv));
        allocate(uh%p(mesh%nv));

        zeta%f =0d0
        do i=1,3
            bih%p(:)%v(i) = 0d0;
            uhq_perp%p(:)%v(i) = 0d0;
            if (nlin) then
                uh%p%v(i)  = mulreal(u%p%v(i),h%f);
            else
                uh%p%v(i)  = u%p%v(i)*1d4
            endif
            ! uh%p%v(i)  = mulreal(u%p%v(i),h%f);
            ! uh%p%v(i)  = u%p%v(i)*h%f;
        enddo

        call vector_ve2tr(uh,uhtr,mesh);
        call div_tr2ve(uhtr,div,mesh);
        !---------------------------------------------
        call ke_energy(u,ke);
            ! phi%f   = g*(h%f + bath%f) + ke%f;
        phi%f   = g*sumreal(h%f, bath%f);
        if (nlin) then
            phi%f   = sumreal(phi%f, ke%f);
            call vector_ve2tr(u,utr,mesh);
            call vort_tr2ve(utr,zeta,mesh);
        endif
        call scalar_ve2tr(phi,phitr,mesh);      
        call scalar_tr2ed(phitr,phied,mesh);
        call grad_tr2ve(phied,phi,gradb,mesh);
        !---------------------------------------------



        !$omp parallel do &
        !$omp default(none) &
        !$omp shared(mesh, u,uv,pq,uperpv) &
        !$omp private(i) &
        !$omp schedule(static)
        do i=1,mesh%nv
            pq(i,:) = mesh%v(i)%p(:);
            uv(i,:) = u%p(i)%v(:);
            uperpv(i,:) = cross_product(pq(i,:),uv(i,:));
        enddo
        !$omp end parallel do

        if (.true.) then
            do i=1,3
                ! uhq_perp%p(:)%v(i) = (zeta%f+f%f)*uperpv(:,i);
                uhq_perp%p(:)%v(i) = sumreal(zeta%f,f%f);
                uhq_perp%p(:)%v(i) = mulreal(uhq_perp%p(:)%v(i),uperpv(:,i));
            enddo
        else
            do i=1,3
                uhq_perp%p(:)%v(i) = 0d0;
            enddo
        endif


        !-------------------------------------------------------------
        if (.true.) then
            do i=1,3
                ! bih%p(:)%v(i) = -1d17*u%p(:)%v(i);
                ! bih%p(:)%v(i) = -1d17*u%p(:)%v(i); ! g5
                ! bih%p(:)%v(i) = -1d13*u%p(:)%v(i);
                ! bih%p(:)%v(i) = -.1d14*u%p(:)%v(i); !g9
                ! bih%p(:)%v(i) = -7.35d14*u%p(:)%v(i); !g8
                ! bih%p(:)%v(i) = -1d16*u%p(:)%v(i); !TC2g7
                ! bih%p(:)%v(i) = -2d17*u%p(:)%v(i); !TC2g6
                ! bih%p(:)%v(i) = -2d18*u%p(:)%v(i); !TC2g5
                ! bih%p(:)%v(i) = -1d19*u%p(:)%v(i); !TC2g4
                ! bih%p(:)%v(i) = -1d20*u%p(:)%v(i); !TC2g3
                ! bih%p(:)%v(i) = -1d21*u%p(:)%v(i); !TC2g2

                ! bih%p(:)%v(i) = 0d0*u%p(:)%v(i); !TC4g6
                ! bih%p(:)%v(i) = -1d12*u%p(:)%v(i); !TC4g6
                ! bih%p(:)%v(i) = -3d14*u%p(:)%v(i); !TC4g6
                ! bih%p(:)%v(i) = -8d14*u%p(:)%v(i); !TC4g6

                ! bih%p(:)%v(i) = -1d17*u%p(:)%v(i); !TC8g6
                ! bih%p(:)%v(i) = -1d18*u%p(:)%v(i); !TC8g5
                ! bih%p(:)%v(i) = -1d19*u%p(:)%v(i); !TC8g4
                ! bih%p(:)%v(i) = -1d20*u%p(:)%v(i); !TC8g3
                ! bih%p(:)%v(i) = -1d22*u%p(:)%v(i); !TC8g2

                ! bih%p(:)%v(i) = -5d15*u%p(:)%v(i); !TC9g7

                ! bih%p(:)%v(i) = -0d0*u%p(:)%v(i); !TC10g6

                ! bih%p(:)%v(i) = -1d22*u%p(:)%v(i); !normg2
                ! bih%p(:)%v(i) = 0d0*u%p(:)%v(i); !g2MODES

                 bih%p(:)%v(i) = 0d0; !TC8g5

            enddo
            call visc2x(bih,lap,mesh);
            call visc2x(lap,bih,mesh);
        endif

    endsubroutine calc_op

    subroutine error_op(u,utr,h,htr,hed,zeta,div,gradb,filef,mesh)
        implicit none
        type(grid_structure), intent(in)      :: mesh
        type(vector_field_cart),intent(in)    :: u
        type(scalar_field),intent(in)         :: h
        type(vector_field_cart),intent(inout) :: gradb,utr
        type(scalar_field),intent(inout)      :: zeta,div,htr,hed
        type(scalar_field)                    :: bath
        type(vector_field_cart)               :: uhq_perp,lap,bih

        type(vector_field_cart)               :: gradu,gradv,uperp

        integer                               :: i
        character(len=100),intent(in)         :: filef
        character(len=100)                    :: filen
        real*8 :: utrnum(mesh%nt,3),utra(mesh%nt,3),errutr(mesh%nt,3), &
                    gradnum(mesh%nv,3),grada(mesh%nv,3),errgrad(mesh%nv,3), &
                    pq(mesh%nv,3),uv(mesh%nv,3),uperpv(mesh%nv,3)

        allocate(bath%f(mesh%nv))
        allocate(uhq_perp%p(mesh%nv));
        allocate(lap%p(mesh%nv));
        allocate(bih%p(mesh%nv));
        allocate(gradu%p(mesh%nt),gradv%p(mesh%nt));
        allocate(uperp%p(mesh%nv));

        bath%f = 0d0
        call scalar_ve2tr(h,htr,mesh);
        call vector_ve2tr(u,utr,mesh);

        ! do i=1,3
        !   utr%p(:)%v(i) = utr%pexact(:)%v(i)
        ! enddo

        call div_tr2ve(utr,div,mesh);
        call vort_tr2ve(utr,zeta,mesh);
        call scalar_tr2ed(htr,hed,mesh);
        ! hed%f = hed%fexact

        call grad_tr2ve(hed,h,gradb,mesh);

        do i=1,3
            pq(:,i) = mesh%v(:)%p(i);
            uv(:,i) = u%p(:)%v(i);
        enddo
        uperpv = cross2d_product(pq,uv);
        do i=1,3
            uperp%p(:)%v(i) = uperpv(:,i)
        enddo

        call grad_vector_ve2tr(u,uperp,gradu,gradv,mesh);

        do i =1,3
            utrnum(:,i) = utr%p(:)%v(i);
            utra(:,i)   = utr%pexact(:)%v(i);
            gradnum(:,i) = g*gradb%p(:)%v(i);
            grada(:,i)   = gradb%pexact(:)%v(i);
        enddo
      
        errutr  = utrnum-utra;
        errgrad = gradnum-grada;

        write(*,*) 'utr:      ', maxval(normm(utrnum)), &
            maxval(normm(utra)), &
            maxval(normm(errutr))

        write(*,*) 'h:        ', maxval(abs(htr%f)), &
            maxval(abs(htr%fexact)), &
            maxval(abs(htr%f-htr%fexact))
        write(*,*) '-------------------------------------------------------------'
        write(*,*) 'div:      ', maxval(abs(div%f)), &
            maxval(abs(div%fexact)), &
            maxval(abs(div%f-div%fexact))
        write(*,*) 'zeta:     ', maxval(abs(zeta%f)), &
            maxval(abs(zeta%fexact)), &
            maxval(abs(zeta%f-zeta%fexact))

        write(*,*) 'gradb:    ', maxval(normm(gradnum)), &
            maxval(normm(grada)), &
            maxval(normm(errgrad))
        write(*,*) '-------------------------------------------------------------'


        write(filen,'(3A)') trim(filef),'gradb.dat';
        open(unit = 100, file = filen(2:100));

        write(filen,'(3A)') trim(filef),'div.dat';
        open(unit = 101, file = filen(2:100));

        write(filen,'(3A)') trim(filef),'zeta.dat';
        open(unit = 102, file = filen(2:100));

        write(filen,'(3A)') trim(filef),'utr.dat';
        open(unit = 103, file = filen(2:100));

        write(filen,'(3A)') trim(filef),'htr.dat';
        open(unit = 104, file = filen(2:100));

        write(filen,'(3A)') trim(filef),'bih.dat';
        open(unit = 105, file = filen(2:100));

        write(filen,'(3A)') trim(filef),'mom.dat';
        open(unit = 106, file = filen(2:100));

        write(filen,'(3A)') trim(filef),'lap.dat';
        open(unit = 107, file = filen(2:100));

        write(filen,'(3A)') trim(filef),'graduvel.dat';
        open(unit = 108, file = filen(2:100));

        do i = 1,mesh%nv
            write(100,*) g*gradb%p(i)%v,gradb%pexact(i)%v
            write(101,*) div%f(i),div%fexact(i)
            write(102,*) zeta%f(i),zeta%fexact(i)
            write(105,*) bih%p(i)%v
            write(107,*) lap%p(i)%v
        enddo

        do i = 1,mesh%nt
            write(103,*) utr%p(i)%v,utr%pexact(i)%v
            write(104,*) htr%f(i),htr%fexact(i)
            write(108,*) gradu%p(i)%v,gradv%p(i)%v
        enddo
      
        do i=100,106
            close(i)
        enddo

    endsubroutine error_op

    subroutine saveoperators(u,h,j,filei,mesh)
        type(grid_structure),intent(in) :: mesh
        type(scalar_field),intent(in) :: h
        type(vector_field_cart),intent(in) :: u
        integer,intent(in) :: j
        integer :: i
        character(len=100),intent(in) :: filei
        character(len=100)            :: filen

        print*, 'Saving',j

        call charsave(filei,filen,'h',j,4)
        open(unit = 110, file = filen);

        call charsave(filei,filen,'u',j,4)
        open(unit = 111, file = filen);

        do i=1,mesh%nv
            write(110,*) h%f(i);
            write(111,*) u%p(i)%v;
        enddo
        do i=110,111
            close(i);
        enddo
    endsubroutine saveoperators

    subroutine charsave(filei,filen,var,j,k)
        character(len=*),intent(in)      :: var
        character(len=100),intent(in)    :: filei
        character(len=100),intent(inout) :: filen
        character(len= 10) :: frmt
        integer, intent(in) :: j,k
        integer :: i,l


        if (j>0) then
            l = floor(log10(j*1d0))+1;
        else
            l = 1;
        endif

        write(frmt,'(A,I1,A)') '(A,I',l,')'

        write(filen,'(5A)') trim(filei(2:100)),'/',trim(var),'/',trim(var);

        do i=1,k-l
            write(filen,'(A,A)') trim(filen),'0';
        enddo
        write(filen,frmt) trim(filen),j;
        write(filen,'(2A)') trim(filen),'.dat';
    endsubroutine charsave

end module swma_operators
