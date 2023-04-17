module init_condC

    !Use main grid data structures
    use constants
    use swmc_operators
    use datastruct, only: &
        grid_structure, &
        scalar_field, &
        vectorinterpol_methods, &
        vector_field_cart

    use smeshpack
    use basic_funcs

    implicit none

    contains

    subroutine inicond(j,u,h,bath,utr,div,grad,ke,zeta,q,uperp,f,fv,mesh,resn)
        implicit none
        type(grid_structure), intent(in) :: mesh
        type(scalar_field), intent(inout)      :: u,h,div,zeta,grad,uperp,f,fv,q,bath,ke
        type(scalar_field)                     :: eta_ed,h_ed
        type(vector_field_cart), intent(inout) :: utr
        type(vector_field_cart) :: uti
        integer,intent(in) :: j
        integer            :: i,m,n
        real*8             :: v1(3),v2(3),pp(3),pu(3),pq(3),lon,lat,r,h0,R0,lonc,latc, &
            ulat,ulon,alpha, &
            dfdx,dfdy,dudx,dudy,dvdx,dvdy,vec(3),up(3), &
            utrf(3)
        real*8             :: uf,ff,uperpf,gradf,hf,etaf,zetaf,qf,bathf,divf,kef, &
            latvec(3),lonvec(3),nr(3),tg(3),u0,b0
        real*8 :: ufv(3)
        real*8,allocatable :: hgg(:)
        real*8 :: dygg,p(3),l1,cc1,cc2,cc3,cc4,beta,lat0, &
            lat1,en,u00,umen,utmp,den,vectmp(3), &
            lat2,l2,clat,e1,e2
        integer :: nygg,jy
        character(len=100) ::resn

        allocate(eta_ed%fexact(mesh%ne),h_ed%f(mesh%ne));
        allocate(uti%p(mesh%nt));


        bath%f = 10d3;


        if (j==2) then ! mountain
            write(resn,*) '../resultCtri/TC2/'

            u0 = 20d0;
            h0 = 5960d0;
            b0 = 2000d0;
            R0 = pi/9d0;
            lonc = -pi*.5d0;
            latc = pi/6d0;

            do i=1,mesh%nt
                pp = mesh%tr(i)%c%p;
                lon = datan2(pp(2),pp(1));
                lat = datan2(pp(3),dsqrt(pp(1)**2+pp(2)**2));

                r   = dsqrt(min(R0**2,(lon-lonc)**2+(lat-latc)**2));      
                bath%f(i)   = b0*(1-r/R0);
                h%f(i) = h0-(rad*omega*u0+u0**2*.5d0)/g*dsin(lat)**2-bath%f(i);

            enddo

            do i =1,mesh%ne
                pu = mesh%ed(i)%c%p;

                lat = datan2(pu(3),dsqrt(pu(1)**2+pu(2)**2));
                f%f(i) = 2*omega*dsin(lat);

                ulon   = u0*dcos(lat);
                ulat = 0d0;

                v1 = ulon*mesh%ed(i)%c%lonvec + ulat*mesh%ed(i)%c%latvec
                v2 = mesh%ed(i)%nr;

                u%f(i) = dot_product(v1,v2);
            enddo

            do i=1,mesh%nv
                pq = mesh%v(i)%p;
                lat = datan2(pq(3),dsqrt(pq(1)**2+pq(2)**2));
                fv%f(i) = 2*omega*dsin(lat);
            enddo

        else if (j==3) then
            write(resn,*) '../resultCtri/TC3/'

            u0=1;
            h0=2d3;
            m = 1;
            n = 1;
            bath%f = 0d0;
            h%f    = 1d3;
            h%fexact = h%f;

            do i=1,mesh%nt
                pp = mesh%tr(i)%c%p;
                lat = datan2(pp(3),dsqrt(pp(1)**2+pp(2)**2));
                lon = datan2(pp(2),pp(1));
                h%f(i) = h0*dcos(m*lon)*dcos(n*lat)**4;

                ulon = -u0*m*dcos(n*lat)**4/dcos(lat)*dsin(lon)*dsin(m*lon);
                ulat = -u0*4*n*dcos(n*lat)**3*dsin(n*lat)*dsin(lon)*dcos(m*lon);

                dudx    = -u0*m*dcos(n*lat)**4/dcos(lat)*(dcos(lon)*dsin(m*lon) + m*dcos(m*lon)*dsin(lon));
                dvdy    = -u0*4*n*(-3*n*dcos(n*lat)**2*dsin(n*lat)**2 + &
                n*dcos(n*lat)**4)*dsin(lon)*dcos(m*lon);

                div%fexact(i) = 1/(rad*dcos(lat))*((dvdy*dcos(lat)-ulat*dsin(lat))+dudx);
                ke%fexact(i) = (ulon**2+ulat**2)*.5d0;
            enddo

            do i=1,mesh%ne
                lon = mesh%ed(i)%c%lon; ! theta
                lat = mesh%ed(i)%c%lat; ! phi

                ulon = -u0*m*dcos(n*lat)**4/dcos(lat)*dsin(lon)*dsin(m*lon);
                ulat = -u0*4*n*dcos(n*lat)**3*dsin(n*lat)*dsin(lon)*dcos(m*lon);

                u%f(i) = dot_product(ulon*mesh%ed(i)%c%lonvec + &
                ulat*mesh%ed(i)%c%latvec,mesh%ed(i)%nr);

                up = cross_product(mesh%ed(i)%c%p, &
                ulon*mesh%ed(i)%c%lonvec + ulat*mesh%ed(i)%c%latvec);

                uperp%fexact(i) = dot_product(up,mesh%ed(i)%nr);

                dfdx    = -h0*m*dsin(m*lon)*dcos(n*lat)**4;
                dfdy    = h0*dcos(m*lon)*(-4*n*dcos(n*lat)**3*dsin(n*lat));

                vec = (1d0/(rad*dcos(lat))*dfdx)*mesh%ed(i)%c%lonvec + &
                        (1d0/rad*dfdy)*mesh%ed(i)%c%latvec;

                grad%fexact(i) = -g*dot_product(vec,mesh%ed(i)%nr);
            enddo

            do i=1,mesh%nv
                lon = mesh%v(i)%lon;
                lat = mesh%v(i)%lat;        

                ulon = -u0*m*dcos(n*lat)**4/dcos(lat)*dsin(lon)*dsin(m*lon);
                ulat = -u0*4*n*dcos(n*lat)**3*dsin(n*lat)*dsin(lon)*dcos(m*lon);

                dudy = -m*(-4*n*dcos(n*lat)**3*dsin(n*lat)/dcos(lat) + &
                dcos(n*lat)**4/dcos(lat)**2*dsin(lat))*dsin(lon)*dsin(m*lon)*u0;
                dvdx = -4*n*dcos(n*lat)**3*dsin(n*lat)*(dcos(lon)*dcos(m*lon) - m*dsin(lon)*dsin(m*lon))*u0;

                zeta%fexact(i) = -1d0/(rad*dcos(lat))*(dcos(lat)*dudy-dsin(lat)*ulon-dvdx);
            enddo

        else if (j==4) then ! Instability
            write(resn,*) '../resultCtri/TC4/';
            alpha = .0d0;

            h%f   = Hm;
            h_ed%f   = Hm;
            do i=1,mesh%nt
                pp = mesh%tr(i)%c%p;
                lat = datan2(pp(3),dsqrt(pp(1)**2+pp(2)**2));
                lon = datan2(pp(2),pp(1));
                lonvec = mesh%tr(i)%c%lonvec
                latvec = mesh%tr(i)%c%latvec

                call TC8(uf,utrf,hf,ff,bathf,uperpf,divf,gradf, &
                    zetaf,qf,etaf,kef,latvec,lonvec,nr,lat,lon,alpha)
                bath%f(i) = hf;
            enddo

            do i =1,mesh%ne
                pu = mesh%ed(i)%c%p;
                lat = datan2(pu(3),dsqrt(pu(1)**2+pu(2)**2));
                lon = datan2(pu(2),pu(1));

                latvec = mesh%ed(i)%c%latvec;
                lonvec = mesh%ed(i)%c%lonvec;
                nr     = mesh%ed(i)%nr;
                tg     = mesh%ed(i)%tg;

                call TC8(uf,utrf,hf,ff,bathf,uperpf,divf,gradf, &
                    zetaf,qf,etaf,kef,latvec,lonvec,nr,lat,lon,alpha)
                ! call p2TC8(uf,ff,uperpf,gradf,hf,etaf,latvec,lonvec,nr,lat,lon,alpha);
                u%f(i)      = uf;
                f%f(i)      = ff;
            enddo

            do i=1,mesh%nv
                pq = mesh%v(i)%p;
                lat = datan2(pq(3),dsqrt(pq(1)**2+pq(2)**2));
                lon = datan2(pq(2),pq(1));

                fv%f(i)  = ff;
            enddo
        else if (j==8) then ! Nonlinear Geostrophy
            write(resn,*) '../resultCtri/TC8/';
            alpha = .0d0;

            do i=1,mesh%nt
                pp = mesh%tr(i)%c%p;
                lat = datan2(pp(3),dsqrt(pp(1)**2+pp(2)**2));
                lon = datan2(pp(2),pp(1));
                lonvec = mesh%tr(i)%c%lonvec
                latvec = mesh%tr(i)%c%latvec

                call TC8(uf,utrf,hf,ff,bathf,uperpf,divf,gradf, &
                    zetaf,qf,etaf,kef,latvec,lonvec,nr,lat,lon,alpha)
                utr%pexact(i)%v = utrf;
                bath%f(i) = bathf;
                h%f(i)   = hf;
                div%fexact(i) = divf;
                ke%fexact(i)  = kef;
            enddo

            do i =1,mesh%ne
                pu = mesh%ed(i)%c%p;
                lat = datan2(pu(3),dsqrt(pu(1)**2+pu(2)**2));
                lon = datan2(pu(2),pu(1));

                latvec = mesh%ed(i)%c%latvec;
                lonvec = mesh%ed(i)%c%lonvec;
                nr     = mesh%ed(i)%nr;
                tg     = mesh%ed(i)%tg;

                call TC8(uf,utrf,hf,ff,bathf,uperpf,divf,gradf, &
                    zetaf,qf,etaf,kef,latvec,lonvec,nr,lat,lon,alpha)
                ! call p2TC8(uf,ff,uperpf,gradf,hf,etaf,latvec,lonvec,nr,lat,lon,alpha);
                u%f(i)      = uf;
                f%f(i)      = ff;
                uperp%fexact(i)  = uperpf;
                grad%fexact(i)   = gradf;
                h_ed%f(i)   = hf;
                eta_ed%fexact(i) = etaf;
            enddo
            call perot_ed2tr(u,uti,mesh)
            call perot_tr2ed(uti,u,mesh)

            do i=1,mesh%nv
                pq = mesh%v(i)%p;
                lat = datan2(pq(3),dsqrt(pq(1)**2+pq(2)**2));
                lon = datan2(pq(2),pq(1));

                call TC8(uf,utrf,hf,ff,bathf,uperpf,divf,gradf,zetaf,qf,etaf,kef,latvec,lonvec,nr,lat,lon,alpha)
                ! call p3TC8(zetaf,qf,ff,lat,lon,alpha)

                fv%f(i)  = ff;
                zeta%fexact(i) = zetaf
                q%fexact(i) = qf;
            enddo
        elseif (j==9) then
            write(resn,*) '../resultCtri/TC9/';

            call TC9mfield(hgg,mesh)

            nygg = 4*FLOOR(SQRT(REAL(mesh%nv)))
            dygg = pi/nygg

            do i =1,mesh%nt
                p = mesh%tr(i)%c%p

                lat = datan2(p(3),dsqrt(p(1)**2+p(2)**2))
                lon = datan2(p(2),p(1))

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
    
            enddo

            lat0 = pi/7d0
            lat1 = pi/2d0-lat0
            en = exp(-4/(lat1 - lat0)**2)
            u00=80
            umen = u00/en
            do i=1,mesh%ne
                p = mesh%ed(i)%c%p

                lat = datan2(p(3),dsqrt(p(1)**2+p(2)**2))
                lon = datan2(p(2),p(1))
                f%f(i) = 2d0*omega*dsin(lat)

                utmp =0

                den = (lat - lat0)*(lat - lat1)
                if (den < 0.0D0) then
                    utmp = umen*exp(1.0D0/den)
                end if
                call convert_vec_sph2cart(utmp, 0d0, p, vectmp)
                u%f(i) = dot_product(vectmp,mesh%ed(i)%nr)

            enddo

            h0 = 120.0D0
            alpha = 1.0D0/3.0D0
            beta = 1.0D0/15.0D0
            lat2 = 0.5D0*piby2
            do i=1,mesh%nt
                p = mesh%tr(i)%c%p

                lat = datan2(p(3),dsqrt(p(1)**2+p(2)**2))
                lon = datan2(p(2),p(1))

                l1 = lon
                l2 = lat
                clat = dcos(l2)

                e1 = EXP(-(l1/alpha)**2)
                e2 = EXP(-((lat2 - l2)/beta)**2)
                h%f(i) = h%f(i)+h0*clat*e1*e2
        

            enddo

            do i=1,mesh%nv
                p = mesh%v(i)%p

                lat = datan2(p(3),dsqrt(p(1)**2+p(2)**2))
                lon = datan2(p(2),p(1))

                fv%f(i) = 2d0*omega*dsin(lat)
            enddo

        elseif(j==10)then
            write(resn,*) '../resultCtri/TC10/';
    
            do i=1,mesh%nt
                pp = mesh%tr(i)%c%p;
                pp = rotation(pp)

                lat = datan2(pp(3),dsqrt(pp(1)**2+pp(2)**2));
                lon = datan2(pp(2),pp(1));
                lonvec = mesh%tr(i)%c%lonvec
                latvec = mesh%tr(i)%c%latvec

                call TC10(hf,uf,ff,dfdy,lat)
                h%f(i)   = hf;

            enddo

            do i =1,mesh%ne
                pu = mesh%ed(i)%c%p;
                pu = rotation(pu)

                lat = datan2(pu(3),dsqrt(pu(1)**2+pu(2)**2));
                lon = datan2(pu(2),pu(1));

                latvec = mesh%ed(i)%c%latvec;
                lonvec = mesh%ed(i)%c%lonvec;

                nr     = mesh%ed(i)%nr;
                tg     = mesh%ed(i)%tg;

                call TC10(hf,uf,ff,dfdy,lat)
                ufv = uf*lonvec

                call convert_vec_sph2cart(uf, 0d0, pu, ufv)
                ufv = Irotation(ufv)
                u%f(i)      = dot_product(ufv,mesh%ed(i)%nr);

                call convert_vec_sph2cart(0d0, dfdy, pu, ufv)
                ufv = Irotation(ufv)
                grad%fexact(i) = dot_product(ufv,nr)

                uperp%fexact(i) = -dot_product(ufv,mesh%ed(i)%tg);
                f%f(i)      = ff;

            enddo

            do i=1,mesh%nv
                pq = mesh%v(i)%p;
                lat = datan2(pq(3),dsqrt(pq(1)**2+pq(2)**2));
                lon = datan2(pq(2),pq(1));

                call TC10(hf,uf,ff,dfdy,lat)

                fv%f(i)  = ff;
            enddo
        endif

        return    
    end subroutine inicond

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
    !TC2
    subroutine p1TC2(h,bath,div,ke,lat,lon)
      implicit none
      real*8,intent(inout) :: h,div,ke,bath
      real*8,intent(in)    :: lat,lon
      real*8               :: r
      real*8,parameter :: u0 = 20d0, &
        h0 = 5960d0, &
        b0 = 2000d0, &
        R0 = pi/9d0, &
        lonc = -pi*.5d0, &
        latc = pi/6d0


        r   = dsqrt(min(R0**2,(lon-lonc)**2+(lat-latc)**2));      
        bath   = b0*(1-r/R0);
        h = h0-(rad*omega*u0+u0**2*.5d0)/g*dsin(lat)**2-bath;
        div = 0d0;
        ke = 0d0;
    endsubroutine p1TC2

    subroutine p2TC2(u,f,uperp,grad,h,eta,latvec,lonvec,nr,lat,lon)
      implicit none
      real*8,intent(inout) :: u,f,uperp,grad,h,eta
      real*8,intent(in)  :: latvec(3),lonvec(3),nr(3),lat,lon
      real*8 :: vec(3),ulon,ulat
      real*8,parameter :: u0 = 20d0, &
        h0 = 5960d0, &
        b0 = 2000d0, &
        R0 = pi/9d0, &
        lonc = -pi*.5d0, &
        latc = pi/6d0

      f = 2*omega*dsin(lat);

      ulon   = u0*dcos(lat);
      ulat = 0d0;

      vec = ulon*lonvec + ulat*latvec;

      u = dot_product(vec,nr);
      uperp = 0d0;
      grad = 0d0;
      h = 0d0;
      eta = 0d0;

    endsubroutine p2TC2

    subroutine p3TC2(zeta,q,f,lat,lon)
      implicit none
      real*8,intent(inout) :: zeta,q,f
      real*8,intent(in)  :: lat,lon
      real*8 :: dudy,dvdx,h,ulon,ulat
      integer :: m=3,n=3
        
      ulon = -m*dcos(n*lat)**4/dcos(lat)*dsin(lon)*dsin(m*lon);
      ulat = -4*n*dcos(n*lat)**3*dsin(n*lat)*dsin(lon)*dsin(m*lon);


      f = 0d0;

      dudy = -m*(-4*n*dcos(n*lat)**3*dsin(n*lat)/dcos(lat) + &
                  dcos(n*lat)**4/dcos(lat)**2*dsin(lat))*dsin(lon)*dsin(m*lon);
      dvdx = -4*n*dcos(n*lat)**3*dsin(n*lat)*(dcos(lon)*dsin(m*lon) + m*dsin(lon)*dcos(m*lon));



      zeta = -1d0/(rad*dcos(lat))*(dcos(lat)*dudy-dsin(lat)*ulon-dvdx);

      h = dcos(m*lon)*dcos(n*lat)**4;

      q = (zeta+f)/h;
      
    endsubroutine p3TC2

    !TC3
    subroutine p1TC3(h,bath,div,ke,lat,lon)
      implicit none
      real*8,intent(inout) :: h,div,ke,bath
      real*8,intent(in)    :: lat,lon
      real*8               :: ulat,ulon,dudx,dvdy
      integer,parameter :: m=1, &
        n=1

      ulon = -m*dcos(n*lat)**4/dcos(lat)*dsin(lon)*dsin(m*lon);
      ulat = -4*n*dcos(n*lat)**3*dsin(n*lat)*dsin(lon)*dsin(m*lon);

      dudx = -m*dcos(n*lat)**4/dcos(lat)*(dcos(lon)*dsin(m*lon) + m*dcos(m*lon)*dsin(lon));

      dvdy = -4*n*(-3*n*dcos(n*lat)**2*dsin(n*lat)**2 + &
        n*dcos(n*lat)**4)*dsin(lon)*dsin(m*lon);


      h  = dcos(m*lon)*dcos(n*lat)**4;
      ke  = (ulon**2+ulat**2)*.5d0;

      div = 1/(rad*dcos(lat))*((dvdy*dcos(lat)-ulat*dsin(lat))+dudx);

      bath = 0d0;

    endsubroutine p1TC3

    subroutine p2TC3(u,f,uperp,grad,h,eta,latvec,lonvec,nr,lat,lon)
      implicit none
      real*8,intent(inout) :: u,f,uperp,grad,h,eta
      real*8,intent(in)  :: latvec(3),lonvec(3),nr(3),lat,lon
      real*8 :: dfdx,dfdy,dvdx,dudy,vec(3),ulon,ulat
      integer,parameter :: m=1, &
        n=1

      ulon = -m*dcos(n*lat)**4/dcos(lat)*dsin(lon)*dsin(m*lon);
      ulat = -4*n*dcos(n*lat)**3*dsin(n*lat)*dsin(lon)*dsin(m*lon);

      f  = 0d0;

      vec = ulon*lonvec + ulat*latvec
      u = dot_product(vec,nr);

      vec = ulon*latvec - ulat*lonvec;
      uperp = dot_product(vec,nr);

      dfdx = -m*dsin(m*lon)*dcos(n*lat)**4;
      dfdy = dcos(m*lon)*(-4*n*dcos(n*lat)**3*dsin(n*lat));

      vec = 1d0/(rad*dcos(lat))*dfdx*lonvec + 1d0/rad*dfdy*latvec;

      grad = -dot_product(vec,nr)*g;

      h  = dcos(m*lon)*dcos(n*lat)**4;

      dudy = -m*(-4*n*dcos(n*lat)**3*dsin(n*lat)/dcos(lat) + &
        dcos(n*lat)**4/dcos(lat)**2*dsin(lat))*dsin(lon)*dsin(m*lon);

      dvdx = -4*n*dcos(n*lat)**3*dsin(n*lat)*(dcos(lon)*dsin(m*lon) + m*dsin(lon)*dcos(m*lon));

      eta = -1d0/(rad*dcos(lat))*(dcos(lat)*dudy-dsin(lat)*ulon-dvdx);

    endsubroutine p2TC3

    subroutine p3TC3(zeta,q,f,lat,lon)
      implicit none
      real*8,intent(inout) :: zeta,q,f
      real*8,intent(in)  :: lat,lon
      real*8 :: dudy,dvdx,h,ulon,ulat
      integer :: m=1,n=1
        
      ulon = -m*dcos(n*lat)**4/dcos(lat)*dsin(lon)*dsin(m*lon);
      ulat = -4*n*dcos(n*lat)**3*dsin(n*lat)*dsin(lon)*dsin(m*lon);


      f = 0d0;

      dudy = -m*(-4*n*dcos(n*lat)**3*dsin(n*lat)/dcos(lat) + &
                  dcos(n*lat)**4/dcos(lat)**2*dsin(lat))*dsin(lon)*dsin(m*lon);
      dvdx = -4*n*dcos(n*lat)**3*dsin(n*lat)*(dcos(lon)*dsin(m*lon) + m*dsin(lon)*dcos(m*lon));



      zeta = -1d0/(rad*dcos(lat))*(dcos(lat)*dudy-dsin(lat)*ulon-dvdx);

      h = dcos(m*lon)*dcos(n*lat)**4;

      q = (zeta+f)/h;
      
    endsubroutine p3TC3

    !TC8--------------------------------------------------------------------
    subroutine p1TC8(h,bath,div,ke,lat,lon,alpha)
      implicit none
      real*8,intent(inout) :: h,div,ke,bath
      real*8,intent(in)    :: lat,lon,alpha
      real*8               :: ulat,ulon
      real*8,parameter     :: u0 = 2d0*pi*rad/(86400d0*12), &
        b0  = 2.94d4/g, &
        cte = 1d0/g*(rad*omega*u0+u0**2*.5d0)

        ulon =  u0*(dcos(lat)*dcos(alpha)+dcos(lon)*dsin(lat)*dsin(alpha));
        ulat = -u0*dsin(lon)*dsin(alpha)


        h = b0 - cte*(-dcos(lon)*dcos(lat)*dsin(alpha)+dsin(lat)*dcos(alpha))**2;
        ke  = (ulon**2+ulat**2)*.5d0;
        div = 0d0;
        bath = 0d0;

    endsubroutine p1TC8

    subroutine p2TC8(u,f,uperp,grad,h,eta,latvec,lonvec,nr,lat,lon,alpha)
      implicit none
      real*8,intent(inout) :: u,f,uperp,grad,h,eta
      real*8,intent(in)  :: latvec(3),lonvec(3),nr(3),lat,lon,alpha
      real*8 :: var,dfdx,dfdy,dvdx,dudy,vec(3),ulon,ulat
      real*8,parameter   :: u0 = 2d0*pi*rad/(86400d0*12), &
        b0  = 2.94d4/g, &
        cte = 1d0/g*(rad*omega*u0+u0**2*.5d0)

        ulon =  u0*(dcos(lat)*dcos(alpha)+dcos(lon)*dsin(lat)*dsin(alpha));
        ulat = -u0*dsin(lon)*dsin(alpha)

        f  = 2d0*omega*(-dcos(lon)*dcos(lat)*dsin(alpha)+dsin(lat)*dcos(alpha));

        vec = ulon*lonvec + ulat*latvec
        u = dot_product(vec,nr);

        vec = ulon*latvec - ulat*lonvec;
        uperp = dot_product(vec,nr);

        var   = 2d0*cte*(dsin(lat)*dcos(alpha)-dcos(lon)*dcos(lat)*dsin(alpha))
        dfdy = var/rad*(dcos(lon)*dsin(lat)*dsin(alpha)+dcos(lat)*dcos(alpha));
        dfdx  = var/(rad*dcos(lat))*dsin(lon)*dcos(lat)*dsin(alpha);

        vec = dfdx*lonvec + dfdy*latvec;

        grad = -dot_product(vec,nr)*g;

        h = b0 - cte*(-dcos(lon)*dcos(lat)*dsin(alpha)+dsin(lat)*dcos(alpha))**2;

        dudy = u0*(dcos(lon)*dcos(lat)*dsin(alpha)-dsin(lat)*dcos(alpha));
        dvdx = -u0*dcos(lon)*dsin(alpha);
        var = -1d0/(rad*dcos(lat))*(dcos(lat)*dudy-dsin(lat)*ulon-dvdx);
        eta = (var+f);

    endsubroutine p2TC8

    subroutine p3TC8(zeta,q,f,lat,lon,alpha)
      implicit none
      real*8,intent(inout) :: zeta,q,f
      real*8,intent(in)  :: lat,lon,alpha
      real*8 :: dudy,dvdx,h,ulon,ulat
      real*8,parameter   :: u0 = 2d0*pi*rad/(86400d0*12), &
        b0  = 2.94d4/g, &
        cte = 1d0/g*(rad*omega*u0+u0**2*.5d0)

        ulon =  u0*(dcos(lat)*dcos(alpha)+dcos(lon)*dsin(lat)*dsin(alpha));
        ulat = -u0*dsin(lon)*dsin(alpha);

        f = 2d0*omega*(-dcos(lon)*dcos(lat)*dsin(alpha)+dsin(lat)*dcos(alpha));

        dudy = u0*(dcos(lon)*dcos(lat)*dsin(alpha)-dsin(lat)*dcos(alpha));
        dvdx = -u0*dcos(lon)*dsin(alpha);
        zeta = -1d0/(rad*dcos(lat))*(dcos(lat)*dudy-dsin(lat)*ulon-dvdx);

        h = b0 - cte*(dcos(lon)*dcos(lat)*dsin(alpha)+dsin(lat)*dcos(alpha))**2;

        q = (zeta+f)/h;
      
    endsubroutine p3TC8

    subroutine TC8(u,uvec,h,f,bath,uperp,div,grad,zeta,q,eta,ke, &
            latvec,lonvec,nr,lat,lon,alpha)
      implicit none
      real*8,intent(inout) :: zeta,f,u,uperp,grad,h,q,eta,bath,div,ke,uvec(3)
      real*8,intent(in)  :: lat,lon,alpha,latvec(3),lonvec(3),nr(3)
      real*8 :: dudy,dvdx,ulon,ulat,var,dfdx,dfdy,vec(3)
      real*8,parameter   :: u0 = 2d0*pi*rad/(86400d0*12), &
        b0  = 2.94d4/g, &
        cte = 1d0/g*(rad*omega*u0+u0**2*.5d0)
  
        ulon =  u0*(dcos(lat)*dcos(alpha)+dcos(lon)*dsin(lat)*dsin(alpha));
        ulat = -u0*dsin(lon)*dsin(alpha);

        f = 2d0*omega*(-dcos(lon)*dcos(lat)*dsin(alpha)+dsin(lat)*dcos(alpha));

        dudy = u0*(dcos(lon)*dcos(lat)*dsin(alpha)-dsin(lat)*dcos(alpha));
        dvdx = -u0*dcos(lon)*dsin(alpha);
        zeta = -1d0/(rad*dcos(lat))*(dcos(lat)*dudy-dsin(lat)*ulon-dvdx);

        h = b0 - cte*(dcos(lon)*dcos(lat)*dsin(alpha)+dsin(lat)*dcos(alpha))**2;

        uvec = ulon*lonvec + ulat*latvec
        u = dot_product(uvec,nr);


        vec = ulon*latvec - ulat*lonvec;
        uperp = dot_product(vec,nr);

        dfdy = -cte/rad*2d0*(dcos(lon)*dcos(lat)*dsin(alpha)+dsin(lat)*dcos(alpha)) &
            *(-dcos(lon)*dsin(lat)*dsin(alpha)+dcos(lat)*dcos(alpha))
        ! dfdy = cte/rad*(dcos(lon)*dcos(lat)*dsin(alpha)-dsin(lat)*dcos(alpha));
        dfdx  = var/(rad*dcos(lat))*dsin(lon)*dcos(lat)*dsin(alpha);

        vec = dfdx*lonvec + dfdy*latvec;
        grad = -dot_product(vec,nr)*g;

        ke  = (ulon**2+ulat**2)*.5d0;
        div = 0d0;
        bath = 0d0;

        eta = (zeta+f);
        q = eta/h;
    endsubroutine TC8
    !-----------------------------------------------------------------------

    !TC10-------------------------------------------------------------------
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
    !-----------------------------------------------------------------------

end module init_condC