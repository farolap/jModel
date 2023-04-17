module init_cond
    !Use main grid data structures
    use constants
    use datastruct
    use smeshpack
    use basic_funcs
    implicit  none

    contains
        subroutine TC2(u,h,bath,f,lon,lat,lonvec,latvec)
            real*8,intent(in) :: lon,lat,lonvec(3),latvec(3)
            real*8,intent(inout) :: u(3),h,bath,f
            real*8 :: C,r,ulon,ulat
            ! integer :: i
            real*8,parameter :: u0 = 2d1, &
                h0 = 5960d0, b0 = 2000d0, R0 = pi/9d0, &
                lonc = -pi*.5d0, latc = pi/6d0

                ulon =  u0*dcos(lat);
                ulat =  0d0;

                C = 1d0/g*(rad*omega*u0+u0**2*.5d0)
                r   = dsqrt(min(R0**2,(lon-lonc)**2+(lat-latc)**2));

                f = 2*omega*dsin(lat);

                bath   = b0*(1d0-r/R0);
                h = h0-C*dsin(lat)**2 - bath;

                u = ulon*lonvec + ulat*latvec;
        endsubroutine TC2

        subroutine TC3(u,h,zeta,div,gradb,momflux,lon,lat,lonvec,latvec)
            real*8,intent(in) :: lon,lat,lonvec(3),latvec(3)
            real*8,intent(inout) :: u(3),h,zeta,div,gradb(3),momflux(3)
            real*8 :: ulon,ulat,dudx,dudy,dvdx,dvdy,dfdx,dfdy,momlon,momlat,gradlon,gradlat
            ! integer :: i
            real*8,parameter :: m=1,n=1

            ulon = -m*dcos(n*lat)**4/dcos(lat)*dsin(lon)*dsin(m*lon);
            ulat = -4*n*dcos(n*lat)**3*dsin(n*lat)*dsin(lon)*dsin(m*lon);

            dudx = -m*dcos(n*lat)**4/dcos(lat)*(dcos(lon)*dsin(m*lon) &
                + m*dcos(m*lon)*dsin(lon));
            dudy = -m*(-4*n*dcos(n*lat)**3*dsin(n*lat)/dcos(lat) + &
                dcos(n*lat)**4/dcos(lat)**2*dsin(lat))*dsin(lon)*dsin(m*lon);
            dvdx = -4*n*dcos(n*lat)**3*dsin(n*lat)*(dcos(lon)*dsin(m*lon) &
                + m*dsin(lon)*dcos(m*lon));
            dvdy = -4*n*(-3*n*dcos(n*lat)**2*dsin(n*lat)**2 + &
                n*dcos(n*lat)**4)*dsin(lon)*dsin(m*lon);

            zeta = -(dcos(lat)*dudy-dsin(lat)*ulon-dvdx)/(rad*dcos(lat))

            momlon = ulon/(rad*dcos(lat))*dudx + ulat/rad*dudy;
            momlat = ulon/(rad*dcos(lat))*dvdx + ulat/rad*dvdy;
            momflux = momlon*lonvec + momlat*latvec;

            div = 1/(rad*dcos(lat))*((dvdy*dcos(lat)-ulat*dsin(lat))+dudx);

            h = dcos(m*lon)*dcos(n*lat)**4;

            dfdx = -m*dsin(m*lon)*dcos(n*lat)**4;
            dfdy = dcos(m*lon)*(-4*n*dcos(n*lat)**3*dsin(n*lat));

            gradlon = 1d0/(rad*dcos(lat))*dfdx;
            gradlat = 1d0/rad*dfdy;

            gradb = (gradlon*lonvec + gradlat*latvec)*g;
            u = ulon*lonvec + ulat*latvec;
        end subroutine TC3

        subroutine TC8(u,h,zeta,div,gradb,momflux,lon,lat,lonvec,latvec)
            implicit none
            real*8,intent(in) :: lon,lat,lonvec(3),latvec(3)
            real*8,intent(inout) :: u(3),h,zeta,div,gradb(3),momflux(3)
            real*8 :: ulon,ulat,dudx,dudy,dvdx,dvdy,momlon,momlat, &
                dfdx,dfdy,f,gradlon,gradlat
            ! integer :: i
            real*8,parameter :: alpha = 0d0*pi/2d0, &
                u0 = 2d0*pi*rad/(86400d0*12), &
                b0 = 2.94d4/g, &
                C = 1d0/g*(rad*omega*u0+u0**2*.5d0);

            f  = 2d0*omega*(-dcos(lon)*dcos(lat)*dsin(alpha) + dsin(lat)*dcos(alpha));

            ulon =  u0*(dcos(lat)*dcos(alpha)+dcos(lon)*dsin(lat)*dsin(alpha));
            ulat =  -u0*dsin(lon)*dsin(alpha);
            u = ulon*lonvec + ulat*latvec;

            dfdx = -C*2d0*(-dcos(lon)*dcos(lat)*dsin(alpha) &
                + dsin(lat)*dcos(alpha))*(dsin(lon)*dcos(lat)*dsin(alpha));
            dfdy = -C*2d0*(-dcos(lon)*dcos(lat)*dsin(alpha) &
                + dsin(lat)*dcos(alpha))*(dcos(lon)*dsin(lat)*dsin(alpha) &
                +dcos(lat)*dcos(alpha));
            gradlon = 1d0/(rad*dcos(lat))*dfdx;
            gradlat = 1d0/rad*dfdy;
            gradb = g*(gradlon*lonvec + gradlat*latvec);

            dudx   = u0*(-dsin(lon)*dsin(lat)*dsin(alpha));
            dudy   = u0*(-dsin(lat)*dcos(alpha)+dcos(lon)*dcos(lat)*dsin(alpha));
            dvdx   = -u0*dcos(lon)*dsin(alpha);
            dvdy   = 0d0;

            momlon = ulon*dudx/(rad*dcos(lat))+ulat*dudy/rad-ulon*ulat*dtan(lat)/rad;
            momlat = ulon*dvdx/(rad*dcos(lat))+ulat*dvdy/rad+ulon**2*dtan(lat)/rad;

            momflux = momlon*lonvec+momlat*latvec
            h = b0 - C*(-dcos(lon)*dcos(lat)*dsin(alpha) + dsin(lat)*dcos(alpha))**2;
            ! ke = (ulon**2+ulat**2)*.5d0;
            zeta = -1d0/(rad*dcos(lat))*(dcos(lat)*dudy-dsin(lat)*ulon-dvdx);
            div = 1/(rad*dcos(lat))*((dvdy*dcos(lat)-ulat*dsin(lat))+dudx);
        endsubroutine TC8

        subroutine TC9(u,h,f,hgg,lat,pp,mesh)
            implicit none
            type(grid_structure),intent(in) :: mesh
            real*8,intent(out) :: u(3),h,f
            real*8,intent(in) :: hgg(:),lat,pp(3)
            real*8 :: l1,den, &
                utmp,vtmp,vectmp(3), &
                beta,cc1,cc2,cc3,cc4, &
                dygg
            real*8,parameter :: u00=80, &
                lat0=pi/7d0
            real*8,parameter :: lat1=pi/2d0-lat0
            real*8,parameter :: en = exp(-4/(lat1 - lat0)**2)
            real*8,parameter :: umen = u00/en
            
            integer :: jy,nygg

            f = 2d0*omega*dsin(lat)
            nygg = 4*FLOOR(SQRT(REAL(mesh%nv)))
            dygg = pi/nygg


            l1 = lat + piby2
            jy = floor(l1/dygg) + 1
            beta = (l1 - (jy - 1)*dygg)/dygg
            if (jy == 1 .or. jy == nygg) then
                ! Linear interpolation
                cc2 = 1.0D0 - beta
                cc3 = beta
                !phi2(if0) = (cc2*hgg(jy) + cc3*hgg(jy+1))*farea(if0,ngrids)
                h=(cc2*hgg(jy) + cc3*hgg(jy+1)) !modif psp
            else
                ! Cubic interpolation
                cc1 = -beta*(beta - 1.0D0)*(beta - 2.0D0)/6.0D0
                cc2 = 0.5D0*(beta + 1.0D0)*(beta - 1.0D0)*(beta - 2.0D0)
                cc3 = -0.5D0*(beta + 1.0D0)*beta*(beta - 2.0D0)
                cc4 = (beta + 1.0D0)*beta*(beta - 1.0D0)/6.0D0
                h = (cc1*hgg(jy-1) + cc2*hgg(jy) + cc3*hgg(jy+1) + cc4*hgg(jy+2))
            endif

            utmp=0._r8
            vtmp=0._r8

            den = (lat - lat0)*(lat - lat1)
            if (den < 0.0D0) then
                utmp = umen*exp(1.0D0/den)
            end if
            call convert_vec_sph2cart(utmp, vtmp, pp, vectmp)
            u = vectmp
        endsubroutine TC9

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

        subroutine TC10(h,u,f,lat,lonv)
            real*8,intent(out) :: h,u(3),f
            real*8,intent(in) :: lat,lonv(3)
            integer,parameter :: k=160
            integer,parameter :: ni=2d0*k+2d0
            real*8 :: F0,C
            real*8,parameter :: h0=1d5/g
            real*8,parameter :: u0=dsqrt(ni*h0*g)

            f = 2d0*omega

            ! if (abs(lat)<pi/2d0) then
            if (lat>0) then
                h = h0*(2d0-dsin(lat)**ni)
                F0 = dcos(lat)/dsin(lat)*f*rad
                C = g*h0*ni*dsin(lat)**(ni-2)*dcos(lat)**2
                u=(-F+dsqrt(F**2+4*C))*.5d0*lonv
            else
                h = 2*h0
                u = 0d0
            endif
        endsubroutine TC10
end module init_cond
