module swmb_operators
    !Use main grid data structures
    use constants
    use datastruct, only: &
        grid_structure, &
        scalar_field, &
        vector_field_cart

    use smeshpack
    ! use aux_funcs

  contains

    subroutine scalar_ve2tr(fve,ftr,mesh)
        implicit none
        type(grid_structure), intent(in)  :: mesh
        type(scalar_field), intent(in)    :: fve
        type(scalar_field), intent(inout) :: ftr
        integer ::i

        !$omp parallel do &
        !$omp default(none) &
        !$omp shared(mesh,fve,ftr) &
        !$omp private(i) &
        !$omp schedule(static)
        do i =1,mesh%nt
            ftr%f(i) = (fve%f(mesh%tr(i)%v(1)) &
                + fve%f(mesh%tr(i)%v(2)) &
                + fve%f(mesh%tr(i)%v(3)))/3d0;
        enddo
        !$omp end parallel do

        ! ftr%f = (fve%f(mesh%tr(:)%v(1)) &
        !     + fve%f(mesh%tr(:)%v(2)) &
        !     + fve%f(mesh%tr(:)%v(3)))/3d0;
    end subroutine scalar_ve2tr

    subroutine scalar_ve2ed(fve,fed,mesh)
        implicit none
        type(grid_structure), intent(in)  :: mesh
        type(scalar_field), intent(in)    :: fve
        type(scalar_field), intent(inout) :: fed
        integer :: i

        !$omp parallel do &
        !$omp default(none) &
        !$omp shared(mesh,fve,fed) &
        !$omp private(i) &
        !$omp schedule(static)
        do i=1,mesh%ne
            fed%f(i) = (fve%f(mesh%ed(i)%v(1)) &
                + fve%f(mesh%ed(i)%v(2)))*.5d0;
        enddo
        !$omp end parallel do

        ! fed%f = (fve%f(mesh%ed(:)%v(1)) &
        !     + fve%f(mesh%ed(:)%v(2)))*.5d0;


    end subroutine scalar_ve2ed

    subroutine vector_ve2tr(fve,ftr,mesh)
        implicit none
        type(grid_structure), intent(in)  :: mesh
        type(vector_field_cart), intent(in)    :: fve
        type(vector_field_cart), intent(inout) :: ftr
        integer :: i

        !$omp parallel do &
        !$omp default(none) &
        !$omp shared(mesh,fve,ftr) &
        !$omp private(i) &
        !$omp schedule(static)
        do i = 1,mesh%nt
            ftr%p(i)%v = (fve%p(mesh%tr(i)%v(1))%v &
                + fve%p(mesh%tr(i)%v(2))%v &
                + fve%p(mesh%tr(i)%v(3))%v)/3d0;
            ftr%p(i)%v(:) = proj_vec_sphere(ftr%p(i)%v(:),mesh%tr(i)%b%p)
        enddo
        !$omp end parallel do

        ! ftr%p(:)%v(i) = (fve%p(mesh%tr(:)%v(1))%v(i) &
        !     + fve%p(mesh%tr(:)%v(2))%v(i) &
        !     + fve%p(mesh%tr(:)%v(3))%v(i))/3d0;

    end subroutine vector_ve2tr

    subroutine grad_ve2ed(fu,grad,mesh)
        implicit none
        type(grid_structure), intent(in)      :: mesh
        type(scalar_field), intent(in)        :: fu
        type(scalar_field), intent(inout)     :: grad
        
        grad%f = (fu%f(mesh%ed(:)%v(2))-fu%f(mesh%ed(:)%v(1)))/(mesh%ed(:)%leng*rad)

    end subroutine grad_ve2ed

    subroutine grad_tr2ed(fu,grad,mesh)
        implicit none
        type(grid_structure), intent(in)      :: mesh
        type(scalar_field), intent(in)        :: fu
        type(scalar_field), intent(inout)     :: grad
        
        grad%f = (fu%f(mesh%ed(:)%sh(2))-fu%f(mesh%ed(:)%sh(1)))/(mesh%edhx(:)%lengc*rad)

    end subroutine grad_tr2ed

    subroutine div_tr2ve(fu,div,mesh)
        implicit none
        type(grid_structure), intent(in)      :: mesh
        type(vector_field_cart), intent(in)   :: fu
        type(scalar_field), intent(inout)     :: div
        integer :: i,j,k,ed,sh
        real*8  :: nrs(3),Ai,le,fint(3),vecd

        div%f = 0d0

        !$omp parallel do &
        !$omp default(none) &
        !$omp shared(mesh, fu, div) &
        !$omp private(i,j,k,ed,sh,nrs,Ai,le,fint,vecd) &
        !$omp schedule(static)
        do i =1,mesh%nv
            Ai = mesh%hx(i)%cdareag*rad;
            do j=1,mesh%v(i)%nnb
                ed = mesh%v(i)%ed(j);
                do k = 1,2            
                    sh   = mesh%ed(ed)%sh(k);
                    le   = mesh%edhx(ed)%lend(k);
                    nrs  = mesh%edhx(ed)%nrs(k,:)*mesh%hx(i)%cord(j,k);

                    fint = fu%p(sh)%v(:);
                    vecd = dot_product(fint,nrs)

                    div%f(i) = div%f(i) + vecd*le;
                enddo
            enddo
            div%f(i) = div%f(i)/Ai
        enddo
        !$omp end parallel do
    end subroutine div_tr2ve

    subroutine vort_tr2ve(fu,zeta,mesh)
        implicit none
        type(grid_structure), intent(in)      :: mesh
        type(vector_field_cart), intent(in)  :: fu
        type(scalar_field), intent(inout)     :: zeta
        integer :: i,j,k,ed,sh
        real*8  :: tgd(3),Ai,le,fint(3),vecd

        zeta%f = 0d0

        !$omp parallel do &
        !$omp default(none) &
        !$omp shared(mesh, fu, zeta) &
        !$omp private(i,j,k,ed,sh,tgd,Ai,le,fint,vecd) &
        !$omp schedule(static)
        do i =1,mesh%nv
            Ai = mesh%hx(i)%cdareag*rad;
            do j=1,mesh%v(i)%nnb
                ed = mesh%v(i)%ed(j);
                do k = 1,2
                    sh   = mesh%ed(ed)%sh(k);
                    le   = mesh%edhx(ed)%lend(k);
                    tgd  = mesh%edhx(ed)%tgs(k,:)*mesh%hx(i)%cord(j,k);

                    fint = fu%p(sh)%v(:);
                    vecd = dot_product(fint,tgd)

                    zeta%f(i) = zeta%f(i) + vecd*le;
                enddo
            enddo
            zeta%f(i) = zeta%f(i)/Ai
        enddo
        !$omp end parallel do
    end subroutine vort_tr2ve

    subroutine grad_ed2tr(he,grad,mesh)
        implicit none
        type(grid_structure), intent(in)       :: mesh
        type(scalar_field), intent(in)         :: he
        type(vector_field_cart), intent(inout) :: grad
        integer :: i,j,ed
        real*8  :: A,le,nes(3),vec(3)

        do i=1,3
            grad%p(:)%v(i) = 0d0;
        enddo

        !$omp parallel do &
        !$omp default(none) &
        !$omp shared(mesh, he, grad) &
        !$omp private(i,j,ed,A,le,nes) &
        !$omp private(vec) &
        !$omp schedule(static)
        do i= 1,mesh%nt
            A    = mesh%tr(i)%areag*rad;
            vec = 0d0;
            do j= 1,3
                ed   = mesh%tr(i)%ed(j);

                le  = mesh%ed(ed)%leng;
                nes = mesh%ed(ed)%nr;

                vec = vec + le*nes*he%f(ed)*mesh%tr(i)%cor(j);
            enddo
            grad%p(i)%v(:) = proj_vec_sphere(vec,mesh%tr(i)%b%p)/A;
        enddo
        !$omp end parallel do
    end subroutine grad_ed2tr

    subroutine ke_tr2ve(u,ke,mesh)
        implicit None
        type(grid_structure), intent(in)        :: mesh
        type(vector_field_cart), intent(in)    :: u
        type(scalar_field), intent(inout)       :: ke
        integer :: i,j,idt
        real*8  :: Ax,Ac
      

        ke%f = 0d0

        !$omp parallel do &
        !$omp default(none) &
        !$omp shared(mesh, u, ke) &
        !$omp private(i,j,idt,Ax,Ac) &
        !$omp schedule(static)
        do i = 1,mesh%nv
            Ax = 0d0;
            do j = 1,mesh%v(i)%nnb
                idt = mesh%v(i)%tr(j);
                Ac = mesh%tr(idt)%areag;

                Ax = Ax +Ac/3d0;
        
                ke%f(i) = ke%f(i) + Ac*norm(u%p(idt)%v)**2/(3d0*2d0);
            enddo
            ke%f(i) = ke%f(i)/Ax;
        enddo
        !$omp end parallel do

    end subroutine ke_tr2ve

    subroutine vectinv_ve2tr(zeta,fv,u,uhq_perp,mesh)
        implicit none
        type(grid_structure),     intent(in)    :: mesh
        type(scalar_field), intent(in)          :: zeta,fv
        type(vector_field_cart), intent(in)     :: u
        type(vector_field_cart), intent(inout)  :: uhq_perp
        real*8 :: uv(mesh%nt,3),uperp(mesh%nt,3),kunit(mesh%nt,3)
        real*8  :: tmp(mesh%nt)
        integer :: i

        uperp = 0d0;
        do i=1,3
            uhq_perp%p(:)%v(i) = 0d0;
            uv(:,i) = u%p(:)%v(i);
            kunit(:,i) = mesh%tr(:)%b%p(i);
        enddo

        uperp = cross2d_product(kunit,uv);

        tmp = ((zeta%f(mesh%tr(:)%v(1))+fv%f(mesh%tr(:)%v(1))) + &
            (zeta%f(mesh%tr(:)%v(2))+fv%f(mesh%tr(:)%v(2))) + &
            (zeta%f(mesh%tr(:)%v(3))+fv%f(mesh%tr(:)%v(3))))/3d0;

        do i=1,3
            uhq_perp%p(:)%v(i) = tmp*uperp(:,i);
        enddo
    end subroutine vectinv_ve2tr

    subroutine coriolisterm(f,u,fuperp,mesh)
        implicit none
        type(grid_structure),     intent(in)    :: mesh
        type(scalar_field), intent(in)          :: f
        type(vector_field_cart), intent(in)     :: u
        type(vector_field_cart), intent(inout)  :: fuperp
        real*8                                  :: uv(mesh%nt,3),uperp(mesh%nt,3),kunit(mesh%nt,3)
        integer :: i

        !$omp parallel do &
        !$omp default(none) &
        !$omp shared(mesh, u,fuperp,uv,kunit,uperp,f) &
        !$omp private(i) &
        !$omp schedule(static)
        do i =1,mesh%nt
            uv(i,:) = u%p(i)%v(:);
            kunit(i,:) = mesh%tr(i)%b%p(:);
            uperp(i,:) = cross_product(kunit(i,:),uv(i,:));
            fuperp%p(i)%v(:) = f%f(i)*uperp(i,:);
        enddo
        !$omp end parallel do

    end subroutine coriolisterm

    subroutine divscalcontrvolume(uh,u,momflux,div,mesh)
        implicit none
        type(grid_structure), intent(in)       :: mesh
        type(vector_field_cart), intent(in)    :: uh,u
        type(vector_field_cart), intent(inout) :: momflux
        type(scalar_field), intent(inout)      :: div
        real*8 :: momu,momv, &
            lonvecqi(3),latvecqi(3),tgs(3),nrs(3),nr(3), &
            pqi(3),ppci(3),pei(3),peppci(3), &
            dei,cor, &
            unri,uhnri,upi(3),uhpi(3),ui,vi
        integer :: i,j,k,vedii,eshi
        real*8  :: Av

        do i=1,3
            momflux%p(:)%v(i) =0d0;
        enddo
        div%f =0d0;

        !$omp parallel do &
        !$omp default(none) &
        !$omp shared(mesh, u,uh, div,momflux) &
        !$omp private(i,j,k,vedii,eshi) &
        !$omp private(pqi,pei,ppci,peppci) &
        !$omp private(Av,lonvecqi,latvecqi,nrs,tgs,nr) &
        !$omp private(unri,uhnri,upi,uhpi,ui,vi,momu,momv) &
        !$omp private(cor,dei) &
        !$omp schedule(static)
        do i =1,mesh%nv
            momu = 0d0
            momv = 0d0
            Av = mesh%hx(i)%cdareag*rad;
            pqi = mesh%v(i)%p
            lonvecqi = mesh%v(i)%c%lonvec
            latvecqi = mesh%v(i)%c%latvec
            do j =1,mesh%v(i)%nnb
                vedii = mesh%v(i)%ed(j)
                pei = mesh%ed(vedii)%c%p
                do k=1,2
                    eshi = mesh%ed(vedii)%sh(k)
                    upi = u%p(eshi)%v
                    uhpi = uh%p(eshi)%v
                    ppci = mesh%tr(eshi)%b%p
                    peppci = ppci+pei; peppci = peppci/norm(peppci)
                    tgs = pei-ppci; tgs = tgs/norm(tgs)
                    nrs = -cross_product(peppci,tgs)
                    nr = pei-pqi; nr = nr/norm(nr)
                    cor = dsign(1d0,sum(nr*nrs))
                    dei = arclen(ppci,pei)
                    unri = cor*sum(upi*nrs)
                    uhnri = cor*sum(uhpi*nrs)

                    ui = sum(upi*lonvecqi)
                    vi = sum(upi*latvecqi)

                    momu = momu + dei*unri*ui
                    momv = momv + dei*unri*vi
                    div%f(i) = dei*unri
                enddo
            enddo
            momflux%p(i)%v = (momu*lonvecqi + &
                momv*latvecqi)/Av
            div%f(i) = div%f(i)/Av
        enddo
        !$omp end parallel do

    end subroutine divscalcontrvolume

    subroutine spherical_correction(u,momflux,mesh)
        type(grid_structure),intent(in) :: mesh
        type(vector_field_cart),intent(in) :: u
        type(vector_field_cart),intent(inout) :: momflux
        real*8 :: ulon,uperp(3),M
        integer ::i

        !$omp parallel do &
        !$omp default(none) &
        !$omp shared(mesh, momflux,u) &
        !$omp private(i,ulon,uperp,M) &
        !$omp schedule(static)
        do i =1,mesh%nt
            ulon = dot_product(u%p(i)%v,mesh%tr(i)%b%lonvec)
            uperp = cross_product(u%p(i)%v,mesh%tr(i)%b%p)
            M = ulon*dtan(mesh%tr(i)%b%lat)/rad
            momflux%p(i)%v = momflux%p(i)%v - M*uperp;
        enddo
        !$omp end parallel do
    endsubroutine spherical_correction

    subroutine visc2x_tr2tr(u,lap,tau,mesh)
        implicit none
        type(grid_structure),     intent(in)   :: mesh
        type(vector_field_cart), intent(in)    :: u
        type(vector_field_cart), intent(inout) :: lap
        real*8,intent(in)                      :: tau
        integer :: i,j,ide,idt(2)
        real*8  :: cor1,cor2,vecr(3),vecr1(3),vecr2(3),A,ui(2),vi(2),nr(3),tg(3)

        do i=1,3
            lap%p(:)%v(i) = 0d0;
            lap%p(:)%v(i) = 0d0;
        enddo

        !$omp parallel do &
        !$omp default(none) &
        !$omp shared(mesh, u, lap) &
        !$omp private(i,j,ide,idt,cor1,cor2,tau,ui,vi) &
        !$omp private(vecr,vecr1,vecr2,nr,tg,A) &
        !$omp schedule(static)
        do i=1,mesh%nt
            vecr = 0d0;
            A = mesh%tr(i)%areag
            do j=1,3
                ide = mesh%tr(i)%ed(j);
                nr = mesh%edhx(ide)%tg
                tg = mesh%edhx(ide)%nr

                idt = mesh%ed(ide)%sh;

                cor1 = dsign(1d0,dot_product(mesh%ed(ide)%c%p-mesh%tr(i)%b%p, &
                    mesh%tr(idt(1))%b%p-mesh%tr(idt(2))%b%p));
                ! cor2 = dsign(1d0,dot_product(mesh%ed(ide)%c%p-mesh%tr(i)%b%p,nr));
                ! print*,dot_product(mesh%ed(ide)%c%p-mesh%tr(i)%b%p,nr)
                ! stop

                ! ui(1) = sum(u%p(idt(1))%v*mesh%tr(idt(1))%b%lonvec)
                ! ui(2) = sum(u%p(idt(1))%v*mesh%tr(idt(1))%b%latvec)

                ! vi(1) = sum(u%p(idt(2))%v*mesh%tr(idt(2))%b%lonvec)
                ! vi(2) = sum(u%p(idt(2))%v*mesh%tr(idt(2))%b%latvec)

                vecr1 = sum(u%p(idt(1))%v*nr)*nr + sum(u%p(idt(1))%v*tg)*tg
                vecr2 = sum(u%p(idt(2))%v*nr)*nr + sum(u%p(idt(2))%v*tg)*tg
                vecr = vecr + (vecr2-vecr1)*cor1 &
                    *mesh%ed(ide)%leng/mesh%edhx(ide)%leng

                ! vecr1 = vecr1 + (ui(1) - vi(1))*cor1 &
                !     *mesh%ed(ide)%leng/mesh%edhx(ide)%leng*cor2 &
                !     *proj_vec_sphere(nr,mesh%tr(i)%b%p)
                ! vecr2 = vecr2 + (ui(2) - vi(2))*cor1 &
                !     *mesh%ed(ide)%leng/mesh%edhx(ide)%leng*cor2 &
                !     *proj_vec_sphere(nr,mesh%tr(i)%b%p);
            enddo
            lap%p(i)%v = lap%p(i)%v + proj_vec_sphere(vecr,mesh%tr(i)%b%p)/(A*rad**2)
        enddo
        !$omp end parallel do
        do i =1,3
            lap%p(:)%v(i) = lap%p(:)%v(i)*tau;
        enddo

    end subroutine visc2x_tr2tr

    ! subroutine visc2x_tr2tr(u,lap,tau,mesh)
    !     implicit none
    !     type(grid_structure),     intent(in)   :: mesh
    !     type(vector_field_cart), intent(in)    :: u
    !     type(vector_field_cart), intent(inout) :: lap
    !     real*8,intent(in)                      :: tau
    !     integer :: i,j,ide,idt(2)
    !     real*8  :: cor,vecr(3),A

    !     do i=1,3
    !         lap%p(:)%v(i) = 0d0;
    !         lap%p(:)%v(i) = 0d0;
    !     enddo

    !     !$omp parallel do &
    !     !$omp default(none) &
    !     !$omp shared(mesh, u, lap) &
    !     !$omp private(i,j,ide,idt,cor,tau) &
    !     !$omp private(vecr,A) &
    !     !$omp schedule(static)
    !     do i=1,mesh%nt
    !         vecr = 0d0;
    !         A = mesh%tr(i)%areag
    !         do j=1,3
    !             ide = mesh%tr(i)%ed(j);

    !             idt = mesh%ed(ide)%sh;
    !             cor = dsign(1d0,dot_product(mesh%ed(ide)%c%p-mesh%tr(i)%b%p,mesh%tr(idt(1))%b%p-mesh%tr(idt(2))%b%p));

    !             vecr = vecr + (u%p(idt(1))%v - u%p(idt(2))%v)*cor*mesh%ed(ide)%leng/mesh%edhx(ide)%leng;
    !         enddo
    !         lap%p(i)%v = lap%p(i)%v + proj_vec_sphere(vecr,mesh%tr(i)%b%p)/(A*rad**2)
    !     enddo
    !     !$omp end parallel do
    !     do i =1,3
    !         lap%p(:)%v(i) = lap%p(:)%v(i)*tau;
    !     enddo

    ! end subroutine visc2x_tr2tr

    subroutine inicond(j,u,h,htr,bath,div,gradb,zeta,ke,f,fv,momfluxc,momfluxv,mesh,resn)
        implicit none
        type(grid_structure), intent(in) :: mesh
        type(scalar_field), intent(inout):: h,htr,div,zeta,f,fv,ke,bath
        type(vector_field_cart), intent(inout) :: u,gradb,momfluxc,momfluxv

        integer,intent(in) :: j
        integer            :: i,m,n
        real*8             :: lon,lat,b0, &
            gradlon,gradlat,alpha,u0,ulon,ulat,C,dfdx,dfdy, &
            dudx,dvdx,dudy,dvdy, &
            xl0(3),yl0(3),momlon,momlat,p(3), &
            lonvec(3),latvec(3), &
            hf,uf,ff,ufv(3)
        real*8             :: h0,R0,lonc,latc,r
        real*8  :: nygg,dygg,l1,beta,cc1,cc2,cc3,cc4,lat0,lat1, &
            en,u00,umen,utmp,den,vectmp(3),lat2,l2,e1,e2,clat
        integer :: jy
        real*8,allocatable :: hgg(:)
        character(len=100) ::resn

        h%fexact = 0d0;
        zeta%fexact  = 0d0;
        div%fexact   = 0d0;
        do i=1,3
            gradb%pexact%v(i) = 0d0;
            u%pexact%v(i) = 0d0;
        enddo

        bath%f = 0d0
        f%f  = 0d0;


        if (j==2) then
            write(resn,*) '../resultB/TC2/'

            u0 = 2d1;
            h0 = 5960d0;
            b0 = 2000d0;
            R0 = pi/9d0;
            lonc = -pi*.5d0;
            latc = pi/6d0;
            C = 1d0/g*(rad*omega*u0+u0**2*.5d0);

            do i=1,mesh%nv
                lon = mesh%v(i)%lon;
                lat = mesh%v(i)%lat;
                fv%f(i) = 2*omega*dsin(lat);

                r   = dsqrt(min(R0**2,(lon-lonc)**2+(lat-latc)**2));
                bath%f(i)   = b0*(1d0-r/R0);
                h%f(i) = h0-C*dsin(lat)**2 - bath%f(i);
            enddo
            do i =1,mesh%nt
                ! lon = mesh%tr(i)%b%lon;
                ! lat = mesh%tr(i)%b%lat;
                ! xl0 = mesh%tr(i)%b%lonvec;
                ! yl0 = mesh%tr(i)%b%latvec;

                p = mesh%tr(i)%b%p
                lat = datan2(p(3),dsqrt(p(1)**2+p(2)**2))
                lon = datan2(p(2),p(1))

                call convert_vec_sph2cart(1d0,0d0,p,xl0)
                yl0 = 0d0


                f%f(i) = 2*omega*dsin(lat);

                ulon =  u0*dcos(lat);
                ulat =  0d0;

                u%p(i)%v = ulon*xl0 + ulat*yl0;
            enddo

        elseif (j==3) then
            write(resn,*) '../resultB/TC3/'

            m=1; n=1;

            do i=1,mesh%nt
                lon = mesh%tr(i)%b%lon;
                lat = mesh%tr(i)%b%lat;

                ulon = -m*dcos(n*lat)**4/dcos(lat)*dsin(lon)*dsin(m*lon);
                ulat = -4*n*dcos(n*lat)**3*dsin(n*lat)*dsin(lon)*dsin(m*lon);

                dudy = -m*(-4*n*dcos(n*lat)**3*dsin(n*lat)/dcos(lat) + &
                            dcos(n*lat)**4/dcos(lat)**2*dsin(lat))*dsin(lon)*dsin(m*lon);
                dvdx = -4*n*dcos(n*lat)**3*dsin(n*lat)*(dcos(lon)*dsin(m*lon) + m*dsin(lon)*dcos(m*lon));


                dfdx = -m*dsin(m*lon)*dcos(n*lat)**4;
                dfdy = dcos(m*lon)*(-4*n*dcos(n*lat)**3*dsin(n*lat));

                gradlon = 1d0/(rad*dcos(lat))*dfdx;
                gradlat = 1d0/rad*dfdy;

                gradb%pexact(i)%v = (gradlon*mesh%tr(i)%b%lonvec + gradlat*mesh%tr(i)%b%latvec)*g;

                u%p(i)%v = ulon*mesh%tr(i)%b%lonvec + ulat*mesh%tr(i)%b%latvec;
            enddo

            do i=1,mesh%nv
                lon = mesh%v(i)%lon;
                lat = mesh%v(i)%lat;

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

                zeta%fexact(i) = -(dcos(lat)*dudy-dsin(lat)*ulon-dvdx)/(rad*dcos(lat))
        
                momlon = ulon/(rad*dcos(lat))*dudx + ulat/rad*dudy;
                momlat = ulon/(rad*dcos(lat))*dvdx + ulat/rad*dvdy;
                momfluxv%pexact(i)%v = momlon*mesh%v(i)%c%lonvec + momlat*mesh%v(i)%c%latvec;

                
                div%fexact(i) = 1/(rad*dcos(lat))*((dvdy*dcos(lat)-ulat*dsin(lat))+dudx);

                h%f(i) = dcos(m*lon)*dcos(n*lat)**4;
            enddo

        elseif (j==4) then
            write(resn,*) '../resultB/TC4/';
            alpha = 0d0*pi/2d0;
            u0 = 2d0*pi*rad/(86400d0*12);
            b0 = 2.94d4/g;

            C = 1d0/g*(rad*omega*u0+u0**2*.5d0);
            do i=1,mesh%nt
                lon = mesh%tr(i)%b%lon;
                lat = mesh%tr(i)%b%lat;

                xl0 = mesh%tr(i)%b%lonvec;
                yl0 = mesh%tr(i)%b%latvec;

                p = mesh%tr(i)%b%p
                lat = datan2(p(3),dsqrt(p(1)**2+p(2)**2))
                lon = datan2(p(2),p(1))

                call convert_vec_sph2cart(1d0,0d0,p,xl0)
                yl0 = 0d0

                f%f(i)  = 2d0*omega*(-dcos(lon)*dcos(lat)*dsin(alpha) + dsin(lat)*dcos(alpha));

                ulon =  u0*(dcos(lat)*dcos(alpha)+dcos(lon)*dsin(lat)*dsin(alpha));
                ulat =  -u0*dsin(lon)*dsin(alpha);

                u%p(i)%v = ulon*xl0 + ulat*yl0;


                dfdx = -C*2d0*(-dcos(lon)*dcos(lat)*dsin(alpha) &
                    + dsin(lat)*dcos(alpha))*(dsin(lon)*dcos(lat)*dsin(alpha));

                dfdy = -C*2d0*(-dcos(lon)*dcos(lat)*dsin(alpha) &
                    + dsin(lat)*dcos(alpha))*(dcos(lon)*dsin(lat)*dsin(alpha) + &
                    dcos(lat)*dcos(alpha));

                gradlon = 1d0/(rad*dcos(lat))*dfdx;
                gradlat = 1d0/rad*dfdy;

                gradb%pexact(i)%v = g*(gradlon*xl0 + gradlat*yl0);

                dudx   = u0*(-dsin(lon)*dsin(lat)*dsin(alpha));
                dudy   = u0*(-dsin(lat)*dcos(alpha)+dcos(lon)*dcos(lat)*dsin(alpha));
                dvdx   = -u0*dcos(lon)*dsin(alpha);
                dvdy   = 0d0;
                momlon = ulon*dudx/(rad*dcos(lat))+ulat*dudy/rad;!-ulon*ulat*dtan(lat)/rad;
                momlat = ulon*dvdx/(rad*dcos(lat))+ulat*dvdy/rad;!+ulon**2*dtan(lat)/rad;

                momfluxc%pexact(i)%v(:) = momlon*xl0+momlat*yl0;
            enddo

            do i=1,mesh%nv
                lat = mesh%v(i)%lat;
                lon = mesh%v(i)%lon;

                xl0 = mesh%v(i)%c%lonvec;
                yl0 = mesh%v(i)%c%latvec;

                fv%f(i)  = 2d0*omega*(-dcos(lon)*dcos(lat)*dsin(alpha) &
                    + dsin(lat)*dcos(alpha));

                bath%f(i) = b0 - C*(-dcos(lon)*dcos(lat)*dsin(alpha) + dsin(lat)*dcos(alpha))**2;

                ulon =  u0*(dcos(lat)*dcos(alpha)+dcos(lon)*dsin(lat)*dsin(alpha));
                ulat =  -u0*dsin(lon)*dsin(alpha);

                ke%fexact(i) = (ulon**2+ulat**2)*.5d0;

                dudx   = u0*(-dsin(lon)*dsin(lat)*dsin(alpha));
                dudy   = u0*(-dsin(lat)*dcos(alpha)+dcos(lon)*dcos(lat)*dsin(alpha));
                dvdx   = -u0*dcos(lon)*dsin(alpha);
                dvdy   = 0d0;
                momlon = ulon*dudx/(rad*dcos(lat))+ulat*dudy/rad-ulon*ulat*dtan(lat)/rad;
                momlat = ulon*dvdx/(rad*dcos(lat))+ulat*dvdy/rad+ulon**2*dtan(lat)/rad;

                zeta%fexact(i) = -1d0/(rad*dcos(lat))*(dcos(lat)*dudy-dsin(lat)*ulon-dvdx);

                div%fexact(i) = 1/(rad*dcos(lat))*((dvdy*dcos(lat)-ulat*dsin(lat))+dudx);

                momfluxv%pexact(i)%v(:) = momlon*xl0+momlat*yl0;

            enddo
            h%f = Hm

            do i=1,3
                u%pexact%v(i) = u%p%v(i);
            enddo
        elseif (j==8) then
            write(resn,*) '../resultB/TC8/';
            alpha = 0d0*pi/2d0;
            u0 = 2d0*pi*rad/(86400d0*12);
            b0 = 2.94d4/g;

            C = 1d0/g*(rad*omega*u0+u0**2*.5d0);
            do i=1,mesh%nt
                lon = mesh%tr(i)%b%lon;
                lat = mesh%tr(i)%b%lat;

                xl0 = mesh%tr(i)%b%lonvec;
                yl0 = mesh%tr(i)%b%latvec;

                p = mesh%tr(i)%b%p
                lat = datan2(p(3),dsqrt(p(1)**2+p(2)**2))
                lon = datan2(p(2),p(1))

                call convert_vec_sph2cart(1d0,0d0,p,xl0)
                !yl0 = 0d0

                f%f(i)  = 2d0*omega*(-dcos(lon)*dcos(lat)*dsin(alpha) + dsin(lat)*dcos(alpha));

                ulon =  u0*(dcos(lat)*dcos(alpha)+dcos(lon)*dsin(lat)*dsin(alpha));
                ulat =  -u0*dsin(lon)*dsin(alpha);

                u%p(i)%v = ulon*xl0 + ulat*yl0;


                dfdx = -C*2d0*(-dcos(lon)*dcos(lat)*dsin(alpha) &
                    + dsin(lat)*dcos(alpha))*(dsin(lon)*dcos(lat)*dsin(alpha));

                dfdy = -C*2d0*(-dcos(lon)*dcos(lat)*dsin(alpha) &
                    + dsin(lat)*dcos(alpha))*(dcos(lon)*dsin(lat)*dsin(alpha) + &
                    dcos(lat)*dcos(alpha));

                gradlon = 1d0/(rad*dcos(lat))*dfdx;
                gradlat = 1d0/rad*dfdy;

                gradb%pexact(i)%v = g*(gradlon*xl0 + gradlat*yl0);

                dudx   = u0*(-dsin(lon)*dsin(lat)*dsin(alpha));
                dudy   = u0*(-dsin(lat)*dcos(alpha)+dcos(lon)*dcos(lat)*dsin(alpha));
                dvdx   = -u0*dcos(lon)*dsin(alpha);
                dvdy   = 0d0;
                momlon = ulon*dudx/(rad*dcos(lat))+ulat*dudy/rad;!-ulon*ulat*dtan(lat)/rad;
                momlat = ulon*dvdx/(rad*dcos(lat))+ulat*dvdy/rad;!+ulon**2*dtan(lat)/rad;

                momfluxc%pexact(i)%v(:) = momlon*xl0+momlat*yl0;
            enddo

            do i=1,mesh%nv
                lat = mesh%v(i)%lat;
                lon = mesh%v(i)%lon;

                xl0 = mesh%v(i)%c%lonvec;
                yl0 = mesh%v(i)%c%latvec;

                fv%f(i)  = 2d0*omega*(-dcos(lon)*dcos(lat)*dsin(alpha) &
                    + dsin(lat)*dcos(alpha));

                h%f(i) = b0 - C*(-dcos(lon)*dcos(lat)*dsin(alpha) + dsin(lat)*dcos(alpha))**2;

                ulon =  u0*(dcos(lat)*dcos(alpha)+dcos(lon)*dsin(lat)*dsin(alpha));
                ulat =  -u0*dsin(lon)*dsin(alpha);

                ke%fexact(i) = (ulon**2+ulat**2)*.5d0;

                dudx   = u0*(-dsin(lon)*dsin(lat)*dsin(alpha));
                dudy   = u0*(-dsin(lat)*dcos(alpha)+dcos(lon)*dcos(lat)*dsin(alpha));
                dvdx   = -u0*dcos(lon)*dsin(alpha);
                dvdy   = 0d0;
                momlon = ulon*dudx/(rad*dcos(lat))+ulat*dudy/rad-ulon*ulat*dtan(lat)/rad;
                momlat = ulon*dvdx/(rad*dcos(lat))+ulat*dvdy/rad+ulon**2*dtan(lat)/rad;

                zeta%fexact(i) = -1d0/(rad*dcos(lat))*(dcos(lat)*dudy-dsin(lat)*ulon-dvdx);

                div%fexact(i) = 1/(rad*dcos(lat))*((dvdy*dcos(lat)-ulat*dsin(lat))+dudx);

                momfluxv%pexact(i)%v(:) = momlon*xl0+momlat*yl0;

            enddo

            h%fexact = h%f;
            do i=1,3
                u%pexact%v(i) = u%p%v(i);
            enddo
        elseif (j==9) then
            write(resn,*) '../resultB/TC9/';

            call TC9mfield(hgg,mesh)

            nygg = 4*FLOOR(SQRT(REAL(mesh%nv)))
            dygg = pi/nygg

            do i =1,mesh%nv
                p = mesh%v(i)%p

                lat = datan2(p(3),dsqrt(p(1)**2+p(2)**2))
                lon = datan2(p(2),p(1))

                fv%f(i) = 2d0*omega*dsin(lat)

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
            do i=1,mesh%nt
                p = mesh%tr(i)%b%p
                lat = datan2(p(3),dsqrt(p(1)**2+p(2)**2))
                lon = datan2(p(2),p(1))

                f%f(i) = 2d0*omega*dsin(lat)
                
                utmp =0

                den = (lat - lat0)*(lat - lat1)
                if (den < 0.0D0) then
                    utmp = umen*dexp(1.0D0/den)
                end if
                call convert_vec_sph2cart(utmp, 0d0, p, vectmp)
                u%p(i)%v = vectmp

                print*,den,umen,dexp(1.0D0/den),utmp
                stop
            enddo
            stop

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

                e1 = EXP(-(l1/alpha)**2)
                e2 = EXP(-((lat2 - l2)/beta)**2)
                h%f(i) = h%f(i)+h0*clat*e1*e2
        

            enddo
    
    
        elseif(j==10)then
            write(resn,*) '../resultB/TC10/';
    
            do i=1,mesh%nt
                p = mesh%tr(i)%b%p;
                p = rotation(p)

                lat = datan2(p(3),dsqrt(p(1)**2+p(2)**2));
                lon = datan2(p(2),p(1));
                lonvec = mesh%tr(i)%b%lonvec
                latvec = mesh%tr(i)%b%latvec

                call TC10(hf,uf,ff,dfdy,lat)
                ufv = uf*lonvec

                call convert_vec_sph2cart(uf, 0d0, p, ufv)
                ufv = Irotation(ufv)
                u%p(i)%v = ufv;
                f%f(i)  = ff;

            enddo

            do i=1,mesh%nv
                p = mesh%v(i)%p
                p = rotation(p)

                lat = datan2(p(3),dsqrt(p(1)**2+p(2)**2))
                lon = datan2(p(2),p(1))

                call TC10(hf,uf,ff,dfdy,lat)

                h%f(i)   = hf;
                fv%f(i)  = ff;
            enddo

        endif

        call scalar_ve2tr(h,htr,mesh);
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


    subroutine calcop(u,h,htr,uh,bath,uhq_perp,zeta,div,gradb,ke,&
        f,fv,momflux,tau,lap,bih,bihm,mesh)
        implicit none
        type(grid_structure), intent(in) :: mesh
        type(scalar_field),intent(in)         :: h,f,fv,bath
        type(vector_field_cart),intent(inout) :: u,gradb,momflux,uh,uhq_perp,lap,bih
        type(scalar_field),intent(inout)      :: htr,zeta,div,ke,bihm
        real*8,intent(in)                     :: tau
        ! logical, intent(in) :: vectinv
        type(scalar_field)                    :: kee,he,phi,phiv,lapm
        type(vector_field_cart)               :: momfluxv,momfluxc,ue
        real*8                                :: taun
        integer                               :: i


        ! --------------------------------------------
        allocate(lapm%f(mesh%nv));
        allocate(phiv%f(mesh%nv),phi%f(mesh%ne))
        allocate(kee%f(mesh%ne),he%f(mesh%ne))
        allocate(ue%p(mesh%ne), &
            momfluxv%p(mesh%nv),momfluxc%p(mesh%nt));
        ! --------------------------------------------

        zeta%f     = 0d0; ke%f       = 0d0;
        bihm%f     = 0d0; lapm%f     = 0d0;
        div%f = 0d0;

        do i=1,3
            momflux%p%v(i)  = 0d0;
            momfluxv%p%v(i)  = 0d0;
            bih%p%v(i)      = 0d0;
            lap%p%v(i)      = 0d0;
            uhq_perp%p%v(i) = 0d0;
        enddo

        if (nlin) then
            call scalar_ve2tr(h,htr,mesh)
        else
            htr%f = 1d4
        endif
        ! call ke_tr2ve(u,ke,mesh);

        if (vectinv)  then
            do i=1,3
                uh%p(:)%v(i) = u%p(:)%v(i)*htr%f;
            enddo
            call div_tr2ve(uh,div,mesh);  

            call ke_tr2ve(u,ke,mesh);
            phiv%f = g*(bath%f + h%f) + ke%f;
            call scalar_ve2ed(phiv,phi,mesh)
            call grad_ed2tr(phi,gradb,mesh);  
  
            call vort_tr2ve(u,zeta,mesh)
            call vectinv_ve2tr(zeta,fv,u,uhq_perp,mesh)

            if (.true.) then
                call visc2x_tr2tr(u,lap,tau,mesh);
                taun = 1d0;
                call visc2x_tr2tr(lap,bih,taun,mesh);
            endif
        else
            do i =1,3
                ! u%p(:)%v(i) = uh%p(:)%v(i)/htr%f;
                u%p(:)%v(i) = mulreal(uh%p(:)%v(i),1/htr%f);
            enddo
            ! phiv%f = g*(bath%f + h%f);
            phiv%f = g*sumreal(bath%f, h%f);

            if (nlin) then
                call divscalcontrvolume(uh,u,momfluxv,div,mesh);
                call vector_ve2tr(momfluxv,momflux,mesh);
            endif
            call div_tr2ve(uh,div,mesh);  

            ! call spherical_correction(u,momflux,mesh)
            call coriolisterm(f,uh,uhq_perp,mesh);

            call scalar_ve2ed(phiv,phi,mesh)
            call grad_ed2tr(phi,gradb,mesh);

            do i =1,3
                ! gradb%p(:)%v(i) = gradb%p(:)%v(i)*htr%f;
                if (nlin) then
                    gradb%p(:)%v(i) = mulreal(gradb%p(:)%v(i),htr%f);
                else
                    gradb%p(:)%v(i) = gradb%p(:)%v(i)*1d4
                endif
            enddo

            if (abs(tau)>1d-14) then
                call visc2x_tr2tr(uh,lap,tau,mesh);
                taun = 1d0;
                call visc2x_tr2tr(lap,bih,taun,mesh);
            endif
        endif

        ! print*,maxval(abs(bih%p%v(1))),maxval(abs(bih%p%v(2))), &
        ! maxval(abs(bih%p%v(3)))
        ! stop


    end subroutine calcop

    subroutine error_op(u,h,zeta,div,ke,gradb,f,momfluxc,momfluxv,filef,mesh)
        implicit none
        type(grid_structure), intent(in)       :: mesh
        type(vector_field_cart),intent(in)     :: u
        type(scalar_field),intent(in)          :: h
        type(vector_field_cart),intent(inout)  :: gradb,momfluxc,momfluxv
        type(scalar_field),intent(inout)       :: zeta,div,ke,f
        type(vector_field_cart)                :: utmp,momflux,uh,lap,bih,gradk,uhq_perp
        type(scalar_field)                     :: htr,he,bath,phi,phiv,bihm
        integer                                :: i
        character(len=100),intent(in)          :: filef
        character(len=100)                     :: filen
        real*8 :: tau,vecgradn(mesh%nt,3),vecgrada(mesh%nt,3),vec(mesh%nv,3),veca(mesh%nv,3), &
            vecm(mesh%nt,3)

        tau = 1d0

        allocate(he%f(mesh%ne),bath%f(mesh%nv));
        allocate(uhq_perp%p(mesh%nt), &
            gradk%p(mesh%nt),    &
            lap%p(mesh%nt),      &
            bih%p(mesh%nt),      &
            bihm%f(mesh%nv),     &
            htr%f(mesh%nt),     &
            momflux%p(mesh%nt), &
            uh%p(mesh%nt),utmp%p(mesh%nt));
        allocate(phi%f(mesh%ne),phiv%f(mesh%nv));


        bath%f = 0d0;
        do i=1,3
            momfluxc%p(:)%v(i) = 0d0;
            utmp%p(:)%v(i) = u%p(:)%v(i);
        enddo

        ! call calcop(utmp,h,htr,uh,bath,uhq_perp,zeta,div,gradb,ke,f,fv,momflux,tau,lap,bih,bihm,mesh)

        phiv%f = g*(bath%f + h%f);
        call coriolisterm(f,u,uhq_perp,mesh);
        call scalar_ve2ed(phiv,phi,mesh);
        call grad_ed2tr(phi,gradb,mesh);
        call scalar_ve2tr(h,htr,mesh);
        do i=1,3
            uh%p(:)%v(i) = u%p(:)%v(i)*htr%f;
        enddo

        call div_tr2ve(uh,div,mesh);
        call vort_tr2ve(u,zeta,mesh);
        call ke_tr2ve(u,ke,mesh);

        call divscalcontrvolume(u,u,momfluxv,div,mesh);
        call div_tr2ve(u,div,mesh);

        call vector_ve2tr(momfluxv,momfluxc,mesh);
        ! call spherical_correction(u,momfluxc,mesh)

        do i=1,3
            vecgradn(:,i) = gradb%p(:)%v(i);
            vecgrada(:,i) = gradb%pexact(:)%v(i);
            vec(:,i) = momfluxv%p(:)%v(i);
            veca(:,i) = momfluxv%pexact(:)%v(i);
            vecm(:,i) = momfluxc%p(:)%v(i)+uhq_perp%p(:)%v(i)+gradb%p(:)%v(i)
        enddo
      
        write(*,*) 'Div:   ',maxval(abs(div%f)),maxval(abs(div%fexact)), &
            maxval(abs(div%f - div%fexact));
        write(*,*) 'Zeta:  ',maxval(abs(zeta%f)),maxval(abs(zeta%fexact)), &
            maxval(abs(zeta%f - zeta%fexact));
        write(*,*) 'Grad:  ',maxval(normm(vecgradn)),maxval(normm(vecgrada)), &
            maxval(normm(vecgradn - vecgrada));
        write(*,*) 'ketr:  ',maxval(abs(ke%f)),maxval(abs(ke%fexact)), &
            maxval(abs(ke%f - ke%fexact));
        write(*,*) 'momv:  ',maxval(normm(vec)),maxval(normm(veca)), &
            maxval(normm(vec - veca));
        write(*,*) 'moment:',maxval(normm(vecm)),0d0,maxval(normm(vecm));
        ! stop


        write(filen,'(3A)') trim(filef),'div.dat';
        open(unit = 100, file = filen(2:100));

        write(filen,'(3A)') trim(filef),'zeta.dat';
        open(unit = 101, file = filen(2:100));

        write(filen,'(3A)') trim(filef),'gradb.dat';
        open(unit = 102, file = filen(2:100));

        write(filen,'(3A)') trim(filef),'ke.dat';
        open(unit = 103, file = filen(2:100));

        write(filen,'(3A)') trim(filef),'momflux.dat';
        open(unit = 104, file = filen(2:100));

        write(filen,'(3A)') trim(filef),'momentum.dat';
        open(unit = 105, file = filen(2:100));

        do i=1,mesh%nv
            write(100,*) div%f(i),div%fexact(i)
            write(101,*) zeta%f(i),zeta%fexact(i)
            write(103,*) ke%f(i),ke%fexact(i)
            write(104,*) momfluxv%p(i)%v(:),momfluxv%pexact(i)%v(:)
        enddo

        do i=1,mesh%nt
            write(102,*) gradb%p(i)%v(:),gradb%pexact(i)%v(:)
            write(105,*) momfluxc%p(i)%v(:)+uhq_perp%p(i)%v(:)+gradb%p(i)%v(:);
        enddo

        do i=100,105
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
        enddo

        do i=1,mesh%nt
          write(111,*) u%p(i)%v;
        enddo
        call FLUSH(110)
        call FLUSH(111)

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

    subroutine savetimesteps(u,h,Ek,zeta,div,time,j)
        type(scalar_field),intent(in) :: h,div,Ek,zeta
        type(vector_field_cart),intent(in) :: u
        real*8,intent(in) :: time
        integer,intent(in) :: j
        real*8 :: hsca,divsca,Eksca,usca


        hsca = maxval(h%f-h%fexact)
        usca = max(maxval(u%p(1)%v),maxval(u%p(2)%v),maxval(u%p(3)%v))
        divsca = maxval(div%f)
        Eksca = sum(Ek%f)
        zetasca = sum(zeta%f)

        write(j,*) time,hsca,usca,divsca,Eksca,zetasca


    end subroutine savetimesteps

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
        enddo

        do i=1,mesh%nt
            read(111,*) u%p(i)%v;
        enddo
        do i=110,111
            close(i);
        enddo
    endsubroutine loadini
end module swmb_operators
