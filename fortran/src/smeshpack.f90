module smeshpack

  use basic_funcs
  use constants
  use datastruct

  implicit none

contains 

    ! subroutine savevars(mesh,fileg)
    !     type(grid_structure),intent(inout) :: mesh
    !     integer :: i
    !     character(len=100) ::

    !     do i =1,2
    !     enddo

    ! endsubroutine savevars

    subroutine calcTopoVar(mesh)
        type(grid_structure),intent(inout) :: mesh

        call calc_tricenter(mesh);
        call edge_mid(mesh);
        call calc_sph_pts(mesh);
        call calc_sph_vecs(mesh)
        call calc_unit_vec(mesh);
        call calc_tri_area(mesh);
        call calc_hex_area(mesh);
        call edgelength(mesh);
        call chgVecDir(mesh);
        call tri_vec_sign(mesh);
        call hex_vec_sign(mesh);
        call calc_unit_vecB(mesh);
        call edgeB_length(mesh);
        call unit_vecB_sign(mesh);
        call tria_weight(mesh);
        call hex_weight(mesh);
        call trisk_hx(mesh);

        ! call perot(mesh);
    end subroutine calcTopoVar

    subroutine calc_sph_pts(mesh)
        implicit none
        type(grid_structure),intent(inout) :: mesh
        real*8 :: pq(mesh%nv,3),pe(mesh%ne,3),pex(mesh%ne,3),pp(mesh%nt,3),ppc(mesh%nt,3), &
            sphpq(mesh%nv,2),sphpe(mesh%ne,2),sphpex(mesh%ne,2), &
            sphpp(mesh%nt,2),sphppc(mesh%nt,2)
        integer :: i
        
        do i =1,3
            pq(:,i)  = mesh%v(:)%p(i);
            pe(:,i)  = mesh%ed(:)%c%p(i);
            pex(:,i) = mesh%edhx(:)%c%p(i);
            pp(:,i)  = mesh%tr(:)%c%p(i);
            ppc(:,i)  = mesh%tr(:)%b%p(i);
        enddo

        sphpq  = spherical_coord(pq);
        sphpe  = spherical_coord(pe);
        sphpex = spherical_coord(pex);
        sphpp  = spherical_coord(pp);
        sphppc = spherical_coord(ppc);

        mesh%v(:)%lon = sphpq(:,1);
        mesh%v(:)%lat = sphpq(:,2);

        mesh%ed(:)%c%lon = sphpe(:,1);
        mesh%ed(:)%c%lat = sphpe(:,2);

        mesh%edhx(:)%c%lon = sphpex(:,1);
        mesh%edhx(:)%c%lat = sphpex(:,2);

        mesh%tr(:)%c%lon = sphpp(:,1);
        mesh%tr(:)%c%lat = sphpp(:,2);

        mesh%tr(:)%b%lon = sphppc(:,1);
        mesh%tr(:)%b%lat = sphppc(:,2);

    end subroutine calc_sph_pts

    subroutine calc_sph_vecs(mesh)
        implicit none
        type(grid_structure),intent(inout) :: mesh
        real*8 :: pq(mesh%nv,3),pe(mesh%ne,3),pex(mesh%ne,3),pp(mesh%nt,3),ppc(mesh%nt,3), &
                sphvecq(mesh%nv,2,3),sphvece(mesh%ne,2,3),sphvecex(mesh%ne,2,3),sphvecp(mesh%nt,2,3), &
                sphvecpc(mesh%nt,2,3)
        integer :: i

        do i=1,3
            pq(:,i)  = mesh%v(:)%p(i);
            pe(:,i)  = mesh%ed(:)%c%p(i);
            pex(:,i) = mesh%edhx(:)%c%p(i);
            pp(:,i)  = mesh%tr(:)%c%p(i);
            ppc(:,i)  = mesh%tr(:)%b%p(i);
        enddo

        sphvecq  = vec2d_sph2cart(pq);
        sphvece  = vec2d_sph2cart(pe);
        sphvecex = vec2d_sph2cart(pex);
        sphvecp  = vec2d_sph2cart(pp);
        sphvecpc = vec2d_sph2cart(ppc);

        
        do i=1,3
            mesh%v(:)%c%lonvec(i)    = sphvecq(:,1,i);
            mesh%v(:)%c%latvec(i)    = sphvecq(:,2,i);
            mesh%ed(:)%c%lonvec(i)   = sphvece(:,1,i);
            mesh%ed(:)%c%latvec(i)   = sphvece(:,2,i);
            mesh%edhx(:)%c%lonvec(i) = sphvecex(:,1,i);
            mesh%edhx(:)%c%latvec(i) = sphvecex(:,2,i);
            mesh%tr(:)%c%lonvec(i) = sphvecp(:,1,i);
            mesh%tr(:)%c%latvec(i) = sphvecp(:,2,i);
            mesh%tr(:)%b%lonvec(i) = sphvecpc(:,1,i);
            mesh%tr(:)%b%latvec(i) = sphvecpc(:,2,i);
        enddo

    end subroutine calc_sph_vecs

    subroutine calc_tricenter(mesh)
        type(grid_structure),intent(inout) :: mesh
        real*8  :: pq(mesh%nv,3),ppc(mesh%nt,3),pp(mesh%nt,3)
        integer :: i,tri(mesh%nt,3)

        do i =1,3
            pq(:,i) = mesh%v(:)%p(i);
            tri(:,i) = mesh%tr(:)%v(i);
        enddo
        ppc = trian_centroid(pq,tri);
        pp  = trian_circum(pq,tri);

        do i =1,3
            mesh%tr(:)%c%p(i) = pp(:,i)
            mesh%tr(:)%b%p(i) = ppc(:,i)
        enddo
    endsubroutine calc_tricenter

    subroutine edge_mid(mesh)
        TYPE(grid_structure) :: mesh
        real*8               :: pq(mesh%nv,3),pp(mesh%nt,3), &
                                pe(mesh%ne,3),pex(mesh%ne,3)
        integer              :: i,sh(mesh%ne,2),edv(mesh%ne,2)

        do i =1,3
            pq(:,i) = mesh%v(:)%p(i);
            pp(:,i) = mesh%tr(:)%c%p(i);
        enddo

        do i =1,2
            sh(:,i) = mesh%ed(:)%sh(i);
            edv(:,i) = mesh%ed(:)%v(i);
        enddo

        pe = edgemidpoint(pq,edv);
        pex  = edgemidpoint(pp,sh);

        do i =1,3
            mesh%edhx(:)%c%p(i) = pex(:,i)
            mesh%ed(:)%c%p(i)   = pe(:,i)
        enddo

    endsubroutine edge_mid

    subroutine calc_unit_vec(mesh)
        type(grid_structure),intent(inout) :: mesh
        real*8 :: pq(mesh%nv,3),pe(mesh%ne,3),pex(mesh%ne,3), pp(mesh%nt,3), &
                tg(mesh%ne,3),nr(mesh%ne,3), &
                tgx(mesh%ne,3),nrx(mesh%ne,3), &
                cor(mesh%ne),corhx(mesh%ne)
        integer :: i,sh(mesh%ne,2),edv(mesh%ne,2)

        do i=1,2
            sh(:,i) = mesh%ed(:)%sh(i);
            edv(:,i) = mesh%ed(:)%v(i);
        enddo

        do i=1,3
            pq(:,i) = mesh%v(:)%p(i);
            pp(:,i) = mesh%tr(:)%c%p(i);
            pe(:,i) = mesh%ed(:)%c%p(i);
            pex(:,i) = mesh%edhx(:)%c%p(i);
        enddo

        tg = edge_tg(pq,edv)
        nr = edge_nr(pe,tg)

        tgx = edge_tg(pp,sh)
        nrx = edge_nr(pex,tgx)

        do i=1,3
            mesh%ed(:)%nr(i) = nr(:,i);
            mesh%ed(:)%tg(i) = tg(:,i);
            mesh%edhx(:)%nr(i) = nrx(:,i);
            mesh%edhx(:)%tg(i) = tgx(:,i);
        enddo

        cor = dsign(1d0,mesh%ed(:)%tg(1)*mesh%edhx(:)%nr(1) &
            +mesh%ed(:)%tg(2)*mesh%edhx(:)%nr(2)+mesh%ed(:)%tg(3)*mesh%edhx(:)%nr(3))
        corhx = dsign(1d0,mesh%ed(:)%nr(1)*mesh%edhx(:)%tg(1) &
            +mesh%ed(:)%nr(2)*mesh%edhx(:)%tg(2)+mesh%ed(:)%nr(3)*mesh%edhx(:)%tg(3))

        do i=1,3
            mesh%edhx(:)%nr(i) = cor*mesh%edhx(:)%nr(i);
            mesh%edhx(:)%tg(i) = corhx*mesh%edhx(:)%tg(i);
        enddo
    endsubroutine calc_unit_vec

    subroutine edgelength(mesh)
        implicit none
        type(grid_structure),intent(inout) :: mesh
        real*8 :: pq1(mesh%ne,3),pq2(mesh%ne,3), &
                pp1(mesh%ne,3),pp2(mesh%ne,3)
        integer :: i,edv(mesh%ne,2),sh(mesh%ne,2)

        do i =1,2
            edv(:,i) = mesh%ed(:)%v(i)
            sh(:,i)  = mesh%ed(:)%sh(i)
        enddo
        do i =1,3
            pq1(:,i) = mesh%v(edv(:,1))%p(i);
            pq2(:,i) = mesh%v(edv(:,2))%p(i);

            pp1(:,i) = mesh%tr(sh(:,1))%c%p(i);
            pp2(:,i) = mesh%tr(sh(:,2))%c%p(i);
        enddo

        mesh%ed(:)%leng   = arclen2d(pq1,pq2);
        mesh%edhx(:)%leng = arclen2d(pp1,pp2);
    end subroutine edgelength

    subroutine chgVecDir(mesh)
        type(grid_structure),intent(inout) :: mesh
        real*8 :: pqe(mesh%ne,2,3),ppe(mesh%ne,2,3), &
            tg(mesh%ne,3),nr(mesh%ne,3), &
            tgx(mesh%ne,3),nrx(mesh%ne,3), &
            dvec(mesh%ne)
        integer :: i,j,sh(mesh%ne,2),edv(mesh%ne,2),vecsign(mesh%ne)

        do i=1,2
            sh(:,i) = mesh%ed(:)%sh(i);
            edv(:,i) = mesh%ed(:)%v(i);
        enddo

        do i=1,3
            do j=1,2
                pqe(:,j,i) = mesh%v(edv(:,j))%p(i);
                ppe(:,j,i) = mesh%tr(sh(:,j))%c%p(i);
            enddo

            nr(:,i) = mesh%ed(:)%nr(i);
            tg(:,i) = mesh%ed(:)%tg(i);
            nrx(:,i) = mesh%edhx(:)%nr(i);
            tgx(:,i) = mesh%edhx(:)%tg(i);
        enddo

        dvec  = 0d0;
        do i=1,3
            dvec  = dvec  + nr(:,i)*(ppe(:,2,i)-ppe(:,1,i));
        enddo

        vecsign = int(dsign(1d0,dvec));
        do i =1,3
            nr(:,i) = vecsign*nr(:,i)
            tg(:,i) = vecsign*tg(:,i)
        enddo

        dvec  = 0d0;
        do i=1,3
            dvec  = dvec  + nr(:,i)*tgx(:,i);
        enddo

        vecsign = int(dsign(1d0,dvec));
        do i =1,3
            nrx(:,i) = vecsign*nrx(:,i)
            tgx(:,i) = vecsign*tgx(:,i)
        enddo

        do i=1,3
            mesh%ed(:)%nr(i)   = nr(:,i);
            mesh%ed(:)%tg(i)   = tg(:,i);
            mesh%edhx(:)%nr(i) = nrx(:,i);
            mesh%edhx(:)%tg(i) = tgx(:,i);
        enddo

    endsubroutine chgVecDir

    subroutine tri_vec_sign(mesh)
        type(grid_structure),intent(inout) :: mesh
        real*8  :: pest(mesh%nt,3),pp(mesh%nt,3),vecref(mesh%nt,3),nr(mesh%nt,3), &
                    dvec(mesh%nt)
        integer :: i,j,etrii(mesh%nt),vecsign(mesh%nt)

        do i =1,3
            etrii = mesh%tr(:)%ed(i)
            do j =1,3
                pest(:,j) = mesh%ed(etrii)%c%p(j);
                pp(:,j) = mesh%tr(:)%c%p(j);
                nr(:,j) = mesh%ed(etrii)%nr(j);
            enddo
            vecref = pest-pp;
            dvec=0d0;
            do j =1,3
                dvec = dvec+vecref(:,j)*nr(:,j)
            enddo
            vecsign = int(dsign(1d0,dvec));
            mesh%tr(:)%cor(i) = vecsign
        enddo
    endsubroutine tri_vec_sign

    subroutine hex_vec_sign(mesh)
        type(grid_structure),intent(inout) :: mesh
        real*8  :: pq(10,3),pex(10,3),vecref(10,3),nr(10,3),dvec(10)
        integer :: i,j,nnb,ved(10),vecsign(10)

        do i =1,mesh%nv
            nnb = mesh%v(i)%nnb;
            ved(1:nnb) = mesh%v(i)%ed;
            do j=1,3
                pq(1:nnb,j) = mesh%v(i)%p(j);
                pex(1:nnb,j) = mesh%edhx(ved(1:nnb))%c%p(j);
                nr(1:nnb,j) = mesh%edhx(ved(1:nnb))%nr(j)
            enddo
            vecref = pex - pq;
            dvec    = sum(vecref*nr,dim=2);
            vecsign = int(dsign(1d0,dvec));
            mesh%hx(i)%cor=vecsign(1:nnb)
        enddo
    endsubroutine hex_vec_sign

    subroutine calc_tri_area(mesh)
        implicit none
        type(grid_structure),intent(inout) :: mesh
        real*8  :: pq(mesh%nv,3),pp(mesh%nt,3),pe(mesh%ne,3),A(mesh%nt,3)
        integer :: i,etri(mesh%nt,3),edv(mesh%ne,2)

        do i=1,3
            pq(:,i) = mesh%v(:)%p(i);
            pp(:,i) = mesh%tr(:)%c%p(i);
            pe(:,i) = mesh%ed(:)%c%p(i);
            etri(:,i) = mesh%tr(:)%ed(i);
        enddo

        do i=1,2
            edv(:,i) = mesh%ed(:)%v(i);
        enddo

        A = triangle2d_area(pe,pq,pp,etri,edv)

        mesh%tr(:)%areag = 0d0;
        do i=1,3
            mesh%tr(:)%ctrve_areag(i) = A(:,i);
            mesh%tr(:)%areag = mesh%tr(:)%areag + A(:,i);
        enddo

    end subroutine calc_tri_area

    subroutine calc_hex_area(mesh)
        implicit none
        type(grid_structure),intent(inout) :: mesh
        real*8                             :: pppe(10,2,3),pqs(10,3)
        integer                            :: i,j,k,nnb,eds(10),trs(10), &
            ted(10,2)

        do i=1,mesh%nv
            nnb = mesh%v(i)%nnb;
            eds(1:nnb) = mesh%v(i)%ed(:);
            trs(1:nnb) = mesh%v(i)%tr(:);
            do j=1,2
                ted(1:nnb,j) = mesh%ed(eds(1:nnb))%sh(j);
            enddo
            do j=1,3
                pqs(1:nnb,j) = mesh%v(i)%p(j);
                do k=1,2
                    pppe(1:nnb,k,j) = mesh%tr(ted(1:nnb,k))%c%p(j);
                enddo
            enddo
            mesh%hx(i)%areag   = hexag_cell_area(pqs(1:nnb,:),pppe(1:nnb,:,:))
            mesh%hx(i)%cdareag = sum(mesh%tr(trs(1:nnb))%areag)/3d0
        enddo
    end subroutine calc_hex_area

    !B-grid special script
    subroutine calc_unit_vecB(mesh)
        implicit none
        type(grid_structure),intent(inout) :: mesh
        real*8  :: pp(mesh%ne,3),pe(mesh%ne,3),nrs(mesh%ne,3),tgs(mesh%ne,3), &
                    pts(mesh%ne,3),normpts(mesh%ne),normtgs(mesh%ne)
        integer :: i,j,etrii(mesh%ne)
        do j=1,3
            pe(:,j) = mesh%ed(:)%c%p(j);
        enddo

        do i=1,2
            etrii = mesh%ed(:)%sh(i)
            do j=1,3
                pp(:,j) = mesh%tr(etrii)%b%p(j);
            enddo
            tgs = pp-pe;
            pts = pp+pe;

            normtgs = normm(tgs);
            normpts = normm(pts);
            do j=1,3
                tgs(:,j) = tgs(:,j)/normtgs;
                pts(:,j) = pts(:,j)/normpts;
            enddo
            nrs = edge_nr(pts,tgs);
            do j=1,3
                mesh%edhx(:)%nrs(i,j) = nrs(:,j);
                mesh%edhx(:)%tgs(i,j) = tgs(:,j);
            enddo
        enddo
        
    end subroutine calc_unit_vecB

    subroutine unit_vecB_sign(mesh)
        type(grid_structure) :: mesh
        real*8 :: pq(6,3),peu(6,3),vecr(6,3),nrs(6,3), &
                scav(6)
        integer :: i,j,k,nnb

        do i =1,mesh%nv
            nnb = mesh%v(i)%nnb;
            do j=1,3
                pq(:,j) = mesh%v(i)%p(j);
                peu(1:nnb,j) = mesh%edhx(mesh%v(i)%ed)%c%p(j);
            enddo
            vecr = peu-pq;
            do j=1,2
                do k=1,3
                    nrs(1:nnb,k) = mesh%edhx(mesh%v(i)%ed)%nrs(j,k)
                enddo
                scav = nrs(:,1)*vecr(:,1) + nrs(:,2)*vecr(:,2) + nrs(:,3)*vecr(:,3);
                mesh%hx(i)%cord(1:nnb,j) = int(dsign(1d0,scav(1:nnb)));
            enddo
        enddo
    endsubroutine unit_vecB_sign

    subroutine edgeB_length(mesh)
        implicit none
        type(grid_structure),intent(inout) :: mesh
        real*8  :: pp(mesh%ne,3),pe(mesh%ne,3)
        integer :: i,j,etrii(mesh%ne)

        do j=1,3
            pe(:,j) = mesh%ed(:)%c%p(j);
        enddo

        do i=1,2
            etrii = mesh%ed(:)%sh(i)
            do j=1,3
                pp(:,j) = mesh%tr(etrii)%b%p(j);
            enddo
            mesh%edhx(:)%lend(i) = arclen2d(pe,pp);
        enddo
    end subroutine edgeB_length

    !TRiSK
    subroutine tria_weight(mesh)
        implicit none
        type(grid_structure),intent(inout) :: mesh
        real*8  :: pq(mesh%nv,3),pp(mesh%nt,3),pe(mesh%ne,3),A(mesh%nt,3)
        integer :: i,etri(mesh%nt,3),tri(mesh%nt,3)

        do i=1,3
            pq(:,i) = mesh%v(:)%p(i);
            pp(:,i) = mesh%tr(:)%c%p(i);
            pe(:,i) = mesh%ed(:)%c%p(i);
            etri(:,i) = mesh%tr(:)%ed(i);
            tri(:,i) = mesh%tr(:)%v(i);
        enddo
        A = triangle2dpq_area(pe,pq,pp,tri,etri)
        do i=1,3
            mesh%tr(:)%trhx_areag(i) = A(:,i);
        enddo
    endsubroutine tria_weight

    subroutine hex_weight(mesh)
        implicit none
        type(grid_structure),intent(inout) :: mesh
        real*8 :: trtrhx_areag(mesh%nt,3)
        real*8 :: pq(mesh%nv,3),pp(mesh%nt,3),pe(mesh%ne,3)
        integer :: i,j,vtr(10),vvc(10),tri(mesh%nt,3),etri(mesh%nt,3),nnb

        do i=1,3
            pq(:,i) = mesh%v(:)%p(i);
            pp(:,i) = mesh%tr(:)%c%p(i);
            pe(:,i) = mesh%edhx(:)%c%p(i);
            etri(:,i) = mesh%tr(:)%ed(i);
            tri(:,i) = mesh%tr(:)%v(i);
        enddo
        trtrhx_areag = triangle2dpq_area(pe,pq,pp,tri,etri);

        !$omp parallel do &
        !$omp default(none) &
        !$omp shared(mesh,tri,trtrhx_areag) &
        !$omp private(i,j,vtr,vvc,nnb) &
        !$omp schedule(static)
        do i=1,mesh%nv
            nnb = mesh%v(i)%nnb;
            vtr(1:nnb) = mesh%v(i)%tr(:);
            vvc(1:nnb) = minloc(abs(tri(vtr(1:nnb),:)-i),2)

            do j=1,nnb
                mesh%hx(i)%trhx_areag(j) = trtrhx_areag(vtr(j),vvc(j));
            enddo
        enddo
        !$omp end parallel do

    endsubroutine hex_weight

    subroutine trisk_hx(mesh)
        implicit none
        type(grid_structure),intent(inout) :: mesh
        real*8 :: A(6,6)
        integer :: i,j,k,nnb

        do i =1,mesh%nv
            A = 0;
            nnb = mesh%v(i)%nnb;
            ! mesh%hx(i)%trhx_areag = mesh%hx(i)%areag/nnb
            do j =1,nnb
                do k =j+1,nnb
                    A(j,k) = sum(mesh%hx(i)%trhx_areag(j:k-1))/mesh%hx(i)%areag;
                    A(k,j) = 1-A(j,k);
                enddo
                A(j,j) = .5d0;
            enddo
            A = A-.5d0;
            mesh%hx(i)%trskwg = A(1:nnb,1:nnb);
        enddo
        
        ! do i=1,mesh%v(10000)%nnb
        !   write(*,*) mesh%hx(10000)%trskwg(i,:)
        ! enddo
        ! write(*,*) mesh%hx(10000)%trhx_areag,mesh%hx(10000)%areag;
        ! stop;

    endsubroutine trisk_hx

    !Harmonic Polynomial
    subroutine harmpol(P,Q,gP,gQ,deg)
        real*8,intent(inout) :: P(:,:,:),Q(:,:,:),gP(:,:,:,:),gQ(:,:,:,:)
        integer,intent(in) :: deg
        real*8 :: xP(size(P,2),size(P,3)),yP(size(P,2),size(P,3)), &
            xQ(size(P,2),size(P,3)),yQ(size(P,2),size(P,3))
        integer :: i

        gP = 0d0;
        gQ = 0d0;
        P = 0d0;
        Q = 0d0;

        gP(1,1,1,1) = 1;
        gQ(1,1,1,2) = 1;

        P(1,1,2) = 1;
        Q(1,2,1) = 1;

        do i=2,deg
            xP = eoshift(P(i-1,:,:),-1,dim=2);
            yP = eoshift(P(i-1,:,:),-1,dim=1);
            xQ = eoshift(Q(i-1,:,:),-1,dim=2);
            yQ = eoshift(Q(i-1,:,:),-1,dim=1);

            gP(i,:,:,1) = i*P(i-1,:,:);
            gP(i,:,:,2) = -i*Q(i-1,:,:);
            gQ(i,:,:,1) = i*Q(i-1,:,:);
            gQ(i,:,:,2) = i*P(i-1,:,:);
            P(i,:,:) = xP-yQ;
            Q(i,:,:) = yP+xQ;
        enddo
    endsubroutine harmpol

    subroutine computepol(pts,P,Q,gP,gQ,Pc,Qc,gPc,gQc)
        real*8,intent(in) :: pts(:,:,:),P(:,:,:),Q(:,:,:), &
            gP(:,:,:,:),gQ(:,:,:,:)
        real*8,intent(inout) :: Pc(size(P,1),size(pts,1),size(pts,2)), &
            Qc(size(P,1),size(pts,1),size(pts,2)), &
            gPc(size(P,1),size(pts,1),size(pts,2),size(pts,3)), &
            gQc(size(P,1),size(pts,1),size(pts,2),size(pts,3))
        real*8 :: polyc(size(pts,1),size(pts,2))
        integer :: i,j,k,l,ni,nj,deg

        deg = size(P,1)
        ni = deg+1;
        nj = ni;
        do i =1,deg
            do j=1,ni
                do k=1,nj
                    polyc = pts(:,:,1)**(k-1)*pts(:,:,2)**(j-1)
                    Pc(i,:,:) = Pc(i,:,:) + P(i,j,k)*polyc;
                    Qc(i,:,:) = Qc(i,:,:) + Q(i,j,k)*polyc;
                    do l=1,2
                        gPc(i,:,:,l) = gPc(i,:,:,l) + gP(i,j,k,l)*polyc;
                        gQc(i,:,:,l) = gQc(i,:,:,l) + gQ(i,j,k,l)*polyc;
                    enddo
                enddo
            enddo
        enddo
    endsubroutine computepol

    subroutine localcoord(p,pr,pl,localvec)
        real*8,intent(in)    :: localvec(:,:,:)
        real*8,intent(in) :: p(:,:,:),pr(:,:)
        real*8,intent(inout) :: pl(:,:,:)
        real*8               :: x(size(pl,1),size(pl,2)),y(size(pl,1),size(pl,2))
        real*8               :: dist(size(p,1),size(p,2)),vec(size(p,1),size(p,2),size(p,3)), &
        vecp(size(p,1),size(p,2),size(p,3))
        integer :: i,j,ni,nj

        ni = size(p,3);
        nj = size(p,2);
        dist = 0d0;
        do i=1,nj
            vec(:,i,:) = p(:,i,:)-pr;
        enddo


        do i=1,ni
            do j=1,nj
                dist(:,j) = dist(:,j) + vec(:,j,i)*pr(:,i)
            enddo
        enddo

        do j=1,nj
            do i=1,ni
                vecp(:,j,i) = p(:,j,i)-dist(:,j)*pr(:,i);
            enddo
        enddo

        x=0d0;
        y=0d0;
        do i =1,ni
            do j =1,nj
                x(:,j) = x(:,j) + (vecp(:,j,i)*localvec(:,i,1));
                y(:,j) = y(:,j) + (vecp(:,j,i)*localvec(:,i,2));
            enddo
        enddo


        pl(:,:,1) = x;
        pl(:,:,2) = y;

        ! write(*,*) localvec(6,:,1)
        ! write(*,*) localvec(6,:,2)
        ! stop;

    endsubroutine localcoord

    subroutine poligavgpol(Pq,Qq,Pp,Qp)
        real*8,intent(in)    :: Pq(:,:,:),Qq(:,:,:)
        real*8,intent(inout) :: Pp(:,:),Qp(:,:)


        Pp(:,:) = Pq(:,:,1) + Pq(:,:,2) + Pq(:,:,3);
        Qp(:,:) = Qq(:,:,1) + Qq(:,:,2) + Qq(:,:,3);
        Pp = Pp/3;
        Qp = Qp/3;
    endsubroutine poligavgpol

    subroutine faceavgpoly(Pq,Qq,Pf,Qf)
        real*8,intent(in)    :: Pq(:,:,:),Qq(:,:,:)
        real*8,intent(inout) :: Pf(:,:,:),Qf(:,:,:)

        Pf(:,:,1) = Pq(:,:,1) + Pq(:,:,2);
        Pf(:,:,2) = Pq(:,:,2) + Pq(:,:,3);
        Pf(:,:,3) = Pq(:,:,3) + Pq(:,:,1);

        Qf(:,:,1) = Qq(:,:,1) + Qq(:,:,2);
        Qf(:,:,2) = Qq(:,:,2) + Qq(:,:,3);
        Qf(:,:,3) = Qq(:,:,3) + Qq(:,:,1);

        Pf = Pf/2; Qf = Qf/2;
    endsubroutine faceavgpoly

    subroutine calcVR(Pwa,Qwa,Pfp,Qfp,x2,gP2p,gQ2p,VR)
        real*8,intent(in)    :: Pfp(:,:,:),Qfp(:,:,:),Pwa(:,:),Qwa(:,:),x2(:,:)
        real*8,intent(in)    :: gP2p(:,:),gQ2p(:,:)
        real*8,intent(inout) :: VR(:,:,:)
        real*8               :: V(size(VR,1),size(VR,2),size(VR,2)),R(size(VR,1),size(VR,2),3), &
            x2p(size(x2,1)),x2f(size(x2,1),size(x2,2)),x2fp(size(x2,1),size(x2,2))
        integer :: i,j

        x2p = (x2(:,1) + x2(:,2) + x2(:,3))/3d0;
        x2f(:,1) = (x2(:,1) + x2(:,2))/2d0;
        x2f(:,2) = (x2(:,2) + x2(:,3))/2d0;
        x2f(:,3) = (x2(:,3) + x2(:,1))/2d0;

        do i=1,size(x2fp,2)
            x2fp(:,i) = x2f(:,i)-x2p;
        enddo

        V(:,1,1) = -x2p;
        V(:,2,1) = 1;
        V(:,1,2) = 1;

        R(:,1,:) = 1;
        R(:,2,:) = 1/4d0*x2fp;

        do i=1,size(Pfp,1)
            V(:,2*i+1,2*i+1) = gP2p(i,:)
            V(:,2*i+2,2*i+2) = gQ2p(i,:)
            do j=1,size(R,3)
                R(:,2*i+1,j) = Pfp(i,:,j)-Pwa(i,:);
                R(:,2*i+2,j) = Qfp(i,:,j)-Qwa(i,:);
            enddo
        enddo
        VR = matrix2dmult(V,R);
    endsubroutine calcVR

    subroutine perot(mesh)
        type(grid_structure),intent(inout) :: mesh
        real*8 :: pqtri(mesh%nt,3,3),pql(mesh%nt,3,2),pp(mesh%nt,3),localvec(mesh%nt,3,2)
        real*8,allocatable :: P(:,:,:),Q(:,:,:),gradP(:,:,:,:),gradQ(:,:,:,:),gP2(:,:,:),gP2p(:,:)
        real*8,allocatable :: Pc(:,:,:),Qc(:,:,:),gradPc(:,:,:,:),gradQc(:,:,:,:),gQ2(:,:,:),gQ2p(:,:)
        real*8,allocatable :: Pf(:,:,:),Qf(:,:,:),Ppc(:,:),Qpc(:,:)
        real*8,allocatable :: Pfp(:,:,:),Qfp(:,:,:),Pwa(:,:),Qwa(:,:)
        real*8             :: A(mesh%nt,3),x2(mesh%nt,3),x2p(mesh%nt)
        integer :: i,j,k
        integer,parameter :: deg=3
        real*8             :: VR(mesh%nt,2*deg+2,3)

        do i =1,3
            do j=1,3
                pqtri(:,i,j) = mesh%v(mesh%tr(:)%v(i))%p(j);
            enddo
            localvec(:,i,1) = mesh%tr(:)%b%lonvec(i);
            localvec(:,i,2) = mesh%tr(:)%b%latvec(i);
            pp(:,i) = mesh%tr(:)%b%p(i);
        enddo

        call localcoord(pqtri,pp,pql,localvec);

        x2 = pql(:,:,1)**2 + pql(:,:,2)**2
        x2p = (x2(:,1) + x2(:,2) + x2(:,3))/3d0;
        ! ------------------------------------------------
        allocate(P(deg,deg+1,deg+1),Q(deg,deg+1,deg+1), &
            gradP(deg,deg+1,deg+1,2),gradQ(deg,deg+1,deg+1,2));
        allocate(Pc(deg,mesh%nt,3),Qc(deg,mesh%nt,3), &
            gradPc(deg,mesh%nt,3,2),gradQc(deg,mesh%nt,3,2));
        allocate(Pf(deg,mesh%nt,3),Qf(deg,mesh%nt,3), &
            Ppc(deg,mesh%nt),Qpc(deg,mesh%nt));
        allocate(Pfp(deg,mesh%nt,3),Qfp(deg,mesh%nt,3), &
            Pwa(deg,mesh%nt),Qwa(deg,mesh%nt));
        allocate(gP2(deg,mesh%nt,3),gQ2(deg,mesh%nt,3), &
            gP2p(deg,mesh%nt),gQ2p(deg,mesh%nt));
        ! ------------------------------------------------

        call computepoln(pql,Pc,Qc,gradPc,gradQc)
        call faceavgpoly(Pc,Qc,Pf,Qf)
        call poligavgpol(Pc,Qc,Ppc,Qpc)

        gP2 = gradPc(:,:,:,1)**2 + gradPc(:,:,:,2)**2;
        gQ2 = gradQc(:,:,:,1)**2 + gradQc(:,:,:,2)**2;

        gP2p = (gP2(:,:,1) + gP2(:,:,2) + gP2(:,:,3))/3d0;
        gQ2p = (gQ2(:,:,1) + gQ2(:,:,2) + gQ2(:,:,3))/3d0;

        do i=1,3
            Pfp(:,:,i) = Pf(:,:,i) -Ppc;
            Qfp(:,:,i) = Qf(:,:,i) -Qpc;
        enddo

        A(:,1) = sphtriarea2d(pp, pqtri(:,1,:), pqtri(:,2,:));
        A(:,2) = sphtriarea2d(pp, pqtri(:,2,:), pqtri(:,3,:));
        A(:,3) = sphtriarea2d(pp, pqtri(:,3,:), pqtri(:,1,:));

        do i =1,3
            do j=1,deg
                Pwa(j,:) = Pwa(j,:) + Pfp(j,:,i)*A(:,i)/mesh%tr(:)%areag;
                Qwa(j,:) = Qwa(j,:) + Qfp(j,:,i)*A(:,i)/mesh%tr(:)%areag;
            enddo
        enddo

        call calcVR(Pwa,Qwa,Pfp,Qfp,x2,gP2p,gQ2p,VR)


        do i=1,mesh%nt
            allocate(mesh%tr(i)%VR(2*deg+2,3));
        enddo

        do k=1,mesh%nt
            do i=1,2*deg+2
                do j=1,3
                    mesh%tr(k)%VR(i,j) = VR(k,i,j)/mesh%tr(k)%areag;
                enddo
            enddo
        enddo

    endsubroutine perot

    subroutine computepoln(pql,P3c,Q3c,gP3c,gQ3c)
        implicit none
        real*8,intent(inout) :: pql(:,:,:),P3c(:,:,:),Q3c(:,:,:),gP3c(:,:,:,:),gQ3c(:,:,:,:)
        real*8,allocatable  :: x(:,:,:),x2(:,:),x2p(:),x2x(:,:,:),x2xp(:,:)
        real*8,allocatable  :: c(:,:)
        real*8,allocatable  :: P2(:,:,:),Q2(:,:,:),gP2(:,:,:,:),gQ2(:,:,:,:)
        real*8,allocatable  :: P2c(:,:,:),Q2c(:,:,:),gP2c(:,:,:,:),gQ2c(:,:,:,:)
        real*8,allocatable  :: P22c(:,:,:),Q22c(:,:,:),P22cp(:,:),Q22cp(:,:)
        integer :: i,ni,nj,nk,nl

        ni = size(gP3c,1);
        nj = size(gP3c,2);
        nk = size(gP3c,3);
        nl = size(gP3c,4);

        allocate(P2(2,3,3),Q2(2,3,3), &
            gP2(2,3,3,2),gQ2(2,3,3,2));

        allocate(P2c(2,nj,nk),Q2c(2,nj,nk), &
            gP2c(2,nj,nk,nl),gQ2c(2,nj,nk,nl))
        allocate(P22c(2,nj,nk),Q22c(2,nj,nk), &
            P22cp(2,nj),Q22cp(2,nj))
        allocate(x(nj,nk,nl),x2(nj,nk),x2p(nj),x2x(nj,nk,nl),c(nj,nl), &
            x2xp(nj,nl))
        ! ------------------------------------
        call harmpol(P2,Q2,gP2,gQ2,2); ! define poly of degree 2
        call computepol(pql,P2,Q2,gP2,gQ2,P2c,Q2c,gP2c,gQ2c) ! Compute degree 2 poly on points

        x = pql;

        P22c = P2c**2;
        Q22c = Q2c**2;
        P22c = 0d0 ! square 2 poly P
        Q22c = 0d0 ! square 2 poly Q
        x2 = 0d0;  ! squared pts
        do i=1,nl
            x2 = x2 + x(:,:,i)**2; 
        enddo

        do i=1,nl ! iterate axis
            x2x(:,:,i) = x2*x(:,:,i); ! numerator of C
        enddo

        P22cp = 0d0 ! centroid square 2 poly P
        Q22cp = 0d0 ! centroid square 2 poly Q
        x2xp =0d0; ! numerator of c
        x2p =0d0; ! 
        do i=1,nk ! iterate vertices
            x2xp = x2xp + x2x(:,i,:)/3d0
            x2p = x2p + x2(:,i)/3d0;
            P22cp = P22cp + P22c(:,:,i)/3d0;
            Q22cp = Q22cp + Q22c(:,:,i)/3d0;
        enddo

        do i=1,nl !iterate axis
            c(:,i) = 1.5d0*x2xp(:,i)/x2p; ! calc c
        enddo

        do i=1,nk
            ! write(*,*) maxval((x(:,:,1)-c(:,1)))
            P3c(ni,:,i) = (x(:,i,1)-c(:,1))*(P2c(2,:,i)-3d0*P22cp(2,:)) &
                - (x(:,i,2)-c(:,2))*(Q2c(2,:,i)-3d0*Q22cp(2,:))
            Q3c(ni,:,i) = (x(:,i,2)-c(:,2))*(P2c(2,:,i)-3d0*P22cp(2,:)) &
                + (x(:,i,1)-c(:,1))*(Q2c(2,:,i)-3d0*Q22cp(2,:))
        enddo

        gP3c(3,:,:,1) = 4d0*P2c(2,:,:)
        gP3c(3,:,:,2) = -4d0*Q2c(2,:,:)
        gQ3c(3,:,:,1) = 4d0*Q2c(2,:,:)
        gQ3c(3,:,:,2) = 4d0*P2c(2,:,:)

        P3c(1:2,:,:) = P2c;
        Q3c(1:2,:,:) = Q2c;
        gP3c(1:2,:,:,:) = gP2c;
        gQ3c(1:2,:,:,:) = gQ2c;
    endsubroutine computepoln

    function getedindexontr(ed, cell, mesh)
        implicit none
        !----------------------------------------------------------------------
        ! getedindexontr
        !   Given an edge number (global index) and a cell number (global)
        !   returns the local index of the edge with respect to cell
        !   returns 0 if edge not in cell
        !   To be used only after mesh is fully structured
        !----------------------------------------------------------------------
        integer, intent(in) :: ed
        integer, intent(in) :: cell
        type(grid_structure), intent(in) :: mesh

        integer:: getedindexontr

        integer:: i

        getedindexontr=0 !index of cell with respect to cell

        do i=1, 3
            !Find the index of edge ed within cell(i)
            if(ed==mesh%tr(cell)%ed(i))then
            getedindexontr=i
            exit
            end if
        end do
            
        return
    end function getedindexontr

    function getedindexonhx(ed, cell, mesh)
        !----------------------------------------------------------------------
        ! getedindexonhx
        !   Given an edge number (global index) and a cell number (global)
        !   returns the local index of the edge with respect to cell
        !   returns 0 if edge not in cell
        !   To be used only after mesh is fully structured
        !----------------------------------------------------------------------
        integer (i4), intent(in) :: ed
        integer (i4), intent(in) :: cell
        type(grid_structure), intent(in) :: mesh

        integer (i4):: getedindexonhx

        integer (i4):: i

        getedindexonhx=0 !index of cell with respect to cell
        do i=1, mesh%v(cell)%nnb
            !Find the index of edge ed within cell(i)
            if(ed==mesh%v(cell)%ed(i))then
                getedindexonhx=i
                exit
            end if
        end do

        return
    end function getedindexonhx
end module smeshpack
