module swmc_operators
    !=============================================================================
    !  Global operators for shallow water model
    !
    ! Pedro da Silva Peixoto (pedrosp@ime.usp.br)
    ! Oct 2018
    !=============================================================================

    !Use main grid data structures
    use datastruct, only: &
        grid_structure, &
        scalar_field, &
        vectorinterpol_methods, &
        vector_field_cart
  
    use constants
    use smeshpack
    use basic_funcs

    contains
        subroutine perot1cor(u,utr,mesh)
            implicit none
            type(grid_structure), intent(in) :: mesh
            type(scalar_field), intent(in):: u
            type(vector_field_cart), intent(inout):: utr
            real*8 :: Ai,dudrk,Pdudrk(3),dei,dek,Punr,nri(3),ui
            integer:: i,j,etrii

            
            do i=1,mesh%nt
                Pdudrk = 0
                Ai = mesh%tr(i)%areag
                do j =1,3
                    etrii = mesh%tr(i)%ed(j)
                    ui = u%f(etrii)
                    nri = mesh%ed(etrii)%nr
                    Punr = sum(utr%p(i)%v*nri)
                    dudrk = ui-Punr
                    dei = mesh%ed(etrii)%leng
                    dek = arclen(mesh%tr(i)%c%p,mesh%ed(etrii)%c%p)
                    Pdudrk = Pdudrk + dudrk*(dei*dek)*nri
                enddo
                Pdudrk = Pdudrk/Ai
                ! print*,Pdudrk
                utr%p(i)%v = utr%p(i)%v + Pdudrk
            enddo
            ! stop
        endsubroutine perot1cor

        subroutine divcorrection(u,h,div,dt,mesh)
            implicit none
            type(grid_structure), intent(in) :: mesh
            type(scalar_field), intent(in):: u,h ! scalar at cells
            type(scalar_field) :: u_ed,udiff! scalar at cells
            type(scalar_field), intent(inout):: div !gradient at edges
            type(vector_field_cart) :: u_tr
            real*8,intent(in) :: dt
            integer:: i

            allocate(u_tr%p(mesh%nt),u_ed%f(mesh%ne))

            call antidiff(h,u,udiff,dt,mesh)
            call perot_ed2tr(u, u_tr, mesh)
            do i=1,mesh%nt
                u_tr%p(i)%v = u_tr%p(i)%v*h%f(i)
            enddo
            call perot_tr2ed(u_tr, u_ed, mesh)
            call div_tr(u_ed,div,mesh)

        endsubroutine divcorrection

        subroutine antidiff(hn,u,udiff,dt,mesh)
            implicit none
            type(grid_structure), intent(in) :: mesh
            type(scalar_field), intent(in):: u,hn ! scalar at cells
            type(scalar_field) :: u_ed,grad,div,ugrad,h_inv! scalar at cells
            type(scalar_field), intent(inout):: udiff !gradient at edges
            type(vector_field_cart) :: u_tr
            real*8,intent(in) :: dt

            allocate(u_tr%p(mesh%nt),u_ed%f(mesh%ne),div%f(mesh%ne),grad%f(mesh%ne))


            h_inv%f = 1d0/(sumreal(abs(hn%f(mesh%ed%sh(1))), &
                abs(hn%f(mesh%ed%sh(2))))+1e-8)
            call perot_ed2tr(u, u_tr, mesh)
            call perot_tr2ed(u_tr, u_ed, mesh)
            call div_ed(u_ed, div, mesh)
            call grad_ed(hn, grad, mesh)

            ugrad%f = mulreal(u_ed%f,grad%f)
            ugrad%f = mulreal(grad%f,mesh%edhx(:)%leng)

            udiff%f = sumreal(div%f,ugrad%f)
            udiff%f = dt*.5d0*mulreal(h_inv%f,udiff%f)
            udiff%f = -mulreal(udiff%f,u_ed%f)

        endsubroutine antidiff

        subroutine div_ed(f_ed,div,mesh)
            implicit none
            type(grid_structure),intent(in) :: mesh
            type(scalar_field),intent(in) :: f_ed
            type(scalar_field),intent(inout) :: div
            real*8  :: pp(3),pu(3),A,le,nr(3),signcor
            integer :: i,j,k,l,sh

            div%f=0d0
            !$omp parallel do &
            !$omp default(none) &
            !$omp shared(mesh, f_ed, div) &
            !$omp private(i,j,k, l, sh,signcor,le,A,pp,pu,nr) &
            !$omp schedule(static)
            do i=1,mesh%ne
                A = 0d0
                do j=1,2
                    sh = mesh%ed(i)%sh(j)
                    A = A + mesh%tr(sh)%areag*rad
                    pp = mesh%tr(sh)%c%p;
                    do k=1, 3
                        l=mesh%tr(sh)%ed(k)
                        pu = mesh%ed(l)%c%p
                        nr = mesh%ed(l)%nr
                        le = mesh%ed(l)%leng
                        signcor = dsign(1d0,dot_product(nr,pu-pp));
                        div%f(i)=div%f(i)+signcor*f_ed%f(l)*le
                    end do
                enddo
                div%f(i)=div%f(i)/A
                end do
                !$omp end parallel do
                return

        endsubroutine div_ed

        !----------------------------------------------------------
        subroutine grad_edhx(f, grad, mesh)
            implicit none
            type(grid_structure), intent(in) :: mesh
            type(scalar_field), intent(in):: f ! scalar at cells
            type(scalar_field), intent(inout):: grad !gradient at edges
            integer(i4):: l,sh(2)
            real*8:: signcor,le,pu(3),ut1(3),pq1(3)

            grad%f = 0

            !$omp parallel do &
            !$omp default(none) &
            !$omp shared(mesh, f, grad) &
            !$omp private(le,pu,sh,ut1,pq1,signcor) &
            !$omp schedule(static)
            do l=1,mesh%ne
                le = mesh%ed(l)%leng*rad;
                pu = mesh%ed(l)%c%p;
                sh = mesh%ed(l)%v;
                ut1 = mesh%ed(l)%tg;
                pq1 = mesh%v(sh(1))%p;
                signcor = -dsign(1d0,dot_product(ut1,pq1-pu));
                grad%f(l) = signcor*(f%f(sh(1))-f%f(sh(2)))/le
            enddo
            !$omp end parallel do
        end subroutine grad_edhx

        subroutine grad_ed(f, grad, mesh)
            implicit none
            !---------------------------------------------------------------
            !Calculate gradient on edges based on scalar at cells
            !---------------------------------------------------------------
            type(grid_structure), intent(in) :: mesh
            type(scalar_field), intent(in):: f ! scalar at cells
            type(scalar_field), intent(inout):: grad !gradient at edges

            integer(i4):: l
            real*8:: de

            grad%f = 0
            !$omp parallel do &
            !$omp default(none) &
            !$omp shared(mesh, f, grad) &
            !$omp private(de) &
            !$omp schedule(static)
            do l=1, mesh%ne
                de = mesh%edhx(l)%leng*rad;
            
                grad%f(l)=-(f%f(mesh%ed(l)%sh(2))-f%f(mesh%ed(l)%sh(1)))
                grad%f(l)=grad%f(l)/de;
            end do
            !$omp end parallel do
            return

        end subroutine grad_ed
        !----------------------------------------------------------

        !Divergent-------------------------------------------------
        subroutine div_tr(uh, div, mesh)
            implicit none
            type(grid_structure), intent(in) :: mesh
            type(scalar_field), intent(in):: uh
            type(scalar_field), intent(inout):: div

            integer(i4):: i, j, l
            real*8 :: signcor,A,le,pp(3),pu(3),nr(3)

            div%f=0d0
            !$omp parallel do &
            !$omp default(none) &
            !$omp shared(mesh, uh, div) &
            !$omp private(j, l, signcor,le,A,pp,pu,nr) &
            !$omp schedule(static)
            do i=1,mesh%nt
                A = mesh%tr(i)%areag*rad;
                pp = mesh%tr(i)%c%p;
                do j=1, 3
                l=mesh%tr(i)%ed(j);
                pu = mesh%ed(l)%c%p;
                nr = mesh%ed(l)%nr
                le = mesh%ed(l)%leng;


                signcor = dsign(1d0,dot_product(nr,pu-pp));
                ! signcor = mesh%tr(i)%cor(j);

                div%f(i)=div%f(i)+signcor*uh%f(l)*le;

                end do
                div%f(i)=div%f(i)/A;
            end do
            !$omp end parallel do
            return

        end subroutine div_tr

        subroutine div_hx(uh, div, mesh)
            implicit none
            type(grid_structure), intent(in)  :: mesh
            type(scalar_field), intent(in)    :: uh
            type(scalar_field), intent(inout) :: div

            integer(i4):: i, j, l
            real*8 :: signcor,A,le,pq(3),pu(3),nr(3)

            div%f=0d0;

            !$omp parallel do &
            !$omp default(none) &
            !$omp shared(mesh, uh, div) &
            !$omp private(j, l, signcor,le,A,pq,pu,nr) &
            !$omp schedule(static)
            do i=1,mesh%nv
                A = mesh%hx(i)%areag*rad;
                pq = mesh%v(i)%p;
                do j=1, mesh%v(i)%nnb
                l  = mesh%v(i)%ed(j);
                pu = mesh%ed(l)%c%p;
                nr = mesh%ed(l)%tg;
                le = mesh%edhx(l)%leng;

                signcor = dsign(1d0,dot_product(nr,pu-pq));

                div%f(i)=div%f(i)+signcor*uh%f(l)*le;
                end do
                div%f(i)=div%f(i)/A;
            end do
            !$omp end parallel do

            return

        end subroutine div_hx
        !----------------------------------------------------------

        !----------------------------------------------------------
        subroutine perot_ed2tr(f_ed, f_tr, mesh)
            implicit none
            type(grid_structure), intent(in)       :: mesh
            type(scalar_field), intent(in)         :: f_ed
            type(vector_field_cart), intent(inout) :: f_tr
            integer  :: i,j,k,lloc(3)
            real*8   :: wee(3),le,Ai,pul(3),pp(3)

            lloc= (/3,1,2/);
            !$omp parallel do &
            !$omp default(none) &
            !$omp shared(mesh, f_ed, f_tr,lloc) &
            !$omp private(i, j, k,le,Ai,pp,pul,wee) &
            !$omp schedule(static)
            do i=1, mesh%nt
                f_tr%p(i)%v=0;
                pp = mesh%tr(i)%c%p;
                Ai = mesh%tr(i)%areag;
                do j=1, 3
                k   = mesh%tr(i)%ed(j);
                pul = mesh%ed(k)%c%p;
                le = mesh%ed(k)%leng;

                if (.true.) then
                    wee = le*arclen(pp,pul)*mesh%ed(k)%nr;
                else
                    wee = 2d0*sphtriarea(pp,mesh%v(mesh%ed(k)%v(1))%p,mesh%v(mesh%ed(k)%v(2))%p)*mesh%ed(k)%nr;
                endif

                ! wee = le*mesh%edhx(k)%leng*mesh%ed(k)%nr*.5d0;

                f_tr%p(i)%v = f_tr%p(i)%v + wee*f_ed%f(k);
                enddo
                f_tr%p(i)%v = f_tr%p(i)%v/Ai
                ! f_tr%p(i)%v = proj_vec_sphere(f_tr%p(i)%v,mesh%tr(i)%c%p)/Ai;
            enddo
            !$omp end parallel do

            if (.false.) then
                call perot1cor(f_ed,f_tr,mesh)
            endif    

            return
        end subroutine perot_ed2tr

        subroutine perot_tr2ed(f_tr,f_ed, mesh)
            implicit none
            type(grid_structure), intent(in)    :: mesh
            type(scalar_field), intent(inout)   :: f_ed ! scalar fields
            type(vector_field_cart), intent(in) :: f_tr

            integer  :: i,j,k,edcelli
            real*8   :: wee(3),de,pu(3),ppl(3),signcor,del(2)

            f_ed%f = 0d0;

            !$omp parallel do &
            !$omp default(none) &
            !$omp shared(mesh, f_ed, f_tr) &
            !$omp private(i, j, k, signcor,edcelli,ppl,pu,wee,de,del) &
            !$omp schedule(static)
            do i=1, mesh%ne
                pu = mesh%ed(i)%c%p;
                de = mesh%edhx(i)%leng;
                ! de =0
                del = (/arclen(mesh%tr(mesh%ed(i)%sh(2))%c%p,pu), &
                    arclen(mesh%tr(mesh%ed(i)%sh(1))%c%p,pu)/)
                do j=1, 2
                k   = mesh%ed(i)%sh(j);
                ppl = mesh%tr(k)%c%p;

                edcelli=getedindexontr(i, k, mesh);
                !  signcor = mesh%tr(k)%nr(edcelli);
                signcor = -dsign(1d0,dot_product(mesh%ed(i)%nr,ppl-pu));

                if (.false.) then
                    wee = del(j)*mesh%ed(i)%nr;
                else
                    wee = arclen(ppl,pu)*mesh%ed(i)%nr;
                    ! de = de+arclen(ppl,pu)
                endif

                f_ed%f(i) = f_ed%f(i) +dot_product(f_tr%p(k)%v,wee);

                enddo
                f_ed%f(i) = f_ed%f(i)/de
            enddo
            !$omp end parallel do

            return
        end subroutine perot_tr2ed
        !----------------------------------------------------------

        !----------------------------------------------------------
        ! subroutine perot_ed2v(f_ed, fv_v,mesh)
        !     implicit none
        !     type(grid_structure), intent(in)       :: mesh
        !     type(scalar_field), intent(in)         :: f_ed ! scalar fields
        !     type(vector_field_cart), intent(inout) :: fv_v
        !     integer  :: i,j,k
        !     real*8   :: wee(3),de,pudl(3),pv(3),pu(3),cor,Ai

        !     !$omp parallel do &
        !     !$omp default(none) &
        !     !$omp shared(mesh, f_ed, fv_v) &
        !     !$omp private(i, j, k,wee,Ai,pudl,pu,de,pv,cor) &
        !     !$omp schedule(static)
        !     do i=1, mesh%nv
        !         fv_v%p(i)%v=0d0;
        !         Ai = mesh%hx(i)%areag;!Ai =0d0;

        !         pv  = mesh%v(i)%p;
        !         do j=1, mesh%v(i)%nnb
        !             k    = mesh%v(i)%ed(j);

        !             pudl = mesh%edhx(k)%c%p;
        !             pu   = mesh%ed(k)%c%p;
        !             de   = mesh%edhx(k)%leng;

        !             cor = dsign(1d0,dot_product(pudl-pv,mesh%edhx(k)%nr))
        !             wee = de*(pudl-pv)*cor;

        !             fv_v%p(i)%v = fv_v%p(i)%v + wee*f_ed%f(k)
        !         enddo
        !         fv_v%p(i)%v = proj_vec_sphere(fv_v%p(i)%v,mesh%v(i)%p)/Ai;
        !         fv_v%p(i)%v = cross_product(pv,fv_v%p(i)%v)

        !     enddo
        !     !$omp end parallel do
        !     return
        ! end subroutine perot_ed2v

        ! subroutine perot_v2ed(fv_v,f_ed, mesh)
        !     implicit none
        !     type(grid_structure), intent(in)  :: mesh
        !     type(scalar_field), intent(inout) :: f_ed ! scalar fields
        !     type(vector_field_cart), intent(in) :: fv_v
        !     integer  :: i,j,k
        !     real*8   :: wee(3),le,pudl(3),pvdl(3)

        !     f_ed%f = 0d0;

        !     !$omp parallel do &
        !     !$omp default(none) &
        !     !$omp shared(mesh, f_ed, fv_v) &
        !     !$omp private(i, j, k,wee,le,pvdl,pudl) &
        !     !$omp schedule(static)
        !     do i=1, mesh%ne
        !         pudl   = mesh%ed(i)%c%p;
        !         le = mesh%ed(i)%leng;
        !         do j=1, 2
        !             k   = mesh%ed(i)%v(j);
        !             pvdl = mesh%v(k)%p;


        !             wee = arclen(pvdl,pudl)*mesh%edhx(i)%nr;
        !             ! wee = arclen(pvdl,pudl)*mesh%edhx(i)%tg;
        !             ! wee = arclen(pvdl,pudl)*mesh%ed(i)%tg;
        !             ! wee = mesh%ed(i)%leng/2d0*mesh%edhx(i)%nr;

        !             f_ed%f(i) = f_ed%f(i) +dot_product(fv_v%p(k)%v,wee);

        !         enddo
        !         f_ed%f(i) = f_ed%f(i)/le
        !     enddo
        !     !$omp end parallel do

        !     return
        ! end subroutine perot_v2ed
        !----------------------------------------------------------

        !----------------------------------------------------------
        subroutine perot_ed2v(f_ed, fv_v,mesh)
            implicit none
            type(grid_structure), intent(in)       :: mesh
            type(scalar_field), intent(in)         :: f_ed ! scalar fields
            type(vector_field_cart), intent(inout) :: fv_v
            integer  :: i,j,k
            real*8   :: wee(3),de,pudl(3),pv(3),pu(3),pp1(3),pp2(3),Ai
      
            !$omp parallel do &
            !$omp default(none) &
            !$omp shared(mesh, f_ed, fv_v) &
            !$omp private(i, j, k,wee,Ai,pudl,pu,de,pv,pp1,pp2) &
            !$omp schedule(static)
            do i=1, mesh%nv
              fv_v%p(i)%v=0d0;
              Ai = mesh%hx(i)%areag;!Ai =0d0;
      
              pv  = mesh%v(i)%p;
              do j=1, mesh%v(i)%nnb        
                k    = mesh%v(i)%ed(j);
      
                pudl = mesh%edhx(k)%c%p;
                pu   = mesh%ed(k)%c%p;
                de   = mesh%edhx(k)%leng;
      
                !wee     = de*cross_product(pudl,arclen(pudl,pv)*mesh%edhx(k)%tg);
                wee     = de*cross_product(pudl,arclen(pudl,pv)*mesh%edhx(k)%nr);
      
                fv_v%p(i)%v = fv_v%p(i)%v + wee*f_ed%f(k)
              enddo
              fv_v%p(i)%v = proj_vec_sphere(fv_v%p(i)%v,mesh%v(i)%p)/Ai;
      
            enddo
            !$omp end parallel do
            return
        end subroutine perot_ed2v

        subroutine perot_v2ed(fv_v,f_ed, mesh)
            implicit none
            type(grid_structure), intent(in)    :: mesh
            type(scalar_field), intent(inout)   :: f_ed ! scalar fields
            type(vector_field_cart), intent(in) :: fv_v
            integer  :: i,j,k
            real*8   :: wee(3),le,pudl(3),pvdl(3),signcor
        
            f_ed%f = 0d0;
        
            !$omp parallel do &
            !$omp default(none) &
            !$omp shared(mesh, f_ed, fv_v) &
            !$omp private(i, j, k, signcor,wee,le,pvdl,pudl) &
            !$omp schedule(static)
            do i=1, mesh%ne
                pudl   = mesh%edhx(i)%c%p;
                le = mesh%ed(i)%leng;
                do j=1, 2
                k   = mesh%ed(i)%v(j);
                pvdl = mesh%v(k)%p;
        
        
                wee = arclen(pvdl,pudl)*mesh%edhx(i)%nr;
                !wee = arclen(pvdl,pudl)*mesh%edhx(i)%tg;
        
                f_ed%f(i) = f_ed%f(i) +dot_product(fv_v%p(k)%v,wee);
        
                enddo
                f_ed%f(i) = f_ed%f(i)/le
            enddo
            !$omp end parallel do
        
            return
        end subroutine perot_v2ed
        !----------------------------------------------------------

        !----------------------------------------------------------
        subroutine coriolis_edhx(u, q_ed, uhq_perp, mesh)
            implicit none
            !---------------------------------------------------------------
            !Calculate Coriolis term
            !---------------------------------------------------------------

            type(grid_structure), intent(in)  :: mesh
            type(scalar_field), intent(in)    :: u, q_ed ! scalar fields
            type(scalar_field), intent(inout) :: uhq_perp !gradient at edges

            !Indexes
            integer(i4):: i,j,k,ks !For node values
            integer(i4):: node !Node index
            integer(i4):: edi !Edge relative to cell i

            integer :: nnb

            !Temporary scalars
            real*8:: les,de,trskw
            real*8:: signcor1,signcor2,signcor
            real*8:: qtmp

            real*8 :: pui(3),ppi(3),pus(3),pqs(3)

            uhq_perp%f=0d0
            !$omp parallel do &
            !$omp default(none) &
            !$omp shared(mesh, u, q_ed,uhq_perp) &
            !$omp private(i,node,edi,de,nnb) &
            !$omp private(j, k,ks,les,trskw, signcor, qtmp) &
            !$omp private(signcor1, signcor2,pui,pqs,pus,ppi) &
            !$omp schedule(static)
            do i=1, mesh%ne
                de = mesh%ed(i)%leng;
                pui = mesh%ed(i)%c%p;
                do j=1,2 !cells sharing edge l
                node = mesh%ed(i)%v(j) ! triangle

                edi  = getedindexonhx(i, node, mesh);
                nnb  = mesh%v(node)%nnb;

                pqs = mesh%v(node)%p;
                ppi = mesh%tr(mesh%v(node)%tr(edi))%c%p;
                signcor1 = -dsign(1d0,dot_product(mesh%ed(i)%nr,ppi-pui));

                do k=1, nnb
                    ks=mesh%v(node)%ed(k);

                    les = mesh%edhx(ks)%leng;
                    pus = mesh%ed(ks)%c%p;
                    signcor2 = dsign(1d0,dot_product(mesh%ed(ks)%tg,pqs-pus));

                    qtmp=0.5d0*(q_ed%f(i)+q_ed%f(ks));
                    signcor=signcor1*signcor2;

                    ! trskw = mesh%hx(node)%trskwg(edi,k);
                    trskw = mesh%hx(node)%trskwg(k,edi);

                    uhq_perp%f(i) = uhq_perp%f(i) + u%f(ks)*qtmp*les*trskw*signcor;
                    ! if (i==7064) then
                    !   write(*,*) j,u%f(ks),trskw,les,signcor1,signcor2,de
                    ! endif
                enddo
                end do
                uhq_perp%f(i)=uhq_perp%f(i)/de;
            end do
            !$omp end parallel do
            return
        end subroutine coriolis_edhx

        subroutine coriolis_ed(u, q_ed, uhq_perp, mesh)
            implicit none
            !---------------------------------------------------------------
            !Calculate Coriolis term
            !---------------------------------------------------------------

            type(grid_structure), intent(in)  :: mesh
            type(scalar_field), intent(in)    :: u, q_ed ! scalar fields
            type(scalar_field), intent(inout) :: uhq_perp !gradient at edges

            !Indexes
            integer(i4):: i,j,k,ks !For node values
            integer(i4):: node !Node index
            integer(i4):: edi !Edge relative to cell i

            integer :: lloce(3,2),llocv(3,2)
            real*8  :: Riv

            !Temporary scalars
            real*8:: les,de,trskw
            real*8:: signcor1,signcor2,signcor
            real*8:: qtmp

            real*8 :: pui(3),pqi(3),pus(3),pps(3)

            lloce(1,:) = (/2,3/);
            lloce(2,:) = (/3,1/);
            lloce(3,:) = (/1,2/);

            llocv(1,:) = (/3,1/);
            llocv(2,:) = (/1,2/);
            llocv(3,:) = (/2,3/);

            
            uhq_perp%f=0d0
            !$omp parallel do &
            !$omp default(none) &
            !$omp shared(mesh, u, q_ed,uhq_perp,lloce,llocv) &
            !$omp private(i,node,edi,de) &
            !$omp private(j, k,ks,les,trskw, signcor, qtmp) &
            !$omp private(signcor1, signcor2,pui,pps,pus,pqi,Riv) &
            !$omp schedule(static)
            do i=1, mesh%ne
                de = mesh%edhx(i)%leng;
                pui = mesh%ed(i)%c%p;
                do j=1,2 !cells sharing edge l
                node=mesh%ed(i)%sh(j) ! triangle
                edi=getedindexontr(i, node, mesh);
                pqi = mesh%v(mesh%tr(node)%v(llocv(edi,2)))%p;
                
                signcor1 = -dsign(1d0,dot_product(mesh%ed(i)%tg,pqi-pui));
                do k=1, 2 ! triangle
                    ks=mesh%tr(node)%ed(lloce(edi,k));
                    
                    les = mesh%ed(ks)%leng;
                    pus = mesh%ed(ks)%c%p;
                    pps = mesh%tr(node)%c%p;
                    signcor2 = dsign(1d0,dot_product(mesh%ed(ks)%nr,pps-pus));          
                    
                    Riv = sum(mesh%tr(node)%trhx_areag(llocv(edi,k:2)));
                    qtmp=0.5d0*(q_ed%f(i)+q_ed%f(ks));
                    signcor=signcor1*signcor2;

                    !  signcor=-signcor1*dsign(1d0,1d0*mesh%tr(node)%nr(lloce(edi,k)));

                    trskw = (Riv/mesh%tr(node)%areag-.5d0);
                    
                    uhq_perp%f(i) = uhq_perp%f(i) + u%f(ks)*qtmp*les*trskw*signcor;
                enddo
                end do
                uhq_perp%f(i)=uhq_perp%f(i)/de;
            end do
            !$omp end parallel do
            return

        end subroutine coriolis_ed
        !----------------------------------------------------------

        !----------------------------------------------------------
        subroutine scalar_ve2ed(fve, fed, mesh)
            implicit none
            !---------------------------------------------------------------
            !Interpolate from triang centers to edges (linear interpolation)
            ! in: ftr - scalar field defined at triangles
            ! out: fed - scalar field defined at edges (must already be allocated)
            !---------------------------------------------------------------
            
            type(grid_structure), intent(in) :: mesh
            type(scalar_field), intent(in):: fve
            type(scalar_field), intent(inout) :: fed
            integer(i4)::  i

            !$omp parallel do &
            !$omp default(none) &
            !$omp shared(mesh,  fed, fve) &
            !$omp private(i) &
            !$omp schedule(static)      
            do i=1,mesh%ne
                fed%f(i)=0.5d0*(fve%f(mesh%ed(i)%v(1))+fve%f(mesh%ed(i)%v(2)))
            end do
            !$omp end parallel do

            return

        end subroutine scalar_ve2ed

        subroutine scalar_tr2ed(ftr, fed, mesh)
            implicit none
            !---------------------------------------------------------------
            !Interpolate from triang centers to edges (linear interpolation)
            ! in: ftr - scalar field defined at triangles
            ! out: fed - scalar field defined at edges (must already be allocated)
            !---------------------------------------------------------------
            
            type(grid_structure), intent(in) :: mesh
            type(scalar_field), intent(in):: ftr
            type(scalar_field), intent(inout) :: fed

            integer(i4)::  l,sh(2)
            real*8:: pp(3),pu(3),p1(3),p2(3),ed_area1,ed_area2

            fed%f = 0d0
            !$omp parallel do &
            !$omp default(none) &
            !$omp shared(mesh,  fed, ftr) &
            !$omp private(l,pp,pu,p1,p2,ed_area1,ed_area2,sh) &
            !$omp schedule(static)

            do l=1,mesh%ne
                pu=mesh%ed(l)%c%p;
                !Calculate PV on edges      
                fed%f(l)=0.5d0*(ftr%f(mesh%ed(l)%sh(1))+ftr%f(mesh%ed(l)%sh(2)))
                
                if (.true.) then
                    sh = mesh%ed(l)%sh;
                    p1=mesh%v(mesh%ed(l)%v(1))%p;
                    p2=mesh%v(mesh%ed(l)%v(2))%p;

                    pp=mesh%tr(sh(1))%c%p;
                    ! ed_area1=sphtriarea(pp, pu, p1)+sphtriarea(pp, pu, p2);
                    ed_area1=sphtriarea(pp, p1, p2);


                    pp=mesh%tr(sh(2))%c%p;
                    ! ed_area2=sphtriarea(pp, pu, p1)+sphtriarea(pp, pu, p2);
                    ed_area2=sphtriarea(pp, p1, p2);

                    ! fed%f(l) = (ed_area1*ftr%f(sh(1))+ed_area2*ftr%f(sh(2)))/(ed_area1+ed_area2);
                    fed%f(l) = (ed_area1*ftr%f(mesh%ed(l)%sh(1))+ &
                        ed_area2*ftr%f(mesh%ed(l)%sh(2)))/(ed_area1+ed_area2);
                    ! fed%f(l) = (ed_area2*ftr%f(mesh%ed(l)%sh(1))+ed_area1*ftr%f(mesh%ed(l)%sh(2)))/(ed_area1+ed_area2);


                endif
                
            end do
            !$omp end parallel do

            return

        end subroutine scalar_tr2ed
        !----------------------------------------------------------

        !----------------------------------------------------------
        subroutine scalar_tr2ve(ftr, fve, mesh)
            implicit none
            !---------------------------------------------------------------
            !Interpolate from triang centers to edges (linear interpolation)
            ! in: ftr - scalar field defined at triangles
            ! out: fed - scalar field defined at edges (must already be allocated)
            !---------------------------------------------------------------
            
            type(grid_structure), intent(in) :: mesh
            type(scalar_field), intent(in):: ftr
            type(scalar_field), intent(inout) :: fve
            real*8  :: inter!,temp(mesh%nt,3)
            integer :: trp(10),pv

            integer(i4)::  l,k

            fve%f = 0
            !$omp parallel do &
            !$omp default(none) &
            !$omp shared(mesh,  ftr, fve) &
            !$omp private(l,inter,pv,trp) &
            !$omp schedule(static)

            do k=1,mesh%nv
                trp = 0
                trp(1:mesh%v(k)%nnb) = mesh%v(k)%tr;
                do l=1,mesh%v(k)%nnb
                pv = minloc(abs(k-mesh%tr(trp(l))%v),1);
                inter = mesh%tr(trp(l))%trhx_areag(pv);
                fve%f(k) = fve%f(k) + inter*ftr%f(trp(l));
                !  fve%f(k) = fve%f(k) + ftr%f(trp(l));
                enddo
                fve%f(k) = fve%f(k)/mesh%hx(k)%areag;
                ! fve%f(k) = fve%f(k)/mesh%v(k)%nnb;
                
            enddo
            !$omp end parallel do

            return

        end subroutine scalar_tr2ve

        subroutine scalar_ve2tr(fve, ftr, mesh)
            implicit none
            !---------------------------------------------------------------
            !Interpolate from triang centers to edges (linear interpolation)
            ! in: ftr - scalar field defined at triangles
            ! out: fed - scalar field defined at edges (must already be allocated)
            !---------------------------------------------------------------
            
            type(grid_structure), intent(in) :: mesh
            type(scalar_field), intent(in):: fve
            type(scalar_field), intent(inout) :: ftr
            real*8  :: A
            integer :: i,j,tri

            ftr%f = 0
            !$omp parallel do &
            !$omp default(none) &
            !$omp shared(mesh,  ftr, fve) &
            !$omp private(i,j,tri,A) &
            !$omp schedule(static)

            do i=1,mesh%nt
                A = mesh%tr(i)%areag
                do j=1,3
                tri = mesh%tr(i)%v(j);
                ftr%f(i) = ftr%f(i)+mesh%tr(i)%trhx_areag(j)*fve%f(tri);
                enddo
                ftr%f(i) = ftr%f(i)/A
            enddo
            !$omp end parallel do

            return

        end subroutine scalar_ve2tr
        !----------------------------------------------------------

        !----------------------------------------------------------
        subroutine apvm(u, q,dt, q_ed, mesh)
            implicit none
            type(grid_structure), intent(in)    :: mesh
            type(scalar_field), intent(in)      :: u,q
            type(scalar_field), intent(inout)   :: q_ed
            type(scalar_field)                  :: dqdx
            real*8      :: pvspar,dt

            integer(i4)::  l,sh(2)

            allocate(dqdx%f(mesh%ne))
            call grad_edhx(q, dqdx, mesh)
            
            
            q_ed%f = 0d0;
            do l=1, mesh%ne
            sh  = mesh%ed(l)%v;


            pvspar = .5d0*u%f(l)*dqdx%f(l);
            q_ed%f(l) = .5d0*sum(q%f(sh))-pvspar*dt;
            end do

            return

        end subroutine apvm
        !----------------------------------------------------------

        !----------------------------------------------------------
        subroutine vort_hx(u, zeta,mesh)
            implicit none

            type(grid_structure), intent(in) :: mesh
            type(scalar_field), intent(in)   :: u
            type(scalar_field), intent(inout):: zeta

            integer(i4):: k, l, ed
            real*8   :: A,signcor,de,pq(3),pu(3)

            zeta%f=0d0

            !$omp parallel do &
            !$omp default(none) &
            !$omp shared(mesh, u, zeta) &
            !$omp private(k,l,A, ed, signcor) &
            !$omp private(de,pq,pu) &
            !$omp schedule(static)
            do k=1, mesh%nv
                A = mesh%hx(k)%areag*rad**2;
                pq = mesh%v(k)%p;

                !loop over triangle edges
                do l=1, mesh%v(k)%nnb
                ed=mesh%v(k)%ed(l);
                pu = mesh%edhx(ed)%c%p;

                de = mesh%edhx(ed)%leng*rad;
                ! signcor=mesh%hx(k)%cor(l);
                signcor=dsign(1d0,dot_product(mesh%ed(ed)%tg,pq-pu));

                zeta%f(k)=zeta%f(k)+u%f(ed)*de*signcor/A;
                
                end do
            end do
            !$omp end parallel do
            return

        end subroutine vort_hx

        subroutine vort_tr(u, zeta, mesh)
            implicit none
            type(grid_structure), intent(in) :: mesh
            type(scalar_field), intent(in):: u
            type(scalar_field), intent(inout):: zeta

            integer(i4):: i, j, l
            real*8 :: signcor,A,le,pp(3),pu(3),nr(3),tg(3)

            zeta%f=0d0
            !$omp parallel do &
            !$omp default(none) &
            !$omp shared(mesh, u, zeta) &
            !$omp private(j, l, signcor,le,A,pp,pu,nr,tg) &
            !$omp schedule(static)
            do i=1,mesh%nt
                A = mesh%tr(i)%areag*rad**2;
                pp = mesh%tr(i)%c%p;
                do j=1, 3
                l=mesh%tr(i)%ed(j);
                pu = mesh%ed(l)%c%p;
                tg = mesh%ed(l)%tg
                nr = mesh%ed(l)%nr
                le = mesh%ed(l)%leng*rad;


                signcor = dsign(1d0,dot_product(nr,pu-pp));

                zeta%f(i)=zeta%f(i)+signcor*u%f(l)*le;
                end do
                zeta%f(i)=zeta%f(i)/A;
            end do
            !$omp end parallel do

            return

        end subroutine vort_tr
        !----------------------------------------------------------

        !--------------------------------------------
        subroutine kinetic_energy_v(u, ke, mesh)
            implicit none
            !---------------------------------------------------------------
            !Calculate kinetic energy at triangles
            ! in: u at cell edges
            ! out: ke at triangles
            !---------------------------------------------------------------
            type(grid_structure), intent(in) :: mesh
            type(scalar_field), intent(in):: u ! velocity at cell edges
            type(scalar_field), intent(inout):: ke !kinetic energy at triangles

            integer(i4):: k, l, ed
            real*8 :: ed_area, A
            real*8 :: pq(3),pp1(3),pp2(3)

            ke%f=0d0
            !$omp parallel do &
            !$omp default(none) &
            !$omp shared(mesh, u, ke) &
            !$omp private(l, ed, A, ed_area,pp1,pp2,pq) &
            !$omp schedule(static)
            do k=1, mesh%nv
            pq = mesh%v(k)%p;
            A = mesh%hx(k)%areag;
            do l=1, mesh%v(k)%nnb
                ed=mesh%v(k)%ed(l);
                pp1 = mesh%tr(mesh%ed(ed)%sh(1))%c%p;
                pp2 = mesh%tr(mesh%ed(ed)%sh(2))%c%p;

                ed_area = mesh%ed(ed)%leng*mesh%edhx(ed)%leng/4d0
                ke%f(k)=ke%f(k)+ed_area*u%f(ed)**2;
            end do
            ke%f(k)=ke%f(k)/A
            end do
            !$omp end parallel do

            return

        end subroutine kinetic_energy_v

        subroutine kinetic_energy_tr(u, ke, mesh)
            implicit none
            !---------------------------------------------------------------
            !Calculate kinetic energy at triangles
            ! in: u at cell edges
            ! out: ke at triangles
            !---------------------------------------------------------------
            type(grid_structure), intent(in) :: mesh
            type(scalar_field), intent(in):: u ! velocity at cell edges
            type(scalar_field), intent(inout):: ke !kinetic energy at triangles

            integer(i4):: k, l, ed
            real*8 :: ed_area, A,pp(1:3),p2(1:3),p3(1:3)
            integer :: vp(2)

            ke%f=0d0
            !$omp parallel do &
            !$omp default(none) &
            !$omp shared(mesh, u, ke) &
            !$omp private(l, ed, A, ed_area, p2, p3,vp,pp) &
            !$omp schedule(static)
            do k=1, mesh%nt
                pp=mesh%tr(k)%c%p;

                !Calculate K energy a la A. Gassman
                A = mesh%tr(k)%areag*rad**2;

                do l=1, 3
                ed=mesh%tr(k)%ed(l);
                if (.true.) then
                    vp = mesh%ed(ed)%v;
                    p2=mesh%v(vp(1))%p;
                    p3=mesh%v(vp(2))%p;
                    ed_area=2d0*sphtriarea(pp, p2, p3)*rad**2;

                else
                    ed_area = mesh%ed(ed)%leng*mesh%edhx(ed)%leng*.5d0*rad**2;
                endif
                    ke%f(k)=ke%f(k)+ed_area*u%f(ed)**2*.5d0;
                end do
                ke%f(k)=ke%f(k)/A
            end do
            !$omp end parallel do

            return

        end subroutine kinetic_energy_tr
        !--------------------------------------------

        !----------------------------------------------------------
        subroutine visc(u,nu,lap,mesh)
            implicit none
            type(grid_structure), intent(in)  :: mesh
            type(scalar_field), intent(in)    :: u
            real*8, intent(in)                :: nu
            type(scalar_field), intent(inout) :: lap
            type(scalar_field)                :: u_perp,div,zeta,q_ed,gradd,gradv
            
            allocate(div%f(mesh%nt),zeta%f(mesh%nv), &
                    gradd%f(mesh%ne),gradv%f(mesh%ne),q_ed%f(mesh%ne), &
                    u_perp%f(mesh%ne));
            q_ed%f = 1d0;

            call coriolis_ed(u, q_ed,u_perp, mesh);

            call div_tr(u, div, mesh);
            call vort_hx(u, zeta, mesh);
            call grad_ed(div, gradd, mesh);
            call grad_edhx(zeta, gradv, mesh);

            lap%f = nu*(gradd%f-gradv%f);

        endsubroutine visc

        subroutine vischx(u,nu,lap,mesh)
            implicit none
            type(grid_structure), intent(in)  :: mesh
            type(scalar_field), intent(in)    :: u
            real*8, intent(in)                :: nu
            type(scalar_field), intent(inout) :: lap
            type(scalar_field)                :: u_perp,div,zeta,q_ed,gradd,gradv

            allocate(div%f(mesh%nv),zeta%f(mesh%nt), &
                    gradd%f(mesh%ne),gradv%f(mesh%ne),q_ed%f(mesh%ne), &
                    u_perp%f(mesh%ne));
            q_ed%f = 1d0;

            ! call coriolis_edhx(u, q_ed,u_perp, mesh);

            call div_hx(u, div, mesh);
            call vort_tr(u, zeta, mesh);
            call grad_edhx(div, gradd, mesh);
            call grad_ed(zeta, gradv, mesh);

            lap%f = nu*(gradd%f-gradv%f);

        endsubroutine vischx
        !----------------------------------------------------------

        !--------------------------------------------
        subroutine calcop(u,h,bath,uh,h_ed,zeta,q,div,uhq_perp,gradb,gradk, &
            ke,nu,lap,bih,f,fv,mesh,dt)

            implicit none
            type(grid_structure), intent(in) :: mesh
            type(scalar_field), intent(in)   :: f,fv,u,h,bath
            type(scalar_field), intent(inout):: div,gradb,gradk,uhq_perp,h_ed,q,zeta, &
                uh,lap,bih
            type(scalar_field)               :: ke,h_ve,q_ed,eta,phi
            type(vector_field_cart)          :: var_tr,uhq_v,varv
            real*8,intent(in) :: dt,nu
            integer :: i
            real*8 :: nu2

            allocate(h_ve%f(mesh%nv),var_tr%p(mesh%nt), &
                    uhq_v%p(mesh%nv),eta%f(mesh%nv), &
                    q_ed%f(mesh%ne),varv%p(mesh%nt), &
                    phi%f(mesh%nt))


            ! ------------------------------------------------------------------------
            h_ed%f = 0d0;
            ke%f     = 0d0;
            lap%f    = 0d0; bih%f    = 0d0;
            q_ed%f   = 0d0; q%f = 0d0;
            eta%f = 0d0; zeta%f = 0d0;
            gradk%f = 0d0;
            ! ------------------------------------------------------------------------

            phi%f = g*(h%f+bath%f)
            if (nlin) then
                !-Kinetic energy--------------------
                if (.false.) then ! kinetic energy
                    call kinetic_energy_tr(u, ke, mesh);
                else             ! Korn, 2017
                    call perot_ed2tr(u, var_tr, mesh);

                    !$omp parallel do &
                    !$omp default(none) &
                    !$omp shared(mesh, ke,var_tr) &
                    !$omp private(i) &
                    !$omp schedule(static)
                    do i =1,mesh%nt
                        ke%f(i) = dsqrt(var_tr%p(i)%v(1)**2+&
                            var_tr%p(i)%v(2)**2+&
                            var_tr%p(i)%v(3)**2);!norm(var_tr%p(i)%v);
                        ke%f(i) = ke%f(i)**2*.5d0;
                    enddo
                    !$omp end parallel do
                endif
                !------------------------------------

                phi%f = phi%f+ke%f
                !-Curl u-----------------------------
                call vort_hx(u, zeta, mesh); ! Curl
                !------------------------------------
                if (.false.) then ! weller u*h
                    call scalar_tr2ed(h, h_ed, mesh)
                    uh%f = u%f*h_ed%f
                elseif (.false.) then
                    call scalar_tr2ed(h, h_ed, mesh)
                    uh%f = u%f*h_ed%f;
                    call perot_ed2tr(uh, var_tr, mesh);
                    call perot_tr2ed(var_tr,uh, mesh);
                else             ! korn u*h
                    call perot_ed2tr(u, var_tr, mesh);
                    do i=1,3
                        var_tr%p(:)%v(i) = var_tr%p(:)%v(i)*h%f;
                    enddo
                    call perot_tr2ed(var_tr,uh, mesh);
                endif
            else
                uh%f = u%f*1d4
                call perot_ed2tr(uh, var_tr, mesh);
                call perot_tr2ed(var_tr,uh, mesh);
            endif
            !-Div--------------------------------
            call div_tr(uh, div, mesh);
            !------------------------------------

            !-grad Phi---------------------------
            if (.true. ) then
                call grad_ed(phi, gradb, mesh);
            else
                phi%f = g*(h%f+bath%f)
                call grad_ed(phi, gradb, mesh);
                call grad_ed(ke, gradk, mesh);
            endif
            !------------------------------------

            !------------------------------------
            call scalar_tr2ve(h, h_ve, mesh)
            q%f = (zeta%f+fv%f)/h_ve%f;
            eta%f = zeta%f+fv%f;
            ! call scalar_tr2ed(h, h_ed, mesh);
            if (.true.) then 
                call scalar_ve2ed(q, q_ed, mesh); ! arith mean
            else
                call apvm(u, q,dt, q_ed, mesh) ! antecipated pot vort
            endif
            !------------------------------------

            !-Coriolis Term----------------------
            if (.false.) then           ! TRSKW
                call coriolis_ed(uh, q_ed, uhq_perp, mesh); ! Nonlinear
            else if (.true.) then      ! Perot
                call perot_ed2v(u, uhq_v,mesh)
                do i =1,3
                    uhq_v%p(:)%v(i) = eta%f*uhq_v%p(:)%v(i);
                enddo
                call perot_v2ed(uhq_v, uhq_perp, mesh);
            else
                call perot_ed2v(uh, uhq_v,mesh)
                do i =1,mesh%nv
                    uhq_v%p(i)%v = q%f(i)*uhq_v%p(i)%v;
                enddo
                call perot_v2ed(uhq_v, uhq_perp, mesh);
            endif
            ! ------------------------------------

            !TRSK
            ! nu = 1d16; call visc(u,nu,lap,mesh);
            ! nu = 1d0;  call visc(lap,nu,bih,mesh);

            ! call visc(u,nu,lap,mesh);
            ! nu2=1d0;
            ! call visc(lap,nu2,bih,mesh);

            ! call visc(bih,nu2,lap,mesh);
            ! call visc(lap,nu2,bih,mesh);

            if (abs(nu)>1d-15) then
                call visc(u,nu,lap,mesh);
                nu2=1d0;
                call visc(lap,nu2,bih,mesh);
            endif

        end subroutine calcop
        !--------------------------------------------


        !--------------------------------------------
        subroutine errorop(u,utr,h,zeta,q,div,uperp,gradb,ke,fv,mesh,dt,file)
            implicit none
            type(grid_structure), intent(in)    :: mesh
            type(scalar_field),   intent(in)    :: fv,u,h
            type(scalar_field)                  :: uh,h_ed,uhq_perp,phi,lap,bih,bath,q_ed,gradk
            type(scalar_field), intent(inout)   :: div,zeta,gradb,uperp,q,ke
            type(vector_field_cart), intent(inout):: utr
            real*8,intent(in)                   :: dt
            real*8                              :: A(mesh%nt),Ae(mesh%ne),Ax(mesh%nv)
            real*8                              :: cen,dTdt,up2l,gr2l,et2l,dv2l,qe2l,grk2l, &
                upmaxl,grmaxl,etmaxl,dvmaxl,qemaxl,grkmaxl,utrmax,utr2l
            integer                             :: i
            character(len=100)                  :: file,files

            type(scalar_field)                  :: momeq;
            real*8                              :: mom2l,mommaxl
            real*8                              :: nu
            real*8                              :: varq(mesh%nv,2),varp(mesh%nt,2), &
                vare(mesh%ne,2),vecp(mesh%nt,3,2)
            type(scalar_field)                  :: f
            type(vector_field_cart) :: up,uv

            nu=1d0;

            allocate(momeq%f(mesh%ne),momeq%fexact(mesh%ne),lap%f(mesh%ne),bih%f(mesh%ne));
            allocate(uh%f(mesh%ne),h_ed%f(mesh%ne));
            allocate(uhq_perp%f(mesh%ne),uhq_perp%fexact(mesh%ne));
            allocate(q_ed%f(mesh%ne),gradk%f(mesh%ne));
            allocate(bath%f(mesh%nt),phi%f(mesh%nt));
            allocate(f%f(mesh%ne));
!            allocate(gradk%f(mesh%ne),f%f(mesh%ne));
            allocate(up%p(mesh%nt),uv%p(mesh%nv));

            A  = mesh%tr(:)%areag;
            Ae = mesh%ed(:)%leng*mesh%edhx(:)%leng;
            Ax = mesh%hx(:)%areag;

            phi%f = g*(bath%f+h%f);
            call calcop(u,h,bath,uh,h_ed,zeta,q,div,uhq_perp,gradb,gradk, &
                ke,nu,lap,bih,f,fv,mesh,dt)
            momeq%f = -(uhq_perp%f-gradb%f-gradk%f);

            call grad_ed(h, gradb, mesh);

            call tenvar(dTdt,uh,gradb,div,phi,mesh);
            
            call perot_ed2v(u, uv,mesh)
            call perot_v2ed(uv, uperp, mesh);
            ! uperp%f = uperp%fexact
            call calc_Corenergies(Cen,u,uperp,mesh)

            call perot_ed2tr(u,up, mesh);
            call perot_tr2ed(up,uh, mesh);

            call div_tr(uh, div, mesh);

            call vort_hx(u, zeta, mesh);
            call grad_ed(h, gradb, mesh);
            
            momeq%fexact = 0d0;

            vare(:,1) = uperp%f; vare(:,2) = uperp%fexact;
            up2l   = l2er(vare,Ae,.true.);
            upmaxl = lmaxer(vare,.false.);


            vare(:,1) = g*gradb%f; vare(:,2) = gradb%fexact;
            gr2l   = l2er(vare,Ae,.true.);
            grmaxl = lmaxer(vare,.false.);

            varp(:,1) = ke%f; varp(:,2) = ke%fexact;
            grk2l   = l2er(varp,A,.true.);
            grkmaxl = lmaxer(varp,.false.);

            varq(:,1) = zeta%f; varq(:,2) = zeta%fexact;
            et2l   = l2er(varq,Ax,.false.)
            etmaxl = lmaxer(varq,.false.)

            varp(:,1) = div%f; varp(:,2) = div%fexact;
            dv2l   = l2er(varp,A,.false.);
            dvmaxl = lmaxer(varp,.false.);

            varq(:,1) = q%f; varq(:,2) = q%fexact;
            qe2l   = l2er(varq,Ax,.false.);
            qemaxl = lmaxer(varq,.false.);

            vare(:,1) = momeq%f; vare(:,2) = momeq%fexact;
            mom2l   = l2er(vare,Ae,.false.);
            mommaxl = lmaxer(vare,.false.);

            do i=1,3
                vecp(:,i,1) = utr%p(:)%v(i);
                vecp(:,i,2) = utr%pexact(:)%v(i);
            enddo
            utr2l = dsqrt(sum(sum((vecp(:,:,1)-vecp(:,:,2))**2,2)*A));
            utrmax = maxval(abs(vecp(:,:,1)-vecp(:,:,2)));

            write(*,'(A)') '---------------------------------------------------------------';
            write(*,'(A)') ' ';
            write(*,'(A,4ES14.7)') 'uperp:',maxval(abs(uperp%f)), &
                maxval(abs(uperp%fexact)),upmaxl,up2l;
            write(*,'(A,4ES14.7)') 'grad: ',g*maxval(abs((gradb%f))), &
                maxval(abs(gradb%fexact)),grmaxl,gr2l;
            write(*,'(A,4ES14.7)') 'ke:   ',maxval(abs((ke%f))), &
                maxval(abs(ke%fexact)),grkmaxl,grk2l;
            write(*,'(A,4ES14.7)') 'div:  ',maxval(abs(div%f)), &
                maxval(abs(div%fexact)),dvmaxl,dv2l;
            write(*,'(A,4ES14.7)') 'zeta: ',maxval(abs(zeta%f)), &
                maxval(abs(zeta%fexact)),etmaxl,et2l;
            write(*,'(A,4ES14.7)') 'q:    ',maxval(abs(q%f)), &
                maxval(abs(q%fexact)),qemaxl,qe2l;
            write(*,'(A,4ES14.7)') 'momeq:',maxval(abs(momeq%f)), &
                maxval(abs(momeq%fexact)),mommaxl,mom2l;
            write(*,'(A,4ES14.7)') 'utr:  ',maxval(normm(vecp(:,:,1))), &
                maxval(normm(vecp(:,:,2))),utrmax,utr2l;
            write(*,'(A)') '---------------------------------------------------------------';
            write(*,'(A)') '-> Energy';
            write(*,'(A,ES14.7)') 'ce:   ',abs(cen);
            write(*,'(A,ES14.7)') 'dTdt: ',abs(dTdt);
            write(*,'(A)') '---------------------------------------------------------------';

            write(files,'(2A)') trim(file),'uperp.dat';
            open(unit = 101, file = files(2:100));

            write(files,'(2A)') trim(file),'gradb.dat';
            open(unit = 102, file = files(2:100));

            write(files,'(2A)') trim(file),'ke.dat';
            open(unit = 103, file = files(2:100));

            write(files,'(2A)') trim(file),'div.dat';
            open(unit = 104, file = files(2:100));

            write(files,'(2A)') trim(file),'zeta.dat';
            open(unit = 105, file = files(2:100));

            write(files,'(2A)') trim(file),'q.dat';
            open(unit = 106, file = files(2:100));

            write(files,'(2A)') trim(file),'momeq.dat';
            open(unit = 107, file = files(2:100));

            write(files,'(2A)') trim(file),'utr.dat';
            open(unit = 108, file = files(2:100));

            do i=1,mesh%nt
                write(103,*) ke%f(i),ke%fexact(i)        
                write(104,*) div%f(i),div%fexact(i)
                write(108,*) utr%p(i)%v(:),utr%pexact(i)%v(:)
            enddo
            do i=1,mesh%ne
                write(101,*) uperp%f(i),uperp%fexact(i)
                write(102,*) gradb%f(i),gradb%fexact(i)
                write(107,*) momeq%f(i),momeq%fexact(i)
            enddo
            do i=1,mesh%nv
                write(105,*) zeta%f(i),zeta%fexact(i)
                write(106,*) q%f(i),q%fexact(i)
            enddo

            do i=101,108
                close(i);
            enddo
        end subroutine errorop

        ! subroutine errorophx(u,h,zeta,q,div,uperp,gradb,ke,fv,mesh,dt,file)
        !     implicit none
        !     type(grid_structure), intent(in)    :: mesh
        !     type(scalar_field),   intent(in)    :: fv,u,h
        !     type(scalar_field)                  :: uh,h_ed,uhq_perp,phi,lap,bih,bath,q_ed
        !     type(scalar_field),   intent(inout) :: div,zeta,gradb,uperp,q,ke
        !     real*8,intent(in)                   :: dt
        !     real*8                              :: A(mesh%nt),Ae(mesh%ne),Ax(mesh%nv)
        !     real*8                              :: cen,dTdt,up2l,gr2l,et2l,dv2l,qe2l,grk2l, &
        !         upmaxl,grmaxl,etmaxl,dvmaxl,qemaxl,grkmaxl
        !     integer                             :: i
        !     character(len=100)                  :: file,files

        !     type(scalar_field)                  :: momeq
        !     real*8                              :: mom2l,mommaxl
        !     real*8                              :: nu
        !     real*8                              :: varq(mesh%nv,2),varp(mesh%nt,2),vare(mesh%ne,2)

        !     nu=1d0;

        !     allocate(momeq%f(mesh%ne),momeq%fexact(mesh%ne),lap%f(mesh%ne),bih%f(mesh%ne));
        !     allocate(uh%f(mesh%ne),h_ed%f(mesh%ne));
        !     allocate(uhq_perp%f(mesh%ne),uhq_perp%fexact(mesh%ne));
        !     allocate(q_ed%f(mesh%ne));
        !     allocate(bath%f(mesh%nv),phi%f(mesh%nv));
            

        !     A  = mesh%tr(:)%areag;
        !     Ae = mesh%ed(:)%leng*mesh%edhx(:)%leng;
        !     Ax = mesh%hx(:)%areag;

        !     phi%f = g*(bath%f+h%f);
        !     call calcophx(u,h,uh,h_ed,bath,zeta,q,div,uhq_perp,gradb,ke,nu,lap,bih,fv,mesh,dt);
        !     momeq%f = -(uhq_perp%f-gradb%f);

        !     phi%f = g*h%f;
        !     call grad_edhx(phi, gradb, mesh);
        !     call tenvar(dTdt,uh,gradb,div,phi,mesh);
        !     call calc_Corenergies(Cen,uh,uhq_perp,mesh)

        !     q_ed%f=1;
        !     call coriolis_edhx(u, q_ed, uperp, mesh);
        !     call div_hx(u, div, mesh);
        !     call vort_tr(u, zeta, mesh);
        !     call kinetic_energy_v(u, ke, mesh);
            
        !     momeq%fexact = 0d0;

        !     vare(:,1) = uperp%f; vare(:,2) = uperp%fexact;
        !     up2l   = l2er(vare,Ae,.false.);
        !     upmaxl = lmaxer(vare,.false.);


        !     vare(:,1) = -gradb%f; vare(:,2) = gradb%fexact;
        !     gr2l   = l2er(vare,Ae,.false.);
        !     grmaxl = lmaxer(vare,.false.);

        !     varq(:,1) = ke%f; varq(:,2) = ke%fexact;
        !     grk2l   = l2er(varq,Ax,.true.);
        !     grkmaxl = lmaxer(varq,.true.);

        !     varp(:,1) = zeta%f; varp(:,2) = zeta%fexact;
        !     et2l   = l2er(varp,A,.true.)
        !     etmaxl = lmaxer(varp,.true.)

        !     varq(:,1) = div%f; varq(:,2) = div%fexact;
        !     dv2l   = l2er(varq,Ax,.false.);
        !     dvmaxl = lmaxer(varq,.false.);

        !     varp(:,1) = q%f; varp(:,2) = q%fexact;
        !     qe2l   = l2er(varp,A,.true.);
        !     qemaxl = lmaxer(varp,.true.);

        !     vare(:,1) = momeq%f; vare(:,2) = momeq%fexact;
        !     mom2l   = l2er(vare,Ae,.false.);
        !     mommaxl = lmaxer(vare,.false.);

        !     write(*,'(A)') '----------------------------------------------------------------------------------------------------';
        !     write(*,'(A)') ' ';
        !     write(*,'(A,4ES14.7)') 'uperp:',maxval(abs(uperp%f)),maxval(abs(uperp%fexact)),upmaxl,up2l;
        !     write(*,'(A,4ES14.7)') 'grad: ',maxval(abs((gradb%f))),maxval(abs(gradb%fexact)),grmaxl,gr2l;
        !     write(*,'(A,4ES14.7)') 'ke:   ',maxval(abs((ke%f))),maxval(abs(ke%fexact)),grkmaxl,grk2l;
        !     write(*,'(A,4ES14.7)') 'div:  ',maxval(abs(div%f)),maxval(abs(div%fexact)),dvmaxl,dv2l;
        !     write(*,'(A,4ES14.7)') 'zeta: ',maxval(abs(zeta%f)),maxval(abs(zeta%fexact)),etmaxl,et2l;
        !     write(*,'(A,4ES14.7)') 'q:    ',maxval(abs(q%f)),maxval(abs(q%fexact)),qemaxl,qe2l;
        !     write(*,'(A,4ES14.7)') 'momeq:',maxval(abs(momeq%f)),maxval(abs(momeq%fexact)),mommaxl,mom2l;
        !     write(*,'(A)') '----------------------------------------------------------------------------------------------------';
        !     write(*,'(A)') '-> Energy';
        !     write(*,'(A,ES14.7)') 'ce:   ',abs(cen);
        !     write(*,'(A,ES14.7)') 'dTdt: ',abs(dTdt);
        !     write(*,'(A)') '----------------------------------------------------------------------------------------------------';

        !     write(files,'(2A)') trim(file),'uhq_perp.dat';
        !     open(unit = 101, file = files(2:100));

        !     write(files,'(2A)') trim(file),'gradb.dat';
        !     open(unit = 102, file = files(2:100));

        !     write(files,'(2A)') trim(file),'ke.dat';
        !     open(unit = 103, file = files(2:100));

        !     write(files,'(2A)') trim(file),'div.dat';
        !     open(unit = 104, file = files(2:100));

        !     write(files,'(2A)') trim(file),'zeta.dat';
        !     open(unit = 105, file = files(2:100));

        !     write(files,'(2A)') trim(file),'q.dat';
        !     open(unit = 106, file = files(2:100));

        !     write(files,'(2A)') trim(file),'momeq.dat';
        !     open(unit = 107, file = files(2:100));

        !     do i=1,mesh%nv
        !         write(103,*) ke%f(i),ke%fexact(i)
        !         write(104,*) div%f(i),div%fexact(i)
        !     enddo
        !     do i=1,mesh%ne
        !         write(101,*) uperp%f(i),uperp%fexact(i)
        !         write(102,*) gradb%f(i),gradb%fexact(i)
        !         write(107,*) momeq%f(i),momeq%fexact(i)
        !     enddo
        !     do i=1,mesh%nt
        !         write(105,*) zeta%f(i),zeta%fexact(i)
        !         write(106,*) q%f(i),q%fexact(i)
        !     enddo

        !     do i=101,107
        !         close(i);
        !     enddo
        ! end subroutine errorophx
        !--------------------------------------------

        function vect_erp1(mesh)
            type(grid_structure), intent(in) :: mesh
            type(scalar_field)               :: vec,vecn
            real*8                           :: vect_erp1(mesh%ne,3)
            type(vector_field_cart)          :: var_tr


            allocate(var_tr%p(mesh%nt),vec%f(mesh%ne),vecn%f(mesh%ne));

            vec%f = 1d0;

            call perot_ed2tr(vec,var_tr, mesh);
            do i =1,mesh%nt
                var_tr%p(i)%v = cross_product(mesh%tr(i)%c%p,var_tr%p(i)%v);
            enddo
            call perot_tr2ed(var_tr,vecn, mesh)

            do i = 1,mesh%ne
                vect_erp1(i,:) = vecn%f(i)*mesh%ed(i)%nr;
            enddo      
        end function vect_erp1

        !--------------------------------------------
        subroutine calc_Corenergies(Cenergy,u,uhq_perp,mesh)
            implicit none
            type(grid_structure), intent(in) :: mesh
            type(scalar_field), intent(in) :: u,uhq_perp
            real*8            :: Ae(mesh%ne)

            real*8, intent(out):: Cenergy

            Ae = mesh%ed(:)%leng*mesh%edhx(:)%leng;
            Cenergy = sum(Ae*uhq_perp%f*u%f);

            return
        end subroutine calc_Corenergies

        !--------------------------------------------
        subroutine tenvar(dTdt,uh,gradb,div,h,mesh)
            implicit none
            type(grid_structure), intent(in) :: mesh
            type(scalar_field), intent(in) :: h,gradb,div,uh
            real*8            :: A(mesh%nt),Ae(mesh%ne),Av(mesh%nv)
            real*8, intent(out):: dTdt

            Av  = mesh%hx(:)%areag;
            A  = mesh%tr(:)%areag;
            Ae = mesh%ed(:)%leng*mesh%edhx(:)%leng;

            dTdt = g*sum(Ae*uh%f*gradb%f)-sum(A*h%f*div%f);
        end subroutine tenvar

        subroutine tenvarhx(dTdt,uh,gradb,div,h,mesh)
            implicit none
            type(grid_structure), intent(in) :: mesh
            type(scalar_field), intent(in) :: h,gradb,div,uh
            real*8            :: A(mesh%nt),Ae(mesh%ne),Av(mesh%nv)
            real*8, intent(out):: dTdt

            Av  = mesh%hx(:)%areag;
            A  = mesh%tr(:)%areag;
            Ae = mesh%ed(:)%leng*mesh%edhx(:)%leng;

            ! dTdt = sum(Ae*u%f*gradb%f*h_ed%f)-sum(Ai*h%f*div%f);
            dTdt = sum(Ae*uh%f*gradb%f)-sum(Av*h%f*div%f);
        end subroutine tenvarhx
        !--------------------------------------------

        !--------------------------------------------
        subroutine saveoperatorstr(u,h,j,filei,mesh)
            type(grid_structure),intent(in) :: mesh
            type(scalar_field),intent(in) :: h,u
            integer,intent(in) :: j
            integer :: i
            character(len=100),intent(in) :: filei
            character(len=100)            :: filen

            print*, 'Saving',j

            call charsave(filei,filen,'h',j,4)
            open(unit = 110, file = filen);

            call charsave(filei,filen,'u',j,4)
            open(unit = 111, file = filen);

            do i=1,mesh%nt
                write(110,*) h%f(i);
            enddo

            do i=1,mesh%ne
                write(111,*) u%f(i)*mesh%ed(i)%nr;
            enddo

            do i=110,111
            close(i);
            enddo
        end subroutine saveoperatorstr

        subroutine saveoperatorsve(pathn,glevel,u,B,l,mesh)
            implicit none
            type(grid_structure), intent(in) :: mesh
            type(scalar_field), intent(in):: B,u
            type(scalar_field)            :: uhq_perp,grad,ke,div,zeta
            integer, intent(in) :: l
            integer             :: i
            character(len=2),intent(in) :: glevel
            character(len=100),intent(in) :: pathn
            character(len=100) :: savefb,savefv,savefg,savefk,savefup,savefdi,savefzt
            
            allocate(uhq_perp%f(mesh%ne),grad%f(mesh%ne),ke%f(mesh%nv),div%f(mesh%nv),zeta%f(mesh%nt));

            call vort_tr(u, zeta, mesh);

            if (l<10) then
                write(savefb,'(A,A,A,A,1I1,A)')  trim(pathn),'snapshots/h/',trim(glevel),'/B000',l,'.dat';
                write(savefv,'(A,A,A,A,1I1,A)')  trim(pathn),'snapshots/v/',trim(glevel),'/v000',l,'.dat';
                write(savefg,'(A,A,A,A,1I1,A)')  trim(pathn),'snapshots/grad/',trim(glevel),'/grad000',l,'.dat';
                write(savefk,'(A,A,A,A,1I1,A)')  trim(pathn),'snapshots/ke/',trim(glevel),'/ke000',l,'.dat';
                write(savefup,'(A,A,A,A,1I1,A)') trim(pathn),'snapshots/uperp/',trim(glevel),'/up000',l,'.dat';
                write(savefdi,'(A,A,A,A,1I1,A)') trim(pathn),'snapshots/div/',trim(glevel),'/div000',l,'.dat';
                write(savefzt,'(A,A,A,A,1I1,A)') trim(pathn),'snapshots/zeta/',trim(glevel),'/zeta000',l,'.dat';
            elseif (l<100) then
                write(savefb,'(A,A,A,A,1I2,A)')  trim(pathn),'snapshots/h/',trim(glevel),'/B00',l,'.dat';
                write(savefv,'(A,A,A,A,1I2,A)')  trim(pathn),'snapshots/v/',trim(glevel),'/v00',l,'.dat';
                write(savefg,'(A,A,A,A,1I2,A)')  trim(pathn),'snapshots/grad/',trim(glevel),'/grad00',l,'.dat';
                write(savefk,'(A,A,A,A,1I2,A)')  trim(pathn),'snapshots/ke/',trim(glevel),'/ke00',l,'.dat';
                write(savefup,'(A,A,A,A,1I2,A)') trim(pathn),'snapshots/uperp/',trim(glevel),'/up00',l,'.dat';
                write(savefdi,'(A,A,A,A,1I2,A)') trim(pathn),'snapshots/div/',trim(glevel),'/div00',l,'.dat';
                write(savefzt,'(A,A,A,A,1I2,A)') trim(pathn),'snapshots/zeta/',trim(glevel),'/zeta00',l,'.dat';

            elseif (l<1000) then
                write(savefb,'(A,A,A,A,1I3,A)')  trim(pathn),'snapshots/h/',trim(glevel),'/B0',l,'.dat';
                write(savefv,'(A,A,A,A,1I3,A)')  trim(pathn),'snapshots/v/',trim(glevel),'/v0',l,'.dat';
                write(savefg,'(A,A,A,A,1I3,A)')  trim(pathn),'snapshots/grad/',trim(glevel),'/grad0',l,'.dat';
                write(savefk,'(A,A,A,A,1I3,A)')  trim(pathn),'snapshots/ke/',trim(glevel),'/ke0',l,'.dat';
                write(savefup,'(A,A,A,A,1I3,A)') trim(pathn),'snapshots/uperp/',trim(glevel),'/up0',l,'.dat';
                write(savefdi,'(A,A,A,A,1I3,A)') trim(pathn),'snapshots/div/',trim(glevel),'/div0',l,'.dat';
                write(savefzt,'(A,A,A,A,1I3,A)') trim(pathn),'snapshots/zeta/',trim(glevel),'/zeta0',l,'.dat';
            else
                write(savefb,'(A,A,A,A,1I4,A)')  trim(pathn),'h/snapshots/',trim(glevel),'/B',l,'.dat';
                write(savefv,'(A,A,A,A,1I4,A)')  trim(pathn),'v/snapshots/',trim(glevel),'/v',l,'.dat';
                write(savefg,'(A,A,A,A,1I4,A)')  trim(pathn),'snapshots/grad/',trim(glevel),'/grad',l,'.dat';
                write(savefk,'(A,A,A,A,1I4,A)')  trim(pathn),'snapshots/ke/',trim(glevel),'/ke',l,'.dat';
                write(savefup,'(A,A,A,A,1I4,A)') trim(pathn),'snapshots/uperp/',trim(glevel),'/up',l,'.dat';
                write(savefdi,'(A,A,A,A,1I4,A)') trim(pathn),'snapshots/div/',trim(glevel),'/div',l,'.dat';
                write(savefzt,'(A,A,A,A,1I4,A)') trim(pathn),'snapshots/zeta/',trim(glevel),'/zeta',l,'.dat';
            endif

            write(*,*) '-->',savefb
            write(*,*) '-->',savefv
            write(*,*) '-->',savefg
            write(*,*) '-->',savefk
            write(*,*) '-->',savefup
            write(*,*) '-->',savefdi
            write(*,*) '-->',savefzt

            
            open(unit=100,file=savefb(2:100));
            open(unit=101,file=savefv(2:100));
            open(unit=102,file=savefg(2:100));
            open(unit=103,file=savefup(2:100));
            open(unit=104,file=savefdi(2:100));
            open(unit=105,file=savefzt(2:100));
            open(unit=106,file=savefk(2:100));

            do i=1,mesh%nv
                write(100,*) B%f(i);
                write(104,*) div%f(i);
                write(106,*) ke%f(i);
            enddo
            do i=1,mesh%ne
                write(101,*) u%f(i);
                write(102,*) grad%f(i);
                write(103,*) uhq_perp%f(i);
            enddo


            do i=1,mesh%nt
                write(105,*) zeta%f(i);
            enddo
            do i = 100,106
                close(i);
            enddo
        end subroutine saveoperatorsve
        !--------------------------------------------

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

end module swmc_operators
