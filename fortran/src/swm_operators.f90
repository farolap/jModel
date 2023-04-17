module swm_operators
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
    vectorinterpol_methods, &
    vector_field_cart

  use smeshpack

  contains
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

      grad%f=-(f%f(mesh%ed(:)%sh(2))-f%f(mesh%ed(:)%sh(1)))/(mesh%edhx(:)%leng*rad);

      return

    end subroutine grad_ed

  subroutine div_tr(uh, div, mesh)
    implicit none

    type(grid_structure), intent(in) :: mesh
    type(scalar_field), intent(in)   :: uh 
    type(scalar_field), intent(inout):: div

    integer(i4):: i, j, l(3)
    real*8 :: A,le(3),pp(3),pu(3,3),nr(3,3),signcor(3)

    div%f=0d0
    !$omp parallel do &
    !$omp default(none) &
    !$omp shared(mesh, uh, div) &
    !$omp private(j, l, signcor,le,A,pp,pu,nr) &
    !$omp schedule(static)
    do i=1,mesh%nt
      A = mesh%tr(i)%areag;
      pp = mesh%tr(i)%c%p;
      l=mesh%tr(i)%ed;

      do j=1,3
        pu(j,:) = mesh%ed(l(j))%c%p;
        nr(j,:) = mesh%ed(l(j))%nr
        signcor(j) =  -dsign(1d0,dot_product(pp-pu(j,:),nr(j,:)))
      enddo
      le      = mesh%ed(l)%leng;

      div%f(i)=sum(signcor*uh%f(l)*le)/(A*rad);

    end do
    !$omp end parallel do

    return

  end subroutine div_tr

  subroutine perot_ed2tr(f_ed, f_tr, mesh)
    implicit none
    type(grid_structure), intent(in)       :: mesh
    type(scalar_field), intent(in)         :: f_ed          ! scalar fields
    type(vector_field_cart), intent(inout) :: f_tr
    integer  :: i,j,k
    real*8   :: wee(3),le,Ai,pul(3),pp(3),v1(3),v2(3)

    !$omp parallel do &
    !$omp default(none) &
    !$omp shared(mesh, f_ed, f_tr) &
    !$omp private(i, j, k,le,Ai,pp,pul,wee,v1,v2) &
    !$omp schedule(static)
    do i=1, mesh%nt
      f_tr%p(i)%v=0;
      pp = mesh%tr(i)%c%p;
      Ai = mesh%tr(i)%areag;
      do j=1, 3
        k   = mesh%tr(i)%ed(j);
        pul = mesh%ed(k)%c%p;
        le = mesh%ed(k)%leng;

        if (.false.) then
          v1 = mesh%v(mesh%ed(k)%v(1))%p;
          v2 = mesh%v(mesh%ed(k)%v(2))%p;
          wee = .5d0*norm(v2-v1)*norm(pp-pul)*mesh%ed(k)%nr;
        else
          wee = le*arclen(pp,pul)*mesh%ed(k)%nr;
        endif

        f_tr%p(i)%v = f_tr%p(i)%v + wee*f_ed%f(k);
      enddo
      f_tr%p(i)%v = f_tr%p(i)%v/Ai
      f_tr%p(i)%v = proj_vec_sphere(f_tr%p(i)%v,mesh%tr(i)%c%p);
    enddo
    !$omp end parallel do
    return
  end subroutine perot_ed2tr

  subroutine perot_tr2ed(f_tr,f_ed, mesh)
    implicit none
    type(grid_structure), intent(in)    :: mesh
    type(scalar_field), intent(inout)   :: f_ed ! scalar fields
    type(vector_field_cart), intent(in) :: f_tr

    integer  :: i,j,k,edcelli
    real*8   :: wee(3),de,pu(3),ppl(3),signcor

    f_ed%f = 0d0;

    !$omp parallel do &
    !$omp default(none) &
    !$omp shared(mesh, f_ed, f_tr) &
    !$omp private(i, j, k, signcor,edcelli,ppl,pu,wee,de) &
    !$omp schedule(static)
    do i=1, mesh%ne
      pu = mesh%ed(i)%c%p;
      de = mesh%edhx(i)%leng;
      do j=1, 2
        k   = mesh%ed(i)%sh(j);
        ppl = mesh%tr(k)%c%p;

        edcelli=getedindexontr(i, k, mesh);
!        signcor = mesh%tr(k)%nr(edcelli);
        signcor = -dsign(1d0,dot_product(mesh%ed(i)%nr,ppl-pu));

        if (.false.) then
          wee = -(ppl-pu)*signcor;
        else
          wee = arclen(ppl,pu)*mesh%ed(i)%nr;
        endif

        f_ed%f(i) = f_ed%f(i) +dot_product(f_tr%p(k)%v,wee);

      enddo
      f_ed%f(i) = f_ed%f(i)/de
    enddo
    !$omp end parallel do

    return
  end subroutine perot_tr2ed

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

        if (.false.) then
          pp1 = mesh%tr(mesh%ed(k)%sh(1))%c%p;
          pp2 = mesh%tr(mesh%ed(k)%sh(2))%c%p;
          wee = norm(pp1-pp2)*norm(pudl-pv)*cross_product(pudl,mesh%edhx(k)%nr);
          Ai  = Ai+norm(pp1-pp2)*norm(pudl-pv)
        else
          wee     = de*cross_product(pudl,arclen(pudl,pv)*mesh%edhx(k)%nr);
        endif

!        write(*,*) dot_product(mesh%edhx(i)%nr,mesh%ed(i)%nr);

        fv_v%p(i)%v = fv_v%p(i)%v + wee*f_ed%f(k)
      enddo
      fv_v%p(i)%v = fv_v%p(i)%v/Ai
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


        if (.false.) then
          signcor = dsign(1d0,dot_product(pvdl-pudl,mesh%ed(i)%tg));
          wee = (pvdl-pudl)*signcor;
        else
          wee = arclen(pvdl,pudl)*mesh%edhx(i)%nr;
        endif

        f_ed%f(i) = f_ed%f(i) +dot_product(fv_v%p(k)%v,wee);

      enddo
      f_ed%f(i) = f_ed%f(i)/le
    enddo
    !$omp end parallel do

    return
  end subroutine perot_v2ed

  subroutine coriolis_edhx(u, q_ed,q, uhq_perp,nlin, mesh)
    implicit none
    type(grid_structure), intent(in) :: mesh
    type(scalar_field), intent(in):: u, q_ed, q ! scalar fields
    type(scalar_field), intent(inout):: uhq_perp !gradient at edges
    integer,intent(in) :: nlin
    integer :: i,j,k,dshi,indeu,eds,vds,indm
    real*8  :: le,wee,des,signcor,pps(3),pus(3),pq(3),pu(3),tv,qtmp,tg(3),nrs(3),Riv

    uhq_perp%f =0d0

    !$omp parallel do &
    !$omp default(none) &
    !$omp shared(mesh, u, q_ed, q, nlin,uhq_perp) &
    !$omp private(i, j, k,le,pu,tg,dshi,indm,pq,indeu,tv) &
    !$omp private(eds,vds,pps,pus,nrs,des,qtmp,Riv,signcor,wee) &
    !$omp schedule(static)
    do i = 1,mesh%ne
      le = mesh%ed(i)%leng;
      pu = mesh%ed(i)%c%p;
      tg = mesh%ed(i)%tg;
      do j = 1,2
        dshi = mesh%ed(i)%v(j);
        indm = mesh%v(dshi)%nnb;
        pq = mesh%v(dshi)%p;
        indeu = minloc(abs(mesh%v(dshi)%ed-i),1);            ! Main Edge  local index
        tv = -dsign(1d0,dot_product(tg,pq-pu));
        do k=1,indm
          eds = mesh%v(dshi)%ed(k);                          ! global index of edge
          vds = mesh%hx(dshi)%v(k);                          ! global index of vertex/triangle center
          pps = mesh%tr(vds)%c%p;                            ! vertex s position of voronoi
          pus = mesh%ed(eds)%c%p;                            ! edge s position
          nrs = mesh%ed(eds)%nr;                             ! edge s position tangent vecto
          des = mesh%edhx(eds)%leng;                         ! voronoi edge length

          if (.false.) then
            qtmp = 0.5d0*(q_ed%f(i)+q_ed%f(eds));
          else
            qtmp = q%f(dshi);
          endif
!          qtmp=1d0;
          Riv = mesh%hx(dshi)%trskwg(indeu,k);
          signcor = dsign(1d0,dot_product(nrs,pps-pus))*tv;
          wee = Riv*signcor;
          uhq_perp%f(i) = uhq_perp%f(i)+wee*des*u%f(eds)*qtmp/le;
        enddo
      enddo
    enddo
    !$omp end parallel do
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
!          signcor=-signcor1*dsign(1d0,1d0*mesh%tr(node)%nr(lloce(edi,k)));

          trskw = (Riv/mesh%tr(node)%areag-.5d0);

          uhq_perp%f(i) = uhq_perp%f(i) + u%f(ks)*qtmp*les*trskw*signcor;
        enddo
      end do
      uhq_perp%f(i)=uhq_perp%f(i)/de;
    end do
    !$omp end parallel do
    return

  end subroutine coriolis_ed

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
      !Calculate PV on edges
      fed%f(i)=0.5d0*(fve%f(mesh%ed(i)%v(1))+fve%f(mesh%ed(i)%v(2)))
    end do
    !$omp end parallel do

    return

  end subroutine scalar_ve2ed

!  subroutine scalar_ve2tr(fve, ftr, mesh)
!    implicit none
!    type(grid_structure), intent(in) :: mesh
!    type(scalar_field), intent(in):: fve
!    type(scalar_field), intent(inout) :: ftr
!    real*8  :: pp(3),arean
!    integer :: i,j,k
!    
!    do i=1,mesh%nt
!      pp = mesh%tr(i)%c%p;
!      do j =1,3
!        k = mesh%tr(i)%v(j);
!!        ed  = mesh%tr(i)%v(modulo(j-1,3));
!!        pu1 =  mesh%ed(ed)%c%p;

!!        ed  = mesh%tr(i)%v(j);
!!        pu2 =  mesh%ed(ed)%c%p;

!!        area = sphtriarea(pp, pu1, pu2);
!        arean = mesh%tr(i)%trhx_areag(j);
!        ftr%f = arean*fve%f(k)
!      enddo
!    enddo
!  
!  end subroutine scalar_ve2tr


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
!        ed_area1=sphtriarea(pp, pu, p1)+sphtriarea(pp, pu, p2);
        ed_area1=sphtriarea(pp, p1, p2);


        pp=mesh%tr(sh(2))%c%p;
!        ed_area2=sphtriarea(pp, pu, p1)+sphtriarea(pp, pu, p2);
        ed_area2=sphtriarea(pp, p1, p2);

!        fed%f(l) = (ed_area1*ftr%f(sh(1))+ed_area2*ftr%f(sh(2)))/(ed_area1+ed_area2);
        fed%f(l) = (ed_area1*ftr%f(mesh%ed(l)%sh(1))+ed_area2*ftr%f(mesh%ed(l)%sh(2)))/(ed_area1+ed_area2);


      endif
      
    end do
    !$omp end parallel do

    return

  end subroutine scalar_tr2ed

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
!        fve%f(k) = fve%f(k) + ftr%f(trp(l));
      enddo
      fve%f(k) = fve%f(k)/mesh%hx(k)%areag;
!      fve%f(k) = fve%f(k)/mesh%v(k)%nnb;
      
    enddo

    return

  end subroutine scalar_tr2ve

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

  subroutine vort_hx(u, zeta,mesh)
    implicit none
    !---------------------------------------------------------------
    !Calculate fields at hexagon
    !   vorticity at hexagon cc based on edge_hx velocities
    !---------------------------------------------------------------
    type(grid_structure), intent(in) :: mesh
    type(scalar_field), intent(in):: u  ! velocity at cell edges
    type(scalar_field), intent(inout):: zeta !absolute vorticity at tr cc

    integer(i4):: k, l, ed,vd
    real*8   :: A,signcor,de,pp(3),pu(3)

    zeta%f=0d0

    !$omp parallel do &
    !$omp default(none) &
    !$omp shared(mesh, u, zeta) &
    !$omp private(k,l,A, ed, vd, signcor) &
    !$omp private(pp,pu,de) &
    !$omp schedule(static)
    do k=1, mesh%nv
      A = mesh%hx(k)%areag*rad**2;

      !loop over triangle edges
      do l=1, mesh%v(k)%nnb
        ed=mesh%v(k)%ed(l);
        vd=mesh%v(k)%tr(l);

        pp = mesh%tr(vd)%c%p;
        pu = mesh%ed(ed)%c%p;
        de = mesh%edhx(ed)%leng*rad;
        signcor=dsign(1d0,dot_product(mesh%ed(ed)%nr,pp-pu));
        zeta%f(k)=zeta%f(k)+u%f(ed)*de*signcor/A;
      end do
    end do
    !$omp end parallel do
    return

  end subroutine vort_hx

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
        ed_area = sphtriarea(pp1,pp2,pq);
        ke%f(k)=ke%f(k)+ed_area*u%f(ed)**2;
      end do
      ke%f(k)=ke%f(k)/A
    end do
    !$omp end parallel do

    return

  end subroutine kinetic_energy_v
!--------------------------------------------

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

  subroutine saveoperators(pathn,glevel,u,B,uhq_perp,grad,ke,div,eta,l,mesh)
    implicit none
    type(grid_structure), intent(in) :: mesh
    type(scalar_field), intent(in):: uhq_perp,grad,ke,div,B,u,eta
    integer, intent(in) :: l
    integer             :: i
    character(len=2),intent(in) :: glevel
    character(len=100),intent(in) :: pathn
    character(len=100) :: savefb,savefv,savefg,savefk,savefup,savefdi,savefzt,savefvt
  
    type(vector_field_cart) :: f_tr

    allocate(f_tr%p(mesh%nt));
    call perot_ed2tr(u, f_tr, mesh);


    if (l<10) then
      write(savefb,'(A,A,A,A,1I1,A)')  trim(pathn),'h/',trim(glevel),'/B000',l,'.dat';
      write(savefv,'(A,A,A,A,1I1,A)')  trim(pathn),'v/',trim(glevel),'/v000',l,'.dat';
      write(savefg,'(A,A,A,A,1I1,A)')  trim(pathn),'grad/',trim(glevel),'/grad000',l,'.dat';
      write(savefk,'(A,A,A,A,1I1,A)')  trim(pathn),'ke/',trim(glevel),'/ke000',l,'.dat';
      write(savefup,'(A,A,A,A,1I1,A)') trim(pathn),'uperp/',trim(glevel),'/up000',l,'.dat';
      write(savefdi,'(A,A,A,A,1I1,A)') trim(pathn),'divv/',trim(glevel),'/div000',l,'.dat';
      write(savefzt,'(A,A,A,A,1I1,A)') trim(pathn),'eta/',trim(glevel),'/eta000',l,'.dat';
      write(savefvt,'(A,A,A,A,1I1,A)') trim(pathn),'vt/',trim(glevel),'/vt000',l,'.dat';

    elseif (l<100) then
      write(savefb,'(A,A,A,A,1I2,A)')  trim(pathn),'h/',trim(glevel),'/B00',l,'.dat';
      write(savefv,'(A,A,A,A,1I2,A)')  trim(pathn),'v/',trim(glevel),'/v00',l,'.dat';
      write(savefg,'(A,A,A,A,1I2,A)')  trim(pathn),'grad/',trim(glevel),'/grad00',l,'.dat';
      write(savefk,'(A,A,A,A,1I2,A)')  trim(pathn),'ke/',trim(glevel),'/ke00',l,'.dat';
      write(savefup,'(A,A,A,A,1I2,A)') trim(pathn),'uperp/',trim(glevel),'/up00',l,'.dat';
      write(savefdi,'(A,A,A,A,1I2,A)') trim(pathn),'divv/',trim(glevel),'/div00',l,'.dat';
      write(savefzt,'(A,A,A,A,1I2,A)') trim(pathn),'eta/',trim(glevel),'/eta00',l,'.dat';
      write(savefvt,'(A,A,A,A,1I2,A)') trim(pathn),'vt/',trim(glevel),'/vt00',l,'.dat';

    elseif (l<1000) then
      write(savefb,'(A,A,A,A,1I3,A)')  trim(pathn),'h/',trim(glevel),'/B0',l,'.dat';
      write(savefv,'(A,A,A,A,1I3,A)')  trim(pathn),'v/',trim(glevel),'/v0',l,'.dat';
      write(savefg,'(A,A,A,A,1I3,A)')  trim(pathn),'grad/',trim(glevel),'/grad0',l,'.dat';
      write(savefk,'(A,A,A,A,1I3,A)')  trim(pathn),'ke/',trim(glevel),'/ke0',l,'.dat';
      write(savefup,'(A,A,A,A,1I3,A)') trim(pathn),'uperp/',trim(glevel),'/up0',l,'.dat';
      write(savefdi,'(A,A,A,A,1I3,A)') trim(pathn),'divv/',trim(glevel),'/div0',l,'.dat';
      write(savefzt,'(A,A,A,A,1I3,A)') trim(pathn),'eta/',trim(glevel),'/eta0',l,'.dat';
      write(savefvt,'(A,A,A,A,1I3,A)') trim(pathn),'vt/',trim(glevel),'/vt0',l,'.dat';
    else
      write(savefb,'(A,A,A,A,1I4,A)')  trim(pathn),'h/',trim(glevel),'/B',l,'.dat';
      write(savefv,'(A,A,A,A,1I4,A)')  trim(pathn),'v/',trim(glevel),'/v',l,'.dat';
      write(savefg,'(A,A,A,A,1I4,A)')  trim(pathn),'grad/',trim(glevel),'/grad',l,'.dat';
      write(savefk,'(A,A,A,A,1I4,A)')  trim(pathn),'ke/',trim(glevel),'/ke',l,'.dat';
      write(savefup,'(A,A,A,A,1I4,A)') trim(pathn),'uperp/',trim(glevel),'/up',l,'.dat';
      write(savefdi,'(A,A,A,A,1I4,A)') trim(pathn),'divv/',trim(glevel),'/div',l,'.dat';
      write(savefzt,'(A,A,A,A,1I4,A)') trim(pathn),'eta/',trim(glevel),'/eta',l,'.dat';
      write(savefvt,'(A,A,A,A,1I4,A)') trim(pathn),'vt/',trim(glevel),'/vt',l,'.dat';
    endif

    write(*,*) '-->',savefb
    write(*,*) '-->',savefv
    write(*,*) '-->',savefg
    write(*,*) '-->',savefk
    write(*,*) '-->',savefup
    write(*,*) '-->',savefdi
    write(*,*) '-->',savefzt
    write(*,*) '-->',savefvt


    
    open(unit=100,file=savefb(2:100));
    open(unit=101,file=savefv(2:100));
    open(unit=102,file=savefg(2:100));
    open(unit=103,file=savefup(2:100));
    open(unit=104,file=savefdi(2:100));
    open(unit=105,file=savefzt(2:100));
    open(unit=106,file=savefk(2:100));
    open(unit=107,file=savefvt(2:100));



    do i=1,mesh%nt
!      f_tr%p(i)%v = proj_vec_sphere(f_tr%p(i)%v,mesh%tr(i)%c%p);

      write(100,*) B%f(i);
      write(104,*) div%f(i);
      write(106,*) ke%f(i);
      write(107,*) f_tr%p(i)%v

    enddo
    do i=1,mesh%ne
      write(101,*) u%f(i);
      write(102,*) grad%f(i);
      write(103,*) uhq_perp%f(i);
    enddo


    do i=1,mesh%nv
      write(105,*) eta%f(i);
    enddo
    do i = 100,107
      close(i);
    enddo
  end subroutine saveoperators

  subroutine inicond(j,u,h,bath,div,grad,ke,zeta,q,uhq_perp,f,fv,mesh,resn,nlin)
    implicit none
    type(grid_structure), intent(in) :: mesh
    type(scalar_field), intent(inout):: u,h,div,zeta,grad,uhq_perp,f,fv,q,bath,ke
    type(scalar_field)               :: eta_ed,h_ed
    type(vector_field_cart)          :: var_tr
    integer,intent(in) :: j,nlin
    integer            :: i,mi,ni
    real*8             :: v1(3),v2(3),pp(3),pu(3),pq(3),lon,lat,r,h0,b0,R0,lonc,latc,hv,var,utheta,uphi,var1,var2,alpha,u0
    character(len=100) ::resn

    bath%f = 10d0

    
    if (j==1) then ! Geostrophic Zonal Flow
      write(resn,*) './result/TC1/snapshots/'
      allocate(eta_ed%fexact(mesh%ne),h_ed%fexact(mesh%ne));

      div%fexact = 0d0;
      f%f  = 10d0;
      fv%f = 10d0;
      do i=1,mesh%nt
        pp = mesh%tr(i)%c%p;
        bath%f(i) = 1d0/pi*cos(pi*(pp(3)+1));
        h%f(i) = 1d0;
        ke%fexact(i) = (pp(1)**2+pp(2)**2)*sin(pi*(pp(3)+1))**2*.5d0;
      enddo
      h%fexact = h%f;

      do i =1,mesh%ne
        pu = mesh%ed(i)%c%p;

        v1(1) = -pu(2)*sin(pi*(pu(3)+1)); ! dudx = 0d0       ; dudy = -sin(theta); dudz = -y*cos(theta)
        v1(2) = pu(1)*sin(pi*(pu(3)+1));  ! dvdx = sin(theta); dvdy = 0d0        ; dvdz = x*cos(theta)
        v1(3) = 0d0;

        v2(1) = mesh%ed(i)%nr(1);
        v2(2) = mesh%ed(i)%nr(2);
        v2(3) = mesh%ed(i)%nr(3);

        u%f(i) = dot_product(v1,v2);
        u%fexact(i) = u%f(i);

!        v1(1) = pu(1)*sin(pi*(pu(3)+1))**2;
!        v1(2) = pu(2)*sin(pi*(pu(3)+1))**2;
!        v1(3) = pi*(pu(1)**2+pu(2)**2)*cos(pi*(pu(3)+1))*sin(pi*(pu(3)+1));
!        gradk%fexact(i) = -dot_product(v1,v2);


        v1(1) = -pu(2)*sin(pi*(pu(3)+1));
        v1(2) = pu(1)*sin(pi*(pu(3)+1));
        v1(3) = 0d0;

        v2(1) = mesh%ed(i)%tg(1);
        v2(2) = mesh%ed(i)%tg(2);
        v2(3) = mesh%ed(i)%tg(3);

        uhq_perp%fexact(i) = dot_product(v1,v2);

        v1(1) = 0d0;
        v1(2) = 0d0;
        v1(3) = -sin(pi*(pu(3)+1));


        v2(1) = mesh%ed(i)%nr(1);
        v2(2) = mesh%ed(i)%nr(2);
        v2(3) = mesh%ed(i)%nr(3);
        grad%fexact(i) = -dot_product(v1,v2);

        h_ed%fexact(i) = 1d0;
        eta_ed%fexact(i) = -pi*(pu(1)**2+pu(2)**2)*cos(pi*(pu(3)+1))+2d0*pq(3)*sin(pi*(pu(3)+1))+10d0;
      enddo

      if (nlin==1) then
        uhq_perp%fexact = uhq_perp%fexact*eta_ed%fexact;
      else
        uhq_perp%fexact = uhq_perp%fexact*10d0;
      endif

      do i=1,mesh%nv
        pq = mesh%v(i)%p;
        zeta%fexact(i) = -pi*(pq(1)**2+pq(2)**2)*cos(pi*(pq(3)+1))+2d0*pq(3)*sin(pi*(pq(3)+1));
        hv = 1d0;
        q%fexact(i) = (zeta%fexact(i)+fv%f(i))/hv;
      enddo

    else if (j==2) then ! mountain
      write(resn,*) './result/TC2/snapshots/'

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
        h%f(i) = h0;

      enddo

      do i =1,mesh%ne
        pu = mesh%ed(i)%c%p;

        lat = datan2(pu(3),dsqrt(pu(1)**2+pu(2)**2));
!        f%f(i) = 86400d0*(2*2*pi/86400d0*sin(lat));
        f%f(i) = 2*omega*dsin(lat);
!        f%f(i) = 10d0;

        v1(1) = -pu(2); ! dudx = 0d0       ; dudy = -sin(theta); dudz = -y*cos(theta)
        v1(2) = pu(1);  ! dvdx = sin(theta); dvdy = 0d0        ; dvdz = x*cos(theta)
        v1(3) = 0d0;

        v2(1) = mesh%ed(i)%nr(1);
        v2(2) = mesh%ed(i)%nr(2);
        v2(3) = mesh%ed(i)%nr(3);

        u%f(i) = dot_product(v1,v2)*2d1;!/norm(v1)*dcos(lat);
      enddo

      do i=1,mesh%nv
        pq = mesh%v(i)%p;
        lat = datan2(pq(3),dsqrt(pq(1)**2+pq(2)**2));
!        fv%f(i) = 86400d0*(2*2*pi/86400d0*sin(lat));
        fv%f(i) = 2*omega*dsin(lat);
      enddo

    else if (j==3) then ! Hollingsworth instability 2
      write(resn,*) './result/TC3/snapshots/'
      f%f = 1d0/mesh%ed(1)%leng;
      fv%f = 1d0/mesh%ed(1)%leng;

!      f%f = 10d0;
!      fv%f = 10d0;

!      h%f = 0d0;
!      h%f = 2.5d0*(1d0/128d0)**2;
!      h%f = 1d1;

      h%f = 1d0;
      do i=1,mesh%nt
        pp = mesh%tr(i)%c%p;
        bath%f(i) = 1d1*f%f(1)/(pi*g)*dcos(pi*(pp(3)+1));

!        lon = datan2(pp(2),pp(1));
!        lat = datan2(pp(3),dsqrt(pp(1)**2+pp(2)**2));
!        if (abs(lon*180/pi)<5 .and. abs(lat*180/pi)<5) then
!          h%f(i) = h%f(i) + h%f(i)/1d3;
!          exit;
!        endif
      enddo

      do i =1,mesh%ne
        pu = mesh%ed(i)%c%p;

        v1(1) = -pu(2);!*sin(pi*(pu(3)+1));
        v1(2) = pu(1);!*sin(pi*(pu(3)+1));
        v1(3) = 0d0;

        v2(1) = mesh%ed(i)%nr(1);
        v2(2) = mesh%ed(i)%nr(2);
        v2(3) = mesh%ed(i)%nr(3);

!        u%f(i) = 50d-1*mesh%ed(1)%leng*dot_product(v1,v2)
        u%f(i) = 1d1*dot_product(v1,v2)*dsin(pi*(pu(3)+1));
      enddo
    else if (j==4) then
      write(resn,*) './result/TC4/snapshots/'

      f%f  = 10d0;
      fv%f = 10d0;
      do i=1,mesh%nt
        pp = mesh%tr(i)%c%p;
        bath%f(i) = 1d0/pi*cos(pi*(pp(3)+1));
        h%f(i) = 1d0;
      enddo
      h%fexact = h%f;

      do i =1,mesh%ne
        pu = mesh%ed(i)%c%p;

        v1(1) = -pu(2)*sin(pi*(pu(3)+1)); ! dudx = 0d0       ; dudy = -sin(theta); dudz = -y*cos(theta)
        v1(2) = pu(1)*sin(pi*(pu(3)+1));  ! dvdx = sin(theta); dvdy = 0d0        ; dvdz = x*cos(theta)
        v1(3) = 0d0;

        v2(1) = mesh%ed(i)%nr(1);
        v2(2) = mesh%ed(i)%nr(2);
        v2(3) = mesh%ed(i)%nr(3);

        u%f(i) = 10d0*dot_product(v1,v2);
      enddo

    else if (j==5) then
      write(resn,*) './result/TC5/snapshots/';

      mi = 1;
      ni = 1;

      f%f  = 10d0;
      fv%f = 10d0;

      do i=1,mesh%nt
        pp = mesh%tr(i)%c%p;
        lon = datan2(pp(2),pp(1)); ! theta
        lat = datan2(pp(3),dsqrt(pp(1)**2+pp(2)**2)); ! phi

        var1   = -4d0*ni**2*(dcos(ni*lat)**4-3d0*dcos(ni*lat)**2*dsin(ni*lat)**2)*dsin(lon)*dcos(mi*lon);
        var2 = -mi*dcos(ni*lat)**4/dcos(lat)*(mi*dsin(lon)*dcos(mi*lon)+dcos(lon)*dsin(mi*lon));
        uphi = -4d0*ni*dcos(ni*lat)**3*dsin(ni*lat)*dsin(lon)*dcos(mi*lon);


        bath%f(i) = 0d0;
        h%f(i) = 1d0;

        div%fexact(i) = h%f(i)/rad*(var1 - (uphi*dtan(lat) + var2/dcos(lat)));

      enddo

      h%fexact = h%f;

      do i =1,mesh%ne
        pu = mesh%ed(i)%c%p;

        lon = datan2(pu(2),pu(1));
        lat = datan2(pu(3),dsqrt(pu(1)**2+pu(2)**2));

!        utheta = mi*dcos(ni*lat)**4/dcos(lat)*dsin(lon)*dsin(mi*lon);!*dcos(lat);
!        uphi = -4d0*ni*dcos(ni*lat)**3*dsin(ni*lat)*dsin(lon)*dcos(mi*lon);!*dcos(lat);

        utheta = dcos(lat)**3*dsin(lon)**2;
        uphi = -4d0*dcos(lat)**3*dsin(lat)*dsin(lon)*dcos(lon);

        call convert_vec_sph2cart(utheta, uphi, pu, v1)


        v2(1) = mesh%ed(i)%nr(1);
        v2(2) = mesh%ed(i)%nr(2);
        v2(3) = mesh%ed(i)%nr(3);

        u%f(i) = dot_product(v1,v2);


        if (nlin==0) then
          uhq_perp%fexact(i) = Hm*dot_product(v1,mesh%ed(i)%tg);
        endif
      enddo
      grad%fexact = 0d0;

    else if (j==6) then
      write(resn,*) './result/TC6/snapshots/';

      f%f  = omega;
      fv%f = omega;

      lonc = 0d0;
      latc = 0d0;
      R0 = .3d0;
      do i=1,mesh%nt
        pp = mesh%tr(i)%c%p;

        lon = datan2(pp(2),pp(1));
        lat = datan2(pp(3),dsqrt(pp(1)**2+pp(2)**2));


        bath%f(i) = exp(-(((lon-lonc)/(2.5*R0))**2+((lat-latc)/R0)**2));
      enddo
      h%f = 1d3;
      div%fexact = 0d0;
!-------------------------------------
      do i =1,mesh%ne
        pu = mesh%ed(i)%c%p;

        lon = datan2(pu(2),pu(1));
        lat = datan2(pu(3),dsqrt(pu(1)**2+pu(2)**2));

        var = -(((lon-lonc)/(2.5*R0))**2+((lat-latc)/R0)**2)
!        uphi = -1d0/(f%f(1)*rad)*(lon-lonc)/(2.5*R0)**2*exp(var);
!        utheta = 1d0/(f%f(1)*rad*dcos(lat))*(lat-latc)/(R0)**2*exp(var);

        uphi = 2d0*g/(f%f(1)*rad)*lat/(R0)**2*exp(var);
        utheta = -2d0*g/(f%f(1)*rad*dcos(lat))*lon/(2.5*R0)**2*exp(var);


        call convert_vec_sph2cart(uphi,utheta, pu, v1)

        v2(1) = mesh%ed(i)%nr(1);
        v2(2) = mesh%ed(i)%nr(2);
        v2(3) = mesh%ed(i)%nr(3);

        u%f(i) = dot_product(v1,v2);
      enddo
!-------------------------------------
    else if (j==7) then ! Barotropically Unstable Jet
      write(resn,*) './result/TC7/snapshots/';

!      do i =1,mesh%nt
!      enddo

      bath%f = 10d0;
      h%f    = 10d0;

!      do i=1,mesh%nt
!        pp = mesh%ed(i)%c%p;
!        lat = datan2(pp(3),dsqrt(pp(1)**2+pp(2)**2)); ! phi
!        h%f(i)    = 80d0/exp(-4d0/(5d0/14d0*pi-pi/7d0)**2)*exp(1d0/((lat-pi/7d0)*(lat-5d0/14d0*pi)))*(pi/7d0+5d0/14d0*pi)/((lat-pi/7d0)**2*(lat-5d0/14d0*pi)**2);
!      enddo

      u%f = 0d0;
      do i=1,mesh%ne
        pu = mesh%ed(i)%c%p;
        lon = datan2(pu(2),pu(1)); ! theta
        lat = datan2(pu(3),dsqrt(pu(1)**2+pu(2)**2)); ! phi
        f%f(i)  = 2*omega*dsin(lat);

        v1(1) = -pu(2);!*sin(pi*(pu(3)+1));
        v1(2) = pu(1);!*sin(pi*(pu(3)+1));
        v1(3) = 0d0;

        v2(1) = mesh%ed(i)%nr(1);
        v2(2) = mesh%ed(i)%nr(2);
        v2(3) = mesh%ed(i)%nr(3);

        if (pi/7d0<=lat .and. lat <= 5d0/14d0*pi) then
          u%f(i) = 80d0/exp(-4d0/(5d0/14d0*pi-pi/7d0)**2)*exp(1d0/((lat-pi/7d0)*(lat-5d0/14d0*pi)))*dot_product(v1,v2)/norm(v1);
        endif
      enddo

      do i=1,mesh%nv
        pq = mesh%v(i)%p;
        lon = datan2(pq(2),pq(1)); ! theta
        lat = datan2(pq(3),dsqrt(pq(1)**2+pq(2)**2)); ! phi
        fv%f(i)  = 2*omega*dsin(lat);
      enddo

    else if (j==8) then ! Nonlinear Geostrophy
      write(resn,*) './result/TC8/snapshots/';
      alpha = 0;
      u0 = 2d0*pi*rad/(86400d0*12);
      b0 = 2.94d4/g;
      allocate(eta_ed%fexact(mesh%ne),h_ed%fexact(mesh%ne),var_tr%p(mesh%nt));

      div%fexact = 0d0;
      do i=1,mesh%nt
        pp = mesh%tr(i)%c%p;
        lat = datan2(pp(3),dsqrt(pp(1)**2+pp(2)**2));
        lon = datan2(pp(2),pp(1));

        bath%f(i) = 0d0;
!        bath%f(i) = b0 - 1d0/g*(rad*omega*u0+u0**2*.5d0)*(dcos(lon)*dcos(lat)*dsin(alpha)+dsin(lat)*dcos(alpha))**2;
        h%f(i) = b0 - 1d0/g*(rad*omega*u0+u0**2*.5d0)*(dcos(lon)*dcos(lat)*dsin(alpha)+dsin(lat)*dcos(alpha))**2;

        uphi   =  u0*(dcos(lat)*dcos(alpha)+dcos(lon)*dsin(lat)*dsin(alpha));
        utheta = -u0*dsin(lon)*dsin(alpha);

        ke%fexact(i) = (uphi**2+utheta**2)*.5d0;

        call convert_vec_sph2cart(uphi, utheta, pp, v1)
        var_tr%p(i)%v = v1;
      enddo
      call perot_tr2ed(var_tr,u, mesh);
      h%fexact = h%f;

      do i =1,mesh%ne
        pu = mesh%ed(i)%c%p;
        lat = datan2(pu(3),dsqrt(pu(1)**2+pu(2)**2));
        lon = datan2(pu(2),pu(1));

        f%f(i)  = 2d0*omega*(-dcos(lon)*dcos(lat)*dsin(alpha)+dsin(lat)*dcos(alpha));

        uphi   =  u0*(dcos(lat)*dcos(alpha)+dcos(lon)*dsin(lat)*dsin(alpha));
        utheta = -u0*dsin(lon)*dsin(alpha);

        call convert_vec_sph2cart(uphi, utheta, pu, v1)

        v2(1) = mesh%ed(i)%nr(1);
        v2(2) = mesh%ed(i)%nr(2);
        v2(3) = mesh%ed(i)%nr(3);

!        u%f(i) = dot_product(v1,v2);
        u%fexact(i) = u%f(i);

        uphi   =  u0*(dcos(lat)*dcos(alpha)+dcos(lon)*dsin(lat)*dsin(alpha));
        utheta = -u0*dsin(lon)*dsin(alpha);

        call convert_vec_sph2cart(-utheta, uphi, pu, v1)
!        v2(1) = mesh%ed(i)%tg(1);
!        v2(2) = mesh%ed(i)%tg(2);
!        v2(3) = mesh%ed(i)%tg(3);

        uhq_perp%fexact(i) = dot_product(v1,v2);

        var    = -2d0/g*(rad*omega*u0+u0**2*.5d0)*(dsin(lat)*dcos(alpha)-dcos(lon)*dcos(lat)*dsin(alpha))
        utheta = var/rad*(dcos(lon)*dsin(lat)*dsin(alpha)+dcos(lat)*dcos(alpha));
        uphi   = var/(rad*dcos(lat))*dsin(lon)*dcos(lat)*dsin(alpha);
        call convert_vec_sph2cart(uphi, utheta, pu, v1)

        v2(1) = mesh%ed(i)%nr(1);
        v2(2) = mesh%ed(i)%nr(2);
        v2(3) = mesh%ed(i)%nr(3);


        grad%fexact(i) = -dot_product(v1,v2)*g;

        h_ed%fexact(i) = b0 - 1d0/g*(rad*omega*u0+u0**2*.5d0)*(dcos(lon)*dcos(lat)*dsin(alpha)+dsin(lat)*dcos(alpha))**2;
!        h_ed%fexact(i) = 1d0;

        uphi = u0*(dcos(lat)*dcos(alpha)+dcos(lon)*dsin(lat)*dsin(alpha));
        var1 = u0*(dcos(lon)*dcos(lat)*dsin(alpha)-dsin(lat)*dcos(alpha));
        var2 = -u0*dcos(lon)*dsin(alpha);
        var = -1d0/(rad*dcos(lat))*(dcos(lat)*var1-dsin(lat)*uphi-var2);
        eta_ed%fexact(i) = (var+f%f(i));
      enddo

      if (nlin==1) then
        uhq_perp%fexact = uhq_perp%fexact*eta_ed%fexact;
      else
        uhq_perp%fexact = uhq_perp%fexact*f%f;
      endif

      do i=1,mesh%nv
        pq = mesh%v(i)%p;
        lat = datan2(pq(3),dsqrt(pq(1)**2+pq(2)**2));
        lon = datan2(pq(2),pq(1));

        fv%f(i)  = 2d0*omega*(-dcos(lon)*dcos(lat)*dsin(alpha)+dsin(lat)*dcos(alpha));

        uphi = u0*(dcos(lat)*dcos(alpha)+dcos(lon)*dsin(lat)*dsin(alpha));
        var1 = u0*(dcos(lon)*dcos(lat)*dsin(alpha)-dsin(lat)*dcos(alpha));
        var2 = -u0*dcos(lon)*dsin(alpha);
        zeta%fexact(i) = -1d0/(rad*dcos(lat))*(dcos(lat)*var1-dsin(lat)*uphi-var2);

        hv = b0 - 1d0/g*(rad*omega*u0+u0**2*.5d0)*(dcos(lon)*dcos(lat)*dsin(alpha)+dsin(lat)*dcos(alpha))**2;
!        hv = 1d0;
        q%fexact(i) = (zeta%fexact(i)+fv%f(i))/hv;
      enddo
    else if (j==9) then ! Nonlinear Geostrophy with constant h
      write(resn,*) './result/TC9/snapshots/';
      u0 = 2d0*pi*rad/(86400d0*12);
      b0 = 2.94d4/g;
      allocate(eta_ed%fexact(mesh%ne),h_ed%fexact(mesh%ne),var_tr%p(mesh%nt));

      div%fexact = 0d0;
      do i=1,mesh%nt
        pp = mesh%tr(i)%c%p;
        lat = datan2(pp(3),dsqrt(pp(1)**2+pp(2)**2));
        lon = datan2(pp(2),pp(1));

        bath%f(i) = b0 - 1d0/g*(rad*omega*u0+u0**2*.5d0)*dsin(lat)**2;
        h%f(i) = 1d0;

        uphi   =  u0*dcos(lat);
        utheta = 0d0;

        ke%fexact(i) = (uphi**2+utheta**2)*.5d0;

        call convert_vec_sph2cart(uphi, utheta, pp, v1)
        var_tr%p(i)%v = v1;
      enddo
      call perot_tr2ed(var_tr,u, mesh);
      h%fexact = h%f;

      do i =1,mesh%ne
        pu = mesh%ed(i)%c%p;
        lat = datan2(pu(3),dsqrt(pu(1)**2+pu(2)**2));
        lon = datan2(pu(2),pu(1));

        f%f(i)  = 2d0*omega*dsin(lat);

        uphi   = u0*dcos(lat);
        utheta = 0d0;

        call convert_vec_sph2cart(uphi, utheta, pu, v1)

        v2(1) = mesh%ed(i)%nr(1);
        v2(2) = mesh%ed(i)%nr(2);
        v2(3) = mesh%ed(i)%nr(3);

!        u%f(i) = dot_product(v1,v2);
        u%fexact(i) = u%f(i);

        uphi   =  u0*dcos(lat);
        utheta = 0d0;

        call convert_vec_sph2cart(-utheta, uphi, pu, v1)

        uhq_perp%fexact(i) = dot_product(v1,v2);


        v2(1) = mesh%ed(i)%nr(1);
        v2(2) = mesh%ed(i)%nr(2);
        v2(3) = mesh%ed(i)%nr(3);

        var    = -2d0/g*(rad*omega*u0+u0**2*.5d0)*dsin(lat)
        utheta = var/rad*dcos(lat);
        uphi   = 0d0;
        call convert_vec_sph2cart(uphi, utheta, pu, v1)

        v2(1) = mesh%ed(i)%nr(1);
        v2(2) = mesh%ed(i)%nr(2);
        v2(3) = mesh%ed(i)%nr(3);


        grad%fexact(i) = -g*dot_product(v1,v2);

        h_ed%fexact(i) = 1d0;

        uphi = u0*dcos(lat);
        var1 = u0*(-dsin(lat));
        var2 = 0d0;
        var = -1d0/(rad*dcos(lat))*(dcos(lat)*var1-dsin(lat)*uphi-var2);
        eta_ed%fexact(i) = (var+f%f(i));
      enddo

      if (nlin==1) then
        uhq_perp%fexact = uhq_perp%fexact*eta_ed%fexact;
      else
        uhq_perp%fexact = uhq_perp%fexact*f%f;
      endif

      do i=1,mesh%nv
        pq = mesh%v(i)%p;
        lat = datan2(pq(3),dsqrt(pq(1)**2+pq(2)**2));
        lon = datan2(pq(2),pq(1));

        fv%f(i)  = 2d0*omega*dsin(lat);

        uphi = u0*dcos(lat);
        var1 = u0*(-dsin(lat));
        var2 = 0d0;
        zeta%fexact(i) = -1d0/(rad*dcos(lat))*(dcos(lat)*var1-dsin(lat)*uphi-var2);

        hv = 1d0;
        q%fexact(i) = (zeta%fexact(i)+fv%f(i))/hv;
      enddo


    else
      write(resn,*) './result/snapshots/'

      f%f  = 10d0;
      fv%f = 10d0;

      call RANDOM_NUMBER(bath%f);
      call RANDOM_NUMBER(h%f);
      h%f=10d0*(h%f+bath%f);

      call RANDOM_NUMBER(u%f)

    endif

    return
  
  end subroutine inicond

  subroutine calcop(u,h,uh,h_ed,phi,zeta,q,div,uhq_perp,gradb,gradk,ke,f,fv,mesh,nlin,dt)
    implicit none
    type(grid_structure), intent(in) :: mesh
    type(scalar_field), intent(in)   :: f,fv,u,h,phi
    type(scalar_field), intent(inout):: div,gradb,uhq_perp,h_ed,q,zeta,uh,gradk
    type(scalar_field)               :: ket,kev,ketv,ke,h_ve,q_ed,eta
    type(vector_field_cart)          :: var_tr,uhq_v,varv
    real*8                           :: alpha
    real*8,intent(in) :: dt
    integer,intent(in) :: nlin
    integer :: i

    allocate(h_ve%f(mesh%nv),var_tr%p(mesh%nt), &
             uhq_v%p(mesh%nv),eta%f(mesh%nv), &
             q_ed%f(mesh%ne),varv%p(mesh%nt));

    allocate(ket%f(mesh%nt),kev%f(mesh%nv),ketv%f(mesh%nt));

!------------------------------------------------------------------------
    ke%f     = 0d0;
    q_ed%f   = f%f;
    q%f      = fv%f;
    gradk%f  = 0d0;
    h_ed%f   = Hm;
    eta%f    = fv%f;

    !-Curl u-----------------------------
    call vort_hx(u, zeta, mesh); ! Curl
    !------------------------------------

    !-grad Phi---------------------------
    call grad_ed(phi, gradb, mesh);
    !------------------------------------


    !-Kinetic energy--------------------
    ket%f  = 0d0;
    kev%f  = 0d0;
    ketv%f = 0d0;

    alpha = 1d0;
    if (.true.) then ! kinetic energy
      call kinetic_energy_tr(u, ket, mesh);
!      call kinetic_energy_v(u, kev, mesh);
!      do i=1,mesh%nt
!        ketv%f(i) = sum(kev%f(mesh%tr(i)%v(:))*mesh%tr(i)%trhx_areag(:))/mesh%tr(i)%areag;
!      enddo      
    else             ! Korn, 2017
      call perot_ed2tr(u, var_tr, mesh);
      do i =1,mesh%nt
        ket%f(i) = norm(var_tr%p(i)%v);
        ket%f(i) = ket%f(i)**2*.5d0;
      enddo
    endif
    ke%f = alpha*ket%f+(1d0-alpha)*ketv%f;
    !------------------------------------

    call scalar_tr2ve(h, h_ve, mesh)
    q%f = (zeta%f+fv%f)/h_ve%f;

    !-Non linear-------------------------
    if (nlin==1) then
      eta%f = zeta%f+fv%f;

      call scalar_tr2ed(h, h_ed, mesh);

      !-grad Kin---------------------------
      call grad_ed(ke, gradk, mesh); ! grad Kin
      !------------------------------------

      if (.true.) then 
        call scalar_ve2ed(q, q_ed, mesh); ! arith mean
      else
        call apvm(u, q,dt, q_ed, mesh) ! antecipated pot vort
      endif
    endif

    !------------------------------------

    !-u*h--------------------------------
    if (.true.) then ! weller u*h
      uh%f = u%f*h_ed%f;
    elseif (.false.) then
      uh%f = u%f*h_ed%f;
      call perot_ed2tr(uh, var_tr, mesh);
      call perot_tr2ed(var_tr,uh, mesh);
    else             ! korn u*h
      call perot_ed2tr(u, var_tr, mesh);
      do i =1,mesh%nt
        if (nlin==1) then
          varv%p(i)%v(:) = var_tr%p(i)%v(:)*h%f(i);
        else
          varv%p(i)%v = var_tr%p(i)%v*Hm;
        endif
      enddo
      call perot_tr2ed(varv,uh, mesh);
    endif

    !------------------------------------

    !-Div--------------------------------
    call div_tr(uh, div, mesh);
    !------------------------------------


    !-Coriolis Term----------------------
    if (.true.) then           ! TRSKW
      call coriolis_ed(uh, q_ed, uhq_perp, mesh); ! Nonlinear
    else if (.false.) then      ! TRSKW-dual
      call coriolis_edhx(uh, q_ed,q, uhq_perp, nlin, mesh); ! Nonlinear
    else if (.true.) then      ! Perot
      call perot_ed2v(u, uhq_v,mesh)
      do i =1,mesh%nv
        uhq_v%p(i)%v = eta%f(i)*uhq_v%p(i)%v;
      enddo
      call perot_v2ed(uhq_v, uhq_perp, mesh);
    else
      call perot_ed2v(uh, uhq_v,mesh)
      do i =1,mesh%nv
        uhq_v%p(i)%v = q%f(i)*uhq_v%p(i)%v;
      enddo
      call perot_v2ed(uhq_v, uhq_perp, mesh);
    endif
    !------------------------------------

  end subroutine calcop

  subroutine ode_rk4(u,h,u_new,h_new,bath,momeq,maseq,f,fv,Forcmo,Forcma,mesh,dt,nlin)
    implicit none

    type(grid_structure), intent(in)    :: mesh
    type(scalar_field),   intent(in)    :: f,fv,bath,Forcmo,Forcma,u,h
    type(scalar_field),   intent(inout) :: u_new,h_new
    real*8,               intent(inout) :: momeq(mesh%ne),maseq(mesh%nt)
    real*8,               intent(in)    :: dt                                           
    integer,              intent(in)    :: nlin
    type(scalar_field)                  :: gradk,ke,uh,div,q,zeta,gradb,uhq_perp,h_ed,phi
    real*8                              :: momeq0(mesh%ne),maseq0(mesh%nt), &
                                           momeq1(mesh%ne),maseq1(mesh%nt), &
                                           momeq2(mesh%ne),maseq2(mesh%nt), &
                                           momeq3(mesh%ne),maseq3(mesh%nt)

    allocate(ke%f(mesh%nt),phi%f(mesh%nt),div%f(mesh%nt), &
             gradb%f(mesh%ne),uh%f(mesh%ne),gradk%f(mesh%ne),uhq_perp%f(mesh%ne),h_ed%f(mesh%ne), &
             q%f(mesh%nv),zeta%f(mesh%nv));

    u_new%f=u%f
    h_new%f=h%f


    phi%f = g*(bath%f+h%f);
    call calcop(u,h,uh,h_ed,phi,zeta,q,div,uhq_perp,gradb,gradk,ke,f,fv,mesh,nlin,dt);
    momeq0 = -(uhq_perp%f-gradb%f-gradk%f);
    maseq0 = -div%f;
    
    
    !First RK step------------------------------------
    u_new%f = u%f + dt * momeq0/ 2d0;
    h_new%f = h%f + dt * maseq0/ 2d0;

    phi%f = g*(bath%f+h_new%f);
    call calcop(u_new,h_new,uh,h_ed,phi,zeta,q,div,uhq_perp,gradb,gradk,ke,f,fv,mesh,nlin,dt);
    momeq1 = -(uhq_perp%f-gradb%f-gradk%f);
    maseq1 = -div%f;
    !-------------------------------------------------

    !Second RK step-----------------------------------
    u_new%f = u%f + dt * momeq1/ 2d0;
    h_new%f = h%f + dt * maseq1/ 2d0;

    phi%f = g*(bath%f+h_new%f);
    call calcop(u_new,h_new,uh,h_ed,phi,zeta,q,div,uhq_perp,gradb,gradk,ke,f,fv,mesh,nlin,dt);
    momeq2 = -(uhq_perp%f-gradb%f-gradk%f);
    maseq2 = -div%f;
    !-------------------------------------------------

    !Third  RK step-----------------------------------
    u_new%f = u%f + dt * momeq2;
    h_new%f = h%f + dt * maseq2;

    phi%f = g*(bath%f+h_new%f);
    call calcop(u_new,h_new,uh,h_ed,phi,zeta,q,div,uhq_perp,gradb,gradk,ke,f,fv,mesh,nlin,dt);
    momeq3 = -(uhq_perp%f-gradb%f-gradk%f);
    maseq3 = -div%f;
    !-------------------------------------------------


    momeq = (momeq0 +2d0*momeq1 +2d0*momeq2 + momeq3)/6d0;
    maseq = (maseq0 +2d0*maseq1 +2d0*maseq2 + maseq3)/6d0;
    u_new%f = u%f + dt*(momeq+Forcmo%f);
    h_new%f = h%f + dt*(maseq+Forcma%f);

    phi%f = g*(bath%f+h_new%f);
    call calcop(u_new,h_new,uh,h_ed,phi,zeta,q,div,uhq_perp,gradb,gradk,ke,f,fv,mesh,nlin,dt);

    return
  end subroutine ode_rk4

  subroutine leapfrog(uo,bo,up,bp,un,bn,bath,dt,fv,Forcmo,Forcma,mesh)
    implicit none
    type(grid_structure), intent(in)    :: mesh
    type(scalar_field),   intent(in)    :: fv,uo,Bo,up,Bp,bath,Forcmo,Forcma
    real*8,intent(in)                   :: dt
    type(scalar_field),   intent(inout) :: un,Bn

    type(scalar_field)                  :: phi,div,ke,q,zeta,B_ve,gradb,gradk,uhq_perp,B_ed,uh,q_ed
    real*8                              :: momeq(mesh%ne),maseq(mesh%nt)

    allocate(phi%f(mesh%nt),div%f(mesh%nt),ke%f(mesh%nt),q%f(mesh%nv),zeta%f(mesh%nv),B_ve%f(mesh%nv), &
             gradb%f(mesh%ne),gradk%f(mesh%ne),uhq_perp%f(mesh%ne),B_ed%f(mesh%ne),uh%f(mesh%ne),q_ed%f(mesh%ne));


    phi%f = g*(Bp%f+bath%f);
    call grad_ed(phi, gradb, mesh);

    call scalar_tr2ve(Bp, B_ve, mesh)
    call scalar_tr2ed(Bp, B_ed, mesh);

    q%f = (zeta%f+fv%f)/B_ve%f;

    call scalar_ve2ed(q, q_ed, mesh); ! arith mean

    uh%f = up%f*B_ed%f;
    call coriolis_ed(uh, q_ed, uhq_perp, mesh);
    call div_tr(uh, div, mesh)

    call kinetic_energy_tr(up, ke, mesh);
    call grad_ed(ke, gradk, mesh);

    momeq = uhq_perp%f-gradb%f-gradk%f;
    maseq = div%f;

    Bn%f = Bo%f - 2d0*dt*(maseq+Forcma%f);
    un%f = uo%f - 2d0*dt*(momeq+Forcmo%f);

!    Bo%f = Bp%f;
!    uo%f = up%f;
!    Bp%f = Bn%f;
!    up%f = un%f;

  end subroutine leapfrog

  subroutine errorop(u,h,uh,h_ed,phi,zeta,q,div,uhq_perp,gradb,ke,f,fv,mesh,nlin,dt,glevel)
    implicit none
    type(grid_structure), intent(in)    :: mesh
    type(scalar_field),   intent(in)    :: f,fv,u,h,phi,ke
    type(scalar_field)                  :: gradk
    type(scalar_field),   intent(inout) :: div,zeta,gradb,uhq_perp,h_ed,q,uh
    real*8,intent(in)                   :: dt
    real*8                              :: A(mesh%nt),Ae(mesh%ne),Ax(mesh%nv)
    real*8                              :: cen,dTdt,up2l,gr2l,et2l,dv2l,qe2l,grk2l
    integer,intent(in)                  :: nlin
    integer                             :: i
    character(len=100)                  :: fileng
    character(len=2),intent(in)         :: glevel

    type(scalar_field)                  :: momeq
    real*8                              :: mom2l

    allocate(gradk%f(mesh%ne),momeq%f(mesh%ne),momeq%fexact(mesh%ne))

    A  = sum(mesh%tr(:)%areag);
    Ae = sum(mesh%ed(:)%leng*mesh%edhx(:)%leng);
    Ax = sum(mesh%hx(:)%areag);

    call calcop(u,h,uh,h_ed,phi,zeta,q,div,uhq_perp,gradb,gradk,ke,f,fv,mesh,nlin,dt);
    call tenvar(dTdt,uh,gradb,div,phi,mesh);
    call calc_Corenergies(Cen,uh,uhq_perp,mesh)

!    phi%f = g*(bath%f+h%f);
    momeq%f = -(uhq_perp%f-gradb%f-gradk%f);
    momeq%fexact = 0d0;


    up2l  = l2er(uhq_perp,Ae)
    gr2l  = l2er(gradb,Ae)
    grk2l = l2er(ke,A)
    et2l  = l2er(zeta,Ax)
    dv2l  = l2er(div,A);
    qe2l  = l2er(q,Ax);

    mom2l  = l2er(momeq,Ae)

    write(*,'(A)') '----------------------------------------------------------------------------------------------------';
    write(*,'(A)') ' ';
    write(*,'(A,F21.16,F21.16,F21.16,F21.16)') 'uhq:  ',maxval(abs(uhq_perp%f)),maxval(abs(uhq_perp%fexact)), &
                                                        maxval(abs(uhq_perp%f-uhq_perp%fexact)),up2l;
    write(*,'(A,F21.16,F21.16,F21.16,F21.16)') 'grad: ',maxval(abs((gradb%f))),maxval(abs(gradb%fexact)), &
                                                        maxval(abs((gradb%f)-gradb%fexact)),gr2l;
    write(*,'(A,F21.16,F21.16,F21.16,F21.16)') 'ke:   ',maxval(abs((ke%f))),maxval(abs(ke%fexact)), &
                                                         maxval(abs((ke%f)-ke%fexact)),grk2l;
    write(*,'(A,F21.16,F21.16,F21.16,F21.16)') 'div:  ',maxval(abs(div%f)),maxval(abs(div%fexact)),maxval(abs(div%f-div%fexact)), &
                                                        dv2l;
    write(*,'(A,4F21.16)') 'zeta: ',maxval(abs(zeta%f)),maxval(abs(zeta%fexact)),maxval(abs(zeta%f-zeta%fexact)),et2l;
    write(*,'(A,4F21.16)') 'q:    ',maxval(abs(q%f)),maxval(abs(q%fexact)),maxval(abs(q%f-q%fexact)),qe2l;
    write(*,'(A,4F21.16)') 'momeq:',maxval(abs(momeq%f)),maxval(abs(momeq%fexact)),maxval(abs(momeq%f-momeq%fexact)),mom2l;
    write(*,'(A)') '----------------------------------------------------------------------------------------------------';
    write(*,'(A)') '-> Energy';
    write(*,'(A,E11.5)') 'ce:   ',abs(cen);
    write(*,'(A,E11.5)') 'dTdt: ',abs(dTdt);
    write(*,'(A)') '----------------------------------------------------------------------------------------------------';

    if (.true.) then
      write(fileng,'(A,A,A)') './result/error/',glevel,'/uhq_perp.dat';
      open(unit = 101, file = fileng);

      write(fileng,'(A,A,A)') './result/error/',glevel,'/gradb.dat';
      open(unit = 102, file = fileng);

      write(fileng,'(A,A,A)') './result/error/',glevel,'/ke.dat';
      open(unit = 103, file = fileng);

      write(fileng,'(A,A,A)') './result/error/',glevel,'/div.dat';
      open(unit = 104, file = fileng);

      write(fileng,'(A,A,A)') './result/error/',glevel,'/zeta.dat';
      open(unit = 105, file = fileng);

      write(fileng,'(A,A,A)') './result/error/',glevel,'/q.dat';
      open(unit = 106, file = fileng);

      write(fileng,'(A,A,A)') './result/error/',glevel,'/momeq.dat';
      open(unit = 107, file = fileng);

      do i=1,mesh%nt
        write(103,*) ke%f(i),ke%fexact(i)
        write(104,*) div%f(i),div%fexact(i)
      enddo
      do i=1,mesh%ne
        write(101,*) uhq_perp%f(i),uhq_perp%fexact(i)
        write(102,*) gradb%f(i),gradb%fexact(i)
        write(107,*) momeq%f(i),momeq%fexact(i)
      enddo
      do i=1,mesh%nv
        write(105,*) zeta%f(i),zeta%fexact(i)
        write(106,*) q%f(i),q%fexact(i)
      enddo

      do i=101,107
        close(i);
      enddo
    endif
  end subroutine errorop

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
    
  end function
  
  subroutine diffusion(u_ed,zeta_hx, lplace, mesh)
    implicit none
    type(grid_structure), intent(in)       :: mesh
    type(scalar_field), intent(in)         :: u_ed,zeta_hx
    type(scalar_field)                     :: div,gdiv, gtzeta_hx
    type(scalar_field), intent(inout)      :: lplace
    type(vector_field_cart)                :: f_tr

    allocate(div%f(mesh%nt),gdiv%f(mesh%ne), gtzeta_hx%f(mesh%ne));
!    allocate(lplace%f(mesh%ne));
    allocate(f_tr%p(mesh%nt));


    call div_tr(u_ed, div, mesh);
    call grad_ed(div, gdiv, mesh);


    call grad_edhx(zeta_hx, gtzeta_hx, mesh)

    lplace%f = gdiv%f-gtzeta_hx%f;

  end subroutine diffusion
end module swm_operators

