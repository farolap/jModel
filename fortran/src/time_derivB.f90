module time_derivB
  use basic_funcs
  use constants
  use datastruct, only: &
    grid_structure, &
    vector_field_cart, &
    scalar_field
  use swmb_operators, only: &
    calcop, &
    scalar_ve2tr

  contains

  subroutine ode_rk4(u,h,u_new,h_new,uh,bath,gradb,div,zeta, &
          momeq,maseq,f,fv,Forcmo,Forcma,mesh,dt,tau,time,j)
    implicit none

    type(grid_structure), intent(in)       :: mesh
    type(scalar_field),   intent(in)       :: f,fv,bath,h,Forcma
    type(scalar_field),   intent(inout)    :: h_new,div,zeta
    type(vector_field_cart), intent(in)    :: Forcmo
    type(vector_field_cart), intent(inout) :: uh,u,u_new,gradb
    real*8,               intent(inout)    :: momeq(mesh%nt,3),maseq(mesh%nv)
    real*8,               intent(in)       :: dt,tau,time
    integer,intent(in) :: j
    type(vector_field_cart)                :: veco,vecn,uhq_perp,momflux,lap,bih
    type(scalar_field)                     :: htr,ke,bihm
    real*8                                 :: momeq0(mesh%nt,3),maseq0(mesh%nv), &
      momeq1(mesh%nt,3),maseq1(mesh%nv), &
      momeq2(mesh%nt,3),maseq2(mesh%nv), &
      momeq3(mesh%nt,3),maseq3(mesh%nv),ui,vi
    integer :: i

    allocate(uhq_perp%p(mesh%nt),momflux%p(mesh%nt), &
      ke%f(mesh%nv), lap%p(mesh%nt),bih%p(mesh%nt), &
      htr%f(mesh%nt), bihm%f(mesh%nv), &
      veco%p(mesh%nt),vecn%p(mesh%nt));
    !--------------------------------------------------------------------

    h_new%f=h%f;

    !First RK step------------------------------------
    call calcop(u,h,htr,uh,bath,uhq_perp,zeta,div,gradb,ke, &
        f,fv,momflux,tau,lap,bih,bihm,mesh)
    if (vectinv) then
      do i=1,3
        veco%p(:)%v(i)=u%p(:)%v(i);
      enddo
    else
      do i=1,3
        veco%p(:)%v(i) = uh%p(:)%v(i);
      enddo
    endif

    do i=1,3
!      momeq0(:,i) = -(uhq_perp%p(:)%v(i) &
!          +gradb%p(:)%v(i)+momflux%p(:)%v(i))+bih%p(:)%v(i);
      momeq0(:,i) = -sumreal(uhq_perp%p(:)%v(i),gradb%p(:)%v(i))
      momeq0(:,i) = sumreal(momeq0(:,i),-momflux%p(:)%v(i))
      momeq0(:,i) = sumreal(momeq0(:,i),bih%p(:)%v(i))
    enddo
    maseq0 = -div%f;


    do i=1,3
!      vecn%p(:)%v(i) = veco%p(:)%v(i) + dt * momeq0(:,i)/2d0;
      vecn%p(:)%v(i) = sumreal(veco%p(:)%v(i), dt * momeq0(:,i)/2d0);
    enddo
!    h_new%f = h%f + dt * maseq0/2d0;
    h_new%f = sumreal(h%f, dt * maseq0/2d0);
    !-------------------------------------------------

    !Second RK step-----------------------------------
    if (vectinv) then
      do i=1,3
        u_new%p(:)%v(i) = vecn%p(:)%v(i);
      enddo
    else
      do i=1,3
        uh%p(:)%v(i) = vecn%p(:)%v(i);
      enddo
    endif

    call calcop(u_new,h_new,htr,uh,bath,uhq_perp,zeta,div,gradb, &
        ke,f,fv,momflux,tau,lap,bih,bihm,mesh);

    do i=1,3
!      momeq1(:,i) = -(uhq_perp%p(:)%v(i) &
!          +gradb%p(:)%v(i)+momflux%p(:)%v(i))+bih%p(:)%v(i);
      momeq1(:,i) = -sumreal(uhq_perp%p(:)%v(i),gradb%p(:)%v(i))
      momeq1(:,i) = sumreal(momeq1(:,i),-momflux%p(:)%v(i))
      momeq1(:,i) = sumreal(momeq1(:,i),bih%p(:)%v(i))
    enddo
    maseq1 = -div%f;

    do i=1,3
!      vecn%p(:)%v(i) = veco%p(:)%v(i) + dt * momeq1(:,i)/2;
      vecn%p(:)%v(i) = sumreal(veco%p(:)%v(i), dt * momeq1(:,i)/2);
    enddo
!    h_new%f = h%f + dt * maseq1/ 2d0;
    h_new%f = sumreal(h%f, dt * maseq1/ 2d0);
    !------------------------------------------------

    !Third  RK step-----------------------------------
    if (vectinv) then
      do i=1,3
        u_new%p(:)%v(i) = vecn%p(:)%v(i);
      enddo
    else
      do i=1,3
        uh%p(:)%v(i) = vecn%p(:)%v(i);
      enddo
    endif
    call calcop(u_new,h_new,htr,uh,bath,&
        uhq_perp,zeta,div,gradb,ke,f,fv,momflux,tau,lap,bih,bihm,mesh);

    do i=1,3
!      momeq2(:,i) = -(uhq_perp%p(:)%v(i) &
!          +gradb%p(:)%v(i)+momflux%p(:)%v(i))+bih%p(:)%v(i);
      momeq2(:,i) = -sumreal(uhq_perp%p(:)%v(i),gradb%p(:)%v(i))
      momeq2(:,i) = sumreal(momeq2(:,i),-momflux%p(:)%v(i))
      momeq2(:,i) = sumreal(momeq2(:,i),bih%p(:)%v(i))
    enddo
    maseq2 = -div%f;

    do i=1,3
!      vecn%p(:)%v(i) = veco%p(:)%v(i) + dt * momeq2(:,i);
      vecn%p(:)%v(i) = sumreal(veco%p(:)%v(i), dt * momeq2(:,i));
    enddo
!    h_new%f = h%f + dt * maseq2;
    h_new%f = sumreal(h%f, dt * maseq2);
    !-------------------------------------------------


    !-------------------------------------------------
    if (vectinv) then
      do i=1,3
        u_new%p(:)%v(i) = vecn%p(:)%v(i);
      enddo
    else
      do i=1,3
        uh%p(:)%v(i) = vecn%p(:)%v(i);
      enddo
    endif
    call calcop(u_new,h_new,htr,uh,bath,uhq_perp,zeta,div,gradb,ke, &
        f,fv,momflux,tau,lap,bih,bihm,mesh);

    do i=1,3
!      momeq3(:,i) = -(uhq_perp%p(:)%v(i) &
!          +gradb%p(:)%v(i)+momflux%p(:)%v(i))+bih%p(:)%v(i);
      momeq3(:,i) = -sumreal(uhq_perp%p(:)%v(i),gradb%p(:)%v(i))
      momeq3(:,i) = sumreal(momeq3(:,i),-momflux%p(:)%v(i))
      momeq3(:,i) = sumreal(momeq3(:,i),bih%p(:)%v(i))
    enddo
    maseq3 = -div%f;!+bihm%f;
    !-------------------------------------------------


    !-------------------------------------------------
    momeq = (momeq0 +2d0*momeq1 +2d0*momeq2 + momeq3)/6d0;
    maseq = (maseq0 +2d0*maseq1 +2d0*maseq2 + maseq3)/6d0;

!    write(*,*) maxval(maseq0), 2d0*maxval(maseq1), &
!        2d0*maxval(maseq2), maxval(maseq3);
!    stop;

    do i=1,3
!      vecn%p(:)%v(i) = veco%p(:)%v(i) + dt * momeq(:,i)+Forcmo%p(:)%v(i);
      vecn%p(:)%v(i) = dt * sumreal(momeq(:,i),Forcmo%p(:)%v(i));
      vecn%p(:)%v(i) = sumreal(vecn%p(:)%v(i), veco%p(:)%v(i));
    enddo
!    h_new%f = h%f + dt*(maseq+Forcma%f);
    h_new%f = dt*sumreal(maseq,Forcma%f);
    h_new%f = sumreal(h_new%f,h%f);
    !-------------------------------------------------
    if (vectinv) then
      do i=1,3
        u_new%p(:)%v(i) = vecn%p(:)%v(i);
      enddo
    else
      do i=1,3
        uh%p(:)%v(i) = vecn%p(:)%v(i);
      enddo
    endif
    call calcop(u_new,h_new,htr,uh,bath,uhq_perp,zeta,div, &
        gradb,ke,f,fv,momflux,tau,lap,bih,bihm,mesh);
    if (saveenergy) then
        call savetimesteps(ke,zeta,time,j,mesh)
    endif

    return
  end subroutine ode_rk4

  subroutine savetimesteps(ke,zeta,time,j,mesh)
      type(grid_structure),intent(in) :: mesh
      type(scalar_field),intent(in) :: ke,zeta
      real*8,intent(in) :: time
      ! real*8 :: kesum(mesh%nt), zetasum(mesh%nv)
      integer :: j

      write(j,*) time/86400d0,sum(ke%f*mesh%tr(:)%areag),sum(zeta%f*mesh%hx(:)%areag)
      call FLUSH(j)

  endsubroutine
end module time_derivB
