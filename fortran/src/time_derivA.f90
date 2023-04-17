module time_derivA
  use constants
  use smeshpack
  use datastruct, only: &
    grid_structure, &
    vector_field_cart, &
    scalar_field
  use swma_operators, only: &
    calc_op

  contains

  subroutine ode_rk4(u,h,u_new,h_new,bath,gradb,momeq,maseq,f,Forcmo,Forcma,mesh,dt,time,j)
    implicit none

    type(grid_structure),      intent(in)         :: mesh
    type(scalar_field),        intent(in)         :: f,bath,h,Forcma
    type(scalar_field),        intent(inout)      :: h_new
    type(vector_field_cart),   intent(in)         :: u,Forcmo
    type(vector_field_cart),   intent(inout)      :: u_new,gradb
    real*8,                    intent(inout)      :: momeq(mesh%nv,3),maseq(mesh%nv)
    real*8,                    intent(in)         :: dt,time
    integer,                   intent(in)         :: j
    type(vector_field_cart)                       :: uhq_perp,lap,bih,utr
    type(scalar_field)                            :: ke,div,zeta
    real*8                                        :: momeq0(mesh%nv,3),maseq0(mesh%nv), &
      momeq1(mesh%nv,3),maseq1(mesh%nv), &
      momeq2(mesh%nv,3),maseq2(mesh%nv), &
      momeq3(mesh%nv,3),maseq3(mesh%nv)
    integer :: i

    allocate(utr%p(mesh%nt));
    allocate(zeta%f(mesh%nv),div%f(mesh%nv),ke%f(mesh%nv));

    allocate(lap%p(mesh%nv));
    allocate(bih%p(mesh%nv));
    allocate(uhq_perp%p(mesh%nv));

    do i=1,3
      u_new%p(:)%v(i)=u%p(:)%v(i);
    enddo
    h_new%f=h%f;

    !First RK step------------------------------------
    call calc_op(u,utr,h,bath,ke,zeta,div,gradb,uhq_perp,lap,bih,f,mesh)
    do i =1,3
!      momeq0(:,i) = -(uhq_perp%p(:)%v(i) + gradb%p(:)%v(i))+bih%p(:)%v(i);
!      u_new%p(:)%v(i) = u%p(:)%v(i)  + dt*momeq0(:,i)/2d0;
      momeq0(:,i) = -sumreal(uhq_perp%p(:)%v(i), gradb%p(:)%v(i));
      momeq0(:,i) = sumreal(momeq0(:,i), bih%p(:)%v(i));
      u_new%p(:)%v(i) = sumreal(u%p(:)%v(i),dt*momeq0(:,i)/2d0);
    enddo
    maseq0 = -div%f;
!    h_new%f = h%f + dt * maseq0/ 2d0;
    h_new%f = sumreal(h%f, dt * maseq0/ 2d0);
    !-------------------------------------------------

    momeq = momeq0
    maseq = maseq0

    !Second RK step-----------------------------------
    call calc_op(u_new,utr,h_new,bath,ke,zeta,div,gradb,uhq_perp,lap,bih,f,mesh)
    do i =1,3
!      momeq1(:,i) = -(uhq_perp%p(:)%v(i) + gradb%p(:)%v(i))+bih%p(:)%v(i);
!      u_new%p(:)%v(i) = u%p(:)%v(i)  + dt*momeq1(:,i)/2d0;
      momeq1(:,i) = -sumreal(uhq_perp%p(:)%v(i), gradb%p(:)%v(i));
      momeq1(:,i) = sumreal(momeq1(:,i), bih%p(:)%v(i));
      u_new%p(:)%v(i) = sumreal(u%p(:)%v(i),dt*momeq1(:,i)/2d0);
    enddo
    maseq1 = -div%f;
!    h_new%f = h%f + dt * maseq1/ 2d0;
    h_new%f = sumreal(h%f, dt * maseq1/ 2d0);

    momeq = momeq1;
    maseq = maseq1;
    
    !-------------------------------------------------

    !Third  RK step-----------------------------------
    call calc_op(u_new,utr,h_new,bath,ke,zeta,div,gradb,uhq_perp,lap,bih,f,mesh)
    do i =1,3
!      momeq2(:,i) = -(uhq_perp%p(:)%v(i) + gradb%p(:)%v(i))+bih%p(:)%v(i);
!      u_new%p(:)%v(i) = u%p(:)%v(i)  + dt*momeq2(:,i)/2d0;
      momeq2(:,i) = -sumreal(uhq_perp%p(:)%v(i), gradb%p(:)%v(i));
      momeq2(:,i) = sumreal(momeq2(:,i), bih%p(:)%v(i));
      u_new%p(:)%v(i) = sumreal(u%p(:)%v(i) , dt*momeq2(:,i));
    enddo
    maseq2 = -div%f;
!    h_new%f = h%f + dt * maseq2;
    h_new%f = sumreal(h%f, dt * maseq2);
    !-------------------------------------------------


    call calc_op(u_new,utr,h_new,bath,ke,zeta,div,gradb,uhq_perp,lap,bih,f,mesh)
    do i =1,3
!      momeq3(:,i) = -(uhq_perp%p(:)%v(i) + gradb%p(:)%v(i))+bih%p(:)%v(i);
      momeq3(:,i) = -sumreal(uhq_perp%p(:)%v(i), gradb%p(:)%v(i));
      momeq3(:,i) = sumreal(momeq3(:,i), bih%p(:)%v(i));
    enddo
    maseq3 = -div%f;
    !-------------------------------------------------
!    momeq = (momeq0 +2d0*momeq1 +2d0*momeq2 + momeq3)/6d0;
!    maseq = (maseq0 +2d0*maseq1 +2d0*maseq2 + maseq3)/6d0;
    do i=1,3
        momeq(:,i) = sumreal(momeq0(:,i),2d0*momeq1(:,i))
        momeq(:,i) = sumreal(momeq(:,i),2d0*momeq2(:,i))
        momeq(:,i) = sumreal(momeq(:,i),momeq3(:,i))
        momeq(:,i) = momeq(:,i)/6d0
    enddo

    maseq = sumreal(maseq0,2d0*maseq1)
    maseq = sumreal(maseq,2d0*maseq2)
    maseq = sumreal(maseq,maseq3)
    maseq = maseq/6d0

    do i =1,3
!      u_new%p(:)%v(i) = u%p(:)%v(i)  + dt*(momeq(:,i)+Forcmo%p(:)%v(i));
      u_new%p(:)%v(i) = sumreal(u%p(:)%v(i), dt*momeq(:,i));
      u_new%p(:)%v(i) = sumreal(u_new%p(:)%v(i), dt*Forcmo%p(:)%v(i));
    enddo
!    h_new%f = h%f + dt*(maseq+Forcma%f);
    h_new%f = sumreal(h%f, dt*maseq);
    h_new%f = sumreal(h_new%f, dt*Forcma%f);

    call calc_op(u_new,utr,h_new,bath,ke,zeta,div,gradb,uhq_perp,lap,bih,f,mesh);
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

        write(j,*) time/86400d0,sum(ke%f*mesh%hx(:)%cdareag),sum(zeta%f*mesh%hx(:)%cdareag)
        call FLUSH(j)

    endsubroutine

end module time_derivA
