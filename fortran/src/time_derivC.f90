module time_derivC
    use basic_funcs
    use constants
    use datastruct, only: &
        grid_structure, &
        vector_field_cart, &
        scalar_field
    use swmc_operators
    contains

        subroutine ode_rk4tr(u,h,u_new,h_new,bath,momeq,maseq,fv,nu,Forcmo,Forcma,mesh,dt)
            implicit none

            type(grid_structure), intent(in) :: mesh
            type(scalar_field),   intent(in) :: fv,bath,Forcmo,Forcma,u,h
            type(scalar_field),   intent(inout) :: u_new,h_new
            real*8,               intent(inout) :: momeq(mesh%ne),maseq(mesh%nt)
            real*8,               intent(in) :: dt,nu
            type(scalar_field) :: ke,uh,div,q,zeta,gradb,uhq_perp, &
                h_ed,phi,lap,bih,f,udiff,gradk

            real*8 :: momeq0(mesh%ne),maseq0(mesh%nt), &
                momeq1(mesh%ne),maseq1(mesh%nt), &
                momeq2(mesh%ne),maseq2(mesh%nt), &
                momeq3(mesh%ne),maseq3(mesh%nt)
            logical,parameter :: MPDATA=.false.
            type(scalar_field) :: momeqsf
            type(vector_field_cart) :: vec_tr

            !----------------------------------------------------------------
            allocate(ke%f(mesh%nt),phi%f(mesh%nt),div%f(mesh%nt), &
                gradb%f(mesh%ne),uh%f(mesh%ne), &
                uhq_perp%f(mesh%ne),h_ed%f(mesh%ne),q%f(mesh%nv), &
                zeta%f(mesh%nv),lap%f(mesh%ne),bih%f(mesh%ne),f%f(mesh%ne), &
                gradk%f(mesh%ne), &
                vec_tr%p(mesh%nt),momeqsf%f(mesh%ne));
            !----------------------------------------------------------------


            u_new%f=u%f
            h_new%f=h%f


            !-------------------------------------------------
            phi%f = g*(bath%f+h%f);
            call calcop(u,h,bath,uh,h_ed,zeta,q,div,uhq_perp, &
                gradb,gradk,ke,nu,lap,bih,f,fv,mesh,dt);
            momeq0 = -(uhq_perp%f-gradb%f) - bih%f;
            maseq0 = -div%f;
            !-------------------------------------------------

            !First RK step------------------------------------
            u_new%f = u%f + dt * momeq0/ 2d0;
            h_new%f = h%f + dt * maseq0/ 2d0;

            phi%f = g*(bath%f+h_new%f);
            call calcop(u_new,h_new,bath,uh,h_ed,zeta,q,div,uhq_perp, &
                gradb,gradk,ke,nu,lap,bih,f,fv,mesh,dt);
            momeq1 = -(uhq_perp%f-gradb%f) - bih%f;
            maseq1 = -div%f;
            !-------------------------------------------------

            !Second RK step-----------------------------------
            u_new%f = u%f + dt * momeq1/ 2d0;
            h_new%f = h%f + dt * maseq1/ 2d0;

            phi%f = g*(bath%f+h_new%f);
            call calcop(u_new,h_new,bath,uh,h_ed,zeta,q,div,uhq_perp, &
                gradb,gradk,ke,nu,lap,bih,f,fv,mesh,dt);
            momeq2 = -(uhq_perp%f-gradb%f) - bih%f;
            maseq2 = -div%f;
            !-------------------------------------------------

            !Third  RK step-----------------------------------
            u_new%f = u%f + dt * momeq2;
            h_new%f = h%f + dt * maseq2;

            phi%f = g*(bath%f+h_new%f);
            call calcop(u_new,h_new,bath,uh,h_ed,zeta,q,div,uhq_perp, &
                gradb,gradk,ke,nu,lap,bih,f,fv,mesh,dt);
            momeq3 = -(uhq_perp%f-gradb%f) - bih%f;
            maseq3 = -div%f;
            !-------------------------------------------------

            momeq = (momeq0 +2d0*momeq1 +2d0*momeq2 + momeq3)/6d0;
            maseq = (maseq0 +2d0*maseq1 +2d0*maseq2 + maseq3)/6d0;

            u_new%f = u%f + dt*(momeq+Forcmo%f);
            h_new%f = h%f + dt*(maseq+Forcma%f);

            phi%f = g*(bath%f+h_new%f);
            call calcop(u_new,h_new,bath,uh,h_ed,zeta,q,div,uhq_perp, &
                gradb,gradk,ke,nu,lap,bih,f,fv,mesh,dt);

            return
    end subroutine ode_rk4tr
    !subroutine semiImplicitAB2(uo,up,un,Bo,Bp,Bn)
    !    implicit none

    !    type(grid_structure), intent(in) :: mesh
    !    type(scalar_field),   intent(in) :: fv,bath,Forcmo,Forcma,u,h
    !    type(scalar_field),   intent(inout) :: u_new,h_new
    !    real*8,               intent(inout) :: momeq(mesh%ne),maseq(mesh%nt)
    !    real*8,               intent(in) :: dt,nu
    !    type(scalar_field) :: ke,uh,div,q,zeta,gradb,uhq_perp, &
    !        h_ed,phi,lap,bih,f,udiff
    !    real*8,parameter :: beta=.55d0,gamm=.55d0
    !    real*8,parameter :: betal = 1-beta,gammal=1-gamm

    !    call calcop(uo,Bo,bath,uh,h_ed,zeta,q,div,uhq_perp, &
    !        gradb,gradk,ke,nu,lap,bih,f,fv,mesh,dt);
    !    Go = uhq_perp%f+gradk%f

    !    call calcop(up,Bp,bath,uh,h_ed,zeta,q,div,uhq_perp, &
    !        gradb,gradk,ke,nu,lap,bih,f,fv,mesh,dt);
    !    Gp = uhq_perp%f+gradk%f

    !    Gl = (1.5d0+.1d0)*Gp-(.5d0-.1d0)*Go

    !    uinter = up%f - dt*(betal*gammal*gradb%f+Gl)
    !    umean = gamm*uinter+gammal*up%f
    !    call calcop(umean,Bp,bath,uh,h_ed,zeta,q,div,uhq_perp, &
    !        gradb,gradk,ke,nu,lap,bih,f,fv,mesh,dt);
    !    hn%f = hp%f - dt*div;
    !    un%f = 


    !endsubroutine

end module time_derivC
