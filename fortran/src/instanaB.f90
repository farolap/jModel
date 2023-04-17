program instanaB
    use swmb_operators

    !Use main grid data structures
    use datastruct, only: &
        grid_structure, &
        scalar_field, &
        vectorinterpol_methods, &
        vector_field_cart

    use smeshpack
    use lmeshpack
    use time_derivB
    ! use init_condB


    implicit none
    type(grid_structure)     :: mesh
    type(scalar_field)       :: bath,Bo,Bp,Btr
    type(vector_field_cart) :: uo,up
    type(vector_field_cart) :: gradb,uh,momfluxc,momfluxv,Forcmo
    type(scalar_field)       :: div,zeta,ke,Forcma,f,fv
    real*8                   :: dt,nu,t
    real*8,allocatable       :: momeq(:,:),maseq(:),ui(:,:),bi(:),erup(:)
    integer                  :: i,j

    real*8 :: alphap,alphan,error,ermin
    real*8,allocatable       :: red(:,:),rtr(:),dist(:)
    real*8 :: equivHeight(10)

    character(len = 100) :: fileg,filen
    character(len = 2) :: glevel

    vectinv=.true.

! -----------------------------------------------------------------------
    glevel = 'g6';
    write(fileg,'(4A)') '../../grid/gridSCVT/',glevel


    write(*,'(2A)') 'Reading Grid: ',fileg
    call loadpts(mesh, fileg);
    call calcTopoVar(mesh);

    !------------------------------------------------------------
    call allocation()
    !------------------------------------------------------------

    allocate(red(mesh%nt,3),rtr(mesh%nv),dist(mesh%nv))

    write(filen,'(3A)') '../resultB/HollingsworthInst/logholls',glevel,'.dat';
    open(unit = 101, file = filen);

    Forcma%f = 0d0;
    Forcmo%p(:)%v(1) = 0d0;
    Forcmo%p(:)%v(2) = 0d0;
    Forcmo%p(:)%v(3) = 0d0;

    dt   = 200d0
    t=0d0
    Bo%f =0d0

    equivHeight = (/1d0,.01d0,.1d0,.5d0,1d0,2d0,3d0,4d0,5d0,10d0/);
    nu=0d0
    ! nu=-5d13
    do i=1,1
        Hm = equivHeight(i)


        
        call inicond(4,uo,Bo,Btr,bath,div,gradb,zeta,ke,f,fv,momfluxc,momfluxv,mesh,filen);

        ! do j =1,mesh%nt
        !     ppi = mesh%tr(j)%c%p;
        !     lat = datan2(ppi(3),dsqrt(ppi(1)**2+ppi(2)**2));
        !     lon = datan2(ppi(2),ppi(1))
        !     dist(j) = dsqrt(lat**2+lon**2);
        ! enddo
        ! k = minloc(dist,1);
        ! Bo%f(k) = Bo%f(k)+Bo%f(k)*1d-3;

        ui(:,1) = uo%p(:)%v(1)
        ui(:,2) = uo%p(:)%v(2)
        ui(:,3) = uo%p(:)%v(3)
        bi = Bo%f

        call ode_rk4(uo,Bo,up,Bp,uh,bath,gradb,div,zeta, &
            momeq,maseq,f,fv,Forcmo,Forcma,mesh,dt,nu,t,201)
        Forcmo%p(:)%v(1) = -(up%p(:)%v(1)-ui(:,1))/dt
        Forcmo%p(:)%v(2) = -(up%p(:)%v(2)-ui(:,2))/dt
        Forcmo%p(:)%v(3) = -(up%p(:)%v(3)-ui(:,3))/dt
        Forcma%f = -(Bp%f-Bi)/dt

        red(:,1) = up%p(:)%v(1)-ui(:,1)
        red(:,2) = up%p(:)%v(2)-ui(:,2)
        red(:,3) = up%p(:)%v(3)-ui(:,3)


        rtr = Bp%f-bi
        alphap = 1d-5/norm(normm(red))
        up%p(:)%v(1) = alphap*red(:,1)+ui(:,1)
        up%p(:)%v(2) = alphap*red(:,2)+ui(:,2)
        up%p(:)%v(3) = alphap*red(:,3)+ui(:,3)

        Bp%f = alphap*rtr+bi

        j=1
        error = 1
        ermin = 1

        do while(error>0)
            call ode_rk4(up,Bp,uo,Bo,uh,bath,gradb,div,zeta, &
                momeq,maseq,f,fv,Forcmo,Forcma,mesh,dt,nu,t,201)
            red(:,1) = uo%p(:)%v(1)-ui(:,1)
            red(:,2) = uo%p(:)%v(2)-ui(:,2)
            red(:,3) = uo%p(:)%v(3)-ui(:,3)
            rtr = Bo%f-bi
            alphan = 1d-5/norm(normm(red))
    
            up%p(:)%v(1) = alphan*red(:,1)+ui(:,1)
            up%p(:)%v(2) = alphan*red(:,2)+ui(:,2)
            up%p(:)%v(3) = alphan*red(:,3)+ui(:,3)

            Bp%f = alphan*rtr+bi
            error =  abs(alphap-alphan)

            write(*,'(I6,A,ES14.7,A,ES14.7,A,ES14.7,A,ES14.7)') j,' |',Hm,' |',&
                error,' |',alphan,' |', dt/log(1d0/alphan)/86400d0;
            ! stop
            write(101,*) j,Hm,error,alphan;
            call flush(101);
            alphap = alphan;
            j=j+1
            if (j>2d4) then
                exit;
            endif


        enddo
        close(101)
        exit
    enddo
    write(filen,'(3A)') '../resultB/HollingsworthInst/eigvalu',glevel,'.dat';
    open(unit = 100, file = filen);
    write(filen,'(3A)') '../resultB/HollingsworthInst/eigvalB',glevel,'.dat';
    open(unit = 101, file = filen);
    
    do i=1,mesh%nv
        write(101,*) Bp%f(i)
    enddo
    do i=1,mesh%nt
        write(100,*) up%p(i)%v
    enddo
    close(100)
    close(101)

    write(*,*) 'FIN!'

    contains
        subroutine allocation()
            allocate(f%f(mesh%nt),fv%f(mesh%nv));
            allocate(bath%f(mesh%nv));
            allocate(Bo%f(mesh%nv),Bp%f(mesh%nv));
            allocate(Btr%f(mesh%nt));
            
            allocate(Bo%fexact(mesh%nv));
            
            allocate(uo%p(mesh%nt), &
                up%p(mesh%nt), &
                uh%p(mesh%nt));
            
            allocate(div%f(mesh%nv),zeta%f(mesh%nv));
            allocate(ke%f(mesh%nv));
            
            allocate(gradb%p(mesh%nt),       &
                momfluxc%p(mesh%nt), &
                momfluxv%p(mesh%nv));
            
            allocate(Forcma%f(mesh%nv));
            
            allocate(maseq(mesh%nv),bi(mesh%nv));

            allocate(uo%pexact(mesh%nt));
            
            allocate(div%fexact(mesh%nv),zeta%fexact(mesh%nv));
            allocate(ke%fexact(mesh%nv));
            
            allocate(Forcmo%p(mesh%nt));
            
            allocate(gradb%pexact(mesh%nt),       &
                momfluxc%pexact(mesh%nt), &
                momfluxv%pexact(mesh%nv));
            allocate(momeq(mesh%nt,3),ui(mesh%nt,3));
            allocate(erup(mesh%nt))
        
        end subroutine allocation
end program instanaB

