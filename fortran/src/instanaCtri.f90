program shall_watC
    use swmc_operators
    use init_condC, only: inicond


    !Use main grid data structures
    use datastruct, only: &
        grid_structure, &
        scalar_field, &
        vectorinterpol_methods

    use smeshpack
    use lmeshpack
    use time_derivC


    implicit none
    type(grid_structure)     :: mesh
    type(scalar_field)       :: bath,Bo,Bp
    type(scalar_field)       :: uo,up,uperp
    type(scalar_field)       :: gradb,Forcmo
    type(scalar_field)       :: div,zeta,q,ke,Forcma,f,fv
    type(vector_field_cart)  :: utr
    real*8                   :: dt,nu,error,ermin,alphap,alphan, &
        equivHeight(10)
    real*8,allocatable       :: momeq(:),maseq(:),ui(:),bi(:),erup(:),red(:),rtr(:), &
        dist(:)
    integer                  :: i,j,k
    real*8                   :: ppi(3),lat,lon

    character(len = 100) :: fileg,filen
    character(len = 2) :: glevel

! -----------------------------------------------------------------------
    glevel = 'g6';
    write(fileg,'(4A)') '../../grid/gridSCVT/',glevel


    write(*,'(2A)') 'Reading Grid: ',fileg

    if (.true.) then
        call loadpts(mesh, fileg);
        call calcTopoVar(mesh)
    else
        call loadall(fileg,mesh)
        call calc_sph_vecs(mesh)
    endif

    call allocation_icos()

    allocate(red(mesh%ne),rtr(mesh%nt),dist(mesh%nt))

    write(filen,'(3A)') '../resultCtri/HollingsworthInst/logholls',glevel,'.dat';
    open(unit = 101, file = filen);

    Forcma%f = 0d0;
    Forcmo%f = 0d0;

    dt   = 400d0

    nu = 6d14
    equivHeight = (/1d0,0d0,0d0,.0d0,0d0,0d0,0d0,0d0,0d0,1d0/);
    do i=1,10
        Hm = equivHeight(i)
        
        call inicond(4,uo,Bo,bath,utr,div,gradb,ke,zeta,q,uperp, &
            f,fv,mesh,filen)

        ! do j =1,mesh%nt
        !     ppi = mesh%tr(j)%c%p;
        !     lat = datan2(ppi(3),dsqrt(ppi(1)**2+ppi(2)**2));
        !     lon = datan2(ppi(2),ppi(1));
        !     dist(j) = dsqrt(lat**2+lon**2);
        ! enddo
        ! k = minloc(dist,1);
        ! Bo%f(k) = Bo%f(k)+Bo%f(k)*1d-3;

        ui = uo%f
        bi = Bo%f

        call ode_rk4tr(uo,Bo,up,Bp,bath,momeq,maseq,fv,nu,Forcmo,Forcma,mesh,dt)
        Forcmo%f = -(up%f-ui)/dt
        Forcma%f = -(Bp%f-Bi)/dt

        red = up%f-ui
        rtr = Bp%f-bi
        alphap = 1d-5/(norm(red))
        up%f = alphap*red+ui
        Bp%f = alphap*rtr+bi

        j=1
        error = 1
        ermin = 1

        do while(error>0)
            call ode_rk4tr(up,Bp,uo,Bo,bath,momeq,maseq,fv,nu,Forcmo,Forcma,mesh,dt)
            red = uo%f-ui
            rtr = Bo%f-bi
            alphan = 1d-5/(norm(red))
    
            up%f = alphan*red+ui
            Bp%f = alphan*rtr+bi
            error =  abs(alphap-alphan)

            write(*,'(I6,A,ES14.7,A,ES14.7,A,ES14.7,A,ES14.7)') j,' |',Hm,' |',&
                error,' |',alphan,' |', dt/log(1d0/alphan)/86400d0;
            ! stop
            write(101,*) j,Hm,error,alphan;
            call flush(101);
            alphap = alphan;
            j=j+1
            if (j>1d5) then
                exit;
            endif


        enddo
        close(101)
        exit
    enddo
    write(filen,'(3A)') '../resultCtri/HollingsworthInst/eigvalu',glevel,'.dat';
    open(unit = 100, file = filen);
    write(filen,'(3A)') '../resultCtri/HollingsworthInst/eigvalB',glevel,'.dat';
    open(unit = 101, file = filen);
    
    do i=1,mesh%nt
        write(101,*) Bp%f(i)
    enddo
    do i=1,mesh%ne
        write(100,*) up%f(i)*mesh%ed(i)%nr
    enddo
    close(100)
    close(101)


    write(*,*) 'FIN!'

    contains
        subroutine allocation_icos()
            allocate(f%f(mesh%ne),fv%f(mesh%nv));
            allocate(bath%f(mesh%nt));
            allocate(Bo%f(mesh%nt),Bp%f(mesh%nt));
            allocate(Bo%fexact(mesh%nt));
            
            allocate(uo%f(mesh%ne),up%f(mesh%ne));
            
            allocate(uperp%f(mesh%ne),uperp%fexact(mesh%ne));
            
            
            allocate(div%f(mesh%nt),zeta%f(mesh%nv),q%f(mesh%nv),q%fexact(mesh%nv));
            allocate(ke%f(mesh%nt));
            
            allocate(gradb%f(mesh%ne));
            
            allocate(Forcma%f(mesh%nt));
            
            allocate(maseq(mesh%nt),bi(mesh%nt));
            !------------------------------------------------------------
            allocate(uo%fexact(mesh%ne));
            
            allocate(div%fexact(mesh%nt),zeta%fexact(mesh%nv));
            allocate(ke%fexact(mesh%nt));
            
            allocate(Forcmo%f(mesh%ne));
            
            allocate(gradb%fexact(mesh%ne));
            allocate(momeq(mesh%ne),ui(mesh%ne));
            allocate(erup(mesh%ne))

            allocate(utr%p(mesh%nt),utr%pexact(mesh%nt))

        end subroutine allocation_icos
end program shall_watC

