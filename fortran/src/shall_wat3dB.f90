program shall_watB
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
    real*8                   :: dt,t,tmax,nu
    real*8,allocatable       :: momeq(:,:),maseq(:),ui(:,:),bi(:),erup(:)
    integer                  :: i,j,k,svf

    character(len = 100) :: fileg,filen,filee,files
    character(len = 2) :: glevel

    vectinv=.true.

    !------------------------------------------------------------
    glevel = 'g8'
    write(fileg,'(4A)') '../../grid/gridSCVT/',glevel


    write(*,'(2A)') 'Reading Grid: ',fileg
    call loadpts(mesh, fileg);
    call calcTopoVar(mesh);

    !------------------------------------------------------------
    call allocation()
    !------------------------------------------------------------

    Forcma%f = 0d0;
    do i=1,3
        Forcmo%p(:)%v(i) = 0d0;
    enddo

    Hm=1d0
    call inicond(8,uo,Bo,Btr,bath,div,gradb,zeta,ke,f,fv,momfluxc,momfluxv,mesh,filen);

    write(filee,*) trim(filen(2:100)),'error/',glevel,'/';
    write(files,*) trim(filen(2:100)),'snapshots/',glevel;

    dt = 50d0;
    ! tmax = 200d0*86400d0
    !tmax = 1500d0*86400d0
    ! tmax = 15d0*86400d0
    tmax = 30d0*3600d0
    t    = 0d0

    k=0;
    svf = int(floor(86400d0/dt));

    j=0;
    if (.true.) then
        call error_op(uo,Bo,zeta,div,ke,gradb,f,momfluxc,momfluxv,filee,mesh);
        stop
    endif
    call saveoperators(uo,Bo,0,files,mesh);
    ! j=4;
    ! t = t + j*86400d0;
    ! call loadini(uo,Bo,j,files,mesh);

    do i=1,3
        ui(:,i) = uo%p(:)%v(i)
        uh%p(:)%v(i) = uo%p(:)%v(i)*Btr%f
    enddo
    bi = Bo%f;



    write(filen,'(3A)') trim(filee),'timeerror.dat';
    open(unit = 200, file = filen(2:100));

    if (saveenergy) then
        write(filen,'(3A)') trim(filee),'energy.dat';
        open(unit = 201, file = filen(2:100));
    endif

    k = 1;
    j=j+1;

    ! nu = -5d-5;
    ! nu = -3.8670268252938655e-05
    ! nu = -1d-6
    ! nu = -3d-4; ! TC8g7
    ! nu = -4d-4; ! TC8g7
    ! nu = -1d-4; ! TC8g6
    ! nu = -1d-4; ! TC8g6
    ! nu = -1d-4; ! TC8g5
    ! nu = -1d15; ! TC8g4
    ! nu = -1d14; ! TC8g4

    ! nu = -3d-4/4; ! TC2g7

    ! nu = -1d-4; ! TC9g7

    ! nu = -1d-4; ! TC9g7
    ! nu = -0d-4; ! TC10g6

    nu = -0d00; ! TC8g5
    do while (t<=tmax)
        call ode_rk4(uo,Bo,up,Bp,uh,bath,gradb,div,zeta, &
            momeq,maseq,f,fv,Forcmo,Forcma,mesh,dt,nu,t,201)
        t=t+dt;

        erup = 0d0;
        do i=1,3
            ! erup = erup + (up%p(:)%v(i)-ui(:,i))**2;
            erup = sumreal(erup,mulreal(up%p(:)%v(i)-ui(:,i),up%p(:)%v(i)-ui(:,i)));
        enddo
        write(*,'(ES14.7,A,ES13.7,A,ES13.7,A,ES13.7,A,I7)') t/86400d0, &
            ' | ',maxval(dsqrt(erup)), &
            ' | ',maxval(abs(Bp%f-bi)), &
            ' | ',maxval(normm(momeq)), &
            ' | ',maxloc(abs(Bp%f-bi));

        write(200,*) t/86400d0,maxval(dsqrt(erup)),maxval(abs(Bp%f-bi));
        call FLUSH(200)

        do i=1,3
            uo%p(:)%v(i) = up%p(:)%v(i);
        enddo
        Bo%f = Bp%f;
        
        if (mod(k,svf) == 0) then
            call saveoperators(uo,Bo,j,files,mesh);
            j=j+1;
        endif
        k = k+1;
    enddo
    close(200)
    close(201)

    call saveoperators(uo,Bo,j,files,mesh);

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

end program shall_watB
