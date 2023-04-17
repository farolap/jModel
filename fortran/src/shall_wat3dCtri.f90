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
    real*8                   :: dt,t,tmax,nu
    real*8,allocatable       :: momeq(:),maseq(:),ui(:),bi(:),erup(:)
    integer                  :: j,k,svf

    character(len = 100) :: fileg,filen,filee,files
    character(len = 2) :: glevel

! -----------------------------------------------------------------------
    glevel = 'g7'
    write(fileg,'(4A)') '../../grid/gridSCVT/',glevel


    write(*,'(2A)') 'Reading Grid: ',fileg

    if (.true.) then
        call loadpts(mesh, fileg);
        call calcTopoVar(mesh)
    else
        call loadall(fileg,mesh)
        call calc_sph_vecs(mesh)
    endif

    call allocation_icos();


    Forcma%f = 0d0;
    Forcmo%f = 0d0;

  
    Hm=1d0
    ! call inicond(8,uo,Bo,bath,div,gradb,ke,zeta,q,uperp,f,fv,mesh,filen)
    call inicond(8,uo,Bo,bath,utr,div,gradb,ke,zeta,q,uperp, &
        f,fv,mesh,filen)

    ui = uo%f;
    bi = Bo%f;

    write(filee,*) trim(filen(2:100)),'error/',glevel,'/';
    write(files,*) trim(filen(2:100)),'snapshots/',glevel;

    if (.false.) then
        call errorop(uo,utr,Bo,zeta,q,div,uperp,gradb,ke,fv,mesh,dt,filee)
        stop
    endif
    call saveoperatorstr(uo,Bo,0,files,mesh);


    dt = 50
    ! dt = 1200d0/2
    ! tmax = 1d3*86400d0
    ! tmax = 5d1*86400d0
    tmax = 15d0*86400d0
    t = 0d0

    ! tmax = 10*dt

    k=0;
    svf = int(floor(86400d0/dt))

    j=0;

    write(filen,'(3A)') trim(filee),'timeerror.dat';
    open(unit = 200, file = filen(2:100));

    k = 1;
    j=j+1;

    nu = 0;
    ! nu = -1d12;
    ! nu = -3.8670268252938655e-05
    do while (t<=tmax)
        ! dt = .1*maxval(mesh%ed(:)%leng)*rad/maxval(uo%f);
        call ode_rk4tr(uo,Bo,up,Bp,bath,momeq,maseq,fv,nu,Forcmo,Forcma,mesh,dt)

        write(*,'(ES14.7,A,ES13.7,A,ES13.7,A,ES13.7,A,ES13.7,A,I7)') t/86400d0, &
            ' | ',maxval(abs(up%f-ui)), &
            ' | ',maxval(abs(Bp%f-bi)), &
            ' | ',maxval(abs(momeq)), &
            ' | ',maxval(abs(maseq)), &
            ' | ',maxloc(abs(Bp%f-bi));

            write(200,*) t/86400d0,maxval(abs(up%f-ui)),maxval(abs(Bp%f-bi));
        call FLUSH()

        uo%f = up%f;
        Bo%f = Bp%f;

        ! call errorop(uo,utr,Bo,zeta,q,div,uperp,gradb,ke,fv,mesh,dt,filee)

        t=t+dt;

        if (mod(k,svf) == 0) then
            call saveoperatorstr(uo,Bo,j,files,mesh);

            ! write(200,*) t
            call FLUSH()
            j=j+1;
        endif

        k = k+1;
    enddo
    close(200)

    call saveoperatorstr(uo,Bo,j,files,mesh);


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

