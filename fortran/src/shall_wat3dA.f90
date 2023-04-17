program shall_watA
    use swma_operators

    !Use main grid data structures
    use datastruct, only: &
        grid_structure, &
        scalar_field, &
        vectorinterpol_methods, &
        vector_field_cart

    use lmeshpack
    use time_derivA
    !   use init_condA


    implicit none
    type(grid_structure)    :: mesh
    type(scalar_field)      :: f,bath,Bo,Bp,Btr,Bed
    type(vector_field_cart) :: uo,up,utr,ued
    type(vector_field_cart) :: gradb,Forcmo
    type(scalar_field)      :: zeta,div,Forcma
    real*8,allocatable      :: ui(:,:),bi(:), &
        momeq(:,:),maseq(:),uerr(:,:)
    integer                 :: i,j,k,svf


    real*8             :: dt,t,tmax


    character(len = 100) :: fileg,filen,filee,files
    character(len = 2) :: glevel

    !--------------------------------------------------
    glevel = 'g6';
    write(fileg,'(4A)') '../../grid/gridSCVT/',glevel


    write(*,'(2A)') 'Reading Grid: ',fileg
    call loadpts(mesh, fileg);
    call calcTopoVar(mesh);

    call allocation();

    do i=1,3
        Forcmo%p(:)%v(i) = 0d0;
    enddo
    Forcma%f = 0d0;

    Hm = 1d0
    call inicond(10,uo,utr,Bo,Btr,Bed,bath,div,gradb,zeta,f,mesh,filen)
    ! call inicond(10,uo,Bo,Btr,bath,div,gradb,zeta, &
    !     f,mesh,filen)
    write(filee,*) trim(filen(2:100)),'error/',glevel,'/';
    write(files,*) trim(filen(2:100)),'snapshots/',glevel,'/';

    if (.false.) then
        call error_op(uo,utr,Bo,Btr,Bed,zeta,div,gradb,filee,mesh)
        stop
    endif

    t=0*86400d0;

    do i=1,3
        ui(:,i) = uo%p(:)%v(i);
    enddo
    bi = Bo%f;

    dt = 50d0;
    tmax = 1.7d0*86400d0
    ! tmax = 30d0*3600d0

    k=0;
    svf = int(floor(86400d0/dt));

    j=0;
    call saveoperators(uo,Bo,0,files,mesh);


    ! t = t + 5*86400d0;
    ! j=j+5*1;
    ! call loadini(uo,Bo,j,files,mesh);

  
    write(filen,*) trim(filee(2:100)),'timeerror.dat'
    open(unit = 200, file = filen(2:100));
    if (saveenergy) then
        write(filen,'(3A)') trim(filee),'energy.dat';
        open(unit = 201, file = filen(2:100));
    endif

    do i =1,3
        up%p(:)%v(i) =uo%p(:)%v(i);
    enddo
    Bp%f =Bo%f;

    j=j+1;
    k=1;
    do while (t<=tmax)
        call ode_rk4(uo,Bo,up,Bp,bath,gradb,momeq,maseq,f,Forcmo,Forcma,mesh,dt,t,201);

        do i=1,3
        uo%p(:)%v(i) = up%p(:)%v(i);
        enddo
        Bo%f = Bp%f;

        do i=1,3
        uerr(:,i) = up%p(:)%v(i)-ui(:,i);
        enddo
        write(*,'(ES14.7,A,ES13.7,A,ES13.7,A,ES13.7)') t/86400d0, &
            ' | ',maxval(normm(uerr)),                          &
            ' | ',maxval(abs(Bp%f-bi)),                         &
            ' | ',maxval(dsqrt(momeq(:,1)**2+momeq(:,2)**2+momeq(:,3)**2))

        write(200,*) t/86400d0,maxval(normm(uerr)), &
                            maxval(abs(Bp%f-bi));
        call FLUSH()

        if (mod(k,svf) == 0) then
            call saveoperators(uo,Bo,j,files,mesh);

            call FLUSH()
            j=j+1;
        endif
        k=k+1;


        t = t+dt;
    enddo
    close(200);
    call saveoperators(uo,Bo,j,files,mesh);
    write(*,*) 'FINISHED!!!!';

    contains
        subroutine allocation()
        allocate(Bo%f(mesh%nv),Btr%f(mesh%nt),Bed%f(mesh%ne));
        allocate(uo%p(mesh%nv),up%p(mesh%nv), &
                utr%p(mesh%nt),ued%p(mesh%ne));
        
        allocate(bi(mesh%nv),ui(mesh%nv,3));
        allocate(maseq(mesh%nv),momeq(mesh%nv,3));
        
        allocate(f%f(mesh%nv));
        
        allocate(bath%f(mesh%nv));
        allocate(div%f(mesh%nv),zeta%f(mesh%nv));
        
        allocate(gradb%p(mesh%nv));
        
        allocate(Bo%fexact(mesh%nv),Btr%fexact(mesh%nt),Bed%fexact(mesh%ne));
        allocate(Bp%f(mesh%nv));
        allocate(uo%pexact(mesh%nv),utr%pexact(mesh%nt));
        
        allocate(div%fexact(mesh%nv),zeta%fexact(mesh%nv));
        allocate(gradb%pexact(mesh%nv));
        
        allocate(Forcmo%p(mesh%nv));
        allocate(Forcma%f(mesh%nv));
        
        allocate(uerr(mesh%nv,3));
        end subroutine allocation

end program shall_watA
