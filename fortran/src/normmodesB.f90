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


    implicit none
    type(grid_structure)     :: mesh
    type(scalar_field)       :: bath,Bo,Bp,Btr
    type(vector_field_cart) :: uo,up
    type(vector_field_cart) :: gradb,uh,Forcmo
    type(scalar_field)       :: div,zeta,ke,Forcma,f,fv
    real*8                   :: dt,t,tmax,nu
    real*8,allocatable       :: momeq(:,:),maseq(:),u(:),v(:)
    integer                  :: i,j,k,svf

    character(len = 100) :: fileg,filen
    character(len = 2) :: glevel

    nlin = .false.
    vectinv = .true.
    !------------------------------------------------------------------------------------
    glevel = 'g2';
    write(fileg,'(4A)') '../../grid/gridSCVT/',glevel


    write(*,'(2A)') 'Reading Grid: ',fileg
    call loadpts(mesh, fileg);
    call calcTopoVar(mesh);

    !------------------------------------------------------------
    call allocation()
    !------------------------------------------------------------
    allocate(u(mesh%nt),v(mesh%nt))

    Forcma%f = 0d0;
    do i=1,3
        Forcmo%p(:)%v(i) = 0d0;
    enddo
    bath%f =0d0

    f%f = 1.4584d-4;!*dsin(mesh%tr(:)%b%lat);
    dt = 1d1;

    write(filen,*) '../resultB/modes/norm',glevel,'.dat'
    open(unit = 100, file = filen(2:100));

    !------------------------------------------------------------------------------------

    nu =-1d-4
    do i =1,mesh%nv+mesh%nt*2

        write(*,*) i,mesh%nv+mesh%nt*2

        Bo%f = 0d0;
        do j=1,mesh%nt
            uo%p(j)%v = 0d0
            uh%p(j)%v = 0d0
        enddo

        if (i<=mesh%nv) then
            Bo%f(i) = 1d0
        elseif (i<=mesh%nv+mesh%nt) then
            uo%p(i-mesh%nv)%v = mesh%tr(i-mesh%nv)%b%lonvec 
            uh%p(i-mesh%nv)%v = mesh%tr(i-mesh%nv)%b%lonvec*1d4
        elseif (i<=mesh%nv+2*mesh%nt) then
            uo%p(i-mesh%nv-mesh%nt)%v = mesh%tr(i-mesh%nv-mesh%nt)%b%latvec
            uh%p(i-mesh%nv-mesh%nt)%v = mesh%tr(i-mesh%nv-mesh%nt)%b%latvec*1d4
        endif
        call ode_rk4(uo,Bo,up,Bp,uh,bath,gradb,div,zeta, &
            momeq,maseq,f,fv,Forcmo,Forcma,mesh,dt,nu,0d0,200)
        do j=1,mesh%nt
            call convert_vec_cart2sph(mesh%tr(j)%b%p, up%p(j)%v , u(j), v(j))
        enddo
        write(100,*) Bp%f,u,v
    enddo
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

            allocate(gradb%p(mesh%nt));

            allocate(Forcma%f(mesh%nv));

            allocate(maseq(mesh%nv));

            allocate(uo%pexact(mesh%nt));

            allocate(div%fexact(mesh%nv),zeta%fexact(mesh%nv));
            allocate(ke%fexact(mesh%nv));

            allocate(Forcmo%p(mesh%nt));

            allocate(gradb%pexact(mesh%nt));
            allocate(momeq(mesh%nt,3));
        end subroutine allocation

end program shall_watB
