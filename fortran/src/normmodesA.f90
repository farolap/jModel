program normmodesA
    use swma_operators

    !Use main grid data structures
    use datastruct, only: &
    grid_structure, &
    scalar_field, &
    vectorinterpol_methods, &
    vector_field_cart

    use lmeshpack
    use time_derivA


    implicit none
    type(grid_structure)    :: mesh
    type(scalar_field)      :: f,bath,Bo,Bp,Btr,Bed
    type(vector_field_cart) :: uo,up,utr,ued
    type(vector_field_cart) :: gradb,Forcmo
    type(scalar_field)      :: Forcma
    real*8,allocatable      :: ui(:,:),bi(:),momeq(:,:),maseq(:),u(:),v(:)
    integer                 :: i,j,k,svf


    real*8             :: dt,t


    character(len = 100) :: fileg,filen
    character(len = 2) :: glevel

    !------------------------------------------------------------------------------------
    nlin = .false.

    glevel = 'g2';
    write(fileg,'(4A)') '../../grid/gridSCVT/',glevel


    write(*,'(2A)') 'Reading Grid: ',fileg
    call loadpts(mesh, fileg);
    call calcTopoVar(mesh);

    call allocation();
    allocate(u(mesh%nv),v(mesh%nv))

    do i=1,3
    Forcmo%p(:)%v(i) = 0d0;
    enddo
    Forcma%f = 0d0;

    bath%f = 0d0;
    f%f = 1.4584d-4;
    dt = 1d1;


    write(filen,*) '../resultA/modes/norm',glevel,'.dat'
    open(unit = 100, file = filen(2:100));

    do i =1,mesh%nv+mesh%nv*2
    write(*,*) i,mesh%nv+mesh%nv*2

    do j =1,mesh%nv
        Bo%f(j) = 0d0;
        uo%p(j)%v = 0d0;
    enddo
    if (i<=mesh%nv) then
        Bo%f(i) = 1d0
    elseif (i<=mesh%nv+mesh%nv) then
        uo%p(i-mesh%nv)%v = mesh%v(i-mesh%nv)%c%lonvec 
    elseif (i<=mesh%nv+2*mesh%nv) then
        uo%p(i-mesh%nv-mesh%nv)%v = mesh%v(i-mesh%nv-mesh%nv)%c%latvec
    endif

    call ode_rk4(uo,Bo,up,Bp,bath,gradb,momeq,maseq,f,Forcmo,Forcma,mesh,dt);
    do j=1,mesh%nv
        call convert_vec_cart2sph(mesh%v(j)%p, up%p(j)%v , u(j), v(j))
    enddo
    write(100,*) Bp%f,u,v

    enddo
    contains
        subroutine allocation()
            allocate(Bo%f(mesh%nv),Btr%f(mesh%nt),Bed%f(mesh%ne));
            allocate(uo%p(mesh%nv),up%p(mesh%nv), &
                   utr%p(mesh%nt),ued%p(mesh%ne));

            allocate(bi(mesh%nv),ui(mesh%nv,3));
            allocate(maseq(mesh%nv),momeq(mesh%nv,3));

            allocate(f%f(mesh%nv));

            allocate(bath%f(mesh%nv));

            allocate(gradb%p(mesh%nv));

            allocate(Bo%fexact(mesh%nv),Btr%fexact(mesh%nt),Bed%fexact(mesh%ne));
            allocate(Bp%f(mesh%nv));
            allocate(uo%pexact(mesh%nv),utr%pexact(mesh%nt));

            allocate(gradb%pexact(mesh%nv));

            allocate(Forcmo%p(mesh%nv));
            allocate(Forcma%f(mesh%nv));

        end subroutine allocation
end program normmodesA
