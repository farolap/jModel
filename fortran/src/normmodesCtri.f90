program adv_eq
    use swmc_operators

    !Use main grid data structures
    use datastruct, only: &
        grid_structure, &
        scalar_field, &
        vectorinterpol_methods

    use smeshpack
    use lmeshpack
    use time_derivC


    implicit none
    type(grid_structure) :: mesh
    type(scalar_field)   :: f,fv,bath,Bo,Bp,uo,up,Forcma,Forcmo

    integer             :: i,j

    real*8,allocatable  :: momeq(:),maseq(:)

    real*8               :: t,dt,H0
    character(len = 100) :: filen,fileg
    character(len = 2) :: glevel

    nlin=.false.
!---------------------------------------------------------------------------------------

    glevel = 'g2';
    write(fileg,'(4A)') '../../grid/gridSCVT/',glevel


    write(*,'(2A)') 'Reading Grid: ',fileg
    call loadpts(mesh, fileg);
    call calcTopoVar(mesh);

    call allocation_icos();
 
    dt=10;

    f%f = 1.4584d-4;
    fv%f = 1.4584d-4;

    bath%f =0d0;
    Forcmo%f = 0d0
    Forcma%f = 0d0

    write(filen,*) '../resultCtri/modes/norm',glevel,'.dat'
    open(unit = 100, file = filen(2:100));

    do i=1,mesh%nt+mesh%ne
        write(*,'(A,2I5,2I5)') 'i:',i,mesh%nt+mesh%ne
        uo%f = 0;
        Bo%f = 0;
        up%f = 0;
        Bp%f = 0;
        if (i<=mesh%nt) then
            Bo%f(i) = 1;
        else
            uo%f(i-mesh%nt) = 1;
        endif

!--------------------------------------------------------------------------------
        call ode_rk4tr(uo,Bo,up,Bp,bath,momeq,maseq,fv,0d0,Forcmo,Forcma,mesh,dt)

!        call ode_rk4(uo,Bo,up,Bp,uh,phi,bath,div,q,eta,uhq_perp,gradb,
!            gradk,B_ed,momeq,maseq,f,fv,nlin,mesh,dt)
!--------------------------------------------------------------------------------

        write(100,*) Bp%f,up%f
    enddo
    close(100)

    contains
        subroutine allocation_icos()
            allocate(f%f(mesh%ne),fv%f(mesh%nv));
            allocate(bath%f(mesh%nt));
            allocate(Bo%f(mesh%nt),Bp%f(mesh%nt));
            allocate(uo%f(mesh%ne),up%f(mesh%ne));
            allocate(Forcma%f(mesh%nt));
            allocate(maseq(mesh%nt));
            allocate(Forcmo%f(mesh%ne));
            allocate(momeq(mesh%ne));

        end subroutine allocation_icos

end program
