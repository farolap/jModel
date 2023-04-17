program construct_grid
    !Use main grid data structures
    use datastruct
    use smeshpack
    use lmeshpack
    use basic_funcs

    implicit none
    type(grid_structure) :: mesh
    ! logical             :: load_pts
    character(len = 100) :: fileg
    character(len = 2)   :: glevel
    ! integer :: i

    glevel = 'g5';
    write(fileg,'(4A)') '../../grid/gridSCVT/',glevel

    ! call loadpts2(mesh,fileg)
    
    call loadpts(mesh, fileg);
    call calcTopoVar(mesh);

    write(*,*) sum(mesh%tr(:)%areag)-4*pi;
    write(*,*) sum(mesh%hx(:)%areag)-4*pi;
    write(*,*) sum(mesh%ed(:)%leng*mesh%edhx(:)%leng)/2-4*pi;

    write(*,*) 'FINISHED!!!'


end program construct_grid