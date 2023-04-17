program shall_watC
  use swmc_operators
  use init_condC, only: inicondhx


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
  type(scalar_field)       :: div,zeta,q,ke,Forcma,f,ft
  real*8                   :: dt,t,tmax,nu
  real*8,allocatable       :: momeq(:),maseq(:),ui(:),bi(:),erup(:)
  integer                  :: i,j,k,svf

  character(len = 100) :: fileg,filen,filee,filen2
  character(len = 2) :: glevel

  ! ---------------------------------------------------------------------------------------------------
  glevel = 'g5';
  write(fileg,'(4A)') '../../grid/gridSCVT/',glevel


  write(*,'(2A)') 'Reading Grid: ',fileg
  call loadpts(mesh, fileg);
  call calcTopoVar(mesh);
  call allocation_icosd();
  
  Forcma%f = 0d0;
  Forcmo%f = 0d0;
  
  call inicondhx(8,uo,Bo,bath,div,gradb,ke,zeta,q,uperp,f,ft,mesh,filen);
  up%fexact = uperp%fexact*f%f;
  ! bo%f =3d3;
  ! call random_number(uo%f); uo%f = 4d1*uo%f;
  ! call random_number(bo%f); bo%f = 3d3*(bo%f-sum(bo%f)/mesh%nv);
  ui = uo%f;
  bi = Bo%f;

  write(filee,*) trim(filen(2:100)),'error/',glevel,'/';

  call errorophx(uo,Bo,zeta,q,div,uperp,gradb,ke,ft,mesh,dt,filee);
  ! stop;

  call saveoperatorsve(filen,glevel,uo,Bo,0,mesh);

  t    = 0d0
  dt   = 2.0d2;
  tmax = 5d0*86400d0;

  write(filen2,'(3A)') trim(filee),'timeerror.dat';
  open(unit = 200, file = filen2(2:100));

  k = 1;
  j=1;
  svf = 173*9999999;

  nu = 0;

  ! nu = -8d12;
  ! nu = -3.8670268252938655e-05
  do while (t<=tmax)
    call ode_rk4hx(uo,Bo,up,Bp,bath,momeq,maseq,ft,nu,Forcmo,Forcma,mesh,dt);

    write(*,'(ES14.7,A,ES13.7,A,ES13.7,A,ES13.7,A,I5)') t/86400d0, &
      ' | ',maxval(abs(up%f-ui)), &
      ' | ',maxval(abs(Bp%f-bi)), &
      ' | ',maxval(abs(momeq)), &
      ' | ',maxloc(abs(Bp%f-bi));

    write(200,*) t/86400d0,maxval(abs(up%f-ui)),maxval(abs(Bp%f-bi));
    call FLUSH()

    uo%f = up%f;
    Bo%f = Bp%f;
    t=t+dt;
    if (mod(k,svf) == 0) then
      call saveoperatorsve(filen,glevel,uo,Bo,j,mesh);

      call FLUSH()
      j=j+1;
    endif


    ! if (t>=200*dt) then
    ! if (t>=15d0*86400) then
      ! exit
    ! endif
    k = k+1;
  enddo
  close(200)

  call saveoperatorsve(filen,glevel,uo,Bo,j,mesh);

  ! call vort_tr(uo,zeta,mesh);
  ! call div_hx(uo,div,mesh);
  
  ! write(filen2,'(3A)') trim(filee),'h.dat';
  ! open(unit = 100, file = filen2(2:100));

  ! write(filen2,'(3A)') trim(filee),'u.dat';
  ! open(unit = 101, file = filen2(2:100));

  ! do i=1,mesh%nv
  !   write(100,*) Bo%f(i),Bi(i)
  !   write(103,*) div%f(i),div%fexact(i)
  ! enddo
  ! do i=1,mesh%nt
  !   write(102,*) zeta%f(i),zeta%fexact(i)
  !   write(103,*) div%f(i),div%fexact(i)
  ! enddo
  ! do i=1,mesh%ne
  !   write(101,*) uo%f(i),ui(i)
  ! enddo

  write(*,*) 'FIN!'
  stop;

  contains
    subroutine allocation_icosd()
      allocate(f%f(mesh%ne),ft%f(mesh%nt));
      allocate(bath%f(mesh%nv));
      allocate(Bo%f(mesh%nv),Bp%f(mesh%nv));
      allocate(Bo%fexact(mesh%nv));
    
      allocate(uo%f(mesh%ne),up%f(mesh%ne));
    
      allocate(uperp%f(mesh%ne),uperp%fexact(mesh%ne));
    
    
      allocate(div%f(mesh%nv),zeta%f(mesh%nt),q%f(mesh%nt),q%fexact(mesh%nt));
      allocate(ke%f(mesh%nv));
    
      allocate(gradb%f(mesh%ne));
    
      allocate(Forcma%f(mesh%nv));
    
      allocate(maseq(mesh%nv),bi(mesh%nv));
      !------------------------------------------------------------
      allocate(uo%fexact(mesh%ne));
    
      allocate(div%fexact(mesh%nv),zeta%fexact(mesh%nt));
      allocate(ke%fexact(mesh%nv));
    
      allocate(Forcmo%f(mesh%ne));
    
      allocate(gradb%fexact(mesh%ne));
      allocate(momeq(mesh%ne),ui(mesh%ne));
      allocate(erup(mesh%ne))
    
    end subroutine allocation_icosd
end program shall_watC

