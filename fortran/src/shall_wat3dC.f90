program shall_watC
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
  type(grid_structure)     :: mesh
  type(scalar_field)       :: bath,Bo,Bp
  type(scalar_field)       :: uo,up,uperp
  type(scalar_field)       :: gradb,Forcmo
  type(scalar_field)       :: div,zeta,q,ke,Forcma,f,fv
  real*8                   :: dt,t,tmax,nu
  real*8,allocatable       :: momeq(:),maseq(:),ui(:),bi(:),erup(:)
  integer                  :: i,j,k,svf

  character(len = 100) :: fileg,filen,filee,filen2
  character(len = 2) :: glevel

!----------------------------------------------------------------------------------------------------------------
  glevel = 'g5';
  write(fileg,'(4A)') '../../grid/gridSCVT/',glevel


  write(*,'(2A)') 'Reading Grid: ',fileg
  call loadpts(mesh, fileg);
  call calcTopoVar(mesh);


  call allocation_icos();

  Forcma%f = 0d0;
  Forcmo%f = 0d0;

  call inicond(8,uo,Bo,bath,div,gradb,ke,zeta,q,uperp,f,fv,mesh,filen);
  ui = uo%f;
  bi = Bo%f;

  write(filee,*) trim(filen(2:100)),'error/',glevel,'/';

  ! call errorop(uo,Bo,zeta,q,div,uperp,gradb,ke,f,fv,mesh,dt,filee);
  call errorop(uo,Bo,zeta,q,div,uperp,gradb,ke,fv,mesh,dt,filee);

  call saveoperators(filen,glevel,uo,Bo,0,mesh);

  ! tmax = 5d1*86400d0;
  t    = 0d0
  dt   = 2.5d2;
  tmax = 15d0*86400d0;

  ! write(filen,'(3A)') trim(filen),'snapshots/';

  write(filen2,'(3A)') trim(filee),'timeerror.dat';
  open(unit = 200, file = filen2(2:100));


  k = 1;
  j=1;
  svf = 173*9999999;

  nu = -5d-5;
  ! nu = -3.8670268252938655e-05
  do while (t<=tmax)
    call ode_rk4(uo,Bo,up,Bp,bath,momeq,maseq,f,fv,nu,Forcmo,Forcma,mesh,dt);
    ! call errorop(up,Bp,zeta,q,div,uperp,gradb,ke,f,fv,mesh,dt,filee);

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
      call saveoperators(filen,glevel,uo,Bo,j,mesh);

!      write(200,*) t
      call FLUSH()
      j=j+1;
    endif


    ! if (t>=10*dt) then
    ! if (t>=5d0*86400) then
      ! exit
    ! endif
    k = k+1;
  enddo
  close(200)

  call saveoperators(filen,glevel,uo,Bo,j,mesh);

  call vort_hx(uo,zeta,mesh);
  call div_tr(uo,div,mesh);
  
  write(filen2,'(3A)') trim(filee),'h.dat';
  open(unit = 100, file = filen2(2:100));

  write(filen2,'(3A)') trim(filee),'u.dat';
  open(unit = 101, file = filen2(2:100));

  do i=1,mesh%nt
    write(100,*) Bo%f(i),Bi(i)
    write(103,*) div%f(i),div%fexact(i)
  enddo
  do i=1,mesh%nv
    write(102,*) zeta%f(i),zeta%fexact(i)
    write(103,*) div%f(i),div%fexact(i)
  enddo
  do i=1,mesh%ne
    write(101,*) uo%f(i),ui(i)
  enddo

  write(*,*) 'FIN!'
  stop;

  contains
    subroutine allocation_icos()
      allocate(f%f(mesh%ne),fv%f(mesh%nv));
      allocate(bath%f(mesh%nt));
      allocate(Bo%f(mesh%nt),Bp%f(mesh%nt));
    
      allocate(Bo%fexact(mesh%nt));
    
      allocate(uo%f(mesh%ne), &
              up%f(mesh%ne));
    
      allocate(uperp%f(mesh%ne),uperp%fexact(mesh%ne));
    
    
      allocate(div%f(mesh%nt),zeta%f(mesh%nv),q%f(mesh%nv),q%fexact(mesh%nv));
      allocate(ke%f(mesh%nt));
    
      allocate(gradb%f(mesh%ne));
    
      allocate(Forcma%f(mesh%ne));
    
      allocate(maseq(mesh%nt),bi(mesh%nt));
      !------------------------------------------------------------
      allocate(uo%fexact(mesh%ne));
    
      allocate(div%fexact(mesh%nt),zeta%fexact(mesh%nv));
      allocate(ke%fexact(mesh%nt));
    
      allocate(Forcmo%f(mesh%ne));
    
      allocate(gradb%fexact(mesh%ne));
      allocate(momeq(mesh%ne),ui(mesh%ne));
      allocate(erup(mesh%ne))
    
    end subroutine allocation_icos
end program shall_watC

