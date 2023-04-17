program shall_watA
  use swma_operators

  !Use main grid data structures
  use datastruct, only: &
    grid_structure, &
    scalar_field, &
    vectorinterpol_methods, &
    vector_field_cart

  use smeshpack
  use time_derivA


  implicit none
  type(grid_structure)     :: mesh
  type(scalar_field)       :: f,bath,Bop,Bo,Bp,Bn,Btr
  type(vector_variable_uv) :: uop,uo,up,un,utr,ued
  type(vector_variable_uv) :: gradb,Forcmo
  type(vector_variable_uv) :: uhq_perp,lap,bih

  type(scalar_field)       :: zeta,div,Forcma
  real*8,allocatable       :: ui(:,:),bi(:),momeq(:,:),maseq(:)
  integer                  :: i,j

  real*8,allocatable       :: red(:,:),rtr(:)
  real*8                   :: u0,alphan,alphap,H0(10),nu(10),error


  real*8             :: dt,t,tmax


  character(len = 100) :: fileg,filen,filee,filen2
  character(len = 2) :: glevel

!----------------------------------------------------------------------------------------------------------------
  glevel = 'g5';
  write(fileg,'(4A)') './grid/',glevel


  write(*,'(2A)') 'Reading Grid: ',fileg
  call loadgrid(mesh, fileg);

  allocate(Bop%f(mesh%nv),Bo%f(mesh%nv),Btr%f(mesh%nt));
  allocate(uop%u(mesh%nv),uo%u(mesh%nv),up%u(mesh%nv),un%u(mesh%nv), &
           utr%u(mesh%nt),ued%u(mesh%ne));
  allocate(uop%v(mesh%nv),uo%v(mesh%nv),up%v(mesh%nv),un%v(mesh%nv), &
           utr%v(mesh%nt),ued%v(mesh%ne));

  allocate(bi(mesh%nv),ui(mesh%nv,2));
  allocate(maseq(mesh%nv),momeq(mesh%nv,2));

  allocate(f%f(mesh%nv));

  allocate(bath%f(mesh%nv));
  allocate(div%f(mesh%nv),zeta%f(mesh%nv));

  allocate(gradb%u(mesh%nv),gradb%v(mesh%nv));

  allocate(Bo%fexact(mesh%nv),Btr%fexact(mesh%nt));
  allocate(Bp%f(mesh%nv));
  allocate(uo%uexact(mesh%nv),utr%uexact(mesh%nt));
  allocate(uo%vexact(mesh%nv),utr%vexact(mesh%nt));

  allocate(div%fexact(mesh%nv),zeta%fexact(mesh%nv));
  allocate(gradb%uexact(mesh%nv),gradb%vexact(mesh%nv));

  allocate(Forcmo%u(mesh%nv),Forcmo%v(mesh%nv));
  allocate(Forcma%f(mesh%nv));

  allocate(red(mesh%nv,2),rtr(mesh%nv));

  Forcma%f = 0d0;
  Forcmo%u = 0d0; Forcmo%v = 0d0;

  dt   = 2d2;

!  nu = (/0d13,5d13,10d13,15d13,20d13,50d13,100d13,200d13,500d13,1000d13/);
!  nu = -4.9d-6
! nu = -(/0d15,1d-13,1d-6,2d-6,2.5d-6,3d-6,3.5d-6,4d-6,4.5d-6,5d-6/);
  nu = 0d0
  H0 = (/1d-3,1d-2,1d-1,1d0,3d0,5d0,10d0,20d0,40d0,60d0/);

  do i =2,10
    if (i<10) then
      write(filen,'(A,A,A,I1,A)') './resultA/TC10/hollsinst/',glevel,'/loghollsvisc',i,'.dat';
    else
      write(filen,'(A,A,A,I2,A)') './resultA/TC10/hollsinst/',glevel,'/loghollsvisc',i,'.dat';
    endif
    open(unit = 101, file = trim(filen));
    
    uo%u = 0d0;
    uo%v = 0d0;
    f%f = 2d0*omega*dsin(mesh%v(:)%lat);

    u0 = 2d0*pi*rad/(86400d0*12);
    do j=1,3
      uo%u = uo%u + u0*dcos(mesh%v(:)%lat)*mesh%v(:)%c%lonvec(j)*mesh%v(:)%c%nr(j);
      uo%v = uo%v + u0*dcos(mesh%v(:)%lat)*mesh%v(:)%c%lonvec(j)*mesh%v(:)%c%tg(j);
    enddo
    bath%f = -(1d0/g*(rad*omega*u0+u0**2*.5d0*1d0))*dsin(mesh%v(:)%lat)**2

    ui(:,1) = uo%u;
    ui(:,2) = uo%v;

    Bo%f = H0(i);
    Bi = Bo%f

    dt   = 2d0**2*200d0;

    call ode_rk4(uo,Bo,up,Bp,bath,gradb,momeq,maseq,f,Forcmo,Forcma,mesh,nu(i),dt)
    Forcmo%u = -(up%u-ui(:,1))/dt;
    Forcmo%v = -(up%v-ui(:,2))/dt;
    Forcma%f = -(Bp%f-Bi)/dt;

    red(:,1) = up%u-ui(:,1);
    red(:,2) = up%v-ui(:,2);

    rtr = Bp%f-bi;
    alphap = 1d-5/norm(normm(red));
    up%u = alphap*red(:,1)+ui(:,1);
    up%v = alphap*red(:,2)+ui(:,2);
    Bp%f = alphap*rtr+bi;

    j=1;
    error = 1;
    do while (error>0)
      call ode_rk4(uo,Bo,up,Bp,bath,gradb,momeq,maseq,f,Forcmo,Forcma,mesh,nu(i),dt)
      red(:,1) = up%u-ui(:,1);
      red(:,2) = up%v-ui(:,2);
      rtr = Bp%f-bi;
      alphan = 1d-5/norm(normm(red));
      uo%u = alphan*red(:,1)+ui(:,1);
      uo%v = alphan*red(:,2)+ui(:,2);
      Bo%f = alphan*rtr+bi;


      error =  abs(alphap-alphan);

      write(*,'(I6,A,ES14.7,A,ES14.7,A,ES14.7,A,ES14.7,A,ES14.7)') j,' |',H0(i),' |',nu(i),' |',error,' |',alphan,' |', &
                                                                   dt/log(1d0/alphan)/86400d0;


      write(101,*) j,nu(i),H0(i),error,alphan;
      call flush(101);
      alphap = alphan;
      j=j+1
      if (j>1d6) then
        exit;
      endif
    enddo
    close(101)
  enddo

  stop;


end program shall_watA
