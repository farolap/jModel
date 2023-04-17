program shall_watB
  use swmb_operators

  !Use main grid data structures
  use datastruct, only: &
    grid_structure, &
    scalar_field, &
    vectorinterpol_methods, &
    vector_field_cart

  use smeshpack
  use time_derivB


  implicit none
  type(grid_structure)     :: mesh
  type(scalar_field)       :: bath,Bo,Bp,Btr
  type(vector_variable_uv) :: uo,up
  type(vector_variable_uv) :: gradb,Forcmo
  type(scalar_field)       :: div,zeta,ke,Forcma,f,fv
  real*8                   :: dt,t
  real*8,allocatable       :: momeq(:,:),maseq(:),ui(:,:),bi(:)
  integer                  :: i,j,k

  real*8,allocatable       :: red(:,:),rtr(:)
  real*8                   :: u0,alphan,alphap,nu(10),H0(10),error

  character(len = 100) :: fileg,filen,filee,filen2
  character(len = 2) :: glevel

!----------------------------------------------------------------------------------------------------------------
  glevel = 'g5';
  write(fileg,'(4A)') './grid/',glevel


  write(*,'(2A)') 'Reading Grid: ',fileg
  call loadgrid(mesh, fileg);

  !------------------------------------------------------------
  allocate(f%f(mesh%nt),fv%f(mesh%nv));
  allocate(bath%f(mesh%nv));
  allocate(Bo%f(mesh%nv),Bp%f(mesh%nv));
  allocate(Btr%f(mesh%nt));
  
  allocate(Bo%fexact(mesh%nv));

  allocate(uo%u(mesh%nt),uo%v(mesh%nt), &
           up%u(mesh%nt),up%v(mesh%nt));

  allocate(div%f(mesh%nv),zeta%f(mesh%nv));
  allocate(ke%f(mesh%nv));

  allocate(gradb%u(mesh%nt),gradb%v(mesh%nt));

  allocate(Forcma%f(mesh%nv));

  allocate(maseq(mesh%nv),bi(mesh%nv));
  !------------------------------------------------------------
  allocate(uo%uexact(mesh%nt),uo%vexact(mesh%nt));

  allocate(div%fexact(mesh%nv),zeta%fexact(mesh%nv));
  allocate(ke%fexact(mesh%nv));

  allocate(Forcmo%u(mesh%nt),Forcmo%v(mesh%nt));

  allocate(gradb%uexact(mesh%nt),gradb%vexact(mesh%nt));
  allocate(momeq(mesh%nt,2),ui(mesh%nt,2));

  allocate(red(mesh%nt,2),rtr(mesh%nv));

  !------------------------------------------------------------

  Forcma%f = 0d0;
  Forcmo%u = 0d0; Forcmo%v = 0d0;

  dt   = 2d2;

  nu = -(/0d0,1d-7,1d-6,2d-6,2.5d-6,3d-6,3.5d-6,4d-6,4.5d-6,5d-6/);
!  nu = 0d0
  H0 = 10d0;!(/1d-3,1d-2,1d-1,1d0,3d0,5d0,10d0,20d0,40d0,60d0/);

  do i =3,10
    if (i<10) then
      write(filen,'(A,A,A,I1,A)') './resultB/TC10/hollsinst/',glevel,'/loghollsvisc',i,'.dat';
    else
      write(filen,'(A,A,A,I2,A)') './resultB/TC10/hollsinst/',glevel,'/loghollsvisc',i,'.dat';
    endif
    open(unit = 101, file = trim(filen));
    
    uo%u = 0d0;
    uo%v = 0d0;
    f%f = 2d0*omega*dsin(mesh%tr(:)%b%lat);
    fv%f = 2d0*omega*dsin(mesh%v(:)%lat);

    u0 = 2d0*pi*rad/(86400d0*12);
    do j=1,3
      uo%u = uo%u + u0*dcos(mesh%tr(:)%b%lat)*mesh%tr(:)%b%lonvec(j)*mesh%tr(:)%b%nr(j);
      uo%v = uo%v + u0*dcos(mesh%tr(:)%b%lat)*mesh%tr(:)%b%lonvec(j)*mesh%tr(:)%b%tg(j);
    enddo

    ui(:,1) = uo%u;
    ui(:,2) = uo%v;

    bath%f = -(1d0/g*(rad*omega*u0+u0**2*.5d0*1d0))*dsin(mesh%v(:)%lat)**2

    Bo%f = H0(i);
    Bi = Bo%f

    dt   = 2d0**2*200d0;

    call ode_rk4(uo,Bo,up,Bp,bath,momeq,maseq,f,fv,Forcmo,Forcma,mesh,dt,nu(i))
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
      call ode_rk4(uo,Bo,up,Bp,bath,momeq,maseq,f,fv,Forcmo,Forcma,mesh,dt,nu(i))
!      call ode_rk4(up,bp,un,Bn,bath,momeq,maseq,f,fv,nu(i),Forcmo,Forcma,mesh,dt,nlin);
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


      write(101,*) j,H0(i),nu(i),error,alphan;
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

end program shall_watB
