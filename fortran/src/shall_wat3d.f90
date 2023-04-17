program adv_eq
  use swmc_operators

  !Use main grid data structures
  use datastruct, only: &
    grid_structure, &
    scalar_field, &
    vectorinterpol_methods, &
    vector_field_cart

  use smeshpack


  implicit none
  type(grid_structure) :: mesh
  type(scalar_field)   :: uo,Bo,un,Bn,bath
  type(scalar_field)   :: div,gradb,gradk,ke,zeta,q,uhq_perp,Forcmo,Forcma,f,fv
  type(scalar_field)   :: uh,B_ed,phi
  integer             :: i
  integer,parameter   :: nlin=1

  real*8             :: dt,t,tmax

  real*8              :: Ti,Tn,Ek,Ep,qi,qn
  real*8,allocatable  :: A(:),Ae(:),Ax(:),ui(:),bi(:),momeq(:),maseq(:)

  character(len = 100) :: fileg,filen
  character(len = 2) :: glevel

!----------------------------------------------------------------------------------------------------------------
  glevel = 'g5';
  write(fileg,'(4A)') './grid/',glevel


  write(*,'(2A)') 'Reading Grid: ',fileg
  call loadgrid(mesh, fileg)
  allocate(A(mesh%nt),Ae(mesh%ne),Ax(mesh%nv));

  A  = mesh%tr(:)%areag;
  Ae = mesh%ed(:)%leng*mesh%edhx(:)%leng;
  Ax = mesh%hx(:)%areag;


  allocate(uo%f(mesh%ne),Bo%f(mesh%nt));

  allocate(bath%f(mesh%nt),div%f(mesh%nt), ke%f(mesh%nt), phi%f(mesh%nt), &
           B_ed%f(mesh%ne),uh%f(mesh%ne),gradb%f(mesh%ne),gradk%f(mesh%ne),uhq_perp%f(mesh%ne),f%f(mesh%ne), &
           zeta%f(mesh%nv),q%f(mesh%nv),fv%f(mesh%nv));

  allocate(momeq(mesh%ne),maseq(mesh%nt));
  allocate(Forcmo%f(mesh%ne),Forcma%f(mesh%nt));

  allocate(uo%fexact(mesh%ne),Bo%fexact(mesh%nt));
  allocate(div%fexact(mesh%nt), &
           ke%fexact(mesh%nt),gradb%fexact(mesh%ne),uhq_perp%fexact(mesh%ne), &
           zeta%fexact(mesh%nv),q%fexact(mesh%nv));

  allocate(ui(mesh%ne),bi(mesh%nt));

  call inicond(8,uo,Bo,bath,div,gradb,ke,zeta,q,uhq_perp,f,fv,mesh,filen,nlin);

!  bath%f = 0d0;
!  CALL RANDOM_NUMBER(uo%f)
!  CALL RANDOM_NUMBER(Bo%f)
!  Bo%f = Bo%f+1d3

  Forcmo%f = 0d0;
  Forcma%f = 0d0;

  phi%f = g*(bath%f+Bo%f);
  t = 0d0;
  dt = 200*2d0;
  tmax = 15d0*86400d0;

  ui = uo%f;
  bi = Bo%f;

  call scalar_tr2ed(Bo, B_ed, mesh);
  Ep = sum(A*g*Bo%f*(Bo%f*.5d0+bath%f));
  Ek = sum(Ae*B_ed%f*uo%f**2*.5d0);
  Ti = Ek+Ep;


  call errorop(uo,Bo,uh,B_ed,phi,zeta,q,div,uhq_perp,gradb,ke,f,fv,mesh,nlin,dt,glevel)
  qi = sum(Ax*q%f)


  i = 1;
  do while (t<=tmax)
    call ode_rk4(uo,bo,un,Bn,bath,momeq,maseq,f,fv,Forcmo,Forcma,mesh,dt,nlin);
    phi%f = g*(Bn%f+bath%f);
    call calcop(un,Bn,uh,B_ed,phi,zeta,q,div,uhq_perp,gradb,gradk,ke,f,fv,mesh,nlin,dt);

    qn = sum(Ax*q%f)
    call scalar_tr2ed(Bn, B_ed, mesh);
    Ep = sum(A*g*Bn%f*(Bn%f*.5d0+bath%f));
    Ek = sum(Ae*B_ed%f*un%f**2*.5d0);
    Tn = Ek+Ep;

    write(*,'(ES14.7,A,ES13.7,A,ES13.7,A,ES13.7,A,ES13.7)') t/86400d0,' | ',maxval(abs(un%f-uo%f)/maxval(ui)),  &
                                                                      ' | ', maxval(abs(Bn%f-Bo%f)/maxval(bi)), &
                                                                      ' | ',abs(Tn-Ti)/Ti,' | ',abs(qn-qi);
!    write(100,*) t/86400d0,maxval(abs(un%f-uo%f)/maxval(ui)),maxval(abs(Bn%f-Bo%f)/maxval(bi)),(Tn-Ti)/Ti,(qn-qi);
!    call FLUSH()

    uo%f = un%f;
    Bo%f = Bn%f;
    t=t+dt

  enddo

  close(100);

end program
