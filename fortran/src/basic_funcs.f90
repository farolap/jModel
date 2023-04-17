module basic_funcs
  use datastruct


  contains

    function sumreal(p1,p2) result(p3)
        implicit none
        real*8,intent(in) :: p1(:),p2(:)
        real*8 :: p3(size(p1,1))
        integer :: i
        ! real*8 :: partial_sum(size(p1,1)),total_sum(size(p1,1))

        ! !$OMP PARALLEL DO REDUCTION(+:p1)
        ! do i = 1, size(p1,1)
        !     p3(i) =  p1(i)+p2(i)
        ! end do
        ! !$OMP END PARALLEL DO

        !$omp parallel do &
        !$omp default(none) &
        !$omp shared(p1,p2,p3) &
        !$omp private(i) &
        !$omp schedule(static)
        do i=1,size(p1,1)
            p3(i) = p1(i)+p2(i)
        enddo
        !$omp end parallel do
    endfunction sumreal

    function mulreal(p1,p2) result(p3)
        implicit none
        real*8,intent(in) :: p1(:),p2(:)
        real*8 :: p3(size(p1,1))
        integer :: i

        !$omp parallel do &
        !$omp default(none) &
        !$omp shared(p1,p2,p3) &
        !$omp private(i) &
        !$omp schedule(static)
        do i=1,size(p1,1)
            p3(i) = p1(i)*p2(i)
        enddo
        !$omp end parallel do
    endfunction mulreal

    function trian_centroid(pq,tri)
      real*8,intent(in)  :: pq(:,:)
      integer,intent(in) :: tri(:,:)
      real*8             :: normpc(size(tri,1)),ppc(size(tri,1),3)
      real*8             :: trian_centroid(size(tri,1),3)
      integer :: i

      ppc = pq(tri(:,1),:) + pq(tri(:,2),:) + pq(tri(:,3),:);
      normpc = normm(ppc);

      do i =1,3
        trian_centroid(:,i) = ppc(:,i)/normpc;
      enddo

      return

    endfunction trian_centroid

    function trian_circum(pq,tri)
      real*8,intent(in)  :: pq(:,:)
      integer,intent(in) :: tri(:,:)
      real*8             :: normp(size(tri,1)),pp(size(tri,1),3)
      real*8             :: trian_circum(size(tri,1),3)
      integer :: i
      
      pp = cross2d_product(pq(tri(:,2),:)-pq(tri(:,1),:), &
                            pq(tri(:,3),:)-pq(tri(:,1),:))
      normp = normm(pp);

      do i =1,3
        trian_circum(:,i) = pp(:,i)/normp;
      enddo
    endfunction trian_circum

    function edgemidpoint(pq,sh)
      real*8,intent(in)  :: pq(:,:)
      integer,intent(in) :: sh(:,:)
      real*8             :: normpe(size(sh,1)),pe(size(sh,1),3), &
          pq1(size(sh,1),3),pq2(size(sh,1),3)
      real*8             :: edgemidpoint(size(sh,1),3)
      integer :: i

      pe = 0;
      pq1 = pq(sh(:,1),:);
      pq2 = pq(sh(:,2),:);
      pe = pq1+pq2;

      normpe = normm(pe);
      do i=1,3
        edgemidpoint(:,i) = pe(:,i)/normpe
      enddo

      return
    endfunction edgemidpoint
    
    function hexag_cell_area(pq,pppe)
      real*8,intent(in)  :: pq(:,:),pppe(:,:,:)
      real*8             :: Ai(size(pq,1)),hexag_cell_area

      Ai = sphtriarea2d(pppe(:,1,:),pppe(:,2,:),pq);
      hexag_cell_area = sum(Ai);
      

    endfunction hexag_cell_area

    function triangle2d_area(pe,pq,pp,etri,edv)
      real*8,intent(in)   :: pe(:,:),pq(:,:),pp(:,:)
      integer,intent(in)  :: edv(:,:),etri(:,:)
      real*8              :: triangle2d_area(size(etri,1),3)
      integer             :: i,etrii(size(etri,1))

      triangle2d_area = 0d0;

      do i =1,3
        etrii = etri(:,i);
        triangle2d_area(:,i) = sphtriarea2d(pe(etrii,:), pq(edv(etrii,1),:), pp) + &
                                sphtriarea2d(pe(etrii,:), pq(edv(etrii,2),:), pp);
      enddo

      return

    endfunction triangle2d_area

    function planartri2d_area(pqt)
      real*8 :: pqt(:,:,:),planartri2d_area(size(pqt,1))
      real*8 :: p1(size(pqt,1),3),p2(size(pqt,1),3),p3(size(pqt,1),3)
      real*8 :: a(size(pqt,1)),b(size(pqt,1)),c(size(pqt,1)), s(size(pqt,1))

      p1 = pqt(:,1,:);
      p2 = pqt(:,2,:);
      p3 = pqt(:,3,:);

      a = normm(p2-p1);
      b = normm(p3-p2);
      c = normm(p1-p3);
      s = (a + b + c) / 2
      planartri2d_area = (s*(s-a)*(s-b)*(s-c)) ** .5d0;
    endfunction planartri2d_area

    function triangle2dpq_area(pe,pq,pp,tri,etri)
      real*8,intent(in)   :: pe(:,:),pq(:,:),pp(:,:)
      integer,intent(in)  :: tri(:,:),etri(:,:)
      real*8              :: triangle2dpq_area(size(etri,1),3)
      integer             :: i,trii(size(tri,1)),lloc(3,2)

      lloc(1,:) = (/3,1/);
      lloc(2,:) = (/1,2/);
      lloc(3,:) = (/2,3/);

      triangle2dpq_area = 0d0;

      do i =1,3
        trii = tri(:,i);
        triangle2dpq_area(:,i) = sphtriarea2d(pq(trii,:), pe(etri(:,lloc(i,1)),:), pp) + &
                                 sphtriarea2d(pq(trii,:), pe(etri(:,lloc(i,2)),:), pp);
      enddo

      return

    endfunction triangle2dpq_area

    function edge_tg(pq,sh)
      real*8,intent(in)   :: pq(:,:)
      integer,intent(in)  :: sh(:,:)
      real*8              :: pq1(size(sh,1),3),pq2(size(sh,1),3), &
                            tg(size(sh,1),3),edge_tg(size(sh,1),3),&
                            normtg(size(sh,1))
      integer :: i

      pq1 = pq(sh(:,1),:);
      pq2 = pq(sh(:,2),:);

      tg = pq2-pq1;
      normtg = normm(tg)
      do i=1,3
        edge_tg(:,i) = tg(:,i)/normtg;
      enddo
      return
    endfunction edge_tg

    function edge_nr(pe,tg)
      real*8,intent(in)   :: pe(:,:),tg(:,:)
      real*8              :: nr(size(tg,1),3),edge_nr(size(tg,1),3)
      integer             :: i


      nr = -cross2d_product(pe,tg);
      do i=1,3
        edge_nr(:,i) = nr(:,i)/normm(nr);
      enddo

    endfunction edge_nr

    function spherical_coord(pc)
      real*8,intent(in)   :: pc(:,:)
      real*8              :: lat(size(pc,1)),lon(size(pc,1)),rade(size(pc,1)), &
                              spherical_coord(size(pc,1),2)


      rade = dsqrt(pc(:,1)**2+pc(:,2)**2);
      lat = datan2(pc(:,3),rade);
      lon = datan2(pc(:,2),pc(:,1));

      spherical_coord(:,1) = lon;
      spherical_coord(:,2) = lat;


    endfunction spherical_coord

    function rotx(theta)
      real*8 :: theta(:)
      real*8 :: rotx(size(theta,1),3,3)

      rotx(:,1,1) = 1; rotx(:,1,2) = 0;           rotx(:,1,3) = 0;
      rotx(:,2,1) = 0; rotx(:,2,2) = dcos(theta); rotx(:,2,3) = -dsin(theta);
      rotx(:,3,1) = 0; rotx(:,3,2) = dsin(theta); rotx(:,3,3) = dcos(theta);
    endfunction rotx

    function roty(theta)
      real*8 :: theta(:)
      real*8 :: roty(size(theta,1),3,3)

      roty(:,1,1) = dcos(theta);  roty(:,1,2) = 0; roty(:,1,3) = dsin(theta);
      roty(:,2,1) = 0;            roty(:,2,2) = 1; roty(:,2,3) = 0;
      roty(:,3,1) = -dsin(theta); roty(:,3,2) = 0; roty(:,3,3) = dcos(theta);
    endfunction roty

    function roty1(theta)
      real*8 :: theta
      real*8 :: roty1(3,3)

      roty1(1,1) = dcos(theta);  roty1(1,2) = 0; roty1(1,3) = dsin(theta);
      roty1(2,1) = 0;            roty1(2,2) = 1; roty1(2,3) = 0;
      roty1(3,1) = -dsin(theta); roty1(3,2) = 0; roty1(3,3) = dcos(theta);
    endfunction roty1

    function rotz(theta)
      real*8 :: theta(:)
      real*8 :: rotz(size(theta,1),3,3)

      rotz(:,1,1) = dcos(theta); rotz(:,1,2) = -dsin(theta); rotz(:,1,3) = 0;
      rotz(:,2,1) = dsin(theta); rotz(:,2,2) =  dcos(theta); rotz(:,2,3) = 0;
      rotz(:,3,1) = 0;           rotz(:,3,2) = 0;            rotz(:,3,3) = 1;
    endfunction rotz

    ! -------------------------------------------------------------
    function cross2d_product(a,b)
        real*8, intent(in):: a(:,:)
        real*8, intent(in):: b(:,:)
        real*8 :: cross2d_product(SIZE(a,1),SIZE(a,2))
        integer :: i

        !$omp parallel do &
        !$omp default(none) &
        !$omp shared(a,b,cross2d_product) &
        !$omp private(i) &
        !$omp schedule(static)
        do i =1,size(cross2d_product,1)
            cross2d_product(i,1) = a(i,2)*b(i,3) - a(i,3)*b(i,2);
            cross2d_product(i,2) = a(i,3)*b(i,1) - a(i,1)*b(i,3);
            cross2d_product(i,3) = a(i,1)*b(i,2) - a(i,2)*b(i,1);
        enddo
        !$omp end parallel do
        ! cross2d_product(:,1) = a(:,2)*b(:,3) - a(:,3)*b(:,2);
        ! cross2d_product(:,2) = a(:,3)*b(:,1) - a(:,1)*b(:,3);
        ! cross2d_product(:,3) = a(:,1)*b(:,2) - a(:,2)*b(:,1);

        return
    end function cross2d_product

    function norm(p)
      !-----------------------------------------
      ! NORM
      ! Calculates the euclidian norm of a vector
      !----------------------------------------------
      real*8, intent(in) :: p(:)
      real*8 :: norm

      norm=dot_product( p, p)
      norm=dsqrt(norm)

      return
    end function norm

    function normm(v)
      real*8  :: v(:,:)
      real*8  :: normm(SIZE(v,1))
      integer :: i

      normm = 0d0;
      do i=1,SIZE(v,2)
        normm = normm+v(:,i)**2
      enddo
      normm = dsqrt(normm);

    end function normm

    function arclen2d(v1,v2)
      real*8  :: v1(:,:),v2(:,:)
      real*8  :: d(SIZE(v1,1))
      real*8  :: arclen2d(SIZE(v1,1))

      d = (v1(:,1)-v2(:,1))**2+(v1(:,2)-v2(:,2))**2+(v1(:,3)-v2(:,3))**2
      
      arclen2d = 2d0 * dasin(dsqrt(d)*.5d0)

      
    end function arclen2d

    function arclen(v1,v2)
        real*8  :: v1(3),v2(3)
        real*8  :: d
        real*8  :: arclen

        ! d = (v1(1)-v2(1))**2+(v1(2)-v2(2))**2+(v1(3)-v2(3))**2    
        ! arclen = 2d0 * dasin(dsqrt(d)*.5d0)

        d = dot_product(v1-v2, v1-v2)
        arclen = 2._r8 * dasin(dsqrt(d)*.5d0)

    end function arclen

    function sphtriarea2d(p1, p2, p3)
      real*8,intent(in) :: p1(:,:),p2(:,:),p3(:,:)
      real*8 :: a(size(p1,1)),b(size(p1,1)),c(size(p1,1)),s(size(p1,1))
      real*8 :: e(size(p1,1)),tmp(size(p1,1))
      real*8 :: sphtriarea2d(size(p1,1))

      !Calculate the sides length's and the semiperimiter
      s=0d0
      a = arclen2d(p1,p2);  !arcdistll(lon(1),lat(1), lon(2),lat(2))
      b = arclen2d(p2,p3);  !arcdistll(lon(2),lat(2), lon(3),lat(3))
      c = arclen2d(p3,p1);  !arcdistll(lon(3),lat(3), lon(1),lat(1))
      s = (a+b+c)/2d0;

      !Calculate spherical triangle excess using L'Huilier's Theorem
      tmp = abs(dtan(s/2d0)*dtan((s-a)/2d0)*dtan((s-b)/2d0)*dtan((s-c)/2d0));
      !Round off error might give almost zero negative numbers => assume zero
      e = 4d0*datan(dsqrt(tmp))
      sphtriarea2d=e

      return
    end function sphtriarea2d

    function vec2d_sph2cart(p)
      real*8, intent(in)  :: p(:,:)
      real*8              :: lat(size(p,1)),lon(size(p,1)), &
                              vec2d_sph2cart(size(p,1),2,3)
      
      lat = datan2(p(:,3),dsqrt(p(:,1)**2+p(:,2)**2));
      lon = datan2(p(:,2),p(:,1));

      vec2d_sph2cart(:,1,1) = -dsin(lon)
      vec2d_sph2cart(:,1,2) = dcos(lon)
      vec2d_sph2cart(:,1,3) = 0d0;

      vec2d_sph2cart(:,2,1) = -dsin(lat)*dcos(lon)
      vec2d_sph2cart(:,2,2) = -dsin(lat)*dsin(lon)
      vec2d_sph2cart(:,2,3) = -(-dcos(lat));
    end function vec2d_sph2cart

    function proj_vec_sphere(v, p)
      !-----------------------------------------------------------
      !  Projects a vector 'v' on the plane tangent to a sphere
      !   Uses the the tangent plane relative to the unit sphere's
      !   point 'p', in cartesian coords
      !-----------------------------------------------------------
      real*8, intent(in), dimension(1:3) :: v
      real*8, intent(in), dimension(1:3) :: p
      real*8, dimension(1:3)  :: proj_vec_sphere
  
      proj_vec_sphere(1:3)=&
           v(1:3)-dot_product(v,p)*p(1:3)/norm(p)
  
      return
    end function proj_vec_sphere

    subroutine convert_vec_sph2cart(vlon, vlat, p, v)
        !---------------------------------------------------------------------
        !	CONVERT_VEC_SPH2CART
        !
        !   Recieves a point p=(x,y,z) and a vector at this point in lon, lat system
        !   in radians (vlon, vlat), ie, 
        !      vlon=West-East direction
        !      vlat=South-North direction
        !   Returns the vector at this point in cartesian reference (v)
        !---------------------------------------------------------------------
        !Point cartesian coords
        real*8, intent(in) :: p(1:3)
        real*8, intent(in) :: vlon
        real*8, intent(in) :: vlat
        !Cartesian vector on point
        real*8, intent(out) :: v(1:3)
        !Auxiliar variables
        real*8:: r
        real*8:: rho

        r=dsqrt(p(1)**2+p(2)**2+p(3)**2)
        rho=dsqrt(p(1)**2+p(2)**2)

        !Case where the point is in the north or south pole
        if(rho==0)then
           v(1)=vlat
           v(2)=vlon
           v(3)=0
           return
        else    
           !The minus sign in vlat is due to the diference between spherical coords and
           !   geographical coords
           v(1)=-vlon*(p(2)/rho) - (vlat*p(1)*p(3))/(rho)
           v(2)=vlon*(p(1)/rho) - (vlat*p(2)*p(3))/(rho)
           v(3)=vlat*rho
        end if

        return    
    end subroutine convert_vec_sph2cart

    subroutine convert_vec_cart2sph(p, v , vlon, vlat)
        !Point cartesian coords
        real(r8), intent(in) :: p(1:3)
        !Cartesian vector on point
        real(r8), intent(in) :: v(1:3)
        !Spherical coord vector on point
        real(r8), intent(out) :: vlat
        real(r8), intent(out) :: vlon
        !Auxiliar variables
        real(r8):: r
        real(r8):: rho
        real(r8):: rvec(1:3)
        real(r8):: latvec(1:3)
        real(r8):: lonvec(1:3)
        real(r8):: zero
        real(r8):: test

        zero=0
        r=dsqrt(p(1)**2+p(2)**2+p(3)**2)
        rho=dsqrt(p(1)**2+p(2)**2)

        !Case where the point is in the north or south pole
        if(rho<10*eps)then
           !print*, "Pole:", v
            vlon=v(2)
            vlat=v(1)
            return
        end if

        rvec=(/p(1),p(2),p(3)/)
        rvec=rvec/r

        latvec=(/-p(1)*p(3),-p(2)*p(3),rho**2/)
        latvec=latvec/rho

        lonvec=(/-p(2), p(1), zero /)
        lonvec=lonvec/rho

        test=dot_product(v,rvec)
        if(abs(test)>10e-5)then
            print *,"CONVERT_VEC_CART2SPH Warning: Vector not tangent to sphere."
            print '(a,3f10.6)',"Vector:",v(1:3)
            print '(a,3f10.6)',"Point:",rvec(1:3)
            print '(a,f10.6)',"Dot Product:", test
            stop
        end if
        vlat=dot_product(v,latvec)
        vlon=dot_product(v,lonvec)

        return    
    end subroutine convert_vec_cart2sph

    function cross_product(a,b)
      real*8, intent(in):: a(3)
      real*8, intent(in):: b(3)
      real*8 :: cross_product(3)
      
      cross_product(1) = a(2)*b(3) - a(3)*b(2);
      cross_product(2) = a(3)*b(1) - a(1)*b(3);
      cross_product(3) = a(1)*b(2) - a(2)*b(1);

      return
    end function cross_product

    function sphtriarea(p1, p2, p3)
      real*8,intent(in) :: p1(3),p2(3),p3(3)
      real*8 :: a,b,c,s
      real*8 :: e,tmp
      real*8 :: sphtriarea

      !Calculate the sides length's and the semiperimiter
      s=0d0
      a = arclen(p1,p2);  !arcdistll(lon(1),lat(1), lon(2),lat(2))
      b = arclen(p2,p3);  !arcdistll(lon(2),lat(2), lon(3),lat(3))
      c = arclen(p3,p1);  !arcdistll(lon(3),lat(3), lon(1),lat(1))
      s = (a+b+c)/2d0;

      !Calculate spherical triangle excess using L'Huilier's Theorem
      tmp = abs(dtan(s/2d0)*dtan((s-a)/2d0)*dtan((s-b)/2d0)*dtan((s-c)/2d0));
      !Round off error might give almost zero negative numbers => assume zero
      e = 4d0*datan(dsqrt(tmp))
      sphtriarea=e

      return
    end function sphtriarea

    function l2er(var,A,normal)
      real*8, intent(in):: var(:,:)
      real*8,intent(in) :: A(:)
      logical,intent(in) :: normal
      real*8            :: S,S2,l2er
  
      l2er=0d0;
  
      S  = sum((var(:,1)-var(:,2))**2*A)/sum(A);
      S2 = sum(var(:,2)**2*A)/sum(A);
  
      if (normal) then
        l2er = dsqrt(S)/dsqrt(S2);
      else
        l2er = dsqrt(S);
      endif
      ! l2er = dsqrt(S);
  
      return
    end function l2er

    function lmaxer(var,normal)
      real*8, intent(in):: var(:,:)
      logical,intent(in) :: normal
      real*8            :: S,S2,lmaxer
  
      lmaxer=0d0;
  
      S  = maxval(abs(var(:,1)-var(:,2)));
      S2 = maxval(abs(var(:,2)));
  
      if (normal) then
        lmaxer = S/S2;
      else
        lmaxer = S;
      endif
      ! lmaxer = S;
  
      return
    end function lmaxer

    function roll2d(var,nbr,axis)
      implicit none
      real*8 :: var(:,:)
      real*8 :: roll2d(size(var,1),size(var,2))
      real*8 :: tmp(size(var,1),size(var,2))
      integer :: nbr,axis,nmax,mve

      tmp =0d0;
      nmax = size(var,nbr);
      mve = nbr+1;

      if (axis==1) then
        tmp(1+nbr:nmax,:) = var(1:nmax-nbr,:);
        tmp(1:nbr,:) = var(nmax-nbr+1:nmax,:);
      elseif (axis==2) then
        tmp(:,1+nbr:nmax) = var(:,1:nmax-nbr);
        tmp(:,1:nbr) = var(:,nmax-nbr+1:nmax);
      endif

      roll2d = tmp
    endfunction roll2d

    subroutine hpolyns(deg,P,Q,gP,gQ)
      implicit none
      integer,intent(in) :: deg
      real*8,intent(inout)  :: P(deg,deg+1,deg+1),Q(deg,deg+1,deg+1),gP(deg,deg+1,deg+1,2),gQ(deg,deg+1,deg+1,2)
      real*8  :: xP(deg+1,deg+1),yP(deg+1,deg+1), xQ(deg+1,deg+1),yQ(deg+1,deg+1)
      integer :: i
      
      P = 0;
      Q = 0;
      gP = 0;
      gQ = 0;

      P(1,1,2) = 1;
      Q(1,2,1) = 1;
      gP(1,1,1,1) = 1;
      gQ(1,1,1,2) = 1;

      do i=2,deg
        xP = roll2d(P(i-1,:,:),1,2);
        yP = roll2d(P(i-1,:,:),1,1);
        xQ = roll2d(Q(i-1,:,:),1,2);
        yQ = roll2d(Q(i-1,:,:),1,1);

        P(i,:,:) = xP-yQ;
        Q(i,:,:) = yP+xQ;

        gP(i,:,:,1) = i*P(i-1,:,:);
        gP(i,:,:,2) = -i*Q(i-1,:,:);

        gQ(i,:,:,1) = i*Q(i-1,:,:);
        gQ(i,:,:,2) = i*P(i-1,:,:);
        
      enddo
      return
    endsubroutine hpolyns

    function matrix2dmult(A,B) result(AB)
      real*8 :: A(:,:,:),B(:,:,:)
      real*8 :: AB(size(A,1),size(A,2),size(B,3))
      integer :: i,j,k

      do i=1,size(A,2)
        do j=1,size(B,3)
          do k=1,size(B,2)
            AB(:,i,j) = AB(:,i,j) + A(:,i,k)*B(:,k,j)
          enddo
        enddo
      enddo
    endfunction matrix2dmult
  end module basic_funcs
