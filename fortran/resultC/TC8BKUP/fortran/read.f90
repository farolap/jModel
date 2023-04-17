program read
    character(len=100) :: fileg,filer1,filer2,filevars(4)
    character(len=2) :: gi
    integer :: i,j,k
    integer :: nt,nv,ne,tmp
    integer,allocatable :: tri(:,:)
    real*8 :: var,varn,vecx(3),tg(3), &
        dudy,ulon
    real*8,allocatable :: latq(:),late(:),latp(:),pq(:,:),pp(:,:),pe(:,:),nr(:,:)
    real*8,parameter :: rad = 6.37122d6 , &
        pi = 3.141592653589793d0
    real*8,parameter :: u0 = 2d0*pi*rad/(86400d0*12)

    gi = 'g2'

    filevars = (/'div  ','ke   ','uperp','zeta '/)
!    write(fileg,'(5A)') '../../../jModel/grid/gridSCVT/',gi,'/',trim(filevars(1)),'.dat'
!    print*, trim(fileg)

    do i=2,8
        write(fileg,'(1A,I1,2A)') '../../../../grid/gridSCVT/g',i,'/', &
            'vertices.dat'
        open(unit=49,file=fileg)

        write(fileg,'(1A,I1,2A)') '../../../../grid/gridSCVT/g',i,'/', &
            'edges.dat'
        open(unit=50,file=fileg)

        write(fileg,'(1A,I1,2A)') '../../../../grid/gridSCVT/g',i,'/', &
            'united.dat'
        open(unit=51,file=fileg)

        write(fileg,'(1A,I1,2A)') '../../../../grid/gridSCVT/g',i,'/', &
            'triangles.dat'
        open(unit=52,file=fileg)

        read(49,*) nv,tmp
        read(50,*) ne,tmp
        read(51,*) ne,tmp
        read(52,*) nt,tmp
        allocate(latq(nv),late(ne),pp(nt,3),pq(nv,3),pe(ne,3),nr(ne,3),tri(nt,3))
        do j =1,nv
            read(49,*) pq(j,:)
            latq(j) = datan2(pq(j,3),dsqrt(pq(j,1)**2+pq(j,2)**2))
        enddo

        do j =1,nv
            read(49,*) pq(j,:)
            read(50,*) pe(j,:)
            read(51,*) nr(j,:),tg
            read(52,*) tri(j,:)

            latq(j) = datan2(pq(j,3),dsqrt(pq(j,1)**2+pq(j,2)**2))
            late(j) = datan2(pe(j,3),dsqrt(pe(j,1)**2+pe(j,2)**2))
        enddo

        do j=1,nt
            print*,j,tri(j,:)
            pp(j,:) = pq(tri(j,1),:)+pq(tri(j,2),:)+pq(tri(j,3),:)
        enddo

        do j=49,52
            close(j)
        enddo
        stop

        do j =1,4
            write(filer1,'(1A,I1,4A)') '../error/g',i,'/', &
                trim(filevars(j)),'.dat'
            write(filer2,'(1A,I1,4A)') './error/g',i,'/', &
                trim(filevars(j)),'.dat'
            open(unit=49,file=filer1)
            open(unit=50,file=filer2)
            print*, trim(filer1)
            do k =1,nv
                if (j==1) then
                    varn=0
                elseif(j==2) then
                    varn=(u0*dcos(latq(j)))**2*.5d0
                elseif(j==3) then
                    call convert_vec_sph2cart(1d0, 0d0, pe(j,:), vecx)
                    varn = u0*dcos(late(j))*sum(vecx*nr(j,:))
                elseif(j==4) then
                    dudy = -u0*dsin(latp(j))
                    ulon = u0*dcos(latp(j))
                    varn = -1d0/(rad*dcos(latp(j)))*(dcos(latp(j))*dudy &
                        -dsin(latp(j))*ulon);
                endif
                read(49,*) var
                write(50,*) var,varn
            enddo
            close(49)
            close(50)
        enddo
    enddo
    contains
        subroutine convert_vec_sph2cart(vlon, vlat, p, v)
            real*8, intent(in) :: p(1:3)
            real*8, intent(in) :: vlon
            real*8, intent(in) :: vlat
            real*8, intent(out) :: v(1:3)
            real*8:: r
            real*8:: rho

            r=dsqrt(p(1)**2+p(2)**2+p(3)**2)
            rho=dsqrt(p(1)**2+p(2)**2)

            if(rho==0)then
               v(1)=vlat
               v(2)=vlon
               v(3)=0
               return
            else    
               v(1)=-vlon*(p(2)/rho) - (vlat*p(1)*p(3))/(rho)
               v(2)=vlon*(p(1)/rho) - (vlat*p(2)*p(3))/(rho)
               v(3)=vlat*rho
            end if

            return    
        end subroutine convert_vec_sph2cart
endprogram read
