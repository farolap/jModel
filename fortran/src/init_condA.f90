module init_condA
    !Use main grid data structures
    use constants
    use datastruct, only: &
    grid_structure, &
        scalar_field, &
        vectorinterpol_methods, &
        vector_field_cart

    use smeshpack
    use basic_funcs
    use init_cond

    contains
        subroutine inicond(j,u,h,htr,bath,div,gradb,zeta, &
                f,mesh,resn)
            implicit none
            type(grid_structure), intent(in) :: mesh
            type(scalar_field), intent(inout):: h,htr,div,zeta,f,bath
            type(vector_field_cart), intent(inout) :: u,gradb

            integer,intent(in) :: j
            integer            :: i
            real*8             :: lon,lat, lonvec(3),latvec(3), &
                uf(3),hf,bathf,ff,gradf(3),zetaf,divf,momfluxf(3), &
                pq(3),up(3)
            real*8,allocatable :: hgg(:)
            character(len=100) ::resn

            h%fexact = 0d0
            zeta%fexact  = 0d0
            div%fexact   = 0d0
            do i=1,3
                gradb%pexact%v(i) = 0d0
                u%pexact%v(i) = 0d0
            enddo
            bath%f = 0d0
            f%f  = 0d0

            if (j==2) then
                write(resn,*) '../resultA/TC2/'
                do i =1,mesh%nv
                    lon = mesh%v(i)%lon
                    lat = mesh%v(i)%lat
                    lonvec = mesh%v(i)%c%lonvec
                    latvec = mesh%v(i)%c%latvec
                    call TC2(uf,hf,bathf,ff,lon,lat,lonvec,latvec)
                    h%f(i) = hf
                    bath%f(i) = bathf
                    u%p(i)%v = uf
                    f%f(i) = ff
                enddo

            elseif (j==3) then
                write(resn,*) '../resultA/TC3/'
                do i =1,mesh%nv
                    lon = mesh%v(i)%lon
                    lat = mesh%v(i)%lat
                    lonvec = mesh%v(i)%c%lonvec
                    latvec = mesh%v(i)%c%latvec
                    call TC3(uf,hf,zetaf,divf,gradf,momfluxf,lon,lat,lonvec,latvec)
                    ! fv%f(i) = ff
                    div%f = divf
                    u%p(i)%v = uf
                    f%f(i) = ff
                    gradb%p(i)%v = gradf
                enddo
            elseif (j==4) then
            elseif (j==5) then
            elseif (j==6) then
            elseif (j==7) then
            elseif (j==8) then
                write(resn,*) '../resultA/TC8/'
                do i =1,mesh%nt
                    lon = mesh%tr(i)%b%lon
                    lat = mesh%tr(i)%b%lat
                    ! lonvec = mesh%tr(i)%b%lonvec
                    ! latvec = mesh%tr(i)%b%latvec
                    call TC8(uf,hf,zetaf,divf,gradf,momfluxf,lon,lat,lonvec,latvec)
                    htr%fexact(i) = hf
                enddo
                do i =1,mesh%nv
                    lon = mesh%v(i)%lon
                    lat = mesh%v(i)%lat
                    lonvec = mesh%v(i)%c%lonvec
                    latvec = mesh%v(i)%c%latvec
                    call TC8(uf,hf,zetaf,divf,gradf,momfluxf,lon,lat,lonvec,latvec)
                    h%f(i) = hf
                    div%fexact(i) = divf
                    zeta%fexact(i) = zetaf
                    u%p(i)%v = uf
                    f%f(i) = ff
                    gradb%pexact(i)%v = gradf
                enddo
            elseif (j==9) then
                write(resn,*) '../resultA/TC9/'
                saveenergy=.true.
                call TC9mfield(hgg,mesh)
                do i =1,mesh%nv
                    lat = mesh%v(i)%lat
                    call TC9(uf,hf,ff,hgg,lat,mesh%v(i)%c%p,mesh)
                    h%f(i) = hf
                    u%p(i)%v = uf
                    f%f(i) = ff
                enddo
            elseif (j==10) then
                write(resn,*) '../resultA/TC10/'
                do i =1,mesh%nv
                    pq = mesh%v(i)%p
                    pq = rotation(pq)

                    lat = datan2(pq(3),dsqrt(pq(1)**2 + pq(2)**2))
                    lonvec = mesh%v(i)%c%lonvec
                    
                    call TC10(hf,uf,ff,lat,lonvec)
                    call convert_vec_sph2cart(norm(uf), 0d0, pq, up)
                    u%p(i)%v = Irotation(up)

                    f%f(i) = ff
                    h%f(i) = hf
                enddo
            endif

        end subroutine inicond

        function rotation(pq) result(pqr)
            real*8 :: pq(3)
            real*8 :: pqr(3)
            real*8 :: Rmat(3,3)
            real*8,parameter :: lat=1*pi/180,&
                lon=3*pi/180

            ! lon = datan2(pq(2),pq(1))
            ! lat = datan2(pq(3),dsqrt(pq(1)**2+pq(2)**2))

            Rmat(1:3,1:3)=0._r8
            Rmat(1,1)=dsin(lat)*dcos(lon)
            Rmat(1,2)=dsin(lat)*dsin(lon)
            Rmat(1,3)=-dcos(lat)
            Rmat(2,1)=-dsin(lon)
            Rmat(2,2)=dcos(lon)
            Rmat(3,1)=dcos(lat)*dcos(lon)
            Rmat(3,2)=dcos(lat)*dsin(lon)
            Rmat(3,3)=dsin(lat)

            ! Rmat = roty1(pi*.5d0)


            pqr = matmul(Rmat,pq)

        end function rotation

        function Irotation(pqr) result(pq)
            real*8 :: pq(3)
            real*8 :: pqr(3)
            real*8 :: Rmat(3,3),RmatT(3,3)
            real*8,parameter :: lat=1*pi/180,&
                lon=3*pi/180

                Rmat(1:3,1:3)=0._r8
                Rmat(1,1)=dsin(lat)*dcos(lon)
                Rmat(1,2)=dsin(lat)*dsin(lon)
                Rmat(1,3)=-dcos(lat)
                Rmat(2,1)=-dsin(lon)
                Rmat(2,2)=dcos(lon)
                Rmat(3,1)=dcos(lat)*dcos(lon)
                Rmat(3,2)=dcos(lat)*dsin(lon)
                Rmat(3,3)=dsin(lat)

            RmatT = transpose(Rmat)
            pq = matmul(RmatT,pqr)

        end function Irotation

end module init_condA