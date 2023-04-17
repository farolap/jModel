module init_condC

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

        subroutine inicond(j,u,h,bath,div,gradb,ke,zeta,q,uperp,f,fv,mesh,resn)
            implicit none
            type(grid_structure), intent(in) :: mesh
            type(scalar_field), intent(inout):: h,div,zeta,f,fv,ke,bath,q
            type(scalar_field), intent(inout) :: u,gradb,uperp

            integer,intent(in) :: j
            integer            :: i
            real*8             :: lon,lat, lonvec(3),latvec(3), &
                uf(3),hf,bathf,ff,gradf(3),zetaf,divf,momfluxf(3), &
                pq(3),ppc(3),up(3),pe(3)
            real*8,allocatable :: hgg(:)
            character(len=100) ::resn

            h%fexact = 0d0
            zeta%fexact = 0d0
            div%fexact  = 0d0
            ke%fexact = 0d0
            gradb%fexact = 0d0
            u%fexact = 0d0
            bath%f = 0d0
            f%f  = 0d0

            if (j==2) then
                write(resn,*) '../resultCtri/TC2/'
                do i =1,mesh%nt
                    lon = mesh%tr(i)%b%lon
                    lat = mesh%tr(i)%b%lat
                    lonvec = mesh%tr(i)%b%lonvec
                    latvec = mesh%tr(i)%b%latvec
                    call TC2(uf,hf,bathf,ff,lon,lat,lonvec,latvec)
                    h%f(i) = hf
                    bath%f(i) = bathf
                enddo
                do i =1,mesh%ne
                    lon = mesh%ed(i)%c%lon
                    lat = mesh%ed(i)%c%lat
                    lonvec = mesh%ed(i)%c%lonvec
                    latvec = mesh%ed(i)%c%latvec
                    call TC2(uf,hf,bathf,ff,lon,lat,lonvec,latvec)
                    u%f(i) = sum(uf*mesh%ed(i)%nr)
                    f%f(i) = ff
                enddo
                do i =1,mesh%nv
                    lon = mesh%v(i)%lon
                    lat = mesh%v(i)%lat
                    lonvec = mesh%v(i)%c%lonvec
                    latvec = mesh%v(i)%c%latvec
                    call TC2(uf,hf,bathf,ff,lon,lat,lonvec,latvec)
                    fv%f(i) = ff
                enddo

            elseif (j==3) then
                write(resn,*) '../resultCtri/TC3/'
                do i =1,mesh%nt
                    lon = mesh%tr(i)%b%lon
                    lat = mesh%tr(i)%b%lat
                    lonvec = mesh%tr(i)%b%lonvec
                    latvec = mesh%tr(i)%b%latvec
                    call TC3(uf,hf,zetaf,divf,gradf,momfluxf,lon,lat,lonvec,latvec)
                    h%f(i) = hf
                    f%f(i) = ff
                    div%fexact(i) = divf
                enddo
                do i =1,mesh%ne
                    lon = mesh%ed(i)%c%lon
                    lat = mesh%ed(i)%c%lat
                    lonvec = mesh%ed(i)%c%lonvec
                    latvec = mesh%ed(i)%c%latvec
                    call TC3(uf,hf,zetaf,divf,gradf,momfluxf,lon,lat,lonvec,latvec)
                    u%f(i) = sum(uf*mesh%ed(i)%nr)
                    uperp%f(i) = sum(uf*mesh%ed(i)%tg)
                    f%f(i) = ff
                    gradb%fexact(i) = sum(gradf*mesh%ed(i)%nr)
                enddo
                do i =1,mesh%nv
                    lon = mesh%v(i)%lon
                    lat = mesh%v(i)%lat
                    lonvec = mesh%v(i)%c%lonvec
                    latvec = mesh%v(i)%c%latvec
                    call TC3(uf,hf,zetaf,divf,gradf,momfluxf,lon,lat,lonvec,latvec)
                    fv%f(i) = ff
                    zeta%fexact(i) = zetaf
                enddo
            elseif (j==4) then
            elseif (j==5) then
            elseif (j==6) then
            elseif (j==7) then
            elseif (j==8) then
                write(resn,*) '../resultCtri/TC8/'
                do i =1,mesh%nt
                    lon = mesh%tr(i)%c%lon
                    lat = mesh%tr(i)%c%lat
                    lonvec = mesh%tr(i)%c%lonvec
                    latvec = mesh%tr(i)%c%latvec
                    call TC8(uf,hf,zetaf,divf,gradf,momfluxf,lon,lat,lonvec,latvec)
                    h%f(i) = hf
                    ke%fexact(i) = sum(uf**2)*.5d0
                    div%fexact(i) = divf
                enddo

                do i =1,mesh%nv
                    lon = mesh%v(i)%lon
                    lat = mesh%v(i)%lat
                    lonvec = mesh%v(i)%c%lonvec
                    latvec = mesh%v(i)%c%latvec
                    call TC8(uf,hf,zetaf,divf,gradf,momfluxf,lon,lat,lonvec,latvec)
                    fv%f(i) = ff
                    zeta%fexact(i) = zetaf
                    q%fexact(i) = (zetaf+ff)/hf
                enddo
                do i =1,mesh%ne
                    lon = mesh%ed(i)%c%lon
                    lat = mesh%ed(i)%c%lat
                    lonvec = mesh%ed(i)%c%lonvec
                    latvec = mesh%ed(i)%c%latvec
                    call TC8(uf,hf,zetaf,divf,gradf,momfluxf,lon,lat,lonvec,latvec)
                    u%f(i) = sum(uf*mesh%ed(i)%nr)
                    uperp%fexact(i) = sum(uf*mesh%ed(i)%tg)
                    f%f(i) = ff
                    gradb%fexact(i) = sum(gradf*mesh%ed(i)%nr)
                enddo
                u%fexact = u%f
                h%fexact = h%f
            elseif (j==9) then
                write(resn,*) '../resultCtri/TC9/'
                saveenergy=.true.
                call TC9mfield(hgg,mesh)
                do i =1,mesh%nt
                    lat = mesh%tr(i)%c%lat
                    call TC9(uf,hf,ff,hgg,lat,mesh%tr(i)%b%p,mesh)
                    h%f(i) = hf
                    f%f(i) = ff
                enddo
                do i =1,mesh%nv
                    lat = mesh%v(i)%lat
                    call TC9(uf,hf,ff,hgg,lat,mesh%v(i)%c%p,mesh)
                    fv%f(i) = ff
                enddo
                do i =1,mesh%ne
                    lat = mesh%ed(i)%c%lat
                    call TC9(uf,hf,ff,hgg,lat,mesh%ed(i)%c%p,mesh)
                    u%f(i) = sum(uf*mesh%ed(i)%nr)
                    uperp%fexact(i) = sum(uf*mesh%ed(i)%tg)
                    f%f(i) = ff
                enddo
            elseif (j==10) then
                write(resn,*) '../resultCtri/TC10/'
                do i =1,mesh%nt
                    ppc = mesh%tr(i)%c%p
                    ppc = rotation(ppc)
                    lat = datan2(ppc(3),dsqrt(ppc(1)**2 + ppc(2)**2))
                    lonvec = mesh%tr(i)%c%lonvec
                    call TC10(hf,uf,ff,lat,lonvec)
                    h%f(i) = hf
                enddo
                do i =1,mesh%ne
                    pe = mesh%ed(i)%c%p
                    pe = rotation(pe)
                    lat = datan2(pe(3),dsqrt(pe(1)**2 + pe(2)**2))

                    lonvec = (mesh%ed(i)%c%lonvec)
                    call TC10(hf,uf,ff,lat,lonvec)
                    call convert_vec_sph2cart(norm(uf), 0d0, pe, up)
                    u%f(i) = sum(Irotation(up)*mesh%ed(i)%nr)
                    uperp%fexact(i) = sum(Irotation(up)*mesh%ed(i)%tg)
                    f%f(i) = ff
                    ! print*,sum(mesh%ed(i)%tg*mesh%ed(i)%nr)
                    ! print*,sum(mesh%edhx(i)%tg*mesh%ed(i)%nr)
                    ! print*,''
                enddo

                do i =1,mesh%nv
                    pq = mesh%v(i)%c%p
                    pq = rotation(pq)
                    lat = datan2(pq(3),dsqrt(pq(1)**2 + pq(2)**2))

                    lonvec = mesh%v(i)%c%lonvec
                    call TC10(hf,uf,ff,lat,lonvec)
                    fv%f(i) = ff
                enddo
            endif

            ! do i =1,mesh%ne
            !     call random_number(lon)
            !     call random_number(lat)
            !     uf = 30*(lon*mesh%ed(i)%c%lonvec + lat*mesh%ed(i)%c%latvec)*dcos(mesh%ed(i)%c%lat)
            !     u%f(i) = sum(uf*mesh%ed(i)%nr)
            !     uperp%fexact(i) = sum(uf*mesh%ed(i)%tg)
            ! enddo

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

end module init_condC
