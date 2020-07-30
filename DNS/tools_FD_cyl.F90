module libraryGKE
  implicit none

  contains

    subroutine rhsXG(x, g, s, Re, r, dr, dz, xrhs, grhs, Nz, Nr, DsDz, DsDr)
      implicit none
      integer             :: i, j, Nz, Nr
      real*8, intent(in)  :: Re, dr, dz
      real*8              :: DxDz, DgDz, DxDr, DgDr
      real*8              :: D2xDz2, D2gDz2, D2xDr2, D2gDr2
      real*8, dimension(2:Nz-1,2:Nr-1)      :: DsDz, DsDr
      real*8, dimension(Nr), intent(in)     :: r
      real*8, dimension(2:Nz-1,2:Nr-1)      :: xrhs, grhs
      real*8, dimension(Nz,Nr), intent(in)  :: x, g, s
      do i=2,Nz-1
        do j=2,Nr-1
        ! First derivatives with respect to z
          DxDz      = (x(i+1,j)-x(i-1,j))/(2d0*dz)
          DgDz      = (g(i+1,j)-g(i-1,j))/(2d0*dz)  !!
          DsDz(i,j) = (s(i+1,j)-s(i-1,j))/(2d0*dz)
        ! First derivatives with respect to r
          DxDr      = (x(i,j+1)-x(i,j-1))/(2d0*dr)
          DgDr      = (g(i,j+1)-g(i,j-1))/(2d0*dr)  !!
          DsDr(i,j) = (s(i,j+1)-s(i,j-1))/(2d0*dr)
        ! Second derivatives with respect to z
          D2xDz2 = (x(i-1,j)-2d0*x(i,j)+x(i+1,j))/(dz**2d0)
          D2gDz2 = (g(i-1,j)-2d0*g(i,j)+g(i+1,j))/(dz**2d0)
        ! Second derivatives with respect to r
          D2xDr2 = (x(i,j-1)-2d0*x(i,j)+x(i,j+1))/(dr**2d0)
          D2gDr2 = (g(i,j-1)-2d0*g(i,j)+g(i,j+1))/(dr**2d0)
        ! Right hand side of the vorticity and angular momentum
          xrhs(i,j) = (DsDz(i,j)*DxDr-DsDr(i,j)*DxDz)/r(j)-(x(i,j)*DsDz(i,j))/(r(j)**2d0)&
                         +2d0*g(i,j)*DgDz/(r(j)**3d0)&
                         +(D2xDz2+D2xDr2+DxDr/r(j)-x(i,j)/r(j)**2d0)/Re
          grhs(i,j) = (DsDz(i,j)*DgDr-DsDr(i,j)*DgDz)/r(j)+(D2gDz2+D2gDr2-DgDr/r(j))/Re
        end do
      end do
    end subroutine rhsXG

    subroutine rhsXG2(x, g, s, Re, r, dr, dz, xrhs, grhs, Nz, Nr, DsDz, DsDr, DgDz, DgDr)
      implicit none
      integer             :: i, j, Nz, Nr
      real*8, intent(in)  :: Re, dr, dz
      real*8              :: DxDz, DxDr
      real*8              :: D2xDz2, D2gDz2, D2xDr2, D2gDr2
      real*8, dimension(2:Nz-1,2:Nr-1)      :: DsDz, DsDr, DgDz, DgDr
      real*8, dimension(Nr), intent(in)     :: r
      real*8, dimension(2:Nz-1,2:Nr-1)      :: xrhs, grhs
      real*8, dimension(Nz,Nr), intent(in)  :: x, g, s
      do i=2,Nz-1
        do j=2,Nr-1
        ! First derivatives with respect to z
          DxDz      = (x(i+1,j)-x(i-1,j))/(2d0*dz)
          DgDz(i,j) = (g(i+1,j)-g(i-1,j))/(2d0*dz)  !!
          DsDz(i,j) = (s(i+1,j)-s(i-1,j))/(2d0*dz)
        ! First derivatives with respect to r
          DxDr      = (x(i,j+1)-x(i,j-1))/(2d0*dr)
          DgDr(i,j)      = (g(i,j+1)-g(i,j-1))/(2d0*dr)  !!
          DsDr(i,j) = (s(i,j+1)-s(i,j-1))/(2d0*dr)
        ! Second derivatives with respect to z
          D2xDz2 = (x(i-1,j)-2d0*x(i,j)+x(i+1,j))/(dz**2d0)
          D2gDz2 = (g(i-1,j)-2d0*g(i,j)+g(i+1,j))/(dz**2d0)
        ! Second derivatives with respect to r
          D2xDr2 = (x(i,j-1)-2d0*x(i,j)+x(i,j+1))/(dr**2d0)
          D2gDr2 = (g(i,j-1)-2d0*g(i,j)+g(i,j+1))/(dr**2d0)
        ! Right hand side of the vorticity and angular momentum
          xrhs(i,j) = (DsDz(i,j)*DxDr-DsDr(i,j)*DxDz)/r(j)-(x(i,j)*DsDz(i,j))/(r(j)**2d0)&
                         +2d0*g(i,j)*DgDz(i,j)/(r(j)**3d0)&
                         +(D2xDz2+D2xDr2+DxDr/r(j)-x(i,j)/r(j)**2d0)/Re
          grhs(i,j) = (DsDz(i,j)*DgDr(i,j)-DsDr(i,j)*DgDz(i,j))/r(j)+(D2gDz2+D2gDr2-DgDr(i,j)/r(j))/Re
        end do
      end do
    end subroutine rhsXG2

    subroutine solve_streamfn(x, s, r, dr, dz, L, D, Nz, Nr, P, Pinv)
      implicit none
      integer                           :: i, j, Nz, Nr, info
      real*8, intent(in)                :: dr, dz
      real*8, dimension(Nr), intent(in) :: r
      real*8, dimension(Nz,Nr)          :: x, s
      real*8, dimension(2:Nz-1,2:Nr-1)  :: D, g
      real*8, dimension(2:Nz-2,2:Nr-1)  :: L
      real*8, dimension(2:Nr-1,2:Nr-1)  :: P, Pinv
    ! Compute right hand side of stream function equation multiplied by dz^2,
    ! NEED TO CHANGE IF BCs are not ZERO DIRICHLET
      do i=2,Nz-1
        do j=2,Nr-1
          s(i,j) = r(j)*x(i,j)*dz**2d0
        end do
      end do
    ! Solve the system A*Psi=rhs with A=L*D*L^T for the interior

      call dgemm('n','t',Nz-2,Nr-2,Nr-2,1.0d0,s(2:Nz-1,2:Nr-1),Nz-2,Pinv,Nr-2,0.0d0,g,Nz-2)

      do j=2,Nr-1
        call dpttrs(Nz-2,1,D(2:Nz-1,j),L(2:Nz-2,j),g(2:Nz-1,j),Nz-2,info)
        if(info.ne.0) then
          print*, 'dpttrs info: ',info, j
        end if
      end do

      call dgemm('n','t',Nz-2,Nr-2,Nr-2,1.0d0,g,Nz-2,P,Nr-2,0.0d0,s(2:Nz-1,2:Nr-1),Nz-2)

    ! The previous multiplication automatically updated the interior values of
    ! the stream function, the boundary conditions are zero in this problem so
    ! no need to update them.
    end subroutine solve_streamfn


    subroutine regSystemMatrices(Nsys,eps,delta,Rasp,f)
    implicit none
    integer,intent(in)    :: Nsys
    real*8, intent(in)    :: eps, delta,Rasp
    real*8, intent(inout) :: f(0:Nsys,1)
    real*8  :: M(0:Nsys,0:Nsys)
    real*8  :: AA,BB
    integer :: info, ipiv(Nsys+1)
!    Define parameters and matrices
    AA = 1/(1d0-Rasp**2d0)
    BB = -Rasp**2d0/(1d0-Rasp**2d0)

    M(0,0) = (1d0-eps)**3d0
    M(0,1) = (1d0-eps)**2d0
    M(0,2) = (1d0-eps)
    M(0,3) = 1
    M(1,0) = (1d0+delta)**3d0
    M(1,1) = (1d0+delta)**2d0
    M(1,2) = (1d0+delta)
    M(1,3) = 1
    M(2,0) = 3d0*(1d0-eps)**2d0
    M(2,1) = 2d0*(1d0-eps)
    M(2,2) = 1
    M(2,3) = 0
    M(3,0) = 3d0*(1d0+delta)**2d0
    M(3,1) = 2d0*(1d0+delta)
    M(3,2) = 1
    M(3,3) = 0

    f(0,1) = 1d0-eps
    f(1,1) = AA*(1d0+delta)+BB/(1d0+delta)
    f(2,1) = 1
    f(3,1) = AA-BB/((1d0+delta)**2d0)
!   Calculate coefficients of cubic spline to regularize analytical
!   solution
    call dgesv(Nsys+1,1,M,Nsys+1,ipiv,f,Nsys+1,info)
    return
    end subroutine regSystemMatrices


    subroutine infBoussinesqBC(vs,ir,r,Nr,Rasp,regOpt)
    implicit none
    integer :: i,Nsys
    parameter (Nsys=3)
    logical, intent(in) :: regOpt
    integer,intent(in)  :: ir, Nr
    real*8, intent(in)  :: r(Nr), Rasp
    real*8, intent(out) :: vs(Nr)
    real*8 :: a,b,c,d,eps,delta,AA,BB
    !NOTE: eps, delta should be intent in as well
    real*8 :: f(0:Nsys,1)
    eps   = 2d-2
    delta = 2d-2
    AA = 1/(1d0-Rasp**2d0)
    BB = -Rasp**2d0/(1d0-Rasp**2d0)
!   Calculate coefficients of cubic splie to regularize analytical
!   solution
    if (regOpt) then
      call regSystemMatrices(Nsys,eps,delta,Rasp,f)
      a = f(0,1)
      b = f(1,1)
      c = f(2,1)
      d = f(3,1)
!     Calculate analytical solution for infinite Boussinesq
      do i=0,Nr
        if (r(i) <= 1d0-eps) then
          vs(i)=r(i)
        elseif (r(i) > 1d0-eps .and. r(i) < 1d0+delta) then
          vs(i)=a*r(i)**3d0+b*r(i)**2d0+c*r(i)+d
        else
          vs(i)=AA*r(i)+BB/r(i)
        endif
      enddo
      call printRegularizationText(AA,BB,a,b,c,d,eps,delta,Rasp)
    else
      do i=0,Nr
        if (r(i) <= 1d0) then
          vs(i)=r(i)
        else
          vs(i)=AA*r(i)+BB/r(i)
        endif
      enddo
      call printAnalyticText(AA,BB,Rasp)
    endif
    return
    end subroutine infBoussinesqBC

    subroutine printRegularizationText(AA,BB,a,b,c,d,eps,delta,Rasp)
    implicit none
    real*8, intent(in) :: AA,BB,a,b,c,d,eps,delta,Rasp
    print*, ''
    print*, '======================================================='
    print*, '=== REGULARIZATION OF ANALYTICAL BOUNDARY CONDITION ==='
    print*, '======================================================='
    print *, ' '
    print *, '--------------------------'
    print *, '---- ODE Coefficients ----'
    print *, '--------------------------'
    print *, ' '
    print *, 'AA:        ', AA
    print *, 'BB:        ', BB
    print *, '--------------------------'
    print *, '--- Cubic Coefficients ---'
    print *, '--------------------------'
    print *, ' '
    print *, 'a:        ', a
    print *, 'b:        ', b
    print *, 'c:        ', c
    print *, 'd:        ', d
    print *, '--------------------------'
    print *, '------ Side Spacing ------'
    print *, '--------------------------'
    print *, ' '
    print *, 'epsilon:  ', eps
    print *, 'delta:    ', delta
    print *, 'A_r  :    ', Rasp
    return
    end subroutine printRegularizationText


    subroutine printAnalyticText(AA,BB,Rasp)
    implicit none
    real*8, intent(in) :: AA,BB,Rasp
    print*, ''
    print*, '====================================='
    print*, '=== ANALYTICAL BOUNDARY CONDITION ==='
    print*, '====================================='
    print *, ' '
    print *, '--------------------------'
    print *, '---- ODE Coefficients ----'
    print *, '--------------------------'
    print *, ' '
    print *, 'AA:        ', AA
    print *, 'BB:        ', BB
    print *, '--------------------------'
    print *, '----- Knife Placement ----'
    print *, '--------------------------'
    print *, ' '
    print *, 'A_r  :    ', Rasp
    return
    end subroutine printAnalyticText

    subroutine BndConds(x, g, s, Bo, w, beta, alpha, time, r, dr, dz,&
                           Nz, Nr, ned, ldiag, mdiag, udiag, ir,vs)
      implicit none
      integer :: i, j, Nz, Nr, ned, ir, info
      real*8  :: Bo, w, beta, alpha, time, dr, dz
      real*8, dimension(Nr)       :: r, vs
      real*8, dimension(Nz,Nr)    :: x, g, s
      real*8, dimension(2:Nr-1)   :: f
      real*8, dimension(3:Nr-1)   :: ldiag
      real*8, dimension(2:Nr-1)   :: mdiag
      real*8, dimension(2:Nr-2)   :: udiag
      real*8, dimension(ir-ned-2) :: a_int, c_int
      real*8, dimension(ir-ned-1) :: b_int
      real*8, dimension(Nr-ir-2)  :: a_ext, c_ext
      real*8, dimension(Nr-ir-1)  :: b_ext
    !---Left and Right Boundaries---!
      !The stream function is zero at the boundaries there is no need to update
      !since we only updated the interior
      x(:,1)  = 0.0d0
      g(:,1)  = 0.0d0
      x(:,Nr) = (0.5d0*s(:,Nr-2)-4d0*s(:,Nr-1))/(r(Nr)*dr**2d0)
      g(:,Nr) = 0.0d0
    !---Bottom Boundary---!
      g(1,2:Nr-1) = 0.0d0
      x(1,2:Nr-1) = (0.5d0*s(3,2:Nr-1)-4d0*s(2,2:Nr-1))/(r(2:Nr-1)*dz**2d0)
    !---Top Boundary, Contaminated Free Surface---!
      x(Nz,2:Nr-1) = (0.5d0*s(Nz-2,2:Nr-1)-4d0*s(Nz-1,2:Nr-1))/(r(2:Nr-1)*dz**2d0)
      if (Bo == 0d0) then
        do i=2,Nr-1
          g(Nz,i) = vs(i)*r(i)
        enddo
      else
        ! Calculate the right hand side of the ang momentum equation
        ! CAREFUL IF BCs ARE NOT ZERO AT THE EDGES
        f(2:ir-ned) = (-2d0*g(Nz-1,2:ir-ned)+0.5d0*g(Nz-2,2:ir-ned))/dz
        f(ir-ned) = f(ir-ned)-Bo*(1/dr**2d0-1/(2d0*r(ir-ned)*dr))*&
                                  (beta+alpha*cos(w*time))*r(ir-ned+1)**2d0
        f(ir-ned+1:ir) = (beta+alpha*cos(w*time))*r(ir-ned+1:ir)**2d0
        f(ir+1:Nr-1) = (-2d0*g(Nz-1,ir+1:Nr-1)+0.5d0*g(Nz-2,ir+1:Nr-1))/dz
        f(ir+1) = f(ir+1)-Bo*(1/dr**2d0+1/(2d0*r(ir+1)*dr))*&
                                        (beta+alpha*cos(w*time))*r(ir)**2d0
        a_int(1:ir-ned-2) = ldiag(3:ir-ned)
        b_int(1:ir-ned-1) = mdiag(2:ir-ned)
        c_int(1:ir-ned-2) = udiag(2:ir-ned-1)

        a_ext(1:Nr-ir-2) = ldiag(ir+2:Nr-1)
        b_ext(1:Nr-ir-1) = mdiag(ir+1:Nr-1)
        c_ext(1:Nr-ir-2) = udiag(ir+1:Nr-2)

        call dgtsv(ir-ned-1,1,a_int,b_int,c_int,f(2:ir-ned),ir-ned-1,info)
        if(info.ne.0) then
          print*, 'dgtsv1 INFO:   ',info
        end if

        call dgtsv(Nr-1-ir,1,a_ext,b_ext,c_ext,f(ir+1:Nr-1),Nr-1-ir,info)
        if(info.ne.0) then
          print*, 'dgtsv2 INFO:   ',info
        end if

        ! Assign values to the angular momentum
        g(Nz,2:ir-ned) = f(2:ir-ned)
        g(Nz,ir-ned+1:ir) = (beta+alpha*cos(w*time))*r(ir-ned+1:ir)**2d0
        g(Nz,ir+1:Nr-1) = f(ir+1:Nr-1)
      endif
    end subroutine BndConds

    subroutine kineticEnergy(Ek, e, g, r, DsDr, DsDz, Nz, Nr, dz, dr)
      implicit none
      integer :: Nz, Nr, j
      real*8  :: dz, dr, Ek
      real*8, dimension(Nr)     :: r
      real*8, dimension(Nz,Nr)  :: g, e
      real*8, dimension(2:Nz-1,2:Nr-1), intent(in) :: DsDz, DsDr
      ! Interior and Top Boundary
      do j=2,Nr-1
        e(2:Nz-1,j) = 0.5d0*(DsDr(:,j)**2d0+g(:,j)**2d0+DsDz(:,j)**2d0)/(r(j)**2d0)
        e(Nz,j) = 0.5d0*(g(Nz,j)/r(j))**2d0
      end do
      ! The left, right and bottom boundaries have kinetic energy zero and are
      ! not imposed because it was initialized to said value.
      ! L2 norm of the kinetic energy
        Ek = sum(sum(e,1))*dz*dr
    end subroutine kineticEnergy

    subroutine observables(Ek, Eg, Ex, ulr, ulv, ulz, ekk, egg, exx, s, g, x, r, DsDr, DsDz, DgDr, DgDz, Nz, Nr, dz, dr)
      implicit none
      integer :: Nz, Nr, j, ilz, ilr
      real*8  :: dz, dr, Ek, Eg, Ex, ulr, ulv, ulz
      real*8, dimension(Nr)     :: r
      real*8, dimension(Nz,Nr)  :: s, g, x, ekk, egg, exx
      real*8, dimension(2:Nz-1,2:Nr-1), intent(in) :: DsDr, DsDz, DgDr, DgDz

      ! -- Kinetic Energy -- !
      ! Interior and Top Boundary
      do j=2,Nr-1
        ekk(2:Nz-1,j) = 0.5d0*(DsDr(:,j)**2d0+g(:,j)**2d0+DsDz(:,j)**2d0)/(r(j)**2d0)
        ekk(Nz,j) = 0.5d0*(g(Nz,j)**2d0+((s(Nz-2,j)-4*s(Nz-1,j))/(2*dz))**2d0)/(r(j)**2d0)
      end do
      ! The left, right and bottom boundaries have kinetic energy zero and are
      ! not imposed because it was initialized to said value. The top boundary
      ! has the derivative of the stream function with respect to z at the
      ! boundary explicitly done with a backward second order finite difference.
      ! L2 norm of the kinetic energy
        Ek = sum(sum(ekk,1))*dz*dr

      ! -- Global Angular Momentum -- !
!!      do j=1,Nr
!!        egg(1:Nz,j) = g(:,j)**2d0
        egg = g**2d0
!!      end do
      ! L2 norm of the angular momentum
        Eg = sum(sum(egg,1))*dz*dr

      ! -- Enstrophy -- !
      do j=2,Nr-1
        exx(2:Nz-1,j) = (DgDz(:,j)/r(j))**2d0+x(:,j)**2d0+(DgDr(:,j)/r(j))**2d0
!        exx(Nz,j) = (DgDz(Nz,j)/r(j))**2d0+x(Nz,j)**2d0+(DgDr(Nz,j)/r(j))**2d0
      end do
        Ex = sum(sum(exx,1))*dz*dr

      ! -- Local Velocities at the point (3*Az/4,3*Ar/4) -- !
      ilr = 3*(Nr-1)/4+1 ! localized r index
      ilz = 3*(Nz-1)/4+1 ! localized z index

      ulr = -DgDz(ilz-1,ilr-1)/r(ilr) ! ilr,ilz -1 because we only have the interior of the derivatavies
      ulv = g(ilz-1,ilr-1)/r(ilr)
      ulz = DgDr(ilz-1,ilr-1)/r(ilr) ! ilr,ilz -1 because we only have the interior of the derivatavies

    end subroutine observables

    subroutine graphs(x, g, s, Re, Bo, beta, alpha, f, w, Hasp, Rasp, Nz, Nr,&
                          ned, dz, dr, dt, time, prefix, ix, init_file)
      implicit none
      integer :: Nz, Nr, ned
      integer :: i, j, init_file, ix
      real*8  :: Re, Bo, beta, alpha, f, w, Hasp, Rasp, dr, dz, dt, time
      real*8, dimension(Nz,Nr) :: x, g, s
      character*128 file_out, prefix
      file_out(1:ix)=prefix(1:ix)
      file_out(ix+1:ix+1)='_'
      write(file_out(ix+2:ix+5),'(i4.4)') init_file
      init_file=init_file+1
      open(unit=10,file=file_out(1:ix+5),form='unformatted')
      write(10) Nz,Nr,ned,dz,dr,dt,time
      write(10) Re,Bo,beta,alpha,f,w,Hasp,Rasp
      write(10) ((s(j,i),j=1,Nz),i=1,Nr),&
                ((x(j,i),j=1,Nz),i=1,Nr),&
                ((g(j,i),j=1,Nz),i=1,Nr)
      close(10)
    end subroutine graphs

end module libraryGKE
