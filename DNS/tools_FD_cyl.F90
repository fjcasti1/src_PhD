module tools_FD_cyl
  implicit none

  contains

    subroutine rhs_wtLt(wt, Lt, sf, Re, r, dr, dz, wt_rhs, Lt_rhs, Nz, Nr, DsfDz, DsfDr, DLtDz, DLtDr)
      implicit none
      integer             :: i, j, Nz, Nr
      real*8, intent(in)  :: Re, dr, dz
      real*8              :: DwtDz, DwtDr
      real*8              :: D2wtDz2, D2LtDz2, D2wtDr2, D2LtDr2
      real*8, dimension(2:Nz-1,2:Nr-1)      :: DsfDz, DsfDr, DLtDz, DLtDr
      real*8, dimension(Nr), intent(in)     :: r
      real*8, dimension(2:Nz-1,2:Nr-1)      :: wt_rhs, Lt_rhs
      real*8, dimension(Nz,Nr), intent(in)  :: wt, Lt, sf
      do i=2,Nz-1
        do j=2,Nr-1
        ! First derivatives with respect to z
          DwtDz      = (wt(i+1,j)-wt(i-1,j))/(2d0*dz)
          DLtDz(i,j) = (Lt(i+1,j)-Lt(i-1,j))/(2d0*dz)  !!
          DsfDz(i,j) = (sf(i+1,j)-sf(i-1,j))/(2d0*dz)
        ! First derivatives with respect to r
          DwtDr      = (wt(i,j+1)-wt(i,j-1))/(2d0*dr)
          DLtDr(i,j) = (Lt(i,j+1)-Lt(i,j-1))/(2d0*dr)  !!
          DsfDr(i,j) = (sf(i,j+1)-sf(i,j-1))/(2d0*dr)
        ! Second derivatives with respect to z
          D2wtDz2 = (wt(i-1,j)-2d0*wt(i,j)+wt(i+1,j))/(dz**2d0)
          D2LtDz2 = (Lt(i-1,j)-2d0*Lt(i,j)+Lt(i+1,j))/(dz**2d0)
        ! Second derivatives with respect to r
          D2wtDr2 = (wt(i,j-1)-2d0*wt(i,j)+wt(i,j+1))/(dr**2d0)
          D2LtDr2 = (Lt(i,j-1)-2d0*Lt(i,j)+Lt(i,j+1))/(dr**2d0)
        ! Right hand side of the vorticity and angular momentum
          wt_rhs(i,j) = (DsfDz(i,j)*DwtDr-DsfDr(i,j)*DwtDz)/r(j)-(wt(i,j)*DsfDz(i,j))/(r(j)**2d0)&
                         +2d0*Lt(i,j)*DLtDz(i,j)/(r(j)**3d0)&
                         +(D2wtDz2+D2wtDr2+DwtDr/r(j)-wt(i,j)/r(j)**2d0)/Re
          Lt_rhs(i,j) = (DsfDz(i,j)*DLtDr(i,j)-DsfDr(i,j)*DLtDz(i,j))/r(j)+(D2LtDz2+D2LtDr2-DLtDr(i,j)/r(j))/Re
        end do
      end do
    end subroutine rhs_wtLt

    subroutine rhs_c(c, sf, Pe, r, dz, dr, Nz, Nr, c_rhs)
      implicit none
      integer              :: j
      integer, intent(in)  :: Nz, Nr
      real*8 , intent(in)  :: Pe, dz, dr
      real*8               :: Dr_cDsfDz, DcDr, D2cDr2
      real*8 , dimension(Nr)                  :: DsfDz
      real*8 , dimension(Nr)    , intent(in)  :: r, c
      real*8 , dimension(Nz,Nr) , intent(in)  :: sf
      real*8 , dimension(2:Nr-1), intent(out) :: c_rhs
      ! Although DsfDz was computed in other subroutines, it is done only for
      ! the interior. Here we need it for the surface. It will be done
      ! internally and not saved since we don't need it elsewhere
      do j=1,Nr
        DsfDz(j)  = (sf(Nz-2,j)-4d0*sf(Nz-1,j))/(2d0*dz)
      end do
      ! Now we use DsfDz to compute Dr_cDsfDz, and later the RHS for the
      ! concentration equation at each r(j)
      do j=2,Nr-1
      ! First derivatives with respect to r
        Dr_cDsfDz = (c(j+1)*DsfDz(j+1)-c(j-1)*DsfDz(j-1))/(2d0*dr)
        DcDr      = (c(j+1)-c(j-1))/(2d0*dr)
      ! Second derivatives with respect to r
        D2cDr2    = (c(j-1)-2d0*c(j)+c(j+1))/(dr**2d0)
      ! Right hand side of the vorticity and angular momentum
        c_rhs(j)  = Dr_cDsfDz/r(j)+(D2cDr2+DcDr/r(j))/Pe
      end do
    end subroutine rhs_c

    subroutine solve_streamfn(wt, sf, r, dz, L, D, Nz, Nr, P, Pinv)
      implicit none
      integer                           :: i, j, Nz, Nr, info
      real*8, intent(in)                :: dz
      real*8, dimension(Nr), intent(in) :: r
      real*8, dimension(Nz,Nr)          :: wt, sf
      real*8, dimension(2:Nz-1,2:Nr-1)  :: D, Lt
      real*8, dimension(2:Nz-2,2:Nr-1)  :: L
      real*8, dimension(2:Nr-1,2:Nr-1)  :: P, Pinv
    ! Compute right hand side of stream function equation multiplied by dz^2,
    ! NEED TO CHANGE IF BCs are not ZERO DIRICHLET
      do i=2,Nz-1
        do j=2,Nr-1
          sf(i,j) = r(j)*wt(i,j)*dz**2d0
        end do
      end do
    ! Solve the system A*Psi=rhs with A=L*D*L^T for the interior

      call dgemm('n','t',Nz-2,Nr-2,Nr-2,1d0,sf(2:Nz-1,2:Nr-1),Nz-2,Pinv,Nr-2,0d0,Lt,Nz-2)

      do j=2,Nr-1
        call dpttrs(Nz-2,1,D(2:Nz-1,j),L(2:Nz-2,j),Lt(2:Nz-1,j),Nz-2,info)
        if(info.ne.0) then
          print*, 'dpttrs info: ',info, j
        end if
      end do

      call dgemm('n','t',Nz-2,Nr-2,Nr-2,1d0,Lt,Nz-2,P,Nr-2,0d0,sf(2:Nz-1,2:Nr-1),Nz-2)

    ! The previous multiplication automatically updated the interior values of
    ! the stream function, the boundary conditions are zero in this problem so
    ! no need to update them.
    end subroutine solve_streamfn
!   subroutine regSystemMatrices(Nsys,eps,delta,Rasp,f)

    subroutine regSystemMatrices(Nsys,eta,eps,delta,f)
      implicit none
      integer,intent(in)    :: Nsys
      real*8, intent(in)    :: eps, delta,eta
      real*8, intent(inout) :: f(0:Nsys,1)
      real*8  :: M(0:Nsys,0:Nsys)
      real*8  :: AA,BB
      integer :: info, ipiv(Nsys+1)
!      Define parameters and matrices
!      AA = 1/(1d0-Rasp**2d0)
!      BB = -Rasp**2d0/(1d0-Rasp**2d0)
      AA = eta**2d0/(eta**2d0-1d0)
      BB = -AA

      M(0,0) = (eta-eps)**3d0
      M(0,1) = (eta-eps)**2d0
      M(0,2) = (eta-eps)
      M(0,3) = 1d0
      M(1,0) = (eta+delta)**3d0
      M(1,1) = (eta+delta)**2d0
      M(1,2) = (eta+delta)
      M(1,3) = 1d0
      M(2,0) = 3d0*(eta-eps)**2d0
      M(2,1) = 2d0*(eta-eps)
      M(2,2) = 1d0
      M(2,3) = 0d0
      M(3,0) = 3d0*(eta+delta)**2d0
      M(3,1) = 2d0*(eta+delta)
      M(3,2) = 1d0
      M(3,3) = 0d0

      f(0,1) = eta-eps
      f(1,1) = AA*(eta+delta)+BB/(eta+delta)
      f(2,1) = 1d0
      f(3,1) = AA-BB/((eta+delta)**2d0)
!     Calculate coefficients of cubic spline to regularize analytical
!     solution
      call dgesv(Nsys+1,1,M,Nsys+1,ipiv,f,Nsys+1,info)
      return
    end subroutine regSystemMatrices


!    subroutine infBoussinesqBC(vs,r,Nr,Rasp,regOpt)
    subroutine infBoussinesqBC(vs,eta,r,Nr,regOpt)
    implicit none
    integer :: i,Nsys
    parameter (Nsys=3)
    logical, intent(in) :: regOpt
    integer,intent(in)  :: Nr
    real*8, intent(in)  :: r(Nr), eta
    real*8, intent(out) :: vs(Nr)
    real*8 :: a,b,c,d,eps,delta,AA,BB
    !NOTE: eps, delta should be intent in as well
    real*8 :: f(0:Nsys,1)
    eps   = 2d-2
    delta = 2d-2
!    AA = 1/(1d0-Rasp**2d0)
!    BB = -Rasp**2d0/(1d0-Rasp**2d0)
    AA = eta**2d0/(eta**2d0-1d0)
    BB = -AA
!   Calculate coefficients of cubic splie to regularize analytical
!   solution
    if (regOpt) then
!      call regSystemMatrices(Nsys,eps,delta,Rasp,f)
      call regSystemMatrices(Nsys,eta,eps,delta,f)
      a = f(0,1)
      b = f(1,1)
      c = f(2,1)
      d = f(3,1)
!     Calculate analytical solution for infinite Boussinesq
      do i=0,Nr
        if (r(i) <= eta-eps) then
          vs(i)=r(i)
        elseif (r(i) > eta-eps .and. r(i) < eta+delta) then
          vs(i)=a*r(i)**3d0+b*r(i)**2d0+c*r(i)+d
        else
          vs(i)=AA*r(i)+BB/r(i)
        endif
      enddo
      call printRegularizationText(AA,BB,a,b,c,d,eps,delta,eta)
    else
      do i=0,Nr
        if (r(i) <= eta) then
          vs(i)=r(i)
        else
          vs(i)=AA*r(i)+BB/r(i)
        endif
      enddo
!      call printAnalyticText(AA,BB,Rasp)
      call printAnalyticText(AA,BB,eta)
    endif
    return
    end subroutine infBoussinesqBC

    subroutine printRegularizationText(AA,BB,a,b,c,d,eps,delta,eta)
    implicit none
    real*8, intent(in) :: AA,BB,a,b,c,d,eps,delta,eta
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
    print *, 'eta  :    ', eta
    return
    end subroutine printRegularizationText


    subroutine printAnalyticText(AA,BB,eta)
    implicit none
    real*8, intent(in) :: AA,BB,eta
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
    print *, 'eta  :    ', eta
    return
    end subroutine printAnalyticText

    subroutine BC_kedgeTop(wt, Lt, sf, Bo, wf, Ro, time, r, dr, dz,&
                           Nz, Nr, ned, ldiag, mdiag, udiag, ir, vs)
      implicit none
      integer :: i, Nz, Nr, ned, ir, info
      real*8  :: Bo, wf, Ro, time, dr, dz
      real*8, dimension(Nr)       :: r, vs
      real*8, dimension(Nz,Nr)    :: wt, Lt, sf
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
      wt(:,1)  = 0d0
      Lt(:,1)  = 0d0
      wt(:,Nr) = (0.5d0*sf(:,Nr-2)-4d0*sf(:,Nr-1))/(r(Nr)*dr**2d0)
      Lt(:,Nr) = 0d0
    !---Bottom Boundary---!
      Lt(1,2:Nr-1) = 0d0
      wt(1,2:Nr-1) = (0.5d0*sf(3,2:Nr-1)-4d0*sf(2,2:Nr-1))/(r(2:Nr-1)*dz**2d0)
    !---Top Boundary, Contaminated Free Surface---!
      wt(Nz,2:Nr-1) = (0.5d0*sf(Nz-2,2:Nr-1)-4d0*sf(Nz-1,2:Nr-1))/(r(2:Nr-1)*dz**2d0)
      if (Bo == 0d0) then
        do i=2,Nr-1
          Lt(Nz,i) = vs(i)*r(i)
        enddo
      else
        ! Calculate the right hand side of the ang momentum equation
        ! CAREFUL IF BCs ARE NOT ZERO AT THE EDGES
        f(2:ir-ned) = (-2d0*Lt(Nz-1,2:ir-ned)+0.5d0*Lt(Nz-2,2:ir-ned))/dz
        f(ir-ned) = f(ir-ned)-Bo*(1/dr**2d0-1/(2d0*r(ir-ned)*dr))*&
                                  (1d0+Ro*dsin(wf*time))*r(ir-ned+1)**2d0
        f(ir-ned+1:ir) = (1d0+Ro*dsin(wf*time))*r(ir-ned+1:ir)**2d0
        f(ir+1:Nr-1) = (-2d0*Lt(Nz-1,ir+1:Nr-1)+0.5d0*Lt(Nz-2,ir+1:Nr-1))/dz
        f(ir+1) = f(ir+1)-Bo*(1/dr**2d0+1/(2d0*r(ir+1)*dr))*&
                                        (1d0+Ro*dsin(wf*time))*r(ir)**2d0
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
        Lt(Nz,2:ir-ned) = f(2:ir-ned)
        Lt(Nz,ir-ned+1:ir) = (1d0+Ro*dsin(wf*time))*r(ir-ned+1:ir)**2d0
        Lt(Nz,ir+1:Nr-1) = f(ir+1:Nr-1)
      endif
    end subroutine BC_kedgeTop

    subroutine solve_concentration(c,sf,Pe,r,dz,dr,dt,Nz,Nr,c_tmp,c_rhs)
      implicit none
      integer :: Nz, Nr
      real*8  :: Pe, dz, dr, dt
      real*8, dimension(Nr) :: r, c, c_tmp, c_rhs
      real*8, dimension(Nz,Nr) :: sf
    !-- Solves the concentration with a predictor-corrector method

      !First RK2 step
      call rhs_c(c,sf,Pe,r,dz,dr,Nz,Nr,c_rhs)
      c_tmp(2:Nr-1) = c(2:Nr-1) + dt*c_rhs(2:Nr-1)
      c_tmp(1)  = c_tmp(2)
      c_tmp(Nr) = c_tmp(Nr-1)

      !Second RK2 step
      call rhs_c(c_tmp,sf,Pe,r,dz,dr,Nz,Nr,c_rhs)
      c(2:Nr-1) = 0.5d0*(c(2:Nr-1) + c_tmp(2:Nr-1) + dt*c_rhs(2:Nr-1))
      c(1)  = c(2)
      c(Nr) = c(Nr-1)
    end subroutine solve_concentration

    subroutine BC_freeSurfTop(wt, Lt, sf, c, Ca, wf, Ro, time, r, dr, dz,&
                           Nz, Nr)
      implicit none
      integer :: i, Nz, Nr
      real*8  :: Ca, wf, Ro, time, dr, dz, den
      real*8, dimension(Nr)       :: r, c, sigma, mu_s, k_s
      real*8, dimension(Nz,Nr)    :: wt, Lt, sf
      real*8  :: DsigmaDr, DmuDr, DmukDr, DsfDz, D2sfDrDz, D3sfDr2Dz

      !The stream function is zero at the boundaries there is no need to update
      !since we only updated the interior

    !--- Left Boundary, symmetry axis ---!
      Lt(:,1)  = 0d0
      wt(:,1)  = 0d0

    !--- Right Boundary, side wall (no-slip) ---!
      Lt(:,Nr) = 0d0
      wt(:,Nr) = (0.5d0*sf(:,Nr-2)-4d0*sf(:,Nr-1))/(r(Nr)*dr**2d0)
      ! NOTE: Second order approximation of D2sfDr2 only with the assumption
      ! that stream function, sf, and its first radial derivative, DsfDr, are
      ! zero. This is true given no-slip boundary conditions

    !--- Bottom Boundary, rotating no-slip ---!
      Lt(1,2:Nr-1) = (1+Ro*dsin(wf*time))*r(2:Nr-1)**2d0
      wt(1,2:Nr-1) = (0.5d0*sf(3,2:Nr-1)-4d0*sf(2,2:Nr-1))/(r(2:Nr-1)*dz**2d0)
      ! NOTE: Second order approximation of D2sfDz2 only with the assumption
      ! that stream function, sf, and its first axial derivative, DsfDz, are
      ! zero. This is true given no-slip boundary conditions

    !--- Top Boundary, Contaminated Free Surface ---!
      call stateEq_surfTension(sigma, c, Nr,'tanh')
      call stateEq_surfShearVisc(mu_s, c, Nr)
!      call stateEq_surfDilatVisc(k_s, mu_s, Nr)
      k_s = 10d0*mu_s
!! ----- Sanity check 1
!      sigma = 0d0
!      mu_s  = 0d0
!      k_s   = 0d0
!! ----- Sanity check 2
!      mu_s  = 0d0
!      k_s   = 0d0
! ----- Sanity check 3
!     k_s   = 0d0
!! ----- Sanity check 4
!
      do i=2,Nr-1
        !-- Condition for angular momentum --!
        ! Denominator first
        den = 3d0*dr/(2d0*dz)+2d0*mu_s(i)/dr+(mu_s(i+1)-mu_s(i-1))/r(i)
        ! Angular momentum
        Lt(Nz,i) = (mu_s(i)*(Lt(Nz,i+1)+Lt(Nz,i-1))/dr &
                  -mu_s(i)*(Lt(Nz,i+1)-Lt(Nz,i-1))/(2d0*r(i)) &
                  +(mu_s(i+1)-mu_s(i-1))*(Lt(Nz,i+1)-Lt(Nz,i-1))/(4d0*dr)&
                  -(Lt(Nz-2,i)-4d0*Lt(Nz-1,i))*dr/(2d0*dz))/den
        !-- Condition for azimuthal vorticity --!
        ! Derivatives needed
        DsigmaDr  = (sigma(i+1)-sigma(i-1))/(2d0*dr)
        DmuDr     = (mu_s(i+1)-mu_s(i-1))/(2d0*dr)
        DmukDr    = (mu_s(i+1)+k_s(i+1)-mu_s(i-1)-k_s(i-1))/(2d0*dr)
        DsfDz     = (sf(Nz-2,i)-4d0*sf(Nz-1,i))/(2d0*dz)
        D2sfDrDz  = (sf(Nz-2,i+1)-4d0*sf(Nz-1,i+1) &
                    -sf(Nz-2,i-1)+4d0*sf(Nz-1,i-1))/(4d0*dr*dr)
        D3sfDr2Dz = (sf(Nz-2,i+1)-4d0*sf(Nz-1,i+1) &
                -2d0*sf(Nz-2,i)  +8d0*sf(Nz-1,i  ) &
                    +sf(Nz-2,i-1)-4d0*sf(Nz-1,i-1))/(2d0*dr**2d0*dr)
        ! Azimuthal vorticity
        wt(Nz,i)  = DsigmaDr/Ca &
                   +(mu_s(i)+k_s(i))*(D2sfDrDz/(r(i)**2d0)-D3sfDr2Dz/r(i)) &
                   -D2sfDrDz*DmukDr/r(i) + 2d0*DsfDz*DmuDr/(r(i)**2d0)

      end do
    end subroutine BC_freeSurfTop

    subroutine stateEq_surfTension(sigma, c, Nr, option)
      implicit none
      character(len=*), intent(in) :: option
      integer, intent(in)  :: Nr
      real*8 , dimension(Nr) , intent(in)    :: c
      real*8 , dimension(Nr) , intent(inout) :: sigma
      real*8 :: sigma0, a0, a1, a2, a3, a4, a5, a6
      if (option == 'tanh' ) then
        ! Equation of State fitting from:
        ! Lopez & Hirsa (2000)
        sigma0 = 66
        a0     = 6.3
        a1     = 6.2
        sigma  = 1 + (a0/sigma0)*dtanh(a1*(1-c))
      elseif (option == 'exp') then
        ! Equation of State fitting from:
        ! Hirsa, Lopez & Miraghaie (2001)
        sigma0 =  72.4d0
        a0     =  1.108d0
        a1     =  32.37d0
        a2     =  20.11d0
        a3     =  97.04d0
        a4     = -45.90d0
        a5     =  sigma0
        a6     = -00.15d0
        sigma  = ((a2+a3*c+a4*c**2d0)/(1d0+dexp(a1*(a0-c))) + &
                      (a5+a6*c**2d0)/(1d0+dexp(a1*(c-a0))))/sigma0
      endif

    end subroutine stateEq_surfTension

    subroutine stateEq_surfShearVisc(mu_s, c, Nr)
      implicit none
      integer, intent(in)  :: Nr
      real*8 , dimension(Nr) , intent(in)    :: c
      real*8 , dimension(Nr) , intent(inout) :: mu_s
      real*8 :: a0, a1
      ! Equation of State fitting from:
      ! Lopez & Hirsa (2000)
      a0   = 1.15d-3
      a1   = 2.6d0
      mu_s = a0*(dexp(a1*c)-1)

    end subroutine stateEq_surfShearVisc

    subroutine kineticEnergy(Ek, e, Lt, r, DsfDr, DsfDz, Nz, Nr, dz, dr)
      implicit none
      integer :: Nz, Nr, j
      real*8  :: dz, dr, Ek
      real*8, dimension(Nr)     :: r
      real*8, dimension(Nz,Nr)  :: Lt, e
      real*8, dimension(2:Nz-1,2:Nr-1), intent(in) :: DsfDz, DsfDr
      ! Interior and Top Boundary
      do j=2,Nr-1
        e(2:Nz-1,j) = 0.5d0*(DsfDr(:,j)**2d0+Lt(:,j)**2d0+DsfDz(:,j)**2d0)/(r(j)**2d0)
        e(Nz,j) = 0.5d0*(Lt(Nz,j)/r(j))**2d0
      end do
      ! The left, right and bottom boundaries have kinetic energy zero and are
      ! not imposed because it was initialized to said value.
      ! L2 norm of the kinetic energy
        Ek = sum(sum(e,1))*dz*dr
    end subroutine kineticEnergy

    subroutine observables(Ek, Eg, Ex, ulr, ulv, ulz, ekk, egg, exx, sf, Lt, wt, r, DsfDr, DsfDz, DLtDr, DLtDz, Nz, Nr, dz, dr)
      implicit none
      integer :: Nz, Nr, j, ilz, ilr
      real*8  :: dz, dr, Ek, Eg, Ex, ulr, ulv, ulz
      real*8, dimension(Nr)     :: r
      real*8, dimension(Nz,Nr)  :: sf, Lt, wt, ekk, egg, exx
      real*8, dimension(2:Nz-1,2:Nr-1), intent(in) :: DsfDr, DsfDz, DLtDr, DLtDz

      ! -- Kinetic Energy -- !
      ! Interior and Top Boundary
      do j=2,Nr-1
        ekk(2:Nz-1,j) = 0.5d0*(DsfDr(:,j)**2d0+Lt(:,j)**2d0+DsfDz(:,j)**2d0)/(r(j)**2d0)
        ekk(Nz,j) = 0.5d0*(Lt(Nz,j)**2d0+((sf(Nz-2,j)-4*sf(Nz-1,j))/(2*dz))**2d0)/(r(j)**2d0)
      end do
      ! The left, right and bottom boundaries have kinetic energy zero and are
      ! not imposed because it was initialized to said value. The top boundary
      ! has the derivative of the stream function with respect to z at the
      ! boundary explicitly done with a backward second order finite difference.
      ! L2 norm of the kinetic energy
        Ek = sum(sum(ekk,1))*dz*dr

      ! -- Global Angular Momentum -- !
        egg = Lt**2d0
      ! L2 norm of the angular momentum
        Eg = sum(sum(egg,1))*dz*dr

      ! -- Enstrophy -- !
      do j=2,Nr-1
        exx(2:Nz-1,j) = (DLtDz(:,j)/r(j))**2d0+wt(:,j)**2d0+(DLtDr(:,j)/r(j))**2d0
      end do
        Ex = sum(sum(exx,1))*dz*dr

      ! -- Local Velocities at the point (3*Az/4,3*Ar/4) -- !
      ilr = 3*(Nr-1)/4+1 ! localized r index
      ilz = 3*(Nz-1)/4+1 ! localized z index

      ulr = -DLtDz(ilz-1,ilr-1)/r(ilr) ! ilr,ilz -1 because we only have the interior of the derivatavies
      ulv = Lt(ilz-1,ilr-1)/r(ilr)
      ulz = DLtDr(ilz-1,ilr-1)/r(ilr) ! ilr,ilz -1 because we only have the interior of the derivatavies

    end subroutine observables

    subroutine graphs_kedgeTop(wt, Lt, sf, Re, Bo, Ro, f, wf, Gama, eta, Nz, Nr,&
                          ned, dz, dr, dt, time, prefix, ix, init_file)
      implicit none
      integer :: Nz, Nr, ned
      integer :: i, j, init_file, ix
      real*8  :: Re, Bo, Ro, f, wf, Gama, eta, dr, dz, dt, time
      real*8, dimension(Nz,Nr) :: wt, Lt, sf
      character*128 file_out, prefix
      file_out(1:ix)=prefix(1:ix)
      file_out(ix+1:ix+1)='_'
      write(file_out(ix+2:ix+5),'(i4.4)') init_file
      init_file=init_file+1
      open(unit=10,file=file_out(1:ix+5),form='unformatted')
      write(10) Nz,Nr,ned,dz,dr,dt,time
      write(10) Re,Bo,Ro,f,wf,Gama,eta
      write(10) ((sf(j,i),j=1,Nz),i=1,Nr),&
                ((wt(j,i),j=1,Nz),i=1,Nr),&
                ((Lt(j,i),j=1,Nz),i=1,Nr)
      close(10)
    end subroutine graphs_kedgeTop

    subroutine graphs_freeSurfTop(wt, Lt, sf, Re, Pe, Ca, Ro, wf, Gama, Nz, Nr,&
                          dz, dr, dt, time, prefix, ix, init_file)
      implicit none
      integer :: Nz, Nr
      integer :: i, j, init_file, ix
      real*8  :: Re, Pe, Ca, Ro, wf, Gama, dr, dz, dt, time
      real*8, dimension(Nz,Nr) :: wt, Lt, sf
      character*128 file_out, prefix
      file_out(1:ix)=prefix(1:ix)
      file_out(ix+1:ix+1)='_'
      write(file_out(ix+2:ix+5),'(i4.4)') init_file
      init_file=init_file+1
      open(unit=10,file=file_out(1:ix+5),form='unformatted')
      write(10) Nz,Nr,dz,dr,dt,time
      write(10) Re,Pe,Ca,Ro,wf,Gama
      write(10) ((sf(j,i),j=1,Nz),i=1,Nr),&
                ((wt(j,i),j=1,Nz),i=1,Nr),&
                ((Lt(j,i),j=1,Nz),i=1,Nr)
      close(10)
    end subroutine graphs_freeSurfTop

end module tools_FD_cyl
