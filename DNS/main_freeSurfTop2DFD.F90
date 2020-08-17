program main_freeSurfaceTop
  use tools_FD_cyl
  implicit none
  logical  :: existe, regOpt
  integer  :: Nz, Nz1, Nr, Nr1, Nsteps, nsaves
  integer  :: i, j, m, ix, irestart
  integer  :: igraph, itseries, ibegin, init_file, info
  real*8   :: Ca, Ca1, Re, Re1, Pe, Pe1, Ro, Ro1, wf, wf1, simTU, NT, NtsT, T
  real*8   :: ALT, RAD, Gama, Gama1, dr, dr1, dz, dz1
  real*8   :: time, oldtime, dt, dt1
  real*8   :: Ek, Eg, Ex, ulr, ulv, ulz
  integer, dimension(:),   allocatable :: ipiv
  real*8,  dimension(:),   allocatable :: ldiag, mdiag, udiag
  real*8,  dimension(:),   allocatable :: eigRe, eigIm
  real*8,  dimension(:),   allocatable :: r, c, c_tmp, c_rhs
  real*8,  dimension(:,:), allocatable :: wt, Lt, sf, wt_tmp, Lt_tmp
  real*8,  dimension(:,:), allocatable :: wt_rhs, Lt_rhs, DsfDz, DsfDr, DLtDz, DLtDr, ekk, egg, exx
  real*8,  dimension(:,:), allocatable :: L, D, B, P, Paux, Pinv, dumeig
  character*128 :: prefix, restart, tfile
!  character*128 :: vfile ! Write surface velocity
  real*8, parameter :: PI = 4.d0*atan(1.0d0)

! Read parameters
  read*, prefix    ! prefix for filenames
  ix=index(prefix//' ',' ')-1
  read*, restart   ! name of restart file
  irestart=index(restart//' ',' ')-1
  read*, Re
  read*, Pe
  read*, Ca
  read*, Ro
  read*, wf
  read*, Gama
  read*, Nr
  read*, Nz
  read*, dt         ! Default dt if wf = 0
  read*, NtsT       ! Number of Time Steps per Period
  read*, NT         ! Number of Periods or Number of Time Units (if wf=0)
  read*, nsaves
  read*, regOpt
  read*, itseries
  read*, init_file
  read*, ibegin     ! controls the start/restart process

! ===========
! = Scaling =
! ===========
! ALT = H/Lc
! RAD = R/Lc
! Characteristic Length: Lc = R
! Gamma = H/R

  ALT = Gama
  RAD = 1d0

  dr = RAD/(Nr-1)
  dz = ALT/(Nz-1)

  if (Ro.gt.0d0.AND.wf.gt.0d0) then  ! Forced system
    T  = 2.d0*PI/wf
    dt = T/NtsT
  elseif (Ro.eq.0d0.AND.wf.gt.0d0) then ! Unforced system, given response frequency
    T  = 2.d0*PI/wf
    dt = T/NtsT
!    Nsteps = CEILING(NT/dt)
  endif
  Nsteps = NT*NtsT
  igraph = Nsteps/nsaves
  simTU  = dt*Nsteps

! Print general case info
  print *, '======================================================='
  print *, '================ INPUT PARAMETERS READ ================'
  print *, '======================================================='
  print *, 'prefix:       ', prefix
  print *, 'restart:      ', restart
  print *, 'Re:           ', Re
  print *, 'Pe:           ', Pe
  print *, 'Ca:           ', Ca
  print *, 'Ro:           ', Ro
  print *, 'wf:           ', wf
  print *, 'Gamma:        ', Gama
  print *, 'Nr:           ', Nr
  print *, 'Nz:           ', Nz
  print *, 'NtsT:         ', NtsT
  print *, 'NT:           ', NT
  print *, 'nsaves:       ', nsaves
  print *, 'regOpt:       ', regOpt
  print *, 'itseries:     ', itseries
  print *, 'init_file:    ', init_file
  print *, 'ibegin:       ', ibegin
  print *, ' '
  print *, 'Calculated:   '
  print *, 'dr:           ', dr
  print *, 'dz:           ', dz
  if (Ro.eq.0d0.AND.wf.gt.0d0) then
    print *, 'Response Period:       ', T
    print *, 'dt:           ', dt
  elseif (Ro.gt.0d0.AND.wf.gt.0d0) then
    print *, 'Forcing Period:       ', T
    print *, 'dt:           ', dt
  else
    print *, 'Period:            No period'
    print *, 'dt:           ', dt
  endif
  print *, 'Nsteps:       ', Nsteps
  print *, 'igraph:       ', igraph
  print *, 'TU simulated: ', simTU
  print *, '======================================================='
  print *, '======================================================='
  print *, ''

  if (Ro.eq.0d0.AND.wf.gt.0d0) then ! Unforced system, given response frequency
    wf = 0
  endif
  print*, 'Initializing Variables'
  call initialize()

  tfile(1:3)='ts_'
  tfile(4:3+ix)=prefix(1:ix)
  inquire(file=tfile(1:3+ix),exist=existe)
  if (.not.existe) then
    open(unit=15,file=tfile(1:3+ix),status='new',form='formatted')
  else
    open(unit=15,file=tfile(1:3+ix),status='old', access='append')
  end if
!  Save surface velocity
!  ! -----  ERASE THIS SECTION ---- !
!  vfile(1:3)='vs_'
!  vfile(4:3+ix)=prefix(1:ix)
!  inquire(file=vfile(1:3+ix),exist=existe)
!  if (.not.existe) then
!    open(unit=25,file=vfile(1:3+ix),status='new',form='formatted')
!  else
!    open(unit=25,file=vfile(1:3+ix),status='old', access='append')
!  end if
!  ! -----  ERASE THIS SECTION ---- !

  oldtime=0d0
  if(ibegin.eq.0) then
    print*,'Starting from static initial conditions'
    wt=0d0
    Lt=0d0
    sf=0d0
  else
    print*,'Reading from restart file'
    open(unit=1,file=restart(1:irestart),status='old',&
        form='unformatted')
    read(1) Nz1,Nr1,dz1,dr1,dt1,oldtime
    read(1) Re1,Pe1,Ca1,Ro1,wf1,Gama1
    read(1) ((sf(j,i),j=1,Nz),i=1,Nr),&
            ((wt(j,i),j=1,Nz),i=1,Nr),&
            ((Lt(j,i),j=1,Nz),i=1,Nr)
    close(1)
  end if
  if(ibegin.eq.2) oldtime=0d0

! Format for writing time-series
  100 format (ES23.15e3,2x,ES23.15e3,2x,ES23.15e3,2x,ES23.15e3,2x,ES23.15e3,2x,ES23.15e3,2x,ES23.15e3,2x,ES23.15e3)
! Format for writing surface velocity
!  109 format (I4,2x,ES23.15e3,2x,ES23.15e3)

! Define first and second derivative operator (together), B, for r direction
! of the stream function PDE.
  do i=3,Nr-2
    B(i,i-1)=1/dr**2.0d0+1/(r(i)*2.0d0*dr)  !Lower diagonal
    B(i,i)=-2/dr**2.0d0                     !Main diagonal
    B(i,i+1)=1/dr**2.0d0-1/(r(i)*2.0d0*dr)  !Upper diagonal
  end do
  !Extremes of the tridiagonal matrix
  B(2,2)=-2/dr**2.0d0
  B(2,3)=1/dr**2.0d0-1/(r(2)*2.0d0*dr)
  B(Nr-1,Nr-1)=-2/dr**2.0d0
  B(Nr-1,Nr-2)=1/dr**2.0d0+1/(r(Nr-1)*2.0d0*dr)

  print*, 'Calculating eigenvalues'
  call dgeev('n','v',Nr-2,B,Nr-2,eigRe,eigIm,dumeig,Nr-2,P,Nr-2,dumeig,&
                                              (Nr-2)*(Nr-2),info)
  if(info.ne.0) then
    print*,'dgeev info=',info
  end if

  print*, 'Calculating Inverse'
  Paux=P
  call dgesv(Nr-2,Nr-2,Paux,Nr-2,ipiv,Pinv,Nr-2,info)
  if(info.ne.0) then
    print*,'dgesv info=',info
  end if

! Set up entries for the Nr-2 tridiagonal matrices such that the streamfunction
! PDE is expressed as A*Psi=rhs
  L=-1d0
  do j=2,Nr-1
    D(2:Nz-1,j)=2.d0-eigRe(j)*dz**2.0d0
    call dpttrf(Nz-2,D(:,j),L(:,j),info)
    if(info.ne.0) then
      print*,'dpttrf info=',info,j
    end if
  end do

  print*, 'Starting time-stepping procedure'
!-----time-stepping procedure-------------------------------------
  do m=1,Nsteps
    time=m*dt+oldtime

    !First RK2 step
    call rhs_wtLt(wt,Lt,sf,Re,r,dr,dz,wt_rhs,Lt_rhs,Nz,Nr,DsfDz,DsfDr,DLtDz,DLtDr)
    wt_tmp(2:Nz-1,2:Nr-1) = wt(2:Nz-1,2:Nr-1) + dt*wt_rhs(2:Nz-1,2:Nr-1)
    Lt_tmp(2:Nz-1,2:Nr-1) = Lt(2:Nz-1,2:Nr-1) + dt*Lt_rhs(2:Nz-1,2:Nr-1)
    call solve_streamfn(wt_tmp,sf,r,dz,L,D,Nz,Nr,P,Pinv)
    call solve_concentration(c,sf,Pe,r,dz,dr,dt,Nz,Nr,c_tmp,c_rhs)
    call BC_freeSurfTop(wt_tmp,Lt_tmp,sf,c,Ca,wf,Ro,time,r,dr,dz,Nz,Nr)

    !Second RK2 step
    call rhs_wtLt(wt_tmp,Lt_tmp,sf,Re,r,dr,dz,wt_rhs,Lt_rhs,Nz,Nr,DsfDz,DsfDr,DLtDz,DLtDr)
    wt(2:Nz-1,2:Nr-1) = 0.5d0*(wt(2:Nz-1,2:Nr-1) + wt_tmp(2:Nz-1,2:Nr-1) + dt*wt_rhs(2:Nz-1,2:Nr-1))
    Lt(2:Nz-1,2:Nr-1) = 0.5d0*(Lt(2:Nz-1,2:Nr-1) + Lt_tmp(2:Nz-1,2:Nr-1) + dt*Lt_rhs(2:Nz-1,2:Nr-1))
    call solve_streamfn(wt,sf,r,dz,L,D,Nz,Nr,P,Pinv)
    call solve_concentration(c,sf,Pe,r,dz,dr,dt,Nz,Nr,c_tmp,c_rhs)
    call BC_freeSurfTop(wt,Lt,sf,c,Ca,wf,Ro,time,r,dr,dz,Nz,Nr)
! ::::::::::::::::::::::::::::::

    call observables(Ek,Eg,Ex,ulr,ulv,ulz,ekk,egg,exx,sf,Lt,wt,r,DsfDr,DsfDz,DLtDr,DLtDz,Nz,Nr,dz,dr)
    ! Outputs
    if (mod(m,igraph).eq.0) then
      call graphs_freeSurfTop(wt,Lt,sf,Re,Pe,Ca,Ro,wf,Gama,Nz,Nr,&
                          dz, dr, dt, time, prefix, ix, init_file)
    end if
    if (mod(m,itseries).eq.0) then
      write(15,100) time, Ek, Eg, Ex, ulr, ulv, ulz
    end if
  end do
  close(15)
  print*, 'Time-stepping procedure FINISHED'

  contains

    subroutine initialize()
      allocate(wt(Nz,Nr))
      allocate(Lt(Nz,Nr))    !Nr and Nz must be number of NODES
      allocate(sf(Nz,Nr))
      allocate( c(Nr))
    ! Define the radial and vertical vectors r and z
      allocate(r(Nr))
        do i=1,Nr
          r(i)=(i-1)*dr
        end do
      allocate(B(2:Nr-1,2:Nr-1))
        B=0d0
      allocate(eigRe(2:Nr-1))
      allocate(eigIm(2:Nr-1))
      allocate(dumeig(2:Nr-1,2:Nr-1))
      allocate(P(2:Nr-1,2:Nr-1))
      allocate(Paux(2:Nr-1,2:Nr-1))
      allocate(Pinv(2:Nr-1,2:Nr-1))
      allocate(ipiv(2:Nr-1))
        Pinv=0d0
        do i=2,Nr-1
          Pinv(i,i)=1d0
        end do
      allocate(D(2:Nz-1,2:Nr-1))
      allocate(L(3:Nz-1,2:Nr-1))

      allocate(ldiag(3:Nr-1))
      allocate(mdiag(2:Nr-1))
      allocate(udiag(2:Nr-2))

      allocate(wt_tmp(Nz,Nr))
      allocate(Lt_tmp(Nz,Nr))
      allocate(wt_rhs(2:Nz-1,2:Nr-1))
      allocate(Lt_rhs(2:Nz-1,2:Nr-1))
      allocate(DsfDz(2:Nz-1,2:Nr-1))
      allocate(DsfDr(2:Nz-1,2:Nr-1))
      allocate(DLtDz(2:Nz-1,2:Nr-1))
      allocate(DLtDr(2:Nz-1,2:Nr-1))
      allocate(ekk(Nz,Nr))
      allocate(egg(Nz,Nr))
      allocate(exx(Nz,Nr))
      ekk = 0d0
      egg = 0d0
      exx = 0d0
    end subroutine initialize

end program main_freeSurfaceTop
