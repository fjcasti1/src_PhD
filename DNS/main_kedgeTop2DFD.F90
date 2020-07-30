program main_kedgeTop2DFD
  use tools_FD_cyl
  implicit none
  logical  :: existe, regOpt
  integer  :: Nz, Nz1, Nr, Nr1, ned, ned1, Nsteps, nsaves
  integer  :: i, j, m, ix, irestart, ir
  integer  :: igraph, itseries, ibegin, init_file, info
  real*8   :: Re, Re1, Bo, Bo1, beta, beta1, alpha, alpha1, f, f1, wf, wf1, simTU, NT, NtsT, T
  real*8   :: Gama, eta, Hasp, Hasp1, Rasp, Rasp1, dr, dr1, dz, dz1
  real*8   :: time, oldtime, dt, dt1, xmax
  real*8   :: Ek, Eg, Ex, ulr, ulv, ulz
  integer, dimension(:),   allocatable :: ipiv
  real*8,  dimension(:),   allocatable :: ldiag, mdiag, udiag
  real*8,  dimension(:),   allocatable :: eigRe, eigIm
  real*8,  dimension(:),   allocatable :: r, vs
  real*8,  dimension(:,:), allocatable :: wt, Lt, sf, wt_tmp, Lt_tmp
  real*8,  dimension(:,:), allocatable :: wt_rhs, Lt_rhs, DsfDz, DsfDr, DLtDz, DLtDr, ekk, egg, exx
  real*8,  dimension(:,:), allocatable :: L, D, B, P, Paux, Pinv, dumeig
  character*128 :: prefix, restart, fileout, tfile
  character*128 :: vfile ! Erase  this line
  character*10 :: TSformat
  real*8, parameter :: PI = 4.d0*atan(1.0d0)

! Read parameters
  read*, prefix    ! prefix for filenames
  ix=index(prefix//' ',' ')-1
  read*, restart   ! name of restart file
  irestart=index(restart//' ',' ')-1
  read*, Bo
  read*, Re
  read*, alpha
  read*, wf
  read*, Gama
  read*, eta
  read*, Nr
  read*, Nz
  read*, ned
  read*, dt         ! Default dt if wf = 0
  read*, NtsT       ! Number of Time Steps per Period
  read*, NT         ! Number of Periods or Number of Time Units (if wf=0)
  read*, nsaves
  read*, regOpt
  read*, itseries
  read*, init_file
  read*, ibegin     ! controls the start/restart process


  Hasp = Gama/eta
  Rasp =  1d0/eta

  beta = 1.d0
  dr= Rasp/(Nr-1)
  dz= Hasp/(Nz-1)

  if (alpha.gt.0d0.AND.wf.gt.0d0) then  ! Forced system
    T  = 2.d0*PI/wf
    dt = T/NtsT
  elseif (alpha.eq.0d0.AND.wf.gt.0d0) then ! Unforced system, given response frequency
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
  print *, 'Bo:           ', Bo
  if (alpha.eq.0d0.AND.wf.gt.0d0) then ! Unforced system, given response frequency
    print *, 'wf:           ', 0
  endif
  print *, 'wf:           ', wf
  print *, 'beta:         ', beta
  print *, 'alpha:        ', alpha
  print *, 'Hasp:         ', Hasp
  print *, 'Rasp:         ', Rasp
  print *, 'Nr:           ', Nr
  print *, 'Nz:           ', Nz
  print *, 'ned:          ', ned
!  if (alpha.eq.0) then
!    print *, 'dt:           ', dt
!  endif
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
  if (alpha.eq.0d0.AND.wf.gt.0d0) then
    print *, 'Response Period:       ', T
    print *, 'dt:           ', dt
  elseif (alpha.gt.0d0.AND.wf.gt.0d0) then
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

  if (alpha.eq.0d0.AND.wf.gt.0d0) then ! Unforced system, given response frequency
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
  ! -----  ERASE THIS SECTION ---- !
  vfile(1:3)='vs_'
  vfile(4:3+ix)=prefix(1:ix)
  inquire(file=vfile(1:3+ix),exist=existe)
  if (.not.existe) then
    open(unit=25,file=vfile(1:3+ix),status='new',form='formatted')
  else
    open(unit=25,file=vfile(1:3+ix),status='old', access='append')
  end if
  ! -----  ERASE THIS SECTION ---- !

  oldtime=0.d0
  if(ibegin.eq.0) then
    print*,'Starting from static initial conditions'
    wt=0d0
    Lt=0d0
    sf=0d0
  else
    print*,'Reading from restart file'
    open(unit=1,file=restart(1:irestart),status='old',&
        form='unformatted')
    read(1) Nz1,Nr1,ned1,dz1,dr1,dt1,oldtime
    read(1) Re1,Bo1,beta1,alpha1,f1,wf1,Hasp1,Rasp1
    read(1) ((sf(j,i),j=1,Nz),i=1,Nr),&
            ((wt(j,i),j=1,Nz),i=1,Nr),&
            ((Lt(j,i),j=1,Nz),i=1,Nr)
    close(1)
  end if
  if(ibegin.eq.2) oldtime=0.d0

! Format for writing time-series
  100 format (ES23.15e3,2x,ES23.15e3,2x,ES23.15e3,2x,ES23.15e3,2x,ES23.15e3,2x,ES23.15e3,2x,ES23.15e3,2x,ES23.15e3)
  109 format (I4,2x,ES23.15e3,2x,ES23.15e3) ! ERASE THIS LINE

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

  if (Bo == 0d0) then
    call infBoussinesqBC(vs,ir,r,Nr,Rasp,regOpt)
    do i=1,Nr
      write(25,109) i, r(i), vs(i)
    enddo
  elseif (Bo > 0d0) then
!   Create matrices to solve the contaminated free surface
    ir=1+0.5d0*(Nr-1)
    !Lower diagonal
    ldiag(3:ir-ned)=Bo*(1/dr**2d0+1/(2d0*r(3:ir-ned)*dr))
    ldiag(ir-ned+1:ir)=0d0
    ldiag(ir+1:Nr-1)=Bo*(1/dr**2d0+1/(2d0*r(ir+1:Nr-1)*dr))
    !Main diagonal
    mdiag(2:ir-ned)=-2d0*Bo/dr**2d0-0.5d0*3/dz
    mdiag(ir-ned+1:ir)=1d0
    mdiag(ir+1:Nr-1)=-2d0*Bo/dr**2d0-0.5d0*3/dz
    !Upper diagonal
    udiag(2:ir-ned)=Bo*(1/dr**2d0-1/(2d0*r(2:ir-ned)*dr))
    udiag(ir-ned+1:ir)=0d0
    udiag(ir+1:Nr-2)=Bo*(1/dr**2d0-1/(2d0*r(ir+1:Nr-2)*dr))
  endif

  print*, 'Starting time-stepping procedure'
!-----time-stepping procedure-------------------------------------
  do m=1,Nsteps
    time=m*dt+oldtime
      !First RK2 step
      !call rhsXG(x,g,s,Re,r,dr,dz,xrhs,grhs,Nz,Nr,DsDz,DsDr)
      call rhsXG2(wt,Lt,sf,Re,r,dr,dz,wt_rhs,Lt_rhs,Nz,Nr,DsfDz,DsfDr,DLtDz,DLtDr)
      wt_tmp(2:Nz-1,2:Nr-1) = wt(2:Nz-1,2:Nr-1) + dt*wt_rhs(2:Nz-1,2:Nr-1)
      Lt_tmp(2:Nz-1,2:Nr-1) = Lt(2:Nz-1,2:Nr-1) + dt*Lt_rhs(2:Nz-1,2:Nr-1)
        !--call solve_streamfn(xtmp,s,r,dr,dz,L,D,Nz,Nr)
      call solve_streamfn(wt_tmp,sf,r,dr,dz,L,D,Nz,Nr,P,Pinv)
      call BndConds(wt_tmp,Lt_tmp,sf,Bo,wf,beta,alpha,time,r,dr,dz,&
                                        Nz,Nr,ned,ldiag,mdiag,udiag,ir,vs)
     !Second RK2 step
      !call rhsXG(xtmp,gtmp,s,Re,r,dr,dz,xrhs,grhs,Nz,Nr,DsDz,DsDr)
      call rhsXG2(wt_tmp,Lt_tmp,sf,Re,r,dr,dz,wt_rhs,Lt_rhs,Nz,Nr,DsfDz,DsfDr,DLtDz,DLtDr)
      wt(2:Nz-1,2:Nr-1) = 0.5d0*(wt(2:Nz-1,2:Nr-1) + wt_tmp(2:Nz-1,2:Nr-1) + dt*wt_rhs(2:Nz-1,2:Nr-1))
      Lt(2:Nz-1,2:Nr-1) = 0.5d0*(Lt(2:Nz-1,2:Nr-1) + Lt_tmp(2:Nz-1,2:Nr-1) + dt*Lt_rhs(2:Nz-1,2:Nr-1))
        !--call solve_streamfn(x,s,r,dr,dz,L,D,Nz,Nr)
      call solve_streamfn(wt,sf,r,dr,dz,L,D,Nz,Nr,P,Pinv)
      call BndConds(wt,Lt,sf,Bo,wf,beta,alpha,time,r,dr,dz,&
                                        Nz,Nr,ned,ldiag,mdiag,udiag,ir,vs)
!!      call kineticEnergy(Ek,ekk,g,r,DsDr,DsDz,Nz,Nr,dz,dr)
      call observables(Ek,Eg,Ex,ulr,ulv,ulz,ekk,egg,exx,sf,Lt,wt,r,DsfDr,DsfDz,DLtDr,DLtDz,Nz,Nr,dz,dr)
    ! Outputs
    if (mod(m,igraph).eq.0) then
      call graphs(wt,Lt,sf,Re,Bo,beta,alpha,f,wf,Hasp,Rasp,Nz,Nr,ned,dz,dr,&
                                                   dt,time,prefix,ix,init_file)
    end if
    if (mod(m,itseries).eq.0) then
      !xmax=maxval(maxval(x,1))
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
      allocate(vs(Nr))
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

end program main_kedgeTop2DFD
