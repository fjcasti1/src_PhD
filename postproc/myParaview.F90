! (nr ,nz ,nn ) dimensions of the read-in solution
! (nrp,nzp,ntp) dimensions of the uniform grid in cylindrical coordinates.
! nr is radius (0,1], nx is the duplication of nr
program myParaview2
  use iso_c_binding, only: c_double, c_int
  implicit none
  integer , parameter :: dp=kind(0d0)
  integer , parameter :: Nf = 12
  integer , parameter :: VTK_UNIT=1000, LOG_UNIT=2000
  real(dp), parameter :: pi=dacos(-1d0)
  integer  :: nr,nx,nz,nn,nrp,nzp,ntp
  integer  :: stage(2)    ! pointers to solutions at N and N-1
  integer  :: i, j, k, nrzt, jff(Nf), i0mode, ix
  real(dp) :: Bo,Re,Ro,wf,gama,reg,dt,tps
  real(dp) :: h1, thc, ths
  real(dp),allocatable,dimension(:,:,:,:) ::  vr,  vt,  vz
  real(dp),allocatable,dimension(:,:,:)   :: pvr, pvt, pvz
  real(dp),allocatable,dimension(:,:,:)   :: pwr, pwt, pwz, pa
  real(dp),allocatable,dimension(:,:)     :: ps
  real(dp),allocatable,dimension(:)       :: r, pr, pz, pth
  real(dp),allocatable,dimension(:,:)     :: dz1
  real(dp),allocatable,dimension(:,:,:)   :: dr1
  logical   :: existe
  character :: prefix*128

  integer, parameter ::  TEST_UNIT = 69

  read(*,*,err=8,end=8) nrp
  read(*,*,err=8,end=8) nzp
  read(*,*,err=8,end=8) ntp
  read(*,*,err=8,end=8) i0mode
  read(*,*,err=8,end=8) (jff(i),i=1,Nf)
  do i=1,5
    read(*,*,err=8,end=8) ! Advance the 5 comments lines
  enddo

  print*,'Generating vtk file.'

!*****************************************************************
!***** START MAIN LOOP *******************************************
!     opens log file storing plotting info
  open(LOG_UNIT,file='ParaviewPlot.log',form='formatted')

!***** SWEEP THROUGH LIST OF FILES AT BOTTOM OF input file ******
!*****   AND PRODUCE SPECIFIED PLOTS ****************************
7 continue

! Read in next file to be plotted or end program if no more files
  read(*,'(a)',err=9,end=9) prefix ! Read from stdin
!                                              another filename
  ix=index(prefix//' ',' ')-1
  if (ix.le.1) goto 9 ! If the line is blank

  inquire(file=prefix(1:ix),exist=existe)
  if (.not.existe) then
     print *, ' '
     print *, 'File ',prefix(1:ix),' does not exist.'
     goto 7
  endif

  print *, ' '
  print*,'filein=',prefix(1:ix)
  write(LOG_UNIT,'(A,A)') 'filein = ',adjustl(prefix(1:ix))

  call reader(prefix)
  nx=2*nr+1
  call deriv(nr,nx,nz,r,dr1,dz1)

! Read in spectral coefficients and translate to uniform physical grids
  call velvort(nr,nx,nz,nn,vr,vt,vz,r,dr1,dz1, &
    nrp,nzp,ntp,pvr,pvt,pvz,ps,pwr,pwt,pwz,i0mode,gama,h1)

! Print description of plots to be made
  print*,'H = ',h1

  write(LOG_UNIT,'(A,f5.3)') 'H = ',h1

! Set up mesh
  do k=0,nrp
    pr(k)=dfloat(k)/dfloat(nrp)
  enddo
  do i=0,nzp
    pz(i)=h1*(dfloat(i)/dfloat(nzp)-0.5d0)
  enddo
  do j=0,ntp
    pth(j)=2d0*pi*dfloat(j)/dfloat(ntp+1)
  enddo

! GENERATE VTK LEGACY FILE
  nrzt=(nrp+1)*(nzp+1)*(ntp+2)
  if (ntp==0) then
    nrzt=(nrp+1)*(nzp+1)
  endif
! Write out VTK standards
  open(VTK_UNIT,file=prefix(1:ix)//'.vtk')
  write(VTK_UNIT,'(A)')'# vtk DataFile Version 3.0'
  write(VTK_UNIT,'(A)')             'Cylindrical geometry'
  write(VTK_UNIT,'(A)')             'ASCII'
  write(VTK_UNIT,'(A)')             'DATASET STRUCTURED_GRID'
  if (ntp>0) then
    write(VTK_UNIT,'(A,3(1x,I3))')    'DIMENSIONS',nrp+1,ntp+2,nzp+1
  else
    write(VTK_UNIT,'(A,3(1x,I3))')    'DIMENSIONS',nrp+1,1,nzp+1
  endif
  write(VTK_UNIT,'(A,1X,I9,1X,A)')  'POINTS',nrzt,'float'
! Write out grid (points in cartesian coordinates)
  print*,'Grid with',nzp+1,' horizontal layers. Writing  grid:'
  do i=0,nzp
    print '(A,$)','.'
    do j=0,ntp
      thc = cos(real(pth(j)))
      ths = sin(real(pth(j)))
      do k=0,nrp
        write(VTK_UNIT,'(3ES25.15e3)') real(pr(k))*thc, real(pr(k))*ths, real(pz(i))
      enddo
    enddo
!  Close the grid
    if (ntp>0) then
      thc = cos(real(pth(0)))
      ths = sin(real(pth(0)))
      do k=0,nrp
        write(VTK_UNIT,'(3ES25.15e3)') real(pr(k))*thc,real(pr(k))*ths,real(pz(i))
      enddo
    endif
  enddo
  print *, ''
! Prepare for data
  write(VTK_UNIT,*)
  write(VTK_UNIT,'(A,1X,I9)')       'POINT_DATA',nrzt

! Writing grid data
  if ( jff(1) > 0 ) call writevtk(pvr,nrp,nzp,ntp,'vr')
  if ( jff(2) > 0 ) call writevtk(pvt,nrp,nzp,ntp,'vt')
  if ( jff(3) > 0 ) call writevtk(pvz,nrp,nzp,ntp,'vz')
  if ( jff(4) > 0 ) call writevtk(pwr,nrp,nzp,ntp,'wr')
  if ( jff(5) > 0 ) call writevtk(pwt,nrp,nzp,ntp,'wt')
  if ( jff(6) > 0 ) call writevtk(pwz,nrp,nzp,ntp,'wz')

  if ( jff(7) > 0 ) then   ! angular momentum (axial component)
    do k=0,ntp
      do j=0,nzp
        do i=0,nrp
          pa(i,j,k)=pr(i)*pvt(i,j,k)
        enddo
      enddo
    enddo
    call writevtk(pa,nrp,nzp,ntp,'rv')
  endif

  if ( jff(8) > 0 ) then    ! kinetic energy
    do k=0,ntp
      do j=0,nzp
        do i=0,nrp
          pa(i,j,k) = pvr(i,j,k)**2d0 + pvt(i,j,k)**2d0 + pvz(i,j,k)**2d0
        enddo
      enddo
    enddo
    call writevtk(pa,nrp,nzp,ntp,'ke')
  endif

  if ( jff(9) > 0 ) then    ! helicity (velocity*vorticity vectors)
    do k=0,ntp
      do j=0,nzp
        do i=0,nrp
          pa(i,j,k) = pvr(i,j,k)*pwr(i,j,k) + pvt(i,j,k)*pwt(i,j,k) + &
                      pvz(i,j,k)*pwz(i,j,k)
        enddo
      enddo
    enddo
    call writevtk(pa,nrp,nzp,ntp,'he')
  endif

  if ( jff(10) > 0 ) then
    do k=0,ntp
      do j=0,nzp
        do i=0,nrp
          pa(i,j,k) = ps(i,j)
        enddo
      enddo
    enddo
    call writevtk(pa,nrp,nzp,ntp,'sf')
  endif

  if ( jff(11) > 0 ) then    ! x-component of velocity
    do k=0,ntp
      thc = cos(pth(k))
      ths = sin(pth(k))
      do j=0,nzp
        do i=0,nrp
          pa(i,j,k) = pvr(i,j,k)*thc - pvt(i,j,k)*ths
        enddo
      enddo
    enddo
    call writevtk(pa,nrp,nzp,ntp,'vx')
  endif

  if ( jff(12) > 0 ) then    ! y-component of velocity
    do k=0,ntp
      thc = cos(pth(k))
      ths = sin(pth(k))
      do j=0,nzp
        do i=0,nrp
          pa(i,j,k) = pvr(i,j,k)*ths + pvt(i,j,k)*thc
        enddo
      enddo
    enddo
    call writevtk(pa,nrp,nzp,ntp,'vy')
  endif
  write(LOG_UNIT,*) ' '
!***** GOING BACK TO START OF MAIN LOOP TO GET ANOTHER COEF. FILE
  goto 7
!***** OUT OF COEF. FILES - CLOSING input file
8 print *, 'Wrong input file; change it and try again.'

9 continue
  close(LOG_UNIT)
  stop

  contains

    subroutine reader(restart)
      implicit none
      integer, parameter :: IN_UNIT = 1000
      character(*), intent(in) :: restart
      if ( allocated(vr) ) call dealloc()
      open(unit=IN_UNIT,file=restart,status='old',form='unformatted')
      read(IN_UNIT) nr,nz,nn,dt,tps,stage(1),stage(2),Bo,Re,Ro,wf,gama,reg
      call alloc()
      ! The fields are read of size nr to be unfolded to size nx, saves space
      call read_field(nr,nz,nn,2,IN_UNIT,vr)
      call read_field(nr,nz,nn,2,IN_UNIT,vt)
      call read_field(nr,nz,nn,2,IN_UNIT,vz)
      close(IN_UNIT)
      return
    end subroutine reader

    subroutine read_field(nr,nz,nn,nstage,in_unit,field)
      use iso_c_binding, only: c_double, c_int
      implicit none
      integer, parameter :: dp=kind(0.d0)
      integer,  intent(in)    :: nr,nz,nn,nstage
      integer,  intent(in)    :: in_unit
      real(dp), intent(inout) :: field(0:nr,0:nz,0:nn-1,nstage)
      integer :: nrk,nzk,nnk,nstagek
      read(in_unit)                                   &
        (                                             &
          (                                           &
            (                                         &
              (                                       &
                field(nrk,nzk,nnk,nstagek), nrk=0,nr  &
              ), nzk=0,nz                             &
            ), nnk=0,nn-1                             &
          ), nstagek=1,nstage                         &
        )
      return
    end subroutine read_field

    subroutine alloc()
      ! Collocation points in r and 1st derivatives
      allocate(r(0:nr),source=0d0)
      allocate(dr1(0:nr,0:nr,2),source=0d0)
      allocate(dz1(0:nz,0:nz),source=0d0)
      ! Velocities in Collocation Points
      allocate(vr(0:nr,0:nz,0:nn-1,2),source=0d0)
      allocate(vt(0:nr,0:nz,0:nn-1,2),source=0d0)
      allocate(vz(0:nr,0:nz,0:nn-1,2),source=0d0)
      ! Projected Variables (in uniform grid)
      allocate(pr (0:nrp),source=0d0)
      allocate(pz (0:nzp),source=0d0)
      allocate(pth(0:ntp),source=0d0)
      allocate(pvr(0:nrp,0:nzp,0:ntp),source=0d0)
      allocate(pvt(0:nrp,0:nzp,0:ntp),source=0d0)
      allocate(pvz(0:nrp,0:nzp,0:ntp),source=0d0)
      allocate(pwr(0:nrp,0:nzp,0:ntp),source=0d0)
      allocate(pwt(0:nrp,0:nzp,0:ntp),source=0d0)
      allocate(pwz(0:nrp,0:nzp,0:ntp),source=0d0)
      allocate(pa (0:nrp,0:nzp,0:ntp),source=0d0)
      allocate(ps (0:nrp,0:nzp),source=0d0)
    end subroutine alloc

    subroutine dealloc()
      ! Collocation points in r and 1st derivatives
      deallocate(r,dr1,dz1)
      ! Velocities in Collocation Points
      deallocate(vr,vt,vz)
      ! Projected Variables (in uniform grid)
      deallocate(pr,pz,pth,pvr,pvt,pvz,pwr,pwt,pwz,pa,ps)
    end subroutine dealloc

end program myParaview2

!=======================================================================

subroutine writevtk(pa,nrp,nzp,ntp,chf)
  implicit none
  integer, parameter :: dp=kind(0d0)
  integer, parameter :: VTK_UNIT=1000, LOG_UNIT=2000
  integer   :: nrp, nzp, ntp, nrzt, i, j ,k
  real(dp)  :: pa(0:nrp,0:nzp,0:ntp), fmin, fmax
  character :: chf*2

! Find min/max for each computed field and print result
  nrzt=(nrp+1)*(nzp+1)*(ntp+1)
  call dminmax(pa,nrzt,fmin,fmax)
  print *, chf, ' min, max: ', fmin, fmax
  write(LOG_UNIT,*) chf, ' min, max: ', fmin, fmax

! vtk legacy file
  print *, 'Writing field:'
  write(VTK_UNIT,'(A)')'SCALARS '//chf//' float'
  write(VTK_UNIT,'(A)')'LOOKUP_TABLE default'
! Scalar function to print
  do i=0,nzp
    print '(A,$)', '.'
    do j=0,ntp
      do k=0,nrp
        write(VTK_UNIT,'(ES25.15e3)') real(pa(k,i,j))
      enddo
    enddo
    ! Close the grid
    if (ntp>0) then
      do k=0,nrp
        write(VTK_UNIT,'(ES25.15e3)') real(pa(k,i,0))
      enddo
    endif
  enddo
  print*

  return
end subroutine writevtk

!=======================================================================

subroutine velvort(nr,nx,nz,nn,vr,vt,vz,r,dr1,dz1, &
    nrp,nzp,ntp,pvr,pvt,pvz,ps,pwr,pwt,pwz,i0mode,gama,h1)

!     Read in the output from evolcrot.e. Output: pvr, pvt, pvz, pdt, ps,
!     pwr, pwt, pwz (velocity components, temperature-deviation,
!     streamfunction [0-mode], and vorticity components in cylindrical
!     coordinates), in physical space.
!     (vr,vt,vz,tem) the read-in solution, r in (0,1]; (r,z) collocation
!     (vre,vte,vze,teme) extended values, r in [-1,1]; (r,z) collocation
!     (vrs,vts,vzs,tems) spectral coefficients in (r,z,theta)
  implicit none
  integer , parameter :: dp=kind(0d0)
  integer , parameter :: RS_UNIT=3000 ! Restart File Unit
  real(dp), parameter :: pi=dacos(-1d0)
  integer   :: nr, nx, nz, nn
  integer   :: nrp, nzp, ntp
  integer   :: in(2), i, j, k, im, iparity
  integer   :: Nrzp, i0mode
  real(dp)  :: h1, alpha, ALT, angle, coef, px, RAD, gama, parity
  real(dp)  :: r(0:nr), dz1(0:nz,0:nz), dr1(0:nr,0:nr,2)
  real(dp)  ::                          wa (0:nr,0:nz,0:nn-1) ! Auxiliary Vorticity
  real(dp)  :: vr (0:nr,0:nz,0:nn-1,2), wr (0:nr,0:nz,0:nn-1)
  real(dp)  :: vt (0:nr,0:nz,0:nn-1,2), wt (0:nr,0:nz,0:nn-1)
  real(dp)  :: vz (0:nr,0:nz,0:nn-1,2), wz (0:nr,0:nz,0:nn-1)
  real(dp)  :: vre(0:nx,0:nz,0:nn-1)  , wre(0:nx,0:nz,0:nn-1) ! Extended Values
  real(dp)  :: vte(0:nx,0:nz,0:nn-1)  , wte(0:nx,0:nz,0:nn-1) ! Extended Values
  real(dp)  :: vze(0:nx,0:nz,0:nn-1)  , wze(0:nx,0:nz,0:nn-1) ! Extended Values
  real(dp)  :: vrs(0:nx,0:nz,0:nn-1)  , wrs(0:nx,0:nz,0:nn-1) ! Spectral Coefs
  real(dp)  :: vts(0:nx,0:nz,0:nn-1)  , wts(0:nx,0:nz,0:nn-1) ! Spectral Coefs
  real(dp)  :: vzs(0:nx,0:nz,0:nn-1)  , wzs(0:nx,0:nz,0:nn-1) ! Spectral Coefs
  real(dp)  :: pvr(0:nrp,0:nzp,0:ntp), pwr(0:nrp,0:nzp,0:ntp)
  real(dp)  :: pvt(0:nrp,0:nzp,0:ntp), pwt(0:nrp,0:nzp,0:ntp)
  real(dp)  :: pvz(0:nrp,0:nzp,0:ntp), pwz(0:nrp,0:nzp,0:ntp)
  real(dp)  :: cx(0:nx), Chr(0:nrp,0:nx)
  real(dp)  :: cz(0:nz), Chz(0:nzp,0:nz)
  real(dp)  :: Ftr(0:ntp,0:nn-1), Fti(0:ntp,0:nn-1)
  real(dp)  :: ps(0:nrp,0:nzp)
  real(dp)  :: ax1(0:nx,0:nx), az1(0:nz,0:nz), aux1(0:nx,0:nz)
  real(dp)  :: aux2(0:nx,0:nzp), aux3(0:nrp,0:nzp,0:nn-1)

  integer, parameter ::  TEST_UNIT = 69

! notice that aspect ratio = radius / height (Benard convection)
! and so h1, RAD and ALT must be appropriately defined
  h1=gama ; nx=2*nr+1
!  RAD=1d0 ; ALT=gama
  RAD=gama ; ALT=1d0 ! IT WAS LIKE THIS!
  in(1) = 1
  in(2) = 2

  print*, 'h1=',h1,gama

! Computing the vorticity wr, wt, wz > --------------------------------
!   Radial component wr=d_theta(v_z)/r-d_z(v_theta)
!   d_theta(v_z)/r

  !do j=0,nz
  !  do i=0,nr
  !    wr(i,j,0)=0d0 ; wr(i,j,nn-1)=0d0
  !  enddo
  !enddo
  wr(:,:,0) = 0d0
  wr(:,:,nn-1) = 0d0

  !do k=1,nn/2-1
  !  do j=0,nz
  !    do i=0,nr
  !      wr(i,j,2*k  ) = dfloat(k)*vz(i,j,2*k  ,in(1))/(RAD*r(i))
  !      wr(i,j,2*k-1) = dfloat(k)*vz(i,j,2*k-1,in(1))/(RAD*r(i))
  !    enddo
  !  enddo
  !enddo
  do k=1,nn/2-1
    do i=0,nr
      wr(i,:,2*k  ) = dfloat(k)*vz(i,:,2*k  ,in(1))/(RAD*r(i))
      wr(i,:,2*k-1) = dfloat(k)*vz(i,:,2*k-1,in(1))/(RAD*r(i))
    enddo
  enddo
! d_z(v_theta)
  alpha=2d0/ALT
  do k=0,nn-1
    call DGEMM('N','T',nr+1,nz+1,nz+1,alpha,vt(0,0,k,in(1)),nr+1,&
     dz1(0,0),nz+1,0d0,wa(0,0,k),nr+1)
  enddo
! wr
  !do k=0,nn-1
  !  do j=0,nz
  !    do i=0,nr
  !      wr(i,j,k) = wr(i,j,k) - wa(i,j,k)
  !    enddo
  !  enddo
  !enddo
  wr = wr - wa

! Azimuthal component wt=d_z(v_r)-d_r(v_z)
! d_z(v_r)
  alpha=2d0/ALT
  do k=0,nn-1
    call DGEMM('N','T',nr+1,nz+1,nz+1,alpha,vr(0,0,k,in(1)),nr+1,&
        dz1(0,0),nz+1,0d0,wt(0,0,k),nr+1)
  enddo
! d_r(v_z)
  call derivr(vz(0,0,0,in(1)),wa(0,0,0),2,1,nr,nz,nn,&
      nr,nz,nn,RAD,dr1)
! wt
  !do k=0,nn-1
  !  do j=0,nz
  !    do i=0,nr
  !      wt(i,j,k) = wt(i,j,k) - wa(i,j,k)
  !    enddo
  !  enddo
  !enddo
  wt = wt - wa

! Axial component wz=[d_r(r*v_theta)-d_theta(v_r)]/r
! d_r(r*v_theta)
  !do k=0,nn-1
  !  do j=0,nz
  !    do i=0,nr
  !      wa(i,j,k) = RAD*r(i)*vt(i,j,k,in(1))
  !    enddo
  !  enddo
  !enddo
  do i=0,nr
    wa(i,:,:) = RAD*r(i)*vt(i,:,:,in(1))
  enddo
  call derivr(wa(0,0,0),wz(0,0,0),2,1,nr,nz,nn,&
      nr,nz,nn,RAD,dr1)
! d_theta(v_r)
  !do j=0,nz
  !  do i=0,nr
  !    wa(i,j,  0) = 0d0
  !    wa(i,j,nn-1) = 0d0
  !  enddo
  !enddo
  wa(:,:,   0) = 0d0
  wa(:,:,nn-1) = 0d0

  !do k=1,nn/2-1
  !  do j=0,nz
  !    do i=0,nr
  !      wa(i,j,2*k  ) = dfloat(k)*vr(i,j,2*k  ,in(1))
  !      wa(i,j,2*k-1) = dfloat(k)*vr(i,j,2*k-1,in(1))
  !    enddo
  !  enddo
  !enddo
  do k=1,nn/2-1
    wa(:,:,2*k  ) = dfloat(k)*vr(:,:,2*k  ,in(1))
    wa(:,:,2*k-1) = dfloat(k)*vr(:,:,2*k-1,in(1))
  enddo
! wz
  !do k=0,nn-1
  !  do j=0,nz
  !    do i=0,nr
  !      wz(i,j,k) = (wz(i,j,k) - wa(i,j,k)) / (RAD*r(i))
  !    enddo
  !  enddo
  !enddo
  do i=0,nr
    wz(i,:,:) = (wz(i,:,:) - wa(i,:,:)) / (RAD*r(i))
  enddo

! </Computing the vorticity wr, wt, wz > -------------------------------

! <SPECTRAL COEFFICIENTS FOR vr, vt, vz, wr, wt, wz> --------------
! Extending (vr,vt,vz,tem,wr,wt,wz) to [-1,1] in x, preserving parity
  !do k=0,nn-1
  !  do j=0,nz
  !    do i=0,nr
  !      vre(i,j,k) = vr(i,j,k,in(1))
  !      vte(i,j,k) = vt(i,j,k,in(1))
  !      vze(i,j,k) = vz(i,j,k,in(1))
  !      wre(i,j,k) = wr(i,j,k)
  !      wte(i,j,k) = wt(i,j,k)
  !      wze(i,j,k) = wz(i,j,k)
  !    enddo
  !  enddo
  !enddo
  vre(0:nr,:,:) = vr(:,:,:,in(1))
  vte(0:nr,:,:) = vt(:,:,:,in(1))
  vze(0:nr,:,:) = vz(:,:,:,in(1))
  wre(0:nr,:,:) = wr(:,:,:)
  wte(0:nr,:,:) = wt(:,:,:)
  wze(0:nr,:,:) = wz(:,:,:)
  !do j=0,nz
  !  do i=0,nr
  !    vre(nx-i,j,0) =-vre(i,j,0)
  !    vte(nx-i,j,0) =-vte(i,j,0)
  !    vze(nx-i,j,0) = vze(i,j,0)
  !    wre(nx-i,j,0) =-wre(i,j,0)
  !    wte(nx-i,j,0) =-wte(i,j,0)
  !    wze(nx-i,j,0) = wze(i,j,0)
  !  enddo
  !enddo
  do i=0,nr
    vre(nx-i,:,0) =-vre(i,:,0)
    vte(nx-i,:,0) =-vte(i,:,0)
    vze(nx-i,:,0) = vze(i,:,0)
    wre(nx-i,:,0) =-wre(i,:,0)
    wte(nx-i,:,0) =-wte(i,:,0)
    wze(nx-i,:,0) = wze(i,:,0)
  enddo
! TODO: VECTORIZE FROM HERE
  iparity=-1 ; parity=dfloat(iparity)
  do k=1,nn-3,2
    do j=0,nz
      do i=0,nr
        vze (nx-i,j,k)   =  parity*vze (i,j,k)
        vre (nx-i,j,k)   = -parity*vre (i,j,k)
        vte (nx-i,j,k)   = -parity*vte (i,j,k)
        wze (nx-i,j,k)   =  parity*wze (i,j,k)
        wre (nx-i,j,k)   = -parity*wre (i,j,k)
        wte (nx-i,j,k)   = -parity*wte (i,j,k)
        vze (nx-i,j,k+1) =  parity*vze (i,j,k+1)
        vre (nx-i,j,k+1) = -parity*vre (i,j,k+1)
        vte (nx-i,j,k+1) = -parity*vte (i,j,k+1)
        wze (nx-i,j,k+1) =  parity*wze (i,j,k+1)
        wre (nx-i,j,k+1) = -parity*wre (i,j,k+1)
        wte (nx-i,j,k+1) = -parity*wte (i,j,k+1)
      enddo
    enddo
    iparity=-iparity ; parity=dfloat(iparity)
  enddo
  do j=0,nz
    do i=0,nr
      vze (nx-i,j,nn-1) =  parity*vze (i,j,nn-1)
      vre (nx-i,j,nn-1) = -parity*vre (i,j,nn-1)
      vte (nx-i,j,nn-1) = -parity*vte (i,j,nn-1)
      wze (nx-i,j,nn-1) =  parity*wze (i,j,nn-1)
      wre (nx-i,j,nn-1) = -parity*wre (i,j,nn-1)
      wte (nx-i,j,nn-1) = -parity*wte (i,j,nn-1)
    enddo
  enddo
! TODO: VECTORIZE UNTIL HERE
! Auxiliary matrices for the computation of the spectral coefficients
  !cx(0)=2d0 ; cx(nx)=2d0
  !do i=1,nx-1
  !  cx(i) = 1d0
  !enddo
  cx = 1d0 ; cx(0) = 2d0 ; cx(nx) = 2d0
  !What is this? !!do k=1,nx,2
  !What is this? !!!  cx(k) = -cx(k)
  !What is this? !!enddo

  !cz(0)=2d0 ; cz(nz)=2d0
  !do j=1,nz-1
  !  cz(j) = 1d0
  !enddo
  cz = 1d0 ; cz(0) = 2d0 ; cz(nz) = 2d0

  do k=0,nz
    do j=0,nz
      az1(j,k) = 2d0*dcos(pi*dfloat(k*j)/dfloat(nz)) / &
          (dfloat(nz)*cz(j)*cz(k))
    enddo
  enddo
  do i=0,nx
    do j=0,nx
      ax1(i,j) = 2d0*dcos(pi*dfloat(i*j)/dfloat(nx))/ &
          (dfloat(nx)*cx(i)*cx(j))
    enddo
  enddo
! Spectral coefficients for (vr,vt,vz,wr,wt,wz) in (r,z)
  do k=0,nn-1
    call DGEMM('N','T',nx+1,nz+1,nz+1,1d0,vre(0,0,k),nx+1,&
        az1(0,0),nz+1,0d0,aux1(0,0),nx+1)
    call DGEMM('N','N',nx+1,nz+1,nx+1,1d0,ax1(0,0),nx+1,&
        aux1(0,0),nx+1,0d0,vrs(0,0,k),nx+1)              ! vrs
    call DGEMM('N','T',nx+1,nz+1,nz+1,1d0,vte(0,0,k),nx+1,&
        az1(0,0),nz+1,0d0,aux1(0,0),nx+1)
    call DGEMM('N','N',nx+1,nz+1,nx+1,1d0,ax1(0,0),nx+1,&
        aux1(0,0),nx+1,0d0,vts(0,0,k),nx+1)              ! vts
    call DGEMM('N','T',nx+1,nz+1,nz+1,1d0,vze(0,0,k),nx+1,&
        az1(0,0),nz+1,0d0,aux1(0,0),nx+1)
    call DGEMM('N','N',nx+1,nz+1,nx+1,1d0,ax1(0,0),nx+1,   &
        aux1(0,0),nx+1,0d0,vzs(0,0,k),nx+1)              ! vz s
    call DGEMM('N','T',nx+1,nz+1,nz+1,1d0,wre(0,0,k),nx+1, &
        az1(0,0),nz+1,0d0,aux1(0,0),nx+1)
    call DGEMM('N','N',nx+1,nz+1,nx+1,1d0,ax1(0,0),nx+1,   &
        aux1(0,0),nx+1,0d0,wrs(0,0,k),nx+1)              ! wr s
    call DGEMM('N','T',nx+1,nz+1,nz+1,1d0,wte(0,0,k),nx+1, &
        az1(0,0),nz+1,0d0,aux1(0,0),nx+1)
    call DGEMM('N','N',nx+1,nz+1,nx+1,1d0,ax1(0,0),nx+1,   &
        aux1(0,0),nx+1,0d0,wts(0,0,k),nx+1)              ! wt s
    call DGEMM('N','T',nx+1,nz+1,nz+1,1d0,wze(0,0,k),nx+1, &
        az1(0,0),nz+1,0d0,aux1(0,0),nx+1)
    call DGEMM('N','N',nx+1,nz+1,nx+1,1d0,ax1(0,0),nx+1,&
        aux1(0,0),nx+1,0d0,wzs(0,0,k),nx+1)              ! wzs
  enddo
! </SPECTRAL COEFFICIENTS FOR vr, vt, vz, wr, wt, wz> -------------
  im=abs(i0mode)
  if (i0mode.lt.0) then
! Azimuthal mode 0 set to zero
    !do i=0,nx
    !  do j=0,nz
    !    vrs(i,j,0)=0d0 ; vts(i,j,0)=0d0 ; vzs(i,j,0)=0d0
    !    wrs(i,j,0)=0d0 ; wts(i,j,0)=0d0 ; wzs(i,j,0)=0d0
    !  enddo
    !enddo
    vrs(:,:,0)=0d0 ; vts(:,:,0)=0d0 ; vzs(:,:,0)=0d0
    wrs(:,:,0)=0d0 ; wts(:,:,0)=0d0 ; wzs(:,:,0)=0d0
  endif
  if (i0mode.eq.0) then
! Non-azimuthal modes set to zero
    !do i=0,nx
    !  do j=0,nz
    !    do k=1,nn-1
    !      vrs(i,j,k)=0d0 ; vts(i,j,k)=0d0 ; vzs(i,j,k)=0d0
    !      wrs(i,j,k)=0d0 ; wts(i,j,k)=0d0 ; wzs(i,j,k)=0d0
    !    enddo
    !  enddo
    !enddo
    vrs=0d0 ; vts=0d0 ; vzs=0d0
    wrs=0d0 ; wts=0d0 ; wzs=0d0
  endif
! Keeping only the i0mode
  if(im*(1000-im) .gt. 0 .and. im .lt. nn/2) then
    !do i=0,nx
    !  do j=0,nz
    !    do k=1,2*im-2
    !      vrs(i,j,k)=0d0 ; vts(i,j,k)=0d0 ; vzs(i,j,k)=0d0
    !      wrs(i,j,k)=0d0 ; wts(i,j,k)=0d0 ; wzs(i,j,k)=0d0
    !    enddo
    !    do k=2*im+1,nn-1
    !      vrs(i,j,k)=0d0 ; vts(i,j,k)=0d0 ; vzs(i,j,k)=0d0
    !      wrs(i,j,k)=0d0 ; wts(i,j,k)=0d0 ; wzs(i,j,k)=0d0
    !    enddo
    !  enddo
    !enddo
    do k=1,2*im-2
      vrs(:,:,k)=0d0 ; vts(:,:,k)=0d0 ; vzs(:,:,k)=0d0
      wrs(:,:,k)=0d0 ; wts(:,:,k)=0d0 ; wzs(:,:,k)=0d0
    enddo
    do k=2*im+1,nn-1
      vrs(:,:,k)=0d0 ; vts(:,:,k)=0d0 ; vzs(:,:,k)=0d0
      wrs(:,:,k)=0d0 ; wts(:,:,k)=0d0 ; wzs(:,:,k)=0d0
    enddo
  endif

! <UNIFORM GRID in (r,theta,z)> ----------------------------------------
! Auxiliary matrices for the transformation from spectral to physical
  do i=0,nrp
    angle=dacos(dfloat(i)/dfloat(nrp))
    do j=0,nx
      Chr(i,j)=dcos(dfloat(j)*angle)
    enddo
  enddo
  do i=0,nzp
    angle=dacos(2d0*dfloat(i)/dfloat(nzp)-1d0)
    do j=0,nz
      Chz(i,j)=dcos(dfloat(j)*angle)
    enddo
  enddo
  !do j=0,ntp
  !  Ftr(j,0)=1d0 ; Fti(j,0)=-1d0
  !enddo
  Ftr(:,0)=1d0 ; Fti(:,0)=-1d0

  px=2d0*pi/(dfloat(ntp+1))
  do k=1,nn/2-1
    do j=0,ntp
      Ftr(j,2*k  ) =-2d0*dsin(px*dfloat(j*k))
      Ftr(j,2*k-1) = 2d0*dcos(px*dfloat(j*k))

      Fti(j,2*k  ) =-Ftr(j,2*k-1)
      Fti(j,2*k-1) = Ftr(j,2*k)
    enddo
  enddo
  do j=0,ntp
    Ftr(j,nn-1) = dcos(px*dfloat(j*nn/2))
    Fti(j,nn-1) =-Ftr(j,nn-1)
  enddo

  Nrzp=(nrp+1)*(nzp+1)
! pvr, radial velocity, equispaced grid in (r,theta,z)
  do k=0,nn-1
    call DGEMM('N','T',nx+1,nzp+1,nz+1,1d0,vrs(0,0,k),nx+1, &
        Chz(0,0),nzp+1,0d0,aux2(0,0),nx+1)
    call DGEMM('N','N',nrp+1,nzp+1,nx+1,1d0,Chr(0,0),nrp+1,  &
        aux2(0,0),nx+1,0d0,aux3(0,0,k),nrp+1)
  enddo
  call DGEMM('N','T',Nrzp,ntp+1,nn,1d0,aux3(0,0,0),Nrzp, &
      Ftr(0,0),ntp+1,0d0,pvr(0,0,0),Nrzp)
! pvt, azimuthal velocity, equispaced grid in (r,theta,z)
  do k=0,nn-1
     call DGEMM('N','T',nx+1,nzp+1,nz+1,1d0,vts(0,0,k),nx+1,&
          Chz(0,0),nzp+1,0d0,aux2(0,0),nx+1)
     call DGEMM('N','N',nrp+1,nzp+1,nx+1,1d0,Chr(0,0),nrp+1,&
          aux2(0,0),nx+1,0d0,aux3(0,0,k),nrp+1)
  enddo
  call DGEMM('N','T',Nrzp,ntp+1,nn,1d0,aux3(0,0,0),Nrzp, &
       Fti(0,0),ntp+1,0d0,pvt(0,0,0),Nrzp)
! pvz, axial velocity, equispaced grid in (r,theta,z)
  do k=0,nn-1
    call DGEMM('N','T',nx+1,nzp+1,nz+1,1d0,vzs(0,0,k),nx+1,&
        Chz(0,0),nzp+1,0d0,aux2(0,0),nx+1)
    call DGEMM('N','N',nrp+1,nzp+1,nx+1,1d0,Chr(0,0),nrp+1,&
        aux2(0,0),nx+1,0d0,aux3(0,0,k),nrp+1)
  enddo
  call DGEMM('N','T',Nrzp,ntp+1,nn,1d0,aux3(0,0,0),Nrzp,&
      Ftr(0,0),ntp+1,0d0,pvz(0,0,0),Nrzp)
! pwr, radial vorticity, equispaced grid in (r,theta,z)
  do k=0,nn-1
    call DGEMM('N','T',nx+1,nzp+1,nz+1,1d0,wrs(0,0,k),nx+1,&
        Chz(0,0),nzp+1,0d0,aux2(0,0),nx+1)
    call DGEMM('N','N',nrp+1,nzp+1,nx+1,1d0,Chr(0,0),nrp+1,  &
        aux2(0,0),nx+1,0d0,aux3(0,0,k),nrp+1)
  enddo
  call DGEMM('N','T',Nrzp,ntp+1,nn,1d0,aux3(0,0,0),Nrzp,      &
      Fti(0,0),ntp+1,0d0,pwr(0,0,0),Nrzp)
! pwt, azimuthal vorticity, equispaced grid in (r,theta,z)
  do k=0,nn-1
    call DGEMM('N','T',nx+1,nzp+1,nz+1,1d0,wts(0,0,k),nx+1, &
        Chz(0,0),nzp+1,0d0,aux2(0,0),nx+1)
    call DGEMM('N','N',nrp+1,nzp+1,nx+1,1d0,Chr(0,0),nrp+1,  &
        aux2(0,0),nx+1,0d0,aux3(0,0,k),nrp+1)
  enddo
  call DGEMM('N','T',Nrzp,ntp+1,nn,1d0,aux3(0,0,0),Nrzp,&
      Ftr(0,0),ntp+1,0d0,pwt(0,0,0),Nrzp)
! pwz, axial vorticity, equispaced grid in (r,theta,z)
  do k=0,nn-1
    call DGEMM('N','T',nx+1,nzp+1,nz+1,1d0,wzs(0,0,k),nx+1,&
        Chz(0,0),nzp+1,0d0,aux2(0,0),nx+1)
    call DGEMM('N','N',nrp+1,nzp+1,nx+1,1d0,Chr(0,0),nrp+1,&
        aux2(0,0),nx+1,0d0,aux3(0,0,k),nrp+1)
  enddo
  call DGEMM('N','T',Nrzp,ntp+1,nn,1d0,aux3(0,0,0),Nrzp, &
      Fti(0,0),ntp+1,0d0,pwz(0,0,0),Nrzp)

! </UNIFORM GRID in (r,theta,z)> ---------------------------------------

!     ps, streamfunction for the 0 Fourier mode, equispaced grid in (r,z)
!     Integramos de abajo arriba (bad if the boundary layer is at the
!     bottom endwall)
!     do i=1,nrp-1
!        ps(i,0)=0d0
!        coef=0.5d0*dfloat(i)/dfloat(nrp*nzp)
!        do j=1,nzp
!           ps(i,j)=ps(i,j-1)-coef*(pvr(i,j-1,0)+pvr(i,j,0))
!        enddo
!     enddo
!     do j=0,nzp
!        ps(0,j)=0d0
!        ps(nrp,j)=0d0
!     enddo
!     Integramos de arriba abajo
  do i=1,nrp-1
    ps(i,nzp)=0d0
    coef=0.5d0*dfloat(i)/dfloat(nrp*nzp)
    do j=nzp,2,-1
      ps(i,j-1)=ps(i,j)+coef*(pvr(i,j-1,0)+pvr(i,j,0))
    enddo
    ps(i,0)=0d0          !se pone a mano por el BL
  enddo
  do j=0,nzp
    ps(0,j)=0d0
    ps(nrp,j)=0d0
  enddo
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++++++++++++ WRITING DEBUG FILE (will erase) +++++++++++++++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  print *, 'alpha =', alpha
  print *, ' '
  print *, 'wa = ', size(wa,1), 'x', size(wa,2), 'x', size(wa,3)
  print *, 'wa = ', nr+1, 'x', nz+1, 'x', nn
  print *, 'in(1) = ', in(1)
  print *, 'in(2) = ', in(2)
  print *, 'vt(0,0,0,in(1)) = ', vt(0,0,0,in(1))
  print *, ' '
  print *, 'nr =', nr
  print *, 'nx =', nx
  print *, 'nz =', nz
  print *, 'nn =', nn
  print *, ' '
!  open(unit=TEST_UNIT,file='testB',status='new',form='unformatted')
!  write(TEST_UNIT) nr,nx,nz,nn,nrp,nzp,ntp
!  write(TEST_UNIT) (r(i),i=0,nr)
!  write(TEST_UNIT) (((dr1(i,j,k),i=0,nr),j=0,nr),k=1,2)
!  write(TEST_UNIT) ((dz1(i,j),i=0,nz),j=0,nz)
!
!  write(TEST_UNIT) (((vr(i,j,k,1),i=0,nr),j=0,nz),k=0,nn-1)
!  write(TEST_UNIT) (((vt(i,j,k,1),i=0,nr),j=0,nz),k=0,nn-1)
!  write(TEST_UNIT) (((vz(i,j,k,1),i=0,nr),j=0,nz),k=0,nn-1)
!  write(TEST_UNIT) (((wr(i,j,k),i=0,nr),j=0,nz),k=0,nn-1)
!  write(TEST_UNIT) (((wt(i,j,k),i=0,nr),j=0,nz),k=0,nn-1)
!  write(TEST_UNIT) (((wz(i,j,k),i=0,nr),j=0,nz),k=0,nn-1)
!
!  write(TEST_UNIT) (((vre(i,j,k),i=0,nx),j=0,nz),k=0,nn-1)
!  write(TEST_UNIT) (((vte(i,j,k),i=0,nx),j=0,nz),k=0,nn-1)
!  write(TEST_UNIT) (((vze(i,j,k),i=0,nx),j=0,nz),k=0,nn-1)
!  write(TEST_UNIT) (((wre(i,j,k),i=0,nx),j=0,nz),k=0,nn-1)
!  write(TEST_UNIT) (((wte(i,j,k),i=0,nx),j=0,nz),k=0,nn-1)
!  write(TEST_UNIT) (((wze(i,j,k),i=0,nx),j=0,nz),k=0,nn-1)
!
!  write(TEST_UNIT) (((vrs(i,j,k),i=0,nx),j=0,nz),k=0,nn-1)
!  write(TEST_UNIT) (((vts(i,j,k),i=0,nx),j=0,nz),k=0,nn-1)
!  write(TEST_UNIT) (((vzs(i,j,k),i=0,nx),j=0,nz),k=0,nn-1)
!  write(TEST_UNIT) (((wrs(i,j,k),i=0,nx),j=0,nz),k=0,nn-1)
!  write(TEST_UNIT) (((wts(i,j,k),i=0,nx),j=0,nz),k=0,nn-1)
!  write(TEST_UNIT) (((wzs(i,j,k),i=0,nx),j=0,nz),k=0,nn-1)
!
!  write(TEST_UNIT) (((pvr(i,j,k),i=0,nrp),j=0,nzp),k=0,ntp)
!  write(TEST_UNIT) (((pvt(i,j,k),i=0,nrp),j=0,nzp),k=0,ntp)
!  write(TEST_UNIT) (((pvz(i,j,k),i=0,nrp),j=0,nzp),k=0,ntp)
!  write(TEST_UNIT) (((pwr(i,j,k),i=0,nrp),j=0,nzp),k=0,ntp)
!  write(TEST_UNIT) (((pwt(i,j,k),i=0,nrp),j=0,nzp),k=0,ntp)
!  write(TEST_UNIT) (((pwz(i,j,k),i=0,nrp),j=0,nzp),k=0,ntp)
!  close(TEST_UNIT)
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  return
end subroutine velvort

!=======================================================================

subroutine dminmax(f,ndim,fmin,fmax)
  implicit none
  integer, parameter :: dp=kind(0d0)
  integer  :: i, ndim
  real(dp) :: f(ndim), fmin, fmax

  fmin=f(1) ; fmax=f(1)
  do i=2,ndim
    fmin = min(fmin,f(i))
    fmax = max(fmax,f(i))
  enddo
  return
end subroutine dminmax

!=======================================================================

subroutine deriv(nr,nx,nz,r,dr1,dz1)
  implicit none
  integer , parameter :: dp=kind(0d0)
  real(dp), parameter :: pi=dacos(-1d0)
  integer  :: i, j, nr, nx, nz
  real(dp) ::  x (0:nx),       z (0:nz)
  real(dp) :: cx (0:nx),      cz (0:nz)
  real(dp) :: dx1(0:nx,0:nx), dz1(0:nz,0:nz)
  real(dp) :: dx2(0:nx,0:nx), dz2(0:nz,0:nz)
  real(dp) :: dr1(0:nr,0:nr,2), dr2(0:nr,0:nr,2), lap(0:nr,0:nr,2)
  real(dp) :: r(0:nr), sign

!***********************************************************************
!*     X                                                               *
!***********************************************************************

  do i=0,nx
    x(i)=dcos(pi*dfloat(i)/dfloat(nx))
  enddo
  do i=0,nr
    r(i)=dcos(pi*dfloat(i)/dfloat(nx))
  enddo
  cx(0)=2d0 ; cx(nx)=2d0
  do i=1,nx-1
    cx(i)=1d0
  enddo

  sign=1d0
  do j=0,nx
    do i=0,nx
      if(i.ne.j) then
        dx1(i,j)=cx(i)*sign/(cx(j)*(x(i)-x(j)))
      endif
      sign=-sign
    enddo
    if(mod(nx,2).ne.0) sign=-sign
  enddo

  dx1(0,0)=(2d0*nx**2d0+1d0)/6d0
  dx1(nx,nx)=-(2d0*nx**2d0+1d0)/6d0
  do i=1,nx-1
    dx1(i,i)=-x(i)/(2d0*(1d0-x(i)**2d0))
  enddo

!***************************

  sign=-1d0
  do j=0,nx
    do i=1,nx-1
      if(i.ne.j) then
        dx2(i,j)=sign*(x(i)**2d0+x(i)*x(j)-2d0)/&
            (cx(j)*(1d0-x(i)**2d0)*(x(i)-x(j))**2d0)
      endif
      sign=-sign
    enddo
    if(mod(nx,2).ne.0) sign=-sign
  enddo

  dx2(0,0)=(nx**4-1d0)/15d0
  dx2(nx,nx)=(nx**4-1d0)/15d0
  do i=1,nx-1
    dx2(i,i)=-((nx**2d0-1d0)*(1d0-x(i)**2d0)+3d0)/&
        (3d0*(1d0-x(i)**2d0)**2d0)
  enddo

  sign=-1d0
  do j=1,nx
    dx2(0,j)=sign*2d0*((2d0*nx**2d0+1d0)*(1d0-x(j))-6d0)/&
         (3d0*cx(j)*(1d0-x(j))**2d0)
    sign=-sign
  enddo

  sign=1d0
  if(mod(nx,2).ne.0) sign=-1d0
  do j=0,nx-1
    dx2(nx,j)=sign*2d0*((2d0*nx**2d0+1d0)*(1d0+x(j))-6d0)/&
         (3d0*cx(j)*(1d0+x(j))**2d0)
    sign=-sign
  enddo

  do i=0,nr
    do j=0,nr
      dr1(i,j,1) = dx1(i,j) - dx1(i,nx-j)
      dr1(i,j,2) = dx1(i,j) + dx1(i,nx-j)

      dr2(i,j,1) = dx2(i,j) - dx2(i,nx-j)
      dr2(i,j,2) = dx2(i,j) + dx2(i,nx-j)

      lap(i,j,1) = dr2(i,j,1) + dr1(i,j,1)/r(i)
      lap(i,j,2) = dr2(i,j,2) + dr1(i,j,2)/r(i)
    enddo
  enddo

!***********************************************************************
!*     Z                                                               *
!***********************************************************************

  do i=0,nz
    z(i)=dcos(pi*dfloat(i)/dfloat(nz))
  enddo

  cz(0)=2d0
  cz(nz)=2d0
  do i=1,nz-1
    cz(i)=1d0
  enddo

  sign=1d0
  do j=0,nz
    do i=0,nz
      if(i.ne.j) then
        dz1(i,j)=cz(i)*sign/(cz(j)*(z(i)-z(j)))
      endif
      sign=-sign
    enddo
    if(mod(nz,2).ne.0) sign=-sign
  enddo

  dz1(0,0)=(2d0*nz**2d0+1d0)/6d0
  dz1(nz,nz)=-(2d0*nz**2d0+1d0)/6d0
  do i=1,nz-1
    dz1(i,i)=-z(i)/(2d0*(1d0-z(i)**2d0))
  enddo

!***************************

  sign=-1d0
  do j=0,nz
    do i=1,nz-1
      if(i.ne.j) then
        dz2(i,j)=sign*(z(i)**2d0+z(i)*z(j)-2d0)/&
            (cz(j)*(1d0-z(i)**2d0)*(z(i)-z(j))**2d0)
      endif
      sign=-sign
    enddo
    if(mod(nz,2).ne.0) sign=-sign
  enddo

  dz2(0,0)=(nz**4-1d0)/15d0
  dz2(nz,nz)=(nz**4-1d0)/15d0
  do i=1,nz-1
    dz2(i,i)=-((nz**2d0-1d0)*(1d0-z(i)**2d0)+3d0)/&
        (3d0*(1d0-z(i)**2d0)**2d0)
  enddo

  sign=-1d0
  do j=1,nz
    dz2(0,j)=sign*2d0*((2d0*nz**2d0+1d0)*(1d0-z(j))-6d0)/&
        (3d0*cz(j)*(1d0-z(j))**2d0)
    sign=-sign
  enddo

  sign=1d0
  if(mod(nz,2).ne.0) sign=-1d0
  do j=0,nz-1
    dz2(nz,j)=sign*2d0*((2d0*nz**2d0+1d0)*(1d0+z(j))-6d0)/&
        (3d0*cz(j)*(1d0+z(j))**2d0)
    sign=-sign
  enddo

  return
end subroutine deriv

!=======================================================================

subroutine derivr(f,f3,in1,in2,nr2,nz2,nn2,nr,nz,nn,RAD,dr1)
! Se hace la derivada radial de cada modo de Fourier, se indica la paridad
! a partir de los valores de k y l.
! Si k=2 y l=1, empieza por modo cero par y sigue n=1 impar.......
! Si k=1 y l=2, empieza por modo cero impar y sigue n=1 par.......
  implicit none
  integer, parameter :: dp=kind(0d0)
  integer  :: nr2, nz2, nn2, nr, nz, nn
  integer  :: in(2),n,j,k,nk,in1,in2
  real(dp) :: dr1(0:nr,0:nr,2), rad
  real(dp) :: alpha, f(0:nr,0:nz,0:nn-1), f3(0:nr,0:nz,0:nn-1)

  alpha=1d0/RAD ; in(1)=in1 ; in(2)=in2

! Hacemos n=0
  call DGEMM('N','N',nr2+1,nz2+1,nr2+1,alpha,dr1(0,0,in(1)),nr+1,&
            f(0,0,0),nr+1,0d0,f3(0,0,0),nr+1)

  j=in(1) ; in(1)=in(2) ; in(2)=j

! Hacemos los otros n
  do n=1,nn2-3,2
    do k=0,1
      nk=n+k
      call DGEMM('N','N',nr2+1,nz2+1,nr2+1,alpha,dr1(0,0,in(1)),nr+1,&
          f(0,0,nk),nr+1,0d0,f3(0,0,nk),nr+1)
    enddo
    j=in(1) ; in(1)=in(2) ; in(2)=j
  enddo

! Hacemos el nn2/2 -1 (solo hay uno)
  nk=nn2-1
  call DGEMM('N','N',nr2+1,nz2+1,nr2+1,alpha,dr1(0,0,in(1)),nr+1,&
      f(0,0,nk),nr+1,0d0,f3(0,0,nk),nr+1)

  return
end subroutine derivr
