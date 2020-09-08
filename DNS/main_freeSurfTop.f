ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c $Id$
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program main_kedgeTop
      implicit none
      integer nr,nx,nz,nmax,nn,nnh,nnhp,nn1,nr1,nz1,nrpz
c     nr, nz, nn all even; nn .ge. 2; nmax=max(nr,nz,nn)
c     nmax is needed in work arrays in spectral_tools.
      parameter(nr=100,nz=100,nn=32,nmax=100)
      parameter(nx=2*nr+1,nnh=nn/2,nnhp=nnh+1,nn1=nn-1)
      parameter(nr1=nr+1,nz1=nz+1,nrpz=nx+nz-2)
      real*8 dx1(0:nx,0:nx)  , dx2(0:nx,0:nx)
      real*8 dz1(0:nz,0:nz)  , dz2(0:nz,0:nz)
      real*8 dr1(0:nr,0:nr,2), dr2(0:nr,0:nr,2), lap(0:nr,0:nr,2)
      real*8 dummyMM(0:nr,0:nr,2)
      real*8 x(0:nx),z(0:nz),r(0:nr),ri(0:nr),riq(0:nr)
C     x(i)=r(i)=cos(pi*i/nx) but with different dimensions.
C     rphys(i)=rad*r(i)=r(i); ri(i)=1d0/r(i)
C     z(i)=cos(pi*i/nz); zphys(i)=alt*(1+z(i))/2=gama*(1+z(i))/2
      real*8 axd(1:nr,1:nr,0:nnhp),vxd(1:nr,0:nnhp),
     &       bxd(1:nr,1:nr,0:nnhp),axn(1:nr,1:nr,0:nnhp),
     &       vxn(1:nr,0:nnhp),bxn(1:nr,1:nr,0:nnhp),
     &       azd(1:nz-1,1:nz-1),vzd(1:nz-1),bzd(1:nz-1,1:nz-1),
     &       azn(1:nz-1,1:nz-1),vzn(1:nz-1),bzn(1:nz-1,1:nz-1)
      real*8 cx01(2),cx03(1:nr,2),cz01,cz02,cz03(1:nz-1),
     &       cz11,cz12,cz13(1:nz-1)
      real*8 tdr(0:nn1,0:nn1),tdi(0:nn1,0:nn1)
      real*8 tir(0:nn1,0:nn1),tii(0:nn1,0:nn1)
      real*8 vr(0:nr,0:nz,0:nn1,2)
      real*8 vt(0:nr,0:nz,0:nn1,2),vz(0:nr,0:nz,0:nn1,2)
      real*8 wkdiag(0:nr,0:nz,0:nn1,2),nlvr(0:nr,0:nz,0:nn1,2)
      real*8 nlvt(0:nr,0:nz,0:nn1,2),nlvz(0:nr,0:nz,0:nn1,2)
      real*8 rotvorz(0:nr,2,0:nn1,2),rotvorr(0:nz,0:nn1,2)
      real*8 dt,pnu,rad,alt
      real*8 Re,Pe,Ca,Ro,wf,c0,gama
      real*8 pert,time
      real*8 h(0:nr,0:nz,0:nn1),ccr(0:nx),ccz(0:nz),h2(0:nr,0:nn1)
      real*8 Ftr(0:nn1,10)  !up to ten values in phase
      real*8 axx(0:nx,10),azz(0:nz,10)  !up to ten values in phase
      real*8 wk(0:nr,0:nz,0:nn1,6)
      integer nsteps,itseries,insec
      integer stage(2),i,j,k,n,m,nrz,ibegin,ix,irestart,init_file
      integer nrpp(10),ntpp(10),nzpp(10),iaxisym,imode,npp
      integer NtsT
      character(128) restart,prefix,file_out
      integer, parameter :: OUT_RESTART_UNIT = 1000
      integer, parameter :: OUT_TS_UNIT      = 1001
c      integer, parameter :: OUT_TSTXT_UNIT   = 2001
      real*8 , parameter :: pi = dacos(-1d0)


      nrz=(nr+1)*(nz+1)

      call readinout(prefix,ix,restart,irestart,r,
     &     Re,Ro,wf,gama,pnu,nsteps,insec,
     &     init_file,iaxisym,ibegin,imode,pert,itseries,npp,nrpp,ntpp,
     &     nzpp,rad,alt,vr,vt,vz,stage,time,file_out,nr,nz,nn,NtsT,dt)
      print*,' '
      print*,'--- Done readinout ---'
      print*,' '


C <INITIALIZATION> ----------------------------------------------------------
      call grid(ccr,ccz,npp,nrpp,ntpp,nzpp,Ftr,axx,azz,nr,nz,nn)
      call deriv(dx1,dz1,dx2,dz2,dr1,dr2,lap,dummyMM,x,z,r,ri,riq,nr,nz)
      call coef(dr1,dz1,cx01,cx03,cz01,cz02,cz03,cz11,cz12,cz13,nr,nz)
c     wkdiag is a workspace, changed on exit; it is reset after call mdiag
      call mdiag(wkdiag,dz2,lap,riq,axd,vxd,bxd,azd,vzd,bzd,axn,
     &     vxn,bxn,azn,vzn,bzn,cx03,cz03,cz13,nr,nz,nn,nmax)
      call transform(tdr,tir,tdi,tii,nn)

C---- Computing the nonlinear terms and pressure bc for the initial
C     solution at N-1
      call tnlrot(stage(2),tdr,tir,tdi,tii,dr1,dz1,ri,vr,vt,vz,
     &     nlvr,nlvt,nlvz,rotvorz,rotvorr,rad,alt,wk,nr,nz,nn)
C </INITIALIZATION> -------------------------------------------------------

C <MAIN LOOP> ----------------------------------------------------------
      do m=1,nsteps

C     Advancing in time
         time=time+dt
c     stage(1) is the solution at N, stage(2) at N-1
         call time_stepper(Re,Ro,wf,dt,pnu,rad,alt,
     &        stage,time,tdr,tir,tdi,tii,wk,dz1,dz2,nr,nz,nn,
     &        dr1,lap,z,r,ri,
     &        axd,vxd,bxd,azd,vzd,bzd,axn,vxn,bxn,azn,vzn,bzn,
     &        cx01,cx03,cz01,cz02,cz03,cz11,cz12,cz13,
     &        vr,vt,vz,nlvr,nlvt,nlvz,rotvorz,rotvorr)
c     stage(1) is the solution at N, stage(2) at N+1
         stage(1)=stage(2) ; stage(2)=3-stage(1)
c     stage(1) is the solution at N+1, stage(2) at N
         n=stage(1)   ! pointer to the just computed solution at time

c     Enforce symmetry
         call enforce_symmetry(iaxisym,vr,vt,vz,nr,nz,nn)

C <Saving results> -----------------------------------------------------
c     write out time series
        if(mod(m,itseries).eq.0) then
c           open(12,file='phasept.'//prefix(1:ix),status='old',
c     &          access='append')
           call tswrite(vr,vt,vz,time,stage,npp,r,axx,azz,Ftr,ccr,
     &                  ccz,h,h2,nr,nz,nn,nx,rad,pnu,Re)
c           close(12)
        endif

c     saving the solution
        if(mod(m,insec).eq.0.OR.m.eq.nsteps) then
          write(file_out(ix+2:ix+5),'(i4.4)') init_file
          init_file=init_file+1
          open(unit=OUT_RESTART_UNIT,file=file_out(1:ix+5),
     &      form='unformatted')
          write(OUT_RESTART_UNIT) nr,nz,nn,dt,time,
     &      stage(1),stage(2),Re,Pe,Ca,Ro,wf,c0,gama
          write(OUT_RESTART_UNIT)
     &      ((((vr(i,j,k,n),i=0,nr),j=0,nz),k=0,nn-1),n=1,2)
          write(OUT_RESTART_UNIT)
     &      ((((vt(i,j,k,n),i=0,nr),j=0,nz),k=0,nn-1),n=1,2)
          write(OUT_RESTART_UNIT)
     &      ((((vz(i,j,k,n),i=0,nr),j=0,nz),k=0,nn-1),n=1,2)
          close(OUT_RESTART_UNIT)
        endif
C </Saving results> ----------------------------------------------------

      enddo
C </MAIN LOOP> ---------------------------------------------------------
      close(OUT_TS_UNIT)
c      close(OUT_TSTXT_UNIT)
      stop
      end program main_kedgeTop

c=======================================================================

      subroutine bcvel(time,r,vr,vt,vz,rad,alt,
     &          nr,nz,nn,Re,Ro,wf,pnu,dz1)
      implicit none
      integer nr,nz,nn,i,j,k
      real*8 r(0:nr),time
      real*8 vr(0:nr,0:nz,0:nn-1),vt(0:nr,0:nz,0:nn-1)
      real*8 vz(0:nr,0:nz,0:nn-1)
      real*8 dz1(0:nz,0:nz)
      real*8 rad,alt,Re,Ro,wf,pnu
      real*8 freg,eps

c     r(i)=cos(pi*i/nx) in (0,1], radial collocation points

c     Regularization factor freg with coefficient eps
      eps = 5d-3

c---- First, no-slip conditions in all boundaries, in all Fourier modes
      do k=0,nn-1
         do i=0,nr
            vr(i, 0,k)=0d0; vt(i, 0,k)=0d0; vz(i, 0,k)=0d0 ! Top
            vr(i,nz,k)=0d0; vt(i,nz,k)=0d0; vz(i,nz,k)=0d0 ! Bottom
         enddo
         do j=1,nz-1
            vr(0,j,k)=0d0; vt(0,j,k)=0d0; vz(0,j,k)=0d0 ! Sidewall
         enddo
      enddo

c---- Now, the non-zero boundary conditions
      ! Bottom
      do i=0,nr
        freg = 1d0-dexp((r(i)**2d0-1)/eps)
        vt(i,nz,0) = -(pnu*Re/rad)*(1d0+Ro*dsin(wf*time))*r(i)*freg
      enddo

      return
      end subroutine bcvel

      subroutine readinout(prefix,ix,restart,irestart,r,
     &     Re,Ro,wf,gama,pnu,nsteps,insec,
     &     init_file,iaxisym,ibegin,imode,pert,itseries,npp,nrp,ntp,
     &     nzp,rad,alt,vr,vt,vz,stage,time,file_out,nr,nz,nn,NtsT,dt)
      implicit none
      integer :: nr,nz,nn,nx,nnh,nnhp,nn1
      integer, parameter :: IN_RESTART_UNIT  = 1000
      integer, parameter :: OUT_TS_UNIT      = 1001
c      integer, parameter :: OUT_TSTXT_UNIT   = 2001
      real*8 , parameter :: pi = dacos(-1d0)
      real*8  r(0:nr),vr(0:nr,0:nz,0:nn-1,2)
      real*8  vt(0:nr,0:nz,0:nn-1,2), vz(0:nr,0:nz,0:nn-1,2)
      real*8  dt ,time ,Re ,Pe ,Ca ,Ro ,wf ,c0 ,gama ,pnu
      real*8  dt1,time1,Re1,Pe1,Ca1,Ro1,wf1,c01,gama1
      real*8  pert,tinicial,rad,alt,unifRand
      integer ix,irestart,nsteps,insec,init_file,iaxisym,ibegin,imode
      integer nrp(10),ntp(10),nzp(10),itseries,nr2,nz2,nn2,stage(2)
      integer i,j,k,n,isig,jj,npp
      integer NtsT,NT,Nsaves
      character*128 restart,prefix,file_out,time_series
c      character*128 txt_time_series
      logical existe

      nx=2*nr+1 ; nnh=nn/2 ; nnhp=nnh+1 ; nn1=nn-1

C <Read input parameters> ---------------------------------------------------
C     Reading from standard input, e.g.: ./evolcrot.e < in_evolcrot &
      read*, prefix    ! Prefix for filenames
      read*, restart   ! Name of restart file
      read*, Re        ! Omega*R^2/nu
      read*, Pe        ! Peclet
      read*, Ca        ! Capillary
      read*, Ro        ! Omega/Omega
      read*, wf        ! Frequency of the sin(Wot) function
      read*, c0        ! Initial concentration
      read*, gama      ! H/R Aspect Ratio
      read*, dt      ! Number of timesteps per period (used to get dt)
      read*, NtsT      ! Number of timesteps per period (used to get dt)
      read*, NT        ! Number of periods
      read*, Nsaves    ! Write Nsaves solutions, used to calculate insec
      read*, itseries  ! Write in time-series-files every itseries time-steps
      read*, init_file ! Number of first output file
      read*, iaxisym   ! m=0 axisym; m=1 Full 3D Sol; m>0 m-Fourier subspace; m<0 rotoreflec.
      do i=1,4
        read *
      enddo
      read*, ibegin    ! Controls the start/restart process
      do i=1,5
        read *
      enddo
      read*, imode     ! Azimuthal mode to be perturbed
      read*, pert      ! Amplitude of the T perturbation

c      call sleep(2)
c      stop
      ix=index(prefix//' ',' ')-1
      irestart=index(restart//' ',' ')-1

c <Calculate some variables>
      nsteps = NtsT*NT        ! Number of time-steps
      insec  = nsteps/Nsaves  ! Write a full soln every insec time-steps

c     Time-Step size
      if (Ro > 0d0 .and. wf > 0d0) then
        dt = (2d0*pi/wf)/NtsT
      elseif (Ro == 0d0 .and. wf > 0d0) then ! Manual input of dt via wf
        dt = wf
      elseif (Ro == 0d0 .and. wf ==0d0) then ! Default dt
        dt = 5d-3
      else
        print*,'(Ro,wf) = (', Ro, wf,'). They should not be negative.'
        stop
      endif

C <Print general case info> -------------------------------------------------
      print*, ''
      print*, '======================================================='
      print*, '================ INPUT PARAMETERS READ ================'
      print*, '======================================================='
      print *, ' '
      print *, 'Prefix:    ', prefix
      if (ibegin .gt. 0)  then
        print *, 'Restart:   ', restart
      else
        print *, 'Restart:   None'
      endif
      print *, 'ibegin:    ', ibegin
      if (ibegin.eq.0) then
        print*,'Start from rest with random perturbation (if specified),
     & set t=0'
      elseif (ibegin.eq.1) then
         print*,'Continue restart solution, keep t'
      elseif (ibegin.eq.2) then
         print*,'Continue restart solution, set t=0'
      elseif (ibegin.eq.3) then
         print*,'Continue restart solution, with random pertrubation,
     & set t=0.'
      elseif (ibegin.eq.-1) then
         print*,'Start from solid body rotation + pert, set t=0.
     & NOT WORKING.'
      else
        print*,'Wrong ibegin value'
        stop
      endif
      print *, ' '
      print *, '-------------------------'
      print *, '-- Physical Parameters --'
      print *, '-------------------------'
      print *, ' '
      print *, 'Re:        ', Re
      print *, 'Ro:        ', Ro
      print *, 'wf:        ', wf
      print *, ' '
      print *, '-------------------------'
      print *, '---- Geometry & Mesh ----'
      print *, '-------------------------'
      print *, ' '
      print *, 'gama:      ', gama
      print *, 'nr:        ', nr
      print *, 'nz:        ', nz
      print *, 'nn:        ', nn
      print *, ' '
      print *, '-------------------------'
      print *, '----- Time Stepping -----'
      print *, '-------------------------'
      print *, ' '
      print *, 'NtsT:      ', NtsT
      print *, 'NT:        ', NT
      print *, ' '
      print *, 'Calculated'
      print *, '----------'
      print *, 'nsteps:    ', nsteps
      print *, 'dt:        ', dt
      print *, ' '
      print *, '-------------------------'
      print *, '--- Writing Solutions ---'
      print *, '-------------------------'
      print *, ' '
      print *, 'Nsaves:    ', Nsaves
      print *, 'itseries:  ', itseries
      print *, 'init_file: ', init_file
      print *, ' '
      print *, 'Calculated'
      print *, '----------'
      print *, 'insec:     ', insec
      print *, ' '
      print *, '-------------------------'
      print *, '----- Perturbations -----'
      print *, '-------------------------'
      print *, ' '
      print *, 'imode:     ', imode
      print *, 'pert:      ', pert

      if (imode.ge.nnh) then
         print*,'The mode perturbed must be less than ',nnh,
     &        ' and it is',imode
         stop
      endif

      print *, ' '
      print *, '-------------------------'
      print *, '---- Simulation Mode ----'
      print *, '-------------------------'
      print *, ' '

      print *, 'iaxisym:    ', iaxisym
      if (iaxisym.eq.0) then
         print*,'Computing into the AXISYMMETRIC subspace'
      elseif (iaxisym.eq.1) then
         print*,'Computing the FULL 3D solution'
      elseif (iaxisym.gt.1) then
         print*,'Computing into the ',iaxisym,' Fourier subspace'
      else
        print*,'Computing into the ',-iaxisym,' rotoreflection subspace'
      endif
      tinicial=0d0

c     Computing at the physical points (r,theta,z)~(nrp,ntp,nzp) where
c     0<=nrp<=100, 0<=ntp<=360, 0<=nzp<=100.
c     There are npp values of vz to be stored in phase.xxx; npp<=10.
      npp=1
      nrp(1)=98 ; ntp(1)=0 ; nzp(1)=50
      print *, ' '
      print *, '-------------------------'
      print *, '----- Probes Placed -----'
      print *, '-------------------------'
      print *, ' '
      print *, 'There are ', npp, 'probes placed'

c     Nondimensional parameters: Reynolds=Omega*R^2/nu, Gamma=H/R;
c        Rossby=omega/Omega
c     Let ell, tau be the space and time scales (for nondimensionalization).
c     Definitions: alt=H/ell, rad=R/ell, pnu=nu*tau/ell^2 (nondimensional).
c     Velocity scale: Omega*R*(tau/ell)=Reynolds*pnu/rad (nondimensional).
c     Nondimensionalization: ell=R, tau=1/Omega
      rad=1d0 ; alt=gama ; pnu=1d0/Re
C </Read input parameters> --------------------------------------------------

C <START/RESTART PROCESS> ---------------------------------------------------
      if(ibegin.eq.0) then
c---- starting from rest
         do n=1,2
            do k=0,nn-1
               do j=0,nz
                  do i=0,nr
                     vr(i,j,k,n)=0.d0
                     vt(i,j,k,n)=0.d0
                     vz(i,j,k,n)=0.d0
                  enddo
               enddo
            enddo
         enddo
         stage(1)=1 ; stage(2)=2
      elseif (ibegin.eq.1 .or. ibegin.eq.2 .or. ibegin.eq.3) then
c---- starting from a previous solution in restart
        inquire(file=restart(1:irestart),exist=existe)
        if (.not.existe) then
          print*,'File ',restart(1:irestart),' does not exist.'
          stop
        endif
        open(IN_RESTART_UNIT,file=restart(1:irestart),
     &    form='unformatted')
        read(IN_RESTART_UNIT) nr2,nz2,nn2,dt1,time1,
     &    stage(1),stage(2),Re1,Pe1,Ca1,Ro1,wf1,c01,gama1
        if (nr2.ne.nr .or. nz2.ne.nz .or. nn2.ne.nn) then
          print*,'incorrect number of spectral modes in ',
     &         restart(1:irestart),' :'
          print*,'nr=',nr2,' nz=',nz2,' nn=',nn2,
     &         '; fix it using ''cambio'' and try again.'
          stop
        endif
        read(IN_RESTART_UNIT)
     &     ((((vr(i,j,k,n),i=0,nr),j=0,nz),k=0,nn1),n=1,2)
        read(IN_RESTART_UNIT)
     &     ((((vt(i,j,k,n),i=0,nr),j=0,nz),k=0,nn1),n=1,2)
        read(IN_RESTART_UNIT)
     &     ((((vz(i,j,k,n),i=0,nr),j=0,nz),k=0,nn1),n=1,2)
        close(IN_RESTART_UNIT)
      elseif(ibegin.eq.-1) then
         do i=0,nr
            r(i)=dcos(pi*dfloat(i)/dfloat(2*nr+1))
            print*,'r=',r(i)
         end do
c---- starting from solid-body rotation
         do n=1,2
            do k=1,nn-1         ! non-axisymmetric modes
               do j=0,nz
                  do i=0,nr
                     vr(i,j,k,n)=0.d0
                     vt(i,j,k,n)=0.d0
                     vz(i,j,k,n)=0.d0
                  enddo
               enddo
            enddo
            do j=0,nz           ! zero Fourier mode (axisymmetric)
               do i=0,nr
                  vr(i,j,0,n)=0.d0
                  vt(i,j,0,n)=-r(i)*(pnu*Re/rad)
                  vz(i,j,0,n)=0.d0
               enddo
            enddo
         enddo
         stage(1)=1 ; stage(2)=2
      endif

      call random_seed()
      if (ibegin.eq.0) then
C     Start from rest with random perturbation in vr
c       The factor of 10 adjust so the perturbation is of the order
c       introduced, since the random number is between 0 and 1
         time=tinicial
         if (imode.eq.0) then
           do j=0,nz
             do i=0,nr
               call random_number(unifRand)
               vr(i,j,0,1)=(2*unifRand-1)*10*pert
               call random_number(unifRand)
               vr(i,j,0,2)=(2*unifRand-1)*10*pert
             enddo
           enddo
         else
           do j=0,nz
             do i=0,nr
               call random_number(unifRand)
               vr(i,j,2*imode-1,1)=(2*unifRand-1)*10*pert
               call random_number(unifRand)
               vr(i,j,2*imode,  1)=(2*unifRand-1)*10*pert
               call random_number(unifRand)
               vr(i,j,2*imode-1,2)=(2*unifRand-1)*10*pert
               call random_number(unifRand)
               vr(i,j,2*imode,  2)=(2*unifRand-1)*10*pert
             enddo
           enddo
         endif
      elseif (ibegin.eq.1) then
C     Continue restart solution, keep t.
         time=time1
      elseif (ibegin.eq.2) then
C     Continue restart solution, set t=0.
         time=tinicial
      elseif (ibegin.eq.3) then
C     Continue restart solution with random perturbation, set t=0.
         time=tinicial
         if (imode.eq.0) then
            do j=0,nz
               do i=0,nr
                  call random_number(unifRand)
                  vr(i,j,0,1) = vr(i,j,0,1) + (2*unifRand-1)*10*pert
                  call random_number(unifRand)
                  vr(i,j,0,2) = vr(i,j,0,2) + (2*unifRand-1)*10*pert
               enddo
            enddo
         else
            do j=0,nz
              do i=0,nr
                call random_number(unifRand)
                vr(i,j,2*imode-1,1) = vr(i,j,2*imode-1,1)
     &                                + (2*unifRand-1)*10*pert
                call random_number(unifRand)
                vr(i,j,2*imode,  1) = vr(i,j,2*imode,  1)
     &                                + (2*unifRand-1)*10*pert
                call random_number(unifRand)
                vr(i,j,2*imode-1,2) = vr(i,j,2*imode-1,2)
     &                                + (2*unifRand-1)*10*pert
                call random_number(unifRand)
                vr(i,j,2*imode,  2) = vr(i,j,2*imode,  2)
     &                                + (2*unifRand-1)*10*pert
              enddo
            enddo
         endif
      endif

C     Computing into the axisymmetric subspace
      if (iaxisym.eq.0) then
         do n=1,2
            do k=1,nn-1
               do j=0,nz
                  do i=0,nr
                     vr(i,j,k,n)=0.d0
                     vt(i,j,k,n)=0.d0
                     vz(i,j,k,n)=0.d0
                  enddo
               enddo
            enddo
         enddo
      endif

C     Computing into the m-Fourier subspace (modes m, 2m, 3m, 4m ...)
      if (abs(iaxisym).gt.1) then
         do n=1,2
            do jj=1,2*abs(iaxisym)-2
               do k=jj,nn-1,2*abs(iaxisym)
                  do j=0,nz
                     do i=0,nr
                        vr(i,j,k,n)=0.d0
                        vt(i,j,k,n)=0.d0
                        vz(i,j,k,n)=0.d0
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif

C     Computing into the m-rotoreflection subspace
      if (iaxisym.lt.0) then
         do n=1,2
c     0-Fourier mode
            do j=0,nz/2-1
               do i=0,nr
                  vr(i,j,0,n)=0.5d0*(vr(i,j,0,n)+vr(i,nz-j,0,n))
                  vr(i,nz-j,0,n)=vr(i,j,0,n)
                  vt(i,j,0,n)=0.5d0*(vt(i,j,0,n)+vt(i,nz-j,0,n))
                  vt(i,nz-j,0,n)=vt(i,j,0,n)
                  vz(i,j,0,n)=0.5d0*(vz(i,j,0,n)-vz(i,nz-j,0,n))
                  vz(i,nz-j,0,n)=-vz(i,j,0,n)
               enddo
            enddo
            do i=0,nr
               vz(i,nz/2,0,n)=0.d0 !mid point
            enddo
c     multiple of abs(iaxisym) Fourier modes
            isig=-1
            do jj=2*abs(iaxisym),nn-1,2*abs(iaxisym)
               do k=jj-1,jj
                  do j=0,nz/2-1
                     do i=0,nr
                        vr(i,j,k,n)=0.5d0*(vr(i,j,k,n)
     &                               +isig*vr(i,nz-j,k,n))
                        vr(i,nz-j,k,n)=isig*vr(i,j,k,n)
                        vt(i,j,k,n)=0.5d0*(vt(i,j,k,n)
     &                               +isig*vt(i,nz-j,k,n))
                        vt(i,nz-j,k,n)=isig*vt(i,j,k,n)
                        vz(i,j,k,n)=0.5d0*(vz(i,j,k,n)
     &                               -isig*vz(i,nz-j,k,n))
                        vz(i,nz-j,k,n)=-isig*vz(i,j,k,n)
                     enddo
                  enddo
                  do i=0,nr     !mid point
                     if (isig.eq.1) then
                        vz(i,nz/2,k,n)=0.d0
                     else
                        vr(i,nz/2,k,n)=0.d0
                        vt(i,nz/2,k,n)=0.d0
                     endif
                  enddo
               enddo
               isig=-isig
            enddo
c     just in case nn is multiple of 2*abs(iaxisym)
            do j=0,nz
               do i=0,nr
                  vr(i,j,nn-1,n)=0.d0
                  vt(i,j,nn-1,n)=0.d0
                  vz(i,j,nn-1,n)=0.d0
               enddo
            enddo
         enddo
      endif

C </START/RESTART PROCESS> --------------------------------------------------

C <Output files> ------------------------------------------------------------
      file_out(1:ix+1)=prefix(1:ix)//'_'

      time_series(1:ix+3)='ts_'//prefix(1:ix)
c      txt_time_series(1:ix+6)='txtts_'//prefix(1:ix)
C     if time_series do not exists it is created
      inquire(file=time_series(1:ix+3),exist=existe)
      if (.not.existe) then
        call init_tsFile(time_series(1:ix+3),nnh,npp,nsteps/itseries)
      endif
      open(OUT_TS_UNIT,file=time_series(1:ix+3),status='old',
     &  form='unformatted',position='append',access='stream')
c     TXT time series
c      inquire(file=txt_time_series(1:ix+6),exist=existe)
c      if (.not.existe) then
c        open(OUT_TSTXT_UNIT,file=txt_time_series(1:ix+6),status='new',
c     &    form='formatted')
c      endif
c      open(OUT_TSTXT_UNIT,file=txt_time_series(1:ix+6),status='old',
c     &  form='formatted',access='append')
C </Output files> -----------------------------------------------------------
      return
      end subroutine readinout

      subroutine init_tsFile(filename,nnh,npp,Nrows)
        implicit none
        integer, parameter  :: fileUnit = 9000
        integer, intent(in) :: nnh, npp, Nrows
        character(128), intent(in) :: filename
        character      :: ts_file_header*256
        character(8)   :: Nrowss
        character(4)   :: Ncols, LastMode, Nprobes, Nmodes
        character, parameter :: NL = NEW_LINE('A')

        write(Nrowss  ,'(I8)') Nrows
        write(Ncols   ,'(I4)') 1+nnh+npp ! time + energy modes + probe fields
        write(LastMode,'(I4)') nnh-1
        write(Nmodes  ,'(I4)') nnh
        write(Nprobes ,'(I4)') npp
        open(fileUnit,file=filename,status='new',
     &    form='unformatted',access='stream')
        ts_file_header = '# FILE HEADER IS 256 BYTES LONG'
     &    // NL //
     &    '# FOLLOWING DATA_TYPE_BITS: 64'
     &    // NL //
     &    '# Nrows='//Nrowss//', Ncols='//Ncols//','
     &    // NL //
     &    '# Nmodes='//Nmodes//', Nprobes='//Nprobes//','
     &    // NL //
     &    '# time Ek_m(0:'//LastMode//') phasept(1:'//Nprobes//')'
     &    // NL
        ts_file_header(255:256) = NL
        write(fileUnit) ts_file_header
        close(fileUnit)
      return
      end subroutine init_tsFile

c========================================================================

      function ran1(idum)
      integer idum,ia,IM,IQ,IR,NTAB,NDIV
      real*8 ran1,AM,EPS,RNMX
      parameter (ia=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     &     NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      integer j,k,iv(NTAB),iy
      save iv,iy
      data iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
         idum=max(-idum,1)
         do 11 j=NTAB+8,1,-1
            k=idum/IQ
            idum=ia*(idum-k*IQ)-IR*k
            if (idum.lt.0) idum=idum+IM
            if (j.le.NTAB) iv(j)=idum
 11      continue
         iy=iv(1)
      endif
      k=idum/IQ
      idum=ia*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      end
C  (C) Copr. 1986-92 Numerical Recipes Software 0!",1'1.:<%:.

c=======================================================================

      subroutine grid(ccr,ccz,npp,nrp,ntp,nzp,Ftr,axx,azz,nr,nz,nn)
      implicit none
      integer nr,nx,nz,nn,nn1,nrpz
      real*8 , parameter :: pi = dacos(-1d0)
      real*8 cx(0:2*nr+1),cz(0:nz),Chr(0:2*nr+1),Chz(0:nz)
      real*8 ax1(0:2*nr+1,0:2*nr+1),az1(0:nz,0:nz)
      real*8 cc(0:2*nr+nz-1),ccr(0:2*nr+1),ccz(0:nz),angle
      real*8 Ftr(0:nn-1,10),axx(0:2*nr+1,10),azz(0:nz,10),px,pnj
      integer i,j,k,m,n,isig,nrp(10),ntp(10),nzp(10),npp
      nx=2*nr+1 ; nn1=nn-1 ; nrpz=nx+nz-2

C <Clenshaw-Curtis weights> -------------------------------------------------
C     Auxiliary matrices for the computation of the Clenshaw-Curtis weights
      cx(0)=2d0 ; cx(nx)=2d0
      do i=1,nx-1
         cx(i)=1d0
      enddo
      cz(0)=2d0 ; cz(nz)=2d0
      do j=1,nz-1
         cz(j)=1d0
      enddo
c     Chebishev collocation -> spectral, axial
      do n=0,nz
         do j=0,nz
            az1(j,n)=2d0*dcos(pi*dfloat(n*j)/dfloat(nz))/
     &           (dfloat(nz)*cz(j)*cz(n))
         enddo
      enddo
c     Chebishev collocation -> spectral, radial
      do i=0,nx
         do m=0,nx
            ax1(i,m)=2d0*dcos(pi*dfloat(i*m)/dfloat(nx))/
     &           (dfloat(nx)*cx(i)*cx(m))
         enddo
      enddo
C     Spectral weights
      cc(0)=1d0 ; cc(1)=0.5d0 ; isig=1
      do i=1,max(nr,nz/2)
         isig=-isig
         cc(2*i)=-1d0/dfloat(4*i*i-1)
         cc(2*i+1)=dfloat(isig*(2*i+1)-1)/dfloat(4*i*(i+1))
      enddo
C     Axial Clenshaw-Curtis weights
      do k=0,nz
         ccz(k)=0.d0
         do j=0,nz/2
            ccz(k)=ccz(k)+az1(2*j,k)*cc(2*j)
         enddo
         ccz(k)=2d0*ccz(k)
      enddo
C     Radial Clenshaw-Curtis weights
      do k=0,nx
         ccr(k)=0.d0
         do j=0,nx
            ccr(k)=ccr(k)+ax1(j,k)*cc(j)
         enddo
      enddo
      do k=0,nr
         ccr(k)=ccr(k)-ccr(nx-k)
      enddo
C </Clenshaw-Curtis weights> ------------------------------------------------

C <AUXILIARY MATRICES> for the transformation from spectral to physical -----
      do j=1,npp
         angle=dacos(1.d-2*dfloat(nrp(j)))
         do i=0,nx
            Chr(i)=dcos(dfloat(i)*angle)
         enddo
         angle=dacos(2.d-2*dfloat(nzp(j))-1d0)
         do i=0,nz
            Chz(i)=dcos(dfloat(i)*angle)
         enddo
         Ftr(0,j)=1d0
         px=2d0*pi*dfloat(ntp(j))/360.d0
         do n=1,nn/2-1
            pnj=px*dfloat(n)
            Ftr(2*n-1,j)=2d0*dcos(pnj) ; Ftr(2*n,j)=-2d0*dsin(pnj)
         enddo
         Ftr(nn-1,j)=dcos(px*dfloat(nn/2))
         call dgemv('T',nx+1,nx+1,1d0,ax1(0,0),nx+1,Chr(0),1,
     &        0.d0,axx(0,j),1)
         call dgemv('T',nz+1,nz+1,1d0,az1(0,0),nz+1,Chz(0),1,
     &        0.d0,azz(0,j),1)
      enddo
C </AUXILIARY MATRICES> -----------------------------------------------------

      return
      end subroutine grid

c=======================================================================

      subroutine explicit_terms(nlvr,nlvt,nlvz,wk,z,stage,nr,nz,nn)
      implicit none
      integer nr,nz,nn,nn1
      integer stage(2),i,j,nrz,nrzt,k,keven,kodd
      real*8  z(0:nz),a
      real*8  nlvr(0:nr,0:nz,0:nn-1,2),wk(0:nr,0:nz,0:nn-1,6)
      real*8  nlvt(0:nr,0:nz,0:nn-1,2),nlvz(0:nr,0:nz,0:nn-1,2)

      nn1=nn-1 ; nrz=(nr+1)*(nz+1) ; nrzt=nrz*nn
c     explicit term in radial Navier-Stokes eq.    -> wk(,,,1)
c     explicit term in azimuthal Navier-Stokes eq. -> wk(,,,2)
c     explicit term in axial Navier-Stokes eq.     -> wk(,,,3)
c     lines with nl: advection term
!$OMP PARALLEL
!$OMP+SHARED(nlvr,nlvt,nlvz,vr,vt,wk,z,r,stage)
!$OMP+PRIVATE(i,j,a)
!$OMP DO
      do j=0,nz
         a=(1d0+z(j))/2d0
         do i=0,nr
            wk(i,j,0,1)=
     &           -2d0*nlvr(i,j,0,stage(1))+nlvr(i,j,0,stage(2))
            wk(i,j,nn1,1)=
     &           -2d0*nlvr(i,j,nn1,stage(1))+nlvr(i,j,nn1,stage(2))

            wk(i,j,0,2)=
     &           -2d0*nlvt(i,j,0,stage(1))+nlvt(i,j,0,stage(2))
            wk(i,j,nn1,2)=
     &           -2d0*nlvt(i,j,nn1,stage(1))+nlvt(i,j,nn1,stage(2))

            wk(i,j,0,3)=-2d0*nlvz(i,j,0,stage(1))+nlvz(i,j,0,stage(2))
            wk(i,j,nn1,3)=
     &           -2d0*nlvz(i,j,nn1,stage(1))+nlvz(i,j,nn1,stage(2))
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL
!$OMP+SHARED(nlvr,nlvt,nlvz,vr,vt,wk,r,stage)
!$OMP+PRIVATE(i,j,k,keven,kodd)
!$OMP DO
      do k=1,nn/2-1
         keven=2*k ; kodd=keven-1
         do j=0,nz
            do i=0,nr
               wk(i,j,kodd,1)=
     &              -2d0*nlvr(i,j,kodd,stage(1))
     &              +nlvr(i,j,kodd,stage(2))
               wk(i,j,keven,1)=
     &              -2d0*nlvr(i,j,keven,stage(1))
     &              +nlvr(i,j,keven,stage(2))

               wk(i,j,kodd,2)=
     &              -2d0*nlvt(i,j,kodd,stage(1))
     &              +nlvt(i,j,kodd,stage(2))
               wk(i,j,keven,2)=
     &              -2d0*nlvt(i,j,keven,stage(1))
     &              +nlvt(i,j,keven,stage(2))

               wk(i,j,kodd,3)=
     &              -2d0*nlvz(i,j,kodd,stage(1))
     &              +nlvz(i,j,kodd,stage(2))
               wk(i,j,keven,3)=
     &              -2d0*nlvz(i,j,keven,stage(1))
     &              +nlvz(i,j,keven,stage(2))
            enddo
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
      return
      end subroutine explicit_terms

c=======================================================================

      subroutine time_stepper(Re,Ro,wf,dt,pnu,rad,alt,
     &     stage,time,tdr,tir,tdi,tii,wk,dz1,dz2,nr,nz,nn,
     &     dr1,lap,z,r,ri,
     &     axd,vxd,bxd,azd,vzd,bzd,axn,vxn,bxn,azn,vzn,bzn,
     &     cx01,cx03,cz01,cz02,cz03,cz11,cz12,cz13,
     &     vr,vt,vz,nlvr,nlvt,nlvz,rotvorz,rotvorr)
c     stage(1) points to N, stage(2) points to N-1 on input, the opposite on output
      implicit none
      integer nr,nz,nn,nnh,nnhp,nn1,nr1,nz1
      integer stage(2),i,j,k,nrz,nrzt,nzt,n
      real*8  dz1(0:nz,0:nz),dr1(0:nr,0:nr,2)
      real*8  dz2(0:nz,0:nz)
      real*8  lap(0:nr,0:nr,2)
      real*8  z(0:nz),r(0:nr),ri(0:nr)
      real*8  axd(1:nr,1:nr,0:nn/2+1),vxd(1:nr,0:nn/2+1),
     &        bxd(1:nr,1:nr,0:nn/2+1),axn(1:nr,1:nr,0:nn/2+1),
     &        vxn(1:nr,0:nn/2+1),bxn(1:nr,1:nr,0:nn/2+1),
     &        azd(1:nz-1,1:nz-1),vzd(1:nz-1),bzd(1:nz-1,1:nz-1),
     &        azn(1:nz-1,1:nz-1),vzn(1:nz-1),bzn(1:nz-1,1:nz-1)
      real*8  cx01(2),cx03(1:nr,2),cz01,cz02,cz03(1:nz-1),
     &        cz11,cz12,cz13(1:nz-1)
      real*8  tdr(0:nn-1,0:nn-1),tdi(0:nn-1,0:nn-1)
      real*8  tir(0:nn-1,0:nn-1),tii(0:nn-1,0:nn-1)
      real*8  vr(0:nr,0:nz,0:nn-1,2)
      real*8  vt(0:nr,0:nz,0:nn-1,2),vz(0:nr,0:nz,0:nn-1,2)
      real*8  nlvr(0:nr,0:nz,0:nn-1,2)
      real*8  nlvt(0:nr,0:nz,0:nn-1,2),nlvz(0:nr,0:nz,0:nn-1,2)
      real*8  rotvorz(0:nr,2,0:nn-1,2),rotvorr(0:nz,0:nn-1,2)
      real*8  wk(0:nr,0:nz,0:nn-1,6)
      real*8  dt,pnu,rad,alt,time,Re,Ro,wf
      real*8  alpha,b,a

      nnh=nn/2 ; nnhp=nnh+1 ; nn1=nn-1 ; nr1=nr+1 ; nz1=nz+1
      nrz=(nr+1)*(nz+1) ; nrzt=nrz*nn ; nzt=(nz+1)*nn

c     hacemos el termino no lineal y la c.c de la presion para el N
      call tnlrot(stage(1),tdr,tir,tdi,tii,dr1,dz1,ri,vr,vt,vz,
     &     nlvr,nlvt,nlvz,rotvorz,rotvorr,rad,alt,wk,nr,nz,nn)

c****************************************************************************
c
c     Calculamos u^/deltat, y la ponemos en wk( , , ,i): 1, vr; 2, vt; 3, vz
c
c     All the explicit terms in Navier-Stokes: advection (nlxx),
c     buoyancy, Coriolis and centrifugal
c
c****************************************************************************
      call explicit_terms(nlvr,nlvt,nlvz,wk,z,stage,nr,nz,nn)

c******************************************************************************
c
c     Calculamos la divergencia de v^/deltat y la ponemos en wk( , , ,4)
c
c******************************************************************************
!$OMP PARALLEL
!$OMP+SHARED(nrzt,wk)
!$OMP+PRIVATE(i)
!$OMP DO
      do i=0,nrzt-1
         wk(i,0,0,4)=0.d0
      enddo
!$OMP END DO
!$OMP END PARALLEL
c     calculamos derivada radial de v^r
      call derivr(wk(0,0,0,1),wk(0,0,0,4),1,2,dr1,rad,nr,nz,nn)
c     a\F1adimos v_r/r(i)-n v_t/r(i)
      alpha=1d0/rad
!$OMP PARALLEL
!$OMP+SHARED(wk,alpha,ri)
!$OMP+PRIVATE(n,j,i)
!$OMP DO
      do j=0,nz
         do i=0,nr
            wk(i,j,0,4)=wk(i,j,0,4)+alpha*ri(i)*wk(i,j,0,1)
            wk(i,j,nn1,4)=wk(i,j,nn1,4)+alpha*ri(i)*wk(i,j,nn1,1)
         enddo
      enddo
!$OMP END DO
!$OMP DO
      do n=1,nnh-1
         do j=0,nz
            do i=0,nr
               wk(i,j,2*n-1,4)=wk(i,j,2*n-1,4)
     &              +alpha*ri(i)*(wk(i,j,2*n-1,1)-
     &              dfloat(n)*wk(i,j,2*n-1,2))
               wk(i,j,2*n,4)=wk(i,j,2*n,4)
     &              +alpha*ri(i)*(wk(i,j,2*n,1)-
     &              dfloat(n)*wk(i,j,2*n,2))
            enddo
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
c     a\F1adimos d_z(v_z)
      alpha=2d0/alt
      do k=0,nn1
         call dgemm('N','T',nr1,nz1,nz1,alpha,wk(0,0,k,3),nr1,
     &        dz1(0,0),nz1,1d0,wk(0,0,k,4),nr1)
      enddo
c******************************************************************************
c
c     Calculamos la presion y la ponemos en wk( , , ,5)
c
c******************************************************************************
c     El segundo miembro esta en wk( , , ,4)
c     Fabricamos el contorno en z=0,1, la derivada z de la presion vale
c     pnu(-2 rotvorz(n)+rotvorz(n-1))
c     -2(nl)_z^{n+}+(nl)_z^{n-1} -2 coriolis_^{n}+coriolis_^{n-1}
c     ESTO ES CERO; condiciones de contorno en las tapas
!$OMP PARALLEL
!$OMP+SHARED(pnu,rotvorz,stage,wk)
!$OMP+PRIVATE(n,i)
!$OMP DO
      do n=0,nn1
         do i=0,nr
            wk(i,0,n,4)=wk(i,0,n,3)
            wk(i,nz,n,4)=wk(i,nz,n,3)
            wk(i,0,n,4)=wk(i,0,n,4)+pnu*(-2d0*rotvorz(i,1,n,stage(1))+
     &           rotvorz(i,1,n,stage(2)))
          wk(i,nz,n,4)=wk(i,nz,n,4)+pnu*(-2d0*rotvorz(i,2,n,stage(1))+
     &           rotvorz(i,2,n,stage(2)))
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
c     condiciones de contorno en la superficie lateral; en r=0, la
c     derivada r de la presion vale pnu(-2 rotvorr(N)+rotvorr(N-1))+
c     -2(nl)_r^{n+}+(nl)_r^{n-1} -2 coriolis_^{n}+coriolis_^{n-1}
!$OMP PARALLEL
!$OMP+SHARED(wk,pnu,rotvorr,stage)
!$OMP+PRIVATE(n,j)
!$OMP DO
      do n=0,nn1
         do j=1,nz-1
            wk(0,j,n,4)=wk(0,j,n,1)
            wk(0,j,n,4)=wk(0,j,n,4)+pnu*(-2d0*rotvorr(j,n,stage(1))+
     &           rotvorr(j,n,stage(2)))
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
      a=0.d0
      call poisp(a,wk(0,0,0,4),wk(0,0,0,5),dz2,lap,axn,vxn,bxn,azn,vzn,
     &     bzn,rad,alt,cx01,cx03,cz01,cz02,cz03,cz11,cz12,cz13,nr,nz,nn)
c     En wk( , , ,5) esta la presion

c*****************************************************************************
c
c     Calculamos u*/deltat=-gradP+u^/deltat+2 v^n/Deltat -v^(n-1)/2Deltatt
c
c*****************************************************************************
c     Calculamos la derivada radial de p y la ponemos en wk( , , ,6); se
c     la restamos a u^ que esta en wk( , , ,1) y la ponemos en wk( , , ,1).
c     Calculamos derivada radial de presion que es par, n=0 par, n=1 impar...
      call derivr(wk(0,0,0,5),wk(0,0,0,6),2,1,dr1,rad,nr,nz,nn)
!$OMP PARALLEL
!$OMP+SHARED(nrzt,wk,vr,stage,dt)
!$OMP+PRIVATE(i)
!$OMP DO
      do i=0,nrzt-1
         wk(i,0,0,1)=wk(i,0,0,1)-wk(i,0,0,6)+
     &        (2d0*vr(i,0,0,stage(1))-0.5d0*vr(i,0,0,stage(2)))/dt
      enddo
!$OMP END DO
!$OMP END PARALLEL
c     Hacemos n/r(i) p
      alpha=1d0/rad
!$OMP PARALLEL
!$OMP+SHARED(wk,ri)
!$OMP+PRIVATE(n,j,i)
!$OMP DO
      do j=0,nz
         do i=0,nr
            wk(i,j,0,6)=0.d0
            wk(i,j,nn1,6)=0.d0
         enddo
      enddo
!$OMP END DO
!$OMP DO
      do n=1,nnh-1
         do j=0,nz
            do i=0,nr
               wk(i,j,2*n-1,6)=ri(i)*wk(i,j,2*n-1,5)*dfloat(n)
               wk(i,j,2*n,6)=ri(i)*wk(i,j,2*n,5)*dfloat(n)
            enddo
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
!$OMP PARALLEL
!$OMP+SHARED(nrzt,wk,rad,vt,stage,dt)
!$OMP+PRIVATE(i)
!$OMP DO
      do i=0,nrzt-1
         wk(i,0,0,2)=wk(i,0,0,2)-wk(i,0,0,6)/rad+
     &       (2d0*vt(i,0,0,stage(1))-0.5d0*vt(i,0,0,stage(2)))/dt
      enddo
!$OMP END DO
!$OMP END PARALLEL
c     hacemos la derivada respecto de z de la presion
      alpha=2d0/alt
      do k=0,nn1
         call dgemm('N','T',nr1,nz1,nz1,alpha,wk(0,0,k,5),nr1,
     &        dz1(0,0),nz1,0.d0,wk(0,0,k,6),nr1)
      enddo
!$OMP PARALLEL
!$OMP+SHARED(nrzt,wk,vz,stage,dt)
!$OMP+PRIVATE(i)
!$OMP DO
      do i=0,nrzt-1
         wk(i,0,0,3)=wk(i,0,0,3)-wk(i,0,0,6)+
     &        (2d0*vz(i,0,0,stage(1))-0.5d0*vz(i,0,0,stage(2)))/dt
      enddo
!$OMP END DO
!$OMP END PARALLEL
c     Ahora tenemos en wk( , , ,i) 1-vr*/deltat, 2-vt*/deltat, 3-vz*/deltat

c**************************************************************************
c
c     Compute u*_l=ur+utheta. Divided by deltat goes to wk( , , ,4)
c     Compute u*_p=ur-utheta. Divided by deltat goes to wk( , , ,6)
c     Resolvemos los poisson de las velocidades
c
c**************************************************************************

c     boundary conditions for the velocity field
      call bcvel(time,r,wk(0,0,0,4),wk(0,0,0,6),wk(0,0,0,3),rad,alt,
     &          nr,nz,nn,Re,Ro,wf,pnu,dz1)
!$OMP PARALLEL
!$OMP+SHARED(alpha,wk)
!$OMP+PRIVATE(n,i)
!$OMP DO
      do n=0,nn1
         do i=0,nr
            alpha=wk(i,0,n,4)
            wk(i,0,n,4)=alpha+wk(i,0,n,6)
            wk(i,0,n,6)=alpha-wk(i,0,n,6)
            alpha=wk(i,nz,n,4)
            wk(i,nz,n,4)=alpha+wk(i,nz,n,6)
            wk(i,nz,n,6)=alpha-wk(i,nz,n,6)
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL
!$OMP+SHARED(alpha,wk)
!$OMP+PRIVATE(n,j)
!$OMP DO
      do n=0,nn1
         do j=1,nz-1
            alpha=wk(0,j,n,4)
            wk(0,j,n,4)=alpha+wk(0,j,n,6)
            wk(0,j,n,6)=alpha-wk(0,j,n,6)
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL

c     inner values for the velocity field
!$OMP PARALLEL
!$OMP+SHARED(wk,pnu)
!$OMP+PRIVATE(n,j,i)
!$OMP DO
      do n=0,nn1
         do j=1,nz-1
            do i=1,nr
               wk(i,j,n,4)=-(wk(i,j,n,1)+wk(i,j,n,2))/pnu
               wk(i,j,n,6)=-(wk(i,j,n,1)-wk(i,j,n,2))/pnu
               wk(i,j,n,3)=-wk(i,j,n,3)/pnu
            enddo
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL

c     solving the Helmholtz equation
      a=1.5d0/(pnu*dt)
      call poisul(a,wk(0,0,0,4),wk(0,0,0,1),dz2,lap,axd,vxd,bxd,azd,
     &     vzd,bzd,rad,alt,nr,nz,nn)
c     En wk( , , ,1) esta u_l=ur+utheta
      a=1.5d0/(pnu*dt)
      call poisup(a,wk(0,0,0,6),wk(0,0,0,2),dz2,lap,axd,vxd,bxd,azd,
     &     vzd,bzd,rad,alt,nr,nz,nn)
c     En wk( , , ,2) esta u_p=ur-utheta
c     Calculamos la velocidades vr,vt
!$OMP PARALLEL
!$OMP+SHARED(vr,stage,wk,nrzt)
!$OMP+PRIVATE(i)
!$OMP DO
      do i=0,nrzt-1
         vr(i,0,0,stage(2))=0.5d0*(wk(i,0,0,1)+wk(i,0,0,2))
         vt(i,0,0,stage(2))=0.5d0*(wk(i,0,0,1)-wk(i,0,0,2))
      enddo
!$OMP END DO
!$OMP END PARALLEL
c     calculamos vz
      a=1.5d0/(pnu*dt)
      call poist(a,wk(0,0,0,3),vz(0,0,0,stage(2)),dz2,lap,axd,vxd,bxd,
     &     azd,vzd,bzd,rad,alt,nr,nz,nn)

c *************************************************************************
c
c                       Paso corrector
c
c**************************************************************************

c**************************************************************************
c
c     Calculamos la divergencia de v*  y la ponemos en wk( , , ,1)
c
c**************************************************************************
!$OMP PARALLEL
!$OMP+SHARED(nrzt,wk)
!$OMP+PRIVATE(i)
!$OMP DO
      do i=0,nrzt-1
         wk(i,0,0,1)=0.d0
      enddo
!$OMP END DO
!$OMP END PARALLEL
c     calculamos derivada radial de v^r
      call derivr(vr(0,0,0,stage(2)),wk(0,0,0,1),1,2,dr1,rad,nr,nz,nn)
c     a\F1adimos v_r/r(i)-n v_t/r(i)
      alpha=1d0/rad
!$OMP PARALLEL
!$OMP+SHARED(wk,alpha,ri,vr,stage)
!$OMP+PRIVATE(j,i)
!$OMP DO
      do j=0,nz
         do i=0,nr
            wk(i,j,0,1)=wk(i,j,0,1)+alpha*ri(i)*vr(i,j,0,stage(2))
            wk(i,j,nn1,1)=wk(i,j,nn1,1)+alpha*ri(i)*vr(i,j,nn1,stage(2))
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL
!$OMP+SHARED(wk,alpha,ri,vr,stage,vt)
!$OMP+PRIVATE(n,j,i)
!$OMP DO
      do n=1,nnh-1
         do j=0,nz
            do i=0,nr
               wk(i,j,2*n-1,1)=wk(i,j,2*n-1,1)
     &              +alpha*ri(i)*(vr(i,j,2*n-1,stage(2))
     &              -dfloat(n)*vt(i,j,2*n-1,stage(2)))

               wk(i,j,2*n,1)=wk(i,j,2*n,1)
     &              +alpha*ri(i)*(vr(i,j,2*n,stage(2))
     &              -dfloat(n)*vt(i,j,2*n,stage(2)))
            enddo
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
c     a\F1adimos d_z(v_z)
      alpha=2d0/alt
      do k=0,nn1
         call dgemm('N','T',nr1,nz1,nz1,alpha,vz(0,0,k,stage(2)),nr1,
     &        dz1(0,0),nz1,1d0,wk(0,0,k,1),nr1)
      enddo

c******************************************************************************
c
c     Calculamos la fi y la ponemos en wk( , , ,2)
c
c******************************************************************************
c     El segundo miembro esta en wk( , , ,4)
c     Fabricamos el contorno
c     en z=0,1, la derivada z de la presion vale 0
c     condiciones de contorno en las tapas
!$OMP PARALLEL
!$OMP+SHARED(wk)
!$OMP+PRIVATE(n,i)
!$OMP DO
      do n=0,nn1
         do i=0,nr
            wk(i,0,n,1)=0.d0
            wk(i,nz,n,1)=0.d0
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
c     condiciones de contorno en la superficie lateral
c     en r=0, la derivada r de la presion vale 0
!$OMP PARALLEL
!$OMP+SHARED(wk)
!$OMP+PRIVATE(n,j)
!$OMP DO
      do n=0,nn1
         do j=1,nz-1
            wk(0,j,n,1)=0.d0
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
      a=0.d0
      call poisp(a,wk(0,0,0,1),wk(0,0,0,2),dz2,lap,axn,vxn,bxn,azn,vzn
     &   ,bzn,rad,alt,cx01,cx03,cz01,cz02,cz03,cz11,cz12,cz13,nr,nz,nn)
c     En wk( , , ,2) esta la  fi

c*****************************************************************************
c
c     Calulamos v^(n+1)=-grad fi+ v*
c
c*****************************************************************************
c     calculamos la derivada radial de fi y la ponemos en wk( , , ,6); se la
c     restamos a u^que esta en vr( , , ,stage(2)) y la ponemos en vr( , , ,stage(2))
c    Calculamos la derivada radial de presion que es par, n=0 par, n=1 impar...
      call derivr(wk(0,0,0,2),wk(0,0,0,6),2,1,dr1,rad,nr,nz,nn)
      b=-1d0
      call daxpy(nrzt,b,wk(0,0,0,6),1,vr(0,0,0,stage(2)),1)
c     Hacemos n/r(i) p
      alpha=1d0/rad
!$OMP PARALLEL
!$OMP+SHARED(wk,ri)
!$OMP+PRIVATE(n,j,i)
!$OMP DO
      do j=0,nz
         do i=0,nr
            wk(i,j,0,6)=0.d0
            wk(i,j,nn1,6)=0.d0
         enddo
      enddo
!$OMP END DO
!$OMP DO
      do n=1,nnh-1
         do j=0,nz
            do i=0,nr
               wk(i,j,2*n-1,6)=ri(i)*wk(i,j,2*n-1,2)*dfloat(n)
               wk(i,j,2*n,6)=ri(i)*wk(i,j,2*n,2)*dfloat(n)
            enddo
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
      b=-alpha
      call daxpy(nrzt,b,wk(0,0,0,6),1,vt(0,0,0,stage(2)),1)
c     hacemos la derivada respecto de z de la presion
      alpha=2d0/alt
      do k=0,nn1
         call dgemm('N','T',nr1,nz1,nz1,alpha,wk(0,0,k,2),nr1,
     &        dz1(0,0),nz1,0.d0,wk(0,0,k,6),nr1)
      enddo
      b=-1d0
      call daxpy(nrzt,b,wk(0,0,0,6),1,vz(0,0,0,stage(2)),1)

      return
      end subroutine time_stepper

c=======================================================================

      subroutine TnlROT(stage,tdr,tir,tdi,tii,dr1,dz1,ri,vr,vt,vz,
     &     nlvr,nlvt,nlvz,rotvorz,rotvorr,rad,alt,wk,nr,nz,nn)
      implicit none
      integer nr,nz,nn,nnh,nn1,nr1,nz1
      integer stage,i,j,k,nrz,nrzt,nzt,n
      real*8  dz1(0:nz,0:nz),dr1(0:nr,0:nr,2),ri(0:nr)
      real*8  tdr(0:nn-1,0:nn-1),tdi(0:nn-1,0:nn-1),tir(0:nn-1,0:nn-1)
      real*8  tii(0:nn-1,0:nn-1)
      real*8  vr(0:nr,0:nz,0:nn-1,2),vt(0:nr,0:nz,0:nn-1,2)
      real*8  vz(0:nr,0:nz,0:nn-1,2)
      real*8  nlvr(0:nr,0:nz,0:nn-1,2),nlvt(0:nr,0:nz,0:nn-1,2)
      real*8  nlvz(0:nr,0:nz,0:nn-1,2)
      real*8  rotvorz(0:nr,2,0:nn-1,2),rotvorr(0:nz,0:nn-1,2)
      real*8  wk(0:nr,0:nz,0:nn-1,6),ep(0:nr,0:nz,0:nn-1)
      real*8  rad,alt,alpha,con

      nnh=nn/2 ; nn1=nn-1 ; nr1=nr+1 ; nz1=nz+1
      nrz=(nr+1)*(nz+1) ; nrzt=nrz*nn ; nzt=(nz+1)*nn
      con=1d0/(dfloat(nn))

c     transformamos Fourier vr
      call dgemm('N','T',nrz,nn,nn,1d0,vr(0,0,0,stage),nrz,tdr(0,0),nn,
     $     0.d0,wk(0,0,0,1),nrz)

c     transformamos Fourier vz
      call dgemm('N','T',nrz,nn,nn,1d0,vz(0,0,0,stage),nrz,tdr(0,0),nn,
     $     0.d0,wk(0,0,0,3),nrz)

c     transformamos Fourier vt
      call dgemm('N','T',nrz,nn,nn,1d0,vt(0,0,0,stage),nrz,tdi(0,0),nn,
     $     0.d0,wk(0,0,0,2),nrz)

c     Definimos una variable que contiene  v_theta/r
!$OMP PARALLEL
!$OMP+SHARED(wk,ri,rad)
!$OMP+PRIVATE(k,j,i)
!$OMP DO
      do k=0,nn1
         do j=0,nz
            do i=0,nr
               wk(i,j,k,4)=wk(i,j,k,2)*ri(i)/rad
            enddo
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
c***********************************************************************
c
c     Hacemos el termino nolineal de la vz
c
c***********************************************************************
c     calculamos la derivada radial de vz que es par, n=0 par, n=1 impar ....
      call derivr(vz(0,0,0,stage),wk(0,0,0,5),2,1,dr1,rad,nr,nz,nn)

c     transformamos fourier la derivada radial de vz
      call dgemm('N','T',nrz,nn,nn,1d0,wk(0,0,0,5),nrz,tdr(0,0),nn,
     $     0.d0,wk(0,0,0,6),nrz)

c     hacemos  vr x partial_r (vz)
!$OMP PARALLEL
!$OMP+SHARED(nrzt,ep,wk)
!$OMP+PRIVATE(i)
!$OMP DO
      do i=0,nrzt-1
         ep(i,0,0)=wk(i,0,0,1)*wk(i,0,0,6)
      enddo
!$OMP END DO
!$OMP END PARALLEL
c     calculamos la derivada acimutal de vz; nn es par
      do j=0,nrz-1
        wk(j,0,0,5)=0.d0
        wk(j,0,nn1,5)=0.d0
      enddo
!$OMP PARALLEL
!$OMP+SHARED(nrz,wk,vz,stage)
!$OMP+PRIVATE(n)
!$OMP DO
      do n=1,nnh-1
         do j=0,nrz-1
            wk(j,0,2*n-1,5)=dfloat(n)*vz(j,0,2*n-1,stage)
            wk(j,0,2*n,  5)=dfloat(n)*vz(j,0,2*n,  stage)
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
c     transformamos fourier la derivada acimutal de vz
      call dgemm('N','T',nrz,nn,nn,1d0,wk(0,0,0,5),nrz,tdi(0,0),nn,
     $     0.d0,wk(0,0,0,6),nrz)

c     hacemos  v_theta/r x partial_theta (vz)
!$OMP PARALLEL
!$OMP+SHARED(nrzt,ep,wk)
!$OMP+PRIVATE(i)
!$OMP DO
      do i=0,nrzt-1
         ep(i,0,0)=ep(i,0,0)+wk(i,0,0,4)*wk(i,0,0,6)
      enddo
!$OMP END DO
!$OMP END PARALLEL
c     hacemos la derivada respecto de z de vz; parcial z de vz(N+1)
      alpha=2d0/alt
      do k=0,nn1
         call dgemm('N','T',nr1,nz1,nz1,alpha,vz(0,0,k,stage),nr1,
     $        dz1(0,0),nz1,0.d0,wk(0,0,k,5),nr1)
      enddo
c     transformamos Fourier la derivada respecto de z de vz
      call dgemm('N','T',nrz,nn,nn,1d0,wk(0,0,0,5),nrz,tdr(0,0),nn,
     $     0.d0,wk(0,0,0,6),nrz)

c     hacemos  v_z x partial_z (vz)
!$OMP PARALLEL
!$OMP+SHARED(nrzt,ep,wk)
!$OMP+PRIVATE(i)
!$OMP DO
      do i=0,nrzt-1
         ep(i,0,0)=ep(i,0,0)+wk(i,0,0,3)*wk(i,0,0,6)
      enddo
!$OMP END DO
!$OMP END PARALLEL
c     antitransformamos Fourier
      call dgemm('N','T',nrz,nn,nn,1d0,ep(0,0,0),nrz,tir,nn,
     $     0.d0,nlvz(0,0,0,stage),nrz)

c***********************************************************************
c
c     Hacemos el termino nolineal de la vr
c
c***********************************************************************

c     calculamos la derivada radial de vr que es impar, n=0 impar, n=1 par ....
      call derivr(vr(0,0,0,stage),wk(0,0,0,5),1,2,dr1,rad,nr,nz,nn)

c     transformamos fourier la derivada radial de vr
      call dgemm('N','T',nrz,nn,nn,1d0,wk(0,0,0,5),nrz,tdr(0,0),nn,
     $     0.d0,wk(0,0,0,6),nrz)

c     hacemos  vr x partial_r (vr)
!$OMP PARALLEL
!$OMP+SHARED(nrzt,ep,wk)
!$OMP+PRIVATE(i)
!$OMP DO
      do i=0,nrzt-1
         ep(i,0,0)=wk(i,0,0,1)*wk(i,0,0,6)
      enddo
!$OMP END DO
!$OMP END PARALLEL
c     calculamos la derivada acimutal de vr; nn es par
!$OMP PARALLEL
!$OMP+SHARED(nrz,wk,vr,stage)
!$OMP+PRIVATE(n,j)
!$OMP DO
      do j=0,nrz-1
        wk(j,0,0,5)=0.d0
        wk(j,0,nn1,5)=0.d0
      enddo
!$OMP END DO
!$OMP DO
      do n=1,nnh-1
         do j=0,nrz-1
            wk(j,0,2*n-1,5)=dfloat(n)*vr(j,0,2*n-1,stage)
            wk(j,0,2*n,  5)=dfloat(n)*vr(j,0,2*n,  stage)
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
c     transformamos Fourier la derivada acimutal de vr
      call dgemm('N','T',nrz,nn,nn,1d0,wk(0,0,0,5),nrz,tdi(0,0),nn,
     $     0.d0,wk(0,0,0,6),nrz)

c     hacemos  v_theta/r x (partial_theta (vr)-v_theta)
!$OMP PARALLEL
!$OMP+SHARED(nrzt,ep,wk)
!$OMP+PRIVATE(i)
!$OMP DO
      do i=0,nrzt-1
         ep(i,0,0)=ep(i,0,0)+wk(i,0,0,4)*(wk(i,0,0,6)-wk(i,0,0,2))
      enddo
!$OMP END DO
!$OMP END PARALLEL
c     hacemos la derivada respecto de z de vr; PARCIAL z DE vr(N+1)
      alpha=2d0/alt
      do k=0,nn1
         call dgemm('N','T',nr1,nz1,nz1,alpha,vr(0,0,k,stage),nr1,
     $        dz1(0,0),nz1,0.d0,wk(0,0,k,5),nr1)
      enddo
c     transformamos Fourier la derivada respecto de z de vr
      call dgemm('N','T',nrz,nn,nn,1d0,wk(0,0,0,5),nrz,tdr(0,0),nn,
     $     0.d0,wk(0,0,0,6),nrz)

c     hacemos  v_z x partial_z (vr)
!$OMP PARALLEL
!$OMP+SHARED(nrzt,ep,wk)
!$OMP+PRIVATE(i)
!$OMP DO
      do i=0,nrzt-1
         ep(i,0,0)=ep(i,0,0)+wk(i,0,0,3)*wk(i,0,0,6)
      enddo
!$OMP END DO
!$OMP END PARALLEL
c     antitransformamos Fourier
      call dgemm('N','T',nrz,nn,nn,1d0,ep(0,0,0),nrz,tir,nn,
     $     0.d0,nlvr(0,0,0,stage),nrz)

c***********************************************************************
c
c     Hacemos el termino nolineal de la v_theta
c
c***********************************************************************

c     derivada radial de v_theta que es impar, n=0 impar, n=1 par ....
      call derivr(vt(0,0,0,stage),wk(0,0,0,5),1,2,dr1,rad,nr,nz,nn)

c     transformamos Fourier la derivada radial de v_theta
      call dgemm('N','T',nrz,nn,nn,1d0,wk(0,0,0,5),nrz,tdi(0,0),nn,
     $     0.d0,wk(0,0,0,6),nrz)

c     hacemos  vr x partial_r (v_theta)
!$OMP PARALLEL
!$OMP+SHARED(nrzt,ep,wk)
!$OMP+PRIVATE(i)
!$OMP DO
      do i=0,nrzt-1
         ep(i,0,0)=wk(i,0,0,1)*wk(i,0,0,6)
      enddo
!$OMP END DO
!$OMP END PARALLEL
c     calculamos la derivada acimutal de v_theta; nn es par
!$OMP PARALLEL
!$OMP+SHARED(nrz,wk,vt,stage)
!$OMP+PRIVATE(n,j)
!$OMP DO
      do j=0,nrz-1
        wk(j,0,0,5)=0.d0
        wk(j,0,nn1,5)=0.d0
      enddo
!$OMP END DO
!$OMP DO
      do n=1,nnh-1
         do j=0,nrz-1
            wk(j,0,2*n-1,5)=-dfloat(n)*vt(j,0,2*n-1,stage)
            wk(j,0,2*n,  5)=-dfloat(n)*vt(j,0,2*n,  stage)
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
c     transformamos Fourier la derivada acimutal de vt
      call dgemm('N','T',nrz,nn,nn,1d0,wk(0,0,0,5),nrz,tdr(0,0),nn,
     $     0.d0,wk(0,0,0,6),nrz)

!$OMP PARALLEL
!$OMP+SHARED(nrzt,ep,wk)
!$OMP+PRIVATE(i)
!$OMP DO
      do i=0,nrzt-1
         ep(i,0,0)=ep(i,0,0)+wk(i,0,0,4)*(wk(i,0,0,6)+wk(i,0,0,1))
      enddo
!$OMP END DO
!$OMP END PARALLEL
c     hacemos la derivada respecto de z de v_theta; PARCIAL z DE v_theta(N+1)
      alpha=2d0/alt
      do k=0,nn1
         call dgemm('N','T',nr1,nz1,nz1,alpha,vt(0,0,k,stage),nr1,
     $        dz1(0,0),nz1,0.d0,wk(0,0,k,5),nr1)
      enddo
c     transformamos Fourier la derivada respecto de z de v_theta
      call dgemm('N','T',nrz,nn,nn,1d0,wk(0,0,0,5),nrz,tdi(0,0),nn,
     $     0.d0,wk(0,0,0,6),nrz)

c     hacemos  v_z x partial_z (v_theta)
!$OMP PARALLEL
!$OMP+SHARED(nrzt,ep,wk)
!$OMP+PRIVATE(i)
!$OMP DO
      do i=0,nrzt-1
         ep(i,0,0)=ep(i,0,0)+wk(i,0,0,3)*wk(i,0,0,6)
      enddo
!$OMP END DO
!$OMP END PARALLEL
c     antitransformamos Fourier
      call dgemm('N','T',nrz,nn,nn,1d0,ep(0,0,0),nrz,tii(0,0),nn,
     $     0.d0,nlvt(0,0,0,stage),nrz)

c*****************************************************************************
c     hacemos el rotacional de la vorticidad que se necesita para la
c     c.c. de la presion
c*****************************************************************************

c     Hacemos d_z(v_r)
      alpha=2d0/alt
      do k=0,nn1
         call dgemm('N','T',nr1,nz1,nz1,alpha,vr(0,0,k,stage),nr1,
     $        dz1(0,0),nz1,0.d0,wk(0,0,k,5),nr1)
      enddo
!$OMP PARALLEL
!$OMP+SHARED(stage,wk,rotvorz,ri,rad)
!$OMP+PRIVATE(k,i)
!$OMP DO
      do k=0,nn1
         do i=0,nr
            rotvorz(i,1,k,stage)=wk(i,0,k,5)*ri(i)/rad
            rotvorz(i,2,k,stage)=wk(i,nz,k,5)*ri(i)/rad
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
c     derivada radial de d_z(v_r) que es impar, n=0 impar, n=1 par ....
      call derivr(wk(0,0,0,5),wk(0,0,0,6),1,2,dr1,rad,nr,nz,nn)
!$OMP PARALLEL
!$OMP+SHARED(stage,wk,rotvorz)
!$OMP+PRIVATE(k,i)
!$OMP DO
      do k=0,nn1
         do i=0,nr
            rotvorz(i,1,k,stage)=rotvorz(i,1,k,stage)+wk(i,0,k,6)
            rotvorz(i,2,k,stage)=rotvorz(i,2,k,stage)+wk(i,nz,k,6)
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
c     Hacemos d_z(v_theta)
      alpha=2d0/alt
      do k=0,nn1
         call dgemm('N','T',nr1,nz1,nz1,alpha,vt(0,0,k,stage),nr1,
     $        dz1(0,0),nz1,0.d0,wk(0,0,k,5),nr1)
      enddo
c     calculamos la derivada acimutal de d_z(v_theta); nn es par
!$OMP PARALLEL
!$OMP+SHARED(nrz,wk)
!$OMP+PRIVATE(n,j)
!$OMP DO
      do j=0,nrz-1
        wk(j,0,0,6)=0.d0
        wk(j,0,nn1,6)=0.d0
      enddo
!$OMP END DO
!$OMP DO
      do n=1,nnh-1
         do j=0,nrz-1
            wk(j,0,2*n-1,6)=-dfloat(n)*wk(j,0,2*n-1,5)
            wk(j,0,2*n,  6)=-dfloat(n)*wk(j,0,2*n,  5)
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL
!$OMP+SHARED(rotvorz,stage,wk,ri,rad)
!$OMP+PRIVATE(k,i)
!$OMP DO
      do k=0,nn1
         do i=0,nr
            rotvorz(i,1,k,stage)=
     &        rotvorz(i,1,k,stage)+wk(i,0,k,6)*ri(i)/rad
            rotvorz(i,2,k,stage)=
     &        rotvorz(i,2,k,stage)+wk(i,nz,k,6)*ri(i)/rad
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
c**********************************************************
c     calculamos la derivada acimutal de v_theta; nn es par
!$OMP PARALLEL
!$OMP+SHARED(nrz,wk,vt,stage)
!$OMP+PRIVATE(n,j)
!$OMP DO
      do j=0,nrz-1
        wk(j,0,0,5)=0.d0
        wk(j,0,nn1,5)=0.d0
      enddo
!$OMP END DO
!$OMP DO
      do n=1,nnh-1
         do j=0,nrz-1
            wk(j,0,2*n-1,5)=-dfloat(n)*vt(j,0,2*n-1,stage)
            wk(j,0,2*n,  5)=-dfloat(n)*vt(j,0,2*n,  stage)
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL
!$OMP+SHARED(rotvorr,stage,wk,rad)
!$OMP+PRIVATE(k,j)
!$OMP DO
      do k=0,nn1
         do j=0,nz
            rotvorr(j,k,stage)=wk(0,j,k,5)/(rad*rad)
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL

c     derivada radial de D_theta(v_theta) que es impar, n=0 impar, n=1 par ....
      call derivr(wk(0,0,0,5),wk(0,0,0,6),1,2,dr1,rad,nr,nz,nn)
!$OMP PARALLEL
!$OMP+SHARED(rotvorr,stage,wk,rad)
!$OMP+PRIVATE(k,j)
!$OMP DO
      do k=0,nn1
         do j=0,nz
            rotvorr(j,k,stage)=rotvorr(j,k,stage)+wk(0,j,k,6)/(rad)
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL

c     Hacemos d_z(v_z)
      alpha=2d0/alt
      do k=0,nn1
         call dgemm('N','T',nr1,nz1,nz1,alpha,vz(0,0,k,stage),nr1,
     $        dz1(0,0),nz1,0.d0,wk(0,0,k,5),nr1)
      enddo
c     derivada radial de d_z(v_z) que es par, n=0 par, n=1 impar ....
      call derivr(wk(0,0,0,5),wk(0,0,0,6),2,1,dr1,rad,nr,nz,nn)
!$OMP PARALLEL
!$OMP+SHARED(rotvorr,stage,wk)
!$OMP+PRIVATE(k,j)
!$OMP DO
      do k=0,nn1
         do j=0,nz
            rotvorr(j,k,stage)=rotvorr(j,k,stage)+wk(0,j,k,6)
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL

      return
      end subroutine TnlROT

c=======================================================================

      subroutine enforce_symmetry(iaxisym,vr,vt,vz,nr,nz,nn)

      implicit none
c     input variables
      integer nr,nz,nn,iaxisym
      real*8 vr(0:nr,0:nz,0:nn-1,2)
      real*8 vt(0:nr,0:nz,0:nn-1,2),vz(0:nr,0:nz,0:nn-1,2)
c     local variables
      integer i,j,k,n,jj,isig

C     Computing into the axisymmetric subspace
         if (iaxisym.eq.0) then
!$OMP PARALLEL
!$OMP+SHARED(vr,vt,vz)
!$OMP+PRIVATE(k,j,i,n)
!$OMP DO
            do k=1,nn-1
               do j=0,nz
                  do i=0,nr
                     do n=1,2
                     vr(i,j,k,n)=0.d0
                     vt(i,j,k,n)=0.d0  ; vz(i,j,k,n)=0.d0
                     enddo
                  enddo
               enddo
            enddo
!$OMP END DO
!$OMP END PARALLEL
         endif
C     Computing into the m-Fourier subspace (modes m, 2m, 3m, 4m ...)
         if (abs(iaxisym).gt.1) then
!$OMP PARALLEL
!$OMP+SHARED(iaxisym,vr,vt,vz)
!$OMP+PRIVATE(jj,k,j,i,n)
!$OMP DO
            do jj=1,2*abs(iaxisym)-2
               do k=jj,nn-1,2*abs(iaxisym)
                  do j=0,nz
                     do i=0,nr
                        do n=1,2
                        vr(i,j,k,n)=0.d0
                        vt(i,j,k,n)=0.d0  ; vz(i,j,k,n)=0.d0
                        enddo
                     enddo
                  enddo
               enddo
            enddo
!$OMP END DO
!$OMP END PARALLEL
         endif
C     Computing into the m-rotoreflection subspace
         if (iaxisym.lt.0) then
            do j=0,nz/2-1
               do i=0,nr
                  vr(i,j,0,n)=0.5d0*(vr(i,j,0,n)+vr(i,nz-j,0,n))
                  vr(i,nz-j,0,n)=vr(i,j,0,n)
                  vt(i,j,0,n)=0.5d0*(vt(i,j,0,n)+vt(i,nz-j,0,n))
                  vt(i,nz-j,0,n)=vt(i,j,0,n)
                  vz(i,j,0,n)=0.5d0*(vz(i,j,0,n)-vz(i,nz-j,0,n))
                  vz(i,nz-j,0,n)=-vz(i,j,0,n)
               enddo
            enddo
            do i=0,nr
               vz(i,nz/2,0,n)=0.d0
            enddo
            isig=-1
            do jj=2*abs(iaxisym),nn-1,2*abs(iaxisym)
               do k=jj-1,jj
                  do j=0,nz/2-1
                     do i=0,nr
                        vr(i,j,k,n)=0.5d0*(vr(i,j,k,n)
     &                               +isig*vr(i,nz-j,k,n))
                        vr(i,nz-j,k,n)=isig*vr(i,j,k,n)
                        vt(i,j,k,n)=0.5d0*(vt(i,j,k,n)
     &                               +isig*vt(i,nz-j,k,n))
                        vt(i,nz-j,k,n)=isig*vt(i,j,k,n)
                        vz(i,j,k,n)=0.5d0*(vz(i,j,k,n)
     &                               -isig*vz(i,nz-j,k,n))
                        vz(i,nz-j,k,n)=-isig*vz(i,j,k,n)
                     enddo
                  enddo
                  do i=0,nr
                     if (isig.eq.1) then
                        vz(i,nz/2,k-1,n)=0.d0  ; vz(i,nz/2,k,n)=0.d0
                     else
                        vr(i,nz/2,k-1,n)=0.d0  ; vr(i,nz/2,k,n)=0.d0
                        vt(i,nz/2,k-1,n)=0.d0  ; vt(i,nz/2,k,n)=0.d0
                     endif
                  enddo
               enddo
               isig=-isig
            enddo
            do j=0,nz
               do i=0,nr
                  vr(i,j,nn-1,n)=0.d0
                  vt(i,j,nn-1,n)=0.d0
                  vz(i,j,nn-1,n)=0.d0
               enddo
            enddo
         endif

      return
      end subroutine enforce_symmetry

c=======================================================================

      subroutine tswrite(vr,vt,vz,time,stage,npp,r,axx,azz,Ftr,ccr,
     &   ccz,h,h2,nr,nz,nn,nx,rad,pnu,Re)
      implicit none
      integer, parameter :: OUT_TS_UNIT      = 1001
c      integer, parameter :: OUT_TSTXT_UNIT   = 2001
      integer nr,nx,nz,nn,nnh,nr1,nz1,stage(2)
      integer i,j,k,n,npp
      real*8 , parameter :: pi = dacos(-1d0)
      real*8  vr(0:nr,0:nz,0:nn-1,2), vt(0:nr,0:nz,0:nn-1,2)
      real*8  vz(0:nr,0:nz,0:nn-1,2), vze(0:nx,0:nz,0:nn-1)
      real*8  h(0:nr,0:nz,0:nn-1),h2(0:nr,0:nn-1)
      real*8  ccr(0:2*nr+1),ccz(0:nz),time,r(0:nr)
      real*8  Ftr(0:nn-1,10),axx(0:2*nr+1,10),azz(0:nz,10)
      real*8  rad,pnu,Re
      real*8  pvz(10),ddot,ekm(0:nn/2)
c     local variables

      nnh=nn/2 ; nr1=nr+1 ; nz1=nz+1

C <PHYSICAL> Computing in a equispaced grid ----------------
c Extending vz to [-1,1] in x, preserving parity
!$OMP PARALLEL
!$OMP+SHARED(vze,vz,stage)
!$OMP+PRIVATE(n,j,i)
!$OMP DO
      do n=0,nn-1
         do j=0,nz
            do i=0,nr
               vze(i,j,n)=vz(i,j,n,stage(1))
            enddo
         enddo
      enddo
!$OMP END DO
!$OMP DO
      do j=0,nz
         do i=0,nr
            vze(nx-i,j,0)=vze(i,j,0)
         enddo
      enddo
!$OMP END DO
!$OMP DO
      do n=1,nn-3,2
         do j=0,nz
            do i=0,nr
               vze(nx-i,j,n)  =dfloat((-1)**n)*vze(i,j,n)
               vze(nx-i,j,n+1)=dfloat((-1)**n)*vze(i,j,n+1)
            enddo
         enddo
      enddo
!$OMP END DO
!$OMP DO
      do j=0,nz
         do i=0,nr
            vze(nx-i,j,nn-1)=dfloat((-1)**(nn-2))*vze(i,j,nn-1)
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
c     Computing vz at the npp physical points (nrp,ntp,nzp)
      do k=1,npp
         pvz(k)=0.d0
         do n=0,nn-1
            do j=0,nz
               do i=0,nx
                  pvz(k)=pvz(k)+axx(i,k)*azz(j,k)*Ftr(n,k)*vze(i,j,n)
               enddo
            enddo
         enddo
      enddo

c     The following code reduces the number of operations, but not
c     significantly, and it is not worth implementing.
c      real*8 auxa(0:nx,0:nz,10),auxb(0:nx,10)
c      do k=1,npp
c         pvz(k)=0.d0
c         do i=0,nx
c            auxb(i,k)=0.d0
c            do j=0,nz
c               auxa(i,j,k)=0.d0
c               do n=0,nn-1
c                  auxa(i,j,k)=auxa(i,j,k)+Ftr(n,k)*vze(i,j,n)
c               enddo
c               auxb(i,k)=auxb(i,k)+azz(j,k)*auxa(i,j,k)
c            enddo
c            pvz(k)=pvz(k)+axx(i,k)*auxb(i,k)
c         enddo
c      enddo

C </PHYSICAL> ----------------------------------------------


C <Write Modal Energies ekm (kinetic and thermal)> -------------
C     computing the kinetic energy in rotating frame
!$OMP PARALLEL
!$OMP+SHARED(r,vr,vt,vz,stage,h,rad,pnu,Re)
!$OMP+PRIVATE(i,j,k)
!$OMP DO
      do j=0,nz
         do i=0,nr
            h(i,j,0)=r(i)*(vr(i,j,0,stage(1))**2
     &           +(vt(i,j,0,stage(1))+r(i)*(pnu*Re/rad) )**2
     &           +vz(i,j,0,stage(1))**2)
         enddo
      enddo
!$OMP END DO
!$OMP DO
      do k=1,nnh-1
         do j=0,nz
            do i=0,nr
               h(i,j,k)=0.5d0*r(i)*
     &              (vr(i,j,2*k-1,stage(1))**2+vr(i,j,2*k,stage(1))**2
     &              +vt(i,j,2*k-1,stage(1))**2+vt(i,j,2*k,stage(1))**2
     &              +vz(i,j,2*k-1,stage(1))**2+vz(i,j,2*k,stage(1))**2)
            enddo
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL
!$OMP+SHARED(h,ccz,h2,ccr,ekm)
!$OMP+PRIVATE(k)
!$OMP DO
      do k=0,nnh-1
         call dgemv('N',nr+1,nz+1,1d0,h(0,0,k),nr+1,ccz(0),1,
     &        0.d0,h2(0,k),1)
         ekm(k)=ddot(nr+1,h2(0,k),1,ccr(0),1)
      enddo

!$OMP END DO
!$OMP END PARALLEL
C </Write Modal Energies and Probes Local Results> ---------------------------
      write(OUT_TS_UNIT) time, (ekm(i),i=0,nnh-1), (pvz(i),i=1,npp)
c  134 format (ES23.15e3,2x,ES23.15e3,2x,ES23.15e3)
c      write(OUT_TSTXT_UNIT,134) time, (ekm(i),i=0,nnh-1),(pvz(i),i=1,npp)

      return
      end subroutine tswrite
C
C
