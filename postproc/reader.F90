!--------------------------- Restart Section ---------------------------
module restart
  use iso_c_binding, only: c_double, c_int
  implicit none
  integer, parameter :: dp=kind(0d0)
  integer  :: nr,nx,nz,nn ! nr is radius (0,1], nx is the duplication of nr
  integer  :: stage(2)    ! pointers to solutions at N and N+1
  real(dp) :: dt,tps
  real(dp) :: Bo,Re,Ro,wf,gama,reg
  real(dp),allocatable,dimension(:,:,:,:) ::  vr, vt, vz
  character(256) :: header
  contains
    subroutine read_restart(restart)
      implicit none
      integer, parameter :: IN_UNIT = 1000
      character(*), intent(in) :: restart
      if ( allocated(vr) ) call dealloc()
      open(unit=IN_UNIT,file=restart,status='old',form='unformatted')
      read(IN_UNIT) nr,nz,nn,dt,tps,stage(1),stage(2),Bo,Re,Ro,wf,gama,reg
      nx = 2*nr+1
      call alloc()
      call read_field(nr,nz,nn,2,IN_UNIT,vr)
      call read_field(nr,nz,nn,2,IN_UNIT,vt)
      call read_field(nr,nz,nn,2,IN_UNIT,vz)
      close(IN_UNIT)
      return
    end subroutine read_restart

    subroutine write_restart(restart)
      implicit none
      integer, parameter :: OUT_UNIT = 1001
      character(*), intent(in) :: restart
      open(unit=OUT_UNIT,file=restart,status='new',form='unformatted')
      write(OUT_UNIT) nr,nz,nn,dt,tps,stage(1),stage(2),Bo,Re,Ro,wf,gama,reg
      call write_field(nr,nz,nn,2,OUT_UNIT,vr)
      call write_field(nr,nz,nn,2,OUT_UNIT,vt)
      call write_field(nr,nz,nn,2,OUT_UNIT,vz)
      close(OUT_UNIT)
      return
    end subroutine write_restart

    subroutine alloc()
      implicit none
      allocate(vr(0:nr,0:nz,0:nn-1,2),source=0d0)
      allocate(vt(0:nr,0:nz,0:nn-1,2),source=0d0)
      allocate(vz(0:nr,0:nz,0:nn-1,2),source=0d0)
      return
    end subroutine alloc

    subroutine dealloc()
      implicit none
      deallocate(vr,vt,vz)
      return
    end subroutine dealloc
end module restart
!------------------------- END Restart Section -------------------------

!------------------------- Time-Series Section -------------------------
module timeseries
  use iso_c_binding, only: c_double, c_int
  implicit none
  integer, parameter :: dp=kind(0d0)
  real(dp),allocatable,dimension(:)   ::  time
  real(dp),allocatable,dimension(:,:) ::  ek, phasept
  integer         :: Nrows, Nprobes, LastMode
  character(256)  :: header
  contains
    subroutine read_timeseries(tsFile)
      implicit none
      integer, parameter :: IN_UNIT = 1000
      character(*), intent(in) :: tsFile
      integer :: i, j
      logical :: existe
      if ( allocated(time) ) call dealloc()
      inquire(file=tsFile,exist=existe)
      if (existe) then
        open(unit=IN_UNIT,file=tsFile,status='old',form='unformatted')!,access='stream')
        read(IN_UNIT) header
        read(IN_UNIT) Nrows, LastMode, Nprobes
        call alloc()
        do i=1,Nrows
          read(IN_UNIT) time(i), (ek(i,j),j=0,LastMode), (phasept(i,j),j=1,Nprobes)
        enddo
        close(IN_UNIT)
      else
        print*,'File ',tsFile,' does not exist.'
      endif
      return
    end subroutine read_timeseries

    subroutine alloc()
      implicit none
      allocate(time   (Nrows)           ,source=0d0)
      allocate(ek     (Nrows,0:LastMode),source=0d0)
      allocate(phasept(Nrows,  Nprobes ),source=0d0)
      return
    end subroutine alloc

    subroutine dealloc()
      implicit none
      deallocate(time,ek,phasept)
      return
    end subroutine dealloc
end module timeseries
!----------------------- END Time-Series Section -----------------------

!-------------------------- Other Subroutines --------------------------

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

subroutine write_field(nr,nz,nn,nstage,out_unit,field)
  use iso_c_binding, only: c_double, c_int
  implicit none
  integer, parameter :: dp=kind(0.d0)
  integer,  intent(in)    :: nr,nz,nn,nstage
  integer,  intent(in)    :: out_unit
  real(dp), intent(inout) :: field(0:nr,0:nz,0:nn-1,nstage)
  integer :: nrk,nzk,nnk,nstagek
  write(out_unit)                                 &
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
end subroutine write_field
