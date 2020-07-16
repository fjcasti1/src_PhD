ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c $Id: spectral_tools_fcyl.f 41 2009-11-19 02:52:24Z marques $
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine coef(dr1,dz1,cx01,cx03,cz01,cz02,cz03,cz11,cz12,cz13,
     &                nr,nz)
      implicit none
      integer nz,nr,i
      real*8 dz1(0:nz,0:nz),dr1(0:nr,0:nr,2)
      real*8 cx01(2),cx03(1:nr,2),cz01,cz02,cz03(1:nz-1),
     $     cz11,cz12,cz13(1:nz-1),d

c     El 2 es para el par y el 1 para el impar
      cx01(2)=1d0/dr1(0,0,2)
      cx01(1)=1d0/dr1(0,0,1)
      do i=1,nr
         cx03(i,2)=-dr1(0,i,2)/dr1(0,0,2)
         cx03(i,1)=-dr1(0,i,1)/dr1(0,0,1)
      enddo

      d=dz1(nz,nz)*dz1(0,0)-dz1(0,nz)*dz1(nz,0)
      cz01=dz1(nz,nz)/d
      cz02=-dz1(0,nz)/d
      do i=1,nz-1
         cz03(i)=(dz1(0,nz)*dz1(nz,i)-dz1(nz,nz)*dz1(0,i))/d
      enddo

      cz11=-dz1(nz,0)/d
      cz12=dz1(0,0)/d
      do i=1,nz-1
         cz13(i)=(dz1(nz,0)*dz1(0,i)-dz1(0,0)*dz1(nz,i))/d
      enddo

      return
      end subroutine coef

c=======================================================================

      subroutine deriv(dx1,dz1,dx2,dz2,dr1,dr2,lap,M,x,z,r,ri,riq,nr,nz)
      implicit none
      integer nr,nx,nz
      real*8 dx1(0:2*nr+1,0:2*nr+1),dz1(0:nz,0:nz),dr1(0:nr,0:nr,2),
     $       dx2(0:2*nr+1,0:2*nr+1),dz2(0:nz,0:nz),dr2(0:nr,0:nr,2),
     $       lap(0:nr,0:nr,2), M(0:nr,0:nr,2)
      real*8 pi,x(0:2*nr+1),cx(0:2*nr+1),z(0:nz),cz(0:nz)
      real*8 r(0:nr), ri(0:nr), riq(0:nr), sign
      integer i,j

      nx=2*nr+1
      pi=dacos(-1d0)

***********************************************************************
*     X                                                               *
***********************************************************************

      do i=0,nx
         x(i)=dcos(pi*dfloat(i)/dfloat(nx))
      enddo
      do i=0,nr
         r(i)=dcos(pi*dfloat(i)/dfloat(nx))
         ri(i)=1d0/r(i)
         riq(i)=ri(i)**2d0
      enddo

      cx(0)=2d0
      cx(nx)=2d0
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

      dx1(0,0)=(2d0*nx*nx+1d0)/6.d0
      dx1(nx,nx)=-(2d0*nx*nx+1d0)/6.d0
      do i=1,nx-1
         dx1(i,i)=-x(i)/(2.0d0*(1d0-x(i)*x(i)))
      enddo

***************************

      sign=-1d0
      do j=0,nx
         do i=1,nx-1
            if(i.ne.j) then
              dx2(i,j)=sign*(x(i)*x(i)+x(i)*x(j)-2d0)/
     $             (cx(j)*(1d0-x(i)**2d0)*(x(i)-x(j))**2)
            endif
            sign=-sign
         enddo
         if(mod(nx,2).ne.0) sign=-sign
      enddo

      dx2(0,0)=(nx**4-1d0)/15.d0
      dx2(nx,nx)=(nx**4-1d0)/15.d0
      do i=1,nx-1
         dx2(i,i)=-((nx*nx-1d0)*(1d0-x(i)**2d0)+3.d0)/
     $        (3.d0*(1d0-x(i)**2d0)**2)
      enddo

      sign=-1d0
      do j=1,nx
         dx2(0,j)=sign*2d0*((2d0*nx*nx+1d0)*(1d0-x(j))-6.d0)/
     $        (3.d0*cx(j)*(1d0-x(j))**2)
         sign=-sign
      enddo

      sign=1d0
      if(mod(nx,2).ne.0) sign=-1d0
      do j=0,nx-1
         dx2(nx,j)=sign*2d0*((2d0*nx*nx+1d0)*(1d0+x(j))-6.d0)/
     $        (3.d0*cx(j)*(1d0+x(j))**2)
         sign=-sign
      enddo

      do j=0,nr
         do i=0,nr
            dr1(i,j,2)=dx1(i,j)+dx1(i,nx-j)
            dr2(i,j,2)=dx2(i,j)+dx2(i,nx-j)
            dr1(i,j,1)=dx1(i,j)-dx1(i,nx-j)
            dr2(i,j,1)=dx2(i,j)-dx2(i,nx-j)
            lap(i,j,2)=dr2(i,j,2)+ri(i)*dr1(i,j,2)
            lap(i,j,1)=dr2(i,j,1)+ri(i)*dr1(i,j,1)
         enddo
      enddo

***********************************************************************
*     Z                                                               *
***********************************************************************

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

      dz1(0,0)=(2d0*nz*nz+1d0)/6.d0
      dz1(nz,nz)=-(2d0*nz*nz+1d0)/6.d0
      do i=1,nz-1
         dz1(i,i)=-z(i)/(2.0d0*(1d0-z(i)**2d0))
      enddo

***************************

      sign=-1d0
      do j=0,nz
         do i=1,nz-1
            if(i.ne.j) then
              dz2(i,j)=sign*(z(i)*z(i)+z(i)*z(j)-2d0)/
     $             (cz(j)*(1d0-z(i)**2d0)*(z(i)-z(j))**2)
            endif
            sign=-sign
         enddo
         if(mod(nz,2).ne.0) sign=-sign
      enddo

      dz2(0,0)=(nz**4-1d0)/15.d0
      dz2(nz,nz)=(nz**4-1d0)/15.d0
      do i=1,nz-1
         dz2(i,i)=-((nz*nz-1d0)*(1d0-z(i)**2d0)+3.d0)/
     $        (3.d0*(1d0-z(i)**2d0)**2)
      enddo

      sign=-1d0
      do j=1,nz
         dz2(0,j)=sign*2d0*((2d0*nz*nz+1d0)*(1d0-z(j))-6.d0)/
     $        (3.d0*cz(j)*(1d0-z(j))**2)
         sign=-sign
      enddo

      sign=1d0
      if(mod(nz,2).ne.0) sign=-1d0
      do j=0,nz-1
         dz2(nz,j)=sign*2d0*((2d0*nz*nz+1d0)*(1d0+z(j))-6.d0)/
     $        (3.d0*cz(j)*(1d0+z(j))**2)
         sign=-sign
      enddo

***********************************************************************
*     Matrix M for Boundary Condition at the top                      *
***********************************************************************

      do j=0,nr
        do i=0,nr
          M(i,j,1) = dr2(i,j,1) + dr1(i,j,1)/r(i) - 1d0/r(i)**2d0
          M(i,j,2) = dr2(i,j,2) + dr1(i,j,2)/r(i) - 1d0/r(i)**2d0
        enddo
      enddo

      return
      end subroutine deriv

c=======================================================================
c Se hace la derivada radial de cada modo de Fourier, se indica la paridad
c a partir de los valores de k y l.
c Si in1=1 y in2=2, empieza por modo cero impar y sigue n=1 par.......
c Si in1=2 y in2=1, empieza por modo cero par y sigue n=1 impar.......

      subroutine derivr(f,f3,in1,in2,dr1,RAD,nr,nz,nn)
      implicit none
      integer nr,nz,nn,nr1,nz1
      REAL*8 dr1(0:nr,0:nr,2)
      REAL*8 RAD,alpha,f(0:nr,0:nz,0:nn-1),f3(0:nr,0:nz,0:nn-1)
      INTEGER in(2),n,j,k,nk,in1,in2

      nr1=nr+1 ; nz1=nz+1
      alpha=1d0/RAD
      in(1)=in1 ; in(2)=in2

c     Hacemos n=0
      call DGEMM('N','N',nr1,nz1,nr1,alpha,dr1(0,0,in(1)),nr1,
     $           f(0,0,0),nr1,0.d0,f3(0,0,0),nr1)
      j=in(1)
      in(1)=in(2)
      in(2)=j

c     Hacemos los otros n
!$OMP PARALLEL
!$OMP+SHARED(alpha,dr1,in,f,f3,nn,nr1,nz1)
!$OMP+PRIVATE(n,k,nk,j)
!$OMP DO
      do n=1,nn-3,2
         do k=0,1
            nk=n+k
            call DGEMM('N','N',nr1,nz1,nr1,alpha,dr1(0,0,in(1)),nr1,
     $           f(0,0,nk),nr1,0.d0,f3(0,0,nk),nr1)
         enddo
         j=in(1)
         in(1)=in(2)
         in(2)=j
      enddo
!$OMP END DO
!$OMP END PARALLEL
c     Hacemos el nn/2 -1 (solo hay uno)
      nk=nn-1
      call DGEMM('N','N',nr1,nz1,nr1,alpha,dr1(0,0,in(1)),nr1,
     $     f(0,0,nk),nr1,0.d0,f3(0,0,nk),nr1)

      return
      end subroutine derivr

c=======================================================================

      subroutine mdiag(work,dz2,lap,riq,axd,vxd,bxd,azd,vzd,bzd,axn,
     &     vxn,bxn,azn,vzn,bzn,cx03,cz03,cz13,nr,nz,nn,nmax)
      implicit none
      integer nr,nz,nn,nmax
      real*8 dz2(0:nz,0:nz),lap(0:nr,0:nr,2),riq(0:nr)
      real*8 dxd(1:nr,1:nr),dxn(1:nr,1:nr),aux1(1:nr)
      real*8 dzd(1:nz-1,1:nz-1),dzn(1:nz-1,1:nz-1),aux2(1:nz-1)
      real*8 axd(1:nr,1:nr,0:nn/2+1),vxd(1:nr,0:nn/2+1),
     $     bxd(1:nr,1:nr,0:nn/2+1),axn(1:nr,1:nr,0:nn/2+1),
     $     vxn(1:nr,0:nn/2+1),bxn(1:nr,1:nr,0:nn/2+1),
     $     azd(1:nz-1,1:nz-1),vzd(1:nz-1),bzd(1:nz-1,1:nz-1),
     $     azn(1:nz-1,1:nz-1),vzn(1:nz-1),bzn(1:nz-1,1:nz-1)
      real*8 cx03(1:nr,2),cz03(1:nz-1),cz13(1:nz-1)
      integer index(nmax), LWORK
      real*8 WORK(1:(nmax+2)*6)
      integer i,j,n,nnq
      integer INFO,IPIV(nmax)

      LWORK=(nmax+2)*6

      do n=0,nn/2+1
         nnq=n*n

****************************************************************************
*     axd  vxd  bxd  r(dirichlet)                                          *
****************************************************************************
         if(mod(n,2).ne.0) then
            do i=1,nr
               do j=1,nr
                  dxd(i,j)=lap(i,j,1)
               enddo
               dxd(i,i)=dxd(i,i)-dfloat(nnq)*riq(i)
            enddo
         else
            do i=1,nr
               do j=1,nr
                  dxd(i,j)=lap(i,j,2)
               enddo
               dxd(i,i)=dxd(i,i)-dfloat(nnq)*riq(i)
            enddo
         endif
         call DGEEV('N','V',nr,dxd,nr,vxd(1,n),aux1,bxd(1,1,n),nr,
     $        axd(1,1,n),nr,WORK,LWORK,INFO)
         do i=1,nr
            index(i)=i
         enddo
         call SHELL(vxd(1,n),nr,index,nmax)
         do i=1,nr
            do j=1,nr
               dxd(i,j)=axd(i,index(nr+1-j),n)
            enddo
         enddo
         do i=1,nr
            do j=1,nr
               axd(i,j,n)=dxd(i,j)
            enddo
         enddo
         do i=1,nr
            aux1(i)=vxd(nr+1-i,n)
         enddo
         do i=1,nr
            vxd(i,n)=aux1(i)
         enddo
         do i=1,nr
            do j=1,nr
               bxd(i,j,n)=0
            enddo
         enddo
         do i=1,nr
            bxd(i,i,n)=1d0
         enddo
         do i=1,nr
            do j=1,nr
               dxd(i,j)=axd(i,j,n)
            enddo
         enddo
         call DGESV(nr,nr,dxd,nr,IPIV,bxd(1,1,n),nr,INFO)

*******************************************************************************
*     axn  vxn  bxn  r(newmann)                                               *
*******************************************************************************
         if(mod(n,2).ne.0) then
            do i=1,nr
               do j=1,nr
                  dxn(i,j)=lap(i,j,1)+lap(i,0,1)*cx03(j,1)
               enddo
               dxn(i,i)=dxn(i,i)-dfloat(nnq)*riq(i)
            enddo
         else
            do i=1,nr
               do j=1,nr
                  dxn(i,j)=lap(i,j,2)+lap(i,0,2)*cx03(j,2)
               enddo
               dxn(i,i)=dxn(i,i)-dfloat(nnq)*riq(i)
            enddo
         endif
         call DGEEV('N','V',nr,dxn,nr,vxn(1,n),aux1,bxn(1,1,n),nr,
     $        axn(1,1,n),nr,WORK,LWORK,INFO)
         do i=1,nr
            index(i)=i
         enddo
         call SHELL(vxn(1,n),nr,index,nmax)
         do i=1,nr
            do j=1,nr
               dxn(i,j)=axn(i,index(nr+1-j),n)
            enddo
         enddo
         do i=1,nr
            do j=1,nr
               axn(i,j,n)=dxn(i,j)
            enddo
         enddo
         do i=1,nr
            aux1(i)=vxn(nr+1-i,n)
         enddo
         do i=1,nr
            vxn(i,n)=aux1(i)
         enddo
         do i=1,nr
            do j=1,nr
               bxn(i,j,n)=0
            enddo
         enddo
         do i=1,nr
            bxn(i,i,n)=1d0
         enddo
         do i=1,nr
            do j=1,nr
               dxn(i,j)=axn(i,j,n)
            enddo
         enddo
         call DGESV(nr,nr,dxn,nr,IPIV,bxn(1,1,n),nr,INFO)
      enddo

***********************************************************************
*     azd  vzd  bzd     z (dirichlet)
***********************************************************************
      do i=1,nz-1
         do j=1,nz-1
            dzd(i,j)=dz2(i,j)
         enddo
      enddo
      call DGEEV('N','V',nz-1,dzd,nz-1,vzd,aux2,bzd,nz-1,azd,nz-1,
     $     WORK,LWORK,INFO)
      do i=1,nz-1
         index(i)=i
      enddo
      call SHELL(vzd,nz-1,index,nmax)
      do i=1,nz-1
         do j=1,nz-1
            dzd(i,j)=azd(i,index(nz-j))
         enddo
      enddo
      do i=1,nz-1
         do j=1,nz-1
            azd(i,j)=dzd(i,j)
         enddo
      enddo
      do i=1,nz-1
         aux2(i)=vzd(nz-i)
      enddo
      do i=1,nz-1
         vzd(i)=aux2(i)
      enddo
      do i=1,nz-1
         do j=1,nz-1
            bzd(i,j)=0
         enddo
      enddo
      do i=1,nz-1
         bzd(i,i)=1d0
      enddo
      do i=1,nz-1
         do j=1,nz-1
            dzd(i,j)=azd(i,j)
         enddo
      enddo
      call DGESV(nz-1,nz-1,dzd,nz-1,IPIV,bzd,nz-1,INFO)

**********************************************************************
*     azn  vzn  bzn            z (newmann)
***********************************************************************
      do i=1,nz-1
         do j=1,nz-1
            dzn(i,j)=dz2(i,j)+dz2(i,0)*cz03(j)+dz2(i,nz)*cz13(j)
         enddo
      enddo
      call DGEEV('N','V',nz-1,dzn,nz-1,vzn,aux2,bzn,nz-1,azn,nz-1,
     $     WORK,LWORK,INFO)
      do i=1,nz-1
         index(i)=i
      enddo
      call SHELL(vzn,nz-1,index,nmax)
      do i=1,nz-1
         do j=1,nz-1
            dzn(i,j)=azn(i,index(nz-j))
         enddo
      enddo
      do i=1,nz-1
         do j=1,nz-1
            azn(i,j)=dzn(i,j)
         enddo
      enddo
      do i=1,nz-1
         aux2(i)=vzn(nz-i)
      enddo
      do i=1,nz-1
         vzn(i)=aux2(i)
      enddo
      do i=1,nz-1
         do j=1,nz-1
            bzn(i,j)=0
         enddo
      enddo
      do i=1,nz-1
         bzn(i,i)=1d0
      enddo
      do i=1,nz-1
         do j=1,nz-1
            dzn(i,j)=azn(i,j)
         enddo
      enddo
      call DGESV(nz-1,nz-1,dzn,nz-1,IPIV,bzn,nz-1,INFO)
      return
      end subroutine mdiag

c=======================================================================

      subroutine shell(arr,n,index,nmax)
      implicit none
      integer nmax
      real*8 arr(nmax),t,aln2i,tiny
      integer index(nmax),i,j,k,l,m,n,lognb2,nnn
      aln2i=1d0/log(2d0)
      tiny=1.d-5
      lognb2=int(alog(float(n))*aln2i+tiny)
      m=n
      do nnn=1,lognb2
         m=m/2
         k=n-m
         do j=1,k
            i=j
 3          continue
            l=i+m
            if(arr(l).lt.arr(i)) then
               t=arr(i)
               k=index(i)
               arr(i)=arr(l)
               index(i)=index(l)
               arr(l)=t
               index(l)=k
               i=i-m
               if(i.ge.1) goto 3
            endif
         enddo
      enddo
      return
      end subroutine shell

c=====================================================================
c     Computes the direct and inverse Fourier transforms for real
c     (u,w,T) and imaginary (v) variables

      subroutine transform(tdr,tir,tdi,tii,nn)
      implicit none
      integer nn,j,n,nd,nd1
      real*8  tdr(0:nn-1,0:nn-1),tdi(0:nn-1,0:nn-1)
      real*8  tir(0:nn-1,0:nn-1),tii(0:nn-1,0:nn-1)
      real*8  con,px,pnj,pi

      pi=dacos(-1d0)
      con=1d0/(dfloat(nn))
      px=2d0*pi*con
      do j=0,nn-1
         tdr(j,0)=1d0
         tir(0,j)=con
         tdi(j,0)=1d0
         tii(0,j)=con
      enddo
      do n=1,nn/2-1
         nd=2*n
         nd1=2*n-1
         do j=0,nn-1
            pnj=px*dfloat(n)*dfloat(j)
            tdr(j,nd1)=2d0*dcos(pnj)
            tdr(j,nd)=-2d0*dsin(pnj)
            tir(nd1,j)=con*dcos(pnj)
            tir(nd,j)=-con*dsin(pnj)
            tdi(j,nd1)=2d0*dsin(pnj)
            tdi(j,nd)=2d0*dcos(pnj)
            tii(nd1,j)=con*dsin(pnj)
            tii(nd,j)=con*dcos(pnj)
         enddo
      enddo
      do j=0,nn-1
         pnj=px*dfloat(nn/2)*dfloat(j)
         tdr(j,nn-1)=dcos(pnj)
         tir(nn-1,j)=con*dcos(pnj)
         tdi(j,nn-1)=dcos(pnj)
         tii(nn-1,j)=con*dcos(pnj)
      enddo
      return
      end subroutine transform


c=======================================================================
c     Poisson solver for the pressure (Neumann boundary conditions on
c     end plates and wall)
      subroutine poisp(a,f,f3,dz2,lap,axn,vxn,bxn,azn,vzn,bzn,RAD,ALT,
     &     cx01,cx03,cz01,cz02,cz03,cz11,cz12,cz13,nr,nz,nn)
      implicit none
      integer nr,nz,nn
      real*8 dz2(0:nz,0:nz),lap(0:nr,0:nr,2)
      real*8 axn(1:nr,1:nr,0:nn/2+1),vxn(1:nr,0:nn/2+1),
     $       bxn(1:nr,1:nr,0:nn/2+1),
     $       azn(1:nz-1,1:nz-1),vzn(1:nz-1),bzn(1:nz-1,1:nz-1)
      real*8 cx01(2),cx03(1:nr,2),cz01,cz02,cz03(1:nz-1),
     $       cz11,cz12,cz13(1:nz-1)
      real*8 RAD,ALT
      real*8 a,f(0:nr,0:nz,0:nn-1),f3(0:nr,0:nz,0:nn-1)
      real*8 t,al2,radq,ar,nmq
      real*8 h(1:nr,1:nz-1),w(1:nr,1:nz-1),aux(1:nr,1:nz-1)
      integer in(2),n,nm,i,j,k,l,s,nk

      ar=RAD/ALT
      radq=RAD*RAD
      al2=ar*ar

c     Rescale Chebyshev radial derivative on the boundary to be used in
c     an interval of length RAD
      do n=0,nn-1
         do j=1,nz-1
            f(0,j,n)=f(0,j,n)*RAD
         enddo
      enddo
c     Rescale Chebyshev axial derivative on the boundary to be used in
c     an interval of length ALT/2d0
      do n=0,nn-1
         do i=0,nr
            f(i,0,n)=f(i,0,n)*ALT/2d0
            f(i,nz,n)=f(i,nz,n)*ALT/2d0
         enddo
      enddo

c=======================================================================
c     n=0 mode: even in radial direction
      in(1)=2
      in(2)=1
      do i=1,nr
         do j=1,nz-1
            h(i,j)=f(i,j,0)*radq-
     $           lap(i,0,in(1))*cx01(in(1))*f(0,j,0)-
     $           4.d0*al2*(dz2(j,0)*(cz01*f(i,0,0)+cz02*f(i,nz,0))+
     $           dz2(j,nz)*(cz11*f(i,0,0)+cz12*f(i,nz,0)))
         enddo
      enddo
      call DGEMM('N','T',nr,nz-1,nz-1,1d0,h(1,1),nr,bzn(1,1),nz-1,0.d0
     $     ,aux(1,1),nr)
      call DGEMM('N','N',nr,nz-1,nr,1d0,bxn(1,1,0),nr,aux(1,1),nr,0.d0
     $     ,w(1,1),nr)
      do i=2,nr
         do s=1,nz-1
            t=1d0/(vxn(i,0)+4.d0*al2*vzn(s))
            w(i,s)=t*w(i,s)
         enddo
      enddo
      do s=2,nz-1
         t=1d0/(vxn(1,0)+4.d0*al2*vzn(s))
         w(1,s)=t*w(1,s)
      enddo
c     Compute w(1,1). Set value in (1,1)
      t=0.d0
      do i=2,nr
         do s=1,nz-1
            t=t+axn(1,i,0)*azn(1,s)*w(i,s)
         enddo
      enddo
      do s=2,nz-1
         t=t+axn(1,1,0)*azn(1,s)*w(1,s)
      enddo
      w(1,1)=(a-t)/(axn(1,1,0)*azn(1,1))

      call DGEMM('N','T',nr,nz-1,nz-1,1d0,w(1,1),nr,azn(1,1),nz-1,0.d0
     $     ,aux(1,1),nr)
      call DGEMM('N','N',nr,nz-1,nr,1d0,axn(1,1,0),nr,aux(1,1),nr,0.d0
     $     ,f3(1,1,0),nr+1)

c      Compute f on the boundaries

c     Wall
      do j=1,nz-1
         f3(0,j,0)=cx01(in(1))*f(0,j,0)
         do k=1,nr
            f3(0,j,0)=f3(0,j,0)+cx03(k,in(1))*f3(k,j,0)
         enddo
      enddo

c     Plates
      do i=0,nr
         f3(i,0,0)=cz01*f(i,0,0)+cz02*f(i,nz,0)
         f3(i,nz,0)=cz11*f(i,0,0)+cz12*f(i,nz,0)
         do l=1,nz-1
            f3(i,0,0)=f3(i,0,0)+cz03(l)*f3(i,l,0)
            f3(i,nz,0)=f3(i,nz,0)+cz13(l)*f3(i,l,0)
         enddo
      enddo
      j=in(1)
      in(1)=in(2)
      in(2)=j

c=======================================================================
c     n!=0 modes
      do n=1,nn-3,2
         nm=(n+1)/2
         nmq=dfloat(nm*nm)
         do k=0,1
            nk=n+k
            do i=1,nr
               do j=1,nz-1
                  h(i,j)=f(i,j,nk)*radq-
     $                 lap(i,0,in(1))*cx01(in(1))*f(0,j,nk) - 4.d0*al2*
     $                 (dz2(j,0)*(cz01*f(i,0,nk)+cz02*f(i,nz,nk))+
     $                 dz2(j,nz)*(cz11*f(i,0,nk)+cz12*f(i,nz,nk)))
               enddo
            enddo

            call DGEMM('N','T',nr,nz-1,nz-1,1d0,h(1,1),nr,
     $           bzn(1,1),nz-1,0.d0,aux(1,1),nr)
            call DGEMM('N','N',nr,nz-1,nr,1d0,bxn(1,1,nm),nr,
     $           aux(1,1),nr,0.d0,w(1,1),nr)

            do i=1,nr
               do s=1,nz-1
                  t=1d0/(vxn(i,nm)+4.d0*al2*vzn(s))
                  w(i,s)=t*w(i,s)
               enddo
            enddo

            call DGEMM('N','T',nr,nz-1,nz-1,1d0,w(1,1),nr,
     $           azn(1,1),nz-1,0.d0,aux(1,1),nr)
            call DGEMM('N','N',nr,nz-1,nr,1d0,axn(1,1,nm),nr,
     $           aux(1,1),nr,0.d0,f3(1,1,nk),nr+1)

c     Compute f on the boundaries

c     Wall
            do j=1,nz-1
               f3(0,j,nk)=cx01(in(1))*f(0,j,nk)
               do i=1,nr
                  f3(0,j,nk)=f3(0,j,nk)+cx03(i,in(1))*f3(i,j,nk)
               enddo
            enddo

c     Plates
            do i=0,nr
               f3(i,0,nk)=cz01*f(i,0,nk)+cz02*f(i,nz,nk)
               f3(i,nz,nk)=cz11*f(i,0,nk)+cz12*f(i,nz,nk)
               do l=1,nz-1
                  f3(i,0,nk)=f3(i,0,nk)+cz03(l)*f3(i,l,nk)
                  f3(i,nz,nk)=f3(i,nz,nk)+cz13(l)*f3(i,l,nk)
               enddo
            enddo

c     end of loop in k
         enddo
         j=in(1)
         in(1)=in(2)
         in(2)=j

c     end of loop in n
      enddo

c=======================================================================
c     n/2-1 mode
      nm=nn/2
      nmq=dfloat(nm*nm)
      nk=nn-1
      do i=1,nr
         do j=1,nz-1
            h(i,j)=f(i,j,nk)*radq-
     $           lap(i,0,in(1))*cx01(in(1))*f(0,j,nk) -
     $           4.d0*al2*(dz2(j,0)*(cz01*f(i,0,nk)+cz02*f(i,nz,nk))+
     $           dz2(j,nz)*(cz11*f(i,0,nk)+cz12*f(i,nz,nk)))
         enddo
      enddo

      call DGEMM('N','T',nr,nz-1,nz-1,1d0,h(1,1),nr,
     $            bzn(1,1),nz-1,0.d0,aux(1,1),nr)
      call DGEMM('N','N',nr,nz-1,nr,1d0,bxn(1,1,nm),nr,
     $            aux(1,1),nr,0.d0,w(1,1),nr)
      do i=1,nr
         do s=1,nz-1
            t=1d0/(vxn(i,nm)+4.d0*al2*vzn(s))
            w(i,s)=t*w(i,s)
         enddo
      enddo
      call DGEMM('N','T',nr,nz-1,nz-1,1d0,w(1,1),nr,
     $     azn(1,1),nz-1,0.d0,aux(1,1),nr)
      call DGEMM('N','N',nr,nz-1,nr,1d0,axn(1,1,nm),nr,
     $     aux(1,1),nr,0.d0,f3(1,1,nk),nr+1)

c     Compute of f on the boundaries

c     Wall
      do j=1,nz-1
         f3(0,j,nk)=cx01(in(1))*f(0,j,nk)
         do k=1,nr
            f3(0,j,nk)=f3(0,j,nk)+cx03(k,in(1))*f3(k,j,nk)
         enddo
      enddo
c     Plates
      do i=0,nr
         f3(i,0,nk)=cz01*f(i,0,nk)+cz02*f(i,nz,nk)
         f3(i,nz,nk)=cz11*f(i,0,nk)+cz12*f(i,nz,nk)
         do k=1,nz-1
            f3(i,0,nk)=f3(i,0,nk)+cz03(k)*f3(i,k,nk)
            f3(i,nz,nk)=f3(i,nz,nk)+cz13(k)*f3(i,k,nk)
         enddo
      enddo

      RETURN
      END subroutine poisp

c=======================================================================
c     Es el poisson de la temperatura para el cilindro vertical.
c     Condiciones de funcion conocida en las tapas i derivada normal en
c     el lateral

      SUBROUTINE POISTA(a,f,f3,dz2,lap,axn,vxn,bxn,azd,vzd,bzd,RAD,ALT,
     &     cx01,cx03,nr,nz,nn)
      implicit none
      integer nr,nz,nn
      real*8 dz2(0:nz,0:nz),lap(0:nr,0:nr,2)
      real*8 axn(1:nr,1:nr,0:nn/2+1),vzd(1:nz-1),bzd(1:nz-1,1:nz-1),
     &     vxn(1:nr,0:nn/2+1),bxn(1:nr,1:nr,0:nn/2+1),azd(1:nz-1,1:nz-1)
      real*8 cx01(2),cx03(1:nr,2),RAD,ALT,t,al2,radq,ar,nmq
      real*8 a,f(0:nr,0:nz,0:nn-1),f3(0:nr,0:nz,0:nn-1)
      real*8 h(1:nr,1:nz-1),w(1:nr,1:nz-1),aux(1:nr,1:nz-1)
      integer in(2),n,nm,i,j,k,s,nk
      ar=RAD/ALT
      radq=RAD*RAD
      al2=ar*ar

c     multiplicamos la derivada radial real en el contorno por RAD para que
c     sea la derivada chebyshev
      do n=0,nn-1
         do j=1,nz-1
            f(0,j,n)=f(0,j,n)*RAD
         enddo
      enddo
c     Hacemos el n=0 que es par en r
      in(1)=2
      in(2)=1
      do i=1,nr
         do j=1,nz-1
            h(i,j)=f(i,j,0)*radq-
     &           lap(i,0,in(1))*cx01(in(1))*f(0,j,0)-
     &           4.d0*al2*(dz2(j,0)*f(i,0,0)+dz2(j,nz)*f(i,nz,0))
         enddo
      enddo
      call DGEMM('N','T',nr,nz-1,nz-1,1d0,h(1,1),nr,
     &     bzd(1,1),nz-1,0.d0,aux(1,1),nr)
      call DGEMM('N','N',nr,nz-1,nr,1d0,bxn(1,1,0),nr,
     &     aux(1,1),nr,0.d0,w(1,1),nr)
      do i=1,nr
         do s=1,nz-1
            t=1d0/(vxn(i,0)+4.d0*al2*vzd(s)-a*radq)
            w(i,s)=t*w(i,s)
         enddo
      enddo

      call DGEMM('N','T',nr,nz-1,nz-1,1d0,w(1,1),nr,
     &     azd(1,1),nz-1,0.d0,aux(1,1),nr)
      call DGEMM('N','N',nr,nz-1,nr,1d0,axn(1,1,0),nr,
     &     aux(1,1),nr,0.d0,f3(1,1,0),nr+1)

C     CONTORN DE LA FUNCIO F

C     TAPES
      do i=0,nr
         f3(i,0,0)=f(i,0,0)
         f3(i,nz,0)=f(i,nz,0)
      enddo

C     LATERALS
      do j=1,nz-1
         f3(0,j,0)=cx01(in(1))*f(0,j,0)
         do k=1,nr
            f3(0,j,0)=f3(0,j,0)+cx03(k,in(1))*f3(k,j,0)
         enddo
      enddo
      j=in(1)
      in(1)=in(2)
      in(2)=j

c     Hacemos los otros n
      do n=1,nn-3,2
         nm=(n+1)/2
         nmq=dfloat(nm*nm)
         do k=0,1
            nk=n+k
            do i=1,nr
               do j=1,nz-1
                  h(i,j)=f(i,j,nk)*radq-
     &                 lap(i,0,in(1))*cx01(in(1))*f(0,j,nk) - 4.d0*al2*
     &                 (dz2(j,0)*f(i,0,nk)+dz2(j,nz)*f(i,nz,nk))
               enddo
            enddo

            call DGEMM('N','T',nr,nz-1,nz-1,1d0,h(1,1),nr,
     &           bzd(1,1),nz-1,0.d0,aux(1,1),nr)
            call DGEMM('N','N',nr,nz-1,nr,1d0,bxn(1,1,nm),nr,
     &           aux(1,1),nr,0.d0,w(1,1),nr)

            do i=1,nr
               do s=1,nz-1
                  t=1d0/(vxn(i,nm)+4.d0*al2*vzd(s)-a*radq)
                  w(i,s)=t*w(i,s)
               enddo
            enddo

            call DGEMM('N','T',nr,nz-1,nz-1,1d0,w(1,1),nr,
     &           azd(1,1),nz-1,0.d0,aux(1,1),nr)
            call DGEMM('N','N',nr,nz-1,nr,1d0,axn(1,1,nm),nr,
     &           aux(1,1),nr,0.d0,f3(1,1,nk),nr+1)

C     CONTORN DE LA FUNCIO F

C     TAPES
            do i=0,nr
               f3(i,0,nk)=f(i,0,nk)
               f3(i,nz,nk)=f(i,nz,nk)
            enddo

C     LATERALS
            do j=1,nz-1
               f3(0,j,nk)=cx01(in(1))*f(0,j,nk)
               do i=1,nr
                  f3(0,j,nk)=f3(0,j,nk)+cx03(i,in(1))*f3(i,j,nk)
               enddo
            enddo

c     cierro el do de k             *
         enddo
         j=in(1)
         in(1)=in(2)
         in(2)=j

c     cierro el do de n             *
      enddo


c     Hacemos el nn/2 -1 (solo hay uno)
      nm=nn/2
      nmq=dfloat(nm*nm)
      nk=nn-1
      do i=1,nr
         do j=1,nz-1
            h(i,j)=f(i,j,nk)*radq-
     &           lap(i,0,in(1))*cx01(in(1))*f(0,j,nk) -
     &           4.d0*al2*(dz2(j,0)*f(i,0,nk)+dz2(j,nz)*f(i,nz,nk))
         enddo
      enddo

      call DGEMM('N','T',nr,nz-1,nz-1,1d0,h(1,1),nr,
     &     bzd(1,1),nz-1,0.d0,aux(1,1),nr)
      call DGEMM('N','N',nr,nz-1,nr,1d0,bxn(1,1,nm),nr,
     &     aux(1,1),nr,0.d0,w(1,1),nr)

      do i=1,nr
         do s=1,nz-1
            t=1d0/(vxn(i,nm)+4.d0*al2*vzd(s)-a*radq)
            w(i,s)=t*w(i,s)
         enddo
      enddo
      call DGEMM('N','T',nr,nz-1,nz-1,1d0,w(1,1),nr,
     &     azd(1,1),nz-1,0.d0,aux(1,1),nr)
      call DGEMM('N','N',nr,nz-1,nr,1d0,axn(1,1,nm),nr,
     &     aux(1,1),nr,0.d0,f3(1,1,nk),nr+1)

C     CONTORN DE LA FUNCIO F

C     TAPES
      do i=0,nr
         f3(i,0,nk)=f(i,0,nk)
         f3(i,nz,nk)=f(i,nz,nk)
      enddo

C     LATERALS
      do j=1,nz-1
         f3(0,j,nk)=cx01(in(1))*f(0,j,nk)
         do k=1,nr
            f3(0,j,nk)=f3(0,j,nk)+cx03(k,in(1))*f3(k,j,nk)
         enddo
      enddo

      return
      end subroutine poista

c=======================================================================

c     Es el poisson de la temperatura para el cilindro vertical, Condiciones
c     de funcion conocida en las tapas y en el lateral

      SUBROUTINE POIST(a,f,f3,dz2,lap,axd,vxd,bxd,azd,vzd,bzd,RAD,ALT,
     &                 nr,nz,nn)
      implicit none
      integer nr,nz,nn
      REAL*8 dz2(0:nz,0:nz),lap(0:nr,0:nr,2)
      REAL*8 axd(1:nr,1:nr,0:nn/2+1),vxd(1:nr,0:nn/2+1),
     &     azd(1:nz-1,1:nz-1),bxd(1:nr,1:nr,0:nn/2+1),
     &     vzd(1:nz-1),bzd(1:nz-1,1:nz-1)
      REAL*8 RAD,ALT,t,al2,radq,ar
      REAL*8 a,f(0:nr,0:nz,0:nn-1),f3(0:nr,0:nz,0:nn-1)
      REAL*8 h(1:nr,1:nz-1),w(1:nr,1:nz-1),aux(1:nr,1:nz-1)
      INTEGER n,nm,i,j,k,s,nk
      ar=RAD/ALT
      radq=RAD*RAD
      al2=ar*ar

c     Hacemos el n=0 que es par en r

      do i=1,nr
         do j=1,nz-1
            h(i,j)=f(i,j,0)*radq-
     $           lap(i,0,2)*f(0,j,0)-
     $           4.d0*al2*(dz2(j,0)*f(i,0,0)+dz2(j,nz)*f(i,nz,0))
         enddo
      enddo

      call DGEMM('N','T',nr,nz-1,nz-1,1d0,h(1,1),nr,
     $     bzd(1,1),nz-1,0.d0,aux(1,1),nr)
      call DGEMM('N','N',nr,nz-1,nr,1d0,bxd(1,1,0),nr,
     $     aux(1,1),nr,0.d0,w(1,1),nr)

      do i=1,nr
         do s=1,nz-1
            t=1d0/(vxd(i,0)+4.d0*al2*vzd(s)-a*radq)
            w(i,s)=t*w(i,s)
         enddo
      enddo

      call DGEMM('N','T',nr,nz-1,nz-1,1d0,w(1,1),nr,
     $     azd(1,1),nz-1,0.d0,aux(1,1),nr)
      call DGEMM('N','N',nr,nz-1,nr,1d0,axd(1,1,0),nr,
     $     aux(1,1),nr,0.d0,f3(1,1,0),nr+1)

C     CONTORN DE LA FUNCIO F

C     TAPES
      do i=0,nr
         f3(i,0,0)=f(i,0,0)
         f3(i,nz,0)=f(i,nz,0)
      enddo

C     LATERALS
      do j=1,nz-1
         f3(0,j,0)=f(0,j,0)
      enddo

c     Hacemos los otros n
      do n=1,nn-3,2
         nm=(n+1)/2
         do k=0,1
            nk=n+k
            do i=1,nr
               do j=1,nz-1
                  h(i,j)=f(i,j,nk)*radq-
     $                 lap(i,0,1)*f(0,j,nk) - 4.d0*al2*
     $                 (dz2(j,0)*f(i,0,nk)+dz2(j,nz)*f(i,nz,nk))
               enddo
            enddo

            call DGEMM('N','T',nr,nz-1,nz-1,1d0,h(1,1),nr,
     $           bzd(1,1),nz-1,0.d0,aux(1,1),nr)
            call DGEMM('N','N',nr,nz-1,nr,1d0,bxd(1,1,nm),nr,
     $           aux(1,1),nr,0.d0,w(1,1),nr)

            do i=1,nr
               do s=1,nz-1
                  t=1d0/(vxd(i,nm)+4.d0*al2*vzd(s)-a*radq)
                  w(i,s)=t*w(i,s)
               enddo
            enddo

            call DGEMM('N','T',nr,nz-1,nz-1,1d0,w(1,1),nr,
     $           azd(1,1),nz-1,0.d0,aux(1,1),nr)
            call DGEMM('N','N',nr,nz-1,nr,1d0,axd(1,1,nm),nr,
     $           aux(1,1),nr,0.d0,f3(1,1,nk),nr+1)

C     CONTORN DE LA FUNCIO F

C     TAPES
            do i=0,nr
               f3(i,0,nk)=f(i,0,nk)
               f3(i,nz,nk)=f(i,nz,nk)
            enddo

C     LATERALS
            do j=1,nz-1
               f3(0,j,nk)=f(0,j,nk)
            enddo

c     cierro el do de k             *
         enddo

c     cierro el do de n             *
      enddo

c     Hacemos el nn/2 -1 (solo hay uno)
      nm=nn/2
      nk=nn-1
      do i=1,nr
         do j=1,nz-1
            h(i,j)=f(i,j,nk)*radq-
     $           lap(i,0,2)*f(0,j,nk) -
     $           4.d0*al2*(dz2(j,0)*f(i,0,nk)+dz2(j,nz)*f(i,nz,nk))
         enddo
      enddo

      call DGEMM('N','T',nr,nz-1,nz-1,1d0,h(1,1),nr,
     $     bzd(1,1),nz-1,0.d0,aux(1,1),nr)
      call DGEMM('N','N',nr,nz-1,nr,1d0,bxd(1,1,nm),nr,
     $     aux(1,1),nr,0.d0,w(1,1),nr)

      do i=1,nr
         do s=1,nz-1
            t=1d0/(vxd(i,nm)+4.d0*al2*vzd(s)-a*radq)
            w(i,s)=t*w(i,s)
         enddo
      enddo

      call DGEMM('N','T',nr,nz-1,nz-1,1d0,w(1,1),nr,
     $     azd(1,1),nz-1,0.d0,aux(1,1),nr)
      call DGEMM('N','N',nr,nz-1,nr,1d0,axd(1,1,nm),nr,
     $     aux(1,1),nr,0.d0,f3(1,1,nk),nr+1)

C     CONTORN DE LA FUNCIO F

C     TAPES
      do i=0,nr
         f3(i,0,nk)=f(i,0,nk)
         f3(i,nz,nk)=f(i,nz,nk)
      enddo

C     LATERALS
      do j=1,nz-1
         f3(0,j,nk)=f(0,j,nk)
      enddo

      RETURN
      END subroutine poist
c     Es el poisson de ur-iutheta para el cilindro vertical,
c     Equivale a hacer ur + utheta (reales e imaginarias). Condiciones
c     de funcion conocida en las tapas y en el lateral

c=======================================================================

      subroutine POISUL(a,f,f3,dz2,lap,axd,vxd,bxd,azd,vzd,bzd,RAD,ALT,
     &                  nr,nz,nn)
      IMPLICIT NONE
      INTEGER nr,nz,nn
      REAL*8 dz2(0:nz,0:nz),lap(0:nr,0:nr,2)
      REAL*8 axd(1:nr,1:nr,0:nn/2+1),vxd(1:nr,0:nn/2+1),
     &     azd(1:nz-1,1:nz-1),bxd(1:nr,1:nr,0:nn/2+1),
     &     vzd(1:nz-1),bzd(1:nz-1,1:nz-1)
      REAL*8 a,f(0:nr,0:nz,0:nn-1),f3(0:nr,0:nz,0:nn-1)
      REAL*8 RAD,ALT,t,al2,radq,ar,nmq
      REAL*8 h(1:nr,1:nz-1),w(1:nr,1:nz-1),aux(1:nr,1:nz-1)
      INTEGER in(2),n,nm,i,j,k,s,nk
      ar=RAD/ALT
      radq=RAD*RAD
      al2=ar*ar

c     Hacemos el n=0 que es impar en r
      in(1)=1
      in(2)=2
      do i=1,nr
         do j=1,nz-1
            h(i,j)=f(i,j,0)*radq-
     $           lap(i,0,in(1))*f(0,j,0)-
     $           4.d0*al2*(dz2(j,0)*f(i,0,0)+dz2(j,nz)*f(i,nz,0))
         enddo
      enddo

      call DGEMM('N','T',nr,nz-1,nz-1,1d0,h(1,1),nr,
     $     bzd(1,1),nz-1,0.d0,aux(1,1),nr)
      call DGEMM('N','N',nr,nz-1,nr,1d0,bxd(1,1,1),nr,
     $     aux(1,1),nr,0.d0,w(1,1),nr)

      do i=1,nr
         do s=1,nz-1
            t=1d0/(vxd(i,1)+4.d0*al2*vzd(s)-a*radq)
            w(i,s)=t*w(i,s)
         enddo
      enddo

      call DGEMM('N','T',nr,nz-1,nz-1,1d0,w(1,1),nr,
     $     azd(1,1),nz-1,0.d0,aux(1,1),nr)
      call DGEMM('N','N',nr,nz-1,nr,1d0,axd(1,1,1),nr,
     $     aux(1,1),nr,0.d0,f3(1,1,0),nr+1)

C     CONTORN DE LA FUNCIO F

C     TAPES
      do i=0,nr
         f3(i,0,0)=f(i,0,0)
         f3(i,nz,0)=f(i,nz,0)
      enddo

C     LATERALS
      do j=1,nz-1
         f3(0,j,0)=f(0,j,0)
      enddo
      j=in(1)
      in(1)=in(2)
      in(2)=j

c     Hacemos los otros n
      do n=1,nn-3,2
         nm=((n+1)/2)-1
         nmq=dfloat(nm*nm)
         do k=0,1
            nk=n+k
            do i=1,nr
               do j=1,nz-1
                  h(i,j)=f(i,j,nk)*radq-
     $                 lap(i,0,in(1))*f(0,j,nk) - 4.d0*al2*
     $                 (dz2(j,0)*f(i,0,nk)+dz2(j,nz)*f(i,nz,nk))
               enddo
            enddo

            call DGEMM('N','T',nr,nz-1,nz-1,1d0,h(1,1),nr,
     $           bzd(1,1),nz-1,0.d0,aux(1,1),nr)
            call DGEMM('N','N',nr,nz-1,nr,1d0,bxd(1,1,nm),nr,
     $           aux(1,1),nr,0.d0,w(1,1),nr)

            do i=1,nr
               do s=1,nz-1
                  t=1d0/(vxd(i,nm)+4.d0*al2*vzd(s)-a*radq)
                  w(i,s)=t*w(i,s)
               enddo
            enddo

            call DGEMM('N','T',nr,nz-1,nz-1,1d0,w(1,1),nr,
     $           azd(1,1),nz-1,0.d0,aux(1,1),nr)
            call DGEMM('N','N',nr,nz-1,nr,1d0,axd(1,1,nm),nr,
     $           aux(1,1),nr,0.d0,f3(1,1,nk),nr+1)

C     CONTORN DE LA FUNCIO F

C     TAPES
            do i=0,nr
               f3(i,0,nk)=f(i,0,nk)
               f3(i,nz,nk)=f(i,nz,nk)
            enddo

C     LATERALS
            do j=1,nz-1
               f3(0,j,nk)=f(0,j,nk)
            enddo

c     cierro el do de k             *
         enddo
         j=in(1)
         in(1)=in(2)
         in(2)=j

c     cierro el do de n             *
      enddo

c     Hacemos el nn/2  (solo hay uno)
      nm=(nn/2)-1
      nmq=dfloat(nm*nm)
      nk=nn-1
      do i=1,nr
         do j=1,nz-1
            h(i,j)=f(i,j,nk)*radq-
     $           lap(i,0,in(1))*f(0,j,nk) -
     $           4.d0*al2*(dz2(j,0)*f(i,0,nk)+dz2(j,nz)*f(i,nz,nk))
         enddo
      enddo

      call DGEMM('N','T',nr,nz-1,nz-1,1d0,h(1,1),nr,
     $     bzd(1,1),nz-1,0.d0,aux(1,1),nr)
      call DGEMM('N','N',nr,nz-1,nr,1d0,bxd(1,1,nm),nr,
     $     aux(1,1),nr,0.d0,w(1,1),nr)

      do i=1,nr
         do s=1,nz-1
            t=1d0/(vxd(i,nm)+4.d0*al2*vzd(s)-a*radq)
            w(i,s)=t*w(i,s)
         enddo
      enddo

      call DGEMM('N','T',nr,nz-1,nz-1,1d0,w(1,1),nr,
     $     azd(1,1),nz-1,0.d0,aux(1,1),nr)
      call DGEMM('N','N',nr,nz-1,nr,1d0,axd(1,1,nm),nr,
     $     aux(1,1),nr,0.d0,f3(1,1,nk),nr+1)

C     CONTORN DE LA FUNCIO F

C     TAPES
      do i=0,nr
         f3(i,0,nk)=f(i,0,nk)
         f3(i,nz,nk)=f(i,nz,nk)
      enddo

C     LATERALS
      do j=1,nz-1
         f3(0,j,nk)=f(0,j,nk)
      enddo

      RETURN
      END subroutine poisul
c     Es el poisson de ur+iutheta para el cilindro vertical,
c     Equivale a hacer ur -utheta (reales e imaginarias). Condiciones
c     de funcion conocida en las tapas y en el lateral

c=======================================================================

      subroutine poisup(a,f,f3,dz2,lap,axd,vxd,bxd,azd,vzd,bzd,RAD,ALT,
     &                  nr,nz,nn)
      implicit none
      integer nr,nz,nn
      real*8 dz2(0:nz,0:nz),lap(0:nr,0:nr,2)
      real*8 axd(1:nr,1:nr,0:nn/2+1),vxd(1:nr,0:nn/2+1),
     &     azd(1:nz-1,1:nz-1),bxd(1:nr,1:nr,0:nn/2+1),
     &     vzd(1:nz-1),bzd(1:nz-1,1:nz-1)
      real*8 RAD,ALT,t,al2,radq,ar,nmq
      real*8 a,f(0:nr,0:nz,0:nn-1),f3(0:nr,0:nz,0:nn-1)
      real*8 h(1:nr,1:nz-1),w(1:nr,1:nz-1),aux(1:nr,1:nz-1)
      integer in(2),n,nm,i,j,k,s,nk

      ar=RAD/ALT
      radq=RAD*RAD
      al2=ar*ar

c     Hacemos el n=0 que es impar en r
      in(1)=1
      in(2)=2
      do i=1,nr
         do j=1,nz-1
            h(i,j)=f(i,j,0)*radq-
     &           lap(i,0,in(1))*f(0,j,0)-
     &           4.d0*al2*(dz2(j,0)*f(i,0,0)+dz2(j,nz)*f(i,nz,0))
         enddo
      enddo

      call DGEMM('N','T',nr,nz-1,nz-1,1d0,h(1,1),nr,
     &     bzd(1,1),nz-1,0.d0,aux(1,1),nr)
      call DGEMM('N','N',nr,nz-1,nr,1d0,bxd(1,1,1),nr,
     &     aux(1,1),nr,0.d0,w(1,1),nr)

      do i=1,nr
         do s=1,nz-1
            t=1d0/(vxd(i,1)+4.d0*al2*vzd(s)-a*radq)
            w(i,s)=t*w(i,s)
         enddo
      enddo

      call DGEMM('N','T',nr,nz-1,nz-1,1d0,w(1,1),nr,
     &     azd(1,1),nz-1,0.d0,aux(1,1),nr)
      call DGEMM('N','N',nr,nz-1,nr,1d0,axd(1,1,1),nr,
     &     aux(1,1),nr,0.d0,f3(1,1,0),nr+1)

C     CONTORN DE LA FUNCIO F

C     TAPES
      do i=0,nr
         f3(i,0,0)=f(i,0,0)
         f3(i,nz,0)=f(i,nz,0)
      enddo

C     LATERALS
      do j=1,nz-1
         f3(0,j,0)=f(0,j,0)
      enddo
      j=in(1)
      in(1)=in(2)
      in(2)=j

c     Hacemos los otros n
      do n=1,nn-3,2
         nm=((n+1)/2)+1
         nmq=dfloat(nm*nm)
         do k=0,1
            nk=n+k
            do i=1,nr
               do j=1,nz-1
                  h(i,j)=f(i,j,nk)*radq-
     &                 lap(i,0,in(1))*f(0,j,nk) - 4.d0*al2*
     &                 (dz2(j,0)*f(i,0,nk)+dz2(j,nz)*f(i,nz,nk))
               enddo
            enddo

            call DGEMM('N','T',nr,nz-1,nz-1,1d0,h(1,1),nr,
     &           bzd(1,1),nz-1,0.d0,aux(1,1),nr)
            call DGEMM('N','N',nr,nz-1,nr,1d0,bxd(1,1,nm),nr,
     &           aux(1,1),nr,0.d0,w(1,1),nr)

            do i=1,nr
               do s=1,nz-1
                  t=1d0/(vxd(i,nm)+4.d0*al2*vzd(s)-a*radq)
                  w(i,s)=t*w(i,s)
               enddo
            enddo

            call DGEMM('N','T',nr,nz-1,nz-1,1d0,w(1,1),nr,
     &           azd(1,1),nz-1,0.d0,aux(1,1),nr)
            call DGEMM('N','N',nr,nz-1,nr,1d0,axd(1,1,nm),nr,
     &           aux(1,1),nr,0.d0,f3(1,1,nk),nr+1)

C     CONTORN DE LA FUNCIO F

C     TAPES
            do i=0,nr
               f3(i,0,nk)=f(i,0,nk)
               f3(i,nz,nk)=f(i,nz,nk)
            enddo

C     LATERALS
            do j=1,nz-1
               f3(0,j,nk)=f(0,j,nk)
            enddo

c     cierro el do de k             *
         enddo
         j=in(1)
         in(1)=in(2)
         in(2)=j

c     cierro el do de n             *
      enddo


c     Hacemos el nn/2 -1 (solo hay uno)
      nm=(nn/2)+1
      nmq=dfloat(nm*nm)
      nk=nn-1
      do i=1,nr
         do j=1,nz-1
            h(i,j)=f(i,j,nk)*radq-
     &           lap(i,0,in(1))*f(0,j,nk) -
     &           4.d0*al2*(dz2(j,0)*f(i,0,nk)+dz2(j,nz)*f(i,nz,nk))
         enddo
      enddo

      call DGEMM('N','T',nr,nz-1,nz-1,1d0,h(1,1),nr,
     &     bzd(1,1),nz-1,0.d0,aux(1,1),nr)
      call DGEMM('N','N',nr,nz-1,nr,1d0,bxd(1,1,nm),nr,
     &     aux(1,1),nr,0.d0,w(1,1),nr)

      do i=1,nr
         do s=1,nz-1
            t=1d0/(vxd(i,nm)+4.d0*al2*vzd(s)-a*radq)
            w(i,s)=t*w(i,s)
         enddo
      enddo

      call DGEMM('N','T',nr,nz-1,nz-1,1d0,w(1,1),nr,
     &     azd(1,1),nz-1,0.d0,aux(1,1),nr)
      call DGEMM('N','N',nr,nz-1,nr,1d0,axd(1,1,nm),nr,
     &     aux(1,1),nr,0.d0,f3(1,1,nk),nr+1)

C     CONTORN DE LA FUNCIO F

C     TAPES
      do i=0,nr
         f3(i,0,nk)=f(i,0,nk)
         f3(i,nz,nk)=f(i,nz,nk)
      enddo

C     LATERALS
      do j=1,nz-1
         f3(0,j,nk)=f(0,j,nk)
      enddo

      RETURN
      END subroutine poisup
