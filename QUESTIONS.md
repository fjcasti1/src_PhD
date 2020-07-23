## QUESTIONS
Collection of questions (so I don't make a mess in my desk and later forget)

### z domain in DNS and VTK
In DNS the z domain is [-Gamma,Gamma] and it is projected to a uniform grid of
domain [-Gamma/2,Gamma/2] in the VTK file. Correct?

### Erase one loop in freeSurfaceTop.f
From:
```Fortran
     do i=1,nr
       ftilde(i,2*k  )   = f(i,0,2*k  )
     enddo
     ! Solve
     call dgetrs('N',nr,1,Mtilde(1,1,1),nr,ipiv1,
&     ftilde(1,2*k  ),nr,info)
```
To:
```Fortran
     ! Solve
     call dgetrs('N',nr,1,Mtilde(1,1,1),nr,ipiv1,
&     f(1,0,2*k  ),nr,info)
```

### Regularization of vt bottom/side-walls in freeSurfTop.f
From:
```Fortran
c   Bottom lid: Uniform rotation with rate Omega=1+Ro*sin(wf*t) rad/s

    do i=0,nr
       vt(i,nz,0)=-r(i)*(pnu*Re/rad)*(1d0+Ro*dsin(wf*tps))
    enddo

c   Side wall: no-slip BC with a gaussian regularization to mantain
c   contuniuity

    mu  = -alt ! Center of the Gaussian
    reg = 1d-2
    do j=1,nz-1
       vt(0,j,0)=-(pnu*Re/rad)*(1d0+Ro*dsin(wf*tps))*
   &               dexp(-(z(j)-mu)**2d0/(2d0*reg**2d0))
    enddo
```
To:
```Fortran
c   Bottom lid: Uniform rotation with rate Omega=1+Ro*sin(wf*t) rad/s

    do i=0,nr
       vt(i,nz,0)=-r(i)*(pnu*Re/rad)*(1d0+Ro*dsin(wf*tps))
    enddo

c   Side wall: no-slip BC with a gaussian regularization to mantain
c   contuniuity

    mu  = -alt ! Center of the Gaussian
    reg = 1d-2
    do j=1,nz-1
       vt(0,j,0)=vt(0,nz,0)*dexp(-(z(j)-mu)**2d0/(2d0*reg**2d0))
    enddo
```
