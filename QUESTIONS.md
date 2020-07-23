## QUESTIONS
Collection of questions (so I don't make a mess in my desk and later forget)

### z domain in DNS and VTK
In DNS the z domain is [-Gamma,Gamma] and it is projected to a uniform grid of
domain [-Gamma/2,Gamma/2] in the VTK file. Correct?

### Erase one loop in freeSurfaceTop.f
From:
```
     do i=1,nr
       ftilde(i,2*k  )   = f(i,0,2*k  )
     enddo
     ! Solve
     call dgetrs('N',nr,1,Mtilde(1,1,1),nr,ipiv1,
&     ftilde(1,2*k  ),nr,info)
```
To:
```
     ! Solve
     call dgetrs('N',nr,1,Mtilde(1,1,1),nr,ipiv1,
&     f(1,0,2*k  ),nr,info)
```