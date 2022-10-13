!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CreateGrid.F90                                 !
!    CONTAINS: subroutine CreateGrid                      !
!                                                         ! 
!    PURPOSE: Compute the indices, grid, grid metrics     !
!     and coefficients for differentiation                !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CreateGrid
      use param                                                 
      use mpih
      implicit none
      integer  :: j,k,n3mo,nclip,i,km,kp
      real :: tstr3, z2dp
      real :: x1,x2,x3,etain,delet
      real :: a33,a33m,a33p


      dx1=rext/real(n1m)
      dx2=rext2/real(n2m)
      dx3=alx3/real(n3m)
      dx1=1.0/dx1
      dx2=1.0/dx2
      dx3=1.0/dx3
      dx1q=dx1*dx1                                                      
      dx2q=dx2*dx2                                                      
      dx3q=dx3*dx3                                                      

!     Create vectors for indices

      do i=1,n1m
        imv(i)=i-1
        ipv(i)=i+1
        if(i.eq.1) imv(i)=n1m
        if(i.eq.n1m) ipv(i)=1
      enddo

!EP   Threaded for CPU affinity

!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(j)
      do j=1,n2m
        jmv(j)=j-1
        jpv(j)=j+1
        if(j.eq.1) jmv(j)=n2m
        if(j.eq.n2m) jpv(j)=1
      enddo
!$OMP  END PARALLEL DO

!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(j)
      do j = 1,n2+1
       jmhv(j) = mod(j,n2m/2+1)
       if(jmhv(j).eq.0) jmhv(j) = n2m/2 + 1
      enddo
!$OMP  END PARALLEL DO

      do k=1,n3m
        kmv(k)=k-1
        kpv(k)=k+1
        if(k.eq.1) kmv(k)=k
        if(k.eq.n3m) kpv(k)=k
      end do

      do k=1,n3m
        kpc(k)=kpv(k)-k
        kmc(k)=k-kmv(k)
        kup(k)=1-kpc(k)
        kum(k)=1-kmc(k)
      enddo

!     Create Coordinates

!     X- Uniform Coordinate

      do  i=1,n1
         x1=real(i-1)/real(n1m)
         tc(i)= rext*x1
        end do
      do i=1,n1m
        tm(i)=(tc(i)+tc(i+1))*0.5d0
      end do

!     Y- Uniform Coordinate
        do  j=1,n2
         x2=real(j-1)/real(n2m)
         rc(j)= rext2*x2
        end do
      do j=1,n2m
        rm(j)=(rc(j)+rc(j+1))*0.5d0
      end do
 
!     Z COORDINATE DEFINITION
!
!
!     UNIFORM GRID
!

      tstr3=tanh(str3)

      if (istr3.eq.0) then
        do k=1,n3
          x3=real(k-1)/real(n3m)
          etaz(k)=alx3*x3
          zc(k)=etaz(k)
        enddo
      endif


!       
!      CLUSTERING AT THE EXTERNAL RADIAL WALL 
!                       and  
!             CLUSTERING AT THE AXIS 
!      

        if (istr3.eq.4) then
         zc(1)=0.0d0
         do k=2,n3
          z2dp=float(2*k-n3-1)/float(n3m)
          zc(k)=(1+tanh(str3*z2dp)/tstr3)*0.5*alx3
          if(zc(k).lt.0.or.zc(k).gt.alx3)then
           write(*,*)'Forza la griglia: ','zc(',k,')=',zc(k)
           stop
          endif
         end do
        end if



      if(istr3.eq.6) then
      nclip = int(str3)
      n3mo = n3+nclip+nclip
      do k=1,n3mo
        etazm(k)=+cos(pi*(float(k)-0.5)/float(n3mo))
      end do
      do k=1,n3
        etaz(k)=etazm(k+nclip)
      end do
      delet = etaz(1)-etaz(n3)
      etain = etaz(1)
      do k=1,n3
        etaz(k)=etaz(k)/(0.5*delet)
      end do
      zc(1) = 0.
      do k=2,n3m
        zc(k) = alx3*(1.-etaz(k))*0.5
      end do
      zc(n3) = alx3
      endif
      
!-----------------------------------------
!
!     STAGGERED COORDINATES AND
!     METRIC QUANTITIES
!
      do k=1,n3m
        zm(k)=(zc(k)+zc(k+1))*0.5d0
        g3rm(k)=(zc(k+1)-zc(k))*dx3
      enddo
      do k=2,n3m
        g3rc(k)=(zc(k+1)-zc(k-1))*dx3*0.5d0
      enddo
      g3rc(1)=(zc(2)-zc(1))*dx3
      g3rc(n3)= (zc(n3)-zc(n3m))*dx3

      do k=1,n3m
        udx3m(k) = dx3/g3rm(k)
        udx3c(k) = dx3/g3rc(k)
      end do
      udx3c(n3) = dx3/g3rc(n3)

!====================================================
!     WRITE GRID INFORMATION
      if(myid.eq.0) then
      open(unit=78,file='axicor.out',status='unknown')
      do k=1,n3
        write(78,345) k,zc(k),zm(k),g3rc(k),g3rm(k)
      end do
      close(78)
 345  format(i4,4(2x,e23.15))

!===================================================
!
!     QUANTITIES FOR DERIVATIVES
!
      open(unit=78,file='fact3.out',status='unknown')
      do k=1,n3m
        write(78,*) k,udx3m(k),udx3c(k)
      end do
        write(78,*) n3,udx3m(n3m),udx3c(n3)
      close(78)

      endif




!     Create Metric and Coefficients for Derivatives


      am3ck(1)=0.d0
      ap3ck(1)=0.d0
      ac3ck(1)=1.d0
      am3ck(n3)=0.d0
      ap3ck(n3)=0.d0
      ac3ck(n3)=1.d0
      do k=2,n3m
      km=k-1
      kp=k+1
      a33=dx3q/g3rc(k)
      a33p=1.d0/g3rm(k)
      a33m=1.d0/g3rm(km)
      ap3ck(k)=a33*a33p
      am3ck(k)=a33*a33m
      ac3ck(k)=-(ap3ck(k)+am3ck(k))
      enddo

      do k=2,n3m-1
      kp=k+1
      km=k-1
      a33=dx3q/g3rm(k)
      a33p= +a33/g3rc(kp)
      a33m= +a33/g3rc(k)
      ap3sk(k)=a33p
      am3sk(k)=a33m
      ac3sk(k)=-(ap3sk(k)+am3sk(k))
      ap3ssk(k)=ap3sk(k)
      am3ssk(k)=am3sk(k)
      ac3ssk(k)=ac3sk(k)
      enddo

      k=1
      kp=k+1
      a33=dx3q/g3rm(k)
      a33p= +a33/g3rc(kp)
      a33m= +a33/g3rc(k)
      ap3sk(k)=a33p
      am3sk(k)=0.d0
      ac3sk(k)=-(a33p+inslws*a33m*2.d0)
      ap3ssk(k)=ap3sk(k)
      am3ssk(k)=am3sk(k)
      ac3ssk(k)=-(a33p+a33m*2.d0)
     
      k=n3m
      kp=k+1
      a33=dx3q/g3rm(k)
      a33p= +a33/g3rc(kp)
      a33m= +a33/g3rc(k)
      am3sk(k)=a33m
      ap3sk(k)=0.d0
      ac3sk(k)=-(a33m+inslwn*a33p*2.d0)
      ap3ssk(k)=ap3sk(k)
      am3ssk(k)=am3sk(k)
      ac3ssk(k)=-(a33m+a33p*2.d0)

      return                                                            
      end                                                               
