      subroutine ReadInputFile
      use param
      implicit none
      character(len=4) :: dummy


      open(unit=15,file='bou.in',status='old')
        read(15,301) dummy
        read(15,*) n1m,n2m,n3m,nsst,nwrit,nread
        read(15,301) dummy
        read(15,*) ntst,tprint,tpin,tmax,walltimemax
        read(15,301) dummy
        read(15,*) alx3,istr3,str3
        read(15,301) dummy
        read(15,*) rext,rext2
        read(15,301) dummy
        read(15,*) ray,pra,dt,resid,cflmax
        read(15,301) dummy
        read(15,*) tsta,starea
        read(15,301) dummy
        read(15,*) inslws,inslwn
        read(15,301) dummy       
        read(15,*) idtv,dtmax,cfllim  
        read(15,301) dummy       
        read(15,*) tframe
301     format(a4)                
      close(15)
      
      ren = dsqrt(ray/pra)
      pec = dsqrt(pra*ray)

      pi=2.d0*dasin(1.d0)                          

      n1=n1m+1 
      n2=n2m+1
      n3=n3m+1

      n1mh=n1m/2+1
      n1mp=n1mh+1
      n2mh=n2m/2+1
      n2mp=n2mh+1

!     Switch from Hours to Seconds

      walltimemax=walltimemax*60.*60.
      return
      end
