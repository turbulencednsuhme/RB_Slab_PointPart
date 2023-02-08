!------------------------------------------------------------
!     update particle position based on velocity
!     find out indices and distribute to slaves
!------------------------------------------------------------
      subroutine UpdatePointParticlePosition

      USE param
      USE pointparticle

      IMPLICIT NONE

      real      :: pos(3),qInt(7)
      real      :: vx_new,vy_new,vz_new
      real      :: vx_old,vy_old,vz_old
      real      :: q_new, q_old
      integer   :: inp,merr

      
      do inp=1,Npointpart

       xpo(inp,1) = xp(inp,1)
       xpo(inp,2) = xp(inp,2)
       xpo(inp,3) = xp(inp,3) 
       xpo(inp,4) = xp(inp,4) 

       vx_new=vxp (inp); vy_new=vyp (inp); vz_new=vzp (inp)
       vx_old=vxpo(inp); vy_old=vypo(inp); vz_old=vzpo(inp)      
       q_new = dtempp(inp); q_old = dtemppo(inp)

       if(ONparticle.eq.0)then
        xp(inp,1) = xp(inp,1) + dt*vx_new
        xp(inp,2) = xp(inp,2) + dt*vy_new
        xp(inp,3) = xp(inp,3) + dt*vz_new
        xp(inp,4) = xp(inp,4) + dt*q_new
       else
        xp(inp,1) = xp(inp,1) + 0.5d0*dt*(3.d0*vx_new - vx_old)
        xp(inp,2) = xp(inp,2) + 0.5d0*dt*(3.d0*vy_new - vy_old)
        xp(inp,3) = xp(inp,3) + 0.5d0*dt*(3.d0*vz_new - vz_old)
        xp(inp,4) = xp(inp,4) + 0.5d0*dt*(3.d0*q_new - q_old)
       end if

!     =================================================================
!               Periodic boundary corrections for particles
!     =================================================================


!      X-boundary check
       if(xp(inp,1).lt.tm(1))   xp(inp,1)=tm(n1m)-(tm(1)-xp(inp,1))
       if(xp(inp,1).gt.tm(n1m)) xp(inp,1)=tm(1)+(xp(inp,1)-tm(n1m))

!      Y-boundary check
       if(xp(inp,2).lt.rm(1))   xp(inp,2)=rm(n2m)-(rm(1)-xp(inp,2))
       if(xp(inp,2).gt.rm(n2m)) xp(inp,2)=rm(1)+(xp(inp,2)-rm(n2m))

!      Z-boundary check
!       if(xp(inp,3).lt.zm(1))   xp(inp,3)=zm(n3m)-(zm(1)-xp(inp,3))
!       if(xp(inp,3).gt.zm(n3m)) xp(inp,3)=zm(1)+(xp(inp,3)-zm(n3m))

!     =================================================================
!                Update old variables by all processors
!     =================================================================

       vxpo(inp)=vxp(inp)
       vypo(inp)=vyp(inp)
       vzpo(inp)=vzp(inp)
       dtemppo(inp)=dtempp(inp)


       kalbo1(inp)=kalb1(inp)
       kalbo2(inp)=kalb2(inp)
       kalbo3(inp)=kalb3(inp)

       qValo1(inp)=qVal1(inp)
       qValo2(inp)=qVal2(inp)
       qValo3(inp)=qVal3(inp)

       vorto1(inp)=vort1(inp)
       vorto2(inp)=vort2(inp)
       vorto3(inp)=vort3(inp)

      end do
     
      return
      end
