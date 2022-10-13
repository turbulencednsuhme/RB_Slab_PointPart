!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: DeallocateArrays.F90                           !
!    CONTAINS: subroutine DeallocateArrays                !
!                                                         ! 
!    PURPOSE: Finalization routine. Deallocates all       !
!     variables used in the code                          !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine DeallocateArrays
      use param
      use mpi_param, only: kstart,kend
      use local_arrays
      use stat_arrays
      use AuxiliaryRoutines
      implicit none


      call DestroyReal1DArray(tc)
      call DestroyReal1DArray(tm)

      call DestroyReal1DArray(rc)
      call DestroyReal1DArray(rm)

      call DestroyReal1DArray(zc)
      call DestroyReal1DArray(zm)

      call DestroyReal1DArray(g3rc)
      call DestroyReal1DArray(g3rm)

      call DestroyReal1DArray(udx3c)
      call DestroyReal1DArray(udx3m)
      call DestroyReal1DArray(etaz)
      call DestroyReal1DArray(etazm)

      call DestroyReal1DArray(ak1)
      call DestroyReal1DArray(ao)

      call DestroyReal1DArray(ak2)
      call DestroyReal1DArray(ap)

      call DestroyReal1DArray(amphk)
      call DestroyReal1DArray(apphk)
      call DestroyReal1DArray(acphk)

      call DestroyReal1DArray(trigx1)

      call DestroyInt1DArray(imv)
      call DestroyInt1DArray(ipv)
      call DestroyInt1DArray(jmv)
      call DestroyInt1DArray(jpv)
      call DestroyInt1DArray(kmv)
      call DestroyInt1DArray(kpv)

      call DestroyInt1DArray(jmhv)

      call DestroyInt1DArray(kmc)
      call DestroyInt1DArray(kpc)
      call DestroyInt1DArray(kup)
      call DestroyInt1DArray(kum)

      call DestroyReal1DArray(ap3j)
      call DestroyReal1DArray(ac3j)
      call DestroyReal1DArray(am3j)

      call DestroyReal1DArray(ap3ck)
      call DestroyReal1DArray(am3ck)
      call DestroyReal1DArray(ac3ck)

      call DestroyReal1DArray(ap3sk)
      call DestroyReal1DArray(am3sk)
      call DestroyReal1DArray(ac3sk)

      call DestroyReal1DArray(ap3ssk)
      call DestroyReal1DArray(am3ssk)
      call DestroyReal1DArray(ac3ssk)

      call DestroyReal2DArray(denbs)
      call DestroyReal2DArray(denbn)
      call DestroyReal2DArray(uwtb)
      call DestroyReal2DArray(lwtb)

      call DestroyReal2DArray(vx_me)
      call DestroyReal2DArray(vy_me)
      call DestroyReal2DArray(vz_me)
      call DestroyReal2DArray(temp_me)

      call DestroyReal2DArray(vx_rms)
      call DestroyReal2DArray(vy_rms)
      call DestroyReal2DArray(vz_rms)
      call DestroyReal2DArray(temp_rms)

      call DestroyReal2DArray(tempvz_me)

      call DestroyReal2DArray(disste)
      call DestroyReal2DArray(dissth)

      call DestroyReal3DArray(vx)
      call DestroyReal3DArray(vy)
      call DestroyReal3DArray(vz)
      call DestroyReal3DArray(pr)
      call DestroyReal3DArray(temp)

      call DestroyReal3DArray(dph)

      call DestroyReal3DArray(rhs)
      call DestroyReal3DArray(dq)
      call DestroyReal3DArray(qcap)
      call DestroyReal3DArray(hro)

      call DestroyReal3DArray(ru1)
      call DestroyReal3DArray(ru2)
      call DestroyReal3DArray(ru3)
      call DestroyReal3DArray(ruro)

      return
      end
