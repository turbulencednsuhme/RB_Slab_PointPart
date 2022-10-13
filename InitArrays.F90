      subroutine InitArrays
      use param
      use mpi_param, only: kstart,kend
      use local_arrays
      use stat_arrays
      use AuxiliaryRoutines
      implicit none

      call AllocateReal1DArray(tc,1,n1)
      call AllocateReal1DArray(tm,1,n1)

      call AllocateReal1DArray(rc,1,n2)
      call AllocateReal1DArray(rm,1,n2)

      call AllocateReal1DArray(zc,1,n3)
      call AllocateReal1DArray(zm,1,n3)

      call AllocateReal1DArray(g3rc,1,n3)
      call AllocateReal1DArray(g3rm,1,n3)

      call AllocateReal1DArray(udx3c,1,n3)
      call AllocateReal1DArray(udx3m,1,n3)
      call AllocateReal1DArray(etaz,1,n3)
      call AllocateReal1DArray(etazm,1,n3+1000)

      call AllocateReal1DArray(ak1,1,n1)
      call AllocateReal1DArray(ao,1,n1)

      call AllocateReal1DArray(ak2,1,n2)
      call AllocateReal1DArray(ap,1,n2)

      call AllocateReal1DArray(amphk,1,n3)
      call AllocateReal1DArray(apphk,1,n3)
      call AllocateReal1DArray(acphk,1,n3)

      call AllocateReal1DArray(trigx1,1,3*n2/2+1)

      call AllocateInt1DArray(imv,1,n1)
      call AllocateInt1DArray(ipv,1,n1)
      call AllocateInt1DArray(jmv,1,n2)
      call AllocateInt1DArray(jpv,1,n2)
      call AllocateInt1DArray(kmv,1,n3)
      call AllocateInt1DArray(kpv,1,n3)

      call AllocateInt1DArray(jmhv,1,n2+1)

      call AllocateInt1DArray(kmc,1,n3)
      call AllocateInt1DArray(kpc,1,n3)
      call AllocateInt1DArray(kup,1,n3)
      call AllocateInt1DArray(kum,1,n3)

      call AllocateReal1DArray(ap3j,1,n2)
      call AllocateReal1DArray(ac3j,1,n2)
      call AllocateReal1DArray(am3j,1,n2)

      call AllocateReal1DArray(ap3ck,1,n3)
      call AllocateReal1DArray(am3ck,1,n3)
      call AllocateReal1DArray(ac3ck,1,n3)

      call AllocateReal1DArray(ap3sk,1,n3)
      call AllocateReal1DArray(am3sk,1,n3)
      call AllocateReal1DArray(ac3sk,1,n3)

      call AllocateReal1DArray(ap3ssk,1,n3)
      call AllocateReal1DArray(am3ssk,1,n3)
      call AllocateReal1DArray(ac3ssk,1,n3)

      call AllocateReal2DArray(denbs,1,n1,1,n2)
      call AllocateReal2DArray(denbn,1,n1,1,n2)
      call AllocateReal2DArray(uwtb,1,n1,1,n2)
      call AllocateReal2DArray(lwtb,1,n1,1,n2)

      call AllocateReal2DArray(vx_me,1,n2m,kstart,kend)
      call AllocateReal2DArray(vy_me,1,n2m,kstart,kend)
      call AllocateReal2DArray(vz_me,1,n2m,kstart,kend)
      call AllocateReal2DArray(temp_me,1,n2m,kstart,kend)

      call AllocateReal2DArray(vx_rms,1,n2m,kstart,kend)
      call AllocateReal2DArray(vy_rms,1,n2m,kstart,kend)
      call AllocateReal2DArray(vz_rms,1,n2m,kstart,kend)
      call AllocateReal2DArray(temp_rms,1,n2m,kstart,kend)

      call AllocateReal2DArray(tempvz_me,1,n2m,kstart,kend)

      call AllocateReal2DArray(disste,1,n2m,kstart,kend)
      call AllocateReal2DArray(dissth,1,n2m,kstart,kend)

      call AllocateReal3DArray(vx,1,n1,1,n2,kstart-1,kend+1)
      call AllocateReal3DArray(vy,1,n1,1,n2,kstart-1,kend+1)
      call AllocateReal3DArray(vz,1,n1,1,n2,kstart-1,kend+1)
      call AllocateReal3DArray(pr,1,n1,1,n2,kstart-1,kend+1)
      call AllocateReal3DArray(temp,1,n1,1,n2,kstart-1,kend+1)

      call AllocateReal3DArray(dph,1,n1,1,n2+1,kstart-1,kend+1)

      call AllocateReal3DArray(rhs,1,n1,1,n2,kstart,kend)
      call AllocateReal3DArray(dq,1,n1,1,n2,kstart,kend)
      call AllocateReal3DArray(qcap,1,n1,1,n2,kstart,kend)
      call AllocateReal3DArray(hro,1,n1,1,n2,kstart,kend)

      call AllocateReal3DArray(ru1,1,n1,1,n2,kstart,kend)
      call AllocateReal3DArray(ru2,1,n1,1,n2,kstart,kend)
      call AllocateReal3DArray(ru3,1,n1,1,n2,kstart,kend)
      call AllocateReal3DArray(ruro,1,n1,1,n2,kstart,kend)
 

      return 
      end 
