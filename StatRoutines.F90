!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: StatRoutines.F90                               !
!    CONTAINS: subroutine CalcStats,WriteStats            !
!                                                         ! 
!    PURPOSE: Calculates and writes out statistics for    !
!     the flow field. All quantities are averaged in the  !
!     two horizontal (homogeneous) directions             !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CalcStats
      use param
      use mpi_param, only: kstart,kend
      use local_arrays, only: vx,vy,vz,temp
      use stat_arrays
      use mpih
      implicit none
      real :: my_vy_rms_vol,my_vx_rms_vol,my_vz_rms_vol
      real :: my_vxvyvz_rms_vol,vy_rms_vol,vx_rms_vol
      real :: vz_rms_vol,vxvyvz_rms_vol
      real :: lvol,volt,usn1m,rradpr
      integer :: i,j,k

      timeint_cdsp = timeint_cdsp + 1

768   format(6e20.10)

      pi = 2.*asin(1.)

      usn1m = 1.0/n1m

      my_vx_rms_vol = 0.0d0
      my_vy_rms_vol = 0.0d0
      my_vz_rms_vol = 0.0d0
      my_vxvyvz_rms_vol = 0.0d0

      do k=kstart,kend
      lvol = g3rm(k)
!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(i,j)
!$OMP& REDUCTION(+:my_vx_rms_vol,my_vy_rms_vol)
!$OMP& REDUCTION(+:my_vz_rms_vol,my_vxvyvz_rms_vol)
         do j=1,n2m
            do i=1,n1m
               vx_me(j,k) = vx_me(j,k) + vx(i,j,k)*usn1m
               vy_me(j,k) = vy_me(j,k) + vy(i,j,k)*usn1m
               vz_me(j,k) = vz_me(j,k) + vz(i,j,k)*usn1m
               temp_me(j,k) = temp_me(j,k) + temp(i,j,k)*usn1m
               vx_rms(j,k) = vx_rms(j,k) + vx(i,j,k)**2*usn1m
               vy_rms(j,k) = vy_rms(j,k) + vy(i,j,k)**2*usn1m
               vz_rms(j,k) = vz_rms(j,k) + vz(i,j,k)**2*usn1m
               temp_rms(j,k) = temp_rms(j,k) + temp(i,j,k)**2*usn1m
               tempvz_me(j,k) = tempvz_me(j,k) + temp(i,j,k)*(vz(i,j,k+1)+vz(i,j,k))*usn1m*0.5

               my_vx_rms_vol = my_vx_rms_vol + lvol*vx(i,j,k)**2
               my_vy_rms_vol = my_vy_rms_vol + lvol*vy(i,j,k)**2
               my_vz_rms_vol = my_vz_rms_vol + lvol*vz(i,j,k)**2
               my_vxvyvz_rms_vol = my_vxvyvz_rms_vol + lvol*( &
                (vx(i+1,j,k)+vx(i,j,k))**2+ &
                (vy(i,j+1,k)+vy(i,j,k))**2+ &
                (vz(i,j,k+1)+vz(i,j,k))**2)*0.25 
            end do
         end do
!$OMP  END PARALLEL DO
      end do

      call MPI_REDUCE(my_vx_rms_vol,vx_rms_vol,1,MDP,MPI_SUM,0, &
       MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(my_vy_rms_vol,vy_rms_vol,1,MDP,MPI_SUM,0, &
       MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(my_vz_rms_vol,vz_rms_vol,1,MDP,MPI_SUM,0, &
       MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(my_vxvyvz_rms_vol,vxvyvz_rms_vol,1,MDP, &
       MPI_SUM,0,MPI_COMM_WORLD,ierr)


      if(myid.eq.0) then
        volt = 1.d0/(real(n1m)*real(n2m)*real(n3m))
        rradpr=dsqrt(ray/pra)
        vx_rms_vol=dsqrt(vx_rms_vol*volt)*rradpr
        vy_rms_vol=dsqrt(vy_rms_vol*volt)*rradpr
        vz_rms_vol=dsqrt(vz_rms_vol*volt)*rradpr
        vxvyvz_rms_vol=dsqrt(vxvyvz_rms_vol*volt)*rradpr
      
       write(94,768) time,vx_rms_vol,vy_rms_vol,vz_rms_vol, &
       vxvyvz_rms_vol
      endif

      return  
      end

!***********************************************************************

      subroutine WriteStats
      use mpih
      use param
      use mpi_param, only: kstart,kend
      use stat_arrays
      use hdf5

      implicit none

      integer hdf_error

      integer(HID_T) :: file_id
      integer(HID_T) :: filespace
      integer(HID_T) :: slabspace
      integer(HID_T) :: memspace

      integer(HID_T) :: dset_vxme
      integer(HID_T) :: dset_vyme
      integer(HID_T) :: dset_vzme
      integer(HID_T) :: dset_tempme

      integer(HID_T) :: dset_vxrms
      integer(HID_T) :: dset_vyrms
      integer(HID_T) :: dset_vzrms
      integer(HID_T) :: dset_temprms
      integer(HID_T) :: dset_tempvzme

      integer(HID_T) :: dset_kindiss
      integer(HID_T) :: dset_thediss

      integer(HSIZE_T) :: dims(2)

      integer(HSIZE_T) :: dims_grid(1)
      integer(HID_T) :: dset_grid
      integer(HID_T) :: dspace_grid

      integer(HID_T) :: plist_id
      integer(HSIZE_T), dimension(2) :: data_count  
      integer(HSSIZE_T), dimension(2) :: data_offset 

      integer :: comm, info
      integer :: ndims

      character*30 filnam1
      character*30 filnamgrid

!RO   Sort out MPI definitions

      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

!RO   Form the name of the file

      filnam1 = 'stafield_data.h5'

!RO   Set offsets and element counts
   
      ndims = 2

      dims(1)=n2m
      dims(2)=n3m

      call h5screate_simple_f(ndims, dims, filespace, hdf_error)

      data_count(1) = n2m
      data_count(2) = kend-kstart+1

      data_offset(1) = 0
      data_offset(2) = kstart-1

!RO   Open file and create dataspace

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id,hdf_error)

      call h5pset_fapl_mpio_f(plist_id, comm, info,hdf_error)

      call h5fcreate_f(filnam1, H5F_ACC_TRUNC_F, file_id, hdf_error, access_prp=plist_id)

      call h5pclose_f(plist_id, hdf_error)

!RO   Create the datasets with default properties.

      call h5dcreate_f(file_id, 'Vx_mean', H5T_NATIVE_DOUBLE, &
                      filespace, dset_vxme, hdf_error)

      call h5dcreate_f(file_id, 'Vy_mean', H5T_NATIVE_DOUBLE, &
                      filespace, dset_vyme, hdf_error)

      call h5dcreate_f(file_id, 'Vz_mean', H5T_NATIVE_DOUBLE, &
                      filespace, dset_vzme, hdf_error)

      call h5dcreate_f(file_id, 'temp_mean', H5T_NATIVE_DOUBLE, &
                      filespace, dset_tempme, hdf_error)

      call h5dcreate_f(file_id, 'Vx_rms', H5T_NATIVE_DOUBLE, &
                      filespace, dset_vxrms, hdf_error)

      call h5dcreate_f(file_id, 'Vy_rms', H5T_NATIVE_DOUBLE, &
                      filespace, dset_vyrms, hdf_error)

      call h5dcreate_f(file_id, 'Vz_rms', H5T_NATIVE_DOUBLE, &
                      filespace, dset_vzrms, hdf_error)

      call h5dcreate_f(file_id, 'temp_rms', H5T_NATIVE_DOUBLE, &
                      filespace, dset_temprms, hdf_error)

      call h5dcreate_f(file_id, 'kindiss_mean', H5T_NATIVE_DOUBLE, &
                      filespace, dset_kindiss, hdf_error)

      call h5dcreate_f(file_id, 'thediss_mean', H5T_NATIVE_DOUBLE, &
                      filespace, dset_thediss, hdf_error)

      call h5dcreate_f(file_id, 'tempvz_mean', H5T_NATIVE_DOUBLE, &
                      filespace, dset_tempvzme, hdf_error)

      call h5sclose_f(filespace, hdf_error)

!RO   Create dataspace in memory

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 


!RO   Select hyperslab  and then write it

!RO   Q1me


      call h5dget_space_f(dset_vxme, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dwrite_f(dset_vxme, H5T_NATIVE_DOUBLE, &
         vx_me(1:n2m,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

!EP   Q2me

      call h5dget_space_f(dset_vyme, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F,&
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,&
                              hdf_error)
       call h5dwrite_f(dset_vyme, H5T_NATIVE_DOUBLE,&
         vy_me(1:n2m,kstart:kend), dims, &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace,& 
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

!EP   Q3me

      call h5dget_space_f(dset_vzme, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F,&
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,&
                              hdf_error)
       call h5dwrite_f(dset_vzme, H5T_NATIVE_DOUBLE,&
         vz_me(1:n2m,kstart:kend), dims, &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace,&
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

!EP   tempme

      call h5dget_space_f(dset_tempme, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F,&
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,&
                              hdf_error)
       call h5dwrite_f(dset_tempme, H5T_NATIVE_DOUBLE,&
         temp_me(1:n2m,kstart:kend), dims, &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

!EP   Q1rms

      call h5dget_space_f(dset_vxrms, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dwrite_f(dset_vxrms, H5T_NATIVE_DOUBLE, &
         vx_rms(1:n2m,kstart:kend), dims, &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

!EP   Q2rms

      call h5dget_space_f(dset_vyrms, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dwrite_f(dset_vyrms, H5T_NATIVE_DOUBLE, &
         vy_rms(1:n2m,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

!EP   Q3rms

      call h5dget_space_f(dset_vzrms, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dwrite_f(dset_vzrms, H5T_NATIVE_DOUBLE, &
         vz_rms(1:n2m,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

!EP   temprms

      call h5dget_space_f(dset_temprms, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dwrite_f(dset_temprms, H5T_NATIVE_DOUBLE, &
         temp_rms(1:n2m,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

!EP   kindiss

      call h5dget_space_f(dset_kindiss, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dwrite_f(dset_kindiss, H5T_NATIVE_DOUBLE, &
         disste(1:n2m,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

!EP   thediss

      call h5dget_space_f(dset_thediss, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F,&
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,&
                              hdf_error)
       call h5dwrite_f(dset_thediss, H5T_NATIVE_DOUBLE,&
         dissth(1:n2m,kstart:kend), dims, &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

!EP   tempvz

      call h5dget_space_f(dset_tempvzme, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dwrite_f(dset_tempvzme, H5T_NATIVE_DOUBLE, &
         tempvz_me(1:n2m,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

!RO   Close properties and datasets

      call h5dclose_f(dset_vxme, hdf_error)
      call h5dclose_f(dset_vyme, hdf_error)
      call h5dclose_f(dset_vzme, hdf_error)
      call h5dclose_f(dset_tempme, hdf_error)

      call h5dclose_f(dset_vxrms, hdf_error)
      call h5dclose_f(dset_vyrms, hdf_error)
      call h5dclose_f(dset_vzrms, hdf_error)
      call h5dclose_f(dset_temprms, hdf_error)

      call h5dclose_f(dset_kindiss, hdf_error)
      call h5dclose_f(dset_thediss, hdf_error)
      call h5dclose_f(dset_tempvzme, hdf_error)

      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

!RO   Write the grid & statistics information
!RO   only if master process

      if (myid.eq.0) then

      ndims=1

      filnamgrid = 'stafield_master.h5'
      call h5fcreate_f(filnamgrid,H5F_ACC_TRUNC_F, file_id, hdf_error)

!RO   Write amount of averages 

      dims_grid(1)=1
      call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)

      call h5dcreate_f(file_id, 'averaging_time', H5T_NATIVE_INTEGER, &
                      dspace_grid, dset_grid, hdf_error)

      call h5dwrite_f(dset_grid, H5T_NATIVE_INTEGER, timeint_cdsp, &
             dims_grid,hdf_error)

      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

!RO   Write Reynolds number

      dims_grid(1)=1
      call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)

      call h5dcreate_f(file_id, 'Ra', H5T_NATIVE_DOUBLE, &
                      dspace_grid, dset_grid, hdf_error)

      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, ray,dims_grid,hdf_error)

      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)
          

!EP   Write Prandtl number

      dims_grid(1)=1
      call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)

      call h5dcreate_f(file_id, 'Pr', H5T_NATIVE_DOUBLE,dspace_grid, dset_grid, hdf_error)

      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, pra,dims_grid,hdf_error)

      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)
     &     

!RO   Write the grid information 

      dims_grid(1)=n2m
      call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)

      call h5dcreate_f(file_id, 'R_cordin', H5T_NATIVE_DOUBLE, &
                      dspace_grid, dset_grid, hdf_error)

      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, rm(1:n2m), &
             dims_grid,hdf_error)

      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      dims_grid(1)=n3m
      call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_id, 'Z_cordin', H5T_NATIVE_DOUBLE, &
                      dspace_grid, dset_grid, hdf_error)

      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, zm(1:n3m), &
              dims_grid, hdf_error)


      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

!RO   Close file

      call h5fclose_f(file_id, hdf_error)

      endif

      return  
      end
 
!***********************************************************************

      subroutine initstst
      use param
      use mpi_param, only: kstart,kend
      use stat_arrays
      use mpih
      use hdf5
      implicit none
      integer :: j,k

      integer hdf_error

      integer(HID_T) :: file_id
      integer(HID_T) :: memspace
      integer(HID_T) :: slabspace

      integer(HID_T) :: dset_vxme
      integer(HID_T) :: dset_vyme
      integer(HID_T) :: dset_vzme
      integer(HID_T) :: dset_tempme

      integer(HID_T) :: dset_vxrms
      integer(HID_T) :: dset_vyrms
      integer(HID_T) :: dset_vzrms
      integer(HID_T) :: dset_temprms
      integer(HID_T) :: dset_tempvzme

      integer(HID_T) :: dset_kindiss
      integer(HID_T) :: dset_thediss

      integer(HSIZE_T) :: dims(2)

      integer(HSIZE_T) :: dims_grid(1)
      integer(HID_T) :: dset_grid
      integer(HID_T) :: dspace_grid

      integer(HID_T) :: plist_id
      integer(HID_T) :: plist_full
      integer(HSIZE_T), dimension(2) :: data_count  
      integer(HSSIZE_T), dimension(2) :: data_offset 

      integer :: comm, info
      integer :: ndims

      character*30 filnam1
      character*30 filnamgrid

!EP   Read or initialize stat arrays

      if(starea.eq.1) then

!RO   Sort out MPI definitions

      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

!RO   Form the name of the file

      filnam1 = 'stafield_data.h5'
      filnamgrid = 'stafield_master.h5'

!RO   Set offsets and element counts
   
      ndims = 2

      dims(1)=n2m
      dims(2)=n3m

      data_count(1)=n2m
      data_count(2)=kend-kstart+1

      data_offset(1) = 0
      data_offset(2) = kstart-1
         
!RO   Open file and create dataspace

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_full, hdf_error)
      call h5pset_fapl_mpio_f(plist_full, comm, info, hdf_error)

      call h5fopen_f(filnam1, H5F_ACC_RDONLY_F, file_id, hdf_error, &
                       access_prp=plist_full)


!RO   Create the datasets with default properties.

      call h5dopen_f(file_id, 'Vx_mean',dset_vxme, hdf_error)
      call h5dopen_f(file_id, 'Vy_mean',dset_vyme, hdf_error)
      call h5dopen_f(file_id, 'Vz_mean',dset_vzme, hdf_error)
      call h5dopen_f(file_id, 'temp_mean',dset_tempme, hdf_error)

      call h5dopen_f(file_id, 'Vx_rms',dset_vxrms, hdf_error)
      call h5dopen_f(file_id, 'Vy_rms',dset_vyrms, hdf_error)
      call h5dopen_f(file_id, 'Vz_rms',dset_vzrms, hdf_error)

      call h5dopen_f(file_id, 'temp_rms',dset_temprms, hdf_error)

      call h5dopen_f(file_id, 'kindiss_mean',dset_kindiss, hdf_error)

      call h5dopen_f(file_id, 'thediss_mean',dset_thediss, hdf_error)

      call h5dopen_f(file_id, 'tempvz_mean',dset_tempvzme, hdf_error)

!RO   Create dataspace in memory

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

!RO   Select hyperslab  and then read it

!RO   Q1me

      call h5dget_space_f(dset_vxme, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dread_f(dset_vxme, H5T_NATIVE_DOUBLE, &
         vx_me(1:n2m,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

!EP   Q2me

      call h5dget_space_f(dset_vyme, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dread_f(dset_vyme, H5T_NATIVE_DOUBLE, &
         vy_me(1:n2m,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

!EP   Q3me

      call h5dget_space_f(dset_vzme, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,&
                              hdf_error)
       call h5dread_f(dset_vzme, H5T_NATIVE_DOUBLE, &
         vz_me(1:n2m,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

!EP   tempme

      call h5dget_space_f(dset_tempme, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dread_f(dset_tempme, H5T_NATIVE_DOUBLE, &
         temp_me(1:n2m,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

!EP   Q1rms

      call h5dget_space_f(dset_vxrms, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,&
                              hdf_error)
       call h5dread_f(dset_vxrms, H5T_NATIVE_DOUBLE, &
         vx_rms(1:n2m,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

!EP   Q2rms

      call h5dget_space_f(dset_vyrms, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dread_f(dset_vyrms, H5T_NATIVE_DOUBLE, &
         vy_rms(1:n2m,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

!EP   Q3rms

      call h5dget_space_f(dset_vzrms, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dread_f(dset_vzrms, H5T_NATIVE_DOUBLE, &
         vz_rms(1:n2m,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

!EP   temprms

      call h5dget_space_f(dset_temprms, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dread_f(dset_temprms, H5T_NATIVE_DOUBLE, &
         temp_rms(1:n2m,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)


!EP   kindiss

      call h5dget_space_f(dset_kindiss, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dread_f(dset_kindiss, H5T_NATIVE_DOUBLE, &
         disste(1:n2m,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

!EP   thediss

      call h5dget_space_f(dset_thediss, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dread_f(dset_thediss, H5T_NATIVE_DOUBLE, &
         dissth(1:n2m,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

!EP   tempvzme

      call h5dget_space_f(dset_tempvzme, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dread_f(dset_tempvzme, H5T_NATIVE_DOUBLE, &
         tempvz_me(1:n2m,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

!RO   Close properties and datasets and file

      call h5dclose_f(dset_vxme, hdf_error)
      call h5dclose_f(dset_vyme, hdf_error)
      call h5dclose_f(dset_vzme, hdf_error)
      call h5dclose_f(dset_tempme, hdf_error)

      call h5dclose_f(dset_vxrms, hdf_error)
      call h5dclose_f(dset_vyrms, hdf_error)
      call h5dclose_f(dset_vzrms, hdf_error)
      call h5dclose_f(dset_temprms, hdf_error)

      call h5dclose_f(dset_kindiss, hdf_error)
      call h5dclose_f(dset_thediss, hdf_error)

      call h5dclose_f(dset_tempvzme, hdf_error)

      call h5sclose_f(memspace, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)


!RO   Read the grid & statistics information

      ndims=1
      dims_grid(1)=1

      call h5fopen_f(filnamgrid, H5F_ACC_RDONLY_F, file_id, hdf_error, &
                       access_prp=plist_full)

      call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)

      call h5dopen_f(file_id, 'averaging_time', dset_grid, hdf_error)

      call h5dread_f(dset_grid, H5T_NATIVE_INTEGER, timeint_cdsp, &
             dims_grid,hdf_error)

      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      call h5pclose_f(plist_full, hdf_error)
      call h5fclose_f(file_id, hdf_error)

!EP starea.eq.0
      else
      
      timeint_cdsp = 0

!EP   Initialize to 0

      do k=kstart,kend
!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(j)
      do j=1,n2m
       vx_me(j,k)    =0.0d0
       vy_me(j,k)    =0.0d0
       vz_me(j,k)    =0.0d0
       temp_me(j,k)  =0.0d0
       vx_rms(j,k)   =0.0d0
       vy_rms(j,k)   =0.0d0
       vz_rms(j,k)   =0.0d0
       temp_rms(j,k) =0.0d0
       disste(j,k) = 0.0d0
       dissth(j,k) = 0.0d0
       tempvz_me(j,k) = 0.0d0
      enddo
!$OMP  END PARALLEL DO 
      enddo

      endif

      return
      end
