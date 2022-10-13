      subroutine Dump3DFiles
      use param
      use mpih
      use mpi_param, only: kstart,kend
      use local_arrays, only: pr,vy,vz,vx,temp
      use hdf5
      implicit none

      integer hdf_error
      integer :: i,j,k,ipp,jpp,kpp,kmm,jmm,imm
      integer :: ip,jp,kp,im,jm,km

      integer(HID_T) :: file_id
      integer(HID_T) :: filespace
      integer(HID_T) :: slabspace
      integer(HID_T) :: memspace

      integer(HID_T) :: dset_vx
      integer(HID_T) :: dset_vy
      integer(HID_T) :: dset_vz
      integer(HID_T) :: dset_ke
      integer(HID_T) :: dset_pr

      integer(HSIZE_T) :: dims(3)

      integer(HID_T) :: plist_id
      integer(HSIZE_T), dimension(3) :: data_count  
      integer(HSSIZE_T), dimension(3) :: data_offset 

      integer :: comm, info
      integer :: ndims
      integer :: mydata
      integer :: my_down, my_up,tag

      real :: tprfi,h11,h22,h33
      real :: h12,h13,h21,h23,h31,h32,grad_ke
      real :: vxcc,vycc,vzcc,fack,facj
      real :: udx1_c, udx1_m, udx2_c, udx2_m, udx3_c, udx3_m
      real, dimension(n2) :: g2rm_t, g2rm_t2
      real, dimension(n3) :: g3rm_t, g3rm_t2
      integer :: itime
      character*30 filnam1,filnam2,filnam3
      character*30 filnam4,filnam5,filnam6
      character*30 filnam7,filnam8,filnam9
      character*30 filnam0
      character*5 ipfi
 
!RO   Form name

      tprfi = 100.
      itime=nint(time*tprfi)
      write(ipfi,82)itime
   82 format(i5.5)

!RO   Sort out MPI definitions

      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

!RO   Form the name of the file

      filnam1 = 'movie/pr'//ipfi//'.h5'
      filnam2 = 'movie/vx_'//ipfi//'.h5'
      filnam3 = 'movie/vy_'//ipfi//'.h5'
      filnam4 = 'movie/vz_'//ipfi//'.h5'
      filnam5 = 'movie/dn_'//ipfi//'.h5'

!RO   Set offsets and element counts
   
      ndims = 3

      dims(1)=n1
      dims(2)=n2
      dims(3)=n3m

      call h5screate_simple_f(ndims, dims, filespace, hdf_error)

      data_count(1) = n1
      data_count(2) = n2
      data_count(3) = kend-kstart+1

      data_offset(1) = 0
      data_offset(2) = 0
      data_offset(3) = kstart-1

!    Dens

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id,hdf_error)

      call h5pset_fapl_mpio_f(plist_id, comm, info,hdf_error)

      call h5fcreate_f(filnam5, H5F_ACC_TRUNC_F, file_id,&
       hdf_error, access_prp=plist_id)

      call h5pclose_f(plist_id, hdf_error)

      call h5dcreate_f(file_id, 'temp', H5T_NATIVE_DOUBLE,&
                      filespace, dset_ke, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

      call h5dget_space_f(dset_ke, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F,&
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,&
                              hdf_error)
       call h5dwrite_f(dset_ke, H5T_NATIVE_DOUBLE,&
         temp(1:n1,1:n2,kstart:kend), dims, &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5dclose_f(dset_ke, hdf_error)

      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

!    Vx

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id,hdf_error)

      call h5pset_fapl_mpio_f(plist_id, comm, info,hdf_error)

      call h5fcreate_f(filnam2, H5F_ACC_TRUNC_F, file_id, hdf_error, access_prp=plist_id)

      call h5pclose_f(plist_id, hdf_error)

      call h5dcreate_f(file_id, 'Vx', H5T_NATIVE_DOUBLE, &
                      filespace,dset_ke, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

      call h5dget_space_f(dset_ke, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,hdf_error)

       call h5dwrite_f(dset_ke, H5T_NATIVE_DOUBLE, &
         vx(1:n1,1:n2,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5dclose_f(dset_ke, hdf_error)

      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

!    Vy

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)

      call h5pset_fapl_mpio_f(plist_id, comm, info,hdf_error)

      call h5fcreate_f(filnam3, H5F_ACC_TRUNC_F, file_id,hdf_error, access_prp=plist_id)

      call h5pclose_f(plist_id, hdf_error)

      call h5dcreate_f(file_id, 'Vy', H5T_NATIVE_DOUBLE,filespace, dset_ke, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

      call h5dget_space_f(dset_ke, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)

      call h5dwrite_f(dset_ke, H5T_NATIVE_DOUBLE, &
         vy(1:n1,1:n2,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5dclose_f(dset_ke, hdf_error)

      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

!    Vz

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id,hdf_error)

      call h5pset_fapl_mpio_f(plist_id, comm, info,hdf_error)

      call h5fcreate_f(filnam4, H5F_ACC_TRUNC_F, file_id,hdf_error, access_prp=plist_id)

      call h5pclose_f(plist_id, hdf_error)

      call h5dcreate_f(file_id, 'Vz', H5T_NATIVE_DOUBLE,filespace,dset_ke, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

      call h5dget_space_f(dset_ke, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,hdf_error)
 
      call h5dwrite_f(dset_ke, H5T_NATIVE_DOUBLE, &
         vz(1:n1,1:n2,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5dclose_f(dset_ke, hdf_error)

      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

      end subroutine Dump3DFiles
      
