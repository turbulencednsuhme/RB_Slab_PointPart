c===============================================
      subroutine block(n, p, irank, istart, iend, blcsz)
      implicit none
      integer,intent(in) :: n,p,irank
      integer,intent(out) :: istart,iend
      integer :: i
      integer,dimension(0:p-1),intent(out) :: blcsz
      
      do i=0,p-1
      blcsz(i) = floor(real((n+p-i-1)/p))
      enddo
      istart = sum(blcsz(0:irank))-blcsz(irank)+1
      iend = istart+blcsz(irank)-1

      end subroutine block
c=================================================           
      subroutine mpi_workdistribution
      use param
      use mpih 
      use mpi_param
      implicit none
      integer :: i
      
      if(.not. allocated(countj)) allocate(countj(0:numtasks-1))
      if(.not. allocated(countjp)) allocate(countjp(0:numtasks-1))
      if(.not. allocated(countk)) allocate(countk(0:numtasks-1))

!EP   For PERIODIC pressure solver
      call block(n2+1, numtasks, myid, jstartp, jendp, countjp)
      djp=jendp-jstartp+1

      call block(n2m, numtasks, myid, jstart, jend, countj)
      dj=jend-jstart+1
      
      call block(n3m, numtasks, myid, kstart, kend, countk)
      dk=kend-kstart+1

#ifdef DEBUG
      write(*,*) "jstart: ",jstart
      write(*,*) "jend: ",jend
      write(*,*) "jstartp: ",jstart
      write(*,*) "jendp: ",jend
      write(*,*) "kstart: ",kstart
      write(*,*) "kend: ",kend
#endif

      if( dj .lt. 1 ) then            
       write(6,*)'process ',myid,' has work load <1 cell in j direction'
       write(6,*)"Check grid dimensions and number of processes"
       
       call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif

      if( dk .lt. 1 ) then            
       write(6,*)'process ',myid,' has work load <1 cell in k direction'
       write(6,*)"Check grid dimensions and number of processes"
       
       call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif
  
      if(.not. allocated(offsetjp)) allocate(offsetjp(0:numtasks-1))
      if(.not. allocated(offsetj)) allocate(offsetj(0:numtasks-1))
      if(.not. allocated(offsetk)) allocate(offsetk(0:numtasks-1))
      
      offsetjp(:)=0
      offsetj(:)=0
      offsetk(:)=0
      do i=1,numtasks-1
        offsetjp(i)= offsetjp(i-1) + countjp(i-1)
        offsetj(i)= offsetj(i-1) + countj(i-1)
        offsetk(i)= offsetk(i-1) + countk(i-1)
      end do
      
      !-------For MPI-IO--------------------------------
      mydata= n2*dk*n1
      mydatam = n2m*dk*n1m

      if(myid .eq. numtasks-1) mydata = n2*(dk+1)*n1
      
      if(.not. allocated(countf)) allocate(countf(0:numtasks-1))
      if(.not. allocated(offsetf)) allocate(offsetf(0:numtasks-1))
       
      call MPI_ALLGATHER(mydata, 1, MPI_INTEGER, countf, 1, MPI_INTEGER,
     & MPI_COMM_WORLD,ierr)
    
      offsetf(:)=0
      do i=1,numtasks-1
        offsetf(i)= offsetf(i-1) + countf(i-1)
      end do
      
      !------------------------------------------------
      
      end subroutine mpi_workdistribution 
c===============================================
      subroutine update_both_ghosts(n1,n2,vx,ks,ke)
      use mpih
      implicit none
      integer, intent(in) :: ks,ke
      real,intent(inout) :: vx(n1,n2,ks-1:ke+1)
      integer,intent(in) :: n1,n2
      integer :: mydata
      integer :: my_down, my_up,tag
      
      mydata= n1*n2
      
      my_down=myid-1
      
      my_up=myid+1

      if(myid .eq. 0) my_down=MPI_PROC_NULL
      if(myid .eq. numtasks-1) my_up=MPI_PROC_NULL

      tag=1
      call MPI_ISEND(vx(1,1,ke), mydata, MDP,
     & my_up,tag,MPI_COMM_WORLD,req(1),ierr)
      
      call MPI_ISEND(vx(1,1,ks), mydata,  MDP,
     & my_down,tag,MPI_COMM_WORLD,req(2), ierr)
     
      call MPI_IRECV(vx(1,1,ks-1), mydata,  MDP, 
     & my_down,tag,MPI_COMM_WORLD,req(3),ierr)
     
      call MPI_IRECV(vx(1,1,ke+1), mydata,  MDP,
     & my_up, tag,MPI_COMM_WORLD,req(4),ierr)
     
      call MPI_Waitall(4,req,status,ierr)

      end subroutine update_both_ghosts
c=========================================
      subroutine update_upper_ghost(n1,n2,vx)
      use mpih
      use mpi_param, only: kstart,kend,dk
      implicit none
      real,intent(inout) :: vx(n1,n2,kstart-1:kend+1)
      integer,intent(in) :: n1,n2
      integer :: mydata
      integer :: my_down, my_up,tag
       
      mydata= n1*n2
      
      my_down= myid-1
      
      my_up= myid+1

      if(myid .eq. 0) my_down= MPI_PROC_NULL
      if(myid .eq. numtasks-1) my_up= MPI_PROC_NULL
     
      tag=1
      
      call MPI_ISEND(vx(1,1,kstart), mydata, MDP,
     & my_down, tag, MPI_COMM_WORLD, req(1), ierr)
      
      call MPI_IRECV(vx(1,1,kend+1), mydata, MDP,
     & my_up,tag, MPI_COMM_WORLD, req(2), ierr)
       
      call MPI_Waitall(2,req,status,ierr)
     
      end subroutine update_upper_ghost
c=========================================
      subroutine update_lower_ghost(n1,n2,vx)
      use mpih
      use mpi_param, only: kstart,kend,dk
      implicit none
      real,intent(inout) :: vx(n1,n2,kstart-1:kend+1)
      integer,intent(in) :: n1,n2
      integer :: mydata
      integer :: my_down, my_up,tag
       
      mydata= n1*n2
      
      my_down= myid-1
      
      my_up= myid+1

      if(myid .eq. 0) my_down= MPI_PROC_NULL
      if(myid .eq. numtasks-1) my_up= MPI_PROC_NULL
      
      
      tag=1
      
      call MPI_ISEND(vx(1,1,kend), mydata,  MDP,
     & my_up, tag, MPI_COMM_WORLD, req(1), ierr)
      
      call MPI_IRECV(vx(1,1,kstart-1), mydata,  MDP,
     & my_down,tag, MPI_COMM_WORLD, req(2), ierr)
       
      call MPI_Waitall(2,req,status,ierr)
     
      end subroutine update_lower_ghost
c================================================
      subroutine mpi_write_continua
      use param
      use mpih
      use mpi_param, only: kstart,kend
      use local_arrays, only: temp,vy,vz,vx
      use hdf5
      implicit none

      integer hdf_error

      integer(HID_T) :: file_id
      integer(HID_T) :: filespace
      integer(HID_T) :: slabspace
      integer(HID_T) :: memspace

      integer(HID_T) :: dset_vx
      integer(HID_T) :: dset_vy
      integer(HID_T) :: dset_vz
      integer(HID_T) :: dset_temp

      integer(HSIZE_T) :: dims(3)

      integer(HID_T) :: plist_id
      integer(HSIZE_T), dimension(3) :: data_count  
      integer(HSSIZE_T), dimension(3) :: data_offset 

      integer :: comm, info
      integer :: ndims

      character*30 filnam1,filnam2,filnam3,filnam4

!RO   Sort out MPI definitions

      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

!RO   Form the name of the file

      filnam1 = 'continua_temp.h5'
      filnam2 = 'continua_vx.h5'
      filnam3 = 'continua_vy.h5'
      filnam4 = 'continua_vz.h5'

!RO   Set offsets and element counts
   
      ndims = 3

      dims(1)=n1
      dims(2)=n2
      dims(3)=n3m

      call h5screate_simple_f(ndims, dims, 
     &                        filespace, hdf_error)

      data_count(1) = n1
      data_count(2) = n2
      data_count(3) = kend-kstart+1

      data_offset(1) = 0
      data_offset(2) = 0
      data_offset(3) = kstart-1

!EP   temp

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id,
     %    hdf_error)

      call h5pset_fapl_mpio_f(plist_id, comm, info,
     &  hdf_error)

      call h5fcreate_f(filnam1, H5F_ACC_TRUNC_F, file_id,
     & hdf_error, access_prp=plist_id)

      call h5pclose_f(plist_id, hdf_error)

      call h5dcreate_f(file_id, 'temp', H5T_NATIVE_DOUBLE,
     &                filespace,
     &                dset_temp, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

      call h5dget_space_f(dset_temp, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F,
     &                      data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,
     &                        hdf_error)
       call h5dwrite_f(dset_temp, H5T_NATIVE_DOUBLE,
     &   temp(1:n1,1:n2,kstart:kend), dims, 
     &   hdf_error, file_space_id = slabspace, mem_space_id = memspace, 
     &   xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5dclose_f(dset_temp, hdf_error)

      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

!EP   vx

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id,
     %    hdf_error)

      call h5pset_fapl_mpio_f(plist_id, comm, info,
     &  hdf_error)

      call h5fcreate_f(filnam2, H5F_ACC_TRUNC_F, file_id,
     & hdf_error, access_prp=plist_id)

      call h5pclose_f(plist_id, hdf_error)

      call h5dcreate_f(file_id, 'Vx', H5T_NATIVE_DOUBLE,
     &                filespace,
     &                dset_vx, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

      call h5dget_space_f(dset_vx, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F,
     &                      data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,
     &                        hdf_error)
       call h5dwrite_f(dset_vx, H5T_NATIVE_DOUBLE,
     &   vx(1:n1,1:n2,kstart:kend), dims, 
     &   hdf_error, file_space_id = slabspace, mem_space_id = memspace, 
     &   xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5dclose_f(dset_vx, hdf_error)

      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

!EP   vy

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id,
     %    hdf_error)

      call h5pset_fapl_mpio_f(plist_id, comm, info,
     &  hdf_error)

      call h5fcreate_f(filnam3, H5F_ACC_TRUNC_F, file_id,
     & hdf_error, access_prp=plist_id)

      call h5pclose_f(plist_id, hdf_error)

      call h5dcreate_f(file_id, 'Vy', H5T_NATIVE_DOUBLE,
     &                filespace,
     &                dset_vy, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

      call h5dget_space_f(dset_vy, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F,
     &                      data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,
     &                        hdf_error)
       call h5dwrite_f(dset_vy, H5T_NATIVE_DOUBLE,
     &   vy(1:n1,1:n2,kstart:kend), dims, 
     &   hdf_error, file_space_id = slabspace, mem_space_id = memspace, 
     &   xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5dclose_f(dset_vy, hdf_error)

      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

!EP   vz

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id,
     %    hdf_error)

      call h5pset_fapl_mpio_f(plist_id, comm, info,
     &  hdf_error)

      call h5fcreate_f(filnam4, H5F_ACC_TRUNC_F, file_id,
     & hdf_error, access_prp=plist_id)

      call h5pclose_f(plist_id, hdf_error)

      call h5dcreate_f(file_id, 'Vz', H5T_NATIVE_DOUBLE,
     &                filespace,
     &                dset_vz, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

      call h5dget_space_f(dset_vz, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F,
     &                      data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,
     &                        hdf_error)
       call h5dwrite_f(dset_vz, H5T_NATIVE_DOUBLE,
     &   vz(1:n1,1:n2,kstart:kend), dims, 
     &   hdf_error, file_space_id = slabspace, mem_space_id = memspace, 
     &   xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5dclose_f(dset_vz, hdf_error)

      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)
      
      if (myid .eq. 0) then
       open(13,file='continua_grid.dat',status='unknown')
       rewind(13)                                                      
       write(13,*) n1,n2,n3
       write(13,*) rext,time
       write(13,*) istr3,str3
       close(13)
      endif

      end subroutine mpi_write_continua
      
c================================================
      
      subroutine mpi_read_continua(n1o,n2o,n3o,ks,ke,intvar,qua)
      use mpih
      use param
      use hdf5
      implicit none
      integer, intent(in) :: ks,ke,n2o,n1o,n3o
      real, dimension(1:n1o,1:n2o,ks-1:ke+1)::qua

      integer hdf_error

      integer(HID_T) :: file_id
      integer(HID_T) :: slabspace
      integer(HID_T) :: memspace

      integer(HID_T) :: dset_qua

      integer(HSIZE_T) :: dims(3)

      integer(HID_T) :: plist_id
      integer(HSIZE_T), dimension(3) :: data_count  
      integer(HSSIZE_T), dimension(3) :: data_offset 

      integer :: comm, info
      integer :: ndims

      integer, intent(in) :: intvar
      character*70 :: filnam1
      character*10 :: dsetname

      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

!EP   Select file and dataset based on intvar

      select case (intvar)
        case (1)
          dsetname = trim('Vx')
          filnam1 = trim('continua_vx.h5')
        case (2)
          dsetname = trim('Vy')
          filnam1 = trim('continua_vy.h5')
        case (3)
          dsetname = trim('Vz')
          filnam1 = trim('continua_vz.h5')
        case (4)
          dsetname = trim('temp')
          filnam1 = trim('continua_temp.h5')
      end select

!RO   Set offsets and element counts
   
      ndims = 3

      dims(1)=n1o
      dims(2)=n2o
      dims(3)=n3o-1


      data_count(1) = n1o
      data_count(2) = n2o
      data_count(3) = ke-ks+1

      data_offset(1) = 0
      data_offset(2) = 0
      data_offset(3) = ks-1


      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id,
     %    hdf_error)


      call h5pset_fapl_mpio_f(plist_id, comm, info,
     &  hdf_error)

      call h5fopen_f(filnam1, H5F_ACC_RDONLY_F, file_id,
     & hdf_error, access_prp=plist_id)

      call h5dopen_f(file_id, dsetname,
     &                dset_qua, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

      call h5dget_space_f(dset_qua, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F,
     &                      data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,
     &                        hdf_error)
       call h5dread_f(dset_qua, H5T_NATIVE_DOUBLE,
     &   qua(1:n1o,1:n2o,ks:ke), dims, 
     &   hdf_error, file_space_id = slabspace, mem_space_id = memspace, 
     &   xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5dclose_f(dset_qua, hdf_error)

      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

      end subroutine mpi_read_continua

