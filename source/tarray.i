c
c
c     ################################################################
c     ## COPYRIGHT (C) 2013 by Xiao Zhu, Pengyu Ren & Jay W. Ponder ##
c     ##                     All Rights Reserved                    ##
c     ################################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  tarray.i  --  temporary array storage for saving calc    ##
c     ##                                                           ##
c     ###############################################################
c
c
c     npair        total number of dipole-dipole within cutoff
c     ta_dipdip    temporary array of dipole-dipole interaction field
c     ta_dip_index temporary index array of dipole-dipole 
c
c
      integer npair
      integer, pointer :: ta_dip_index(:,:)
      real*8, pointer :: ta_dipdip(:,:)
      common /tarray/ npair,ta_dipdip,ta_dip_index
