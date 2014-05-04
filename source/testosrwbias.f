c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2014  by  Jay William Ponder  ##
c     ##           Jayvee Abella, Pengyu Ren           ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  program testosrwbias  --  test osrw derivatives           ##
c     ##                                                            ##
c     ################################################################
c
c
c     "testosrwbias" computes the analytical and numerical
c     gradient vectors of the potential energy function with respect
c     to Cartesian coordinates and lambda for a given lambda using osrw
c
c
      program testosrwbias
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bath.i'
      include 'deriv.i'
      include 'energi.i'
      include 'inform.i'
      include 'inter.i'
      include 'iounit.i'
      include 'mutant.i'
      include 'osrwi.i'
      include 'solute.i'
      include 'usage.i'
      integer i,j,next
      real*8 etot,f,f0,eps,eps0,old,energy,DTenergy
      real*8 DTu,DTu0
      real*8 totnorm,ntotnorm,rms,nrms
      real*8 oldl
      real*8 f1, f2, f3
      real*8 ndudl, nd2udl2
      real*8 t0,t1,tlength, told
      real*8 omp_get_wtime
      real*8, allocatable :: denorm(:)
      real*8, allocatable :: ndenorm(:)
      real*8, allocatable :: detot(:,:)
      real*8, allocatable :: ndetot(:,:)
      
      logical exist,query
      logical doanalyt,donumer,dofull,doendpoints
      logical dosweep, dotime
      character*1 answer
      character*1 axis(3)
      character*120 record
      character*120 string
      data axis  / 'X','Y','Z' /
c
c
c     set up the structure and mechanics calculation
c
      call initial
      call getxyz
      call mechanic
      
c
c    only works if there was previous histogram file
c
      if (.not.restartosrw) then
        print *, "Need histogram file"
        call fatal
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (detot(3,n))
      
      
c
c     decide whether to do timing compare
c
      dotime = .false.
      call nextarg (answer,exist)
      if (.not. exist) then
         write (iout,9)
   9    format (/,' Do timing test [N] :  ',$)
         read (input,10)  record
   10    format (a120)
         next = 1
         call gettext (record,answer,next)
      end if
      call upcase (answer)
      if (answer .eq. 'Y')  dotime = .true. 
      
      if (dotime) then
        lambda = 0.8d0
        print *, ""
        print *, "Lambda =", lambda
        print *, "Running 100 iterations of osrw"
        propagateLambda = .false.
        osrwon = .true.
        t0=omp_get_wtime()
        do i = 1,100
          call osrw
        end do
        t1=omp_get_wtime()
        tlength = t1 - t0
        told = tlength
        print *, "osrw 100x:", tlength
        print *, "Running 100 iterations of normal TINKER"
        osrwon = .false.
        lambda = 1.0d0
        t0=omp_get_wtime()
        do i = 1,100
          call gradient (etot, detot)
        end do
        t1=omp_get_wtime()
        tlength = t1 - t0
        print *, "normal 100x:", tlength
        tlength = told / tlength
        print *, "Code is slower by:", tlength
        call fatal
      end if    
      
c
c     decide whether to do an analytical gradient calculation
c
      doanalyt = .true.
      call nextarg (answer,exist)
      if (.not. exist) then
         write (iout,15)
   15    format (/,' Compute the Analytical Gradient Vector [Y] :  ',$)
         read (input,20)  record
   20    format (a120)
         next = 1
         call gettext (record,answer,next)
      end if
      call upcase (answer)
      if (answer .eq. 'N')  doanalyt = .false.
c
c     decide whether to do a numerical gradient calculation
c
      donumer = .true.
      call nextarg (answer,exist)
      if (.not. exist) then
         write (iout,30)
   30    format (/,' Compute the Numerical Gradient Vector [Y] :   ',$)
         read (input,40)  record
   40    format (a120)
         next = 1
         call gettext (record,answer,next)
      end if
      call upcase (answer)
      if (answer .eq. 'N')  donumer = .false.
c
c     get the stepsize for numerical gradient calculation
c
      if (donumer) then
         eps = -1.0d0
         eps0 = 1.0d-5
         query = .true.
         call nextarg (string,exist)
         if (exist) then
            read (string,*,err=50,end=50)  eps
            query = .false.
         end if
   50    continue
         if (query) then
            write (iout,60)  eps0
   60       format (/,' Enter a Numerical Stepsize [',d7.1,
     &                 ' Ang] :  ',$)
            read (input,70,err=50)  eps
   70       format (f20.0)
         end if
         if (eps .le. 0.0d0)  eps = eps0
      end if
c
c     decide whether to output results by gradient component
c
      dofull = .true.
      if (n .gt. 100) then
         dofull = .false.
         call nextarg (answer,exist)
         if (.not. exist) then
            write (iout,80)
   80       format (/,' Output Breakdown by Gradient Component',
     &                 ' [N] :  ',$)
            read (input,90)  record
   90       format (a120)
            next = 1
            call gettext (record,answer,next)
         end if
         call upcase (answer)
         if (answer .eq. 'Y')  dofull = .true.
      end if

        

c
c     compute the analytical gradient components
c
      propagateLambda = .false.
      if (doanalyt) then
c         call gradient2 (etot,detot)
         call osrw
         DTu0 = esum
      end if
c
c     print the total dedl of the system
c
      if (doanalyt) then
         if (digits .ge. 8) then
            write (iout,100)  dudl
  100       format (/,' Total dedl :',8x,f20.8,' Kcal/mole')
         else if (digits .ge. 6) then
            write (iout,110)  dudl
  110       format (/,' Total dedl :',8x,f18.6,' Kcal/mole')
         else
            write (iout,120)  dudl
  120       format (/,' Total dedl :',8x,f16.4,' Kcal/mole')
         end if
c
c     print the lambda related variables over individual components
c
     
         print *, 'Lambda related terms'
         print *, ''
         print *, 'lam',lambda
         print *, 'e',DTu0
         print *, 'dudl',dudl
         if (lambda .eq. 1.0d0 .or. lambda .eq. 0.0d0) then
           print *, "Derivatives undefined at lambda endpoints"
           print *, "Make sure 0 < lambda < 1"
           call fatal
         end if
         
      end if
      

c     print a header for the gradients of individual potentials
c
      if (dofull) then
         write (iout,200)
  200    format (/,' Cartesian Gradient Breakdown by Individual',
     &                 ' Components :')
         
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (ndetot(3,n))
      
c     check numerical dudl, d2udl2
      if (donumer) then
        oldl = lambda
        
        lambda = lambda - 0.5d0*eps
        call osrw
        DTu0 = esum
        
        lambda = lambda + eps
        call osrw
        DTu = esum
        ndudl = (DTu - DTu0) / eps
              
        lambda = oldl      
              
      end if
      
      if (doanalyt) then
        print *, 'Analytic'
        call osrw
        DTu0 = esum
        print *, 'dudl',dudl
c        print *, 'd2udl2',d2udl2
        print *, ''
      end if
      if (donumer) then
        print *, 'Numeric'
        print *, 'dudl',ndudl
c        print *, 'd2udl2',nd2udl2
        print *, ''
      end if
      
c
c     get the Cartesian component two-sided numerical gradients

      do i = 1, n
         if (donumer .and. use(i)) then
            do j = 1, 3
c      numerical gradient calculation (spatial)
               if (j .eq. 1) then
                  old = x(i)
                  x(i) = x(i) - 0.5d0*eps
               else if (j .eq. 2) then
                  old = y(i)
                  y(i) = y(i) - 0.5d0*eps
               else if (j .eq. 3) then
                  old = z(i)
                  z(i) = z(i) - 0.5d0*eps
               end if
c              DTenergy test
               call osrw
               DTu0 = esum
               
               if (j .eq. 1) then
                  x(i) = x(i) + eps
               else if (j .eq. 2) then
                  y(i) = y(i) + eps
               else if (j .eq. 3) then
                  z(i) = z(i) + eps
               end if

c             DTenergy test
               call osrw
               DTu = esum
               ndetot(j,i) = (DTu - DTu0) / eps
               
               if (j .eq. 1) then
                  x(i) = old
               else if (j .eq. 2) then
                  y(i) = old
               else if (j .eq. 3) then
                  z(i) = old
               end if
               
            end do
         end if
c
c     print analytical gradients of each energy term for each atom

         if (dofull .and. use(i)) then
            call osrw
            DTu0 = esum
            do j = 1, 3
               if (doanalyt) then
                  print *, 'Analytic',i,axis(j)
                  print *, 'desum', desum(j,i)
c                  print *, 'd2udl(r)',d2udl(j,i)
                  print *, ''
               end if
c
c     print numerical gradients of each energy term for each atom
c
               if (donumer) then
                  print *, 'Numeric',i,axis(j)
                  print *, 'detot', ndetot(j,i)
c                  print *, 'd2udl(r)',nd2udl(j,i)
                  print *, ''
               end if
            end do
         end if
      end do

      
c
c     perform dynamic allocation of some local arrays
c
      allocate (denorm(n))
      allocate (ndenorm(n))
c
c     print the total gradient components for each atom
c
      write (iout,300)
  300 format (/,' gradient of dE/dl :')
      if (digits .ge. 8) then
         write (iout,310)
  310    format (/,2x,'Type',4x,'Atom',10x,'dE/dXdl',11x,'dE/dYdl',
     &              11x,'dE/dZdl',11x,'Norm',/)
      else if (digits .ge. 6) then
         write (iout,320)
  320    format (/,2x,'Type',6x,'Atom',11x,'dE/dXdl',9x,'dE/dYdl',
     &              9x,'dE/dZdl',11x,'Norm',/)
      else
         write (iout,330)
  330    format (/,2x,'Type',6x,'Atom',14x,'dE/dXdl',7x,'dE/dYdl',
     &              7x,'dE/dZdl',10x,'Norm',/)
      end if
      totnorm = 0.0d0
      ntotnorm = 0.0d0
      do i = 1, n
         if (doanalyt .and. use(i)) then
            f1 = desum(1,i)
            f2 = desum(2,i)
            f3 = desum(3,i)
            denorm(i) = f1**2 + f2**2 +
     &                       f3**2
            totnorm = totnorm + denorm(i)
            denorm(i) = sqrt(denorm(i))
            if (digits .ge. 8) then
               write (iout,340)  i,f1,f2,f3,denorm(i)
  340          format (' Anlyt',i8,1x,3f16.8,f16.8)
            else if (digits .ge. 6) then
               write (iout,350)  i,f1,f2,f3,denorm(i)
  350          format (' Anlyt',2x,i8,3x,3f14.6,2x,f14.6)
            else
               write (iout,360)  i,f1,f2,f3,denorm(i)
  360          format (' Anlyt',2x,i8,7x,3f12.4,2x,f12.4)
            end if
         end if
         if (donumer .and. use(i)) then
            ndenorm(i) = ndetot(1,i)**2 + ndetot(2,i)**2 +
     &                        ndetot(3,i)**2
            ntotnorm = ntotnorm + ndenorm(i)
            ndenorm(i) = sqrt(ndenorm(i))
            if (digits .ge. 8) then
               write (iout,370)  i,(ndetot(j,i),j=1,3),ndenorm(i)
  370          format (' Numer',i8,1x,3f16.8,f16.8)
            else if (digits .ge. 6) then
               write (iout,380)  i,(ndetot(j,i),j=1,3),ndenorm(i)
  380          format (' Numer',2x,i8,3x,3f14.6,2x,f14.6)
            else
               write (iout,390)  i,(ndetot(j,i),j=1,3),ndenorm(i)
  390          format (' Numer',2x,i8,7x,3f12.4,2x,f12.4)
            end if
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (detot)
      deallocate (ndetot)
      deallocate (denorm)
      deallocate (ndenorm)
c
c     print the total norm for the analytical gradient
c
      write (iout,400)
  400 format (/,' Total Gradient Norm and RMS Gradient per Atom :')
      if (doanalyt) then
         totnorm = sqrt(totnorm)
         if (digits .ge. 8) then
            write (iout,410)  totnorm
  410       format (/,' Anlyt',6x,'Total Gradient Norm Value',6x,f20.8)
         else if (digits .ge. 6) then
            write (iout,420)  totnorm
  420       format (/,' Anlyt',6x,'Total Gradient Norm Value',6x,f18.6)
         else
            write (iout,430)  totnorm
  430       format (/,' Anlyt',6x,'Total Gradient Norm Value',6x,f16.4)
         end if
      end if
c
c     print the total norm for the numerical gradient
c
      if (donumer) then
         ntotnorm = sqrt(ntotnorm)
         if (digits .ge. 8) then
            write (iout,440)  ntotnorm
  440       format (' Numer',6x,'Total Gradient Norm Value',6x,f20.8)
         else if (digits .ge. 6) then
            write (iout,450)  ntotnorm
  450       format (' Numer',6x,'Total Gradient Norm Value',6x,f18.6)
         else
            write (iout,460)  ntotnorm
  460       format (' Numer',6x,'Total Gradient Norm Value',6x,f16.4)
         end if
      end if
c
c     print the rms per atom norm for the analytical gradient
c
      if (doanalyt) then
         rms = totnorm / sqrt(dble(nuse))
         if (digits .ge. 8) then
            write (iout,470)  rms
  470       format (/,' Anlyt',6x,'RMS Gradient over All Atoms',
     &                 4x,f20.8)
         else if (digits .ge. 6) then
            write (iout,480)  rms
  480       format (/,' Anlyt',6x,'RMS Gradient over All Atoms',
     &                 4x,f18.6)
         else
            write (iout,490)  rms
  490       format (/,' Anlyt',6x,'RMS Gradient over All Atoms',
     &                 4x,f16.4)
         end if
      end if
c
c     print the rms per atom norm for the numerical gradient
c
      if (donumer) then
         nrms = ntotnorm / sqrt(dble(nuse))
         if (digits .ge. 8) then
            write (iout,500)  nrms
  500       format (' Numer',6x,'RMS Gradient over All Atoms',4x,f20.8)
         else if (digits .ge. 6) then
            write (iout,510)  nrms
  510       format (' Numer',6x,'RMS Gradient over All Atoms',4x,f18.6)
         else
            write (iout,520)  nrms
  520       format (' Numer',6x,'RMS Gradient over All Atoms',4x,f16.4)
         end if
      end if
c
c     perform any final tasks before program exit
c
      call final
      end

