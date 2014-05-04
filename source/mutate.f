c
c
c     ##############################################################
c     ##  COPYRIGHT (C) 2009 by Chuanjie Wu & Jay William Ponder  ##
c     ##                   All Rights Reserved                    ##
c     ##############################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine mutate  --  set parameters for hybrid system  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "mutate" constructs the hybrid hamiltonian for a specified
c     initial state, final state and mutation parameter "lambda"
c
c

      subroutine mutate
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'inform.i'
      include 'iounit.i'
      include 'katoms.i'
      include 'keys.i'
      include 'mutant.i'
      include 'osrwi.i'
      integer i,j,k,ihyb
      integer it0,it1,next
      integer list(20)
      character*20 keyword
      character*120 record
      character*120 string
c
c
c     set defaults for lambda and soft core vdw parameters
c
      lambda = 1.0d0
      vlambda = 1.0d0
      elambda = 1.0d0
      scexp = 5.0d0
      scalpha = 0.7d0
      
c
c     JRA more variables
c      
      scexpv = 1.0d0
      scalphav = 0.05d0
      scexpm = 1.0d0
      scalpham = 2.0d0
      scexpp = 2.0d0
      polstart = 0.75d0
      polend = 1.0d0
      osrwon = .false.
      doNoLigandCondensed = .false.
      doVaporElec = .false.
      use_soft = .false.
      isrelative = .false.
      rel_ligA = .false.
      debugjv = .false.
      envenvon = .true.
      dorealsplit = .false.
      restartosrw = .false.
      
c
c     JRA osrw variables
c
      tempsystem = 298.15d0
      dtosrw = 0.001d0
      propagateLambda = .false.
      biasMag = 0.002d0
      lambdaBins = 201
      FLambdaBins = 401
      dL = 1.0d0 / (lambdaBins - 1)
      dL_2 = dL / 2.0d0
      dFL = 2.0
      dFL_2 = dFL / 2.0d0
      countInterval = 10
      energyCount = 0
      biasCutoff = 5
      totalEnergy = 0.0d0
      totalFreeEnergy = 0.0d0
      minFLambda = -(dFL*FLambdaBins) / 2.0d0
      minLambda = -dL_2
      maxFLambda = minFLambda + FLambdaBins*dFL
      fLambdaUpdates = 0
      fLambdaPrintInterval = 100
      printFrequency = 100
      halfThetaVelocity = 0.0d0
      theta = 0.0d0
      totalCounts = 0
      thetaFriction = 1.0E-19
      thetaMass = 1.0E-18
      hisverbose = .false.
      
      allocate (recursionKernel(lambdaBins,FLambdaBins))
      do i = 1, lambdaBins
        FLambda(i) = 0.0d0
        do j = 1, FLambdaBins
          recursionKernel(i,j) = 0
        end do
      end do
      
c
c     zero number of hybrid atoms and hybrid atom list
c
      nmut = 0
      nmut1 = 0
      nmut2 = 0
      do i = 1, n
         mut(i) = .false.
c         imut(i) = 0
c         imut1(i) = 0
c         imut2(i) = 0
      end do
      do i = 1, 20
         list(i) = 0
      end do
c
c     search keywords for free energy perturbation options
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:7) .eq. 'LAMBDA ') then
            string = record(next:120)
            read (string,*,err=20)  lambda
            vlambda = lambda
            elambda = lambda
            theta = asin(sqrt(lambda))
         else if (keyword(1:11) .eq. 'VDW-LAMBDA ') then
            string = record(next:120)
            read (string,*,err=20)  vlambda
         else if (keyword(1:11) .eq. 'ELE-LAMBDA ') then
            string = record(next:120) 
            read (string,*,err=20)  elambda
         else if (keyword(1:7) .eq. 'MUTATE ') then
            string = record(next:120)
            read (string,*,err=20)  ihyb,it0,it1
            nmut = nmut + 1
            imut(nmut) = ihyb
            mut(ihyb) = .true.
            type0(nmut) = it0
            type1(nmut) = it1
            class0(nmut) = atmcls(it0)
            class1(nmut) = atmcls(it1)
c
c       JRA - changed way in which LIGAND is read
c             (should be minor, includes printing)
c
         else if (keyword(1:7) .eq. 'LIGAND ') then
            string = record(next:120)
            read (string,*,err=10,end=10)  (list(k),k=1,20)
   10       continue
            k = 1
            do while (list(k) .ne. 0) 
               if (list(k) .gt. 0) then
                  j = list(k)
                  nmut = nmut + 1
                  imut(nmut) = j
                  mut(j) = .true.
                  type0(nmut) = 0
                  type1(nmut) = type(j)
                  class0(nmut) = 0
                  class1(nmut) = class(j)
                  k = k + 1
                  write(iout,15) j,name(j)
   15             format(1x,'Atoms in Ligand: ',i5,a5)
               else
                  do j = abs(list(k)), abs(list(k+1))
                     nmut = nmut + 1
                     imut(nmut) = j
                     mut(j) = .true.
                     type0(nmut) = 0
                     type1(nmut) = type(j)
                     class0(nmut) = 0
                     class1(nmut) = class(j)
                     write(iout,15) j,name(j)
                  end do
                  k = k + 2
               end if
            end do
c
c        JRA additional accepted keywords
c
         else if (keyword(1:8) .eq. 'LIGAND1 ') then
            string = record(next:120)
            read (string,*,err=16,end=16)  (list(k),k=1,20)
   16       continue
            k = 1
            do while (list(k) .ne. 0) 
               if (list(k) .gt. 0) then
                  j = list(k)
                  nmut1 = nmut1 + 1
                  imut1(nmut1) = j
                  mut(j) = .true.
                  type0(nmut1) = 0
                  type1(nmut1) = type(j)
                  class0(nmut1) = 0
                  class1(nmut1) = class(j)
                  k = k + 1
                  write(iout,17) j,name(j)
   17             format(1x,'Atoms in Ligand1: ',i5,a5)
               else
                  do j = abs(list(k)), abs(list(k+1))
                     nmut1 = nmut1 + 1
                     imut1(nmut1) = j
                     mut(j) = .true.
                     type0(nmut1) = 0
                     type1(nmut1) = type(j)
                     class0(nmut1) = 0
                     class1(nmut1) = class(j)
                     write(iout,17) j,name(j)
                  end do
                  k = k + 2
               end if
            end do
            nmut = nmut + nmut1
          else if (keyword(1:8) .eq. 'LIGAND2 ') then
            string = record(next:120)
            read (string,*,err=18,end=18)  (list(k),k=1,20)
   18       continue
            k = 1
            do while (list(k) .ne. 0) 
               if (list(k) .gt. 0) then
                  j = list(k)
                  nmut2 = nmut2 + 1
                  imut2(nmut2) = j
                  mut(j) = .true.
                  type0(nmut2) = 0
                  type1(nmut2) = type(j)
                  class0(nmut2) = 0
                  class1(nmut2) = class(j)
                  k = k + 1
                  write(iout,19) j,name(j)
   19             format(1x,'Atoms in Ligand2: ',i5,a5)
               else
                  do j = abs(list(k)), abs(list(k+1))
                     nmut2 = nmut2 + 1
                     imut2(nmut2) = j
                     mut(j) = .true.
                     type0(nmut2) = 0
                     type1(nmut2) = type(j)
                     class0(nmut2) = 0
                     class1(nmut2) = class(j)
                     write(iout,19) j,name(j)
                  end do
                  k = k + 2
               end if
            end do
            nmut = nmut + nmut2
            isrelative = .true.
         else if(keyword(1:14) .eq. 'OSRW-ABSOLUTE ') then
            osrwon = .true.
            isrelative = .false.
         else if(keyword(1:14) .eq. 'OSRW-RELATIVE ') then
            osrwon = .true.
            isrelative = .true.
         else if(keyword(1:12) .eq. 'HIS-VERBOSE ') then
            hisverbose = .true.
         else if (keyword(1:12) .eq. 'VDW-SCALPHA ') then
            string = record(next:120) 
            read (string,*,err=20)  scalphav
         else if (keyword(1:12) .eq. 'ELE-SCALPHA ') then
            string = record(next:120) 
            read (string,*,err=20)  scalpham
         else if (keyword(1:10) .eq. 'VDW-SCEXP ') then
            string = record(next:120) 
            read (string,*,err=20)  scexpv
         else if (keyword(1:10) .eq. 'ELE-SCEXP ') then
            string = record(next:120) 
            read (string,*,err=20)  scexpm
         else if (keyword(1:10) .eq. 'POL-SCEXP ') then
            string = record(next:120) 
            read (string,*,err=20)  scexpp
         else if (keyword(1:9) .eq. 'POLSTART ') then
            string = record(next:120) 
            read (string,*,err=20)  polstart
         else if (keyword(1:7) .eq. 'POLEND ') then
            string = record(next:120) 
            read (string,*,err=20)  polend
         else if(keyword(1:20) .eq. 'DONOLIGANDCONDENSED ') then
            doNoLigandCondensed = .true.
         else if(keyword(1:12) .eq. 'DOVAPORELEC ') then
            doVaporElec = .true.
         else if (keyword(1:8) .eq. 'BIASMAG ') then
            string = record(next:120) 
            read (string,*,err=20)  biasMag
         else if (keyword(1:9) .eq. 'FLAMPRINT ') then
            string = record(next:120) 
            read (string,*,err=20)  fLambdaPrintInterval
         else if (keyword(1:10) .eq. 'THETAMASS ') then
            string = record(next:120) 
            read (string,*,err=20)  thetaMass
         else if (keyword(1:14) .eq. 'THETAFRICTION ') then
            string = record(next:120) 
            read (string,*,err=20)  thetaFriction  
         end if
   20    continue
      end do
      
c
c     JRA use one lambda for osrw simulation
c
      if (elambda .eq. vlambda) then
        lambda = elambda
      end if
      if (osrwon .and. elambda .ne. vlambda) then
        print *, "elambda != vlambda, using elambda"
        lambda = elambda
      end if
      
c
c     JRA when osrwon = .false., code is in original state
c         else, print parameters
c
      if (.not. osrwon) then
        elambda = 1.0d0
        vlambda = 1.0d0
        lambda = 1.0d0
        nmut = 0
      	do i = 1, n
          mut(i) = .false.
      	end do
      	use_soft = .false.
      else
        use_soft = .true.
        print *, "Softcore Parameters"
        print *, "vdw_scalpha", scalphav
        print *, "vdw scexp", scexpv
        print *, ""
        print *, "ele scalpha", scalpham
        print *, "ele scexp", scexpm
        print *, ""
        print *, "pol scexp", scexpp
        print *, "pol start", polstart
        print *, "pol end", polend
        print *, ""
        print *, "DONOLIGANDCONDENSED", doNoLigandCondensed
        print *, "DOVAPORELEC", doVaporElec
        print *, ""
        print *, "biasMag", biasMag
        print *, ""
        
c     JRA - restart osrw if there exists histogram file
        call gethis
      end if
      
c
c     scale electrostatic parameter values based on lambda
c          JRA not used with osrw simulation
      if (.not. osrwon) then
      if (elambda.ge.0.0d0 .and. elambda.lt.1.0d0)  call altelec
      end if
c

c
c     write the status of the current free energy perturbation step
c
      if (nmut.ne.0 .and. .not.silent) then
         write (iout,30)  vlambda
   30    format (/,' Free Energy Perturbation :',f15.3,
     &              ' Lambda for van der Waals')
         write (iout,40)  elambda
   40    format (' Free Energy Perturbation :',f15.3,
     &              ' Lambda for Electrostatics')
      end if
      return
      end
      
c
c
c     ###################################################
c     ##            COPYRIGHT (C)  2013  by            ##
c     ##         Jayvee Abella and Pengyu Ren          ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine gethis  --  get histogram structure             ##
c     ##                                                             ##
c     #################################################################
c
c
c     "gethis" asks for a histogram file to restart OSRW
c
c
      subroutine gethis
      implicit none
      include 'files.i'
      include 'inform.i'
      include 'iounit.i'
      include 'osrwi.i'
      include 'output.i'
      integer ihis
      integer freeunit
      logical exist
      character*120 hisfile
c
c     try to restart using prior recursionKernel
c
      hisfile = filename(1:leng)//'.his'
      call version (hisfile,'old')
      inquire (file=hisfile,exist=exist)
      if (exist) then
         print *, "Restarting from latest histogram file"
         ihis = freeunit ()
         open (unit=ihis,file=hisfile,status='old')
         rewind (unit=ihis)
         call readhis (ihis)
         close (unit=ihis)


      end if
      return
      end      
      
      
c
c
c     ###################################################
c     ##            COPYRIGHT (C)  2013  by            ##
c     ##         Jayvee Abella and Pengyu Ren          ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine readhis  --  input of OSRW history           ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "readhis"  sets recursionKernel to previous history
c
c
      subroutine readhis (ihis)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'files.i'
      include 'inform.i'
      include 'iounit.i'
      include 'mutant.i'
      include 'osrwi.i'
      include 'titles.i'
      integer i,j,k,m
      integer ihis,nmax
      integer next,size
      integer first,last
      integer nexttext
      integer trimtext
      real*8  dummy, updateFLambda
      logical exist,opened
      logical quit,reorder
      logical clash
      character*120 hisfile
      character*401 record
      character*401 string
      character*3   lc
      character*25  fstr

      deallocate(recursionKernel)
      
      quit = .true.
      read(ihis,*) lambda, halfThetaVelocity, minFLambda, FLambdaBins
      
      allocate(recursionKernel(lambdaBins,FLambdaBins))
      
      lc = '401'
      fstr = '('//lc//'i4)'
      do i = 1, lambdaBins
c        do k = 1, FLambdaBins
c          write (ihis,90) recursionKernel(i,k)
         read (ihis,fstr,err=70,end=70)
     &             (recursionKernel(i,j),j=1,FLambdaBins)
   60    continue
      end do
      quit = .false.
   70 continue  
   
      maxFLambda = minFLambda + dFL*FLambdaBins
   
      elambda = lambda
      vlambda = lambda
   
c
c     an error occurred in reading the histogram file
c
      if (quit) then
         print *, "ERROR: READING IN HISTOGRAM FILE"
         print *, ""
         call fatal
      end if

      energyCount = 0
      do i = 1, lambdaBins
        do j = 1, FLambdaBins
          energyCount = energyCount + recursionKernel(i,j)
        end do
      end do
      
      
c    Reset theta variable
      theta = asin(sqrt(lambda))
      
      print *, "lambda:",lambda
      print *, "halfThetaVelocity:",halfThetaVelocity
      print *, "minFLambda, maxFlambda", minFlambda, maxFLambda
      print *, "FLambdaBins, dFL", FLambdaBins, dFL
      print *, ""
      restartosrw = .true.
      dummy = updateFLambda(.true.)

      return
      end      

c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine altelec  --  mutated electrostatic parameters  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "altelec" constructs the mutated electrostatic parameters
c     based on the lambda mutation parameter "elmd"
c
c
      subroutine altelec
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'charge.i'
      include 'mpole.i'
      include 'mutant.i'
      include 'polar.i'
      include 'potent.i'
      integer i,j,k
c
c
c     set electrostatic parameters for partial charge models
c
      if (use_charge) then
         do i = 1, nion
            if (mut(i)) then
               pchg(i) = pchg(i) * elambda
            end if
         end do
      end if
c
c     set electrostatic parameters for polarizable multipole models
c
      if (use_mpole) then
         do i = 1, npole
            k = ipole(i)
            if (mut(k)) then
               do j = 1, 13
                  pole (j,i) = pole(j,i) * elambda
               end do
            end if
         end do
         do i = 1, npolar
            if (mut(i)) then
               polarity(i) = polarity(i) * elambda
            end if
         end do
      end if
      return
      end
