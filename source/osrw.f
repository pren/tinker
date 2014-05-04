c
c
c     ###################################################
c     ##            COPYRIGHT (C)  2013  by            ##
c     ##         Jayvee Abella and Pengyu Ren          ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine osrw - an implementation of osrw algorithm  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "osrw" calls the subroutines to calculate the potential
c     energy terms and sums up to form the total energy
c     and adds and maintains OSRW biasing gaussians
c     Based off code from Force Field X written in Java
c     Array indicies are calculated as if 0 is first element
c       so whenever arrays are accessed, the indicies are incremented
c
c
      subroutine osrw
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bath.i'
      include 'bound.i'
      include 'cutoff.i'
      include 'deriv.i'
      include 'energi.i'
      include 'ewald.i'
      include 'inform.i'
      include 'iounit.i'
      include 'mutant.i'
      include 'osrwi.i'
      include 'polpot.i'
      include 'potent.i'
      include 'rigid.i'
      include 'vdwpot.i'
      
      integer lambdaBin, FLambdaBin
      integer binForLambda, binForFLambda
      integer lcenter, lcount, FLcenter
      integer iL, iFL, i ,j
      real*8  dte, DTenergy
      real*8  biasEnergy
      real*8  dGdLambda, dGdFLambda
      real*8  ls2, FLs2, deltaL, deltaL2, mirrorFactor
      real*8  deltaFL, deltaFL2, weight, bias
      real*8  dEdU, updateFLambda, currentFreeEnergy
      real*8  freeEnergy
      logical toPrint

c
c     In RESPA, OSRW is propagated with slowly varying terms
c        (Set in respa.f)
c      
      if (propagateLambda) then
        energyCount = energyCount + 1
      end if
      
c      
c     Way to speeden code, not used or tested
c      if (propagateLambda .and. energyCount .gt. 0) then
c            
c            Metadynamics grid counts (every 'countInterval' steps).
c             
c        if (mod(energyCount,countInterval) .eq. 0) then
c          poleps = 0.00000001
c        else 
c          poleps = 0.01
c        end if
c      end if
c      poleps = 0.01
      
      
c
c     Get potential energy of dual topology system
c       (or absolute system)
c      
      dte = DTenergy(lambda)
      biasEnergy = 0.0d0
c
c     Get corresponding bin for current lambda
c       and calculated dudl states
c      
      call checkRecursionKernelSize(dudl)     
      lambdaBin = binForLambda(lambda)
      FLambdaBin = binForFLambda(dudl)
      
      
      
c    
c     Calculate recursion kernel G(L, F_L) and its derivatives with respect
c     to L and F_L.
c              
      dGdLambda = 0.0d0
      dGdFLambda = 0.0d0
      ls2 = (2.0d0 * dL) * (2.0d0 * dL)
      FLs2 = (2.0d0 * dFL) * (2.0d0 * dFL)
      iL = -biasCutoff
      do while (iL .le. biasCutoff)
        lcenter = lambdaBin + iL
        deltaL = lambda - (lcenter * dL)
        deltaL2 = deltaL*deltaL
c
c       mirror conditions for recursion kernel counts
c         (Bias is added outside the lambda boundaries to
c           prevent the lambda particle from being stuck at the 
c           endpoints)
c
        lcount = lcenter
        mirrorFactor = 1.0d0
        if (lcount .eq. 0 .or. lcount .eq. lambdaBins - 1) then
          mirrorFactor = 2.0d0
        else if (lcount .lt. 0) then
          lcount = -lcount
        else if (lcount .gt. lambdaBins - 1) then
c         number of bins past the last bin
          lcount = lcount - (lambdaBins - 1)
c         mirror bin
          lcount = -lcount + (lambdaBins - 1)
        end if
        iFL = -biasCutoff
        do while (iFL .le. biasCutoff)
          FLcenter = FLambdaBin + iFL
c
c         If either of the following FL edge conditions are true, then
c         there are no counts and we continue
c
          if (FLcenter .lt. 0 .or. FLcenter .ge. FLambdaBins) then
            iFL = iFL + 1
            cycle
          end if
          deltaFL = dudl - (minFLambda + FLcenter*dFL + dFL_2)
          deltaFL2 = deltaFL*deltaFL
          weight = mirrorFactor*recursionKernel(lcount+1,FLcenter+1)
          bias = weight*biasMag
     &           *exp(-deltaL2/(2.0d0*ls2))
     &           *exp(-deltaFL2/(2.0d0*FLs2))
          biasEnergy = biasEnergy + bias
          dGdLambda = dGdLambda - (deltaL / ls2 * bias)
          dGdFLambda = dGdFLambda - (deltaFL / FLs2 * bias)
          iFL = iFL + 1
        end do
        iL = iL + 1
      end do
            
c      if (mod(energyCount,10) .eq. 0) then
c      print *, "dudlc", dudl,d2udl2,dte
c      print *, "dGdL",dGdLambda,dGdFLambda,biasEnergy
c      end if
      
c
c     Lambda gradient due to recursion kernel G(L,F_L)
c        (and save unbiased dudl)
c      
      dEdU = dudl      
      dudl = dudl + dGdLambda + dGdFLambda*d2udl2
      
c
c     Cartesian coordinate gradient due to recursion kernel G(L,F_L)
c
      do i = 1,n
        do j = 1,3
          desum(j,i) = desum(j,i) + dGdFLambda*d2udl(j,i)
        end do
      end do
      
      if (propagatelambda .and. energyCount .gt. 0) then
c
c       update free energy F(L) every ~10 steps
c
        if (mod(energyCount,10) .eq. 0) then
          fLambdaUpdates = fLambdaUpdates + 1
          toPrint = mod(fLambdaUpdates,fLambdaPrintInterval) .eq. 0
          totalFreeEnergy = updateFLambda(toPrint)
        end if
c
c       Save histogram state whenever dyn files are created
c
c        if (mod(energyCount,printFrequency) .eq. 0) then
        if (mod(energyCount,iwrite) .eq. 0) then
          call prthis
          print *, "Wrote to histogram"
          print *, ""
        end if
      end if
      
c
c     compute the energy and gradient for the recursion slave
c     at F(L) using interpolation
c
      freeEnergy = currentFreeEnergy()
      biasEnergy = biasEnergy + freeEnergy
      
      if (verbose) then
c      if (.true.) then
      print *, "Bias Energy:", biasEnergy, freeEnergy
      print *, "OSRW Potential:", dte+biasEnergy, "(kcal/mol)"
      end if
      
      if (propagateLambda .and. energyCount .gt. 0) then
c            
c       Metadynamics grid counts (every 'countInterval' steps).
c           (10 steps)       
c
        if (mod(energyCount,countInterval) .eq. 0) then
c          call checkRecursionKernelSize(dEdU)
          recursionKernel(lambdaBin+1,FLambdaBin+1) = 
     &           recursionKernel(lambdaBin+1,FLambdaBin+1) + 1
        end if
    
      end if
      
      if (propagateLambda) then
        call langevin
      else
c        equilibration???  
      end if

c   
c     Final output of osrw
c     esum, desum, dudl, d2udl2, d2udl  
c 
      totalEnergy = dte + biasEnergy
      esum = totalEnergy
      return 
      end
      

c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine - checkRecursionKernelSize                       ##
c     ##                                                              ##
c     ##################################################################
c
c      Update RecursionKernel size if dudl falls outside of range
c             maxFLambda, minFLambda
c

      subroutine checkRecursionKernelSize (dudlambda)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'mutant.i'
      include 'osrwi.i'
      integer newFLambdaBins, i, offset, j
      real*8  dudlambda
      real*8  origDeltaG, updateFLambda
      real*8, allocatable :: rkCopy (:,:)
      
      
      if (dudlambda .gt. maxFLambda) then
        allocate (rkCopy(lambdaBins,FLambdabins))
        print *, "dudlambda > maxFLambda",dudlambda,maxFLambda
        origDeltaG = updateFLambda(.false.)
        
        newFLambdaBins = FLambdaBins
        do while (minFLambda + newFLambdaBins*dFL .lt. dudLambda)
          newFLambdaBins = newFLambdaBins + 100
        end do
        do i = 1, lambdaBins
          do j = 1, FLambdaBins
            rkCopy(i,j) = recursionKernel(i,j)
          end do
        end do
        deallocate(recursionKernel)
        allocate(recursionKernel(lambdaBins,newFLambdaBins))
c       initialize to 0
        do i = 1, lambdaBins
          do j = 1, newFLambdaBins
            recursionKernel(i,j) = 0
          end do
        end do
        do i = 1, lambdaBins
          do j = 1, FLambdaBins
            recursionKernel(i,j) = rkCopy(i,j)
          end do
        end do
        FLambdaBins = newFLambdaBins
        maxFLambda = minFLambda + dFL*FLambdaBins
        if (.not. origDeltaG .eq. updateFLambda(.false.)) then
          print *, "origDeltaG != updateFLambda"
          print *, origDeltaG, updateFLambda(.false.)
          call fatal
        end if
        deallocate(rkCopy)
      end if
      
      if (dudlambda .lt. minFLambda) then
        allocate (rkCopy(lambdaBins,FLambdabins))
        print *, "dudlambda < minFLambda",dudlambda,minFLambda
        origDeltaG = updateFLambda(.false.)
        
        offset = 100
        do while (dudLambda .lt. minFLambda - offset*dFL)
          offset = offset + 100
        end do
        newFLambdaBins = FLambdaBins + offset
        do i = 1, lambdaBins
          do j = 1, FLambdaBins
            rkCopy(i,j) = recursionKernel(i,j)
          end do
        end do
        deallocate(recursionKernel)
        allocate(recursionKernel(lambdaBins,newFLambdaBins))
c       initialize to 0
        do i = 1, lambdaBins
          do j = 1, newFLambdaBins
            recursionKernel(i,j) = 0
          end do
        end do
        do i = 1, lambdaBins
          do j = 1, FLambdaBins
            recursionKernel(i,j+offset) = rkCopy(i,j)
          end do
        end do
        FLambdaBins = newFLambdaBins
        minFLambda = minFLambda - offset*dFL
        if (.not. origDeltaG .eq. updateFLambda(.false.)) then
          print *, "origDeltaG != updateFLambda"
          print *, origDeltaG, updateFLambda(.false.)
          call fatal
        end if
        deallocate(rkCopy)
      end if
      
      return
      end
      

c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine - langevin                                       ##
c     ##                                                              ##
c     ##################################################################
c
c      Use langevin dynamics to update lambda particle position
c

      subroutine langevin
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bath.i'
      include 'math.i'
      include 'mutant.i'
      include 'osrwi.i'
      include 'units.i'
      real*8  R
      real*8  randomConvert
      real*8  rt2, randomForce, dEdLi, sinTheta
      real*8  randomConvert2, normal
      
c
c     Random force conversion to kcal/mol/A
c
      randomConvert = sqrt(4.184d0) / 10E9
      randomConvert2 = randomConvert*randomConvert
      
      rt2 = 2.0d0*gasconst*tempsystem*thetaFriction / dtosrw
      randomForce = sqrt(rt2)*(normal()) / randomConvert
c
c     du/d(theta) = dudl * dl/d(theta)
c
      dEdLi = -dudl * sin(2.0d0 * theta)
      
      
      halfThetaVelocity = (halfThetaVelocity 
     &                 * (2.0d0 * thetaMass - thetaFriction * dtosrw)
     &          + randomConvert2 * 2.0d0 * dtosrw 
     &          * (dEdLi + randomForce))
     &           / (2.0d0 * thetaMass + thetaFriction * dtosrw)

      theta = theta + dtosrw * halfThetaVelocity
      
      if (theta .gt. pi) then
        theta = theta - 2.0d0*pi
      else if (theta .le. -pi) then
        theta = theta + 2.0d0*pi
      end if

      sinTheta = sin(theta)
      lambda = sinTheta*sinTheta
      
      return
      end
      
c
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  function - currentFreeEnergy                                ##
c     ##                                                              ##
c     ##################################################################
c
c
c     compute the energy and gradient for the recursion slave
c     at F(L) using interpolation
c
c
c

      function currentFreeEnergy
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'mutant.i'
      include 'osrwi.i'
      integer iL0,iL1
      real*8  currentFreeEnergy, biasEnergy
      real*8  L0,L1,FL0,FL1,deltaFL
      logical done
      
      biasEnergy = 0.0d0
      iL0 = 0
      do while(iL0 .lt. lambdaBins - 1)
        iL1 = iL0 + 1
c
c       Find bin centers and values for interpolation/extrapolation 
c        points
c
        L0 = iL0*dL
        L1 = L0+ dL
        FL0 = FLambda(iL0+1)
        FL1 = FLambda(iL1+1)
        deltaFL = FL1-FL0
c            
c       If the lambda is less than or equal to the upper limit, this is
c       the final interval. Set the upper limit to L, compute the partial
c       derivative and break.
c                   
        done = .false.
        if (lambda .le. L1) then
          done = .true.
          L1 = lambda
        end if
c
c       Upper limit - lower limit of the integral of the 
c        extrapolation/interpolation
c
        biasEnergy = biasEnergy + (FL0 * L1 
     &                          + deltaFL * L1*(0.5d0 * L1 - L0) / dL)
        biasEnergy = biasEnergy - (FL0 * L0 
     &                          + deltaFL * L0 * (-0.5d0 * L0) / dL)
        if (done) then
c
c         Compute the gradient dF(L)/dL at L
c
          dudl = dudl - (FL0 + (L1 - L0) * deltaFL / dL)
          exit
        end if
        iL0 = iL0 + 1
      end do
      
      currentFreeEnergy = -biasEnergy
      
      return
      end
      
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  function - updateFLambda                                    ##
c     ##                                                              ##
c     ##################################################################
c
c
c

      function updateFLambda (printlambda)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bath.i'
      include 'mutant.i'
      include 'osrwi.i'
      include 'units.i'
      integer iL, jFL, ulFL, llFL, count, lambdaCount
      real*8  updateFLambda
      real*8  freeEnergy, lla, ula
      real*8  ensembleAverageFLambda, partitionFunction
      real*8  currentFLambda, evaluateKernel
      real*8  weight, delta, deltaFreeEnergy, llL, ulL
      logical printlambda
      
      if (printlambda) then
        print *, "Count LambdaBins F_LambdaBins <F_L> dG G"
      end if
      
      freeEnergy = 0.0d0
      totalCounts = 0
      tempsystem = 298.15d0
      
      iL = 0
      do while(iL .lt. lambdaBins)
        ulFL = -1
        llFL = -1
        
c       find the smallest FL bin
        jFL = 0
        do while (jFL .lt. FLambdaBins)
          count = recursionKernel(iL+1,jFL+1)
          if (count .gt. 0) then
            llFL = jFL
            exit
          end if
          jFL = jFL + 1
        end do

c       find the largest FL bin
        jFL = FLambdaBins - 1
        do while (jFL .ge. 0)
          count = recursionKernel(iL+1,jFL+1)
          if (count .gt. 0) then
            ulFL = jFL
            exit
          end if
          jFL = jFL - 1
        end do
        
        lambdaCount = 0
c   The FL range sample for lambda bin [iL*dL .. (iL+1)*dL]
        lla = 0.0d0
        ula = 0.0d0
        if (ulFL .eq. -1) then
          FLambda(iL+1) = 0.0d0
        else
          ensembleAverageFLambda = 0.0d0
          partitionFunction = 0.0d0
          jFL = llFL
          do while (jFL .le. ulFL)
            currentFLambda = minFLambda + jFL*dFL + dFL_2
            weight = exp(evaluateKernel(iL,jFL) /(gasconst*tempsystem))
            ensembleAverageFLambda = ensembleAverageFLambda
     &                              + currentFLambda*weight
            partitionFunction = partitionFunction + weight
            lambdaCount = lambdaCount + recursionKernel(iL+1,jFL+1)
            jFL = jFL + 1
          end do
          FLambda(iL+1) = ensembleAverageFLambda / partitionFunction
          lla = minFLambda + llFL*dFL
          ula = minFLambda + (ulFL+1)*dFL
        end if
        
c       The first and last lambda bins are half size        
        delta = dL
        if (iL .eq. 0 .or. iL .eq. lambdaBins - 1) then
          delta = dL_2
        end if
        deltaFreeEnergy = FLambda(iL+1)*delta
        freeEnergy = freeEnergy + deltaFreeEnergy
        totalCounts = totalCounts + lambdaCount

      
      if (printlambda) then
        llL = iL*dL - dL_2
        ulL = llL + dL
        if (llL .lt. 0.0d0) llL = 0.0d0
        if (ulL .gt. 1.0d0) ulL = 1.0d0
        write (*,40)  lambdaCount,llL,ulL,lla,ula,
     &                      FLambda(iL+1),deltaFreeEnergy,freeEnergy
   40   format (i5,2f7.3,2f10.1,3f10.2)

      end if
      
      iL = iL + 1
      end do
      
      updateFLambda = freeEnergy
      
      
      print *, "The free energy is ", freeEnergy, "kcal/mol"
      print *, "from count number",totalCounts
      print *, "current lam:", lambda,dudl
      print *, ""
      
      return
      end
      
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  function - evaluateKernel                                   ##
c     ##                                                              ##
c     ##################################################################
c
c
c

      function evaluateKernel (cLambdai, cF_Lambda)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'mutant.i'
      include 'osrwi.i'
      integer cLambdai, cF_Lambda
      integer iL, Lcenter, lcount, jFL, FLcenter
      real*8  evaluateKernel
      real*8  sum, vL, vFL, Ls2, FLs2
      real*8  deltaL, deltaL2, mirrorFactor, deltaFL, deltaFL2
      real*8  weight, ener
      
c
c     Compute the value of L and FL for the center of the current bin
c      
      vL = cLambdai*dL
      vFL = minFLambda + cF_Lambda*dFL + dFL_2
c
c     Set the variances for the Gaussian bias
c
      Ls2 = 2.0d0 * dL * 2.0d0 * dL
      FLs2 = 2.0d0 * dFL * 2.0d0 * dFL
      sum = 0.0d0
      
      iL = -biasCutoff
      do while (iL .le. biasCutoff)
        Lcenter = cLambdai + iL
        deltaL = vL - Lcenter*dL
        deltaL2 = deltaL*deltaL
c       mirror condition for lambda counts
        lcount = Lcenter
        mirrorFactor = 1.0d0
        if (lcount .eq. 0 .or. lcount .eq. lambdaBins - 1) then
c
c         The width of the first and last bins is dLambda_2, so the
c          mirror condition is to double their counts
c
          mirrorFactor = 2.0d0
        else if (lcount .lt. 0) then
          lcount = -lcount
        else if (lcount .gt. lambdaBins - 1) then
c         number of bins past the last bin
          lcount = lcount - (lambdaBins - 1)
c         mirror bin
          lcount = -lcount + (lambdaBins - 1)
        end if
        
        jFL = -biasCutoff
        do while (jFL .le. biasCutoff)
          FLcenter = cF_Lambda + jFL
c
c         For FLambda outside the count matrix the weight is 0
c          so we continue
c
          if (FLcenter .lt. 0 .or. FLcenter .ge. FLambdaBins) then
            jFL = jFL + 1
            cycle
          end if
          deltaFL = vFL - (minFLambda + FLcenter * dFL + dFL_2)
          deltaFL2 = deltaFL * deltaFL
          weight = mirrorFactor *recursionKernel(lcount+1,FLcenter+1)
          if (weight .gt. 0) then
            ener = weight * biasMag * exp(-deltaL2 / (2.0d0 * Ls2))
     &                       * exp(-deltaFL2 / (2.0d0 * FLs2))
            sum = sum + ener
          end if
          jFL = jFL + 1
        end do
        iL = iL + 1
      end do
      
      evaluateKernel = sum
      
      return
      end
      
c
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  function - binForLambda                                     ##
c     ##                                                              ##
c     ##################################################################
c
c
c

      function binForLambda (lambdai)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'mutant.i'
      include 'osrwi.i'
      integer binForLambda
      real*8  lambdai
      
      binForLambda = int( (lambdai - minLambda) / dL )
      if (binForLambda .lt. 0) binForLambda = 0
      if (binForLambda .ge. lambdaBins) then
        binForLambda = lambdaBins - 1
      end if
      
      return
      end
      
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  function - binForFLambda                                    ##
c     ##                                                              ##
c     ##################################################################
c
c
c

      function binForFLambda (dudlambda)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'mutant.i'
      include 'osrwi.i'
      integer binForFLambda
      real*8  dudlambda
      
      binForFLambda = int( (dudlambda - minFLambda) / dFL )
            
      if (binForFLambda .eq. FLambdaBins) then
        binForFLambda = FLambdaBins - 1
      end if
      
      if (binForFLambda .ge. FLambdaBins) then
        print *, "binForFLambda >= FLambdaBins",binForFLambda
        binForFLambda = FLambdaBins - 1
c        call fatal
      else if (binForFLambda .lt. 0) then
        print *, "binForFLambda < 0", binForFLambda
        binForFLambda = 0
c        call fatal
      end if
      
      return
      end

c
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine prthis  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "prthis"
c
c
      subroutine prthis
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'files.i'
      include 'inform.i'
      include 'mutant.i'
      include 'osrwi.i'
      include 'output.i'
      include 'titles.i'
      integer ihis,freeunit
      integer i,k
      integer size,crdsiz
      real*8 crdmin,crdmax
      logical opened, exist
      character*2 atmc
      character*2 crdc
      character*2 digc
      character*3 lc
      character*25 fstr
      character*120 xyzfile
      character*120 hisfile
c
c
c     write the current coordinates to a file after an error
c
      if (hisverbose) then
      ihis = freeunit ()
      hisfile = filename(1:leng)//'.his'
      call version (hisfile,'new')
      open (unit=ihis,file=hisfile,status='new')
      
      else if (.not. hisverbose) then
      ihis = freeunit ()
      hisfile = filename(1:leng)//'.his'
      inquire (file=hisfile,exist=exist)
      if (exist) then
         open (unit=ihis,file=hisfile,status='old')
         rewind (unit=ihis)
      else
         open (unit=ihis,file=hisfile,status='new')
      end if
      
      end if
      
      write(ihis,*) lambda, halfThetaVelocity, minFLambda, FLambdaBins
c
c     write out the coordinate line for each atom
c
      lc = '401'
      fstr = '('//lc//'i4)'
      do i = 1, lambdaBins
c        do k = 1, FLambdaBins
c          write (ihis,90) recursionKernel(i,k)
         write (ihis,fstr)  (recursionKernel(i,k),k=1,FLambdaBins)
c        end do
      end do    
    
c   90 format (//FLambdaBins//I6)   
      
      
      close (unit=ihis)
      return
      end
