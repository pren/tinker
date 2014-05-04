c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2014  by  Jay William Ponder  ##
c     ## Michael Schneiders, Jayvee Abella, Pengyu Ren ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  osrwi.i --  parameters for osrw implementation         ##
c     ##                                                         ##
c     #############################################################
c
c
c     dtosrw           time step of MD simulation in picoseconds
c     tempsystem       temperature of MD simulation in Kelvin
c     propagateLambda  flag to indicate movement of lambda particle
c     biasMag          height of biasing gaussian
c     lambdaBins       number of lambda bins in recursionKernel
c     FLambdaBins      number of flambda bins in recursionKernel
c     dL               width of lambda bin
c     dFL              width of flambda bin
c     dL_2             half width of lambda bin
c     dFL_2            half width of flambda bin
c     countInterval    steps between update of recursionKernel
c     energyCount      number of times osrw energy was called
c     biasCutoff       number bins from center with nonzero gaussian
c     totalEnergy      biased potential osrw energy
c     totalFreeEnergy  current free energy est from recursionKernel
c     minLambda        minimum value of lambda
c     minFLambda       minimum value of flambda
c     maxFLambda       maximum value of flambda
c     FLambda          ensemble average of flambda
c     fLambdaUpdates   number of times FLambda is updated
c     fLambdaPrintInterval steps*10 between recursionKernel summary
c     printFrequency   NOT USED
c     halfThetaVelocity  related to lambda particle velocity
c     theta              related to lambda particle value
c     totalCounts      total recursionKernel counts
c     thetaFriction      related to lambda particle friction
c     thetaMass          related to lambda particle mass
c     hisverbose       flag to print multiple histogram files
c     recursionKernel  histogram that saves (lambda,dudl) states
c     restartosrw    flag to indicate that osrw was restarted


      integer, pointer :: recursionKernel(:,:)
      integer FLambdaBins, lambdaBins
      integer energyCount, countInterval, biasCutoff
      integer fLambdaUpdates, fLambdaPrintInterval, printFrequency
      integer totalCounts
      real*8  dL, dFL, dFL_2, biasMag, dL_2
      real*8  minLambda, minFLambda, totalEnergy, totalFreeEnergy
      real*8  maxFLambda
      real*8  FLambda
      real*8  dtosrw, tempsystem, halfThetaVelocity, theta
      real*8  thetaMass, thetaFriction
      logical propagateLambda, hisverbose, restartosrw
     
      common /osrwi/ dL, dFL, dFL_2, biasMag, dL_2,
     &                minLambda, minFLambda, totalEnergy, 
     &                totalFreeEnergy,
     &                maxFLambda, FLambda(201),
     &                dtosrw, tempsystem,
     &                halfThetaVelocity, theta,
     &                thetaMass, thetaFriction,
     &                recursionKernel,
     &                FLambdaBins, lambdaBins, energyCount, 
     &                countInterval, biasCutoff,
     &                fLambdaUpdates, fLambdaPrintInterval,
     &                printFrequency, totalCounts,
     &                propagateLambda, hisverbose, restartosrw