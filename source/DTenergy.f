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
c     ##  function DTenergy  --  energy of an alchemical system  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "energy" calls the subroutines to calculate the potential
c     energy terms and sums up to form the total energy along with
c     derivatives wrt Cartesian and lambda variables
c
c
      function DTenergy (lam)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cutoff.i'
      include 'deriv.i'
      include 'energi.i'
      include 'ewald.i'
      include 'inter.i'
      include 'iounit.i'
      include 'mpole.i'
      include 'mutant.i'
      include 'polar.i'
      include 'potent.i'
      include 'rigid.i'
      include 'vdwpot.i'
      include 'virial.i'
      integer i,j,xyz_type,k
      real*8 lam
      real*8 DTenergy
      real*8 cutoff
      real*8 aewaldold
      real*8 DTesum
      real*8 etot
      real*8 plam,elam
      real*8 plamn, dplamn, d2plamn
      real*8 elamn, delamn, d2elamn
      real*8 plamn2, dplamn2, d2plamn2,plam2
      real*8 mpolecutold, vdwcutold
      real*8 totalbond, totalvdw, totalem
      real*8 totalpol, totalereal
      real*8 solverecip, solvpol, solvereal
      real*8 ligem, ligpol, restr_ener
      real*8 total2em
      real*8 total2pol, total2ereal
      real*8 lig2em, lig2pol
      real*8 dummy1, dummy2, dummy3
      real*8 dummy4, dummy5
      real*8 polarLambdaScale,virold(3,3)
      real*8, allocatable :: detot(:,:)
      logical use_polarold
      
      allocate (detot(3,n))

c
c     Zero out energy and derivative terms
c
      DTesum = 0.0d0
      dudl = 0.0d0
      d2udl2 = 0.0d0
      do i = 1,n
        do j = 1,3
          detot(j,i) = 0.0d0
          desum(j,i) = 0.0d0
          d2udl(j,i) = 0.0d0
        end do
      end do
      
c
c     Assert 0 <= lam <= 1, cant testlambdagrad for lam = 0,1
c
      if (lam .gt. 1.0d0) lam = 1.0d0
      if (lam .lt. 0.0d0) lam = 0.0d0
      
      do i = 1, 3
        do j = 1, 3
          virold(j,i) = 0.0d0
        end do
      end do
c     Use one lambda value
      elambda = lam
      vlambda = lam
         
c ********************************************************************
c     System energy with only one mutating entity
c ********************************************************************
c     totalbond + totalvdw + lam*totalrecip + (soft real - fxn of lam)
c      + (lam*totalpol - if pol >= polstart)
c      + (1-lam)*restraint (if pro-lig)
c
c     (1-lam)*solvrecip + solvreal (nonsoft interactions)
c       +  solvep (if pol <= polstart)
c       + (1-lam)*solvep (if pol >= polstart)
c
c     (1-lam)*ligem
c     + ligpol (if lam <= polstart)
c     + (1-lam)*ligpol (if lam >= polstart)

      if (.not. isrelative) then
      
c     Based on current lambda, calculate lambda^n
      polarLambdaScale = 1.0d0 / (polend - polstart)
      plam = (lam - polstart) * polarLambdaScale
      elam = lam
      
      if (scexpp .ge. 2.0d0) then
        plamn = plam**scexpp
        dplamn = scexpp*plam**(scexpp-1.0d0)
        d2plamn = scexpp*(scexpp-1.0d0)*plam**(scexpp-2.0d0)
      else if (scexpp .ge. 0.0d0) then
        plamn = plam
        dplamn = 1.0d0
        d2plamn = 0.0d0
      end if
      
      if (scexpm .ge. 2.0d0) then
        elamn = elam**scexpm
        delamn = scexpm*elam**(scexpm-1.0d0)
        d2elamn = scexpm*(scexpm-1.0d0)*elam**(scexpm-2.0d0)
      else if (scexpp .ge. 0.0d0) then
        elamn = elam
        delamn = 1.0d0
        d2elamn = 0.0d0
      end if
      
      
c     chain rule
      dplamn = dplamn*polarLambdaScale
      d2plamn = d2plamn*polarLambdaScale*polarLambdaScale
      

c --------------------------------------------------------------------
c     get full system energy
c      
c     totalbond + totalvdw + lam*totalrecip + (soft real - fxn of lam)
c      + (lam*totalpol - if pol >= polstart)
c      + (1-lam)*restraint (if pro-lig)
      
c     turn off polar calculations if lam <= polstart
c     scale recip and pol by elamn, plamn
      use_polarold = use_polar
      if (lam .le. polstart) use_polar = .false.
      dorealsplit = .true.
      call gradient2 (etot,detot)
      use_polar = use_polarold
      
c     em - ereal = erecip + eself
      totalbond = esum - ev - ep - em
      totalvdw = ev
      totalem = em
      totalereal = ereal
      totalpol = ep
      
c    Bonded energy and softcore vdw energies are calculated
c      for both systems in one call of ehal1.f (inside esum)
      DTesum = esum - ep - (em - ereal)
      DTesum = DTesum + elamn*(em - ereal)
      dudl = dedlv + dedlm + delamn*(em - ereal)
      d2udl2 = d2edl2v + d2edl2m + d2elamn*(em-ereal)
      do i = 1,n
        do j = 1,3
          detot(j,i) = desum(j,i) - dep(j,i) - (dem(j,i)-dereal(j,i))
          detot(j,i) = detot(j,i) + elamn*(dem(j,i)-dereal(j,i))
          d2udl(j,i) = d2edlv(j,i)+ d2edlg(j,i) + 
     &                            delamn*(dem(j,i)-dereal(j,i))
        end do
      end do

c      if (lam .le. polstart) then
c       do nothing
      if (lam .gt. polstart .and. lam .le. polend) then
        DTesum = DTesum + plamn*ep
        dudl = dudl + dplamn*ep
        d2udl2 = d2udl2 + d2plamn*ep
        do i = 1,n
          do j = 1,3
            detot(j,i) = detot(j,i) + plamn*dep(j,i)
            d2udl(j,i) = d2udl(j,i) + dplamn*dep(j,i)
          end do
        end do
      end if
      
c      group restraint
      if (use_geom) then
        DTesum = DTesum - eg
        DTesum = DTesum + (1-lam)*eg
        dudl = dudl - lam*eg
        d2udl2 = d2udl2 - eg
        do i = 1,n
          do j = 1,3
            detot(j,i) = detot(j,i) - deg(j,i)
            detot(j,i) = detot(j,i) + (1-lam)*deg(j,i)
            d2udl(j,i) = d2udl(j,i) - lam*deg(j,i)
          end do
        end do
      end if
c --------------------------------------------------------------------

                  
c --------------------------------------------------------------------
c     get solvent only energy
c
c     (1-lam)*solvrecip + solvreal (nonsoft interactions)
c       +  solvep (if pol <= polstart)
c       + (1-lam)*solvep (if pol >= polstart)
        
      xyz_type = 1
      call changeXYZ(xyz_type)
      use_polarold = use_polar

      osrwon = .false.
      dorealsplit = .false.
      call eleconly 
      osrwon = .true.
      use_polar = use_polarold
      call restoreXYZ(xyz_type)

      solverecip = em - ereal
      solvpol = ep
      solvereal = ereal
      DTesum = DTesum + (1-elamn)*(em - ereal) + ereal
      dudl = dudl - delamn*(em - ereal)
      d2udl2 = d2udl2 - d2elamn*(em - ereal)
      do i = 1,n
        do j = 1,3
          detot(j,i) = detot(j,i) + (1-elamn)*(dem(j,i)-dereal(j,i))
     &                            + dereal(j,i)
          d2udl(j,i) = d2udl(j,i) - delamn*(dem(j,i)-dereal(j,i))
        end do
      end do

      if (lam .le. polstart) then
        DTesum = DTesum + ep
        do i = 1,n
          do j = 1,3
            detot(j,i) = detot(j,i) + dep(j,i)
          end do
        end do
      else if (lam .le. polend) then
        DTesum = DTesum + (1-plamn)*ep
        dudl = dudl - dplamn*ep
        d2udl2 = d2udl2  - d2plamn*ep
        do i = 1,n
          do j = 1,3
            detot(j,i) = detot(j,i) + (1-plamn)*dep(j,i)
            d2udl(j,i) = d2udl(j,i)-dplamn*dep(j,i)
          end do
        end do
      end if
c --------------------------------------------------------------------


            
c --------------------------------------------------------------------
c     get energy for ligand only in vacuum
c
c     (1-lam)*ligem
c     + ligpol (if lam <= polstart)
c     + (1-lam)*ligpol (if lam >= polstart)

      if (doVaporElec) then
      xyz_type = 2
      call changeXYZ(xyz_type)
c     turn off ewald to just run double loop
      aewaldold = aewald
c      aewald = 0.0d0
      osrwon = .false.
      use_mlist = .false.
      use_list = .false.
      use_ewald = .false.
      use_bounds = .false.
      mpolecutold = mpolecut
      vdwcutold = vdwcut
      mpolecut = 1.0d12
      vdwcut = 1.0d12
      call eleconly 
      vdwcut = vdwcutold
      mpolecut = mpolecutold
      use_bounds = .true.
      use_mlist = .true.
      use_ewald = .true.
      use_list = .true.
      osrwon = .true.
      aewald = aewaldold

      ligem = em
      ligpol = ep
      DTesum = DTesum + (1-elamn)*em
      dudl = dudl - delamn*em
      d2udl2 = d2udl2 - d2elamn*em
      do i = 1,n
        do j = 1,3
          detot(j,i) = detot(j,i) + (1-elamn)*dem(j,i)
          d2udl(j,i) = d2udl(j,i) - delamn*dem(j,i)
        end do
      end do
      call restoreXYZ(xyz_type)
      
      if (lam .le. polstart) then
        DTesum = DTesum + ep
        do i = 1,n
          do j = 1,3
            detot(j,i) = detot(j,i) + dep(j,i)
          end do
        end do
      else if (lam .le. polend) then
        DTesum = DTesum  + (1-plamn)*ep
        dudl = dudl - dplamn*ep
        d2udl2 = d2udl2 - d2plamn*ep
        do i = 1,n
          do j = 1,3
            detot(j,i) = detot(j,i) + (1-plamn)*dep(j,i)
            d2udl(j,i) = d2udl(j,i)-dplamn*dep(j,i)
          end do
        end do
      end if   
      end if
c --------------------------------------------------------------------
      
c     debug
c      print *, "DTdudl", dudl
c      print *, "reg dudl",dum1,dum2
c      if (dudl .gt. 500d0) then
c        print *, "warning gt",dudl,lam
c        call prterr
c        call fatal
c      else if (dudl .lt. -500d0) then
c        print *, "warning lt",dudl,lam
c        call prterr
c        call fatal
c      end if

      if (debugjv) then
        print *, "totalbond", totalbond
        print *, "totalvdw", totalvdw
        if (lam .eq. 0.0d0) then
          print *, "em0", solvereal+solverecip+ligem
          print *, "pol0", solvpol+ligpol
        else if (lam .eq. 1.0d0) then
          print *, "em1", totalem-totalereal+solvereal
          print *, "pol1", totalpol
        end if
      end if
c ********************************************************************
c     End Absolute
c ********************************************************************       
      
c ********************************************************************
c     System with two mutating entities
c ********************************************************************   
c     totalbond (one call of gradient)
c     totalvdw (one call of ehal1.f - mutating systems
c                         do not interact with each other)
c     ( = Uenv-A_VDW(1-lam) + Uenv-B_VDW(lam) + Uenv_VDW + UA + UB )
c     UA_real(1-lam) + Uenv-A_real(1-lam) + Uenv_real
c     UB_real(lam) + Uenv-B(lam)
c     UA_em_vac(lam) + UB_em_vac(1-lam)
c     
c     0 < lam < 0.25
c        lam1*(Uenv-A_pol) + (1-lam1)*(Uenv_pol + UA_pol) + UB_pol
c     0.25 < lam < 0.75
c        Uenv_pol + UA_pol + UB_pol
c     0.75 < lam < 1
c        lam2*(Uenv-B_pol) + (1-lam2)*(Uenv_pol + UB_pol) + UA_pol
c     (1-lam)*total1recip + lam*total2recip
c     
c     Also includes implicit restraint bond
      else if (isrelative) then
      
c     Based on current lambda, calculate lambda^n
      polarLambdaScale = 1.0d0 / (0.25d0)
      plam = (0.25d0 - lam) * polarLambdaScale
      plam2 = (lam - 0.75d0) * polarLambdaScale
      elam = lam
      
      if (scexpp .ge. 2.0d0) then
        plamn = plam**scexpp
        dplamn = -scexpp*plam**(scexpp-1.0d0)
        d2plamn = scexpp*(scexpp-1.0d0)*plam**(scexpp-2.0d0)
        plamn2 = plam2**scexpp
        dplamn2 = scexpp*plam2**(scexpp-1.0d0)
        d2plamn2 = scexpp*(scexpp-1.0d0)*plam2**(scexpp-2.0d0)
      else if (scexpp .ge. 0.0d0) then
        plamn = plam
        dplamn = -1.0d0
        d2plamn = 0.0d0
        plamn = plam2
        dplamn = 1.0d0
        d2plamn = 0.0d0
      end if
      
      if (scexpm .ge. 2.0d0) then
        elamn = elam**scexpm
        delamn = scexpm*elam**(scexpm-1.0d0)
        d2elamn = scexpm*(scexpm-1.0d0)*elam**(scexpm-2.0d0)
      else if (scexpp .ge. 0.0d0) then
        elamn = elam
        delamn = 1.0d0
        d2elamn = 0.0d0
      end if
      
      
c     chain rule
      dplamn = dplamn*polarLambdaScale
      d2plamn = d2plamn*polarLambdaScale*polarLambdaScale
      dplamn2 = dplamn2*polarLambdaScale
      d2plamn2 = d2plamn2*polarLambdaScale*polarLambdaScale    
      
c --------------------------------------------------------------------
c     get energies for total system 1
c
c     totalbond + totalvdw + (1-lam)*total1recip
c      + (UA_real(1-lam) + Uenv-A_real(1-lam))
c      + lam1*total1pol (if lam <= .25)

      rel_ligA = .true.
      use_polarold = use_polar
      xyz_type = 5
      if (lam .gt. 0.25d0) use_polar = .false.
      call changeXYZ(xyz_type) 

c     calculate energy
      dorealsplit = .true.
      call gradient2 (etot,detot)
      use_polar = use_polarold

c     restore xyz
      call restoreXYZ(xyz_type)
      
c     update DTesum
       totalbond = esum - ev - em - ep
       totalvdw = ev
       totalereal = ereal
       totalem = em
       totalpol = ep
       restr_ener = eg
c    Bonded energy and softcore vdw energies are calculated
c      for both systems in one call of ehal1.f
       DTesum = esum - ep - (em - ereal) + (1-elamn)*(em-ereal)
       dudl = dedlv + dedlm - delamn*(em-ereal)
       d2udl2 = d2edl2v + d2edl2m - d2elamn*(em-ereal)
       
      do i = 1,n
        do j = 1,3
          detot(j,i) = desum(j,i) - dep(j,i) - (dem(j,i)-dereal(j,i))
     &                              + (1-elamn)*(dem(j,i)-dereal(j,i))
          d2udl(j,i) = d2edlv(j,i)+ d2edlg(j,i) 
     &            - delamn*(dem(j,i)-dereal(j,i))
        end do
      end do
      if (lam .le. 0.25d0) then
        DTesum = DTesum + (plamn)*ep
        dudl = dudl + dplamn*ep
        d2udl2 = d2udl2 + d2plamn*ep
        do i = 1,n
          do j = 1,3
            detot(j,i) = detot(j,i) + (plamn)*dep(j,i) 
            d2udl(j,i) = d2udl(j,i) + dplamn*dep(j,i)
          end do
        end do
      end if
      
      dummy1 = DTesum
c      print *, esum, ep, em, ereal
      
c      if (use_geom) then
c        print *, "Geometric restraint not yet implemented with"
c        print *, "osrw-relative simulations"
c        call fatal
c      end if
c --------------------------------------------------------------------
      

c --------------------------------------------------------------------
c     get energies for total system 2
c
c     lam*total2recip + 
c     (Uenv-B_real(lam) + UB_real(lam))
c     + lam2*total2pol (if lam2 >= .75)
 
      rel_ligA = .false.
      xyz_type = 6
      call changeXYZ(xyz_type) 

c     calculate energy
      use_polarold = use_polar
      if (lam .lt. 0.75d0) use_polar = .false.
      dorealsplit = .true.
      call eleconly
      use_polar = use_polarold
      
c     restore xyz
      call restoreXYZ(xyz_type)  
c     update DTesum
      total2ereal = ereal
      total2em = em
      total2pol = ep
      DTesum = DTesum + elamn*(em-ereal) + ereal
      dudl = dudl + delamn*(em-ereal) + dedlm
      d2udl2 = d2udl2 + d2elamn*(em-ereal) + d2edl2m
      do i = 1,n
        do j = 1,3
          detot(j,i) = detot(j,i) + (elamn)*(dem(j,i)-dereal(j,i))
     &                       + dereal(j,i)
          d2udl(j,i) = d2udl(j,i) + (delamn)*(dem(j,i)-dereal(j,i))
     &                       + d2edlg(j,i)
        end do
      end do
      if (lam .ge. 0.75d0) then
        DTesum = DTesum + (plamn2)*ep
        dudl = dudl + dplamn2*ep
        d2udl2 = d2udl2 + d2plamn2*ep
        do i = 1,n
          do j = 1,3
            detot(j,i) = detot(j,i) + (plamn2)*dep(j,i)
            d2udl(j,i) = d2udl(j,i) + dplamn2*dep(j,i)
          end do
        end do
      end if
      
      dummy2 = DTesum - dummy1
c --------------------------------------------------------------------
      
      
c --------------------------------------------------------------------
c    get energy for environment
c
c    solvreal (Uenv_real)
c    + (1-lam)*ep (if lam <= .25)
c    + (1-lam2)*ep (if lam >= .75)
c    + ep (else)

c     changexyz
      xyz_type = 1
      call changeXYZ(xyz_type)
c     Condensed phase SCF w/o the ligand present
c     For DualTopologyEnergy calculations it can be turned off.
c      if (.not. doNoLigandCondensed) then
c        use_polar = .false.
c      end if
c     calculate energy
      osrwon = .false.
      dorealsplit = .false.
      call eleconly
      osrwon = .true.
c     restore xyz  
      call restoreXYZ(xyz_type) 
c     update DTesum 
      
      solvpol = ep
      solvereal = ereal
      
      DTesum = DTesum + ereal
      do i = 1,n
        do j = 1,3
          detot(j,i) = detot(j,i) + dereal(j,i)
        end do
      end do
      if (lam .le. 0.25d0) then
        DTesum = DTesum + (1-plamn)*ep
        dudl = dudl - dplamn*ep
        d2udl2 = d2udl2 - d2plamn*ep
        do i = 1,n
          do j = 1,3
            detot(j,i) = detot(j,i) + (1-plamn)*dep(j,i)
            d2udl(j,i) = d2udl(j,i) - dplamn*dep(j,i)
          end do
        end do
      else if (lam .ge. 0.75d0) then
        DTesum = DTesum + (1-plamn2)*ep
        dudl = dudl - dplamn2*ep
        d2udl2 = d2udl2 - d2plamn2*ep
        do i = 1,n
          do j = 1,3
            detot(j,i) = detot(j,i) + (1-plamn2)*dep(j,i)
            d2udl(j,i) = d2udl(j,i) - dplamn2*dep(j,i)
          end do
        end do
      else
        DTesum = DTesum + ep
        do i = 1,n
          do j = 1,3
            detot(j,i) = detot(j,i) + dep(j,i)
          end do
        end do
      end if
      
      dummy3 = DTesum - dummy1 - dummy2
c --------------------------------------------------------------------

c --------------------------------------------------------------------
c     get energies for ligand 1 only  
c
c     lam*ligem
c     + (1-lam)*ligep (if lam <= .25)
c     + ligep (else)

      xyz_type = 3
      call changeXYZ(xyz_type) 
c     turn off ewald to just run double loop
      aewaldold = aewald
c      aewald = 0.0d0
      osrwon = .false.
      use_mlist = .false.
      use_ewald = .false.
      use_bounds = .false.
      use_list = .false.
      mpolecutold = mpolecut
      vdwcutold = vdwcut
      mpolecut = 1.0d12
      vdwcut = 1.0d12
      call eleconly 
      vdwcut = vdwcutold
      mpolecut = mpolecutold
      use_bounds = .true.
      use_mlist = .true.
      use_ewald = .true.
      use_list = .true.
      osrwon = .true.
      aewald = aewaldold

      ligem = em
      ligpol = ep
      DTesum = DTesum + (elamn)*em
      dudl = dudl + delamn*em
      d2udl2 = d2udl2 + d2elamn*em
      do i = 1,n
        do j = 1,3
          detot(j,i) = detot(j,i) + (elamn)*dem(j,i)
          d2udl(j,i) = d2udl(j,i) + delamn*dem(j,i)
        end do
      end do
      call restoreXYZ(xyz_type)
      
      if (lam .le. 0.25d0) then
        DTesum = DTesum  + (1-plamn)*ep
        dudl = dudl - dplamn*ep
        d2udl2 = d2udl2 - d2plamn*ep
        do i = 1,n
          do j = 1,3
            detot(j,i) = detot(j,i) + (1-plamn)*dep(j,i)
            d2udl(j,i) = d2udl(j,i)-dplamn*dep(j,i)
          end do
        end do
      else 
        DTesum = DTesum + ep
        do i = 1,n
          do j = 1,3
            detot(j,i) = detot(j,i) + dep(j,i)
          end do
        end do
      end if   
      
      dummy4 = DTesum-dummy1-dummy2-dummy3
c --------------------------------------------------------------------


c --------------------------------------------------------------------
c     get energies for ligand 2 only  
c
c     (1-lam)*lig2em
c     + (1-lam2)*lig2pol (if lam >= .75)
c     + lig2pol (else)

      xyz_type = 4
      call changeXYZ(xyz_type) 
c     turn off ewald to just run double loop
      aewaldold = aewald
c      aewald = 0.0d0
      osrwon = .false.
      use_mlist = .false.
      use_ewald = .false.
      use_bounds = .false.
      use_list = .false.
      mpolecutold = mpolecut
      vdwcutold = vdwcut
      mpolecut = 1.0d12
      vdwcut = 1.0d12
      call eleconly 
      vdwcut = vdwcutold
      mpolecut = mpolecutold
      use_bounds = .true.
      use_mlist = .true.
      use_ewald = .true.
      use_list = .true.
      osrwon = .true.
      aewald = aewaldold

      lig2em = em
      lig2pol = ep
      DTesum = DTesum + (1-elamn)*em
      dudl = dudl - delamn*em
      d2udl2 = d2udl2 - d2elamn*em
      do i = 1,n
        do j = 1,3
          detot(j,i) = detot(j,i) + (1-elamn)*dem(j,i)
          d2udl(j,i) = d2udl(j,i) - delamn*dem(j,i)
        end do
      end do
      call restoreXYZ(xyz_type)
      
      if (lam .ge. 0.75d0) then
        DTesum = DTesum  + (1-plamn2)*ep
        dudl = dudl - dplamn2*ep
        d2udl2 = d2udl2 - d2plamn2*ep
        do i = 1,n
          do j = 1,3
            detot(j,i) = detot(j,i) + (1-plamn2)*dep(j,i)
            d2udl(j,i) = d2udl(j,i)-dplamn2*dep(j,i)
          end do
        end do
      else 
        DTesum = DTesum + ep
        do i = 1,n
          do j = 1,3
            detot(j,i) = detot(j,i) + dep(j,i)
          end do
        end do
      end if  
      
      dummy5 = DTesum-dummy1-dummy2-dummy3-dummy4
c --------------------------------------------------------------------
      
      if (debugjv) then
        print *, "totalbond", totalbond
        print *, "totalvdw", totalvdw
        print *, "solpol", solvpol
        print *, "recipA ligApol", totalem-totalereal, ligpol
        print *, "recipB ligBpol", total2em-total2ereal, lig2pol
        print *, "Restraint term", restr_ener
        print *, dummy1,dummy2,dummy3
        print *, dummy4,dummy5
        if (lam .eq. 0.0d0) then
          print *, "emA", totalem-totalereal+ligem+solvereal
          print *, "ligBonly", lig2em
          print *, "polA", totalpol+ligpol
        else if (lam .eq. 1.0d0) then
          print *, "emB", total2em-total2ereal+lig2em+solvereal
          print *, "ligAonly",ligem
          print *, "polB", total2pol+lig2pol
        end if
      end if
      
      
      end if 
c ********************************************************************
c     End Relative
c ********************************************************************

c   
c     Final output of DTenergy
c     DTenergy, desum, dudl, d2udl2, d2udl  
c 
      do i = 1, 3
        do j = 1, 3
          vir(j,i) = virold(j,i)
        end do
      end do
      
      do i = 1,n
        do j = 1,3
          desum(j,i) = detot(j,i)
        end do
      end do
      
      esum = DTesum
      DTenergy = DTesum
      
c
c     check for an illegal value for the total energy
c
      if (isnan(DTesum)) then
         write (iout,50)
   50    format (/,' DTENERGY  --  Illegal Value for the Total',
     &              ' Potential Energy')
         call fatal
      end if 
          
      deallocate(detot)
      
      return
      end
      

c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine changeXYZ - change which atoms interact         ##
c     ##                                                              ##
c     ##################################################################
c
c
c

      subroutine changeXYZ (xyz_type)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'mpole.i'
      include 'mutant.i'
      include 'polar.i'
      integer i,j,k
      integer xyz_type
      real*8, allocatable :: todelete(:)

      allocate (todelete(n))
c      
c     copy original (total) information into old
c
      do i = 1,n
        do j = 1, 13
          poleold(j,i) = pole(j,i)
        end do
        polarityold(i) = polarity(i)
      end do
      
c     1 - solv only, 2 - lig only
c     3 - lig1 only, 4 - lig2 only
c     5 - tot1 only, 6 - tot2 only      

c      
c     decide which atoms to delete
c
      do i = 1,n
         todelete(i) = 0.0d0
      end do
      if (xyz_type .eq. 1) then
         if (isrelative) then
           do i = 1,nmut1
              todelete(imut1(i)) = 1.0d0
           end do
           do i = 1,nmut2
              todelete(imut2(i)) = 1.0d0
           end do
         else 
           do i = 1,nmut
              todelete(imut(i)) = 1.0d0
           end do
         end if
      else if (xyz_type .eq. 2) then
         do i = 1,n
           todelete(i) = 1.0d0
         end do
         do i = 1,nmut
           todelete(imut(i)) = 0.0d0
         end do
      else if (xyz_type .eq. 3) then
         do i = 1,n
           todelete(i) = 1.0d0
         end do
         do i = 1,nmut1
           todelete(imut1(i)) = 0.0d0
         end do
      else if (xyz_type .eq. 4) then
         do i = 1,n
           todelete(i) = 1.0d0
         end do
         do i = 1,nmut2
           todelete(imut2(i)) = 0.0d0
         end do
      else if (xyz_type .eq. 5) then
         do i = 1,nmut2
            todelete(imut2(i)) = 1.0d0
         end do
      else if (xyz_type .eq. 6) then
         do i = 1,nmut1
            todelete(imut1(i)) = 1.0d0
         end do
      end if

c      
c     Zero out electrostatic and polarization parameters
c      
         do i = 1, npole
            k = ipole(i)
            if (todelete(k) .eq. 1.0d0) then
               do j = 1, 13
                   pole(j,i) = 0.0d0
               end do
            end if
         end do
         do i = 1, npole
           k = ipole(i)
           if (todelete(k) .eq. 1.0d0) then
               polarity(i) = 0.0d0
            end if
         end do
  
c
c    Reduce loops to save time in vacuum calculation
c  
      if (xyz_type .eq. 2) then
        npole = nmut
        npolar = nmut
      else if (xyz_type .eq. 3) then
        npole = nmut
        npolar = nmut
      else if (xyz_type .eq. 4) then
        npole = nmut
        npolar = nmut
      end if
      
      deallocate (todelete)
      return
      end      
      
      

c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine restoreXYZ - restore to total system            ##
c     ##                                                              ##
c     ##################################################################
c
c
c

      subroutine restoreXYZ (xyz_type)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'mpole.i'
      include 'mutant.i'
      include 'polar.i'
      integer i,j,k,xyz_type
      
      if (xyz_type .eq. 2) then
        npole = n
        npolar = n
      else if (xyz_type .eq. 3) then
        npole = n
        npolar = n
      else if (xyz_type .eq. 4) then
        npole = n
        npolar = n
      end if
      
      do i = 1,n
        do j = 1, 13
          pole(j,i) = poleold(j,i)
        end do
        polarity(i) = polarityold(i)
      end do
      
      return
      end
      
      
      
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine eleconly  --  find energy & gradient components  ##
c     ##                      of electrostatic and polarization only  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "eleconly" calls subroutines to calculate the potential energy
c     and first derivatives with respect to Cartesian coordinates and lambda
c     for electrostatic and polarization only
c
c
      subroutine eleconly 
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cutoff.i'
      include 'deriv.i'
      include 'energi.i'
      include 'inter.i'
      include 'iounit.i'
      include 'potent.i'
      include 'rigid.i'
      include 'vdwpot.i'
      include 'virial.i'
      integer i,j
      real*8 cutoff
      

c
c
c     zero out each of the potential energy components
c      
      em = 0.0d0
      ep = 0.0d0
      
c
c     zero out each of the first derivative components
c
      do i = 1, n
         do j = 1, 3
            dem(j,i) = 0.0d0
            dep(j,i) = 0.0d0
         end do
      end do
c
c     zero out the virial and the intermolecular energy
c
      do i = 1, 3
         do j = 1, 3
            vir(j,i) = 0.0d0
         end do
      end do
      einter = 0.0d0

c
c     maintain any periodic boundary conditions
c
      if (use_bounds .and. .not.use_rigid)  call bounds
c
c     update the pairwise interaction neighbor lists
c
      if (use_list)  call nblist
c
c     remove any previous use of the replicates method
c
      cutoff = 0.0d0
      call replica (cutoff)
c
c     many implicit solvation models require Born radii
c
      if (use_born)  call born
c
c     alter bond and torsion constants for pisystem
c
      if (use_orbit)  call picalc

c
c     call the electrostatic energy and gradient routines
c
      if (use_charge)  call echarge1
      if (use_chgdpl)  call echgdpl1
      if (use_dipole)  call edipole1
      if (use_mpole .or. use_polar)  call empole1
      if (use_rxnfld)  call erxnfld1
      
      if (.false.) then
        print *, "elec"
        print *, "ev,em,ep"
        print *, ev,em,ep
      end if

      
c
c     check for an illegal value for the total energy
c
      if (isnan(esum)) then
         write (iout,10)
   10    format (/,' ELECONLY  --  Illegal Value for the Total',
     &              ' Potential Energy')
         call fatal
      end if
      return
      end      
      
      

c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine gradient2  --  find energy & gradient components ##
c     ##                                  (same as gradient)          ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "gradient2" is the same is gradient without the redirecting statement
c         to osrw
c
c
      subroutine gradient2 (energy,derivs)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cutoff.i'
      include 'deriv.i'
      include 'energi.i'
      include 'inter.i'
      include 'iounit.i'
      include 'potent.i'
      include 'rigid.i'
      include 'vdwpot.i'
      include 'virial.i'
      integer i,j
      real*8 energy,cutoff
      real*8 derivs(3,*)
      
c
c
c     zero out each of the potential energy components
c
      eb = 0.0d0
      ea = 0.0d0
      eba = 0.0d0
      eub = 0.0d0
      eaa = 0.0d0
      eopb = 0.0d0
      eopd = 0.0d0
      eid = 0.0d0
      eit = 0.0d0
      et = 0.0d0
      ept = 0.0d0
      ebt = 0.0d0
      ett = 0.0d0
      ev = 0.0d0
      ec = 0.0d0
      ecd = 0.0d0
      ed = 0.0d0
      em = 0.0d0
      ep = 0.0d0
      er = 0.0d0
      es = 0.0d0
      elf = 0.0d0
      eg = 0.0d0
      ex = 0.0d0
c
c     zero out each of the first derivative components
c
      do i = 1, n
         do j = 1, 3
            deb(j,i) = 0.0d0
            dea(j,i) = 0.0d0
            deba(j,i) = 0.0d0
            deub(j,i) = 0.0d0
            deaa(j,i) = 0.0d0
            deopb(j,i) = 0.0d0
            deopd(j,i) = 0.0d0
            deid(j,i) = 0.0d0
            deit(j,i) = 0.0d0
            det(j,i) = 0.0d0
            dept(j,i) = 0.0d0
            debt(j,i) = 0.0d0
            dett(j,i) = 0.0d0
            dev(j,i) = 0.0d0
            dec(j,i) = 0.0d0
            decd(j,i) = 0.0d0
            ded(j,i) = 0.0d0
            dem(j,i) = 0.0d0
            dep(j,i) = 0.0d0
            der(j,i) = 0.0d0
            des(j,i) = 0.0d0
            delf(j,i) = 0.0d0
            deg(j,i) = 0.0d0
            dex(j,i) = 0.0d0
         end do
      end do
c
c     zero out the virial and the intermolecular energy
c
      do i = 1, 3
         do j = 1, 3
            vir(j,i) = 0.0d0
         end do
      end do
      einter = 0.0d0
c
c     maintain any periodic boundary conditions
c
      if (use_bounds .and. .not.use_rigid)  call bounds
c
c     update the pairwise interaction neighbor lists
c
      if (use_list)  call nblist
c
c     remove any previous use of the replicates method
c
      cutoff = 0.0d0
      call replica (cutoff)
c
c     many implicit solvation models require Born radii
c
      if (use_born)  call born
c
c     alter bond and torsion constants for pisystem
c
      if (use_orbit)  call picalc
c
c     call the local geometry energy and gradient routines
c
      if (use_bond)  call ebond1
      if (use_angle)  call eangle1
      if (use_strbnd)  call estrbnd1
      if (use_urey)  call eurey1
      if (use_angang)  call eangang1
      if (use_opbend)  call eopbend1
      if (use_opdist)  call eopdist1
      if (use_improp)  call eimprop1
      if (use_imptor)  call eimptor1
      if (use_tors)  call etors1
      if (use_pitors)  call epitors1
      if (use_strtor)  call estrtor1
      if (use_tortor)  call etortor1
c
c     call the van der Waals energy and gradient routines
c
      if (use_vdw) then
         if (vdwtyp .eq. 'LENNARD-JONES')  call elj1
         if (vdwtyp .eq. 'BUCKINGHAM')  call ebuck1
         if (vdwtyp .eq. 'MM3-HBOND')  call emm3hb1
         if (vdwtyp .eq. 'BUFFERED-14-7')  call ehal1
         if (vdwtyp .eq. 'GAUSSIAN')  call egauss1
      end if
c
c     call the electrostatic energy and gradient routines
c
      if (use_charge)  call echarge1
      if (use_chgdpl)  call echgdpl1
      if (use_dipole)  call edipole1
      if (use_mpole .or. use_polar)  call empole1
      if (use_rxnfld)  call erxnfld1
c
c     call any miscellaneous energy and gradient routines
c
      if (use_solv)  call esolv1
      if (use_metal)  call emetal1
      if (use_geom)  call egeom1
      if (use_extra)  call extra1
      
c
c     sum up to get the total energy and first derivatives
c
      esum = eb + ea + eba + eub + eaa + eopb + eopd + eid + eit
     &          + et + ept + ebt + ett + ev + ec + ecd + ed + em
     &          + ep + er + es + elf + eg + ex
      energy = esum
            
      do i = 1, n
         do j = 1, 3
            desum(j,i) = deb(j,i) + dea(j,i) + deba(j,i)
     &                      + deub(j,i) + deaa(j,i) + deopb(j,i)
     &                      + deopd(j,i) + deid(j,i) + deit(j,i)
     &                      + det(j,i) + dept(j,i) + debt(j,i)
     &                      + dett(j,i) + dev(j,i) + dec(j,i)
     &                      + decd(j,i) + ded(j,i) + dem(j,i)
     &                      + dep(j,i) + der(j,i) + des(j,i)
     &                      + delf(j,i) + deg(j,i) + dex(j,i)
            derivs(j,i) = desum(j,i)
         end do
      end do
            
c
c     check for an illegal value for the total energy
c
      if (isnan(esum)) then
         write (iout,10)
   10    format (/,' GRADIENT2  --  Illegal Value for the Total',
     &              ' Potential Energy')
         call fatal
      end if
      return
      end
