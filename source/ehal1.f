c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ehal1  --  buffered 14-7 energy & derivatives  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ehal1" calculates the buffered 14-7 van der Waals energy and
c     its first derivatives with respect to Cartesian coordinates
c
c
      subroutine ehal1
      implicit none
      real*8 elrc,vlrc
      include 'cutoff.i'
      include 'energi.i'
      include 'vdwpot.i'
      include 'virial.i'
c
c
c     choose the method for summing over pairwise interactions
c
      if (use_lights) then
         call ehal1b
      else if (use_vlist) then
         call ehal1c
      else
         call ehal1a
      end if
c
c     apply long range van der Waals correction if desired
c
      if (use_vcorr) then
         call evcorr1 (elrc,vlrc)
         ev = ev + elrc
         vir(1,1) = vir(1,1) + vlrc
         vir(2,2) = vir(2,2) + vlrc
         vir(3,3) = vir(3,3) + vlrc
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine ehal1a  --  double loop buffer 14-7 vdw derivs  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "ehal1a" calculates the buffered 14-7 van der Waals energy and
c     its first derivatives with respect to Cartesian coordinates
c     using a pairwise double loop
c
c
      subroutine ehal1a
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'couple.i'
      include 'deriv.i'
      include 'energi.i'
      include 'group.i'
      include 'inter.i'
      include 'molcul.i'
      include 'mutant.i'
      include 'shunt.i'
      include 'usage.i'
      include 'vdw.i'
      include 'vdwpot.i'
      include 'virial.i'
      integer i,j,k
      integer ii,iv,it
      integer kk,kv,kt
      integer, allocatable :: iv14(:)
      real*8 e,de,eps,rdn
      real*8 fgrp,rv,rv7
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 redi,rediv
      real*8 redk,redkv
      real*8 dedx,dedy,dedz
      real*8 rho,rho6,rho7
      real*8 tau,tau7,scal
      real*8 s1,s2,t1,t2
      real*8 dt1drho,dt2drho
      real*8 dtau,gtau
      real*8 taper,dtaper
      real*8 rik,rik2,rik3
      real*8 rik4,rik5
      real*8 rik6,rik7
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8, allocatable :: xred(:)
      real*8, allocatable :: yred(:)
      real*8, allocatable :: zred(:)
      real*8, allocatable :: vscale(:)
      logical proceed,usei
      logical muti,mutk
      character*6 mode
c
c
c     zero out the van der Waals energy and first derivatives
c
      ev = 0.0d0
      do i = 1, n
         dev(1,i) = 0.0d0
         dev(2,i) = 0.0d0
         dev(3,i) = 0.0d0
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (iv14(n))
      allocate (xred(n))
      allocate (yred(n))
      allocate (zred(n))
      allocate (vscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         vscale(i) = 1.0d0
         iv14(i) = 0
      end do
c
c     set the coefficients for the switching function
c
      mode = 'VDW'
      call switch (mode)
c
c     apply any reduction factor to the atomic coordinates
c
      do k = 1, nvdw
         i = ivdw(k)
         iv = ired(i)
         rdn = kred(i)
         xred(i) = rdn*(x(i)-x(iv)) + x(iv)
         yred(i) = rdn*(y(i)-y(iv)) + y(iv)
         zred(i) = rdn*(z(i)-z(iv)) + z(iv)
      end do
c
c     find van der Waals energy and derivatives via double loop
c
      do ii = 1, nvdw-1
         i = ivdw(ii)
         iv = ired(i)
         redi = kred(i)
         rediv = 1.0d0 - redi
         it = jvdw(i)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         usei = (use(i) .or. use(iv))
         muti = mut(i)
c
c     set interaction scaling coefficients for connected atoms
c
         do j = 1, n12(i)
            vscale(i12(j,i)) = v2scale
         end do
         do j = 1, n13(i)
            vscale(i13(j,i)) = v3scale
         end do
         do j = 1, n14(i)
            vscale(i14(j,i)) = v4scale
            iv14(i14(j,i)) = i
         end do
         do j = 1, n15(i)
            vscale(i15(j,i)) = v5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kk = ii+1, nvdw
            k = ivdw(kk)
            kv = ired(k)
            mutk = mut(k)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k) .or. use(kv))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               kt = jvdw(k)
               xr = xi - xred(k)
               yr = yi - yred(k)
               zr = zi - zred(k)
               call image (xr,yr,zr)
               rik2 = xr*xr + yr*yr + zr*zr
c
c     check for an interaction distance less than the cutoff
c
               if (rik2 .le. off2) then
                  rik = sqrt(rik2)
                  rv = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (iv14(k) .eq. i) then
                     rv = radmin4(kt,it)
                     eps = epsilon4(kt,it)
                  end if
                  eps = eps * vscale(k)
c
c     get the energy and gradient, via soft core if necessary
c
                  if ((muti .and. .not.mutk) .or.
     &                (mutk .and. .not.muti)) then
                     rho = rik / rv
                     rho6 = rho**6
                     rho7 = rho6 * rho
                     eps = eps * vlambda**scexp
                     scal = scalpha * (1.0d0-vlambda)**2
                     s1 = 1.0d0 / (scal+(rho+dhal)**7)
                     s2 = 1.0d0 / (scal+rho7+ghal)
                     t1 = (1.0d0+dhal)**7 * s1
                     t2 = (1.0d0+ghal) * s2
                     dt1drho = -7.0d0*(rho+dhal)**6 * t1 * s1
                     dt2drho = -7.0d0*rho6 * t2 * s2
                     e = eps * t1 * (t2-2.0d0)
                     de = eps * (dt1drho*(t2-2.0d0)+t1*dt2drho) / rv
                  else
                     rv7 = rv**7
                     rik6 = rik2**3
                     rik7 = rik6 * rik
                     rho = rik7 + ghal*rv7
                     tau = (dhal+1.0d0) / (rik + dhal*rv)
                     tau7 = tau**7
                     dtau = tau / (dhal+1.0d0)
                     gtau = eps*tau7*rik6*(ghal+1.0d0)*(rv7/rho)**2
                     e = eps*tau7*rv7*((ghal+1.0d0)*rv7/rho-2.0d0)
                     de = -7.0d0 * (dtau*e+gtau)
                  end if
c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. cut2) then
                     rik3 = rik2 * rik
                     rik4 = rik2 * rik2
                     rik5 = rik2 * rik3
                     taper = c5*rik5 + c4*rik4 + c3*rik3
     &                          + c2*rik2 + c1*rik + c0
                     dtaper = 5.0d0*c5*rik4 + 4.0d0*c4*rik3
     &                           + 3.0d0*c3*rik2 + 2.0d0*c2*rik + c1
                     de = e*dtaper + de*taper
                     e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     e = e * fgrp
                     de = de * fgrp
                  end if
c
c     find the chain rule terms for derivative components
c
                  de = de / rik
                  dedx = de * xr
                  dedy = de * yr
                  dedz = de * zr
c
c     increment the total van der Waals energy and derivatives
c
                  ev = ev + e
                  if (i .eq. iv) then
                     dev(1,i) = dev(1,i) + dedx
                     dev(2,i) = dev(2,i) + dedy
                     dev(3,i) = dev(3,i) + dedz
                  else
                     dev(1,i) = dev(1,i) + dedx*redi
                     dev(2,i) = dev(2,i) + dedy*redi
                     dev(3,i) = dev(3,i) + dedz*redi
                     dev(1,iv) = dev(1,iv) + dedx*rediv
                     dev(2,iv) = dev(2,iv) + dedy*rediv
                     dev(3,iv) = dev(3,iv) + dedz*rediv
                  end if
                  if (k .eq. kv) then
                     dev(1,k) = dev(1,k) - dedx
                     dev(2,k) = dev(2,k) - dedy
                     dev(3,k) = dev(3,k) - dedz
                  else
                     redk = kred(k)
                     redkv = 1.0d0 - redk
                     dev(1,k) = dev(1,k) - dedx*redk
                     dev(2,k) = dev(2,k) - dedy*redk
                     dev(3,k) = dev(3,k) - dedz*redk
                     dev(1,kv) = dev(1,kv) - dedx*redkv
                     dev(2,kv) = dev(2,kv) - dedy*redkv
                     dev(3,kv) = dev(3,kv) - dedz*redkv
                  end if
c
c     increment the internal virial tensor components
c
                  vxx = xr * dedx
                  vyx = yr * dedx
                  vzx = zr * dedx
                  vyy = yr * dedy
                  vzy = zr * dedy
                  vzz = zr * dedz
                  vir(1,1) = vir(1,1) + vxx
                  vir(2,1) = vir(2,1) + vyx
                  vir(3,1) = vir(3,1) + vzx
                  vir(1,2) = vir(1,2) + vyx
                  vir(2,2) = vir(2,2) + vyy
                  vir(3,2) = vir(3,2) + vzy
                  vir(1,3) = vir(1,3) + vzx
                  vir(2,3) = vir(2,3) + vzy
                  vir(3,3) = vir(3,3) + vzz
c
c     increment the total intermolecular energy
c
                  if (molcule(i) .ne. molcule(k)) then
                     einter = einter + e
                  end if
               end if
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(i)
            vscale(i12(j,i)) = 1.0d0
         end do
         do j = 1, n13(i)
            vscale(i13(j,i)) = 1.0d0
         end do
         do j = 1, n14(i)
            vscale(i14(j,i)) = 1.0d0
         end do
         do j = 1, n15(i)
            vscale(i15(j,i)) = 1.0d0
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate interaction energy with other unit cells
c
      do ii = 1, nvdw
         i = ivdw(ii)
         iv = ired(i)
         redi = kred(i)
         rediv = 1.0d0 - redi
         it = jvdw(i)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         usei = (use(i) .or. use(iv))
         muti = mut(i)
c
c     set interaction scaling coefficients for connected atoms
c
         do j = 1, n12(i)
            vscale(i12(j,i)) = v2scale
         end do
         do j = 1, n13(i)
            vscale(i13(j,i)) = v3scale
         end do
         do j = 1, n14(i)
            vscale(i14(j,i)) = v4scale
            iv14(i14(j,i)) = i
         end do
         do j = 1, n15(i)
            vscale(i15(j,i)) = v5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kk = ii, nvdw
            k = ivdw(kk)
            kv = ired(k)
            mutk = mut(k)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k) .or. use(kv))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               kt = jvdw(k)
               do j = 1, ncell
                  xr = xi - xred(k)
                  yr = yi - yred(k)
                  zr = zi - zred(k)
                  call imager (xr,yr,zr,j)
                  rik2 = xr*xr + yr*yr + zr*zr
c
c     check for an interaction distance less than the cutoff
c
                  if (rik2 .le. off2) then
                     rik = sqrt(rik2)
                     rv = radmin(kt,it)
                     eps = epsilon(kt,it)
                     if (use_polymer) then
                        if (rik2 .le. polycut2) then
                           if (iv14(k) .eq. i) then
                              rv = radmin4(kt,it)
                              eps = epsilon4(kt,it)
                           end if
                           eps = eps * vscale(k)
                        end if
                     end if
c
c     get the energy and gradient, via soft core if necessary
c
                     if ((muti .and. .not.mutk) .or.
     &                   (mutk .and. .not.muti)) then
                        rho = rik / rv
                        rho6 = rho**6
                        rho7 = rho6 * rho
                        eps = eps * vlambda**scexp
                        scal = scalpha * (1.0d0-vlambda)**2
                        s1 = 1.0d0 / (scal+(rho+dhal)**7)
                        s2 = 1.0d0 / (scal+rho7+ghal)
                        t1 = (1.0d0+dhal)**7 * s1
                        t2 = (1.0d0+ghal) * s2
                        dt1drho = -7.0d0*(rho+dhal)**6 * t1 * s1
                        dt2drho = -7.0d0*rho6 * t2 * s2
                        e = eps * t1 * (t2-2.0d0)
                        de = eps * (dt1drho*(t2-2.0d0)+t1*dt2drho) / rv
                     else
                        rv7 = rv**7
                        rik6 = rik2**3
                        rik7 = rik6 * rik
                        rho = rik7 + ghal*rv7
                        tau = (dhal+1.0d0) / (rik + dhal*rv)
                        tau7 = tau**7
                        dtau = tau / (dhal+1.0d0)
                        gtau = eps*tau7*rik6*(ghal+1.0d0)*(rv7/rho)**2
                        e = eps*tau7*rv7*((ghal+1.0d0)*rv7/rho-2.0d0)
                        de = -7.0d0 * (dtau*e+gtau)
                     end if
c
c     use energy switching if near the cutoff distance
c
                     if (rik2 .gt. cut2) then
                        rik3 = rik2 * rik
                        rik4 = rik2 * rik2
                        rik5 = rik2 * rik3
                        taper = c5*rik5 + c4*rik4 + c3*rik3
     &                             + c2*rik2 + c1*rik + c0
                        dtaper = 5.0d0*c5*rik4 + 4.0d0*c4*rik3
     &                              + 3.0d0*c3*rik2 + 2.0d0*c2*rik + c1
                        de = e*dtaper + de*taper
                        e = e * taper
                     end if
c
c     scale the interaction based on its group membership
c
                     if (use_group) then
                        e = e * fgrp
                        de = de * fgrp
                     end if
c
c     find the chain rule terms for derivative components
c
                     de = de / rik
                     dedx = de * xr
                     dedy = de * yr
                     dedz = de * zr
c
c     increment the total van der Waals energy and derivatives
c
                     if (i .eq. k)  e = 0.5d0 * e
                     ev = ev + e
                     if (i .eq. iv) then
                        dev(1,i) = dev(1,i) + dedx
                        dev(2,i) = dev(2,i) + dedy
                        dev(3,i) = dev(3,i) + dedz
                     else
                        dev(1,i) = dev(1,i) + dedx*redi
                        dev(2,i) = dev(2,i) + dedy*redi
                        dev(3,i) = dev(3,i) + dedz*redi
                        dev(1,iv) = dev(1,iv) + dedx*rediv
                        dev(2,iv) = dev(2,iv) + dedy*rediv
                        dev(3,iv) = dev(3,iv) + dedz*rediv
                     end if
                     if (i .ne. k) then
                        if (k .eq. kv) then
                           dev(1,k) = dev(1,k) - dedx
                           dev(2,k) = dev(2,k) - dedy
                           dev(3,k) = dev(3,k) - dedz
                        else
                           redk = kred(k)
                           redkv = 1.0d0 - redk
                           dev(1,k) = dev(1,k) - dedx*redk
                           dev(2,k) = dev(2,k) - dedy*redk
                           dev(3,k) = dev(3,k) - dedz*redk
                           dev(1,kv) = dev(1,kv) - dedx*redkv
                           dev(2,kv) = dev(2,kv) - dedy*redkv
                           dev(3,kv) = dev(3,kv) - dedz*redkv
                        end if
                     end if
c
c     increment the internal virial tensor components
c
                     vxx = xr * dedx
                     vyx = yr * dedx
                     vzx = zr * dedx
                     vyy = yr * dedy
                     vzy = zr * dedy
                     vzz = zr * dedz
                     vir(1,1) = vir(1,1) + vxx
                     vir(2,1) = vir(2,1) + vyx
                     vir(3,1) = vir(3,1) + vzx
                     vir(1,2) = vir(1,2) + vyx
                     vir(2,2) = vir(2,2) + vyy
                     vir(3,2) = vir(3,2) + vzy
                     vir(1,3) = vir(1,3) + vzx
                     vir(2,3) = vir(2,3) + vzy
                     vir(3,3) = vir(3,3) + vzz
c
c     increment the total intermolecular energy
c
                     einter = einter + e
                  end if
               end do
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(i)
            vscale(i12(j,i)) = 1.0d0
         end do
         do j = 1, n13(i)
            vscale(i13(j,i)) = 1.0d0
         end do
         do j = 1, n14(i)
            vscale(i14(j,i)) = 1.0d0
         end do
         do j = 1, n15(i)
            vscale(i15(j,i)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (iv14)
      deallocate (xred)
      deallocate (yred)
      deallocate (zred)
      deallocate (vscale)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine ehal1b  --  buffered 14-7 vdw derivs via lights  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "ehal1b" calculates the buffered 14-7 van der Waals energy and
c     its first derivatives with respect to Cartesian coordinates
c     using the method of lights
c
c
      subroutine ehal1b
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bound.i'
      include 'boxes.i'
      include 'cell.i'
      include 'couple.i'
      include 'deriv.i'
      include 'energi.i'
      include 'group.i'
      include 'inter.i'
      include 'iounit.i'
      include 'light.i'
      include 'molcul.i'
      include 'mutant.i'
      include 'shunt.i'
      include 'usage.i'
      include 'vdw.i'
      include 'vdwpot.i'
      include 'virial.i'
      integer i,j,k
      integer ii,iv,it
      integer kk,kv,kt
      integer kgy,kgz
      integer start,stop
      integer, allocatable :: iv14(:)
      real*8 e,de,eps,rdn
      real*8 fgrp,rv,rv7
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 redi,rediv
      real*8 redk,redkv
      real*8 dedx,dedy,dedz
      real*8 rho,rho6,rho7
      real*8 tau,tau7,scal
      real*8 s1,s2,t1,t2
      real*8 dt1drho,dt2drho
      real*8 dtau,gtau
      real*8 taper,dtaper
      real*8 rik,rik2,rik3
      real*8 rik4,rik5
      real*8 rik6,rik7
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8, allocatable :: xred(:)
      real*8, allocatable :: yred(:)
      real*8, allocatable :: zred(:)
      real*8, allocatable :: vscale(:)
      real*8, allocatable :: xsort(:)
      real*8, allocatable :: ysort(:)
      real*8, allocatable :: zsort(:)
      logical proceed,usei
      logical muti,mutk
      logical prime,repeat
      character*6 mode
c
c
c     zero out the van der Waals energy and first derivatives
c
      ev = 0.0d0
      do i = 1, n
         dev(1,i) = 0.0d0
         dev(2,i) = 0.0d0
         dev(3,i) = 0.0d0
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (iv14(n))
      allocate (xred(n))
      allocate (yred(n))
      allocate (zred(n))
      allocate (vscale(n))
      allocate (xsort(8*n))
      allocate (ysort(8*n))
      allocate (zsort(8*n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         vscale(i) = 1.0d0
         iv14(i) = 0
      end do
c
c     set the coefficients for the switching function
c
      mode = 'VDW'
      call switch (mode)
c
c     apply any reduction factor to the atomic coordinates
c
      do j = 1, nvdw
         i = ivdw(j)
         iv = ired(i)
         rdn = kred(i)
         xred(j) = rdn*(x(i)-x(iv)) + x(iv)
         yred(j) = rdn*(y(i)-y(iv)) + y(iv)
         zred(j) = rdn*(z(i)-z(iv)) + z(iv)
      end do
c
c     transfer the interaction site coordinates to sorting arrays
c
      do i = 1, nvdw
         xsort(i) = xred(i)
         ysort(i) = yred(i)
         zsort(i) = zred(i)
      end do
c
c     use the method of lights to generate neighbors
c
      call lights (off,nvdw,xsort,ysort,zsort)
c
c     loop over all atoms computing the interactions
c
      do ii = 1, nvdw
         i = ivdw(ii)
         iv = ired(i)
         redi = kred(i)
         rediv = 1.0d0 - redi
         it = jvdw(i)
         xi = xsort(rgx(ii))
         yi = ysort(rgy(ii))
         zi = zsort(rgz(ii))
         usei = (use(i) .or. use(iv))
         muti = mut(i)
c
c     set interaction scaling coefficients for connected atoms
c
         do j = 1, n12(i)
            vscale(i12(j,i)) = v2scale
         end do
         do j = 1, n13(i)
            vscale(i13(j,i)) = v3scale
         end do
         do j = 1, n14(i)
            vscale(i14(j,i)) = v4scale
            iv14(i14(j,i)) = i
         end do
         do j = 1, n15(i)
            vscale(i15(j,i)) = v5scale
         end do
c
c     loop over method of lights neighbors of current atom
c
         if (kbx(ii) .le. kex(ii)) then
            repeat = .false.
            start = kbx(ii) + 1
            stop = kex(ii)
         else
            repeat = .true.
            start = 1
            stop = kex(ii)
         end if
   10    continue
         do j = start, stop
            kk = locx(j)
            kgy = rgy(kk)
            if (kby(ii) .le. key(ii)) then
               if (kgy.lt.kby(ii) .or. kgy.gt.key(ii))  goto 20
            else
               if (kgy.lt.kby(ii) .and. kgy.gt.key(ii))  goto 20
            end if
            kgz = rgz(kk)
            if (kbz(ii) .le. kez(ii)) then
               if (kgz.lt.kbz(ii) .or. kgz.gt.kez(ii))  goto 20
            else
               if (kgz.lt.kbz(ii) .and. kgz.gt.kez(ii))  goto 20
            end if
            k = ivdw(kk-((kk-1)/nvdw)*nvdw)
            kv = ired(k)
            mutk = mut(k)
            prime = (kk .le. nvdw)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k) .or. use(kv))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               kt = jvdw(k)
               xr = xi - xsort(j)
               yr = yi - ysort(kgy)
               zr = zi - zsort(kgz)
               if (use_bounds) then
                  if (abs(xr) .gt. xcell2)  xr = xr - sign(xcell,xr)
                  if (abs(yr) .gt. ycell2)  yr = yr - sign(ycell,yr)
                  if (abs(zr) .gt. zcell2)  zr = zr - sign(zcell,zr)
                  if (monoclinic) then
                     xr = xr + zr*beta_cos
                     zr = zr * beta_sin
                  else if (triclinic) then
                     xr = xr + yr*gamma_cos + zr*beta_cos
                     yr = yr*gamma_sin + zr*beta_term
                     zr = zr * gamma_term
                  end if
               end if
               rik2 = xr*xr + yr*yr + zr*zr
c
c     check for an interaction distance less than the cutoff
c
               if (rik2 .le. off2) then
                  rik = sqrt(rik2)
                  rv = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (prime) then
                     if (iv14(k) .eq. i) then
                        rv = radmin4(kt,it)
                        eps = epsilon4(kt,it)
                     end if
                     eps = eps * vscale(k)
                  end if
c
c     get the energy and gradient, via soft core if necessary
c
                  if ((muti .and. .not.mutk) .or.
     &                (mutk .and. .not.muti)) then
                     rho = rik / rv
                     rho6 = rho**6
                     rho7 = rho6 * rho
                     eps = eps * vlambda**scexp
                     scal = scalpha * (1.0d0-vlambda)**2
                     s1 = 1.0d0 / (scal+(rho+dhal)**7)
                     s2 = 1.0d0 / (scal+rho7+ghal)
                     t1 = (1.0d0+dhal)**7 * s1
                     t2 = (1.0d0+ghal) * s2
                     dt1drho = -7.0d0*(rho+dhal)**6 * t1 * s1
                     dt2drho = -7.0d0*rho6 * t2 * s2
                     e = eps * t1 * (t2-2.0d0)
                     de = eps * (dt1drho*(t2-2.0d0)+t1*dt2drho) / rv
                  else
                     rv7 = rv**7
                     rik6 = rik2**3
                     rik7 = rik6 * rik
                     rho = rik7 + ghal*rv7
                     tau = (dhal+1.0d0) / (rik + dhal*rv)
                     tau7 = tau**7
                     dtau = tau / (dhal+1.0d0)
                     gtau = eps*tau7*rik6*(ghal+1.0d0)*(rv7/rho)**2
                     e = eps*tau7*rv7*((ghal+1.0d0)*rv7/rho-2.0d0)
                     de = -7.0d0 * (dtau*e+gtau)
                  end if
c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. cut2) then
                     rik3 = rik2 * rik
                     rik4 = rik2 * rik2
                     rik5 = rik2 * rik3
                     taper = c5*rik5 + c4*rik4 + c3*rik3
     &                          + c2*rik2 + c1*rik + c0
                     dtaper = 5.0d0*c5*rik4 + 4.0d0*c4*rik3
     &                           + 3.0d0*c3*rik2 + 2.0d0*c2*rik + c1
                     de = e*dtaper + de*taper
                     e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     e = e * fgrp
                     de = de * fgrp
                  end if
c
c     find the chain rule terms for derivative components
c
                  de = de / rik
                  dedx = de * xr
                  dedy = de * yr
                  dedz = de * zr
c
c     increment the total van der Waals energy and derivatives
c
                  ev = ev + e
                  if (i .eq. iv) then
                     dev(1,i) = dev(1,i) + dedx
                     dev(2,i) = dev(2,i) + dedy
                     dev(3,i) = dev(3,i) + dedz
                  else
                     dev(1,i) = dev(1,i) + dedx*redi
                     dev(2,i) = dev(2,i) + dedy*redi
                     dev(3,i) = dev(3,i) + dedz*redi
                     dev(1,iv) = dev(1,iv) + dedx*rediv
                     dev(2,iv) = dev(2,iv) + dedy*rediv
                     dev(3,iv) = dev(3,iv) + dedz*rediv
                  end if
                  if (k .eq. kv) then
                     dev(1,k) = dev(1,k) - dedx
                     dev(2,k) = dev(2,k) - dedy
                     dev(3,k) = dev(3,k) - dedz
                  else
                     redk = kred(k)
                     redkv = 1.0d0 - redk
                     dev(1,k) = dev(1,k) - dedx*redk
                     dev(2,k) = dev(2,k) - dedy*redk
                     dev(3,k) = dev(3,k) - dedz*redk
                     dev(1,kv) = dev(1,kv) - dedx*redkv
                     dev(2,kv) = dev(2,kv) - dedy*redkv
                     dev(3,kv) = dev(3,kv) - dedz*redkv
                  end if
c
c     increment the internal virial tensor components
c
                  vxx = xr * dedx
                  vyx = yr * dedx
                  vzx = zr * dedx
                  vyy = yr * dedy
                  vzy = zr * dedy
                  vzz = zr * dedz
                  vir(1,1) = vir(1,1) + vxx
                  vir(2,1) = vir(2,1) + vyx
                  vir(3,1) = vir(3,1) + vzx
                  vir(1,2) = vir(1,2) + vyx
                  vir(2,2) = vir(2,2) + vyy
                  vir(3,2) = vir(3,2) + vzy
                  vir(1,3) = vir(1,3) + vzx
                  vir(2,3) = vir(2,3) + vzy
                  vir(3,3) = vir(3,3) + vzz
c
c     increment the total intermolecular energy
c
                  if (.not.prime .or. molcule(i).ne.molcule(k)) then
                     einter = einter + e
                  end if
               end if
            end if
   20       continue
         end do
         if (repeat) then
            repeat = .false.
            start = kbx(ii) + 1
            stop = nlight
            goto 10
         end if
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(i)
            vscale(i12(j,i)) = 1.0d0
         end do
         do j = 1, n13(i)
            vscale(i13(j,i)) = 1.0d0
         end do
         do j = 1, n14(i)
            vscale(i14(j,i)) = 1.0d0
         end do
         do j = 1, n15(i)
            vscale(i15(j,i)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (iv14)
      deallocate (xred)
      deallocate (yred)
      deallocate (zred)
      deallocate (vscale)
      deallocate (xsort)
      deallocate (ysort)
      deallocate (zsort)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ehal1c  --  buffered 14-7 vdw derivs via list  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ehal1c" calculates the buffered 14-7 van der Waals energy and
c     its first derivatives with respect to Cartesian coordinates
c     using a pairwise neighbor list
c
c
      subroutine ehal1c
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bound.i'
      include 'couple.i'
      include 'deriv.i'
      include 'energi.i'
      include 'group.i'
      include 'inter.i'
      include 'molcul.i'
      include 'mutant.i'
      include 'neigh.i'
      include 'shunt.i'
      include 'usage.i'
      include 'vdw.i'
      include 'vdwpot.i'
      include 'virial.i'
      integer i,j,k,indexmut
      integer ii,iv,it
      integer kk,kv,kt
      integer, allocatable :: iv14(:)
      real*8 e,de,eps,rdn
      real*8 fgrp,rv,rv7
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 redi,rediv
      real*8 redk,redkv
      real*8 dedx,dedy,dedz
      real*8 rho,rho6,rho7
      real*8 tau,tau7,scal
      real*8 s1,s2,t1,t2
      real*8 dt1drho,dt2drho
      real*8 dtau,gtau
      real*8 taper,dtaper
      real*8 rik,rik2,rik3
      real*8 rik4,rik5
      real*8 rik6,rik7
c     JRA additional variables
      real*8 vlamij, dedlsignv
      real*8 g1,g2,g3
      real*8 t11,t12,t13,t21,t22,t23
      real*8 dedxdl
      real*8 drhodr
      real*8 l_n,dl_ndl,dscaldl,dt1dl
      real*8 dleftdl, dt2dl, dedl_temp
      real*8 d2l_ndl2, d2scaldl2, d2t1dl2
      real*8 d2leftdl2, d2t2dl2, d2edl2_temp
      real*8 d2t1dldrho, d2leftdldrho, d2t2dldrho
      real*8 dedldrho, dedldr
      real*8 dedldx, dedldy, dedldz
      real*8 dedlvt, d2edl2vt
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8 evt,eintert
      real*8 virt(3,3)
      real*8, allocatable :: xred(:)
      real*8, allocatable :: yred(:)
      real*8, allocatable :: zred(:)
      real*8, allocatable :: vscale(:)
      real*8, allocatable :: devt(:,:)
      real*8, allocatable :: d2edlvt(:,:)
      logical proceed,usei
      logical muti,mutk,soft
      logical ifrom1,kfrom1,ifrom2,kfrom2
      character*6 mode
c
c
c     zero out the van der Waals energy and first derivatives
c
      ev = 0.0d0
      do i = 1, n
         dev(1,i) = 0.0d0
         dev(2,i) = 0.0d0
         dev(3,i) = 0.0d0
      end do
      
c
c     zero out vdw lambda derivatives JRA
c
      dedlv = 0.0d0
      d2edl2v = 0.0d0
      do i = 1, n
         d2edlv(1,i) = 0.0d0
         d2edlv(2,i) = 0.0d0
         d2edlv(3,i) = 0.0d0
      end do      
      
c
c     perform dynamic allocation of some local arrays
c         JRA additional array
c
      allocate (iv14(n))
      allocate (xred(n))
      allocate (yred(n))
      allocate (zred(n))
      allocate (vscale(n))
      allocate (devt(3,n))
      allocate (d2edlvt(3,n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         vscale(i) = 1.0d0
         iv14(i) = 0
      end do
c
c     set the coefficients for the switching function
c
      mode = 'VDW'
      call switch (mode)
c
c     apply any reduction factor to the atomic coordinates
c
      do k = 1, nvdw
         i = ivdw(k)
         iv = ired(i)
         rdn = kred(i)
         xred(i) = rdn*(x(i)-x(iv)) + x(iv)
         yred(i) = rdn*(y(i)-y(iv)) + y(iv)
         zred(i) = rdn*(z(i)-z(iv)) + z(iv)
      end do
c
c     transfer global to local copies for OpenMP calculation
c
      evt = ev
      eintert = einter
      do i = 1, n
         devt(1,i) = dev(1,i)
         devt(2,i) = dev(2,i)
         devt(3,i) = dev(3,i)
      end do
      do i = 1, 3
         virt(1,i) = vir(1,i)
         virt(2,i) = vir(2,i)
         virt(3,i) = vir(3,i)
      end do 
      dedlvt = dedlv
      d2edl2vt = d2edl2v
      do i = 1, n
         d2edlvt(1,i) = d2edlv(1,i)
         d2edlvt(2,i) = d2edlv(2,i)
         d2edlvt(3,i) = d2edlv(3,i)
      end do
c
c     set OpenMP directives for the major loop structure
c         JRA added additional variables
c
!$OMP PARALLEL default(private) shared(nvdw,ivdw,ired,kred,
!$OMP& jvdw,xred,yred,zred,use,nvlst,vlst,n12,n13,n14,n15,
!$OMP& i12,i13,i14,i15,v2scale,v3scale,v4scale,v5scale,
!$OMP& use_group,off2,radmin,epsilon,radmin4,epsilon4,ghal,dhal,
!$OMP& osrwon,isrelative,vlambda,nmut1,nmut2,
!$OMP& cut2,scalphav,scexpv,mut,c0,c1,c2,c3,c4,c5,molcule)
!$OMP& firstprivate(vscale,iv14) shared(evt,devt,virt,eintert,
!$OMP& dedlvt,d2edl2vt,d2edlvt)
!$OMP DO reduction(+:evt,devt,virt,eintert,
!$OMP& dedlvt,d2edl2vt,d2edlvt) schedule(guided)
c
c     find van der Waals energy and derivatives via neighbor list
c
      do ii = 1, nvdw
         i = ivdw(ii)
         iv = ired(i)
         redi = kred(i)
         rediv = 1.0d0 - redi
         it = jvdw(i)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         usei = (use(i) .or. use(iv))
         muti = mut(i)
c
c     set interaction scaling coefficients for connected atoms
c
         do j = 1, n12(i)
            vscale(i12(j,i)) = v2scale
         end do
         do j = 1, n13(i)
            vscale(i13(j,i)) = v3scale
         end do
         do j = 1, n14(i)
            vscale(i14(j,i)) = v4scale
            iv14(i14(j,i)) = i
         end do
         do j = 1, n15(i)
            vscale(i15(j,i)) = v5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kk = 1, nvlst(ii)
            k = ivdw(vlst(kk,ii))
            kv = ired(k)
            mutk = mut(k)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k) .or. use(kv))
            
c         JRA - Ligand 1 and ligand 2 do not interact with each other
c         ifrom1 is referring to i, kfrom1 referring to k
            ifrom1 = .false.
            kfrom1 = .false.
            ifrom2 = .false.
            kfrom2 = .false.
            if (osrwon .and .isrelative) then
c          JRA - assumes mut atoms are on top and in order (1 then 2)
              ifrom1 = i .gt. 0 .and. i .le. nmut1
              kfrom1 = k .gt. 0 .and. k .le. nmut1
              ifrom2 = i .gt. nmut1+1 .and. i .le. nmut1+nmut2
              kfrom2 = k .gt. nmut1+1 .and. k .le. nmut1+nmut2
              if ((ifrom1 .and. kfrom2) .or.
     &            (ifrom2 .and. kfrom1)) then
                proceed = .false.
              end if
            end if            
            
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               kt = jvdw(k)
               xr = xi - xred(k)
               yr = yi - yred(k)
               zr = zi - zred(k)
               call image (xr,yr,zr)
               rik2 = xr*xr + yr*yr + zr*zr
c
c     check for an interaction distance less than the cutoff
c
               if (rik2 .le. off2) then
                  rik = sqrt(rik2)
                  rv = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (iv14(k) .eq. i) then
                     rv = radmin4(kt,it)
                     eps = epsilon4(kt,it)
                  end if
                  eps = eps * vscale(k)
c
c     get the energy and gradient, via soft core if necessary
c           JRA altered gradient calculation slightly
c
                  l_n = 1.0d0
                  soft = (muti .and. .not.mutk) .or.
     &                (mutk .and. .not.muti)
c         JRA - If mutating entity is from system 1, use (1-lambda)
                  vlamij = 1.0d0
                  dedlsignv = 1.0d0
                  if (osrwon .and. soft) then
                    vlamij = vlambda
                    if (isrelative .and. (ifrom1 .or. kfrom1)) then
                       vlamij = 1.0d0 - vlambda
                       dedlsignv = -1.0d0                     
                    end if
                  end if
                  if (soft) then
                     rho = rik / rv
                     rho6 = rho**6
                     rho7 = rho6 * rho
                     
                     l_n = vlamij**scexpv
c                     eps = eps * vlambda**scexpv
                     scal = scalphav * (1.0d0-vlamij)**2
                     s1 = 1.0d0 / (scal+(rho+dhal)**7)
                     s2 = 1.0d0 / (scal+rho7+ghal)
                     t1 = (1.0d0+dhal)**7 * s1
                     t2 = (1.0d0+ghal) * s2
                     dt1drho = -7.0d0*(rho+dhal)**6 * t1 * s1
                     dt2drho = -7.0d0*rho6 * t2 * s2
c                     e = eps * t1 * (t2-2.0d0)
c                     de = eps * (dt1drho*(t2-2.0d0)+t1*dt2drho) / rv
                     drhodr = 1.0d0/rv
                     e = eps*l_n*t1*(t2-2.0d0)
                     de = eps*l_n*(dt1drho*(t2-2.0d0)+t1*dt2drho)*drhodr
                  else
                     rv7 = rv**7
                     rik6 = rik2**3
                     rik7 = rik6 * rik
                     rho = rik7 + ghal*rv7
                     tau = (dhal+1.0d0) / (rik + dhal*rv)
                     tau7 = tau**7
                     dtau = tau / (dhal+1.0d0)
                     gtau = eps*tau7*rik6*(ghal+1.0d0)*(rv7/rho)**2
                     e = eps*tau7*rv7*((ghal+1.0d0)*rv7/rho-2.0d0)
                     de = -7.0d0 * (dtau*e+gtau)
                  end if
c
c     use energy switching if near the cutoff distance
c          JRA change to keep taper and dtaper
c
                  taper = 1.0d0
                  dtaper = 0.0d0
                  if (rik2 .gt. cut2) then
                     rik3 = rik2 * rik
                     rik4 = rik2 * rik2
                     rik5 = rik2 * rik3
                     taper = c5*rik5 + c4*rik4 + c3*rik3
     &                          + c2*rik2 + c1*rik + c0
                     dtaper = 5.0d0*c5*rik4 + 4.0d0*c4*rik3
     &                           + 3.0d0*c3*rik2 + 2.0d0*c2*rik + c1
c                     de = e*dtaper + de*taper
c                     e = e * taper
                  end if
                  de = e*dtaper + de*taper
                  e = e * taper
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     e = e * fgrp
                     de = de * fgrp
                  end if
                  
c
c     JRA calculate derivatives wrt lambda
c        
                  if (soft) then
                     dl_ndl = 1.0d0*dedlsignv
                     if (scexpv .gt. 1.5d0) then
                       dl_ndl = dedlsignv*scexpv*vlamij**(scexpv-1.0d0)
                     end if
                     dscaldl = -dedlsignv*2.0d0*scalphav*(1-vlamij)
                     dt1dl = -dscaldl*t1*s1
                     dleftdl = l_n*dt1dl + dl_ndl*t1
                     dt2dl = -dscaldl*t2*s2
                     g1 = dl_ndl * t1 * (t2 - 2.0d0)
                     g2 = l_n * dt1dl * (t2 - 2.0d0)
                     g3 = l_n * t1 * dt2dl
                     dedlvt = dedlvt + eps * (g1+g2+g3) * taper
c                     dedl_temp = l_n*t1*dt2dl + dleftdl*(t2-2.0d0)
                     dedl_temp = g1+g2+g3
                     if (use_group) then
                        dedl_temp = dedl_temp * fgrp
                     end if
c                     dedlvt = dedlvt + eps*dedl_temp*taper
                     
                     if (scexpv-2.0d0 .lt. 0.0d0) then
                        d2l_ndl2 = 0.0d0
                     else
                        d2l_ndl2 = scexpv*(scexpv-1.0d0)
     &                                  *vlamij**(scexpv-2.0d0)
                     end if 
                     d2scaldl2 = 2.0d0*scalphav
                     d2t1dl2 = -s1*(2.0d0*dt1dl*dscaldl + t1*d2scaldl2)
                     d2leftdl2 = l_n*d2t1dl2 + dl_ndl*dt1dl
     &                          + dl_ndl*dt1dl + d2l_ndl2*t1
                     d2t2dl2 = -s2*(2.0d0*dt2dl*dscaldl + t2*d2scaldl2)
                     d2edl2_temp = l_n*t1*d2t2dl2 + dleftdl*dt2dl
     &                         + dleftdl*dt2dl + d2leftdl2*(t2-2.0d0)
                     if (use_group) then
                        d2edl2_temp = d2edl2_temp * fgrp
                     end if
                     d2edl2vt = d2edl2vt + eps*d2edl2_temp*taper
                     
c                     d2t1dldrho = -14.0d0*dt1dl*s1*(rho+dhal)**6
c                     d2leftdldrho = l_n*d2t1dldrho + dl_ndl*dt1drho
c                     d2t2dldrho = -14.0d0*dt2dl*s2*rho6
c                     dedldrho = l_n*(t1*d2t2dldrho + dt1drho*dt2dl)
c     &                   + d2leftdldrho*(t2-2.0d0) + dleftdl*dt2drho
c                     dedldr = eps * dedldrho * drhodr
                     t11  = -dl_ndl*(t2-2.0d0)*(-dt1drho)*drhodr
                     t12 = -l_n * dt2dl * (-dt1drho)*drhodr
                     t13 = 2.0 * l_n * (t2-2.0d0) * (-dt1drho)*drhodr
     &                                         * dscaldl * s1
                     t21 = -dl_ndl * t1 * (-dt2drho)*drhodr
                     t22 = -l_n * dt1dl * (-dt2drho)*drhodr
                     t23 = 2.0 * l_n * t1 * (-dt2drho)*drhodr 
     &                                         * dscaldl * s2
                     dedldr = eps * (t11+t12+t13+t21+t22+t23)
                     if (use_group) then
                        dedldr = dedldr * fgrp
                     end if
                     dedldr = eps*dedl_temp*dtaper + dedldr*taper
                     
                     dedldx = dedldr*xr/rik
                     dedldy = dedldr*yr/rik
                     dedldz = dedldr*zr/rik
                     
                     if (i .eq. iv) then
                        d2edlvt(1,i) = d2edlvt(1,i) + dedldx
                        d2edlvt(2,i) = d2edlvt(2,i) + dedldy
                        d2edlvt(3,i) = d2edlvt(3,i) + dedldz
                     else
                        d2edlvt(1,i) = d2edlvt(1,i) + dedldx*redi
                        d2edlvt(2,i) = d2edlvt(2,i) + dedldy*redi
                        d2edlvt(3,i) = d2edlvt(3,i) + dedldz*redi
                        d2edlvt(1,iv) = d2edlvt(1,iv) + dedldx*rediv
                        d2edlvt(2,iv) = d2edlvt(2,iv) + dedldy*rediv
                        d2edlvt(3,iv) = d2edlvt(3,iv) + dedldz*rediv
                     end if
                     if (k .eq. kv) then
                        d2edlvt(1,k) = d2edlvt(1,k) - dedldx
                        d2edlvt(2,k) = d2edlvt(2,k) - dedldy
                        d2edlvt(3,k) = d2edlvt(3,k) - dedldz
                     else
                        redk = kred(k)
                        redkv = 1.0d0 - redk
                        d2edlvt(1,k) = d2edlvt(1,k) - dedldx*redk
                        d2edlvt(2,k) = d2edlvt(2,k) - dedldy*redk
                        d2edlvt(3,k) = d2edlvt(3,k) - dedldz*redk
                        d2edlvt(1,kv) = d2edlvt(1,kv) - dedldx*redkv
                        d2edlvt(2,kv) = d2edlvt(2,kv) - dedldy*redkv
                        d2edlvt(3,kv) = d2edlvt(3,kv) - dedldz*redkv
                     end if
                     
                  end if     

c
c     find the chain rule terms for derivative components
c
                  de = de / rik
                  dedx = de * xr
                  dedy = de * yr
                  dedz = de * zr
c
c     increment the total van der Waals energy and derivatives
c
                  evt = evt + e
                  if (i .eq. iv) then
                     devt(1,i) = devt(1,i) + dedx
                     devt(2,i) = devt(2,i) + dedy
                     devt(3,i) = devt(3,i) + dedz
                  else
                     devt(1,i) = devt(1,i) + dedx*redi
                     devt(2,i) = devt(2,i) + dedy*redi
                     devt(3,i) = devt(3,i) + dedz*redi
                     devt(1,iv) = devt(1,iv) + dedx*rediv
                     devt(2,iv) = devt(2,iv) + dedy*rediv
                     devt(3,iv) = devt(3,iv) + dedz*rediv
                  end if
                  if (k .eq. kv) then
                     devt(1,k) = devt(1,k) - dedx
                     devt(2,k) = devt(2,k) - dedy
                     devt(3,k) = devt(3,k) - dedz
                  else
                     redk = kred(k)
                     redkv = 1.0d0 - redk
                     devt(1,k) = devt(1,k) - dedx*redk
                     devt(2,k) = devt(2,k) - dedy*redk
                     devt(3,k) = devt(3,k) - dedz*redk
                     devt(1,kv) = devt(1,kv) - dedx*redkv
                     devt(2,kv) = devt(2,kv) - dedy*redkv
                     devt(3,kv) = devt(3,kv) - dedz*redkv
                  end if
c
c     increment the internal virial tensor components
c
                  vxx = xr * dedx
                  vyx = yr * dedx
                  vzx = zr * dedx
                  vyy = yr * dedy
                  vzy = zr * dedy
                  vzz = zr * dedz
                  virt(1,1) = virt(1,1) + vxx
                  virt(2,1) = virt(2,1) + vyx
                  virt(3,1) = virt(3,1) + vzx
                  virt(1,2) = virt(1,2) + vyx
                  virt(2,2) = virt(2,2) + vyy
                  virt(3,2) = virt(3,2) + vzy
                  virt(1,3) = virt(1,3) + vzx
                  virt(2,3) = virt(2,3) + vzy
                  virt(3,3) = virt(3,3) + vzz
c
c     increment the total intermolecular energy
c
                  if (molcule(i) .ne. molcule(k)) then
                     eintert = eintert + e
                  end if
               end if
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(i)
            vscale(i12(j,i)) = 1.0d0
         end do
         do j = 1, n13(i)
            vscale(i13(j,i)) = 1.0d0
         end do
         do j = 1, n14(i)
            vscale(i14(j,i)) = 1.0d0
         end do
         do j = 1, n15(i)
            vscale(i15(j,i)) = 1.0d0
         end do
      end do
c
c     end OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
c
c     transfer local to global copies for OpenMP calculation
c
      ev = evt
      einter = eintert
      do i = 1, n
         dev(1,i) = devt(1,i)
         dev(2,i) = devt(2,i)
         dev(3,i) = devt(3,i)
      end do
      do i = 1, 3
         vir(1,i) = virt(1,i)
         vir(2,i) = virt(2,i)
         vir(3,i) = virt(3,i)
      end do
c
c   JRA transfer for lambda derivatives
c
      dedlv = dedlvt
      d2edl2v = d2edl2vt
      do i = 1, n
         d2edlv(1,i) = d2edlvt(1,i)
         d2edlv(2,i) = d2edlvt(2,i)
         d2edlv(3,i) = d2edlvt(3,i)
      end do
c
c     perform deallocation of some local arrays
c        JRA additional array
c
      deallocate (iv14)
      deallocate (xred)
      deallocate (yred)
      deallocate (zred)
      deallocate (vscale)
      deallocate (devt)
      deallocate (d2edlvt)
      return
      end
