c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2014  by  Jay William Ponder  ##
c     ## Michael Schneiders, Jayvee Abella, Pengyu Ren ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  mutant.i  --  hybrid atoms for free energy perturbation  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     lambda     generic weighting between initial and final states
c     vlambda    state weighting value for electrostatic potentials
c     elambda    state weighting value for van der Waals potentials
c     scexp
c     scalpha
c     eslvt
c     eslut
c     nmut       number of atoms mutated from initial to final state
c     imut       atomic sites differing in initial and final state
c     type0      atom type of each atom in the initial state system
c     class0     atom class of each atom in the initial state system
c     type1      atom type of each atom in the final state system
c     class1     atom class of each atom in the final state system
c     mut        true if an atom is to be mutated, false otherwise
c     eupdated   flag to mark updated energy value in BAR method
c
c     softRadCount   counts mutated radius for vdw correction 
c     scexpv         softcore vdw exponent parameter
c     scexpm         softcore empole exponent parameter
c     scalphav       softcore vdw alpha parameter
c     scalpham       softcore empole alpha parameter
c     scexpp         softcore polarization exponent parameter
c     polstart       starting lambda value to begin polarization
c     polend         ending lambda value to end polarization
c     osrwon      flag to signal calculation of osrw/softcore code
c     doNoLigandCondensed   flag to signal calculation of env energy
c     doVaporElec    flag to signal calculation of lig ener in vacuum
c     use_soft       flag to signal calculation of softcore in empole
c     isrelative     flag for calc. of relative vs absolute trans
c     nmut1          number of mutating atoms in system 1
c     nmut2          number of mutating atoms in system 2
c     imut1          atomic sites differing in system 1
c     imut2          atomic sites differing in system 2
c     poleold        original pole values
c     polarityold    original polarity values
c     rel_ligA       flag to use 1-lambda in elec calc
c     debugjv        flag to print energy components in testlambdagrad
c     envenvon       flag to calculate pro-pro or env-env real space...
c                    ...if false, only calculates soft interactions
c     dorealsplit    flag to split real calculation since polar...
c                        ...doesnt use soft
c
c     dedlv          d(ev) / dlambda
c     d2edl2v        d2(ev) / d2lambda
c     d2edlv         spatial gradient of dedlv
c     dedlm          d(em) / dlambda
c     d2edl2m        d2(em) / d2lambda
c     d2edlg         spatial gradient of dedlm
c     d2edlt         NOT USED (related to torque)
c     ereal          real portion of em
c     dereal         real portion of dem
c     dudl           d(e) / dlambda
c     d2udl2         d2(e) / d2lambda
c     d2udl          spatial gradient of dudl
c

      integer nmut,imut,nmut1,nmut2
      integer imut1,imut2
      integer type0,class0
      integer type1,class1
      integer softRadCount
      real*8 lambda
      real*8 vlambda,elambda
      real*8 scexp,scalpha
      real*8 eslvt,eslut
      real*8 scexpv,scexpm,scalphav,scalpham
      real*8 scexpp,polstart,polend
      real*8 dedlv,d2edl2v,d2edlv
      real*8 dedlm,d2edl2m,d2edlg,d2edlt
      real*8 ereal,dereal
      real*8 ereciprel,dereciprel
      real*8 dudl,d2udl2,d2udl
      real*8  poleold,polarityold
      logical mut,eupdated,osrwon,use_soft
      logical doNoLigandCondensed, doVaporElec
      logical isrelative, rel_ligA, debugjv, envenvon, dorealsplit
      common /mutant/ lambda,vlambda,elambda,scexp,scalpha,eslvt,eslut,
     &                scexpv,scexpm,scalphav,scalpham,scexpp,
     &                polstart,polend,dedlv,d2edl2v,d2edlv(3,maxatm),
     &                dedlm,d2edl2m,d2edlg(3,maxatm),d2edlt(3,maxatm),
     &                ereal,dereal(3,maxatm),
     &                dudl,d2udl2,d2udl(3,maxatm),
     &                poleold(13,maxatm),polarityold(maxatm),
     &                nmut,imut(maxatm),nmut1,nmut2,imut1(maxatm),
     &                imut2(maxatm),type0(maxatm),class0(maxatm),
     &                type1(maxatm),class1(maxatm),
     &                softRadCount(maxatm),mut(maxatm),
     &                eupdated,osrwon,use_soft,doNoLigandCondensed,
     &                doVaporElec, isrelative, rel_ligA, debugjv,
     &                envenvon, dorealsplit
