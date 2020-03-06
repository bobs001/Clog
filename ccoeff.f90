! 
! ROUTINE: bps_ccoeff_ab_mass(nni, ep, mp, zp, ia, ib, betab, zb, mb, nb, &
!            a_ab, a_ab_sing, a_ab_reg, a_ab_qm)
! 
! Assume a plasma composed of several species b, each separately in 
! thermal equilibrium with themselves but not necessarily with each 
! other[1]. This routine returns several useful components of the 
! corresponding C-coefficients introduced in Note [2] below (BPS).
! 
! UNITS: A_{pb} has units of [MeV/micron] (subject to change in updates)
! 
! THE PHYSICS:
! The various subsystems b will exchange coulomb energy and they will
! eventually equilibrate to a common temperature. The C-coefficients 
! introduced in Ref. [2] encode this coulomb energy exchange, exactly to 
! leading and next-to-leading orders in the plasma coupling constant g.
! See Refs. [3,4,5] for more details. For a weakly coupled plasma (g << 1), 
! the BPS calculation is essentially exact, and the error is O(g). Physical 
! properties of interest, such as the stopping power dE/dx and the temperature 
! equilibration rate between plasma species, can be obtained directly from 
! the C-coefficients. 
!
! USAGE:
! Since electrons are thousands of times lighter than ions, separate
! electron and ion temperatures may develop, Te and Ti. The output is
! organized into electron contributions and total ion contributions
!(sum over all ions). The plasma must be weakly coupled.
!
! INPUT: nni, ep, zp, mp, ia, ib, betab, zb, mb, nb
!
! Describe the incident projectile and the background plasma. 
!
! projectile input quantities:
! ep : classical kinetic energy of the projectile [keV]
! zp : charge of the projectile in units of Z_p [dimensionless]
! mp : mass of the projectile [keV], i.e. mp = mp[grams]*c^2
!
! plasma input quantities:
! nni  : Number of total plasma species = number ion species + 1
! zb   : Charges of the plasma species. By convention zp(1) is the 
!      : electron plasma component. [dimensionless, Array]
! betab: Inverse temperatures of the plasma components. For an
!        electron-ion plasma, set betab(1)=1/T_e and all other
!        values of the array to 1/T_I [keV^-1].
! mb   : Masses of the plasma species [keV]. 
! nb   : Number densities of the plasma species [cm^-3]. 
! ia   : First plasma species [usually the projectile]
! ib   : Second plasma species
!
! OUTPUT: a_ab, a_ab_sing, a_ab_reg, a_ab_qm
!
! Each plasma component b makes a linear contribution A_b to the total 
! C-coefficient, i.e. A = sum_b A_b [5]. Each A_b in turn can be be 
! decomposed into a classical-quantum or electron-ion contributions.
!
! classical electron  : ac_e
! classical ion       : ac_i [sum over all ions]
! classical total     : ac_tot = ac_e + ac_i 
! quantum   electron  : aq_e
! quantum   ion       : aq_i [sum over all ions]
! quantum   total     : a_tot = aq_e + aq_i
! total     electric  : a_e = ac_e + aq_e
! total     ion       : a_i = ac_i + aq_i
! total               ! a_tot = a_e + a_i
!
! NOTES:
! [1] The temperatures T_b may therefore all differ. By convention
!     I take b=1 for the electron component of the plasma. A very
!     useful and interesting parameter regime is the one in which 
!     the ions have a common temperature T_I and the electron have
!     a temperature T_e, usually with T_e =/= T_I.  See also USAGE
!     and note [3] below.
!
! [2] BPS paper
!     L. Brown, D. Preston, and R. Singleton~Jr., 
!     "Charged Particle Motion in a Highly Ionized Plasma",
!     Physics Reports, 410 (2005) 237
!     [arXiv:physics/0501084]
!
! [3] The code employs rationalized cgs units in which the dimensionless
!     plasma coupling parameter is defined by g = e^2 kappa/(4Pi*T); in 
!     these units the Debye wavenumber is determined by kappa^2 = e^2 n/T 
!     and the plasma frequency by omega^2 = e^2 n/m. A weakly coupled 
!     plasma is one for which g << 1, i.e. a plasma with thermal kinetic 
!     energy (of order the temperature T) dominates the coulomb potential 
!     energy (for particles separated by a Debye length). In the more 
!     common non-rationalized cgs units, we define g = e^2 kappa/T, with 
!     kappa^2 = 4 Pi e^2 n/T and omega^2 = 4 Pi e^2 n/m. 
!
! [4] For coulomb energy exchange processes, the leading and next-to-
!     leading order terms in the plasma coupling g are proportional to 
!     g^2*ln(g) and g^2, respectively. That is to say, for a property 
!     denoted by F, one can expand F in powers of g in the form:
!
!       F(g,eta) = A(eta)*g^2*ln(g) + B(eta)*g^2 + O(g^3) 
!                 = A(eta)*g^2*[ ln(C(eta)*g)+O(g) ],
!
!     where eta is the dimensionless quantum parameter (eta <<1 means
!     extreme classical scattering while eta >>1 means the extreme quantum 
!     limit). The relative error of BPS is therefore O(g). At the center of 
!     the sun g=0.4, and so the error of Ref. [1] is only of order 4% in 
!     this case. For the processes of charged stopping power and electron-
!     ion temperature equilibration, Ref. [1] calculates the corresponding
!     functions A(eta) and B(eta) exactly, including all orders in the 
!     two-body quantum scattering parameter eta = e^2/4Pi*hbar*v_thermal 
!     (this means that BPS gives the correct interpolation between the 
!     classical and quantum regimes, exact to leading and next-to-leading 
!     order). The O(g^3) terms physically correspond to 3-body correlations 
!     within the plasma, and for a sufficiently weak plasma these are 
!     negligible. For strongly coupled plasmas (g >> 1), all terms in a 
!     g-expansion are important and the BPS calculation is not applicable.
!
! [5] It makes sense to talk about separate linear contribution A_b
!     contributing *from* a given plasma component b only for a weakly 
!     coupled plasma. More exactly, A=sum_b A_b holds true only up to 
!     leading and next-to-leading order in the plasma coupling g. This is 
!     the order to which Ref. [2] calculates all quantities, and therefore 
!     BPS works to a consistent order g^2*in g + g^2.
!
!
! vp = projectile speed [cm/s]
!
!            [ 2 Ep ]
!    = c sqrt[ ---- ]  where mp and Ep are in keV and
!            [ mpc2 ]  c is the speed of light in cm/s
!
!
!        ep^2 kb^2
! c1  =  --------- = 2 zp**2 * Be * kb^2 * a0  [keV/cm]
!           4 Pi                               [MeV/micron]
!
!       [ betab mbc2 ]^1/2  
! c2 =  |----------- |      vp   [dimensionless]
!       [   2 Pi     ]      
!        
! where
!          e^2 
! Be  =  --------- = 13.606E-3   [keV]
!        8 Pi a0                 Bohr radius: a0=5.29E-9 cm
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! main driver for C-coefficient for general quantum and electron-mass regimes
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      SUBROUTINE bps_ccoeff_ab_mass(nni, ep, mp, zp, ia, ib, betab, zb, mb, nb, &
            a_ab, a_ab_sing, a_ab_reg, a_ab_qm)
      USE physvars
      USE mathvars    
        IMPLICIT NONE                                             ! Plasma:
        INTEGER,                            INTENT(IN)  :: nni    !  number of ions
        REAL,                               INTENT(IN)  :: ep     !  energy input [keV]
        REAL,                               INTENT(IN)  :: mp     !  mass [keV]
        REAL,                               INTENT(IN)  :: zp     !  charge
        INTEGER,                            INTENT(IN)  :: ia     !  
        INTEGER,                            INTENT(IN)  :: ib     !  
        REAL,    DIMENSION(1:nni+1),        INTENT(IN)  :: betab  !  temp array [1/keV]
        REAL,    DIMENSION(1:nni+1),        INTENT(IN)  :: mb     !  mass array [keV]
        REAL,    DIMENSION(1:nni+1),        INTENT(IN)  :: nb     !  density [1/cc]
        REAL,    DIMENSION(1:nni+1),        INTENT(IN)  :: zb     !  charge array
                                                                  !
                                                                  ! C-coeffs [MeV/micron]
        REAL,                               INTENT(OUT) :: a_ab
        REAL,                               INTENT(OUT) :: a_ab_sing
        REAL,                               INTENT(OUT) :: a_ab_reg
        REAL,                               INTENT(OUT) :: a_ab_qm

        REAL,    DIMENSION(1:nni+1)  :: mpb, mbpb, kb2, ab
        REAL                         :: vp, zp2, k, k2, kd, kd2, a, b, eta
        REAL                         :: ac_r, ac_s, aq, c1, c2

        REAL, PARAMETER              :: EPS_SMALL_E=2.E-4
        REAL, PARAMETER              :: EPS_SMALL_E_SING=2.E-4
        REAL, PARAMETER              :: EPS_SMALL_E_REG=2.E-4
!
! initialize components of C-coefficients
!
        kb2=8*PI*A0CM*BEKEV*zb*zb*nb*betab
        kd2 = SUM(kb2)                ! [1/cm^2]
        kd  = SQRT(kd2)               ! [1/cm]
        k2  = kb2(1)                  ! [1/cm^2]
        k   = SQRT(k2)                ! [1/cm]   k = k_e
!
! Loop over charged plasma species
!
        mpb = mp*mb/(mp+mb)            ! [keV]
        mbpb= mb/mpb                   ! [dimensionless]
        vp =CC*SQRT(2*ep/mp)           ! [cm/s]
        zp2=zp**2                      ! [dimensionless]
                                       ! ab=(1/2) betab(ib)*mbc2(ib)*vp2/CC2
        ab  =0.5*betab*mb*vp*vp/CC2    ! [dimensionless] 
        IF (zb(ib) .NE. 0.) THEN
        a  =ab(ib)
        b  =-Log(2*betab(ib)*BEKEV*ABS(zp*zb(ib))*k*A0CM*mbpb(ib) )-2*GAMMA+2
        eta=ABS(zp*zb(ib))*2.1870E8/vp ! defined with projectile velocity vp
        c1=2*zp2*BEKEV*kb2(ib)*A0CM    ! [keV/cm] c1 = e_p^2 kappa_b^2/(4 Pi)
        c1=c1*1.E-7                    ! [MeV/micron]  
        c2=SQRT(a/PI)                  ! [dimensionless] 
                                       ! c2=SQRT(betab(ib)*mb(ib)/TWOPI)*vp/CC 
!
! C_{ab}-classical-singular 
!
        CALL c_sing_mass(a,b,ac_s) 
        a_ab_sing=c1*c2*ac_s
!
! C_{ab}-classical-regular 
!
        CALL c_reg_mass(nni,ia,ib,vp,k2,kb2,betab,mb,ac_r)
        a_ab_reg=c1*ac_r
!
! C_{ab}-quantum
!
        CALL c_quantum_mass(ia,ib,a,eta,aq) ! eta = dimensionless quantum param.
        a_ab_qm=c1*c2*aq
!
! C_{ab}-total
!
        a_ab=a_ab_sing + a_ab_reg + a_ab_qm
        ENDIF
      END SUBROUTINE bps_ccoeff_ab_mass
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Assembles the matrix C_{ab} of the C-coefficients.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      SUBROUTINE bps_ccoeff_ab_matrix(nni, ep, betab, zb, mb, nb,    &
        c_ab, c_ab_sing, c_ab_reg, c_ab_qm, c_tot, c_i, c_e, cc_tot, &
        cc_i, cc_e, cq_tot, cq_i, cq_e, cc_s_i, cc_s_e, cc_r_i, cc_r_e)
      USE physvars
      USE mathvars    
        IMPLICIT NONE                                             ! Plasma:
        INTEGER,                            INTENT(IN)  :: nni    !  number of ions
        REAL,                               INTENT(IN)  :: ep     !  energy
        REAL,    DIMENSION(1:nni+1),        INTENT(IN)  :: betab  !  temp array [1/keV]
        REAL,    DIMENSION(1:nni+1),        INTENT(IN)  :: zb     !  charge array
        REAL,    DIMENSION(1:nni+1),        INTENT(IN)  :: mb     !  mass array [keV]
        REAL,    DIMENSION(1:nni+1),        INTENT(IN)  :: nb     !  density [1/cc]
                                                                  !
                                                                  ! C-coeffs [MeV/micron]
        REAL,    DIMENSION(1:nni+1,1:nni+1),INTENT(OUT) :: c_ab
        REAL,    DIMENSION(1:nni+1,1:nni+1),INTENT(OUT) :: c_ab_sing
        REAL,    DIMENSION(1:nni+1,1:nni+1),INTENT(OUT) :: c_ab_reg
        REAL,    DIMENSION(1:nni+1,1:nni+1),INTENT(OUT) :: c_ab_qm
        REAL,    DIMENSION(1:nni+1),        INTENT(OUT) :: c_tot
        REAL,    DIMENSION(1:nni+1),        INTENT(OUT) :: c_i
        REAL,    DIMENSION(1:nni+1),        INTENT(OUT) :: c_e
        REAL,    DIMENSION(1:nni+1),        INTENT(OUT) :: cc_tot
        REAL,    DIMENSION(1:nni+1),        INTENT(OUT) :: cc_i
        REAL,    DIMENSION(1:nni+1),        INTENT(OUT) :: cc_e
        REAL,    DIMENSION(1:nni+1),        INTENT(OUT) :: cq_tot
        REAL,    DIMENSION(1:nni+1),        INTENT(OUT) :: cq_i
        REAL,    DIMENSION(1:nni+1),        INTENT(OUT) :: cq_e
        REAL,    DIMENSION(1:nni+1),        INTENT(OUT) :: cc_s_i
        REAL,    DIMENSION(1:nni+1),        INTENT(OUT) :: cc_s_e
        REAL,    DIMENSION(1:nni+1),        INTENT(OUT) :: cc_r_i
        REAL,    DIMENSION(1:nni+1),        INTENT(OUT) :: cc_r_e

        REAL    :: cab, cab_sing, cab_reg, cab_qm
        REAL    :: mp, zp
        INTEGER :: ia, ib

        c_i   = 0
        cc_s_i= 0
        cc_r_i= 0
        cc_i  = 0
        cq_i  = 0
        DO ia=1,nni+1
          mp=mb(ia)
          zp=zb(ia)
          DO ib=1,nni+1
            CALL bps_ccoeff_ab_mass(nni, ep, mp, zp, ia, ib, betab, zb, mb, nb, &
            cab, cab_sing, cab_reg, cab_qm) !*! change to bps_acoeff_ab_mass
            c_ab(ia,ib)     =cab
            c_ab_sing(ia,ib)=cab_sing
            c_ab_reg(ia,ib) =cab_reg
            c_ab_qm(ia,ib)  =cab_qm 
            IF (ib == 1) THEN
               c_e(ia)   = cab 
               cc_s_e(ia)= cab_sing
               cc_r_e(ia)= cab_reg
               cc_e(ia)  = cab_sing + cab_reg
               cq_e(ia)  = cab_qm
            ELSE
               c_i(ia)   = c_i(ia)    + cab 
               cc_s_i(ia)= cc_s_i(ia) + cab_sing
               cc_r_i(ia)= cc_r_i(ia) + cab_reg
               cc_i(ia)  = cc_i(ia)   + cab_sing + cab_reg
               cq_i(ia)  = cq_i(ia)   + cab_qm
            ENDIF
          ENDDO
          c_tot(ia) = c_e(ia)  + c_i(ia)
          cc_tot(ia)= cc_e(ia) + cc_i(ia)
          cq_tot(ia)= cq_e(ia) + cq_i(ia)
        ENDDO
      END SUBROUTINE bps_ccoeff_ab_matrix


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Returns C_{p I} = \sum_i C_{p i} for backward compatibility
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      SUBROUTINE bps_ccoeff_ei_mass(nni, ep, zp, mp, betab, zb, mb, nb, &
            c_tot, c_i, c_e, cc_tot, cc_i, cc_e, cq_tot, cq_i, cq_e, &
            cc_s_i, cc_s_e, cc_r_i, cc_r_e)
      USE physvars
      USE mathvars    
      USE controlvars  
        IMPLICIT NONE                                      ! Plasma:
        INTEGER,                     INTENT(IN)  :: nni    !  number of ions
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: betab  !  temp array [1/keV]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: mb     !  mass array [keV]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: nb     !  density [1/cc]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: zb     !  charge array
                                                           !
                                                           ! Projectile  
        REAL,                        INTENT(IN)  :: ep     !  projectile energy [keV]
        REAL,                        INTENT(IN)  :: mp     !  projectile mass   [keV]
        REAL,                        INTENT(IN)  :: zp     !  projectile charge
                                                           !
                                                           ! C-coeffs [MeV/micron]
        REAL,                        INTENT(OUT) :: c_tot  !  electron + ion
        REAL,                        INTENT(OUT) :: c_i    !  ion contribution
        REAL,                        INTENT(OUT) :: c_e    !  electron contribution
        REAL,                        INTENT(OUT) :: cc_tot !  classical
        REAL,                        INTENT(OUT) :: cc_i   !  classical
        REAL,                        INTENT(OUT) :: cc_e   !  classical
        REAL,                        INTENT(OUT) :: cq_tot !  quantum
        REAL,                        INTENT(OUT) :: cq_i   !  quantum
        REAL,                        INTENT(OUT) :: cq_e   !  quantum
        REAL,                        INTENT(OUT) :: cc_s_i
        REAL,                        INTENT(OUT) :: cc_s_e 
        REAL,                        INTENT(OUT) :: cc_r_i
        REAL,                        INTENT(OUT) :: cc_r_e

        REAL     :: cdum, cc_s, cc_r, cq
        INTEGER  :: ia, ib, nnb
!
! initialize components of C-coefficients
!
        c_tot =0  ! electron + ion
        c_i   =0  ! ion contribution
        c_e   =0  ! electron contribution
        cc_tot=0  ! classical total
        cc_e  =0  ! classical electron
        cc_i  =0  ! classical ion
        cq_tot=0  ! quantum total
        cq_e  =0  ! quantum electron
        cq_i  =0  ! quantum ion
        cc_s_i=0 
        cc_s_e=0 
        cc_r_i=0
        cc_r_e=0

        NNB = nni+1                 ! number of ions + electrons
        ia=1
        DO ib=1,nni+1
        IF (zb(ib) .NE. 0.) THEN
            CALL bps_ccoeff_ab_mass(nni, ep, mp, zp, ia, ib, betab, zb, mb, nb, &
            cdum, cc_s, cc_r, cq)
            CALL x_collect(ib, NNB, cc_s, cc_r, cq,       &
            c_tot, c_i, c_e, cc_tot, cc_i, cc_e, cq_tot,  &
            cq_i, cq_e, cc_s_i, cc_s_e, cc_r_i, cc_r_e)
        ENDIF
        ENDDO
      END SUBROUTINE bps_ccoeff_ei_mass
      
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! singular contribution for non-zero electron mass
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      SUBROUTINE c_sing_mass(a, b, cc_s)
        REAL,    INTENT(IN)  :: a
        REAL,    INTENT(IN)  :: b
        REAL,    INTENT(OUT) :: cc_s
        REAL                 :: u0, u1, du, u, um
        INTEGER, PARAMETER :: NS=1000 ! integration regions: must be even
        REAL,    PARAMETER :: UPM=0.7745966692E0 ! parameters for Gaussian Quad
        REAL,    PARAMETER :: W13=0.5555555556E0, W2=0.8888888889E0
           ac_s=0
           u0=0
           u1=1
           du=(u1-u0)/NS
           u=u0-du
           DO iu=1,NS,2 ! Gaussian quadrature
              u=u+2.E0*du
              cc_s=ac_s+W2*dcab_sing(u,a,b)
              um=u-du*UPM
              cc_s=ac_s+W13*dcab_sing(um,a,b)
              um=u+du*UPM
              cc_s=ac_s+W13*dcab_sing(um,a,b)
           ENDDO
           cc_s=cc_s*du
      END SUBROUTINE c_sing_mass
!
      FUNCTION dcab_sing(u, a, b)
        IMPLICIT NONE
        REAL,        INTENT(IN)  :: u        ! [dimensionless]
        REAL,        INTENT(IN)  :: a        ! [dimensionless] 
                                             ! a=(1/2)*beta*mpc2*vp^2/C^2
        REAL,        INTENT(IN)  :: b        ! [dimensionless]
        REAL                     :: dcab_sing! [dimensionless]
        dcab_sing=SQRT(u)*EXP(-a*u)*(-LOG(u/(1-u)) + b)
      END FUNCTION dcab_sing

!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! regular contribution for non-zero electron mass
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      SUBROUTINE c_reg_mass(nni, ia, ib, vp, k2, kb2, betab, mb, cc_r)
      USE physvars
        IMPLICIT NONE
        INTEGER,                     INTENT(IN)  :: nni 
        INTEGER,                     INTENT(IN)  :: ia
        INTEGER,                     INTENT(IN)  :: ib
        REAL,                        INTENT(IN)  :: vp
        REAL,                        INTENT(IN)  :: k2
        REAL,                        INTENT(IN)  :: kb2
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: betab
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: mb
        REAL,                        INTENT(OUT) :: cc_r
        REAL,    DIMENSION(1:nni+1)  :: ab
!       INTEGER, PARAMETER :: NR=10 ! integration regions: must be even
        INTEGER, PARAMETER :: NR=100 ! integration regions: must be even
        REAL,    PARAMETER :: UPM=0.7745966692E0 ! parameters for Gaussian Quad
        REAL,    PARAMETER :: W13=0.5555555556E0, W2=0.8888888889E0
        REAL               :: u0, u1, du, u, um, d_cab_reg
        INTEGER            :: iu
        ab=SQRT(0.5*betab*mb)*vp/CC
        cc_r=0
        u0=0.0
        u1=1.
!       u1=MIN(1.,5/(ab(ib)**2)) ! support can lie << 1
        du=(u1-u0)/NR
        u=u0-du
        DO iu=1,NR,2 ! Gaussian quadrature
           u=u+2.*du
           cc_r=cc_r+W2*d_cab_reg(u,vp,ia,ib,nni,k2,kb2,betab,mb)
           um=u-du*UPM
           cc_r=cc_r+W13*d_cab_reg(um,vp,ia,ib,nni,k2,kb2,betab,mb)
           um=u+du*UPM
           cc_r=cc_r+W13*d_cab_reg(um,vp,ia,ib,nni,k2,kb2,betab,mb)
        ENDDO
        cc_r=cc_r*du
      END SUBROUTINE c_reg_mass
      
      FUNCTION d_cab_reg(u, vp, ia, ib, nni, k2, kb2, betab, mb)
      USE mathvars
      USE physvars
        IMPLICIT NONE
        REAL,                        INTENT(IN)  :: u      ! [dimensionless]
        REAL,                        INTENT(IN)  :: vp     ! Projectile velocity [cm/s]
        INTEGER,                     INTENT(IN)  :: ia     ! Species number
        INTEGER,                     INTENT(IN)  :: ib     ! Species number
        INTEGER,                     INTENT(IN)  :: nni    ! Number of ion species
        REAL,                        INTENT(IN)  :: k2     ! Wave-number squared [1/cm^2]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: kb2    ! Debye wavenumber squared [1/cm^2]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: betab  ! Temperature array [1/keV]
        REAL,    DIMENSION(1:nni+1), INTENT(IN)  :: mb     ! Mass array [keV]
        REAL                                     :: d_cab_reg! [dimensionless]
        REAL,    DIMENSION(1:nni+1) :: kbar2b, ab, ab2
        REAL                        :: fr, fi, fabs, farg, h, r_ib
        REAL                        :: kcb, bm_ic, bm_ib, a_ic, a_ib, ex, au
        INTEGER                     :: ic
        ab=SQRT(0.5*betab*mb)*vp/CC
        ab2=ab*ab
        kbar2b=kb2/k2
        CALL frfi(u,nni,kbar2b,ab,fr,fi,fabs,farg)
        h=2*(fr*farg + fi*LOG(fabs))*u
        r_ib=0
        bm_ib=betab(ib)*mb(ib)
        a_ib =ab(ib)*ab(ib)
        DO ic=1,nni+1
           kcb=kb2(ic)/k2
           bm_ic=betab(ic)*mb(ic)
           a_ic =ab(ic)*ab(ic)
           IF (ic == ib) THEN
              ex=1.
           ELSE
              au=(a_ic-a_ib)/u !*! is this correct for C^ll?
              ex=EXP(-au)
           ENDIF
           r_ib=r_ib + kcb*SQRT(bm_ic/bm_ib)*ex
        ENDDO      
        r_ib=1./r_ib
!
        d_cab_reg=-r_ib*h/TWOPI
      END FUNCTION d_cab_reg

!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! quantum contribution for non-zero electron mass
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      SUBROUTINE c_quantum_mass(ia, ib, a, eta, aq)
        IMPLICIT NONE
        INTEGER, INTENT(IN)  :: ia    ! species index
        INTEGER, INTENT(IN)  :: ib    ! species index
        REAL,    INTENT(IN)  :: a     ! [dimensionless] (1/2) betab mb vp^2
        REAL,    INTENT(IN)  :: eta   ! [dimensionless] ep eb/4pi hbar vp
        REAL,    INTENT(OUT) :: aq 
        REAL               :: u0, u1, du, u, um
        INTEGER, PARAMETER :: NQ=1000            ! integration regions quantum : must be even
        REAL,    PARAMETER :: UPM=0.7745966692E0 ! parameters for Gaussian Quad
        REAL,    PARAMETER :: W13=0.5555555556E0, W2=0.8888888889E0
        REAL    :: dcq
        INTEGER :: iu
        aq=0
        u0=0.
        aq=0
        IF (ib == ia) THEN
           u0=0
           u1=4./SQRT(a)
        ELSE
           u0=1-10./SQRT(a)
           u0=MAX(0.,u0)  
           u1=1+10./SQRT(a)
        ENDIF
        du=(u1-u0)/NQ
        u=u0-du
        DO iu=1,NQ,2 ! Gaussian quadrature
           u=u+2.E0*du
           aq=aq+W2*dcq(u,a,eta)
           um=u-du*UPM
           aq=aq+W13*dcq(um,a,eta)
           um=u+du*UPM
           aq=aq+W13*dcq(um,a,eta)
        ENDDO
        aq=aq*du
      END SUBROUTINE c_quantum_mass

      FUNCTION dcq(u, a, eta)
      USE physvars
        IMPLICIT NONE
        REAL,                        INTENT(IN)  :: u          ! [dimensionless]
        REAL,                        INTENT(IN)  :: a          ! [dimensionless]
        REAL,                        INTENT(IN)  :: eta        ! [dimensionless]
        REAL                                     :: dcq  ! [dimensionless]
        REAL            :: repsi, au, eu, au2, ap, am, psilog, ch, sh, csh
        eu=eta/u 
        psilog=repsi(eu) - LOG(eu)
        au =2*a*u
        au2=a*u*u
        ap = au-au2-a
        am =-au-au2-a
        ch =0.5*(EXP(ap)+EXP(am))
        sh =0.5*(EXP(ap)-EXP(am))
        csh=2*(ch - sh/au)/au
        dcq=-psilog*csh
      END FUNCTION dcq


