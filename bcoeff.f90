! 
! ROUTINE: bps_bcoeff_ab_mass(nni, ep, mp, zp, ia, ib, betab, zb, mb, nb, &
!            b_ab, b_ab_sing, b_ab_reg, b_ab_qm)
! 
! Assume a plasma composed of several species b, each separately in 
! thermal equilibrium with themselves but not necessarily with each 
! other[1]. This routine returns several useful components of the 
! corresponding C-coefficients introduced in Note [2] below (BPS).
! 
! UNITS: B_{pb} has units of **[MeV/micron] (subject to change in updates)
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
! classical electron  : bc_e
! classical ion       : bc_i [sum over all ions]
! classical total     : bc_tot = bc_e + bc_i 
! quantum   electron  : bq_e
! quantum   ion       : bq_i [sum over all ions]
! quantum   total     : b_tot = bq_e + bq_i
! total     electric  : b_e = bc_e + bq_e
! total     ion       : b_i = bc_i + bq_i
! total               ! b_tot = b_e + b_i
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
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! main driver for B-coefficient for general quantum and electron-mass regimes
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      SUBROUTINE bps_bcoeff_ab_mass(nni, ep, mp, zp, ia, ib, betab, zb, mb, nb, &
            b_ab, b_ab_sing, b_ab_reg, b_ab_qm)
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
        REAL,                               INTENT(OUT) :: b_ab
        REAL,                               INTENT(OUT) :: b_ab_sing
        REAL,                               INTENT(OUT) :: b_ab_reg
        REAL,                               INTENT(OUT) :: b_ab_qm

        REAL :: c_ab, a_ab, vp
        REAL :: c_ab_sing, a_ab_sing
        REAL :: c_ab_reg, a_ab_reg
        REAL :: c_ab_qm, a_ab_qm

        vp = CC*SQRT(2*ep/mp)                    
        
        CALL bps_ccoeff_ab_mass(nni, ep, mp, zp, ia, ib, betab, zb, mb, nb, &
             c_ab, c_ab_sing, c_ab_reg, c_ab_qm)
        CALL bps_acoeff_ab_mass(nni, ep, mp, zp, ia, ib, betab, zb, mb, nb, &
             a_ab, a_ab_sing, a_ab_reg, a_ab_qm)
        b_ab = c_ab - CC*a_ab/betab(ib)/vp/1000.
        b_ab_sing = c_ab_sing - CC*c_ab_sing/betab(ib)/vp/1000.
        b_ab_reg = c_ab_reg - CC*c_ab_reg/betab(ib)/vp/1000.
        b_ab_qm = c_ab_qm - CC*c_ab_qm/betab(ib)/vp/1000.
        !*! need to form B from C and A.
        
        
        
      END SUBROUTINE bps_bcoeff_ab_mass
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Assembles the matrix B_{ab} of the B-coefficients.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      SUBROUTINE bps_bcoeff_ab_matrix(nni, ep, betab, zb, mb, nb,    &
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
            CALL bps_bcoeff_ab_mass(nni, ep, mp, zp, ia, ib, betab, zb, mb, nb, &
            cab, cab_sing, cab_reg, cab_qm) 
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
      END SUBROUTINE bps_bcoeff_ab_matrix

!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Returns B_{p I} = \sum_i B_{p i} for backward compatibility
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      SUBROUTINE bps_bcoeff_ei_mass(nni, ep, zp, mp, betab, zb, mb, nb, &
            b_tot, b_i, b_e, bc_tot, bc_i, bc_e, bq_tot, bq_i, bq_e, &
            bc_s_i, bc_s_e, bc_r_i, bc_r_e)
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
                                                           ! B-coeffs [MeV/micron]**
        REAL,                        INTENT(OUT) :: b_tot  !  electron + ion
        REAL,                        INTENT(OUT) :: b_i    !  ion contribution
        REAL,                        INTENT(OUT) :: b_e    !  electron contribution
        REAL,                        INTENT(OUT) :: bc_tot !  classical
        REAL,                        INTENT(OUT) :: bc_i   !  classical
        REAL,                        INTENT(OUT) :: bc_e   !  classical
        REAL,                        INTENT(OUT) :: bq_tot !  quantum
        REAL,                        INTENT(OUT) :: bq_i   !  quantum
        REAL,                        INTENT(OUT) :: bq_e   !  quantum
        REAL,                        INTENT(OUT) :: bc_s_i
        REAL,                        INTENT(OUT) :: bc_s_e 
        REAL,                        INTENT(OUT) :: bc_r_i
        REAL,                        INTENT(OUT) :: bc_r_e

        REAL     :: cdum, cc_s, cc_r, cq, adum, ac_s, ac_r, aq, bc_s, bc_r, bq, vp
        INTEGER  :: ia, ib, nnb
!
! initialize components of C-coefficients
!
        b_tot =0  ! electron + ion
        b_i   =0  ! ion contribution
        b_e   =0  ! electron contribution
        bc_tot=0  ! classical total
        bc_e  =0  ! classical electron
        bc_i  =0  ! classical ion
        bq_tot=0  ! quantum total
        bq_e  =0  ! quantum electron
        bq_i  =0  ! quantum ion
        bc_s_i=0 
        bc_s_e=0 
        bc_r_i=0
        bc_r_e=0

        NNB = nni+1                 ! number of ions + electrons
        ia=1
        DO ib=1,nni+1
        IF (zb(ib) .NE. 0.) THEN
            CALL bps_acoeff_ab_mass(nni, ep, mp, zp, ia, ib, betab, zb, mb, nb, &
                 adum, ac_s, ac_r, aq)
            CALL bps_ccoeff_ab_mass(nni, ep, mp, zp, ia, ib, betab, zb, mb, nb, &
                 cdum, cc_s, cc_r, cq)

            vp = CC*SQRT(2*ep/mp)            
            bc_s = cc_s - CC*ac_s/betab(ib)/vp/100.  !C - A/beta vp
            bc_r = cc_r - CC*ac_r/betab(ib)/vp/100.
            bq = cq - CC*aq/betab(ib)/vp/1000.
            
            CALL x_collect(ib, NNB, bc_s, bc_r, bq,       &
            b_tot, b_i, b_e, bc_tot, bc_i, bc_e, bq_tot,  &
            bq_i, bq_e, bc_s_i, bc_s_e, bc_r_i, bc_r_e)
        ENDIF
        ENDDO
      END SUBROUTINE bps_bcoeff_ei_mass
      
! vp = projectile speed [cm/s]
!
!            [ 2 Ep ]
!    = c sqrt[ ---- ]  where mp and Ep are in keV and
!            [ mpc2 ]  c is the speed of light in cm/s
