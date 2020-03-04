      PROGRAM dedx
!
! This program finds the range of a charged partcile in a plasma, where the
! charged particle is created with threshold energy Ep. 
!
! nts = number of time steps
!
      USE allocatablevars
      USE bpsvars
      USE physvars
      USE mathvars
        IMPLICIT NONE
        REAL    :: dedx_tot,  dedx_i,  dedx_e
        REAL    :: dedxc_tot, dedxc_i, dedxc_e
        REAL    :: dedxq_tot, dedxq_i, dedxq_e
        REAL    :: dedx_a_tot, dedx_a_i, dedx_a_e, dedxc_a_tot, dedxc_a_i, dedxc_a_e
        REAL    :: dedxq_a_tot, dedxq_a_i, dedxq_a_e, dedxc_a_s_i, dedxc_a_s_e
        REAL    :: dedxc_a_r_i, dedxc_a_r_e

        INTEGER :: j, nit

        REAL    :: te, ti, ne, ep, mp, zp, epp
        INTEGER :: nni

        ! * 3/3/2020
        ! Compare dE/dx as calculated by this code and using a numerical
        ! derivative of A. dedx.f90 vs acoeff.f90
        ! Ep=3.5MeV has good agreement.
        ! Ep=0.9MeV has bad agreement
        ! The lower velocity behavior is wrong
        !
      
!
! number of iterations
        nit=100
!
! A: alpha particle projectile
!
        ep=3540.    ! Projectile energy  [keV]
        mp=4*MPKEV  ! Projectile mass    [keV]
        zp=2.       ! Projectile charge  [e]
!
!
! DT plasma with alpha particle projectile
!
        CALL define_plasma_dt(te,ti,ne,nni)
!
! plot the regular and singular contributions
!
        OPEN  (1, FILE='plot_acoeff.out')
        CALL write_output(ep,mp,zp,te,ti,ne,nni,betab,zb,mb,nb)
!
! evolution
!
        WRITE(6,'(A)') '#'
        WRITE(6,'(A, 10X,A7, 10X,A6, 15X,A6, 15X,A8)') '#','E [MeV]','dedx_e', 'dedx_I', 'dedx_tot'
        WRITE(6,'(A)') '#'
        WRITE(1,'(A)') '#'
        WRITE(1,'(A, 10X,A7, 10X,A6, 15X,A6, 15X,A8)') '#','E [MeV]','dedx_e', 'dedx_I', 'dedx_tot'
        WRITE(1,'(A)') '#'
        epp=ep ! [keV] projectile energy         
        j = 26
        IF (epp .EQ. 0) epp=1.E-5
        CALL dedx_bps(nni, epp, zp, mp, betab, zb, mb, nb,   &
             dedx_tot, dedx_i, dedx_e, dedxc_tot, dedxc_i, & 
             dedxc_e, dedxq_tot, dedxq_i, dedxq_e) ! [MeV/micron] with epp/1000. in MeV
        
        CALL acoeff_dedx_bps(nni,ep,zp,mp,betab,zb,mb,nb,  &
            dedx_a_tot, dedx_a_i, dedx_a_e, dedxc_a_tot, dedxc_a_i, dedxc_a_e, & 
            dedxq_a_tot, dedxq_a_i, dedxq_a_e, dedxc_a_s_i, dedxc_a_s_e,       &
            dedxc_a_r_i, dedxc_a_r_e)

! total
        WRITE (6,'(I6,E17.8,6E22.13)') j, epp/1000., dedx_e, dedx_i, dedx_tot
        WRITE (6,'(I6,E17.8,6E22.13)') j, epp/1000., dedx_a_e, dedx_a_i, dedx_a_tot        
        WRITE (1,'(I6,E17.8,6E22.13)') j, epp/1000., dedx_e, dedx_i, dedx_tot
        WRITE (1,'(I6,E17.8,6E22.13)') j, epp/1000., dedx_a_e, dedx_a_i, dedx_a_tot

! classical
        WRITE (6,'(I6,E17.8,6E22.13)') j, epp/1000., dedxc_e, dedxc_i, dedxc_tot
        WRITE (6,'(I6,E17.8,6E22.13)') j, epp/1000., dedxc_a_e, dedxc_a_i, dedxc_a_tot        
        WRITE (1,'(I6,E17.8,6E22.13)') j, epp/1000., dedxc_e, dedxc_i, dedxc_tot
        WRITE (1,'(I6,E17.8,6E22.13)') j, epp/1000., dedxc_a_e, dedxc_a_i, dedxc_a_tot

! quantum
        WRITE (6,'(I6,E17.8,6E22.13)') j, epp/1000., dedxq_e, dedxq_i, dedxq_tot
        WRITE (6,'(I6,E17.8,6E22.13)') j, epp/1000., dedxq_a_e, dedxq_a_i, dedxq_a_tot        
        WRITE (1,'(I6,E17.8,6E22.13)') j, epp/1000., dedxq_e, dedxq_i, dedxq_tot
        WRITE (1,'(I6,E17.8,6E22.13)') j, epp/1000., dedxq_a_e, dedxq_a_i, dedxq_a_tot
        
        
        CLOSE (1)
        END PROGRAM dedx

    SUBROUTINE define_plasma_dt(te, ti, ne, nni)
    USE allocatablevars
    USE physvars
      IMPLICIT NONE
      REAL                                :: te     ! [keV]
      REAL                                :: ti     ! [keV]
      REAL                                :: ne     ! [cm^-3]
      INTEGER                             :: nni    ! number of ion species
!     REAL,    DIMENSION(:), ALLOCATABLE  :: betab  ! [1/keV]
!     REAL,    DIMENSION(:), ALLOCATABLE  :: mb     ! [keV]
!     REAL,    DIMENSION(:), ALLOCATABLE  :: nb     ! [cm^-3]
!     REAL,    DIMENSION(:), ALLOCATABLE  :: zb     ! [e]
!     REAL,    DIMENSION(:), ALLOCATABLE  :: gb, etab, mpb, mbpb

      te=1.       ! Electron temperature       [keV]
      ti=1.       ! Ion temperature            [keV]
      ne=2.E26    ! electron number density    [1/cc]    
      nni=2       ! number of ion species
      
      ALLOCATE(betab(1:nni+1),zb(1:nni+1),mb(nni+1),nb(1:nni+1))   ! allocatablevars
      ALLOCATE(gb(1:nni+1),etab(1:nni+1),mpb(nni+1),mbpb(1:nni+1)) ! allocatablevars

      zb(1)=-1.    ! Species charges
      zb(2)=+1.    ! 
      zb(3)=+1.    !
      mb(1)=MEKEV  ! Species masses [keV]
      mb(2)=2*MPKEV!
      mb(3)=3*MPKEV!
!
! Construct density and temperature arrays
!
      nb(1)=1.                          ! ONLY FOR EQUIMOLAR DT
      nb(2:nni+1)=1./(zb(2:nni+1)*nni)  ! charge neutrality
      nb=nb*ne                          ! number density array [cm^-3]
      betab(1)=1./te                    ! inverse temp array   [keV^-1]
      betab(2:nni+1)=1./ti              !
    END SUBROUTINE define_plasma_dt


    SUBROUTINE write_output(ep, mp, zp, te, ti, ne, nni, betab, zb, mb, nb)
    USE physvars
      IMPLICIT NONE
      REAL,                        INTENT(IN) :: ep     ! [keV]
      REAL,                        INTENT(IN) :: mp     ! [keV]
      REAL,                        INTENT(IN) :: zp     ! [e]
      REAL,                        INTENT(IN) :: te     ! [keV]
      REAL,                        INTENT(IN) :: ti     ! [keV]
      REAL,                        INTENT(IN) :: ne     ! [cm^-3]
      INTEGER,                     INTENT(IN) :: nni    ! number of ion species
      REAL,    DIMENSION(1:nni+1), INTENT(IN) :: betab  ! [1/keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN) :: zb     ! [e]
      REAL,    DIMENSION(1:nni+1), INTENT(IN) :: mb     ! [keV]
      REAL,    DIMENSION(1:nni+1), INTENT(IN) :: nb     ! [cm^-3]


      REAL  :: etae, ge, gi, ze
      REAL,    DIMENSION(1:nni+1) :: gb 

      REAL  :: vp
      REAL,    DIMENSION(1:nni+1) :: ab, etab
!
! write header
!
      WRITE(6,'(A)') '#'
      WRITE(6,'(A, 3X,A21)') '#','Projectile Parameters'
      WRITE(6,'(A, 4X,A13, X,D12.4)') '#','charge      :', zp
      WRITE(6,'(A, 4X,A13, X,D12.4)') '#','mass   [amu]:', mp/AMUKEV
      WRITE(6,'(A, 4X,A13, X,D12.4)') '#','mass   [keV]:', mp
      WRITE(6,'(A, 4X,A13, X,D12.4)') '#','energy [keV]:', ep
      WRITE(6,'(A)') '#'
      WRITE(6,'(A, 3X,A17)') '#','Plasma Parameters'
      WRITE(6,'(A, 4X,A28, X,D12.4, X,D12.4)')  '#','electron and ion temp [keV]:', te, ti
      WRITE(6,'(A, 4X,A28, X,D12.4, X,5D12.4)') '#','number density nb   [cm^-3]:', nb
      WRITE(6,'(A, 4X,A28, X,D12.4, X,5D12.4)') '#','mass array mb         [amu]:', mb/AMUKEV
      WRITE(6,'(A, 4X,A28, X,D12.4, X,5D12.4)') '#','mass array mb         [keV]:', mb
      WRITE(6,'(A, 4X,A28, X,D12.4, X,5D12.4)') '#','charge array zb            :', zb

      WRITE(1,'(A)') '#'  
      WRITE(1,'(A)') '#'
      WRITE(1,'(A, 3X,A17)') '#','Plasma Parameters'
      WRITE(1,'(A, 4X,A28, X,D12.4, X,D12.4)')  '#','electron and ion temp [keV]:', te, ti
      WRITE(1,'(A, 4X,A28, X,D12.4, X,5D12.4)') '#','number density nb   [cm^-3]:', nb
      WRITE(1,'(A, 4X,A28, X,D12.4, X,5D12.4)') '#','mass array mb         [amu]:', mb/AMUKEV
      WRITE(1,'(A, 4X,A28, X,D12.4, X,5D12.4)') '#','mass array mb         [keV]:', mb
      WRITE(1,'(A, 4X,A28, X,D12.4, X,5D12.4)') '#','charge array zb            :', zb
      WRITE(1,'(A)') '#'
      WRITE(1,'(A)') '#'
      WRITE(1,'(A, 3X,A21)') '#','Projectile Parameters'
      WRITE(1,'(A, 4X,A13, X,D12.4)') '#','charge      :', zp
      WRITE(1,'(A, 4X,A13, X,D12.4)') '#','mass   [amu]:', mp/AMUKEV
      WRITE(1,'(A, 4X,A13, X,D12.4)') '#','mass   [keV]:', mp
      WRITE(1,'(A, 4X,A13, X,D12.4)') '#','energy [keV]:', ep
!
! Print plasma parameters
!

      CALL param(nni, betab, nb, gb, ge, gi, etae, ze)
      WRITE(6,'(A1, 2X,A9, 4X,A7, 4X,A7, 8X,A2, 7X,A2, 10X,A4, 9X,A9)') &
        '#','ne[cm^-3]','Te[keV]', 'Ti[keV]','ge','gi','etae','ze/2**1.5'
      WRITE(6,'(A1,11D12.4)') '#',ne, te, ti, ge, gi, etae, ze/2**1.5
      WRITE(1,'(A1, 2X,A9, 4X,A7, 4X,A7, 8X,A2, 7X,A2, 10X,A4, 9X,A9)') &
        '#','ne[cm^-3]','Te[keV]', 'Ti[keV]','ge','gi','etae','ze/2**1.5'
      WRITE(1,'(A1,11D12.4)') '#',ne, te, ti, ge, gi, etae, ze/2**1.5

      vp  = CC*SQRT(2*ep/mp)       ! [cm/s]
      ab  =0.5*betab*mb*vp*vp/CC2  ! [dimensionless] (1/2) betab(ib)*mb(ib)*vp2/CC2
      etab=ABS(zp*zb)*2.1870E8/vp  ! [dimensionless]

      WRITE(6,'(A7, 3D12.4)') '# ab  =', ab
      WRITE(6,'(A7, 3D12.4)') '# etab=', etab
      WRITE(1,'(A7, 3D12.4)') '# ab  =', ab
      WRITE(1,'(A7, 3D12.4)') '# etab=', etab

!      WRITE(6,'(A)') "#"
!      WRITE(6,'(A5, 2X,A7, 5X,A5, 7X,A3, 9X,A3, 9X,A6, 6X,A4, 8X,A4, 8X,A6, 6X,A4, 8X,A4)') &
!        "#  iu", "ep[keV]", "a_tot", "a_e", "a_i", "ac_tot", "ac_e_lim", "ac_i", &
!        "aq_tot", "aq_e", "aq_i"
!      WRITE(1,'(A)') "#"
!      WRITE(1,'(A5, 2X,A7, 5X,A5, 7X,A3, 9X,A3, 9X,A6, 6X,A4, 8X,A4, 8X,A6, 6X,A4, 8X,A4)') &
!        "#  iu", "ep[keV]", "a_tot", "a_e", "a_i", "ac_tot", "ac_e", "ac_i", "aq_tot", "aq_e","aq_i"

    END SUBROUTINE write_output

