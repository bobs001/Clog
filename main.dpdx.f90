      PROGRAM dedx
!
! This program calculates the C-coefficients of BPS
!
! nts = number of time steps
!
      USE allocatablevars
      USE bpsvars
      USE physvars
      USE mathvars
      IMPLICIT NONE

      INTEGER :: j, nit      
!     REAL    ::  p_tot, p_i, p_e, pc_tot, pc_i, pc_e
!     REAL    ::  pq_tot, pq_i, pq_e, pc_s_i, pc_s_e, pc_r_i, pc_r_e

      INTEGER :: nni        
      REAL    :: te, ti, ne, ep, mp, zp, epp, de, scale
      REAL    :: dedx_a_tot, dedx_a_i, dedx_a_e, dedxc_a_tot, dedxc_a_i, dedxc_a_e
      REAL    :: dedxq_a_tot, dedxq_a_i, dedxq_a_e, dedxc_a_s_i, dedxc_a_s_e, dedxc_a_r_i
      REAL    :: dedxc_a_r_e, c_e, c_i, c_tot, cc_e, cc_i, cc_tot, cq_e, cq_i, cq_tot           
      REAL    :: cc_s_i, cc_s_e, cc_r_i, cc_r_e
      REAL    :: dpdx_e, dpdx_i, dpdx_tot, dpdxc_e, dpdxc_i, dpdxc_tot, dpdxq_e, dpdxq_i
      REAL    :: dpdxq_tot, vp, fact
      real    :: dedx_tot, dedx_i, dedx_e, dedxc_tot, dedxc_i, dedxc_e, dedxq_tot, dedxq_i, dedxq_e
      REAL    :: p_tot, pc_tot, pq_tot, p_e, p_i, pc_e, pc_i, pq_e, pq_i
      
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
      OPEN(1, FILE='dpdx_1.out')  ! v \cdot dP/dx
      OPEN(2, FILE='dpdx_2.out')  !
      OPEN(3, FILE='dpdx_3.out')  !
      OPEN(4, FILE='dpdx_4.out')  !      

      CALL write_output(ep,mp,zp,te,ti,ne,nni,betab,zb,mb,nb)
!
! evolution
!
      de=ep/nit
      epp=0
      scale = 1 ! dE/dx keV/cm
      DO j=0,nit
         epp=j*de
         IF (epp .EQ. 0) epp=de/2.0

         ! plot dE/dx and C^ll
         CALL acoeff_dedx_bps(nni, scale, epp, zp, mp, betab, zb, mb, nb,        &
              dedx_a_tot, dedx_a_i, dedx_a_e, dedxc_a_tot, dedxc_a_i, dedxc_a_e, & 
              dedxq_a_tot, dedxq_a_i, dedxq_a_e, dedxc_a_s_i, dedxc_a_s_e,       &
              dedxc_a_r_i, dedxc_a_r_e)
         WRITE(1,'(I6,E17.8,9E22.13)') j, epp, dedx_a_e, dedx_a_i, dedx_a_tot, &
              dedxc_a_e, dedxc_a_i, dedxc_a_tot, dedxq_a_e, dedxq_a_i, dedxq_a_tot
         
         CALL bps_ccoeff_ei_mass(nni, scale, epp, zp, mp, betab, zb, mb, nb, &
              c_tot, c_i, c_e, cc_tot, cc_i, cc_e, cq_tot, cq_i, cq_e,       &
              cc_s_i, cc_s_e, cc_r_i, cc_r_e)
         WRITE(2,'(I6,E17.8,9E22.13)') j, epp, c_e, c_i, c_tot, cc_e, cc_i, cc_tot, &
              cq_e, cq_i, cq_tot
         
         ! plot v*dP/dx using acoeff_dedx_bps
         vp = CC*SQRT(2*epp/mp)
         fact = CC*CC/(mp*vp)
         dpdx_e = dedx_a_e + fact*c_e
         dpdx_i = dedx_a_i + fact*c_i
         dpdx_tot = dedx_a_tot + fact*c_tot
         dpdxc_e = dedxc_a_e + fact*cc_e
         dpdxc_i = dedxc_a_i + fact*cc_i
         dpdxc_tot = dedxc_a_tot + fact*cc_tot
         dpdxq_e = dedxq_a_e + fact*cq_e
         dpdxq_i = dedxq_a_i + fact*cq_i
         dpdxq_tot  = dedxq_a_tot + fact*cq_tot
         WRITE(3,'(I6,E17.8,9E22.13)') j, epp, dpdx_e, dpdx_i, dpdx_tot, dpdxc_e, dpdxc_i, &
              dpdxc_tot, dpdxq_e, dpdxq_i, dpdxq_tot
         !
         !WRITE(6,'(I6,E17.8,9E22.13)') j, epp, mp*(vp/CC)**2*betab(1), dedx_a_e/(fact*c_e), dedx_a_i/(fact*c_i), &
         !     dedx_a_tot/(fact*c_tot)
         
         ! plot v*dP/dx using dedx_bps
         CALL dedx_bps(nni, scale, epp, zp, mp, betab, zb, mb, nb, &
              dedx_tot, dedx_i, dedx_e, dedxc_tot, dedxc_i, & 
              dedxc_e, dedxq_tot, dedxq_i, dedxq_e)
         ! p_tot \equiv v^k dP^k/dx = dE/dx + C^ll/m v
         vp = CC*SQRT(2*ep/mp)            
         p_tot = dedx_tot + c_tot*CC*CC/(mp*vp)
         pc_tot = dedxc_tot + cc_tot*CC*CC/(mp*vp)
         pq_tot = dedxq_tot + cq_tot*CC*CC/(mp*vp)
         p_e = dedx_e + c_e*CC*CC/(mp*vp)
         p_i = dedx_i + c_i*CC*CC/(mp*vp)
         pc_e = dedxc_e + cc_e*CC*CC/(mp*vp)
         pc_i = dedxc_i + cc_i*CC*CC/(mp*vp)
         pq_e = dedxq_e + cq_e*CC*CC/(mp*vp)
         pq_i = dedxq_i + cq_i*CC*CC/(mp*vp)
         WRITE(4,'(I6,E17.8,9E22.13)') j, epp, p_e, p_i, p_tot, pc_e, pc_i, &
              pc_tot, pq_e, pq_i, pq_tot
         
      ENDDO
        
      CLOSE(1)
      CLOSE(2)
      CLOSE(3)
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
      INTEGER :: i
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

      DO i=1,1
         WRITE(i,'(A)') '#'  
         WRITE(i,'(A)') '#'
         WRITE(i,'(A, 3X,A17)') '#','Plasma Parameters'
         WRITE(i,'(A, 4X,A28, X,D12.4, X,D12.4)')  '#','electron and ion temp [keV]:', te, ti
         WRITE(i,'(A, 4X,A28, X,D12.4, X,5D12.4)') '#','number density nb   [cm^-3]:', nb
         WRITE(i,'(A, 4X,A28, X,D12.4, X,5D12.4)') '#','mass array mb         [amu]:', mb/AMUKEV
         WRITE(i,'(A, 4X,A28, X,D12.4, X,5D12.4)') '#','mass array mb         [keV]:', mb
         WRITE(i,'(A, 4X,A28, X,D12.4, X,5D12.4)') '#','charge array zb            :', zb
         WRITE(i,'(A)') '#'
         WRITE(i,'(A)') '#'
         WRITE(i,'(A, 3X,A21)') '#','Projectile Parameters'
         WRITE(i,'(A, 4X,A13, X,D12.4)') '#','charge      :', zp
         WRITE(i,'(A, 4X,A13, X,D12.4)') '#','mass   [amu]:', mp/AMUKEV
         WRITE(i,'(A, 4X,A13, X,D12.4)') '#','mass   [keV]:', mp
         WRITE(i,'(A, 4X,A13, X,D12.4)') '#','energy [keV]:', ep
      ENDDO
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

