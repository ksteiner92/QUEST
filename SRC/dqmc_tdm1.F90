module DQMC_TDM1
#include "dqmc_include.h"

  use DQMC_UTIL
  use DQMC_STRUCT
!  use LAPACK_MOD
!  use BLAS_MOD
  use DQMC_GTAU

  implicit none 
  
  !
  ! This module is designed for the computation of time dependent 
  ! measurement (TDM), which requires the unequal time Green's 
  ! function Gtau (up and down).
  !
  ! Two measurements are taken in this module: Green's function (G) and Chi
  ! (chi) function.  Both are stored in a 3 dimensional array. 
  ! The first dimension is for space, the second is for time, and the 
  ! last dimension is for bins.
  !


  type tdmarray

     integer              :: n
     integer              :: nclass
     integer              :: nk
     integer              :: np

     integer, pointer     :: D(:,:)
     integer, pointer     :: F(:)

     integer, pointer     :: phase(:,:)
     complex(wp), pointer :: ftk(:,:) 
     complex(wp), pointer :: ftw(:,:) 

     real(wp),    pointer :: values(:,:,:) 
     complex(wp), pointer :: valuesk(:,:,:) 

     real(wp),    pointer :: tlink(:,:) 

     character(label_len), pointer :: clabel(:)

  end type tdmarray

  ! Index of the array varaiables
  integer, parameter  :: NTDMARRAY = 10
  integer, parameter  :: IGFUN = 1
  integer, parameter  :: IGFUP = 2
  integer, parameter  :: IGFDN = 3
  integer, parameter  :: ISPXX = 4
  integer, parameter  :: ISPZZ = 5
  integer, parameter  :: IDENS = 6
  integer, parameter  :: IPAIR = 7
  integer, parameter  :: ICOND = 8
  integer, parameter  :: IFSUP = 9
  integer, parameter  :: IFSDN = 10

  ! Index of the array varaiables
  character(len=12), parameter :: &
            pname(NTDMARRAY) = (/ &
                  "Gfun        ", &
                  "Gfun up     ", &
                  "Gfun dn     ", &
                  "SxSx        ", &
                  "SzSz        ", &
                  "Den-Den     ", &
                  "S-wave      ", &
                  "Conductivity", &
                  "GfunSelf up ", &
                  "GfunSelf dn "/)

  type TDM1
     integer  :: L
     integer  :: nbin
     integer  :: avg
     integer  :: err
     integer  :: idx

     integer  :: tmp
     integer  :: cnt

     type(IndexListPtr), dimension(:), allocatable :: classToIdxPair
     logical  :: selfen = .false.

     logical  :: compute=.false.

     real(wp) :: dtau
     real(wp), pointer :: sgn(:) 
     type(tdmarray), pointer :: properties(:) 

     ! Fourier transform matrix for bosonic and fermionic fields
     complex(wp), pointer :: ftwfer(:,:) 
     complex(wp), pointer :: ftwbos(:,:) 

  end type TDM1
  
contains

 !--------------------------------------------------------------------!
  
  subroutine DQMC_TDM1_Init(L, dtau, T1, nBin, S, Gwrap)
    use DQMC_Geom_Wrap
    !
    ! Purpose
    ! =======
    !    This subroutine initializes TDM1. 
    !
    ! Arguments
    ! =========
    !
    type(TDM1), intent(inout) :: T1      ! time dependent measurement
    integer, intent(in)       :: L       ! No of time slice
    integer, intent(in)       :: nBin    ! No of Bins
    real(wp), intent(in)      :: dtau
    type(Struct), intent(in)  :: S
    type(GeomWrap), intent(in):: Gwrap
    type(IndexList), pointer  :: next

    ! ... local variables ...
    integer     :: i, j, k

    ! ... Executable ...

    T1%L      =  L
    T1%dtau   =  dtau
    T1%nbin   =  nBin

    T1%tmp    =  nBin + 1
    T1%avg    =  nBin + 1
    T1%err    =  nBin + 2
    T1%idx    =  1

    T1%compute = .true.
    
    ! Allocate storages
    allocate(T1%sgn(nBin+2))
    ! initialize values
    T1%sgn   = ZERO

    ! generate the array of linked lists which maps for every index pair class onto
    ! all corresponding index pairs of this class
    allocate(T1%classToIdxPair(S%nclass))
    do i = 1, S%nSite
        do j = 1, S%nSite
            k = S%D(i,j)
            if (.not. associated(T1%classToIdxPair(k)%ptr)) then
                allocate(T1%classToIdxPair(k)%ptr)
                T1%classToIdxPair(k)%ptr%i = i
                T1%classToIdxPair(k)%ptr%j = j
                T1%classToIdxPair(k)%ptr%n = 1
                T1%classToIdxPair(k)%ptr%last => T1%classToIdxPair(k)%ptr
            else
                allocate(T1%classToIdxPair(k)%ptr%last%next)
                T1%classToIdxPair(k)%ptr%last%next%i = i
                T1%classToIdxPair(k)%ptr%last%next%j = j
                T1%classToIdxPair(k)%ptr%n = T1%classToIdxPair(k)%ptr%n + 1
                T1%classToIdxPair(k)%ptr%last => T1%classToIdxPair(k)%ptr%last%next
            end if
        end do
    end do

    !do i = 1, S%nclass
    !    write(*,*) "class", i, T1%classToIdxPair(i)%ptr%n, "elements:"
    !    next => T1%classToIdxPair(i)%ptr
    !    do
    !        write(*,*), next%i, next%j
    !        if (.not. associated(next%next)) then
    !            exit
    !        end if
    !        next => next%next
    !    end do
    !end do

    call  DQMC_TDM1_InitFTw(T1)
    allocate(T1%properties(NTDMARRAY))
    do i = 1, NTDMARRAY
       if ((i .eq. IFSDN .or. i .eq. IFSUP) .and. .not. T1%selfen) cycle
       call DQMC_TDM1_InitProp(T1, S, Gwrap, i)
    enddo
    
  end subroutine DQMC_TDM1_Init

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM1_InitFTw(T1)
    !
    ! Purpose
    ! =======
    ! Fill the matricess of fourier coefficients necessary to transform
    ! from imaginary time to imaginary frequency. Two matrices are
    ! needed: ftwbos for bosonic fields (spin, charge...) and 
    ! ftwfer for fermionic fields (particles)
    !
    ! Arguments
    ! =========
    !
    type(TDM1), intent(inout) :: T1      ! TDM to be freed

    integer :: iw, it, L
    real(wp) :: x, pi

    ! ... Executable ...

    if (.not.T1%compute) return

    L = T1%L
    pi = acos(-1.0_wp)

    allocate(T1%ftwbos(0:L-1,0:L-1))
    allocate(T1%ftwfer(0:L-1,0:L-1))

    do iw = 0, L-1
       do it = 0, L-1
          x = 2*it*iw*pi/L
          T1%ftwbos(iw,it) = T1%dtau * cmplx(cos(x),sin(x), kind=wp)
       enddo
    enddo

    do iw = 0, L-1
       do it = 0, L-1
          x = 2*it*(iw+0.5_wp)*pi/L
          T1%ftwfer(iw,it) = T1%dtau * cmplx(cos(x),sin(x), kind=wp)
       enddo
    enddo

  end subroutine DQMC_TDM1_InitFTw

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM1_InitProp(T1, S, Gwrap, iprop)
    use DQMC_Geom_Wrap
    !
    ! Purpose
    ! =======
    ! Initialize contents of tdmarray for iprop
    !
    ! Arguments
    ! =========
    type(TDM1), intent(inout) :: T1
    type(Struct), intent(in)  :: S
    type(GeomWrap), intent(in):: Gwrap
    integer, intent(in)       :: iprop

    integer :: nk, np, npp, nclass

    select case (iprop)

       case(IGFUN, IGFUP, IGFDN)

          nk     = Gwrap%RecipLattice%nclass_k
          nclass = S%nClass
          np     = Gwrap%lattice%natom
          npp    = (np*(np+1))/2

          nullify(T1%properties(iprop)%tlink)
          T1%properties(iprop)%n      =  S%nSite
          T1%properties(iprop)%nclass =  nclass
          T1%properties(iprop)%D      => S%D
          T1%properties(iprop)%F      => S%F
          T1%properties(iprop)%nk     =  nk
          T1%properties(iprop)%np     =  np
          T1%properties(iprop)%ftk    => Gwrap%RecipLattice%FourierC
          T1%properties(iprop)%ftw    => T1%ftwfer
          T1%properties(iprop)%phase  => S%gf_phase
          T1%properties(iprop)%clabel  => S%clabel
          allocate(T1%properties(iprop)%values(nclass,0:T1%L-1,T1%err))
          allocate(T1%properties(iprop)%valuesk(nk*npp,0:T1%L-1,T1%err))

       case(ISPXX, ISPZZ, IDENS, IPAIR)

          nk = Gwrap%GammaLattice%nclass_k
          nclass = S%nClass
          np     = Gwrap%lattice%natom
          npp    = (np*(np+1))/2

          nullify(T1%properties(iprop)%tlink)
          T1%properties(iprop)%n      =  S%nSite
          T1%properties(iprop)%nclass =  S%nClass
          T1%properties(iprop)%D      => S%D
          T1%properties(iprop)%F      => S%F
          T1%properties(iprop)%nk     =  Gwrap%GammaLattice%nclass_k
          T1%properties(iprop)%np     =  np
          T1%properties(iprop)%ftk    => Gwrap%GammaLattice%FourierC
          T1%properties(iprop)%ftw    => T1%ftwbos
          T1%properties(iprop)%phase  => S%chi_phase
          T1%properties(iprop)%clabel  => S%clabel
          allocate(T1%properties(iprop)%values(nclass,0:T1%L-1,T1%err))
          allocate(T1%properties(iprop)%valuesk(nk*npp,0:T1%L-1,T1%err))

       case(IFSUP, IFSDN)

          T1%properties(iprop)%tlink  => Gwrap%Hamilt%Uv
          T1%properties(iprop)%n      =  S%nSite
          T1%properties(iprop)%nclass =  S%nClass
          T1%properties(iprop)%D      => S%D
          T1%properties(iprop)%F      => S%F
          T1%properties(iprop)%nk     =  Gwrap%GammaLattice%nclass_k
          T1%properties(iprop)%np     =  Gwrap%lattice%natom
          nullify(T1%properties(iprop)%ftk)
          T1%properties(iprop)%ftw    => T1%ftwbos
          T1%properties(iprop)%phase  => S%chi_phase
          T1%properties(iprop)%clabel  => S%clabel
          allocate(T1%properties(iprop)%values(S%nClass,0:T1%L-1,T1%err))
          nullify(T1%properties(iprop)%valuesk)

       case(ICOND)

          nclass = 1

          ! "n" stores the number of sites in one primitive cell
          T1%properties(iprop)%n      =  Gwrap%lattice%natom
          ! I am reusing "nk" to store the number of cells in the cluster.
          T1%properties(iprop)%nk     =   S%nSite/Gwrap%lattice%natom
          ! I am setting nclass to 1. nclass corresponds, in this
          ! routine, to the length of the array in which the property
          ! is stored. By setting it to 1, I am assuming that only <jx jx>
          ! is desired. One may want to compute the entire tensor in which
          ! case this number should be increased accordingly. Or one may
          ! want to split the conductivity in different contributions e.g.
          ! in a multilayer system it is interesting to look at the contribution
          ! that each layer gives to the conductivity. As I said, for now,
          ! this is set to one : only one of the diagonal component of the
          ! tensor is computed.
          T1%properties(iprop)%nclass =  nclass
          ! Don't have much of a use for some of the pointer since we are
          ! here interested only in the k=0 component of the jj correlation.
          ! We can reuse such pointers to store the necessary info though.
          ! Here I set D to contain the list of sites making up the primitive
          ! links
          T1%properties(iprop)%D      => Gwrap%hamilt%plink
          ! And store the number of primitive links in np
          T1%properties(iprop)%np     =  Gwrap%hamilt%nplink
          ! There are noareal classes so F=1.
          allocate(T1%properties(iprop)%F(nclass))
          T1%properties(iprop)%F(1) = 1
          ! No use for F at this point or ftk at this point
          nullify(T1%properties(iprop)%ftk)
          T1%properties(iprop)%ftw     => T1%ftwbos
          ! Not sure this is needed. Probably not since current is bosonic.
          T1%properties(iprop)%phase   => S%chi_phase
          ! This would be needed if nclass is larger then one. I am assigning
          ! to it an empty string in this case.
          allocate(T1%properties(iprop)%clabel(nclass))
          T1%properties(iprop)%clabel(1)=" "
          ! Since number of class is 1
          allocate(T1%properties(iprop)%values(nclass,0:T1%L-1,T1%err))
          ! We do not need FT storage
          nullify(T1%properties(iprop)%valuesk)
          ! Amazingly, I realized I had no type_real pointer in this class.
          ! I had to add one just to store the tlink matrix. No biggie but
          ! inelegant.
          T1%properties(iprop)%tlink => Gwrap%hamilt%tlink

     end select

     T1%properties(iprop)%values  = 0.0_wp
     if(associated(T1%properties(iprop)%valuesk)) &
         T1%properties(iprop)%valuesk = 0.0_wp

  end subroutine DQMC_TDM1_InitProp

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM1_Free(T1)
    !
    ! Purpose
    ! =======
    !    This subroutine frees TDM1.
    !
    ! Arguments
    ! =========
    !
    type(TDM1), intent(inout) :: T1
    type(IndexList), pointer :: current, next
    integer :: i

    ! ... Executable ...

    if (.not.T1%compute) return
    do i = 1, NTDMARRAY
       if (associated(T1%properties(i)%values))  deallocate(T1%properties(i)%values)
       if (associated(T1%properties(i)%valuesk)) deallocate(T1%properties(i)%valuesk)
       nullify(T1%properties(i)%D)
       nullify(T1%properties(i)%F)
       nullify(T1%properties(i)%ftk)
       nullify(T1%properties(i)%ftw)
    enddo

    deallocate(T1%ftwbos)
    deallocate(T1%ftwfer)
    deallocate(T1%properties)

    do i = 1, size(T1%classToIdxPair)
        current => T1%classToIdxPair(i)%ptr
        if (.not. associated(current)) exit
        next => current%next
        do
            deallocate(current)
            if (.not. associated(next)) exit
            current => next
            if (.not. associated(current)) exit
            next => current%next
        enddo
    end do
    deallocate(T1%classToIdxPair)

  end subroutine DQMC_TDM1_Free

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM1_Meas(T1, tau)
    !
    ! Purpose
    ! =======
    !    This subroutine fills the bin. It assumes that tau%A_up and,
    !    when necessary, tau%A_dn have been filled.
    !
    ! Arguments
    ! =========
    !
    type(TDM1), intent(inout)   :: t1
    type(Gtau), intent(inout)   :: tau
    
    ! ... Local var ...
    integer  :: i, k, m, L, cnt, dt, i0, it, j0, jt, dtau, iprop
    real(wp) :: sgn, factor
    real(wp),pointer :: up0t(:,:)
    real(wp),pointer :: upt0(:,:)
    real(wp),pointer :: dn0t(:,:)
    real(wp),pointer :: dnt0(:,:)
    real(wp),pointer :: up00(:,:)
    real(wp),pointer :: uptt(:,:)
    real(wp),pointer :: dn00(:,:)
    real(wp),pointer :: dntt(:,:)
    real(wp), pointer :: values(:,:,:)

    if (.not.T1%compute) return
    ! ... executable ...
    cnt = 0
    L   = tau%L
    k   = mod(tau%north,2)
    m   = (tau%north-k) / 2

    upt0 => tau%upt0
    up0t => tau%up0t
    dnt0 => tau%dnt0
    dn0t => tau%dn0t
    up00 => tau%up00
    uptt => tau%uptt
    dn00 => tau%dn00
    dntt => tau%dntt

    blocks: do i0 = 1, tau%nb
       do dtau = 0, tau%nb-1
          it = mod(i0+dtau-1,tau%nb) + 1

          ! Stored value
          call DQMC_Gtau_DumpA(tau, TAU_UP, it, i0)
          if (tau%comp_dn .or. .not.tau%neg_u) &
             call DQMC_Gtau_DumpA(tau, TAU_DN, it, i0)
          
          jt = tau%it_up; j0 = tau%i0_up
          call DQMC_TDM1_Compute(T1, upt0, up0t, dnt0, dn0t, up00, uptt, dn00, dntt, jt, j0)

          ! Decrement index tau%it. If north is even do only north/2-1 decrements.
          do dt = 1, m-1+k
             call DQMC_change_gtau_time(tau, TPLUS, TAU_UP)
             if (tau%comp_dn) then
                call DQMC_change_gtau_time(tau, TPLUS, TAU_DN)
             elseif (.not.tau%neg_u) then
                call DQMC_Gtau_CopyUp(tau)
             endif

             jt = tau%it_up; j0 = tau%i0_up
             call DQMC_TDM1_Compute(T1, upt0, up0t, dnt0, dn0t, up00, uptt, dn00, dntt, jt, j0)
          enddo

          if (m .gt. 0) then
             call DQMC_Gtau_DumpA(tau, TAU_UP, it, i0)
             if (tau%comp_dn .or. .not.tau%neg_u) &
                call DQMC_Gtau_DumpA(tau, TAU_DN, it, i0)
          endif
          ! Increment index tau%it
          do dt = 1, m
             call DQMC_change_gtau_time(tau, TMINUS, TAU_UP)
             if (tau%comp_dn) then
                call DQMC_change_gtau_time(tau, TMINUS, TAU_DN)
             elseif (.not.tau%neg_u) then
                call DQMC_Gtau_CopyUp(tau)
             endif
             jt = tau%it_up; j0 = tau%i0_up
             call DQMC_TDM1_Compute(T1, upt0, up0t, dnt0, dn0t, up00, uptt, dn00, dntt, jt, j0)
                
          enddo

       enddo

       cnt = cnt + 1
    enddo blocks
 
    if (i0 .ne. tau%nb+1) then
       write(*,*) "Up and down time slices are mismatched. Stop"
       stop
    endif

    sgn = tau%sgnup * tau%sgndn
    do iprop = 1, NTDMARRAY
       if ((iprop .eq. IFSDN .or. iprop .eq. IFSUP) .and. .not. T1%selfen) cycle
       values => T1%properties(iprop)%values
       !$OMP PARALLEL DO SCHEDULE(STATIC), PRIVATE(it, factor)
       do i = 1, T1%properties(iprop)%nClass
           do it = 0, L-1
             factor = sgn/(T1%properties(iprop)%F(i)*cnt)
             values(i,it,T1%idx)   = values(i,it,T1%idx)   + factor*values(i,it,T1%tmp)
           end do
       end do
       !$OMP END PARALLEL DO
       values(:,:,T1%tmp)   = ZERO
    enddo

    T1%sgn(T1%idx) =  T1%sgn(T1%idx) + sgn
    T1%cnt = T1%cnt + 1

  end subroutine DQMC_TDM1_Meas
  
  subroutine DQMC_TDM1_Meas_Para(T1, pT1, tau, ptau)
    !
    ! Purpose
    ! =======
    !    This subroutine fills the bin. It assumes that tau%A_up and,
    !    when necessary, tau%A_dn have been filled.
    !
    ! Arguments
    ! =========
    !
    type(TDM1), intent(inout)   :: t1
    type(TDM1), intent(inout)   :: pT1(:)
    type(Gtau), intent(inout)   :: tau
    type(Gtau), intent(inout)   :: ptau(:)
    
    ! ... Local var ...
    integer  :: i, k, m, n, L, cnt, dt, i0, it, j0, jt, iprop
    real(wp) :: sgn, factor

    real(wp), pointer :: values(:,:,:)

    if (.not.T1%compute) return
    ! ... executable ...
    cnt = 0
    L   = tau%L
    k   = mod(tau%north,2)
    m   = (tau%north-k) / 2

    !$omp parallel private(it, i0, jt, j0)
    !$omp do
    do i = 1, tau%nb * tau%nb
       it = mod(i,tau%nb)+ 1
       i0 = (i-1)/tau%nb + 1
          ! Stored value
          call DQMC_Gtau_DumpA(ptau(i), TAU_UP, it, i0)
          if (tau%comp_dn .or. .not.tau%neg_u) &
             call DQMC_Gtau_DumpA(ptau(i), TAU_DN, it, i0)          
          jt = ptau(i)%it_up; j0 = ptau(i)%i0_up
          call DQMC_TDM1_Compute(pT1(i), ptau(i)%upt0, ptau(i)%up0t, ptau(i)%dnt0, &
             ptau(i)%dn0t, ptau(i)%up00, ptau(i)%uptt, ptau(i)%dn00, ptau(i)%dntt, jt, j0)

          ! Decrement index tau%it. If north is even do only north/2-1 decrements.
          do dt = 1, m-1+k
             call DQMC_change_gtau_time(ptau(i), TPLUS, TAU_UP)
             if (tau%comp_dn) then
                call DQMC_change_gtau_time(ptau(i), TPLUS, TAU_DN)
             elseif (.not.tau%neg_u) then
                call DQMC_Gtau_CopyUp(ptau(i))
             endif
             jt = ptau(i)%it_up; j0 = ptau(i)%i0_up
             call DQMC_TDM1_Compute(pT1(i), ptau(i)%upt0, ptau(i)%up0t, ptau(i)%dnt0, &
                ptau(i)%dn0t, ptau(i)%up00, ptau(i)%uptt, ptau(i)%dn00, ptau(i)%dntt, jt, j0)
          enddo

          if (m .gt. 0) then
             call DQMC_Gtau_DumpA(ptau(i), TAU_UP, it, i0)
             if (tau%comp_dn .or. .not.tau%neg_u) &
                call DQMC_Gtau_DumpA(ptau(i), TAU_DN, it, i0)
          endif
          ! Increment index tau%it
          do dt = 1, m
             call DQMC_change_gtau_time(ptau(i), TMINUS, TAU_UP)
             if (tau%comp_dn) then
                call DQMC_change_gtau_time(ptau(i), TMINUS, TAU_DN)
             elseif (.not.tau%neg_u) then
                call DQMC_Gtau_CopyUp(ptau(i))
             endif
             jt = ptau(i)%it_up; j0 = ptau(i)%i0_up
             call DQMC_TDM1_Compute(pT1(i), ptau(i)%upt0, ptau(i)%up0t, ptau(i)%dnt0, &
                ptau(i)%dn0t, ptau(i)%up00, ptau(i)%uptt, ptau(i)%dn00, ptau(i)%dntt, jt, j0)
          enddo
    enddo
    !$omp end do
    !$omp end parallel
    cnt = cnt + tau%nb

    ! Sum up the multi-threads value
    do n = 1, tau%nb * tau%nb
       do iprop = 1, NTDMARRAY
          if ((iprop .eq. IFSDN .or. iprop .eq. IFSUP) .and. .not. T1%selfen) cycle
          do i = 1, T1%properties(iprop)%nclass
             do it = 0, L-1
                   T1%properties(iprop)%values(i,it,T1%tmp) = T1%properties(iprop)%values(i,it,T1%tmp)+&
                      pT1(n)%properties(iprop)%values(i,it,T1%tmp)
             enddo
          enddo
       enddo
    enddo
 
    sgn = tau%sgnup * tau%sgndn
    do iprop = 1, NTDMARRAY
       if ((iprop .eq. IFSDN .or. iprop .eq. IFSUP) .and. .not. T1%selfen) cycle
       values => T1%properties(iprop)%values
       do it = 0, L-1
          do i = 1, T1%properties(iprop)%nClass
             factor = sgn/(T1%properties(iprop)%F(i)*cnt)
             values(i,it,T1%idx)   = values(i,it,T1%idx)   + factor*values(i,it,T1%tmp)
          end do
       end do
       values(:,:,T1%tmp)   = ZERO
       do n = 1, tau%nb*tau%nb
          values => pT1(n)%properties(iprop)%values
          values(:,:,T1%tmp) = ZERO
       enddo
    enddo

    T1%sgn(T1%idx) =  T1%sgn(T1%idx) + sgn
    T1%cnt = T1%cnt + 1

  end subroutine DQMC_TDM1_Meas_Para

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM1_Compute(T1, upt0, up0t, dnt0, dn0t, up00, uptt, dn00, dntt, it, i0)
    !
    ! Purpose
    ! =======
    !    This subroutine assembles the time dependent properties
    !    starting from the 1-body Green's function
    !
    ! Arguments
    ! =========
    !
    type(TDM1), intent(inout)    :: T1
    real(wp), intent(in)         :: up0t(:,:), upt0(:,:)
    real(wp), intent(in)         :: dnt0(:,:), dn0t(:,:)
    real(wp), intent(in)         :: up00(:,:), uptt(:,:)
    real(wp), intent(in)         :: dn00(:,:), dntt(:,:)
    integer, intent(in)          :: it, i0

    ! ... Local scalar ...

    integer  :: nclass, i, j, a, b, k, dt, dt1, dt2, iprop, dij
    real(wp), pointer :: value1(:)
    real(wp), pointer :: value2(:)
    real(wp) :: factor

    integer :: nsite, ic, jc, l1s, l1e, l2s, l2e
    real(wp) :: tl1up2up, tl1dn2dn, tl1up2dn, tl1dn2up, dxl1, dxl2, dxl1l2, Fabjt0
    integer  :: dims(3), indices(3)
    type(IndexList), pointer :: next

    ! ... Executable ...
    if (.not.T1%compute) return

    dt = it - i0
    if (dt .gt. 0) then
       ! it > i0
       dt1  =  dt
       dt2  =  T1%L - dt
       factor  = 0.25d0
    elseif (dt .lt. 0) then
       ! it < i0
       dt1  =  dt + T1%L
       dt2  = -dt
       factor = -0.25d0
    else
       dt1 = 0
       dt2 = 0
       factor = 0.5d0
    endif

    nclass = size(T1%classToIdxPair)
    dims(1:3) = T1%properties(IFSDN)%nk

    !$OMP PARALLEL DO SCHEDULE(STATIC), PRIVATE(next, i, j, value1, value2, indices, a, b, iprop)
    do k = 1, nclass
        do iprop = 1, NTDMARRAY
            if ((iprop .eq. IFSDN .or. iprop .eq. IFSUP) .and. .not. T1%selfen) cycle
            value1  => T1%properties(iprop)%values(:, dt1, T1%tmp)
            value2  => T1%properties(iprop)%values(:, dt2, T1%tmp)
            next => T1%classToIdxPair(k)%ptr
            do
                i = next%i
                j = next%j
                select case(iprop)
                        case (IGFUN)
                            value1(k)  = value1(k) + factor*(upt0(i,j) + dnt0(i,j))
                            if (dt .ne. 0) then
                                value2(k)  = value2(k) - factor*(up0t(i,j) + dn0t(i,j))
                            end if
                        case (IGFUP)
                            value1(k)  = value1(k) + 2*factor*upt0(i,j)
                            if (dt .ne. 0) then
                                value2(k)  = value2(k) - 2*factor*up0t(i,j)
                            end if
                        case (IGFDN)
                            value1(k)  = value1(k) + 2*factor*dnt0(i,j)
                            if (dt .ne. 0) then
                                value2(k)  = value2(k) - 2*factor*dn0t(i,j)
                            end if
                        case (ISPXX)
                            if (dt .ne. 0) then
                                value1(k)  = value1(k) - (up0t(j,i)*dnt0(i,j) &
                                        + dn0t(j,i)*upt0(i,j))/2
                                value2(k)  = value2(k) - (up0t(i,j)*dnt0(j,i) &
                                    + dn0t(i,j)*upt0(j,i))/2
                            else
                                value1(k)  = value1(k) - up0t(j,i)*dnt0(i,j) &
                                        - dn0t(j,i)*upt0(i,j)
                            end if
                        case (ISPZZ)
                            if (dt .ne. 0) then
                                value1(k)  = value1(k) - (up0t(j,i)*upt0(i,j) &
                                        + dn0t(j,i)*dnt0(i,j) - (uptt(i,i)-dntt(i,i))*(up00(j,j)-dn00(j,j)) )*0.5_wp
                                value2(k)  = value2(k) - (dn0t(i,j)*dnt0(j,i) &
                                    + up0t(i,j)*upt0(j,i) - (uptt(j,j)-dntt(j,j))*(up00(i,i)-dn00(i,i)) )*0.5_wp
                            else
                                value1(k)  = value1(k) - (up0t(j,i)*upt0(i,j) &
                                        + dn0t(j,i)*dnt0(i,j) - (uptt(i,i)-dntt(i,i))*(up00(j,j)-dn00(j,j)) )
                            end if
                        case (IDENS)
                            if (dt .ne. 0) then
                                value1(k)  = value1(k) - ( up0t(j,i)*upt0(i,j) &
                                        + dn0t(j,i)*dnt0(i,j) - (uptt(i,i)+dntt(i,i))*(up00(j,j)+dn00(j,j)) &
                                        - 2.0_wp*(1.0_wp - uptt(i,i) - up00(j,j) ) &
                                        - 2.0_wp*(1.0_wp - dntt(i,i) - dn00(j,j) ) )*0.5_wp
                                value2(k)  = value2(k) - (dn0t(i,j)*dnt0(j,i) &
                                        + up0t(i,j)*upt0(j,i) - (uptt(j,j)+dntt(j,j))*(up00(i,i)+dn00(i,i)) &
                                        - 2.0_wp*(1.0_wp - up00(i,i) - uptt(j,j) ) &
                                        - 2.0_wp*(1.0_wp - dn00(i,i) - dntt(j,j) ) )*0.5_wp
                            else
                                value1(k)  = value1(k) - ( up0t(j,i)*upt0(i,j) &
                                        + dn0t(j,i)*dnt0(i,j) - (uptt(i,i)+dntt(i,i))*(up00(j,j)+dn00(j,j)) &
                                        - 2.0_wp*(1.0_wp - uptt(i,i) - up00(j,j) ) &
                                        - 2.0_wp*(1.0_wp - dntt(i,i) - dn00(j,j) ) )
                            end if

                        case(IPAIR)
                            if (dt .ne. 0) then
                                value1(k)  = value1(k) + (upt0(i,j)*dnt0(i,j) &
                                        + dn0t(j,i)*up0t(j,i))/2
                                value2(k)  = value2(k) + (upt0(j,i)*dnt0(j,i) &
                                        + dn0t(i,j)*up0t(i,j))/2
                            else
                                value1(k)  = value1(k) + upt0(i,j)*dnt0(i,j) &
                                        + dn0t(j,i)*up0t(j,i)
                            end if
                        case (IFSUP)
                            if (i .eq. j) then
                                dij = 1
                            else
                                dij = 0
                            end if
                            if (dt .ne. 0) then
                                value1(k)  = value1(k) - T1%properties(IFSUP)%tlink(0,0)**2 * ((1.0 - dntt(i,i)) * upt0(i,j) * (1.0 - dn00(j,j)) + &
                                        (dij - dn0t(j,i)) * dnt0(i,j) * upt0(i,j)) * 0.5
                                value2(k)  = value2(k) - T1%properties(IFSUP)%tlink(0,0)**2 * ((1.0 - dntt(j,j)) * upt0(j,i) * (1.0 - dn00(i,i)) + &
                                        (dij - dn0t(i,j)) * dnt0(j,i) * upt0(j,i)) * 0.5
                            else
                                value1(k)  = value1(k) - T1%properties(IFSUP)%tlink(0,0)**2 * (1.0 - dntt(i,i)) * upt0(i,j) * (1.0 - dn00(j,j)) + &
                                        (dij - dn0t(j,i)) * dnt0(i,j) * upt0(i,j)
                            end if
                        case (IFSDN)
                            if (i .eq. j) then
                                dij = 1
                            else
                                dij = 0
                            end if
                            if (dt .ne. 0) then
                                value1(k)  = value1(k) - T1%properties(IFSUP)%tlink(0,0)**2 * ((1.0 - uptt(i,i)) * dnt0(i,j) * (1.0 - up00(j,j)) + &
                                        (dij - up0t(j,i)) * upt0(i,j) * dnt0(i,j)) * 0.5
                                value2(k)  = value2(k) - T1%properties(IFSUP)%tlink(0,0)**2 * ((1.0 - uptt(j,j)) * dnt0(j,i) * (1.0 - up00(i,i)) + &
                                        (dij - up0t(i,j)) * upt0(j,i) * dnt0(j,i)) * 0.5
                            else
                                value1(k)  = value1(k) - T1%properties(IFSUP)%tlink(0,0)**2 * (1.0 - uptt(i,i)) * dnt0(i,j) * (1.0 - up00(j,j)) + &
                                        (dij - up0t(j,i)) * upt0(i,j) * dnt0(i,j)
                            end if
                end select
                if (.not. associated(next%next)) then
                    exit
                end if
                next => next%next
            end do
        end do
    end do
    !$OMP END PARALLEL DO

    if (dt .ne. 0) then

       ! Uniform (q=0) current structure factor. This quantity can be used to get approximated dc conductivity.
       value1  => T1%properties(ICOND)%values(:, dt1, T1%tmp)
       value2  => T1%properties(ICOND)%values(:, dt2, T1%tmp)

       ! J-J correlation is following Simone's trick:
       
       ! n:  the number of sites in one primitive cell
       ! nk: the number of cells in the cluster
       ! D:  the list of sites making up the primitive links
       ! np: the number of all links, instead of primitive links

       ! Double Loops over all links:

       nsite  = T1%properties(ICOND)%n*T1%properties(ICOND)%nk

       ! Call nk the number of unit cells, than the first nk entries of plink
       ! are the translation of the first primitive link. Likewise from nk+1 to 2nk
       ! you have all the links obtained by translation of the second link and so on

       do i = 0,  T1%properties(ICOND)%np - 1
          do j = 0,  T1%properties(ICOND)%np - 1

             ! find out if this type of link contribute to conductivity
             dxl1  = T1%properties(ICOND)%tlink(3, i)
             dxl2  = T1%properties(ICOND)%tlink(3, j)

             if (abs(dxl1)<1.0E-5 .or. abs(dxl2)<1.0E-5) then
                value1(1) = value1(1)
                value2(1) = value2(1)
             else

                ! x-coordinate (for JxJx component of J-J tensor) difference
                dxl1l2 = dxl1*dxl2*0.5_wp/nsite
                tl1up2up = T1%properties(ICOND)%tlink(1, i)*T1%properties(ICOND)%tlink(1, j)
                tl1dn2dn = T1%properties(ICOND)%tlink(2, i)*T1%properties(ICOND)%tlink(2, j)
                tl1up2dn = T1%properties(ICOND)%tlink(1, i)*T1%properties(ICOND)%tlink(2, j)
                tl1dn2up = T1%properties(ICOND)%tlink(2, i)*T1%properties(ICOND)%tlink(1, j)

                do ic = 0,  T1%properties(ICOND)%nk - 1
                  do jc = 0,  T1%properties(ICOND)%nk - 1
                     l1s = T1%properties(ICOND)%D(1, i*T1%properties(ICOND)%nk+ic) + 1
                     l1e = T1%properties(ICOND)%D(2, i*T1%properties(ICOND)%nk+ic) + 1
                     l2s = T1%properties(ICOND)%D(1, j*T1%properties(ICOND)%nk+jc) + 1
                     l2e = T1%properties(ICOND)%D(2, j*T1%properties(ICOND)%nk+jc) + 1

                     ! jx-jx for links l and l':

                     value1(1) = value1(1) - dxl1l2*tl1up2up*((uptt(l1e,l1s)-uptt(l1s,l1e))*(up00(l2e,l2s)-up00(l2s,l2e)) ) &
                                           - dxl1l2*tl1up2dn*((uptt(l1e,l1s)-uptt(l1s,l1e))*(dn00(l2e,l2s)-dn00(l2s,l2e)) ) &
                                           - dxl1l2*tl1dn2up*((dntt(l1e,l1s)-dntt(l1s,l1e))*(up00(l2e,l2s)-up00(l2s,l2e)) ) &
                                           - dxl1l2*tl1dn2dn*((dntt(l1e,l1s)-dntt(l1s,l1e))*(dn00(l2e,l2s)-dn00(l2s,l2e)) )
  
                     ! cross terms
                     value1(1) = value1(1) + dxl1l2*tl1up2up*( upt0(l1e,l2s)*up0t(l2e,l1s) + upt0(l1s,l2e)*up0t(l2s,l1e) &
                                           - upt0(l1e,l2e)*up0t(l2s,l1s) - upt0(l1s,l2s)*up0t(l2e,l1e) )
                     value1(1) = value1(1) + dxl1l2*tl1dn2dn*( dnt0(l1e,l2s)*dn0t(l2e,l1s) + dnt0(l1s,l2e)*dn0t(l2s,l1e) &
                                           - dnt0(l1e,l2e)*dn0t(l2s,l1s) - dnt0(l1s,l2s)*dn0t(l2e,l1e) )

                     ! value2: see the rules at the beginning of routine
                     value2(1) = value2(1) - dxl1l2*tl1up2up*((up00(l1e,l1s)-up00(l1s,l1e))*(uptt(l2e,l2s)-uptt(l2s,l2e)) ) &
                                           - dxl1l2*tl1up2dn*((up00(l1e,l1s)-up00(l1s,l1e))*(dntt(l2e,l2s)-dntt(l2s,l2e)) ) &
                                           - dxl1l2*tl1dn2up*((dn00(l1e,l1s)-dn00(l1s,l1e))*(uptt(l2e,l2s)-uptt(l2s,l2e)) ) &
                                           - dxl1l2*tl1dn2dn*((dn00(l1e,l1s)-dn00(l1s,l1e))*(dntt(l2e,l2s)-dntt(l2s,l2e)) )

                     ! cross terms
                     value2(1) = value2(1) + dxl1l2*tl1up2up*( up0t(l1e,l2s)*upt0(l2e,l1s) + up0t(l1s,l2e)*upt0(l2s,l1e) &
                                           - up0t(l1e,l2e)*upt0(l2s,l1s) - up0t(l1s,l2s)*upt0(l2e,l1e) )
                     value2(1) = value2(1) + dxl1l2*tl1dn2dn*( dn0t(l1e,l2s)*dnt0(l2e,l1s) + dn0t(l1s,l2e)*dnt0(l2s,l1e) &
                                           - dn0t(l1e,l2e)*dnt0(l2s,l1s) - dn0t(l1s,l2s)*dnt0(l2e,l1e) )
                  end do
                 end do
               endif
             end do
           end do

    else
       value1  => T1%properties(ICOND)%values(:, dt1, T1%tmp)

       nsite  = T1%properties(ICOND)%n*T1%properties(ICOND)%nk

       do i = 0,  T1%properties(ICOND)%np - 1
          do j = 0,  T1%properties(ICOND)%np - 1

             ! find out if this type of link contribute to conductivity
             dxl1  = T1%properties(ICOND)%tlink(3, i)
             dxl2  = T1%properties(ICOND)%tlink(3, j)
             if (abs(dxl1)<1.0E-5 .or. abs(dxl2)<1.0E-5) then
                value1(1) = value1(1)
             else

             dxl1l2 = dxl1*dxl2/nsite
             tl1up2up = T1%properties(ICOND)%tlink(1, i)*T1%properties(ICOND)%tlink(1, j)
             tl1dn2dn = T1%properties(ICOND)%tlink(2, i)*T1%properties(ICOND)%tlink(2, j)
             tl1up2dn = T1%properties(ICOND)%tlink(1, i)*T1%properties(ICOND)%tlink(2, j)
             tl1dn2up = T1%properties(ICOND)%tlink(2, i)*T1%properties(ICOND)%tlink(1, j)

             do ic = 0,  T1%properties(ICOND)%nk - 1
               do jc = 0,  T1%properties(ICOND)%nk - 1
                  l1s = T1%properties(ICOND)%D(1, i*T1%properties(ICOND)%nk+ic) + 1
                  l1e = T1%properties(ICOND)%D(2, i*T1%properties(ICOND)%nk+ic) + 1
                  l2s = T1%properties(ICOND)%D(1, j*T1%properties(ICOND)%nk+jc) + 1
                  l2e = T1%properties(ICOND)%D(2, j*T1%properties(ICOND)%nk+jc) + 1

                  ! jx-jx for links l and l':
                  ! include x-coordinate (for JxJx component of J-J tensor) difference
  
                  value1(1) = value1(1) - dxl1l2*tl1up2up*((uptt(l1e,l1s)-uptt(l1s,l1e))*(up00(l2e,l2s)-up00(l2s,l2e)) ) &
                                     - dxl1l2*tl1up2dn*((uptt(l1e,l1s)-uptt(l1s,l1e))*(dn00(l2e,l2s)-dn00(l2s,l2e)) ) &
                                     - dxl1l2*tl1dn2up*((dntt(l1e,l1s)-dntt(l1s,l1e))*(up00(l2e,l2s)-up00(l2s,l2e)) ) &
                                     - dxl1l2*tl1dn2dn*((dntt(l1e,l1s)-dntt(l1s,l1e))*(dn00(l2e,l2s)-dn00(l2s,l2e)) )

                  ! cross terms
                  value1(1) = value1(1) + dxl1l2*tl1up2up*( upt0(l1e,l2s)*up0t(l2e,l1s) + upt0(l1s,l2e)*up0t(l2s,l1e) &
                                     - upt0(l1e,l2e)*up0t(l2s,l1s) - upt0(l1s,l2s)*up0t(l2e,l1e) )
                  value1(1) = value1(1) + dxl1l2*tl1dn2dn*( dnt0(l1e,l2s)*dn0t(l2e,l1s) + dnt0(l1s,l2e)*dn0t(l2s,l1e) &
                                - dnt0(l1e,l2e)*dn0t(l2s,l1s) - dnt0(l1s,l2s)*dn0t(l2e,l1e) )
             end do
            end do
           endif
        end do
       end do


    endif

  end subroutine DQMC_TDM1_Compute

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM1_Avg(T1)
    !
    ! Purpose
    ! =======
    !    This subroutine average properties in a bin and
    !    increment the bin count (idx).
    !
    ! Arguments
    ! =========
    !
    type(TDM1), intent(inout) :: T1                 ! T1

    ! ... local scalar ...
    integer  :: nl, idx, i
    real(wp) :: factor

    ! ... Executable ...
    if (.not.T1%compute) return
    idx    = T1%idx
    factor = ONE/T1%cnt

    ! Compute average on Green's function
    do i = 1, NTDMARRAY
       if ((i .eq. IFSDN .or. i .eq. IFSUP) .and. .not. T1%selfen) cycle
       nl = T1%properties(i)%nClass * T1%L
       call dscal(nl, factor, T1%properties(i)%values(:,0,idx), 1)
    enddo

    T1%sgn(idx) = T1%sgn(idx)*factor
    T1%cnt = 0
    T1%idx = T1%idx + 1

  end subroutine DQMC_TDM1_Avg

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM1_GetErr(T1)
    use dqmc_mpi
    !
    ! Purpose
    ! =======
    !    This subroutine compute the error in tdm using the jackknife
    !
    ! Arguments
    ! =========
    !
    type(TDM1), intent(inout) :: T1                 ! T1

    ! ... local scalar ...
    integer   :: i, j, iprop
    integer   :: nproc, n, avg, err, mpi_err
    real(wp)  :: sum_sgn, sgn(T1%nBin), y(T1%nBin), data(T1%nBin)
    real(wp)  :: average, error

    if (.not.T1%compute) return
    ! ... Executable ...
    nproc  = qmc_sim%size
    n      = T1%nbin
    avg    = T1%avg
    err    = T1%err

    if (nproc .eq. 1) then

       data = T1%sgn(1:n)
       call DQMC_JackKnife(n, T1%sgn(avg), T1%sgn(err), data , &
            y, sgn, sum_sgn)

       do iprop = 1, NTDMARRAY
          if ((iprop .eq. IFSDN .or. iprop .eq. IFSUP) .and. .not. T1%selfen) cycle
          !$OMP PARALLEL DO SCHEDULE(STATIC), PRIVATE(j, data, average, error), FIRSTPRIVATE(y)
          do i = 1, T1%properties(iprop)%nClass
             do j = 0, T1%L-1
                data =  T1%properties(iprop)%values(i, j, 1:n)
                call DQMC_SignJackKnife(n, average, error, data, y, sgn, sum_sgn)
                T1%properties(iprop)%values(i, j, avg) = average
                T1%properties(iprop)%values(i, j, err) = error
             enddo
          end do
          !$OMP END PARALLEL DO
       enddo

    else

       mpi_err = 0

#      ifdef _QMC_MPI
          
          !Average sign
          call mpi_allreduce(T1%sgn(1), T1%sgn(avg), 1, mpi_double, &
             mpi_sum, mpi_comm_world, mpi_err)

          !Average properties
          do iprop = 1, NTDMARRAY
             if ((iprop .eq. IFSDN .or. iprop .eq. IFSUP) .and. .not. T1%selfen) cycle
             binptr => T1%properties(iprop)%values(:,:,1)
             aveptr => T1%properties(iprop)%value(:,:,avg)
             n = T1%properties(iprop)%nClass * T1%L
             call mpi_allreduce(binptr, aveptr, n, mpi_double, &
                mpi_sum, mpi_comm_world, mpi_err)
          enddo

          !Compute average over n-1 processors
          do iprop = 1, NTDMARRAY
             if ((iprop .eq. IFSDN .or. iprop .eq. IFSUP) .and. .not. T1%selfen) cycle
             binptr => T1%properties(iprop)%values(:,:,1)
             aveptr => T1%properties(iprop)%values(:,:,avg)
             binptr = (aveptr - binptr) / dble(nproc - 1)
          enddo
          T1%sgn(1)   = (T1%sgn(avg) - T1%sgn(1)) / dble(nproc - 1)

          !Store average amongst all processors
          do iprop = 1, NTDMARRAY
             if ((iprop .eq. IFSDN .or. iprop .eq. IFSUP) .and. .not. T1%selfen) cycle
             aveptr => T1%properties(iprop)%values(:,:,avg)
             aveptr =  aveptr / T1%sgn(avg) 
          enddo
          T1%sgn(:,avg)     = T1%sgn(:,avg) / dble(nproc)

          !Store jackknife in the processor bin
          do iprop = 1, NTDMARRAY
             if ((iprop .eq. IFSDN .or. iprop .eq. IFSUP) .and. .not. T1%selfen) cycle
             binptr => T1%properties(iprop)%values(:,:,1)
             binptr =  binptr / T1%sgn(1) 
          enddo

          !Compute error
          do iprop = 1, NTDMARRAY
             if ((iprop .eq. IFSDN .or. iprop .eq. IFSUP) .and. .not. T1%selfen) cycle
             binptr => T1%properties(iprop)%values(:,:,1)
             errptr => T1%properties(iprop)%values(:,:,err)
             n = T1%properties(iprop)%nClass * T1%L
             call mpi_allreduce(binptr**2, errptr, n, mpi_double, &
                 mpi_sum, mpi_comm_world, mpi_err)
             errptr = errptr / dble(nproc) - aveptr**2 
             errptr = sqrt(errptr * dble(nproc-1))
          enddo

#      endif

    endif

  end subroutine DQMC_TDM1_GetErr

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM1_Print(T1, OPT)
    use dqmc_mpi
    !
    ! Purpose
    ! =======
    !    This subroutine prints properties to file
    !
    ! Arguments
    ! =========
    !
    type(TDM1), intent(in)   :: T1                 ! T1
    integer, intent(in)      :: OPT

    integer             :: i, j, iprop
    real(wp)            :: tmp(T1%L, 2)
    character(len=15)   :: label(T1%L)
    character(len=slen) :: title

    ! ... Executable ...
    if (.not.T1%compute) return

    if (qmc_sim%rank .ne. 0) return

    do j = 1, T1%L
       write(label(j),'(f15.8)') (j-1)*T1%dtau
    enddo

    do iprop = 1, NTDMARRAY
       if ((iprop .eq. IFSDN .or. iprop .eq. IFSUP) .and. .not. T1%selfen) cycle
       do i = 1, T1%properties(iprop)%nclass
          do j = 0, T1%L-1
             tmp(j+1, 1:2) = T1%properties(iprop)%values(i, j, T1%avg:T1%err)
          enddo
          title=pname(iprop)//" "//trim(adjustl(T1%properties(iprop)%clabel(i)))
          call DQMC_Print_Array(0, T1%L , title, label, tmp(:, 1:1), tmp(:, 2:2), OPT)
          write(OPT,'(1x)')
       enddo
    enddo

  end subroutine DQMC_TDM1_Print

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM1_GetKFT(T1)

    type(tdm1), intent(inout) :: T1

    integer :: ip, it, n, nclass, np, nk, ibin
    integer, pointer :: class(:,:)
    complex(wp), pointer :: wgtftk(:,:)
    integer, pointer :: phase(:,:)

    real(wp), pointer :: value(:)
    complex(wp), pointer :: valuet(:)

    if (.not.T1%compute) return
 
    !Loop over properties to Fourier transform
    do ip = 1, NTDMARRAY

       if (.not.associated(T1%properties(ip)%valuesk)) cycle

       ! Aliases
       n        =  T1%properties(ip)%n
       nclass   =  T1%properties(ip)%nclass
       np       =  T1%properties(ip)%np
       nk       =  T1%properties(ip)%nk
       class    => T1%properties(ip)%D
       wgtftk   => T1%properties(ip)%ftk
       phase    => T1%properties(ip)%phase

       !Fourier transform each bin and average
       do ibin = T1%avg, 1, -1

          ! More aliases
          do it = 0, T1%L-1
             value  =>  T1%properties(ip)%values(:,it,ibin)
             valuet =>  T1%properties(ip)%valuesk(:,it,ibin)
             call dqmc_getFTk(value, n, nclass, class, np, nk, wgtftk, phase, valuet)
          enddo

       enddo ! Loop over bins

    enddo ! Loop over properties

  end subroutine DQMC_TDM1_GetKFT

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM1_GetErrKFT(T1)

    use DQMC_MPI

    type(tdm1), intent(inout) :: T1

    integer :: ip, it, n, nproc, i

    complex(wp), pointer  :: average(:)
    complex(wp), pointer  :: binval(:) 
    complex(wp), pointer  :: error(:)  
    complex(wp), pointer  :: temp(:)   
 
    !Loop over properties to Fourier transform
    nproc = qmc_sim%size

    if (.not.T1%compute) return

    if (nproc .eq. 1) then

       do ip = 1, NTDMARRAY

          if (.not.associated(T1%properties(ip)%valuesk)) cycle

          do it = 0, T1%L-1

             average  => T1%properties(ip)%valuesk(:,it,T1%avg)
             error    => T1%properties(ip)%valuesk(:,it,T1%err)

             !Fourier transform each bin and average
             do i = 1, T1%nbin

                binval => T1%properties(ip)%valuesk(:,it,i)
                error  = error  +  cmplx((real(average-binval))**2,(aimag(average-binval))**2, kind=wp)

             enddo ! Loop over bins

             error  = (T1%nbin-1)*error/T1%nbin
             error = cmplx(sqrt(real(error)),sqrt(aimag(error)), kind=wp)
  
          enddo

       enddo ! Loop over properties

    else

       do ip = 1, NTDMARRAY

          if (.not.associated(T1%properties(ip)%valuesk)) cycle
    
          n = T1%properties(ip)%nk * T1%properties(ip)%np
          allocate(temp(n))
          
          do it = 0, T1%L-1
             average  => T1%properties(ip)%valuesk(:,it,T1%avg)
             error    => T1%properties(ip)%valuesk(:,it,T1%err)
             binval   => T1%properties(ip)%valuesk(:,it,1)
             temp     =  cmplx((real(average-binval))**2,(aimag(average-binval))**2, kind=wp)
#            ifdef _QMC_MPI
             call mpi_allreduce(temp, error, n, mpi_double, mpi_sum, mpi_comm_world, i)
#            endif
             error  = (nproc-1)*error/nproc
             error = cmplx(sqrt(real(error)),sqrt(aimag(error)), kind=wp)
          enddo

          deallocate(temp)

       enddo ! Loop over properties

    endif


  end subroutine DQMC_TDM1_GetErrKFT

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM1_PrintKFT(T1, OPT)
    use dqmc_mpi
    !
    ! Purpose
    ! =======
    !    This subroutine prints properties to file
    !
    ! Arguments
    ! =========
    !
    type(TDM1), intent(in)   :: T1                 ! T1
    integer, intent(in)      :: OPT

    integer             :: i, j, k, ip, jp, iprop
    integer             :: np, npp
    complex(wp)         :: tmp(T1%L, 2)
    character(len=15)   :: label(T1%L)
    character(len=60) :: title

    ! ... Executable ...
    if (.not.T1%compute) return

    if (qmc_sim%rank .ne. 0) return

    do j = 1, T1%L
       write(label(j),'(f15.8)') (j-1)*T1%dtau
       label(j) = adjustl(label(j))
    enddo

    do iprop = 1, NTDMARRAY
       if (.not.associated(T1%properties(iprop)%valuesk)) cycle
       np = T1%properties(iprop)%np
       npp = (np*(np+1))/2
       do k = 1, T1%properties(iprop)%nk
          i = (k-1)*npp
          do ip = 1, np
             do jp = ip, np
                i = i + 1
                do j = 0, T1%L-1
                   tmp(j+1, 1:2) = T1%properties(iprop)%valuesk(i, j, T1%avg:T1%err)
                enddo
                write(title,'(A,i3,A,i3,A,i3,A)') 'k=',k,'   pair=',ip,',',jp
                title=pname(iprop)//" "//trim(adjustl(title))
                call DQMC_Print_Array(0, T1%L , title, label, tmp(:, 1:1), tmp(:, 2:2), OPT)
                write(OPT,'(1x)')
             enddo
          enddo
       enddo
    enddo

  end subroutine DQMC_TDM1_PrintKFT

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM1_ImprovedSelfEnergy(T1, tau, OPT)
    type(TDM1), intent(in)    :: T1
    type(gtau), intent(inout) :: tau
    integer, intent(in)       :: OPT

    integer :: L, n, a, b, i, j, ab, t, nclass, indices(3), dims(3)
    !complex(wp), allocatable :: work(:)
    real(wp), allocatable :: S(:,:)
    !complex(wp), allocatable :: G(:,:), GS(:,:), S(:,:)
    !integer, allocatable     :: ipiv(:)
    !integer     :: info
    real(wp)            :: tmp(T1%L, 2)
    character(len=15)   :: label(T1%L)
    character(len=slen) :: title

    L      =  T1%L
    n     =  T1%properties(IGFUN)%n
    nclass = T1%properties(IFSUP)%nk

    !allocate(work(n))
    !allocate(ipiv(n))
    !allocate(G(n,n))
    !allocate(GS(n,n))
    !allocate(S(0:L-1, nclass))

    !dims(1:3) = n

!    !$OMP PARALLEL DO SCHEDULE(STATIC), PRIVATE(a, b, i, j, ab, indices)
!    do t = 0, L - 1
!        write(label(t + 1),'(f15.8)') t * T1%dtau
!        S(t,:) = 0.0
!        do a = 1, n
!            do b = 1, n
!                indices(2) = b
!                ab = T1%properties(IFSUP)%D(a,b)
!                do i = 1, n
!                    indices(1) = i
!                    do j = 1, n
!                        indices(3) = j
!                        S(t,ab) = S(t,ab) + 0.5_wp * (T1%properties(IFSUP)%tlink(j,b) + T1%properties(IFSUP)%tlink(b,j)) * &
!                                T1%properties(IFSUP)%values(DQMC_TDM1_GetUniqueIndexOfTuple(indices, dims, 3), t, T1%avg) / &
!                                T1%properties(IGFUP)%values(T1%properties(IGFUP)%D(a,i), t, T1%avg)
!                    end do
!                end do
!            end do
!        end do
!    end do
!    !$OMP END PARALLEL DO

!    !$OMP PARALLEL DO SCHEDULE(STATIC), PRIVATE(a, b, i, j, ab, indices)
!    do t = 0, L - 1
!        write(label(t + 1),'(f15.8)') t * T1%dtau
!        S(t,:) = 0.0
!        do a = 1, n
!            do b = 1, n
!                indices(2) = b
!                ab = T1%properties(IFSUP)%D(a,b)
!                do i = 1, n
!                    indices(1) = i
!                    do j = 1, n
!                        indices(3) = j
!                        S(t,ab) = S(t,ab) + 0.5_wp * (T1%properties(IFSUP)%tlink(j,b) + T1%properties(IFSUP)%tlink(b,j)) * &
!                                T1%properties(IFSUP)%values(DQMC_TDM1_GetUniqueIndexOfTuple(indices, dims, 3), t, T1%avg) / &
!                                T1%properties(IGFUP)%values(T1%properties(IGFUP)%D(a,i), t, T1%avg)
!                    end do
!                end do
!            end do
!        end do
!    end do
!    !$OMP END PARALLEL DO

    !do t = 0, L - 1
    !   write(label(t + 1),'(f15.8)') t * T1%dtau
       ! Fill matrix. Note that G is complex symmetric. Not hermitian.
    !   do i = 1, n
    !       do j = 1, n
    !           G(i,j) = T1%properties(IGFUP)%values(T1%properties(IGFUP)%D(i,j), t, T1%avg)
    !           GS(i,j) = T1%properties(IFSUP)%values(T1%properties(IFSUP)%D(i,j), t, T1%avg)
    !       end do
    !   end do
       ! solve the linear equation system G*Sigma = B, where G is the Green's function
       ! Sigma the self energy matrix and B GSigma calculated through the correlation
       ! function F = -<T c_a(tau)c^+_b(tau')c^+_j(tau')c_j(tau') >
    !   call zgetrf(n, n, G(1:n,1:n), n, ipiv(1:n), info)
    !   call zgetrs('N', n, n, G(1:n,1:n), n, ipiv(1:n), GS(1:n, 1:n), n, info)
    !   S(t,:) = 0.0
    !   do i = 1, n
    !       do j = 1, n
    !           S(t,T1%properties(IGFUP)%D(i,j)) = GS(i,j)
    !       end do
    !   end do
    !enddo

    do i = 1, nclass
        do j = 0, L - 1
            tmp(j+1,1) = S(j, i)
            tmp(j+1,2) = 0.0
        enddo
        title = "SelfEnergy up "//trim(adjustl(T1%properties(IFSUP)%clabel(i)))
        call DQMC_Print_Array(0, L , title, label, tmp(:, 1:1), tmp(:, 2:2), OPT)
        write(OPT,'(1x)')
    enddo

    !deallocate(work)
    !deallocate(ipiv)
    !deallocate(G)
    !deallocate(GS)
    deallocate(S)

  end subroutine DQMC_TDM1_ImprovedSelfEnergy

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM1_SelfEnergy(T1, tau, OPT)

    use DQMC_MPI

    type(TDM1), intent(in)    :: T1
    type(gtau), intent(inout) :: tau
    integer, intent(in)       :: OPT

    real(wp),    allocatable  :: g0tau(:,:), tdmg0(:,:)
    complex(wp), allocatable  :: tdmg0k(:,:), tdmgk(:,:)
    complex(wp), allocatable  :: tdmg0kw(:,:,:), tdmgkw(:,:,:)
    complex(wp), allocatable  :: avgSE(:,:,:), errSE(:,:,:), binSE(:,:,:)

    integer :: i, j, k, h, m
    integer :: L, n, nclass, np, nk, npp
    integer,     pointer :: class(:,:)
    integer,     pointer :: ph(:,:)   
    complex(wp), pointer :: ftk(:,:)  
    complex(wp), pointer :: ftw(:,:)  

    complex(wp)          :: tmp(T1%L,2)
    character(len=15)    :: label(T1%L)
    character(len=50)    :: title

    integer, parameter  :: gflist(3) = (/IGFUN, IGFUP, IGFDN/)
    integer, parameter  :: splist(3) = (/  0, TAU_UP, TAU_DN/)
    real(wp), parameter :: pi = 3.1415926535898

    if (.not.T1%compute) return

    L      =  T1%L
    n      =  T1%properties(IGFUN)%n
    nclass =  T1%properties(IGFUN)%nclass
    np     =  T1%properties(IGFUN)%np
    nk     =  T1%properties(IGFUN)%nk
    class  => T1%properties(IGFUN)%D
    ftk    => T1%properties(IGFUN)%ftk
    ftw    => T1%properties(IGFUN)%ftw
    ph     => T1%properties(IGFUN)%phase

    npp = np*(np+1)/2

    do j = 1, L
       write(label(j),'(f15.8)') (2*j-1)*pi/(T1%dtau*L)
       label(j) = adjustl(label(j))
    enddo

    !non-interacting green's function
    allocate(g0tau(n,n))
    allocate(tdmg0(nclass,0:L-1))
    allocate(tdmg0k(nk*npp,0:L-1))

    !temporary storage
    allocate(tdmgk(npp,0:L-1))

    !frequency dependent green's functions
    allocate(tdmg0kw(np,np,0:L-1))
    allocate(tdmgkw(np,np,0:L-1))

    !self energy
    allocate(avgSE(np,np,0:L-1))
    allocate(errSE(np,np,0:L-1))
    allocate(binSE(np,np,0:L-1))

    do h = 1, 3

       ! Get G for the non-interacting system
       tdmg0 = 0.0_wp
       do m = 0, L-1
          call dqmc_Gtau_GetG0(n, tau, splist(h), m, g0tau)
          do i = 1, n
             do j = 1, n
                k = class(i,j)
                tdmg0(k,m) = tdmg0(k,m) + g0tau(i,j) 
             enddo
          enddo
       enddo
       do k = 1, T1%properties(h)%nClass
          tdmg0(k,:) = tdmg0(k,:) / T1%properties(IGFUN)%F(k)
       enddo

       ! Get G in the k-space
       do m = 0, L-1
          call dqmc_getFTk(tdmg0(:,m), n, nclass, class, np, nk, ftk, ph, tdmg0k(:,m))
       enddo

       do k = 1, nk
 
          i = (k-1) * npp + 1
          j = k * npp

          ! Transform G0 from tau to iwn
          tdmgk = tdmg0k(i:j,0:L-1)
          call convert_to_iwn(tdmgk, tdmg0kw)
          call invertG(tdmg0kw)

          ! Transform G from tau to iwn for average
          tdmgk = T1%properties(gflist(h))%valuesk(i:j,0:L-1,T1%avg)
          call convert_to_iwn(tdmgk, tdmgkw)
          call invertG(tdmgkw)

          !Compute average self-energy
          avgSE = tdmg0kw - tdmgkw

          errSE = ZERO
          do m = 0, T1%nbin-1

             ! Transform G from tau to iwn for bin "m"
             tdmgk = T1%properties(gflist(h))%valuesk(i:j,0:L-1,m+1)
             call convert_to_iwn(tdmgk, tdmgkw)
             call invertG(tdmgkw)

             ! Compute self-energy for bin
             binSE = tdmg0kw - tdmgkw
             if (qmc_sim%size .eq. 1) &
                errSE = errSE + cmplx((real(binSE-avgSE))**2,(aimag(binSE-avgSE))**2, kind=wp)
  
          enddo 

          if (qmc_sim%size .gt. 1) then

             m = qmc_sim%size
             ! Reuse tdmgkw for temporary storage
             tdmgkw =  cmplx((real(binSE-avgSE))**2,(aimag(binSE-avgSE))**2, kind=wp)
#            ifdef _QMC_MPI
             call mpi_allreduce(tdmgkw, errSE, n, mpi_double_complex, mpi_sum, mpi_comm_world, i)
#            endif

          endif

          errSE = cmplx(sqrt(real(errSE)),sqrt(aimag(errSE)), kind=wp) * sqrt(dble(m-1)/m)

          ! Take care of printing
          if (qmc_sim%rank .eq. 0) then
             do i = 1, np
                do j = 1, np
                   tmp(1:L,1) = avgSE(i,j,0:L-1)
                   tmp(1:L,2) = errSE(i,j,0:L-1)
                   write(title,'(A,i3)') trim(pname(gflist(h)))//" SelfEn k=", k
                   write(title,'(A,i3,A,i3)') trim(adjustl(title))//'   pair=',i,',',j
                   call DQMC_Print_Array(0, L , title, label, tmp(:, 1:1), tmp(:, 2:2), OPT)
                   write(OPT,'(1x)')
                enddo
             enddo
          endif

       enddo

    enddo

    contains

       subroutine convert_to_iwn(tdmgtau, tdmgw)
          complex(wp), intent(in)  :: tdmgtau(npp, 0:L-1)
          complex(wp), intent(out) :: tdmgw(np, np, 0:L-1)
          ! Local variables
          complex(wp) :: valuetl(0:L-1), valuewl(0:L-1)
          integer     :: ipl, jpl, ijpl
          complex(wp), parameter :: unum=(1.0_wp,0.0_wp), nil=(0.0_wp,0.0_wp)
          ijpl = 0
          do ipl = 1, np
             do jpl = ipl, np
                ijpl = ijpl + 1
                valuetl = tdmgtau(ijpl, 0:L-1)
                if (ipl .eq. jpl) valuetl(0) = valuetl(0) - 0.5_wp
                call zgemv('N', L, L, unum, ftw, L, valuetl, 1, nil, valuewl, 1)
                tdmgw(ipl, jpl, 0:L-1) = valuewl
                if (ipl .ne. jpl) then
                   valuetl = conjg(tdmgtau(ijpl, 0:L-1))
                   call zgemv('N', L, L, unum, ftw, L, valuetl, 1, nil, valuewl, 1)
                   tdmgw(jpl, ipl, 0:L-1) = valuewl
                endif
             enddo
          enddo
       end subroutine convert_to_iwn

     
       subroutine invertG(tdmgw)
          complex(wp), target, intent(inout) :: tdmgw(np, np, 0:L-1)
          ! Local variables
          complex(wp) :: work(np)
          complex(wp), pointer :: gw(:,:)
          integer     :: ipiv(np)
          integer     :: iwl, info
          do iwl = 0, L-1
             ! Fill matrix. Note that G is complex symmetric. Not hermitian.
             gw => tdmgw(1:np, 1:np, iwl)
             call zgetrf(np, np, gw, np, ipiv, info)
             call zgetri(np, gw, np, ipiv, work, np, info)
          enddo
       end subroutine invertG


  end subroutine DQMC_TDM1_SelfEnergy

  !--------------------------------------------------------------------!

end module DQMC_TDM1
