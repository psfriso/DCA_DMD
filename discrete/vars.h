!
! CommLine & I/O
!
   integer, parameter :: NFILES = 5
 
   integer unit_i, unit_o,unit_ener, unit_traj, unit_dca
 
   type(commLineOption) :: files(NFILES) = (/&
   commLineOption("-i",    "param",          "formatted",   "old",     "Settings"),&
   commLineOption("-pdbin", "raw PDB file",  "formatted",   "old",     "Initial PDB"),&
 !  commLineOption("-ener", "energy",         "formatted",   "unknown", "Energies"),&
   commLineOption("-pdbtarg",  "target.pdb", "formatted",   "old",     "Target PDB"),&
   commLineOption("-dca", "dca_results.dat", "formatted",   "old",     "full List Coev Pairs"),&
   commLineOption("-o",    "log",            "formatted",   "unknown", "Calculation Log")/)

  INTEGER, PARAMETER :: DBL=SELECTED_REAL_KIND(p=13)! Double precision real number definition processor independent
  INTEGER, PARAMETER :: SGL=SELECTED_REAL_KIND(p=6) ! Single Precision real number definition processor independent
 !
 real(dp):: tacum,tacrect,taccorr,tacact      
 ! Coordinates and velocities 
   type(point), allocatable :: rsp(:)! Single precision version of r for binary I/O
   type(pointDP), allocatable :: rcoord(:)!
   type(point) :: rcm ! C of M  
         
! Structure
   character(len=4), allocatable :: res(:), atp(:)
   character(len=4), allocatable :: atoms(:)
   character(len=1), allocatable :: chain(:)
   integer, allocatable :: rnum(:)
   integer natom
   !
   type(simulRT) :: refRp 
   type(simulRT) :: rp 
   !        
   ! Pedro
   REAL*8 :: error
   real*8, dimension(3,3) :: U               ! Rotation Matrix
   real*8, dimension(3) :: center1, center2  ! Center of initial and target structures
   !
   
 ! Potentials & energies
   real xmassa, calcEkin
   
! Interaction pairs

! Collisions
   integer ibloc
! Time
   real*8 tinit, tsetup, tfin
!
   integer i,j, ioerr
!
! Old setup
   type(struc) ::  str
!
! PEDRO
   real(dp), allocatable :: dist_INI(:,:)
!   character(len=50) :: fmtsteppot='(2I4,I2,x,I1,4(2(f7.1,X),L1,X),L1,X,f9.4,f8.3)'
   INTEGER :: nwritten=0
   integer  ,allocatable :: semilla(:)     
   integer :: cont
   integer :: npairs
   integer :: k
   type(resultsTable), allocatable, dimension(:) :: summary
   real :: initialener=0._SGL
   integer :: todos
   integer, allocatable, dimension(:,:) :: typeint
   integer :: positives=0
   integer :: negatives=0
   integer , allocatable, dimension(:,:) :: cv 
   type( dcaout ),dimension(:), allocatable :: listaDca  
   integer :: allpairs
integer :: nL=0
real :: fence
type(newmodels), allocatable,dimension(:) :: genModels
integer :: acceptedModels
