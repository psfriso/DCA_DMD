 MODULE paramSet
 use intlist
 use nrtypes
! Input param
!
   real, save :: &
      temp = 300., &
      tmin = 1.e-22*1.e15
  integer, save :: &
      nbloc = 1000, &
      seed = 2381
  real(dp), save :: &
      tsnap = 1000., &
      tcorr = 100., &
      tact = 15., &
      trect = 100.
      
 real :: goener=0.300    ! 0.5
  real :: sigmago=0.1500   ! 0.15 
  real :: sigmabond=0.0200
  real :: sigmaPbond=0.0500 !
  real(dp) :: rcut1= 200.00         ! 200
INTEGER :: nskip=5
real(dp)::rcontacts=50.0d0
real :: mindistance=5.5000000 !5.5
real :: dReq=0.99950000 !0.999
integer :: inequilibrium=40 !50
real :: numericzero=0.0000001
real :: filterbydistance = 6.50  !6.5
integer :: minaadist=5
real :: distanceInitialState=10.0  !9.5
real,PARAMETER ::    rvdw = 2.0 
real,PARAMETER ::    rhc  = 2.0
real,PARAMETER ::    xm   = 0.012
real ::    xsum = 1.0/xm !0.5 * (1./xm + 1./xm)
real ::    distHC=(rhc + rhc)**2
real :: maxDistanceTocontact= 6.5 ! 6.5+5.5 =12 A
real :: repulsiveMin=4.0
real :: repulsiveMax=6.0
LOGICAL :: globalMin=.TRUE.
real :: largerdist=1.10_sp  !1.10
CONTAINS
!


!===============================================================================
 subroutine readInputParamSet (unit_i)
   integer, intent(IN) :: unit_i
!   
   namelist /input/ seed,rcut1,nskip,goener,sigmago,dReq,repulsiveMax,mindistance,globalMin
   !
   read (unit_i, INPUT)
   ! checking 
   if (TRECT.gt.TSNAP) TRECT = TSNAP
   if (TCORR.gt.TRECT) TCORR = TRECT
   if (TACT.gt.TCORR)  TACT = TCORR

 end subroutine readInputParamSet
!===============================================================================
 subroutine writeInputParamSet (unit_o)
   integer, intent(IN) :: unit_o
   ! Copiem l'arxiu de parametres. Pendent format
   write (unit_o, *)
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, '(" | CALCULATION PARAMETERS                                   |")')
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, '(" | Simulation settings                                      |")')
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, '(" | Simulation Time (ps) (Nbloc x TSnap)       |",f12.3," |")') NBLOC * TSNAP / 1.e3
   write (unit_o, '(" | Output structure (fs)             | TSnap  |",f12.3," |")') TSNAP 
   write (unit_o, '(" | Re-scoring target (fs)            | Trect  |",f12.3," |")') TRECT
   write (unit_o, '(" | Update velocities (fs)            | Tcorr  |",f12.3," |")') TCORR
   write (unit_o, '(" | Update Lists, collision times (fs)| Tact   |",f12.3," |")') TACT
   write (unit_o, '(" | Min. accepted colision time (fs)  | TMin   |",f12.8," |")') TMIN   
   write (unit_o, '(" | Temperature (K)                   | Temp   |",6X,f6.2, " |")') TEMP
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, '(" | Other                                                    |")')  
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, '(" | Random generator seed                      |",7X,i5  " |")') seed
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, '(" | Strcutural & Evolution                                   |")')  
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, '(" | Energy of interaction (kcal/mol)  | GOENER |",f12.3," |")') goener
   write (unit_o, '(" | Well Amplitude, non bonded        | sigmago|",f12.3," |")') sigmago
   write (unit_o, '(" | Well Amplitude,     bonded      | sigmabond|",f12.3," |")') sigmabond
   write (unit_o, '(" | Well Amplitude, pseudo-bonded  | sigmaPbond|",f12.3," |")') sigmaPbond
   write (unit_o, '(" | Radius^2 of interaction non bond  | rcut1  |",f12.3," |")') rcut1
   write (unit_o, '(" | Number of skipping frames         | nskip  |",I12," |")') nskip
   write (unit_o, '(" | Radius^2 of interaction contact | rcontacts|",f12.3," |")') rcontacts
   write (unit_o, '(" | Minimum distance force contact | midistance|",f12.3," |")') mindistance
   write (unit_o, '(" | Factor reducing max distance      |  dReq  |",f12.8," |")') dReq
   write (unit_o, '(" | Minimum distance force contact | midistance|",f12.3," |")') mindistance
   write (unit_o, '(" | # frames in equilibrium     | inequilibrium|",I12," |")') inequilibrium
   write (unit_o, '(" | Numeric Zero                  | numericzero|",f12.3," |")') numericzero
   write (unit_o, '(" | distance min to simulate | filterbydistance|",f12.3," |")') filterbydistance
   write (unit_o, '(" | distance in seq to simulate    | minaadist |",I12," |")') minaadist
   write (unit_o, '(" | Initial State  dist(A)|distanceInitialState|",f12.3," |")') distanceInitialState
   write (unit_o, '(" | Minimum distance force contact | midistance|",f12.3," |")') mindistance
   write (unit_o, '(" | van deer Waals Radius          | rvdw      |",f12.3," |")') rvdw
   write (unit_o, '(" | Radius of Hardcore             | rhc       |",f12.3," |")') rhc
   write (unit_o, '(" | Particles mass                 | xm        |",f12.3," |")') xm
   write (unit_o, '(" | Max dist allowed to OK|MaxDistanceToContact|",f12.3," |")') MaxDistanceToContact
   write (unit_o, '(" | NonContact Well MIN         | repulsiveMin |",f12.3," |")') repulsiveMin
   write (unit_o, '(" | NonContact Well MAX         | repulsiveMAX |",f12.3," |")') repulsiveMax
   write (unit_o, '(" | How min distance is computed  |GlobalMin   |",L12," |")') globalmin
   write (unit_o, '(" | Factor enlarger initial dist  | largerdist |",f12.3," |")') largerdist
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, '(" -------PS 27 MAR 2014 --------------------------------------")')

      
 end subroutine writeInputParamSet
!===============================================================================
 END MODULE paramSet