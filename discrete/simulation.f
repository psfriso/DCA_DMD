!     
! File:   simulation.f
! Author: psfriso
!
! Created on April 1, 2014, 11:13 AM
!
Module simulation

    Use intlist
    use geometryDP
    use constants
    use nrtypes
    use steppotentials
    use evolution
    !
    type simulrt
        type(stepPotInt) , allocatable :: stepPts(:,:)
        type(pointDP), allocatable :: r(:)
        type(pointDP), allocatable :: v(:)
        real(sp) :: ekin
        real(sp) :: epotgo
        real(sp) :: ekin0
        real(sp) :: temp
        type(intpList), allocatable :: blist(:)
        type(overlapscheck), allocatable :: clist(:)
        type(intpList), allocatable :: nblist(:)
        type(intpList), allocatable :: Flist(:)
        logical, allocatable :: toUpdate(:)
        real(dp) :: temps
        real(dp) , allocatable :: tpart(:)
        integer , allocatable :: ipart(:)
        integer :: iev
        integer :: ierr
        type(ROC), allocatable :: performance(:)
        integer :: ibloc
        integer :: cont
        integer :: seed
        integer :: natom
    end type simulrt
    !
    CONTAINS
!============================================================================
    subroutine initializeSimulrt(str)
        !
        Use intlist
        use geometryDP
        use constants
        use nrtypes
        use steppotentials
  !
        implicit none
        !
        type(simulrt), intent(inout) :: str
        !
        str % stepPts =stepPotInt(0,0,(/stepPot(0.00000,0.00000),stepPot(0.00000,0.00000),&
        stepPot(0.00000,0.00000),stepPot(0.00000,0.00000)/),.FALSE.,0.0000)
        
        str % r = pointDP(0._dp,0._dp,0._dp)
        str % v = pointDP(0._dp,0._dp,0._dp)
        str % nblist % nats =  0
        !str % nblist % idata(:)  = intData(0 , 0 , str % stepPts ,0._sp,0._dp,0._sp) )
        str % clist % nats = 0
       ! str % clist % vecinos % patnum = 0
       ! str % clist % vecinos % simp = 0
       ! str % clist % vecinos % dmin = 0._sp
        str % toUpdate =.FALSE.
        str % blist % nats =  0
        !str % blist % idata(:)  =intData(0,0,str % stepPts,0._sp,0._dp,0._sp) )
        str % Flist % nats =  0
        !str % Flist % idata(:)  =intData(0,0,str % stepPts ,0._sp,0._dp,0._sp) )
            !
        str % ekin = 0._sp
        str % epotgo = 0._sp
        str % temp = 0._sp
        str % temps = 0._dp
        str % tpart = 0._dp
        str % ipart = 0
        str % iev = 0
        str % ierr = 0
        !
        str % performance = ROC(0,0,0,0,0._sp,0,0,0._sp,0._sp)
        str % ibloc = 0
        str % cont  = 0
        str % seed  = 0
    END subroutine initializeSimulrt
!============================================================================
!============================================================================
    subroutine deallocateSimulrt(natom,str)
        !
        implicit none
        !
        integer, intent(IN) :: natom
        type(simulrt), intent(inout) :: str
        !
        integer :: i, ioerr
        ioerr=0
        !
              !
       ! DO i=1,natom
       !     call deallocateintPList(4, ioerr,str % blist(i))
       !     IF (ioerr .gt. 0) call errorAllocMem(ioerr, 'Bonded List')
       !     call deallocateintPList(INT(0.6*natom), ioerr,str % nblist(i))
       !     IF (ioerr .gt. 0) call errorAllocMem(ioerr, 'Non Bonded List')
       !     call  deallocateCList(natom, ioerr,str % clist(i))
       !     IF (ioerr .gt. 0) call errorAllocMem(ioerr, 'contacts List')
       !     call deallocateintPList(1, ioerr,str % Flist(i))
       !     IF(ioerr .gt. 0) call errorAllocMem(ioerr, 'force List')
       ! END DO
        !
        deallocate( &
        str % stepPts,&
        str % r,&
        str % v, & 
        str % nblist, &
        str % clist, &
        str % toUpdate, &
        str % blist, &
        str % Flist, &
        str % ipart, &
        str % tpart, &
        str % performance, stat=ioerr)
        IF (ioerr .gt. 0) call errorAllocmem(ioerr, 'allocateSimulRT')
 
    END subroutine deallocateSimulrt
!============================================================================
!============================================================================
    function allocateSimulrt(natom,todos) result (str)
        !
        implicit none
        !
        integer, intent(IN) :: natom
        integer, intent(IN) :: todos
        type(simulrt) :: str
        !
        integer :: i, ioerr
        ioerr=0
        !
        allocate( &
        str % stepPts(natom,natom),&
        str % r(natom), &
        str % v(natom), & 
        str % nblist(natom), &
        str % clist(natom), &
        str % toUpdate(natom), &
        str % blist(natom), &
        str % Flist(natom), &
        str % ipart(natom), &
        str % tpart(natom), &
        str % performance(todos) , stat=ioerr)
        IF (ioerr .gt. 0) call errorAllocmem(ioerr, 'allocateSimulRT')
        !
        DO i=1,natom
            str % blist(i) = allocateintPList(4, ioerr)
            IF (ioerr .gt. 0) call errorAllocMem(ioerr, 'Bonded List')
            str % nblist(i) = allocateintPList(INT(natom), ioerr)
            IF (ioerr .gt. 0) call errorAllocMem(ioerr, 'Non Bonded List')
            str % clist(i) = allocateCList(INT(natom), ioerr)
            IF (ioerr .gt. 0) call errorAllocMem(ioerr, 'contacts List')
            str % Flist(i) = allocateintPList(1, ioerr)
            IF(ioerr .gt. 0) call errorAllocMem(ioerr, 'force List')
        END DO
        !
        str % ekin = 0._sp
        str % epotgo = 0._sp
        str % temp = 0._sp
        str % temps = 0._dp
        str % ipart = 0
        str % iev = 0
        str % ierr = 0
        str % seed = 0
        !
        str % performance = ROC(0,0,0,0,0._sp,0,0,0._sp,0._sp)
        str % ibloc = 0
        str % cont  = 0
    END FUNCTION allocateSimulrt
!============================================================================
    function cloneSimulrt(str ) result( newstr )
        type(simulRT), intent(IN) :: str
        type(simulRT) :: newStr
        !
        newstr=str
    end function cloneSimulrt
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
    subroutine printHelloWorld(unit_o)
        !
        implicit none
        !
        integer unit_o
        !
        write(unit_o,*) "Hello World"
    end subroutine printHelloWorld
    !
!============================================================================
    SUBROUTINE blistThem(natom,xsum,dist_ini,stepPts,blist)
       !
        use constants
        use stepPotentials
        use intlist
        use nrtypes
       !
       Implicit None
       !
       integer, intent(in) :: natom
       real(dp), intent(in), dimension(natom,natom) :: dist_ini
       type(stepPotInt), intent(inout),dimension(natom,natom) :: stepPts
       real, intent(in) :: xsum
       type(intplist), dimension(natom), intent(inout) :: blist
       !
       integer :: i,j
       integer :: contarBondedInt
       contarBondedInt=0
       !
       do j = 2,natom
         do i = 1, j - 1
             if (stepPts(i, j) % tipInt .eq. COVB .or.stepPts(i, j) % tipInt .eq. COVF  ) then
                 stepPts(i, j) % active = .false.
                 stepPts(j, i) % active = .false.
! j > i only   
                 IF(dist_INI(i,j) >15.0 ) THEN
                     stepPts(i, j) % tipInt=0
                     stepPts(j, i) % tipInt=0
                     contarBondedInt=contarBondedInt+1
                 ELSE
                 blist(i) % nats = blist(i) % nats + 1
                 blist(i) % iData(blist(i) % nats) = intData(j, blist(j) % nats, stepPts(i, j), xsum, 1.e15, 0.)
                 contarBondedInt=contarBondedInt+1
                 ENDIF
              endif
           enddo
      enddo

        
    END SUBROUTINE blistThem
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!===================================================================    
    SUBROUTINE modifyStepPts(natom,p1,p2,dist_ini,stepPts)
        !
        use steppotentials
        use constants
        use paramset
        use geometryDP
        use nrtypes
        !
        IMPLICIT NONE
        !
        integer, intent(IN) :: p1
        integer, intent(IN) :: p2
        integer, intent(IN) :: natom
        real(dp) , dimension(natom,natom) :: dist_ini
        type(stepPotInt), intent(inout),dimension(natom,natom) :: stepPts
        !
         stepPts(p1, p2)%tipint=SPRING
         stepPts(p1, p2)%dmin=2*rhc
         stepPts(p1, p2)%step(1)%r=2.05*rvdw
         stepPts(p1, p2)%step(2)%r=largerdist*dist_ini(p1,p2)
         stepPts(p2, p1) = stepPts(p1, p2)
        !  
         
    END SUBROUTINE modifyStepPts
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!=========================================================================== 
  SUBROUTINE flistThem(natom,xsum,stepPts,flist)
      !
      use constants
      use steppotentials
      use intlist
      !
      IMPLICIT NONE
      !
      integer , intent(in) :: natom
      real , intent(in) :: xsum
      type(stepPotInt), intent(inout),dimension(natom,natom) :: stepPts
      type(intplist), dimension(natom), intent(inout) :: flist
      !
      integer i,j
      !
      Flist%nats=0      
       do j = 2,natom
        do i = 1, j - 1
            if (stepPts(i, j) % tipInt .eq. SPRING ) then
                stepPts(i, j) % active = .false.
                stepPts(j, i) % active = .false.
                Flist(i) % nats =  Flist(i) % nats + 1
                Flist(i) % iData(Flist(i) % nats) = intData(j, Flist(j) % nats, stepPts(i, j), xsum, 1.e15, 0.)                
             endif
          enddo
     enddo
     !
  END SUBROUTINE flistTHem
!============================================================================
end module simulation
