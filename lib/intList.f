 MODULE intList
 use Constants
 use geometryDP
   
  type stepPot
        real r
        real e
    end type stepPot

    type stepPotInt
        integer nstep, tipInt
        type (stepPot) step(MAXSTEPS)
        logical active
        real dmin
    end type stepPotInt

    type stepPotIntDef
        character(len = 4) :: id
        integer nstep, stepTip
        real r(MAXSTEPS), e(MAXSTEPS)
    end type stepPotIntDef
    
       type intDataTop
        integer :: id, tipInt
        real :: rref, eref
    end type intDataTop    
    
    type stepPotIntDefList
        integer npots
        type(stepPotIntDef) list(MAXSPOTDEF)
    end type stepPotIntDefList

 type intData
    integer patnum, simp ! patnum: num atom de la parella, simp: index en llista de la parella
    type(stepPotInt) :: stepPt
    real :: xsum
    real*8 :: timp
    real :: deltak    
 end type intData
 
 type intpList
    integer nats
    type(intData), pointer :: iData(:)
 end type intpList
 

  type closeones
     integer patnum, simp
     real dmin
 endtype closeones
 
 type overlapscheck
     integer nats
     type(closeOnes), pointer :: vecinos(:)
 end type overlapscheck
 
 type newmodels
     integer newones
     integer indx
     type(pointDP), pointer :: rgen(:)
 endtype newmodels
 
 CONTAINS

 function allocateintPList(natom, ioerr) result (pl)
 integer, intent(IN) :: natom
 integer, intent(OUT) :: ioerr
 type (intpList) pl
  allocate (pl%iData(natom), stat=ioerr)
  pl%nats=0
 end function allocateintPList
 !
 subroutine DEallocateintPList(natom, ioerr,pl)
 integer, intent(IN) :: natom
 integer, intent(OUT) :: ioerr
 type (intpList), intent(inout) :: pl
  if ( associated(pl%iData) )deallocate (pl%iData, stat=ioerr)
  !pl%nats=0
 end subroutine deallocateintPList
 
 !
 
 function allocateCList(natom, ioerr) result (pl)
 integer, intent(IN) :: natom
 integer, intent(OUT) :: ioerr
 type (overlapscheck) pl
  allocate (pl%vecinos(natom), stat=ioerr)
  pl%nats=0
 end function allocateCList

 subroutine deallocateCList(natom, ioerr , pl)
 integer, intent(IN) :: natom
 integer, intent(OUT) :: ioerr
 type (overlapscheck) , intent(inout) :: pl
  if (associated(pl%vecinos)) deallocate (pl%vecinos, stat=ioerr)
  !pl%nats=0
 end subroutine deallocateCList
 
 function allocateNewmodels(natom,ioerr) result (pl)
     integer, intent(in) :: natom
     integer, intent(out):: ioerr
     type(newmodels) pl
     allocate ( pl%rgen(natom),stat=ioerr )
     pl%newones=0
     pl%indx =0
 end function allocateNewmodels 
 
 
 END MODULE intList
