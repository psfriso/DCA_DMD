
!===============================================================================   
 subroutine readCommLine (files, unit_i, unit_o, NFILES)
 USE commLine
  integer, intent(IN) :: NFILES
  integer, intent(OUT) :: unit_i, unit_o
  type(commLineOption), intent(INOUT) :: files(NFILES) 
!
  call inputArgs(files)
  unit_i = openFn(files, '-i')
  if (fileName(files,'-o').ne.'log') then
        unit_o = openFn(files, '-o')
  else
        unit_o = 6
  endif
 end subroutine readCommLine
!===============================================================================   
 subroutine programHeader(unit_o)
   integer unit_o
#ifdef OMP
   integer omp_get_max_threads
#endif
   write (unit_o, *) "================================================="
   write (unit_o, *) "=                                               ="
   write (unit_o, *) "=            Move Evolution                     ="
   write (unit_o, *) "=                                               ="
   write (unit_o, *) "= P. Sfriso A. Emperador, JL. Gelpi, M.Orozco   =" 
   write (unit_o, *) "=                                               ="
   write (unit_o, *) "=                  (c) 2014                     =" 
   write (unit_o, *) "================================================="
   write (unit_o, *)
#ifdef OMP
   write (unit_o, '(" Running OpenMP version on ", i3, " processors")') omp_get_max_threads()
   write (unit_o, *)
#endif
 end subroutine programHeader
!===============================================================================   
 subroutine errorAllocmem (ioerr, text)
  integer ioerr
  character(len=*) text
  if (ioerr.ne.0) then
    write (0, '("Error allocating memory",a30)') text
    stop 1
  endif
 end subroutine errorAllocmem
!===============================================================================   
 subroutine writeSnapshot (unit_traj, ibloc, r, atom, res, chain, rnum, natom)
 use geometryDP
  integer, intent(IN) :: natom, unit_traj, ibloc
  type(pointDP), intent(IN) :: r(natom)
  character(len=4), intent(IN) :: atom(natom), res(natom)
  character(len=1), intent(IN) :: chain(natom)
  integer, intent(IN) :: rnum(natom)
  integer i, j
     write (unit_traj, '("MODEL",8X,I4)') ibloc
     do i = 1,natom
      j=i
      write (unit_traj, '("ATOM",2X,I5,2X,A3,1X,A3,1X,A1,1X,I3,4X,3F8.3)') i, atom(j), res(j), chain(j), rnum(j), r(i)
     enddo
     write (unit_traj, '("ENDMDL")')
 end subroutine writeSnapshot
!===============================================================================   

!!===============================================================================   
!! call writeEnergies(unit_ener, unit_o, temps, epotgo, ekin0, natom)
!subroutine writeEnergies(unit_ener, unit_o, temps, epotgo, ekin0, natom,dist,dist2,error)
! !   
!    use constants
!    integer, intent(IN) :: unit_ener, unit_o, natom
!    real, intent(IN) ::  epotgo, ekin0,dist,dist2
!    real*8, intent(IN) :: temps
!    real, intent(in) ::error
!    !   
!   ! write (unit_ener, *)temps * 1.e-3, epotgo, ekin0/1.5/natom/8.314,dist
!    write (unit_o, *) real(temps * 1.e-3), epotgo, dist,dist2,error
!    
!end subroutine writeEnergies
!!===============================================================================   
 SUBROUTINE writePairList(unit_o,natom,contarBondedInt,stepPts,inicont)
      use stepPotentials
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: unit_o
      integer, intent(in) :: inicont
      INTEGER, INTENT(IN) :: natom
      INTEGER, INTENT(IN) :: contarBondedInt
      TYPE( stepPotInt), INTENT (IN) :: stepPts(natom,natom)
      !
    write (unit_o, '(" Initial Pairs list")')
    write (unit_o, '(" ==================")')
    write (unit_o, '(" Total:              ",i9)') natom * (natom - 1) / 2
    write (unit_o, '(" Bonded:             ",i9)') contarBondedInt
    write (unit_o, '(" Initial contacts:   ",i9)') inicont
    write (unit_o, '(" Non Bonded:         ",i9)') count(stepPts % active)
    write (unit_o, '("   GO 1 SQ WELL      ",i9)') count(stepPts % nstep .eq. 2 .and.stepPts%active)
    !
    write (unit_o, *)
    write (unit_o, '(" System setup completed")')
    write (unit_o, *)
  END SUBROUTINE  writePairList
 !=================================================================================
subroutine writeOut( unit_o, temps, epotgo, dist,dist2,error)
 !   
    integer, intent(IN) :: unit_o
    real, intent(IN) ::  epotgo,dist,dist2
    real*8, intent(IN) :: temps
    real, intent(in) ::error
  !
    write (unit_o, '(f10.3,X,f15.3,x,3(f10.3,x))') real(temps * 1.e-3), epotgo, dist,dist2,error
    
end subroutine writeOut
