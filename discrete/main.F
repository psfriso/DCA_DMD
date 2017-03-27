! Basado en discrete0.26
!
!  ./discrete -pdbin 1usg.pdb -i dmdtest.in -dca lista_dca.dat
!
program discrete
    use commLine
    use geometry
    use geometryDP
    use stepPotentials
    use intList
    use paramSet
    use ls_rmsd
    use Structure
    use txtUtils
    use nr
    use evolution
    use simulation
    !
    !
    IMPLICIT NONE
    ! 
#include "vars.h"

 
  call readCommLine(files, unit_i, unit_o, NFILES)
  call readInputParamSet(unit_i)
 !
 ! SetUp
 ! 
  str = loadStructurePDB(openFn(files, '-pdbin'))
!
  call random_seed(size=i)
  allocate(semilla(i))
  semilla = (/ (seed, j=1, i) /)    ! Obviously, there are better choices!
  call random_seed(put=semilla)
  deallocate(semilla)
 !
 !
!  open(unit=99,file='traj.pdb')
 ! unit_traj=60
 !
 natom=str%natoms
 !
 allocate(typeint(natom,natom),stat = ioerr)
 call errorAllocmem(ioerr, 'NAtom') 
 !
 typeint=0
 call topology(natom, typeint)
 ! 
 call printHelloWorld(unit_o)
 ! 
 allocate (atoms(natom), rnum(natom), res(natom), atp(natom), chain(natom),stat = ioerr)
 call errorAllocmem(ioerr, 'NAtom') 
 !
 atoms=str%ats(:) %atomId
 rnum= str %ats(:) %resnum
 res= str %ats(:) %resIdByAtom
 chain= str %ats(:) %chainId
 atp= str %ats(:) %atType
 !
 ! Coordinates
 allocate (rsp(natom), rcoord(natom), stat = ioerr)
 call errorAllocmem(ioerr, 'Coordinates')
 !
 rsp=str % ats(:) % r
 !
    do i = 1, natom
        rcoord(i) = SPtoDP(rsp(i))
    enddo
    !
      rcm = calcCM(natom, rcoord, xm)
    !
    do i = 1, natom
        rcoord(i) = rcoord(i) - SPtoDP(rcm)
    enddo

    !
     allocate(dist_INI(natom,natom))
     dist_INI=0.0
!
    DO j=1,natom
        DO i=1,natom
            dist_INI(i,j)=calcDistDP(rcoord(i),rcoord(j))
        END DO
    END DO
!
 deallocate(rsp,atp)   
 !
 !END SETUP
 !
    call cpu_time(tinit)
    call programHeader(unit_o)
    call printFileSummary(files, unit_o)
    !
    call writeInputParamSet(unit_o)
    !
    ! Todos los atomos son CA es una asignacion a una lista entera
    xmassa = xm*natom !sum(xm(1:natom))
    !
    unit_dca  = openFn(files, '-dca' )  
    !
    allpairs=cuantosDCAcuentan( unit_dca )
    write(*,*)"all pairs", allpairs
    !
    allocate( listadca(allpairs),stat = ioerr)
    call errorAllocmem(ioerr, 'NAtom') 
    !
    rewind(unit_dca)
  !
    listadca=dcaout(0, 0 , 0  , 0._sp,0,0._sp,.FALSE.)
  !
    call gen_listadca(unit_dca,allpairs,listadca)
  !  
    call indexx_sp(-listadca%dcascore,listadca%indx)
   !
   close(unit_dca)
   !
   call fill_listadca(allpairs,natom,dist_ini,listaDCA,nL)
   !
   todos= cuantosCuentan(natom,dist_ini)
   !
   npairs=count( listadca( : )%simulate )
   !
   !
   allocate( cv(2, npairs ),summary( npairs ) , stat = ioerr)
   call errorAllocmem(ioerr, 'lista de deformaciones')
   !
   cont=1
   DO i=1,allpairs
       IF ( listadca( listadca(i)%indx )%simulate ) then
           cv(1,cont ) =listadca( listadca(i)%indx )%p1
           cv(2,cont ) =listadca( listadca(i)%indx )%p2
           cont = cont +1
       END IF
   END DO
   
   allocate(genModels(npairs),stat=ioerr)
   call errorAllocMem(ioerr, 'generated models')
  !
   do i = 1, npairs
       genModels(i) = allocateNewmodels(natom, ioerr)
       call errorAllocMem(ioerr, 'newModels')
   enddo
  !
  do i=1,npairs
      genModels(i)%newones=i
      genModels(i)%indx=-1
      DO j=1,natom
          genModels(i)%rgen(j)=pointDP(0._dp,0._dp,0._dp)
      END DO
      !
  END DO
    
    refRp=allocateSimulrt(natom,todos)
    rp=allocateSimulrt(natom,todos)
    call initializeSimulrt(refrp)
    call initializeSimulrt(rp)
  
    ! Estructura de la proteina
    do j = 2, natom
       do i = 1, j-1
             refRp%stepPts(i, j) = getStepPotFromDEF(dist_INI(i,j),typeint(i,j))
             refRp%stepPts(j,i) = refRp%stepPts(i,j)
       enddo
   enddo
   call fill_performance(natom,todos,allpairs,dist_ini,listadca,refRp%performance)
   !
   DO i = 1, todos
       refRp%performance(i)%distance=REAL(dist_INI(refRp%performance(i)%p1,refRp%performance(i)%p2))
   END DO
   !
   call blistThem(natom,xsum,dist_ini,refRp%stepPts,refRp%blist)
  !
   refRp%r=pointDP(0._dp,0._dp,0._dp)
   DO i=1,natom
       refRp%r(i) = rcoord(i)
   END DO
  !
   call activateClist(natom, refRp%r, distHC, rcontacts,refRp%clist)
   call activateStepPot(refRp%stepPts,dist_INI ,rcut1,natom, refRp%nblist, xsum)
   !
  call cpu_time(tsetup)
  !
  call calcEpot(natom, refRp%r, refRp%stepPts, refRp%epotgo )
  !
  summary=resultsTable(0,0,0._SGL,0._SGL,0._SGL,0._dp,0,.False.)
  open (80, file='report.dat')
  !
  refRp % temp = TEMP
  refRp % seed = SEED
  refRp % natom = natom
  ! 
  call thermalize( natom, refRp%temp, refRp%v, xm)
  refRp % ekin0 = refRp % ekin
  !
   positives=sum(refRp%performance(:)%isDCA)
   negatives=todos-positives
   !
   !    A escribir todo...
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, '(" | Simulation Parameters                                    |")')
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, '(" | Number of Particles                  | natom   |",I9,"|")') natom
   write (unit_o, '(" | Number of DCA Pairs read             | allpairs|",I9,"|")') allpairs
   write (unit_o, '(" | Number of DCA Pairs matching initial |         |",I9,"|")') sum(listadca(:)%isinitial) 
   write (unit_o, '(" | Number of DCA pairs of maximum MCC   | nL      |",I9,"|")') nL
   write (unit_o, '(" | Number of DCA Pairs to simulate      | npairs  |",I9,"|")') npairs
   write (unit_o, '(" | Number of all Pairs                  | todos   |",I9,"|")') todos
   write (unit_o, '(" | Number of true positives             |positives|",I9,"|")') positives
   write (unit_o, '(" | Number of true negatives             |negatives|",I9,"|")') negatives
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, *)
   !
   !
   cont=0
   cont=count(refRp%stepPts%tipInt ==COVB ) +count(refRp%stepPts%tipInt ==COVF )
   call writePairList(unit_o,natom,cont/2,refRp%stepPts,sum(refRp%clist(:)%nats))
   !
   deallocate( listadca, stat = ioerr )
   call errorAllocmem(ioerr, 'lista dca 2')
   !
   refRp % ekin0 = calcEkin(refRp%v, xm, natom)
   !
  !-----------------------------------------------------------------------------
  npairsloop:  DO k=1,npairs
  !
 
  initialener = refRp % epotgo
  call initializeSimulrt( rp )
  rp = cloneSimulrt(refrp)
  call modifyStepPts(natom, cv(1,k), cv(2,k), dist_ini, rp%stepPts)
  call flistThem(natom,xsum,rp%stepPts,rp%flist)
  rp%cont=0
  rp%temps=0._dp
 !   
    write(unit_o,*) " Contact between ", cv(1,k),cv(2,k), " pair ", k, "out of", npairs 
  !----------------------------------------------------------------------------
    tsnaploop: do
        rp %ibloc=rp%ibloc+1
        tacum = 0.
       
        do while (tacum .lt. TSNAP)
            call activateClist(natom, rp%r, distHC, rcontacts,rp%clist)
            call activateStepPotIter( rp )
            call actualizeForces(natom,cv(1,k),cv(2,k),rp%stepPts,rp%flist)
            tacrect = 0._dp
  !----------------------------------------------------------------------------
            do while (tacrect .lt. TRECT)
                taccorr = 0.
                rp % iev = 0
                rp % ierr = 0
  !----------------------------------------------------------------------------
                do while (taccorr .lt. TCORR)
                    call springForce( rp )
                    call colisioBond( rp )
                    call colisioNonBond( rp )
                    call DMDIntLoop(rp , tacact )
                    taccorr = taccorr + tacact
                    rp % ekin0 = calcEkin(rp%v, xm, natom)
             !       write(unit_o,*)"TEMP ",rp(k)%ekin0/1.5/natom/8.314 
                 !   write(*,*)rp(k) % ekin0 ,omp_get_thread_num()
                END DO
  !----------------------------------------------------------------------------
              tacrect = tacrect + taccorr
  !----------------------------------------------------------------------------
            END DO
            tacum= tacum+ tacrect
        END DO
  !---------------------------------------------------------------------------- 
        IF (abs(  rp %stepPts(cv(1,k), cv(2,k))% step(2)%r -mindistance)< numericzero ) THEN
             rp %cont=rp %cont+1
        END iF
 !
        center1=0.0_dp
        center2=0.0_dp
        U=0.0_dp
        call rmsd(natom, rp%r,rcoord, 1, U, center2, center1, error,.FALSE.)
        call pdbsuperpos(natom, U , center1, center2,rp%r)
        !
        IF ( rp%ibloc > NBLOC .or. rp%cont > inequilibrium ) EXIT tsnaploop      
      !  IF (globalMin)  call findDistMin(todos,natom,rp%r,rp%performance )
        call findDistMin(todos,natom,rp%r,rp%performance )
        call thermalize( natom, rp%temp, rp%v, xm)
        call calcEpot(natom, rp%r, rp%stepPts, rp%epotgo)
       ! call writeout(unit_o,rp(k)%temps,rp(k)%epotgo,&
       !                real(rp(k)%stepPts(cv(1,k), cv(2,k))% step(2)%r),& 
       !      real(calcDistDP(rp(k)%r(cv(1,k)),rp(k)%r(cv(2,k)))),real(error))
        ! call writeSnapShot(99, ibloc,rp(k)%r , atoms, res, chain, rnum, natom)
             rp % ekin0 = calcEkin(rp%v, xm, natom)
        !     write(unit_o,*)"out TEMP ",rp(k)%ekin0/1.5/natom/8.314 
       enddo tsnaploop
    !----------------TSNAP Loop-----------------------------------------------------
      !  IF (.not.globalMin)  call findDistMinFinal(todos,natom,rp%r,rp%performance )
        call findDistMinFinalMix(todos,natom,rp%r,rp%performance )
        !
        summary(k)%drms=error
        summary(k)%denergy=rp%epotgo-initialener
        summary(k)%p1=cv(1,k)
        summary(k)%p2=cv(2,k)
        summary(k)%dist=real(calcDistDP(rp%r(cv(1,k)),rp%r(cv(2,k))))-&
                              real(rp%stepPts(cv(1,k), cv(2,k))% step(2)%r)
        summary(k)%area=Real(AUC(todos,natom,positives,negatives, rp%performance))
        write(80,'(I4,2X,I4,2X,2X,3(f10.3,X),f8.5,I5,2X,L1)') summary(k)
        genmodels(k)%rgen=rp%r
   !
   END DO npairsloop
   !------------------------ npairs loop----------------------------------------
   close(80)
    !     
    fence=0._sp
    fence=calcFence(npairs,summary)
  !
    genmodels%indx=summary%indx
  !
    call writeitdown(npairs,natom,unit_o,summary,fence,atoms,res,chain,rnum,genModels)
  !
    call cpu_time(tfin)
    write (unit_o, *)
    write (unit_o, '(" T I M I N G S")')
    write (unit_o, '(" ============= ")')
    write (unit_o, '(" Setup: " f10.2," s")') tsetup - tinit
    write (unit_o, '(" Traj:  ",f10.2," min ")') (tfin - tsetup)/60.
    write (unit_o, *)
    
end program discrete
!============================================================================= 
