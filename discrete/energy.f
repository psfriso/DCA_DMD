 
! ===============================================
 SUBROUTINE activateClist(natom, r, distHC, rcontacts,clist)
        USE geometryDP
        USE intList
        !
        IMPLICIT NONE
        !
        TYPE(overlapscheck), DIMENSION(natom), INTENT(INOUT) :: clist
        INTEGER, INTENT(IN) :: natom
        TYPE (pointDP) , DIMENSION(natom), INTENT(IN) :: r
        REAL*8, INTENT(IN) :: rcontacts
        REAL, INTENT(IN) :: distHC
        !
        integer :: i, j
        !
        clist%nats = 0
        do i=1,natom
            clist(i)%vecinos=closeOnes(0,0,0.0)
        end do
        !
        ! vamos a generar la lista de vecinos para sacar posibles overlaps
        !
        do j = 3, natom
            do i = 1, j - 2 ! el j-1 ya lo tiene en cuenta el bonded term
                IF(calcDist2DP(r(i),r(j)) < rcontacts ) THEN
                    clist(i) % nats = clist(i) % nats + 1
                    clist(j) % nats = clist(j) % nats + 1
                    !
                    clist(i) % vecinos(clist(i) % nats ) = closeOnes(j,clist(i) % nats,distHC)
                    clist(j) % vecinos(clist(j) % nats ) = closeOnes(i,clist(j) % nats,distHC)
                end if
            end do
        end do
    END SUBROUTINE activateClist

!===============================================================================
    !===============================================================================
   subroutine activateStepPotIter( rpl )
   use stepPotentials
   use geometryDP
   use intList
   use paramSet
   use simulation
   use nrtypes
   !
   implicit none
   !
   type(simulRT), intent(inout) :: rpl
   !
   !integer, intent(IN) :: natom
   !TYPE (pointDP), INTENT(IN), Dimension(natom) :: r
   !type(stepPotInt), intent(INOUT) :: stepPts(natom,natom)
   !type(intpList), intent(INOUT) :: nblist(natom)
   !real, intent(IN) :: xsum
   integer i, j
   real*8 :: rij_INI
!   type(pointDP) rj
!   logical :: IsActive
!
   where (rpl%stepPts%active)
      rpl%stepPts%active = .false.
   end where
   
   rpl%nblist%nats = 0
  ! write(*,*) "here 0"

  do j = 2,rpl%natom
    !  rj = r(j) ! intentem millorar cache hits a r(i)
      do i = 1,j-1
         rij_INI = calcDist2DP(rpl%r(i),rpl%r(j))
          
           isactive1:  IF(rij_INI< rcut1) THEN
                rpl%stepPts(i,j)%active = .true.
                rpl%nblist(i)%nats = rpl%nblist(i)%nats + 1
                rpl%nblist(j)%nats = rpl%nblist(j)%nats + 1
                rpl%nblist(i)%iData(rpl%nblist(i)%nats) = intData(j, rpl%nblist(j)%nats, rpl%stepPts(i,j), xsum, 1.e15, 0.)
                rpl%nblist(j)%iData(rpl%nblist(j)%nats) = intData(i, rpl%nblist(i)%nats, rpl%stepPts(i,j), xsum, 1.e15, 0.)
             END IF isactive1
   
    
      enddo
  enddo
  
  end subroutine activateStepPotIter
!===============================================================================
!===============================================================================
   subroutine activateStepPot(stepPts,dist_INI,rcut2GO2,natom, nblist, xsum)
   use stepPotentials
   use geometryDP
   use intList
   integer, intent(IN) :: natom
   real*8, intent(IN) :: rcut2GO2
   REAL*8, INTENT (IN) :: dist_INI(natom,natom)
   type(stepPotInt), intent(INOUT) :: stepPts(natom,natom)
   type(intpList), intent(INOUT) :: nblist(natom)
   real, intent(IN) :: xsum
   integer i, j
   real*8 :: rij_INI
!   type(pointDP) rj
!   logical :: IsActive
!
   where (stepPts%active)
      stepPts%active = .false.
   end where
   
   nblist%nats = 0
  ! write(*,*) "here 0"

  do j = 2,natom
    !  rj = r(j) ! intentem millorar cache hits a r(i)
      do i = 1,j-1
         rij_INI = dist_INI(i,j)**2
          
           isactive1:  IF(rij_INI< rcut2GO2) THEN
                stepPts(i,j)%active = .true.
                nblist(i)%nats = nblist(i)%nats + 1
                nblist(j)%nats = nblist(j)%nats + 1
                nblist(i)%iData(nblist(i)%nats) = intData(j, nblist(j)%nats, stepPts(i,j), xsum, 1.e15, 0.)
                nblist(j)%iData(nblist(j)%nats) = intData(i, nblist(i)%nats, stepPts(i,j), xsum, 1.e15, 0.)
             END IF isactive1
   
    
      enddo
  enddo
  
 end subroutine activateStepPot
!===============================================================================
 
! subroutine thermalize2(seed, iev, natom, TEMP, xmassa, v, xm, ekin)
! use geometryDP
! use ran_mod
!   integer, intent(IN) :: iev, seed, natom
!   real, intent(IN) ::  xmassa, temp, xm
!   type(pointDP), intent(OUT) :: v(natom)
!   real, intent(OUT) :: ekin
!   integer i, kk,NegSeed
!   real calcEkin
!   type(pointDP) vcm
!   real*8 fi, sto,aux
!
!   aux=0.0
!   negSeed=seed
!   fi=ran1(negSeed)
  ! fi=ran1(seed)
   !   write(*,*)"iev desde thermalize", iev
!   vcm = pointDP(0., 0., 0.)
!   do i = 1,natom
!      kk = (seed + 3 * i + 1 + iev)
!!   write(*,*)"kk desde thermalize", kk
!      fi=ran1(kk)
!      v(i)%x = fi - 0.5
!!   write(*,*)"kk desde thermalize", kk
!       kk = (seed + 3 * i + 1 + iev)
!      fi=ran1(kk)
!      v(i)%y = fi - 0.5
!!   write(*,*)"kk desde thermalize", kk
!      kk = (seed + 3 * i + 1 + iev)
!      fi=ran1(kk)
!      v(i)%z = fi - 0.5
!!
!      vcm = vcm + xm*1.d0 * v(i)
!   enddo
!   vcm = (1.d0/xmassa) * vcm
! !  write(*,*)" vel CM", vcm
!   do i = 1,natom
!      v(i) = v(i) - vcm
!   enddo
!   ! sto=sqrt( 1.5 *natom*TEMP * 0.001987 / calcEkin(...)  )
!   sto = sqrt(1.5 * natom * TEMP * 8.314/ calcEkin(v, xm, natom))
!   do i = 1,natom
!      v(i) =  sto * v(i)
! !     write(*,*) v(i)

!   enddo
   
  ! call rnd_gauss ( aux, xm, temp)

  ! write(*,*)" velocidades", v(1)%x, aux
!   ekin = 1.5 * natom * TEMP ! falta Kb No cal recalcularla 
   
! end subroutine thermalize2
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!=======================================================================
 
 subroutine thermalize( natom, TEMP, v, xm)
     
   use geometryDP
   
   integer, intent(IN) :: natom
   real, intent(IN) ::  temp, xm
   type(pointDP), intent(OUT) :: v(natom)
 !  real, intent(OUT) :: ekin
!   integer i, kk,NegSeed
   real calcEkin
   type(pointDP) vcm
   real*8 fi, sto!,aux
   integer :: i
   !
   vcm = pointDP(0., 0., 0.)
   !
    do i = 1,natom
        call rnd_gauss ( fi, xm, temp)
        v(i)%x = fi
        !
        call rnd_gauss ( fi, xm, temp)
        v(i)%y = fi
        !
        call rnd_gauss ( fi, xm, temp)
        v(i)%z = fi
        !
!
      vcm = vcm + xm*1.d0 * v(i)
   enddo
!
      do i = 1,natom
        v(i) = v(i) - vcm
     enddo
 !
  sto = 0
  sto = sqrt(1.5 * natom * TEMP * 8.314/ calcEkin(v, xm, natom))
   do i = 1,natom
      v(i) =  sto * v(i)
 !     write(*,*) v(i)

   enddo    
     
 end subroutine thermalize
 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!==========================================================
 !
 
!========================================================================
 subroutine rnd_gauss ( random, xm, T)
    !
      use nrtypes
    !
    implicit none
  
    !
   ! INTEGER, PARAMETER :: DBL=SELECTED_REAL_KIND(p=13)! Double precision real number definition processor independent
  !  INTEGER, PARAMETER :: SGL=SELECTED_REAL_KIND(p=6) ! Single Precision real number definition processor independent
!
    real(sp) , INTENT(IN) :: xm,T
    real(dp) , intent ( out ) :: random
 !   real(DBL) , parameter :: R = 8.314_DBL
    real(dp) :: R
    real(sp) :: conversionSItoAFs
!     	
    real(dp) :: rnd1 , rnd2, std_dev!, pi
 !   real, parameter :: conversionSItoAFs=0.00001_SGL
    R = 8.314_dp
    conversionSItoAFs=0.00001_sp
!    
    std_dev = sqrt((T*R) / xm)
 !   pi = 4.0d0 * atan ( 1.0d0 )
    call random_number ( rnd1 )
    call random_number ( rnd2 )
    random = conversionSItoAFs*std_dev * sqrt ( -2.0d0 * log ( rnd1 ) ) * cos ( 2.0d0 * pi * rnd2 )
 end subroutine rnd_gauss
 
 
!========================================================================
   subroutine calcEpot(natom, r, stepPts, epotgo)
        use geometryDP
        use stepPotentials
        use constants

        integer, intent(IN) :: natom
        type(pointDP), intent(IN) :: r(natom)
        type(stepPotInt), intent(IN) :: stepPts(natom, natom)
        real, intent(OUT) :: epotgo
        real epot(MAXTIPINT)
        real dist
        integer i, j, k
        !PENDENT TREBALLAR AMB NBLIST
        epotgo = 0.
        epot = 0.
        do j = 2, natom
            do i = 1, j - 1
                if (stepPts(i, j) % active) then
                    dist = sqrt(real(calcDist2DP(r(i), r(j))))
                    k = stepPts(i, j) % nstep
                    do while ((k .gt. 1).and.dist .lt. stepPts(i, j) % step(k) % r)
                        epot(stepPts(i, j) % tipInt) = epot(stepPts(i, j) % tipInt) - stepPts(i, j) % step(k) % e / FACTE
                        k = k - 1
                    enddo
                    if (dist .lt. stepPts(i, j) % step(k) % r) &
                    epot(stepPts(i, j) % tipInt) = epot(stepPts(i, j) % tipInt) - stepPts(i, j) % step(k) % e / FACTE
                endif
            enddo
       enddo
       epotgo = epot(SSEC)
   end subroutine calcEpot
    !========================================================================
 pure function  calcEkin (v, xm, natom) result (ekin)
 use geometryDP
   integer, intent(IN) :: natom
   real ekin
   type(pointDP), intent(IN) :: v(natom)
   real, intent(IN) :: xm
   real*8, parameter :: a2 = 1.d-20*1.d30 ! 10^7 /4.184
   integer i
   ekin = 0.
   do i = 1,natom
      ekin = ekin + 0.5 * xm* real(A2 * dotDP(v(i), v(i)))
   enddo
 end function calcEkin
!===============================================================================

