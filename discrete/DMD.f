!!
!! File:   DMD.f
!! Author: gelpi
!!
!! Created on 14 de marzo de 2012, 16:54
!!
subroutine dmdIntloop(rpl,tacact)
use intList
use geometryDP
use paramSet
use simulation
use nrtypes
!
implicit none
   real(dp), intent(INOUT) :: tacact
   type(simulRT), intent(inout) :: rpl
!
   real(dp)  ::tevent, tevent1
   integer  ::mem1, mem2, npair1, npair2, i,j
!
   tacact = 0._dp
   tevent=0._dp
   tevent1=0._dp
   rpl%toUpdate = .true.
! evolucio temporal
!------------------------------------------------------------------------------
   do while (tacact .lt. TACT)
! busca quina es la propera colisio
       do i = 1, rpl%natom
            if (rpl % toUpdate(i)) then
                rpl % tpart(i) = 1.d15               
                do j = 1, rpl % nblist(i)%nats
                    if (rpl % nblist(i)%iData(j)%timp.lt.rpl % tpart(i)) then
                        rpl%tpart(i) = rpl%nblist(i)%iData(j)%timp
                        rpl%ipart(i) = j ! ipart recull index NO num d'atom
                    endif
                    rpl%toUpdate(i) = .false.
                enddo  
            endif
        enddo
        tevent = 1.e15
        call nextCol(mem1, mem2, npair1, npair2, tevent, rpl%ipart, rpl%tpart, rpl%natom, rpl%nblist)
        !
        tevent1 = tevent - rpl%temps
        rpl%temps = tevent
        tacact = tacact + tevent1
        rpl%iev = rpl%iev + 1
        ! translacio i variacio dels temps
        do i = 1, rpl%natom ! mantenim la versio inline degut a la barreja de tipus de real !!
            rpl%r(i)%x = rpl%r(i)%x + rpl%v(i)%x * tevent1
            rpl%r(i)%y = rpl%r(i)%y + rpl%v(i)%y * tevent1
            rpl%r(i)%z = rpl%r(i)%z + rpl%v(i)%z * tevent1
        enddo
        call updateV(rpl%r, rpl%v, rpl%nblist(mem1)%iData(npair1)%deltak, xm, xsum, mem1, mem2, rpl%natom)
!        ! Assegurem que trespassem la barrera en la direcciÃ³ correcta
        rpl%r(mem1)%x = rpl%r(mem1)%x + rpl%v(mem1)%x * tevent1/100.
        rpl%r(mem1)%y = rpl%r(mem1)%y + rpl%v(mem1)%y * tevent1/100.
        rpl%r(mem1)%z = rpl%r(mem1)%z + rpl%v(mem1)%z * tevent1/100.
        rpl%r(mem2)%x = rpl%r(mem2)%x + rpl%v(mem2)%x * tevent1/100.
        rpl%r(mem2)%y = rpl%r(mem2)%y + rpl%v(mem2)%y * tevent1/100.
        rpl%r(mem2)%z = rpl%r(mem2)%z + rpl%v(mem2)%z * tevent1/100.
        !
        ! ara actualitza els temps de colisio per a les dues particules que han xocat
        ! anulem la colisio per a que no es repeteixi
        rpl%nblist(mem1)%iData(npair1)%timp = 1.d15
        rpl%nblist(mem2)%iData(npair2)%timp = 1.d15
        call updateXocPart(mem1, rpl%nblist, rpl%temps, rpl%r, rpl%v, TMIN, rpl%natom, rpl%toUpdate)
        call updateXocPart(mem2, rpl%nblist, rpl%temps, rpl%r, rpl%v, TMIN, rpl%natom, rpl%toUpdate)
    enddo
!    deallocate (tpart, ipart)
    ! end do while (tacact.lt.tact)------------------------------------------------
end subroutine dmdIntloop

!
! File:   DMD.f
! Author: gelpi
!
! Created on 14 de marzo de 2012, 16:54
!
!subroutine dmdIntloop(srt, temps, taccorr )
!    use intList
!    use geometry
!    use Simulation
!    use constants
!    use paramSet
!!
!    type(simulRT), intent(INOUT) :: srt
!    real, intent(inout) ::temps,taccorr
!!    real*8, intent(IN) :: TACT
!!    real, intent(IN) :: TMIN
!    !
!   ! real, intent (IN) ::
!    real*8 tevent, tevent1, tevent100, tacact
!    integer mem1, mem2, npair1, npair2, i
!    !
!    tacact = 0.
!    srt % toUpdate = .true.
!    ! evolucio temporal
!    !------------------------------------------------------------------------------
!    do while (tacact .lt. TACT)
!        ! update colision lists
!        do i = 1, srt % natom
!            if (srt % toUpdate(i)) then
!                srt % ipart(i) = minloc(srt % nblist(i) % iData(1:srt % nblist(i) % nats) % timp, 1)
!!Patch to avoid error when no minimum is found, all partners below timp?, to check
!                if (srt%ipart(i).gt.0) srt % tpart(i) = srt % nblist(i) % iData(srt % ipart(i)) % timp
!            endif
!        enddo
!        srt % toUpdate = .false.
!        ! seeks for the mext collision
!        call nextCol(mem1, mem2, npair1, npair2, tevent, srt % ipart, srt % tpart, srt % natom, srt % nblist)
!        !
!        tevent1 = tevent -  temps
!        tevent100 = tevent1/100.
!        temps = tevent
!        tacact = tacact + tevent1
!        srt % iev = srt % iev + 1
!        ! integracio v i r
!        do i = 1, srt % natom
!           srt % r(i) = srt % r(i) + tevent1 * srt % v(i)
!        enddo
!        !call updateV(srt % r, srt % v, srt % nblist(mem1) % iData(npair1) % deltak, &
!        !     xm,xsum, mem1, mem2, srt % natom)
!        !     
!        call updateV(srt %r, srt %v, srt %nblist(mem1)%iData(npair1)%deltak, xm, srt %nblist(mem1)%iData(npair1)%xsum, &
!        mem1, mem2, srt %natom)     
!             
!        ! Assegurem que trespassem la barrera en la direcció correcta
!        srt % r(mem1) % x = srt % r(mem1) % x + srt % v(mem1) % x * tevent100
!        srt % r(mem1) % y = srt % r(mem1) % y + srt % v(mem1) % y * tevent100
!        srt % r(mem1) % z = srt % r(mem1) % z + srt % v(mem1) % z * tevent100
!        srt % r(mem2) % x = srt % r(mem2) % x + srt % v(mem2) % x * tevent100
!        srt % r(mem2) % y = srt % r(mem2) % y + srt % v(mem2) % y * tevent100
!        srt % r(mem2) % z = srt % r(mem2) % z + srt % v(mem2) % z * tevent100
!        !
!        ! ara actualitza els temps de colisio per a les dues particules que han xocat
!        ! anulem la colisio per a que no es repeteixi
!        srt % nblist(mem1) % iData(npair1) % timp = 0.00!TMIN*0.01 
!        srt % nblist(mem2) % iData(npair2) % timp = 0.00!TMIN*0.01
!        call updateXocPart(mem1, srt % nblist, temps, srt % r, srt % v, TMIN, srt % natom, srt % toUpdate)
!        call updateXocPart(mem2, srt % nblist, temps, srt % r, srt % v, TMIN, srt % natom, srt % toUpdate)
!    enddo
!    ! end do while (tacact.lt.tact)------------------------------------------------
!    taccorr = taccorr + tacact                        
!end subroutine dmdIntloop
