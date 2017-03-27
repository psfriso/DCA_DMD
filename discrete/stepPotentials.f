
MODULE stepPotentials
    use Constants
    use intList    
    
CONTAINS
    !===========================================================================
 subroutine topology(natom,typeint)
        !
        implicit none
        !
        integer, intent(in) :: natom
        integer, intent(inout), dimension(natom,natom) :: typeint
        !
        integer i, j
        !
     !   write(*,*) "topology", natom, size(typeint)
        !
        DO j=2,natom
           DO i=1,j-1
              IF(abs(j-i)== 1) THEN
                  typeint(i,j)=COVB
              ELSE
                  IF (abs(j-i) .le. 3) THEN
                     typeint(i,j)=COVF
                  ELSE
                     typeint(i,j) = SSEC
                  END IF
              END IF
              typeint(j,i) =typeint(i,j)
          END DO
        END DO
    end subroutine topology
 !
    !===============================================================================
    function getStepPotFromDEF(rrefINI ,tipINT) result (st)
        !
        Use paramSet
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: tipINT
        type(stepPotInt) :: st
        real*8, intent(IN) :: rrefINI
        !
        
        st=stepPotInt(0,0,(/stepPot(0.00000,0.00000),stepPot(0.00000,0.00000),&
        stepPot(0.00000,0.00000),stepPot(0.00000,0.00000)/),.FALSE.,0.0000)
        
!       
         st % nstep = 2
         !
         SELECT CASE (tipINT)
             CASE (COVB)
               ! Covalent Bonds (i,i+1)
                   st % step(1) = stepPot( rrefINI * DBLE((1-sigmabond)), -GOENER * FACTE)
                   st % step(2) = stepPot( rrefINI * DBLE((1+sigmabond)), +GOENER * FACTE)
             CASE (SSEC)
               ! Non Bonded interactions
               IF (rrefINI < 15.0 ) THEN
                   ! Pozo con un solo minimo
                   st % step(1) = stepPot( rrefINI * DBLE((1-sigmago)), -GOENER * FACTE)
                   st % step(2) = stepPot( rrefINI * DBLE((1+sigmago)), +GOENER * FACTE)
               ELSE
                   st % step(1) = stepPot( repulsiveMin , +GOENER * FACTE)
                   st % step(2) = stepPot( repulsiveMax , -GOENER * FACTE)
               END IF
             CASE (COVF)
               ! Pseudo Bonds
                   st % step(1) = stepPot( rrefINI * DBLE((1-sigmaPbond)), -GOENER * FACTE)
                   st % step(2) = stepPot( rrefINI * DBLE((1+sigmaPbond)), +GOENER * FACTE)
             CASE (SPRING)
                 CONTINUE
         END SELECT
         
        st % tipInt = tipINT
        st % active = .FALSE.
        st % dmin = disthc

    end function getStepPotFromDef
  !===============================================================================

    END MODULE stepPotentials
