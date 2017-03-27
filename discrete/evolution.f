!     
! File:   evo.f
! Author: psfriso
!
! Created on March 11, 2014, 1:37 PM
!

MODULE evolution
!
    use geometryDP
    use nrtypes
!
    IMPLICIT NONE
!
    PRIVATE
!
    PUBLIC resultsTable,dcaout, ROC, area ,cuantosCuentan,AUC,findDistMin,&
           cuantosDCAcuentan,gen_listadca,fill_listadca,fill_performance,&
           calcFence,findDistMinFinal,writeitdown,actualizeForces,findDistMinFinalMix
    !
    type resultsTable
        integer p1
        integer p2
      !  character(len=20) filename
        real drms
        real denergy
        real dist
        real(sp) area
        integer indx
        logical goodenough
    end type resultsTable
    !
    type ROC
        integer p1
        integer p2
        integer indx
        integer isDCA
        real distance
        integer tp
        integer fp
        real tpr
        real fpr
    end type ROC
    !
    type dcaout
        integer p1
        integer p2
        integer indx
        real dcascore
        integer isinitial
        real(dp) mcc
        logical simulate
    endtype dcaout
    !
 !   type ilist
 !      integer coev
 !       integer cont
 !   END TYPE ilist
    !
    CONTAINS
 !=========================================================================
    
    subroutine actualizeForces(natom,p1,p2,stepPts,flist)
        !
        use paramSet
        !
        implicit none
        !
        integer , intent(in) :: natom
        integer , intent(in) :: p1
        integer , intent(in) :: p2
        type(stepPotInt), intent(inout) :: stepPts(natom,natom)
        type(intpList), intent(inout) :: flist(natom)
        !
        stepPts(p1, p2)% step(2)%r=Max(stepPts(p1, p2)% step(2)%r*dReq,mindistance)
        stepPts(p2, p1)= stepPts( p1 , p2)
        Flist(p1) % iData(Flist(p1) % nats) = intData(p2,Flist(p2) % nats, stepPts(p1, p2),xsum, 1.e15, 0.)
      end subroutine actualizeForces
 !=========================================================================
      subroutine gen_listadca(unit_dca,allpairs,listadca)
        !
        use paramSet
        !
        IMPLICIT NONE
        !
        integer, intent( IN ) :: unit_dca
        integer, intent( IN ) :: allpairs
        type( dcaout ) , dimension(allpairs), intent(INOUT) :: listadca
        !
        integer :: cont
        integer :: p1
        integer :: p2
        integer :: ioerr
        real :: aux
        real :: scoredca
        !
        cont=1
        p1=0
        p2=0
        scoredca=0
        !
        DO
            READ(unit_dca,*,IOSTAT=ioerr)p1,p2,aux,scoredca
            IF ( ioerr /= 0 )EXIT
            IF (abs(p1-p2).ge.minaadist )then
               listadca(cont)%p1 =p1
               listadca(cont)%p2 =p2
               listadca(cont)%dcascore=scoredca
               cont = cont +1
            END IF
        END DO
        !
    END SUBROUTINE gen_listadca
 !=========================================================================
        !
    integer function cuantosDCAcuentan( unit_dca )
    !
      use paramSet
    !
      Implicit None
    !
      integer, intent( IN ) :: unit_dca
    !
      integer :: p1,p2,ioerr
      real :: aux
    !  
      cuantosDCAcuentan=0
 
      DO 
          READ(unit_dca,*,IOSTAT=ioerr)p1,p2,aux,aux
          IF ( ioerr /= 0 )EXIT
          IF (abs(p1-p2).ge.minaadist )then
              cuantosDCAcuentan = cuantosDCAcuentan +1
          END IF
      END DO
      !
      END FUNCTION cuantosDCAcuentan
 !=========================================================================
 function area ( x1 , x2 , y1 , y2)
    !
    !  Por el metodo de los trapecios
    !
    Implicit none
    !
       real(dp) :: area         
    !
       real(dp), intent(in) :: x1
       real(dp), intent(in) :: x2
       real(dp), intent(in) :: y1
       real(dp), intent(in) :: y2
             !
       area=0._dp
       area = 0.5_dp * (y2 + y1) * ( x2 - x1 )   
             !
end function area
!=========================================================================
    integer function cuantosCuentan(natom,dist_ini)
    !
    use paramSet
    !
    implicit none
    !
    integer, intent(in) :: natom
    real*8, dimension(natom,natom), intent(in) :: dist_ini
    !
    integer i,j
    !
     cuantosCuentan=0
    DO j=2,natom
        DO i=1,j-1
            IF ( j -i .ge. minaadist ) THEN
                IF ( dist_INI(i,j) > filterbydistance ) THEN
                    cuantosCuentan=cuantosCuentan+1
                END IF
            END IF
        END DO 
    END DO
    !
    END FUNCTION cuantosCuentan
! ==========================================================================   
 double precision function AUC(todos,natom,positives,negatives, performance)
        !
        use geometryDP
        use nr
        !
        IMPLICIT NONE
        !
        integer, intent(in) :: todos
        integer, intent(in) :: natom
        integer, intent(in) :: positives
        integer, intent(in) :: negatives
        type(ROC), intent(INOUT), dimension(todos) :: performance
        !
        integer i
        !
     !   2) ahora a ordenar los datos por distancia
          performance%indx=0
    !      
         call indexx_sp(performance%distance,performance%indx) 
    !
    !   3) es momento de calcular los true positive respecto coevolucion
    !      para cada elemento ordenado (por distancias) de la lista anterior
    !      fijarse el acumulado de dca contacts hasta dicho i, incluido
    !
    !   4) los false positive son los no true positive (TP + FP = i ) para 
    !      un i cualquiera      
    !
    !   3) y 4) lo calculamos a la vez      
    !
    !
         DO i = 1, todos
             performance( performance(i)%indx ) % tp = sum( performance(performance(1:i)%indx )%isDCA)
             performance( performance(i)%indx ) % fp = i-performance( performance(i)%indx ) % tp
         END DO
!         
    !     
    !   5) Calcular los tpr ( tp / p ) p
    !                   fpr ( fp / n)  n 
    !   
         DO i=1,todos
             performance( performance(i)%indx ) % tpr =real(performance( performance(i)%indx ) % tp) /real(positives)
             performance( performance(i)%indx ) % fpr =real(performance( performance(i)%indx ) % fp) /real(negatives)
       !      write(45,*)performance( performance(i)%indx ) % fpr,performance( performance(i)%indx ) % tpr
         END DO
  
   ! WRITE(unit_o,'(I4,I4,I10,X,I4,f8.3,I10,I10,2f8.5)')performance(performance%indx )
   !
   !   6) Calcular el AUC (area under the curve) del par { fpr , tpr }
   !      usaremos el metodo de los trapecios. 
   !      
   !      6.1) Calcular el area entre dos puntos consecutivos.
   !      6.2) Sumar todas las areas
   !      

         
  !
        AUC=0._dp
         DO i=1,todos-1
         
             AUC= AUC + area( dble(performance(performance(i)%indx )%fpr),&
                                    dble(performance(performance(i+1)%indx )%fpr),&
                                    dble(performance(performance(i)%indx )%tpr),&
                                    dble(performance(performance(i+1)%indx )%tpr) )
                                    
            if(.FALSE.) write(*,*) i, AUC
         END DO
         !
         
    !     write(*,*) "Area dp ", areadp
    !  
     !
         
       END FUNCTION AUC
 !==================================================================
 !
      SUBROUTINE findDistMin(todos,natom,r,performance )
           !
           use geometryDP
           !
           IMPLICIT NONE
           !
           integer, intent(IN) :: todos
           integer, intent(IN) :: natom
           type(pointDP), intent(IN), dimension(natom) :: r
           type(ROC), intent(INOUT), dimension(todos) :: performance
           !
           integer :: i
           !
           
            
    !    1) primero tenemos que saber la distancia minima por residuo 
    !        entre la esctructura inicial y final. Es decir si fue en 
    !        algun momento un contacto.
    !       
       !    
       DO i = 1, todos
            performance(i)%distance= MIN(REAL(calcdistDP( r(performance(i)%p1), r(performance(i)%p2) ) ),&
                                              performance(i)%distance)
       END DO

     END SUBROUTINE findDistMin
!=============================================================================
 
 subroutine fill_listadca(allpairs,natom,dist_ini,listaDCA,maximoindex) 
     !
     use paramSet
     use nrtypes
     !
     IMPLICIT NONE
     !
     integer, intent( IN ) :: allpairs
     integer, intent( IN ) :: natom
     real(dp), dimension(natom,natom) , INTENT(IN) :: dist_ini
     type( dcaout ) , dimension(allpairs), intent(INOUT) :: listadca
     integer , intent(OUT) :: maximoindex
     !
     real(dp) :: tp,tn,fn,fp
     real(dp):: denominador
     integer,dimension(1) :: maximo
     !
     integer :: i,j
     integer :: p1=0
     integer :: p2=0
     integer :: nL=0
     !
     maximo=0
     !  Cuantos pares de residuos hay que tener en cuenta
     !  dado un PDB inicial. Excluir residuos consecutivos
     !  hasta i+, y distancia inicial < tol
     !
     DO i=1,allPairs
        p1 = listadca(i)%p1
        p2 = listadca(i)%p2
        IF ( dist_ini(p1,p2) .le. distanceInitialState ) THEN
            listadca(i)%isinitial=1
        END IF
     END DO
    !
    !DO i=1,25*natom
    DO i=1,allpairs
    !    
      tp=0._dp
      tn=0._dp
      fn=0._dp
      fp=0._dp
    !
      tp=count(listadca(  listadca( 1:i )%indx  )%isinitial==1)
      fp=count(listadca(  listadca( 1:i )%indx  )%isinitial==0)
      tn=count(listadca(  listadca( (i+1) : allpairs  )%indx  )%isinitial==0)
      fn=count(listadca(  listadca( (i+1) : allpairs  )%indx  )%isinitial==1)
    !
    !  Vamos a calcular el MCC para cada par de coev ordenado 
    !  MCC == matthews correlation coefficient  
    !  
      denominador=0._dp
      denominador=(tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)
      denominador=dsqrt( denominador )
    !
      listadca( listadca( i )%indx  )%mcc = (tp*tn - fp*fn) /denominador
     END DO
    !
    maximo= maxloc( listadca( listadca( : )%indx  )%mcc  )
    !
    !  Escribir MCC
    IF ( .FALSE. ) THEN
        OPEN( unit=14 , file = 'mcc.dat')
        DO i = 1, 25*natom
            write(14,*), i, listadca( listadca( i )%indx  )%mcc
        END DO
        CLOSE(14)
    END IF
    !
    nL=maximo(1)
    !
    !   Decidir cuales de todos los contactos de Coev hay que simular  
    !
    maximoindex=nL
    !
    listadca%simulate=.FALSE.
    DO i=1,nL
        !
        p1 = listadca( listadca( i )%indx  )%p1
        p2 = listadca( listadca( i )%indx  )%p2
        !
        IF( dist_ini(p1,p2) > filterbydistance  ) THEN
            !
            listadca( listadca( i )%indx  )%simulate=.TRUE.
        END IF
        !
    END DO   
    !
 END SUBROUTINE fill_listadca
!=============================================================================
   
 subroutine fill_performance(natom,todos,allpairs,dist_ini,listadca,performance)
     !
     use paramSet
     use nrtypes
     !
     IMPLICIT NONE
     !
     integer , intent(in) :: natom
     integer , intent(in) :: todos
     integer , intent(in) :: allpairs
     real(dp), dimension(natom,natom), intent(in) :: dist_ini
     type( dcaout ),dimension(allpairs), intent(IN) :: listaDca  
     type(ROC), dimension(todos),intent(INOUT) :: performance
     !
     integer :: cont
     integer :: i,j
     integer :: p1=0
     integer :: p2=0
     !
     performance%isDCA=0
     performance%p1=-1
     performance%p2=-1
     performance%distance=0._sp
     !
     cont=1
     DO j=2,natom
          DO i=1,j-1
              IF ( j -i .ge. minaadist ) THEN
                  IF ( dist_INI(i,j) > filterbydistance ) THEN
                      performance( cont )%p1 = i
                      performance( cont )%p2 = j
                      cont = cont+1
                  END IF
              END IF
          END DO 
      END DO
     !
     cont=0
     DO i=1,todos
         !
         p1 = performance( i )%p1
         p2 = performance( i )%p2
         !
         performance( i )%distance= dist_ini( p1,p2)
         !
         DO j=1,allpairs
             IF ( listadca( j )%simulate ) THEN
                IF ( p1 == listadca( j )%p1 .and. p2 == listadca( j )%p2) then
                     performance( i )%isdca=1
                     cont=cont+1
                 END IF
             END IF
         END DO
      END DO
    !
  END SUBROUTINE fill_performance
!=============================================================================
  real function calcFence(npairs,summary)
  !
  use nrtypes
  use nr
  !
  implicit none
  !
  integer, intent(IN) :: npairs
  type(resultsTable),dimension(npairs),intent(INOUT) :: summary
  !
  integer :: q1indx,q3indx
  real(sp) :: q1
  real(sp) :: q3
  real(sp) :: iqr
   !
  calcFence=0._sp
    !
    !  Vamos a seleccionar los modelos que mejor se ajustan a la coevolucion
    !  
    !  1) ordenar de forma decreciente la lista de resultados
    
     summary%indx=0
     call indexx_sp(summary%area,summary%indx)
    !
    !
    !  2) Seleccionar los valores en Q1 y Q3
    !  Q1 = AUC( en el indice que deja por debajo 25/100 valores ) 
    !  Q3 = AUC( en el indice que deja por debajo 75/100 valores ) 
    !
     q1indx=0
     q3indx=0
    !
     q1indx = Ceiling( 0.25_sp * npairs)
     q3indx = Floor( 0.75_sp * npairs)
    !
     q1=0._sp
     q3=0._sp
    !
     q1 = summary( summary( q1indx )%indx )%area
     q3 = summary( summary( q3indx )%indx )%area
    !
     IQR=0
     IQR=q3-q1
    !
     calcFence = 0._sp
     calcFence = q3 + 1.5_sp * IQR
    ! 
    !
 END FUNCTION calcFence 
!=============================================================================
 !==================================================================
 !
   SUBROUTINE findDistMinFinalMix(todos,natom,r,performance )
           !
           use geometryDP
           !
           IMPLICIT NONE
           !
           integer, intent(IN) :: todos
           integer, intent(IN) :: natom
           type(pointDP), intent(IN), dimension(natom) :: r
           type(ROC), intent(INOUT), dimension(todos) :: performance
           !
           integer :: i
           !
           
            
    !    1) primero tenemos que saber la distancia minima por residuo 
    !        entre la esctructura inicial y final. Es decir si fue en 
    !        algun momento un contacto.
    !       
       !    
       DO i = 1, todos
            performance(i)%distance= 0.5*calcdistDP(r(performance(i)%p1),r(performance(i)%p2) )+ 0.5*performance(i)%distance
       END DO

     END SUBROUTINE findDistMinFinalMix
!=============================================================================
 !==================================================================
 !
   SUBROUTINE findDistMinFinal(todos,natom,r,performance )
           !
           use geometryDP
           !
           IMPLICIT NONE
           !
           integer, intent(IN) :: todos
           integer, intent(IN) :: natom
           type(pointDP), intent(IN), dimension(natom) :: r
           type(ROC), intent(INOUT), dimension(todos) :: performance
           !
           integer :: i
           !
           
            
    !    1) primero tenemos que saber la distancia minima por residuo 
    !        entre la esctructura inicial y final. Es decir si fue en 
    !        algun momento un contacto.
    !       
       !    
       DO i = 1, todos
            performance(i)%distance= MIN(REAL(calcdistDP( r(performance(i)%p1), r(performance(i)%p2) ) ),&
                                              performance(i)%distance)
       END DO

     END SUBROUTINE findDistMinFinal
!=============================================================================
 
     SUBROUTINE writeitdown(npairs,natom,unit_o,summary,fence,atoms,res,chain,rnum,genModels)
         !
         use paramSet
         !
         IMPLICIT NONE
         !
         integer, intent(in) :: npairs
         integer, intent(in) :: unit_o
         !integer, intent(in) :: unit_models
         !integer, intent(in) :: unit_summary
         integer, intent(in) :: natom
         type(resultstable), intent(in) , dimension(npairs) :: summary
         real, intent(in) :: fence
         character(len=4), intent(IN) :: atoms(natom), res(natom)
         character(len=1), intent(IN) :: chain(natom)
         integer, intent(IN) :: rnum(natom)
         type(newmodels), intent(in),dimension(npairs) :: genModels
         !
         integer :: acceptedModels , cont , i , unit_models, unit_summary
         !
         open(unit=81,file='GeneratedModels.pdb')
         open(unit=82,file='models_report.dat')
         !
         unit_models=81
         unit_summary=82
         !
         acceptedModels=0
         !
         Write(unit_summary,*) " Fence ", fence 
         !
         DO i=1,npairs
             IF ( summary( summary( i )%indx )%area .GE. fence .and.&
                  summary( summary( i )%indx )%dist .LT. maxDistanceTocontact) THEN
                    acceptedModels=acceptedModels+1
                    write(unit_summary,'(I4,I4,2X,I4,2X,2X,3(f10.3,X),f8.5,I5,2X,L1)') acceptedModels, summary( summary( i )%indx )
             END IF
         END DO
        !
        cont=0
        DO i=1,npairs
          IF ( summary( summary( i )%indx )%area .GE. fence .and.&
               summary( summary( i )%indx )%dist .LT. maxDistanceTocontact) THEN
                 cont=cont+1
                 call writeSnapShot(unit_models, cont, genmodels( genmodels(i)%indx  )%rgen, atoms,&
                                res, chain, rnum, natom)
          END IF
        END DO
        !
        Write(unit_o,*) " Escritos ", acceptedModels , " modelos"
    !
    !  vamos a guardar los mejores ? 20 modelos en el caso que que ninguna trayectoria
    !  supere fence
    !
       cont=0
       IF (acceptedModels < 2 ) Then
           DO i=1,25 
                 IF ( summary( summary( i )%indx )%area .GE. fence .and.&
                 summary( summary( i )%indx )%dist .LT. maxDistanceTocontact) THEN
                    cont=cont+1
                    call writeSnapShot(unit_models, cont, genmodels( genmodels(i)%indx  )%rgen, &
                                    atoms, res, chain, rnum, natom)
                    Write(unit_o,*) " NO AUC > FENCE: Best ",cont ," models written instead"
                 END IF
           END DO
       !
       END IF
       !
       close(82)
       close(81)
       !
    END SUBROUTINE writeitdown
!=============================================================================  
END MODULE evolution
!
! OLD CHECK (in main)
 !   IF ( allpairs /=  (natom-minaadist+1)*(natom-minaadist)/2) THEN
  !      WRITE( unit_o, * ) " DCA file do not match with number of residues in PDB file"
 !       STOP
 !   END IF
  !
