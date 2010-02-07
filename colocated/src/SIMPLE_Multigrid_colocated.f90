SUBROUTINE SIMPLE_Multigrid_colocated(ppOUT,uuOUT  ,vvOUT  ,wwOUT,   & !OUT
                                      ppIN, uuIN_co,vvIN_co,wwIN_co, & !IN
                                      NXm,NYm,NZm,                   & !IN
                                      NXm2,NYm2,NZm2,                &
                                      NXm3,NYm3,NZm3,                &
                                      NXm4,NYm4,NZm4,                & 
                                      NXm5,NYm5,NZm5 )                 !IN
! skip right boundary cell, no residuls are computed at NX 
USE COMDAT
INTEGER(4)   ,INTENT(IN)           :: NXm,NYm,NZm
INTEGER(4)   ,INTENT(IN)           :: NXm2,NYm2,NZm2
INTEGER(4)   ,INTENT(IN)           :: NXm3,NYm3,NZm3
INTEGER(4)   ,INTENT(IN)           :: NXm4,NYm4,NZm4
INTEGER(4)   ,INTENT(IN)           :: NXm5,NYm5,NZm5
REAL(8)    ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1) ,INTENT(IN)    :: ppIN
REAL(8)    ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1) ,INTENT(OUT)   :: ppOUT
REAL(8)    ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1) ,INTENT(IN)    :: uuIN_co,wwIN_co,vvIN_co
REAL(8)    ,DIMENSION(0:NXm   ,0:NYm+1  ,0:NZm+1)    :: uuIN
REAL(8)    ,DIMENSION(0:NXm+1 ,0:NYm+1  ,0:NZm  )    :: wwIN
REAL(8)    ,DIMENSION(0:NXm+1 ,0:NYm    ,0:NZm+1)    :: vvIN
REAL(8)    ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1) ,INTENT(OUT)   :: uuOUT,wwOUT,vvOUT
REAL(8)    ,DIMENSION(1:NXm ,1:NYm ,1:NZm )     :: kp
REAL(8)    ,DIMENSION(1:NXm2,1:NYm2,1:NZm2)     :: kp2
REAL(8)    ,DIMENSION(1:NXm3,1:NYm3,1:NZm3)     :: kp3
REAL(8)    ,DIMENSION(1:NXm4,1:NYm4,1:NZm4)     :: kp4
REAL(8)    ,DIMENSION(1:NXm5,1:NYm5,1:NZm5)     :: kp5
REAL(8)    ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1)    :: PPtemp
REAL(8)    ,DIMENSION(1:NXm ,1:NYm ,1:NZm )    :: Residual, Resi, Matrix
REAL(8)    ,DIMENSION(1:NXm2,1:NYm2,1:NZm2)    :: Resi2
REAL(8)    ,DIMENSION(1:NXm3,1:NYm3,1:NZm3)    :: Resi3
REAL(8)    ,DIMENSION(1:NXm4,1:NYm4,1:NZm4)    :: Resi4
REAL(8)    ,DIMENSION(1:NXm5,1:NYm5,1:NZm5)    :: Resi5
REAL(8)    ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1)    :: ppp,PPc
REAL(8)    ,DIMENSION(0:NXm2+1,0:NYm2+1 ,0:NZm2+1)  :: ppp2
REAL(8)    ,DIMENSION(0:NXm3+1,0:NYm3+1 ,0:NZm3+1)  :: ppp3
REAL(8)    ,DIMENSION(0:NXm4+1,0:NYm4+1 ,0:NZm4+1)  :: ppp4
REAL(8)    ,DIMENSION(0:NXm5+1,0:NYm5+1 ,0:NZm5+1)  :: ppp5
REAL(8)    ,DIMENSION(0:NXm   ,0:NYm+1  ,0:NZm+1)    :: uutemp
REAL(8)    ,DIMENSION(0:NXm+1 ,0:NYm+1  ,0:NZm  )    :: wwtemp
REAL(8)    ,DIMENSION(0:NXm+1 ,0:NYm    ,0:NZm+1)    :: vvtemp

INTEGER(4)             :: idum


uuOUT=0;wwOUT=0;vvOUT=0
ppOUT=0;PPc=0
ppp=0;ppp2=0;ppp3=0;ppp4=0;ppp5=0
Resi=0;Resi2=0;Resi3=0;Resi4=0;Resi5=0

uutemp=uuIN;wwtemp=wwIN;vvtemp=vvIN

PPtemp=0
PPtemp(1:NXm,1:NYm,1:NZm)=ppIN(1:NXm,1:NYm,1:NZm)



!!!!!! moving backward one cell produce perfect results

DO i=1,NXm-1
DO m=1,NYm
DO j=1,NZm

uutemp(i,m,j) = 0.5*( uuIN_co(i,m,j) + uuIN_co(i+1,m,j) ) 

ENDDO;ENDDO;ENDDO

DO i=1,NXm
DO m=1,NYm-1
DO j=1,NZm

vvtemp(i,m,j) = 0.5*( vvIN_co(i,m,j) + vvIN_co(i,m+1,j) ) 

ENDDO;ENDDO;ENDDO

DO i=1,NXm
DO m=1,NYm
DO j=1,NZm-1

wwtemp(i,m,j) = 0.5*( wwIN_co(i,m,j) + wwIN_co(i,m,j+1) ) 

ENDDO;ENDDO;ENDDO

 CALL       BC_wall   (uutemp,vvtemp,wwtemp,NXm,NYm,NZm,3)  !both wall
 CALL       VelocityBC(uutemp,vvtemp,wwtemp,NXm,NYm,NZm)

uuIN=uutemp;vvIN=vvtemp;wwIN=wwtemp

IF((ANH==1).OR.(ANH==2)) THEN ! NonHydrostatic without Hydrostatic Correction and Adjustable NonHydrostatic.
        ! find kp
        !---------------------------------------------------
        CALL INIT_KP_All(kp, NXm, NYm, NZm)        
        IF (Multilevel>1) CALL INIT_KP_All(kp2,NXm2,NYm2,NZm2)
        IF (Multilevel>2) CALL INIT_KP_All(kp3,NXm3,NYm3,NZm3)         
        IF (Multilevel>3) CALL INIT_KP_All(kp4,NXm4,NYm4,NZm4)
        IF (Multilevel>4) CALL INIT_KP_All(kp5,NXm5,NYm5,NZm5)
        !---------------------------------------------------

   
    LOOP_1: DO iSOR=1,ii_SOR 
!************************************************************************************
    ! Find Residual
   !(============================================================================
  If (SurfaceCorrection==0) THEN
        DO m=1,NYm;DO i=1,NXm;DO j=NZm,NZm
           wwtemp(i,m,j)=wwtemp(i,m,j-1) - (uutemp(i,m,j)-uutemp(i-1,m,j))*dz/dx   &
                                         - (vvtemp(i,m,j)-vvtemp(i,m-1,j))*dz/dy
        ENDDO;ENDDO;ENDDO 
  ENDIF

    Residual = 0
        DO m=1,NYm;DO j=1,NZm;DO i=1,NXm !i=1,nx
  
        Residual(i,m,j) = (   (1./dx)*(uutemp(i,m,j)-uutemp(i-1,m,j))   &
                             +(1./dy)*(vvtemp(i,m,j)-vvtemp(i,m-1,j))   &
                             +(1./dz)*(wwtemp(i,m,j)-wwtemp(i,m,j-1))  )
        ENDDO;ENDDO;ENDDO
   !(============================================================================  

        err=MAXVAL(ABS(Residual)) 

        IF (iSOR==1) THEN
            WRITE(18,*) err, iSOR
            WRITE(*,*) err, iSOR
        ENDIF 
 
 !************************************************************************************
 IF (ISOR_ErrMax==211) THEN

        IF ((err<=errmax).and.(iSOR>=ii_SOR_Min)) then
            WRITE(18,*) err, iSOR
            WRITE(*,*) err, iSOR
            EXIT LOOP_1
        ENDIF
 
 ELSEIF (ISOR_ErrMax==210) THEN

        IF (iSOR>=ii_SOR_Min) then
            WRITE(18,*) err, iSOR
            WRITE(*,*) err, iSOR
            EXIT LOOP_1
        ENDIF 
 
 ELSEIF (ISOR_ErrMax==201) THEN

        IF ((err<=errmax)) then
            WRITE(18,*) err, iSOR
            WRITE(*,*) err, iSOR
            EXIT LOOP_1
        ENDIF 
 
 ELSEIF (ISOR_ErrMax==200) THEN

        IF ((err<=errmax).or.(iSOR>=ii_SOR_Min)) then
            WRITE(18,*) err, iSOR
            WRITE(*,*) err, iSOR
            EXIT LOOP_1
        ENDIF
 ELSE
    WRITE(*,*) 'check ISOR_ErrMax';STOP
 ENDIF
!************************************************************************************ 
 
        Resi=Residual*den_w/dt
!ppp=0;ppp2=0;ppp3=0;ppp4=0;ppp5=0


        !Fine Grid (
        DO idum = 1, i_multi1
         CALL POISSON_Revised_Multigrid(ppp,                      & !INOUT
                                        Resi,kp,NXm,NYm,NZm)        !IN
        ENDDO
        !)
        !-----------------------------------------------
        CALL MatrixUpdate(Matrix,              & !OUT
                          ppp,kp,NXm,NYm,NZm)    !IN 
        Resi = Resi - Matrix                                               
        PPtemp=PPtemp+ppp; PPc=PPc+ppp
        !-----------------------------------------------

IF(Multilevel==1 .or. err<errLevel1) GOTO 915
        
        !Coarse2 Grid (
        CALL AveMatrix2Small(Resi2,NXm2,NYm2,NZm2,Resi,NXm,NYm,NZm)
        
        CALL AveMatrix2SmallP(ppp2,NXm2,NYm2,NZm2,ppp,NXm,NYm,NZm)
        DO idum = 1, i_multi2
         CALL POISSON_Revised_Multigrid(ppp2,                          & !INOUT
                                        Resi2,kp2,NXm2,NYm2,NZm2)        !IN
        ENDDO
        CALL CloneMatrix2Large(ppp,NXm,NYm,NZm,ppp2,NXm2,NYm2,NZm2)
        !)

        !Fine Grid (
        DO idum = 1, i_multi1
         CALL POISSON_Revised_Multigrid(ppp,                      & !INOUT
                                        Resi,kp,NXm,NYm,NZm)        !IN
        ENDDO
        !)
        !-----------------------------------------------
        CALL MatrixUpdate(Matrix,              & !OUT
                          ppp,kp,NXm,NYm,NZm)    !IN 
        Resi = Resi - Matrix                                               
        PPtemp=PPtemp+ppp; PPc=PPc+ppp
        !-----------------------------------------------

IF(Multilevel==2 .or. err<errLevel2) GOTO 915

        !Coarse3 Grid (
        CALL AveMatrix2Small(Resi2,NXm2,NYm2,NZm2,Resi,NXm,NYm,NZm)
        CALL AveMatrix2Small(Resi3,NXm3,NYm3,NZm3,Resi2,NXm2,NYm2,NZm2)
        
        CALL AveMatrix2SmallP(ppp3,NXm3,NYm3,NZm3,ppp2,NXm2,NYm2,NZm2)
        DO idum = 1, i_multi3
         CALL POISSON_Revised_Multigrid(ppp3,                          & !INOUT
                                        Resi3,kp3,NXm3,NYm3,NZm3)        !IN
        ENDDO
        CALL CloneMatrix2Large(ppp2,NXm2,NYm2,NZm2,ppp3,NXm3,NYm3,NZm3)
        !)

        !Coarse2 Grid (

        !CALL AveMatrix2Small(Resi2,NXm2,NYm2,NZm2,Resi,NXm,NYm,NZm)
        !CALL AveMatrix2SmallP(ppp2,NXm2,NYm2,NZm2,ppp,NXm,NYm,NZm)
        DO idum = 1, i_multi2
         CALL POISSON_Revised_Multigrid(ppp2,                          & !INOUT
                                       Resi2,kp2,NXm2,NYm2,NZm2)        !IN
        ENDDO
        CALL CloneMatrix2Large(ppp,NXm,NYm,NZm,ppp2,NXm2,NYm2,NZm2)
        !)

        !Fine Grid (
        DO idum = 1, i_multi1
         CALL POISSON_Revised_Multigrid(ppp,                      & !INOUT
                                        Resi,kp,NXm,NYm,NZm)        !IN
        ENDDO
        !)
        !-----------------------------------------------
        CALL MatrixUpdate(Matrix,              & !OUT
                          ppp,kp,NXm,NYm,NZm)    !IN 
        Resi = Resi - Matrix                                               
        PPtemp=PPtemp+ppp; PPc=PPc+ppp
        !-----------------------------------------------

IF(Multilevel==3 .or. err<errLevel3) GOTO 915

        !Coarse4 Grid (
        CALL AveMatrix2Small(Resi2,NXm2,NYm2,NZm2,Resi,NXm,NYm,NZm)
        CALL AveMatrix2Small(Resi3,NXm3,NYm3,NZm3,Resi2,NXm2,NYm2,NZm2)
        CALL AveMatrix2Small(Resi4,NXm4,NYm4,NZm4,Resi3,NXm3,NYm3,NZm3)
        
        CALL AveMatrix2SmallP(ppp4,NXm4,NYm4,NZm4,ppp3,NXm3,NYm3,NZm3)
        DO idum = 1, i_multi4
         CALL POISSON_Revised_Multigrid(ppp4,                          & !INOUT
                                       Resi4,kp4,NXm4,NYm4,NZm4)        !IN
        ENDDO
        CALL CloneMatrix2Large(ppp3,NXm3,NYm3,NZm3,ppp4,NXm4,NYm4,NZm4)
        !)
        
        !Coarse3 Grid (
        !CALL AveMatrix2Small(Resi3,NXm3,NYm3,NZm3,Resi2,NXm2,NYm2,NZm2)
        !CALL AveMatrix2SmallP(ppp3,NXm3,NYm3,NZm3,ppp2,NXm2,NYm2,NZm2)
        DO idum = 1, i_multi3
         CALL POISSON_Revised_Multigrid(ppp3,                          & !INOUT
                                        Resi3,kp3,NXm3,NYm3,NZm3)        !IN
        ENDDO
        CALL CloneMatrix2Large(ppp2,NXm2,NYm2,NZm2,ppp3,NXm3,NYm3,NZm3)
        !)

        !Coarse2 Grid (
        !CALL AveMatrix2Small(Resi2,NXm2,NYm2,NZm2,Resi,NXm,NYm,NZm)
        !CALL AveMatrix2SmallP(ppp2,NXm2,NYm2,NZm2,ppp,NXm,NYm,NZm)
        DO idum = 1, i_multi2
        CALL POISSON_Revised_Multigrid(ppp2,                          & !INOUT
                                       Resi2,kp2,NXm2,NYm2,NZm2)        !IN
        ENDDO
        CALL CloneMatrix2Large(ppp,NXm,NYm,NZm,ppp2,NXm2,NYm2,NZm2)
        !)

        !Fine Grid (
        DO idum = 1, i_multi1
         CALL POISSON_Revised_Multigrid(ppp,                      & !INOUT
                                        Resi,kp,NXm,NYm,NZm)        !IN
        ENDDO
        !)
        !-----------------------------------------------
        CALL MatrixUpdate(Matrix,              & !OUT
                          ppp,kp,NXm,NYm,NZm)    !IN 
        Resi = Resi - Matrix                                               
        PPtemp=PPtemp+ppp; PPc=PPc+ppp
        !-----------------------------------------------        

IF(Multilevel==4 .or. err<errLevel4) GOTO 915

        !Coarse5 Grid (
        CALL AveMatrix2Small(Resi2,NXm2,NYm2,NZm2,Resi,NXm,NYm,NZm)
        CALL AveMatrix2Small(Resi3,NXm3,NYm3,NZm3,Resi2,NXm2,NYm2,NZm2)
        CALL AveMatrix2Small(Resi4,NXm4,NYm4,NZm4,Resi3,NXm3,NYm3,NZm3)       
        CALL AveMatrix2Small(Resi5,NXm5,NYm5,NZm5,Resi4,NXm4,NYm4,NZm4)
        
        CALL AveMatrix2SmallP(ppp5,NXm5,NYm5,NZm5,ppp4,NXm4,NYm4,NZm4)
        DO idum = 1, i_multi5
         CALL POISSON_Revised_Multigrid(ppp5,                          & !INOUT
                                        Resi5,kp5,NXm5,NYm5,NZm5)        !IN
        ENDDO
        CALL CloneMatrix2Large(ppp4,NXm4,NYm4,NZm4,ppp5,NXm5,NYm5,NZm5)
        !)

        !Coarse4 Grid (
        !CALL AveMatrix2Small(Resi4,NXm4,NYm4,NZm4,Resi3,NXm3,NYm3,NZm3)
        !CALL AveMatrix2SmallP(ppp4,NXm4,NYm4,NZm4,ppp3,NXm3,NYm3,NZm3)
        DO idum = 1, i_multi4
         CALL POISSON_Revised_Multigrid(ppp4,                          & !INOUT
                                        Resi4,kp4,NXm4,NYm4,NZm4)        !IN
        ENDDO
        CALL CloneMatrix2Large(ppp3,NXm3,NYm3,NZm3,ppp4,NXm4,NYm4,NZm4)
        !)
        
        !Coarse3 Grid (
        !CALL AveMatrix2Small(Resi3,NXm3,NYm3,NZm3,Resi2,NXm2,NYm2,NZm2)
        !CALL AveMatrix2SmallP(ppp3,NXm3,NYm3,NZm3,ppp2,NXm2,NYm2,NZm2)
        DO idum = 1, i_multi3
         CALL POISSON_Revised_Multigrid(ppp3,                          & !INOUT
                                        Resi3,kp3,NXm3,NYm3,NZm3)        !IN
        ENDDO
        CALL CloneMatrix2Large(ppp2,NXm2,NYm2,NZm2,ppp3,NXm3,NYm3,NZm3)
        !)        

        !Coarse2 Grid (
        DO idum = 1, i_multi2
         CALL POISSON_Revised_Multigrid(ppp2,                          & !INOUT
                                        Resi2,kp2,NXm2,NYm2,NZm2)        !IN
        ENDDO
        CALL CloneMatrix2Large(ppp,NXm,NYm,NZm,ppp2,NXm2,NYm2,NZm2)
        !)
              
        !Fine Grid (
        DO idum = 1, i_multi1
         CALL POISSON_Revised_Multigrid(ppp,                      & !INOUT
                                        Resi,kp,NXm,NYm,NZm)        !IN
        ENDDO
        !)

        !-----------------------------------------------
        !CALL MatrixUpdate(Matrix,              & !OUT
        !                  ppp,kp,NXm,NYm,NZm)    !IN 
        !Resi = Resi - Matrix                                               
        PPtemp=PPtemp+ppp; PPc=PPc+ppp
        !-----------------------------------------------
        !Residual = Resi*dt/den_w

915     CONTINUE
   !)============================================================================ 
        !(4) Correct uu,vv,ww
        DO m=1,NYm;DO i=1,NXm-1;DO j=1,NZm
           uutemp(i,m,j)=uuIN(i,m,j) -dt/dx*(PPc(i+1,m,j)-PPc(i,m,j))/den_w
        ENDDO;ENDDO;ENDDO
        DO m=1,NYm-1;DO i=1,NXm;DO j=1,NZm
           vvtemp(i,m,j)=vvIN(i,m,j) -dt/dy*(PPc(i,m+1,j)-PPc(i,m,j))/den_w
        ENDDO;ENDDO;ENDDO        
        DO m=1,NYm;DO i=1,NXm;DO j=1,NZm-1 
           wwtemp(i,m,j)=wwIN(i,m,j)-dt/dz*(PPc(i,m,j+1)-PPc(i,m,j))/den_w
        ENDDO;ENDDO;ENDDO
  If (SurfaceCorrection==0) THEN
        DO m=1,NYm;DO i=1,NXm;DO j=NZm,NZm
           wwtemp(i,m,j)=wwtemp(i,m,j-1) - (uutemp(i,m,j)-uutemp(i-1,m,j))*dz/dx   &
                                         - (vvtemp(i,m,j)-vvtemp(i,m-1,j))*dz/dy
        ENDDO;ENDDO;ENDDO 
  ELSE 
        DO m=1,NYm;DO i=1,NXm;DO j=NZm,NZm 
           wwtemp(i,m,j)=wwIN(i,m,j)-dt/dz*(PPc(i,m,j+1)-PPc(i,m,j))/den_w
        ENDDO;ENDDO;ENDDO 
  ENDIF      
   !)============================================================================ 
    ENDDO LOOP_1

    
    
    IF( ANH == 2) THEN !Nonhydrostatic with Hydrostatic Correction -> Adjustable NonHydrostatic

        DO j=1,NZm !-1
        DO m=1,NYm
        DO i=1,NXm   
             wwtemp(i,m,j)   = wwtemp(i,m,j-1) -  (uutemp(i,m,j)-uutemp(i-1,m,j))*dz/dx   &
                                               -  (vvtemp(i,m,j)-vvtemp(i,m-1,j))*dz/dy
        ENDDO;ENDDO;ENDDO
        
    ENDIF

ELSEIF (ANH == 0) THEN ! Fully Hydrostatic

        DO j=1,NZm !-1
        DO m=1,NYm
        DO i=1,NXm   
             wwtemp(i,m,j)   = wwtemp(i,m,j-1) -  (uutemp(i,m,j)-uutemp(i-1,m,j))*dz/dx   &
                                               -  (vvtemp(i,m,j)-vvtemp(i,m-1,j))*dz/dy
        ENDDO;ENDDO;ENDDO

ENDIF    

! correct colocated

PPc(0,:,:)     = PPc(1,:,:)
PPc(NXm+1,:,:) = PPc(NXm,:,:)
PPc(:,0,:)     = PPc(:,1,:)
PPc(:,NYm+1,:) = PPc(:,NYm,:)
PPc(:,:,0)     = PPc(:,:,1)
PPc(:,:,NZm+1) = 0

        DO m=1,NYm;DO i=1,NXm;DO j=1,NZm
           uuOUT(i,m,j)=uuIN_co(i,m,j) -dt/dx*0.5*(PPc(i+1,m,j)-PPc(i-1,m,j))/den_w

           vvOUT(i,m,j)=vvIN_co(i,m,j) -dt/dy*0.5*(PPc(i,m+1,j)-PPc(i,m-1,j))/den_w

           wwOUT(i,m,j)=wwIN_co(i,m,j) -dt/dz*0.5*(PPc(i,m,j+1)-PPc(i,m,j-1))/den_w
        ENDDO;ENDDO;ENDDO

 CALL       BC_wall_colocated   (uuOUT,vvOUT,wwOUT,NXm,NYm,NZm,3) 
 CALL       VelocityBC_colocated(uuOUT,vvOUT,wwOUT,NXm,NYm,NZm) 

ppOUT=PPtemp
!uuOUT=uutemp
!vvOUT=vvtemp
!wwOUT=wwtemp

CONTINUE

ENDSUBROUTINE SIMPLE_Multigrid_colocated