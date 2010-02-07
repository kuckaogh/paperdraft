SUBROUTINE VELOCITYONLY_colocated(uu_OUT,vv_OUT,ww_OUT,pp_OUT,           & !OUT
                         uu_old,vv_old,ww_old,                           & !IN old 
                         uu_m,  vv_m,  ww_m,  pp_m, volsum_m, den_m,     & !IN m
                         NXm,NYm,NZm,wall,pp_solve)                        !IN
!_________________________________________________________________________________________
USE COMDAT
INTEGER(4)   ,INTENT(IN)           :: NXm,NYm,NZm,wall
INTEGER(4)                                     ,INTENT(IN)    :: pp_solve
REAL(8)    ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1) ,INTENT(IN)    :: pp_m
REAL(8)    ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1) ,INTENT(IN)    :: den_m
REAL(8)    ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1) ,INTENT(IN)    :: uu_old,uu_m
REAL(8)    ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1) ,INTENT(IN)    :: ww_old,vv_m
REAL(8)    ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1) ,INTENT(IN)    :: vv_old,ww_m
REAL(8)    ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1) ,INTENT(OUT)   :: uu_OUT
REAL(8)    ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1) ,INTENT(OUT)   :: vv_OUT
REAL(8)    ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1) ,INTENT(OUT)   :: ww_OUT
REAL(8)    ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1) ,INTENT(OUT)   :: pp_OUT
REAL(8)    ,DIMENSION(0:NXm+1,0:NYm+1)         ,INTENT(IN)    :: volsum_m

REAL(8)    ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1)    :: uustar,vvstar,wwstar
!REAL(8)    ,DIMENSION(0:NXm  ,0:NYm  ,0:NZm+1)    :: BaCX,BaCY
REAL(8)    ,DIMENSION(1:NXm  ,1:NYm  ,2)         :: BaT
REAL(8)    ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1)    :: voll_m
REAL(8)    ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1)    :: ppIN
INTEGER(4)             :: NXm2,NYm2,NZm2
INTEGER(4)             :: NXm3,NYm3,NZm3
INTEGER(4)             :: NXm4,NYm4,NZm4
INTEGER(4)             :: NXm5,NYm5,NZm5


!-----------------------------------------------------(
 CALL ASSIGNVOL(voll_m, volsum_m,NXm ,NYm ,NZm)
!-----------------------------------------------------) 

!-----------------------------------------------------(
 CALL BaroTropic_colocated(BaT,                       &  !OUT                                 
                   volsum_m,NXm ,NYm)                     !IN
!-----------------------------------------------------)                   
                   
!-------------------------------------------------------------------------(
  CALL UPREDICT_colocated(uustar,vvstar,wwstar,                      & !OUT
                          uu_old,vv_old,ww_old,uu_m,vv_m,ww_m,pp_m,  & !IN
                          BaT,den_m,NXm,NYm,NZm,wall)                  !IN
!-------------------------------------------------------------------------)

!-------------------------------------------------------------(
IF(pp_solve>0) THEN

        ppIN=pp_m

        NXm2=FLOOR(REAL(NXm)/2); NYm2=MAX(FLOOR(REAL(NYm)/2),1); NZm2=FLOOR(REAL(NZm)/2)
        NXm3=FLOOR(REAL(NXm2)/2);NYm3=MAX(FLOOR(REAL(NYm2)/2),1);NZm3=FLOOR(REAL(NZm2)/2)
        NXm4=FLOOR(REAL(NXm3)/2);NYm4=MAX(FLOOR(REAL(NYm3)/2),1);NZm4=FLOOR(REAL(NZm3)/2)
        NXm5=FLOOR(REAL(NXm4)/2);NYm5=MAX(FLOOR(REAL(NYm4)/2),1);NZm5=FLOOR(REAL(NZm4)/2)

!DO i=0,NXm
!DO m=0,NYm+1
!DO j=0,NZm+1
!
!uuI(i,m,j) = 0.5*( uu_OUT(i,m,j) + uu_OUT(i+1,m,j) ) 
!
!ENDDO;ENDDO;ENDDO
!
!DO i=0,NXm+1
!DO m=0,NYm
!DO j=0,NZm+1
!
!vvI(i,m,j) = 0.5*( vv_OUT(i,m,j) + vv_OUT(i,m+1,j) ) 
!
!ENDDO;ENDDO;ENDDO
!
!DO i=0,NXm+1
!DO m=0,NYm+1
!DO j=0,NZm
!
!wwI(i,m,j) = 0.5*( ww_OUT(i,m,j) + ww_OUT(i,m,j+1) ) 
!
!ENDDO;ENDDO;ENDDO

  CALL SIMPLE_Multigrid_colocated(pp_OUT, uu_OUT, vv_OUT, ww_OUT,        & !OUT
                                  ppIN,   uustar, vvstar, wwstar,        & !IN
                                  NXm,NYm,NZm,                           & !IN
                                  NXm2,NYm2,NZm2,                        &
                                  NXm3,NYm3,NZm3,                        &
                                  NXm4,NYm4,NZm4,                        & 
                                  NXm5,NYm5,NZm5 )
ELSE
 pp_OUT=pp_m
ENDIF
!-------------------------------------------------------------)
  

 CALL       BC_wall_colocated   (uu_OUT,vv_OUT,ww_OUT,NXm,NYm,NZm,wall) 
 CALL       VelocityBC_colocated(uu_OUT,vv_OUT,ww_OUT,NXm,NYm,NZm) 

ENDSUBROUTINE VELOCITYONLY_colocated