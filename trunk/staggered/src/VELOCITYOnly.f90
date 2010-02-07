SUBROUTINE VELOCITYONLY(uu_OUT,vv_OUT,ww_OUT,pp_OUT,                  & !OUT
                         uu_old,vv_old,ww_old,                         & !IN old 
                         uu_m,  vv_m,  ww_m,  pp_m, volsum_m, den_m,   & !IN m
                         pp_solve)                            !IN
!_________________________________________________________________________________________
USE COMDAT
INTEGER(4)                                  ,INTENT(IN)    :: pp_solve
REAL(8)    ,DIMENSION(0:nx+1,0:NY+1,0:NZ+1) ,INTENT(IN)    :: pp_m
REAL(8)    ,DIMENSION(0:nx+1,0:NY+1,0:NZ+1) ,INTENT(IN)    :: den_m
REAL(8)    ,DIMENSION(0:NX  ,0:NY+1,0:NZ+1) ,INTENT(IN)    :: uu_old,uu_m
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ  ) ,INTENT(IN)    :: ww_old,vv_m
REAL(8)    ,DIMENSION(0:NX+1,0:NY  ,0:NZ+1) ,INTENT(IN)    :: vv_old,ww_m
REAL(8)    ,DIMENSION(0:NX  ,0:NY+1,0:NZ+1) ,INTENT(OUT)   :: uu_OUT
REAL(8)    ,DIMENSION(0:NX+1,0:NY  ,0:NZ+1) ,INTENT(OUT)   :: vv_OUT
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ  ) ,INTENT(OUT)   :: ww_OUT
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1) ,INTENT(OUT)   :: pp_OUT
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1)        ,INTENT(IN)    :: volsum_m

REAL(8)    ,DIMENSION(0:NX  ,0:NY+1,0:NZ+1)    :: uustar
REAL(8)    ,DIMENSION(0:NX+1,0:NY  ,0:NZ+1)    :: vvstar
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ  )    :: wwstar
REAL(8)    ,DIMENSION(0:NX  ,0:NY  ,0:NZ+1)    :: BaCX,BaCY
REAL(8)    ,DIMENSION(1:NX  ,1:NY  ,2)         :: BaT
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1)    :: voll_m
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1)    :: ppIN



!-----------------------------------------------------(
 CALL ASSIGNVOL_R(voll_m, volsum_m)
!-----------------------------------------------------) 

!-----------------------------------------------------(
 CALL BaroTropic_R(BaT,                       &  !OUT                                 
                   volsum_m)                     !IN
 CALL BaroClinic_R(BaCX,BaCY,                 &  !OUT                                 
                   den_m)              !IN  ! Be careful that topcell is not input into velocity subroutine

!-----------------------------------------------------)                   
                   
!-------------------------------------------------------------------------(
uustar=0;vvstar=0;wwstar=0


 CALL UPREDICT_CENTRAL(uustar,vvstar,wwstar,                      & !OUT
                       uu_old,vv_old,ww_old,uu_m,vv_m,ww_m,pp_m,  & !IN
                       BaCX,BaCY,BaT,den_m)               !IN

!-------------------------------------------------------------------------)

!-------------------------------------------------------------(
IF(pp_solve>0) THEN
 ppIN=pp_m
 IF(pp_reset) ppIN=0  ! 05121902 
  CALL SIMPLE_R_S_Multigrid_incre(pp_OUT, uu_OUT, vv_OUT, ww_OUT,  & !OUT
                            ppIN,   uustar, vvstar, wwstar)             !IN
ELSE
 pp_OUT=pp_m
ENDIF
!-------------------------------------------------------------)
  


ENDSUBROUTINE VELOCITYONLY