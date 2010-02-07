SUBROUTINE TimeIntegration_colocated(uu  ,vv  ,ww,  pp  ,volsum  ,den  ,ff  ,ffX  ,ffY  ,ffZ,   & !OUT
                           uu_0,vv_0,ww_0,pp_0,volsum_0,den_0,ff_0,ffX_0,ffY_0,ffZ_0, &
                           NXm ,NYm ,NZm, wall)   !IN

USE COMDAT
IMPLICIT NONE !{
INTEGER(4)   ,INTENT(IN)           :: NXm,NYm,NZm,wall
REAL(8)  ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1) ,INTENT(IN)  :: pp_0 
REAL(8)  ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1) ,INTENT(IN)  :: den_0 
REAL(8)  ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1) ,INTENT(IN)  :: uu_0,ww_0,vv_0 
REAL(8)  ,DIMENSION(0:NXm+1,0:NYm+1)         ,INTENT(IN)  :: volsum_0 
REAL(8)  ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1) ,INTENT(IN)  :: ff_0,ffX_0,ffY_0,ffZ_0
REAL(8)  ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1) ,INTENT(OUT) :: uu,vv,ww
REAL(8)  ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1) ,INTENT(OUT) :: den
REAL(8)  ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1) ,INTENT(OUT) :: pp
REAL(8)  ,DIMENSION(0:NXm+1,0:NYm+1)         ,INTENT(OUT) :: volsum
REAL(8)  ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1) ,INTENT(OUT) :: ff,ffX,ffY,ffZ


REAL(8)  ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1)  :: voll
REAL(8)  ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1)  :: uu_m1,uu_m2,uu_m3,uu_m4,uu_m5,uu_I
REAL(8)  ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1)  :: vv_m1,vv_m2,vv_m3,vv_m4,vv_m5,vv_I
REAL(8)  ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1)  :: ww_m1,ww_m2,ww_m3,ww_m4,ww_m5,ww_I
REAL(8)  ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1)  :: pp_m1,pp_m2,pp_m3,pp_m4,pp_m5 !,pp_I
REAL(8)  ,DIMENSION(0:NXm+1,0:NYm+1)          :: volsum_m1,volsum_m2,volsum_m3,volsum_m4,volsum_m5,volsum_I
REAL(8)  ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1)  :: den_m1,den_m2,den_m3,den_m4,den_m5,den_I
REAL(8)  ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1)  :: ff_m1,ffX_m1,ffY_m1,ffZ_m1
REAL(8)  ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1)  :: ff_m2,ffX_m2,ffY_m2,ffZ_m2
REAL(8)  ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1)  :: ff_m3,ffX_m3,ffY_m3,ffZ_m3
REAL(8)  ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1)  :: ff_I,ffX_I,ffY_I,ffZ_I


dt  = dt_o 

 SELECTCASE(TIME_ORDER)
 CASE(1) ! Euler
  STOP

 CASE(2) ! Predictor-Corrector 
  STOP

 CASE(3)
!---------------------------------------------------------( 

!CFL=2
!Ruuth and Spiteri  4 stage 3rd order optimal SSPRK
!
!U_1    =     U^n           + 1/2 dt F(U^n),
!U_2    =     U_1           + 1/2 dt F(U_1),
!U_3    = 2/3 U^n + 1/3 U_2 + 1/6 dt F(U_2),
!U^n+1  =     U_3           + 1/2 dt F(U_3),
 !============================================================================ 

 !U 1
dt=0.5*dt_o
 !============================================================================
  CALL VELOCITYONLY_colocated(uu_m1 ,vv_m1 ,ww_m1 ,pp_m1,                     & !OUT
                     uu_0  ,vv_0  ,ww_0  ,                           & !IN old 
                     uu_0  ,vv_0  ,ww_0  ,pp_0  ,volsum_0  ,den_0  , & !IN m
                     NXm,NYm,NZm,wall,pp_solve_1)
 !============================================================================                     
  CALL HEIGHTVOL_colocated(volsum_m1,voll,uu_0,vv_0,volsum_0,NXm,NYm,NZm)
  IF (height_smoothing==1) CALL HEIGHT_SMOOTH(volsum_m1,NXm,NYm)
 !============================================================================ 
  CALL TRANSPORT_CIP_colocated(ff_m1,  ffX_m1,  ffY_m1,  ffZ_m1,        & ! OUT
                         ff_0  , ffX_0  , ffY_0  , ffZ_0  ,       & ! IN
                        uu_0  ,  vv_0  ,  ww_0,NXm,NYm,NZm)          ! IN m
  CALL DENSITYSAL_R(den_m1,ff_m1,NXm,NYm,NZm)
 !============================================================================ 
 !U 2
dt=0.5*dt_o 
 !============================================================================                    
  CALL VELOCITYONLY_colocated(uu_m2, vv_m2, ww_m2, pp_m2,                     & !OUT
                     uu_m1 ,vv_m1 ,ww_m1 ,                           & !IN old
                     uu_m1 ,vv_m1 ,ww_m1 ,pp_m1, volsum_m1, den_m1,  & !IN m
                     NXm,NYm,NZm,wall,pp_solve_2)
 !============================================================================                     
  CALL HEIGHTVOL_colocated(volsum_m2,voll,uu_m1,vv_m1,volsum_m1,NXm,NYm,NZm)         !IN
  IF (height_smoothing==1) CALL HEIGHT_SMOOTH(volsum_m2,NXm,NYm)
 !============================================================================ 
  CALL TRANSPORT_CIP_colocated(ff_m2,  ffX_m2,  ffY_m2,  ffZ_m2,        & ! OUT 
                         ff_m1,  ffX_m1,  ffY_m1,  ffZ_m1,        & ! IN
                        uu_m1,  vv_m1,  ww_m1,NXm,NYm,NZm)          ! IN m
  CALL DENSITYSAL_R(den_m2,ff_m2,NXm,NYm,NZm)
 !============================================================================ 
 !U 3
dt=1./6*dt_o
 uu_I     = 2./3*uu_0     + 1./3*uu_m2
 vv_I     = 2./3*vv_0     + 1./3*vv_m2
 ww_I     = 2./3*ww_0     + 1./3*ww_m2 
 den_I    = 2./3*den_0    + 1./3*den_m2   
 volsum_I = 2./3*volsum_0 + 1./3*volsum_m2
 ff_I     = 2./3*ff_0     + 1./3*ff_m2      
 ffX_I    = 2./3*ffX_0    + 1./3*ffX_m2
 ffY_I    = 2./3*ffY_0    + 1./3*ffY_m2
 ffZ_I    = 2./3*ffZ_0    + 1./3*ffZ_m2
 !============================================================================                     
  CALL VELOCITYONLY_colocated(uu_m3, vv_m3, ww_m3, pp_m3,                     & !OUT
                     uu_I  ,vv_I  ,ww_I,                            & !IN old 
                     uu_m2, vv_m2, ww_m2, pp_m2, volsum_m2, den_m2,  & !IN m
                     NXm,NYm,NZm,wall,pp_solve_3)                     
 !============================================================================                     
  CALL HEIGHTVOL_colocated(volsum_m3,voll,uu_m2,vv_m2,volsum_I,NXm,NYm,NZm)
  IF (height_smoothing==1) CALL HEIGHT_SMOOTH(volsum_m3,NXm,NYm)
 !============================================================================ 
  CALL TRANSPORT_CIP_colocated(ff_m3,  ffX_m3,  ffY_m3,  ffZ_m3,        & ! OUT 
                         ff_I  , ffX_I  , ffY_I  , ffZ_I  ,       & ! IN
                        uu_m2,  vv_m2,  ww_m2,NXm,NYm,NZm)          ! IN m
  CALL DENSITYSAL_R(den_m3,ff_m3,NXm,NYm,NZm)
 !============================================================================ 
 !U 4 - final
 dt=0.5*dt_o
 !============================================================================                     
  CALL VELOCITYONLY_colocated(uu  ,   vv,     ww     , pp,                     & !OUT
                    uu_m3,  vv_m3,  ww_m3  ,                           & !IN old 
                    uu_m3,  vv_m3,  ww_m3, pp_m3, volsum_m3, den_m3,  & !IN m
                    NXm,NYm,NZm,wall,1)                     
 !============================================================================                     
  CALL HEIGHTVOL_colocated(volsum,voll,uu_m3,vv_m3,volsum_m3,NXm,NYm,NZm)
  IF (height_smoothing==1) CALL HEIGHT_SMOOTH(volsum,NXm,NYm)
 !============================================================================ 
  CALL TRANSPORT_CIP_colocated(ff,    ffX,    ffY,    ffZ,           & ! OUT !Recycle ff_m1
                         ff_m3, ffX_m3, ffY_m3, ffZ_m3,        & ! IN
                        uu_m3,  vv_m3,  ww_m3,NXm,NYm,NZm)           ! IN m
  CALL DENSITYSAL_R(den,ff,NXm,NYm,NZm)
 !============================================================================


!---------------------------------------------------------)
  
 CASE(4) ! RK4 System
!---------------------------------------------------------(

 !K 1
 !============================================================================   
 uu_m1     = uu_0
 vv_m1     = vv_0
 ww_m1     = ww_0    
 pp_m1     = pp_0    
 den_m1    = den_0   
 volsum_m1 = volsum_0

                   
 !K 2
 dt=0.5*dt_o 
 !============================================================================                    
  CALL VELOCITYONLY_colocated(uu_m2, vv_m2, ww_m2, pp_m2,                     & !OUT
                     uu_0  ,vv_0  ,ww_0  ,                           & !IN old
                     uu_m1, vv_m1, ww_m1, pp_m1, volsum_m1, den_m1,  & !IN m
                     NXm,NYm,NZm,wall,pp_solve_2)
 !============================================================================                     
  CALL HEIGHTVOL_colocated(volsum_m2,voll,uu_m1,vv_m1,volsum_0,NXm,NYm,NZm)
  IF (height_smoothing==1) CALL HEIGHT_SMOOTH(volsum_m2,NXm,NYm)
 !============================================================================ 
  CALL TRANSPORT_CIP_colocated(ff_m1,  ffX_m1,  ffY_m1,  ffZ_m1,        & ! OUT !Recycle ff_m1
                         ff_0  , ffX_0  , ffY_0  , ffZ_0  ,       & ! IN
                        uu_m1,  vv_m1,  ww_m1,NXm,NYm,NZm)          ! IN m
  CALL DENSITYSAL_R(den_m2,ff_m1,NXm,NYm,NZm)
 !============================================================================ 
                     
 !K 3
 dt=0.5*dt_o
 !============================================================================                     
  CALL VELOCITYONLY_colocated(uu_m3, vv_m3, ww_m3, pp_m3,                     & !OUT
                     uu_0  ,vv_0  ,ww_0  ,                           & !IN old 
                     uu_m2, vv_m2, ww_m2, pp_m2, volsum_m2, den_m2,  & !IN m
                     NXm,NYm,NZm,wall,pp_solve_3)                     
 !============================================================================                     
  CALL HEIGHTVOL_colocated(volsum_m3,voll,uu_m2,vv_m2,volsum_0,NXm,NYm,NZm)
  IF (height_smoothing==1) CALL HEIGHT_SMOOTH(volsum_m3,NXm,NYm)
 !============================================================================ 
  CALL TRANSPORT_CIP_colocated(ff_m1,  ffX_m1,  ffY_m1,  ffZ_m1,        & ! OUT !Recycle ff_m1
                         ff_0  , ffX_0  , ffY_0  , ffZ_0  ,       & ! IN
                        uu_m2,  vv_m2,  ww_m2,NXm,NYm,NZm)          ! IN m
  CALL DENSITYSAL_R(den_m3,ff_m1,NXm,NYm,NZm)
 !============================================================================ 
                     
 !K 4
 dt=dt_o
 !============================================================================                     
  CALL VELOCITYONLY_colocated(uu_m4, vv_m4, ww_m4, pp_m4,                     & !OUT
                     uu_0  ,vv_0  ,ww_0  ,                           & !IN old 
                     uu_m3, vv_m3, ww_m3, pp_m3, volsum_m3, den_m3,  & !IN m
                     NXm,NYm,NZm,wall,pp_solve_4)                     
 !============================================================================                     
  CALL HEIGHTVOL_colocated(volsum_m4,voll,uu_m3,vv_m3,volsum_0,NXm,NYm,NZm)
  IF (height_smoothing==1) CALL HEIGHT_SMOOTH(volsum_m4,NXm,NYm)
 !============================================================================ 
  CALL TRANSPORT_CIP_colocated(ff_m1,  ffX_m1,  ffY_m1,  ffZ_m1,        & ! OUT !Recycle ff_m1
                         ff_0  , ffX_0  , ffY_0  , ffZ_0  ,       & ! IN
                        uu_m3,  vv_m3,  ww_m3,NXm,NYm,NZm)          ! IN m
  CALL DENSITYSAL_R(den_m4,ff_m1,NXm,NYm,NZm)
 !============================================================================
 
 !K 5 - final
dt=dt_o
 !============================================================================ 
uu_m5     =1./6*(uu_m1     + 2*uu_m2     + 2*uu_m3     + uu_m4)
vv_m5     =1./6*(vv_m1     + 2*vv_m2     + 2*vv_m3     + vv_m4)
ww_m5     =1./6*(ww_m1     + 2*ww_m2     + 2*ww_m3     + ww_m4)
pp_m5     =1./6*(pp_m1     + 2*pp_m2     + 2*pp_m3     + pp_m4)
den_m5    =1./6*(den_m1    + 2*den_m2    + 2*den_m3    + den_m4)
volsum_m5 =1./6*(volsum_m1 + 2*volsum_m2 + 2*volsum_m3 + volsum_m4)
 !============================================================================                     
  CALL VELOCITYONLY_colocated(uu    ,vv    ,ww,    pp,                        & !OUT
                     uu_0  ,vv_0  ,ww_0  ,                           & !IN old 
                     uu_m5, vv_m5, ww_m5, pp_m5, volsum_m5, den_m5,  & !IN m
                     NXm,NYm,NZm,wall,1)                     
 !============================================================================                     
  CALL HEIGHTVOL_colocated(volsum,voll,uu_m5,vv_m5,volsum_0,NXm,NYm,NZm)
  IF (height_smoothing==1) CALL HEIGHT_SMOOTH(volsum,NXm,NYm)
 !============================================================================ 
  CALL TRANSPORT_CIP_colocated(ff    ,  ffX   ,  ffY   ,  ffZ   ,        & ! OUT
                         ff_0  , ffX_0  , ffY_0  , ffZ_0  ,        & ! IN
                        uu_m5,  vv_m5,  ww_m5,NXm,NYm,NZm)             ! IN m
  CALL DENSITYSAL_R(den,ff,NXm,NYm,NZm)
 !============================================================================
 
  CASE DEFAULT
   WRITE(*,*) ' Time_Order Error! '; STOP   
!---------------------------------------------------------)
 ENDSELECT
dt=dt_o


ENDSUBROUTINE TimeIntegration_colocated