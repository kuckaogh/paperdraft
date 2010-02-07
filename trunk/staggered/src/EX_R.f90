SUBROUTINE EX_R

USE COMDAT
IMPLICIT NONE !{


REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1)    :: ff,ffX,ffY,ffZ   ,ff_ori
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1)    :: ff_0,ffX_0,ffY_0,ffZ_0
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1)    :: den,den_0  
REAL(8)    ,DIMENSION(0:NX  ,0:NY+1,0:NZ+1)    :: uu_0
REAL(8)    ,DIMENSION(0:NX  ,0:NY+1,0:NZ+1)    :: uu 
REAL(8)    ,DIMENSION(0:NX+1,0:NY  ,0:NZ+1)    :: vv_0
REAL(8)    ,DIMENSION(0:NX+1,0:NY  ,0:NZ+1)    :: vv 
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ  )    :: ww_0
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ  )    :: ww 
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1)    :: pp
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1)    :: pp_0
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1)    :: voll,voll_0
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1)           :: volsum,volsum_0 




REAL(8)                                        :: totalsalin_ori=0, totalsalin=0
INTEGER :: time_i=0


REAL(8)    ,DIMENSION(0:NX  ,0:NY+1,0:NZ+1)    :: uu_m1,uu_m2,uu_m3,uu_m4,uu_m5,uu_I
REAL(8)    ,DIMENSION(0:NX+1,0:NY  ,0:NZ+1)    :: vv_m1,vv_m2,vv_m3,vv_m4,vv_m5,vv_I
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ  )    :: ww_m1,ww_m2,ww_m3,ww_m4,ww_m5,ww_I
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1)    :: pp_m1,pp_m2,pp_m3,pp_m4,pp_m5,pp_I
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1)           :: volsum_m1,volsum_m2,volsum_m3,volsum_m4,volsum_m5,volsum_I
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1)    :: den_m1,den_m2,den_m3,den_m4,den_m5,den_I
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1)    :: ff_m1,ffX_m1,ffY_m1,ffZ_m1

!
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1)    :: ff_m2,ffX_m2,ffY_m2,ffZ_m2
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1)    :: ff_m3,ffX_m3,ffY_m3,ffZ_m3
!REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1)    :: ff_m4,ffX_m4,ffY_m4,ffZ_m4
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1)    :: ff_I,ffX_I,ffY_I,ffZ_I




ff=0;ffX=0;ffY=0;ffZ=0
ff_0=0;ffX_0=0;ffY_0=0;ffZ_0=0


den=0;den_0=0 

uu_0=0
uu=0 
vv_0=0
vv=0 
ww_0=0
ww=0 
pp=0
voll=0;voll_0=0 
volsum=0;volsum_0=0   


uu_Max_CL_T = 0; vv_Max_CL_T = 0; ww_Max_CL_T = 0; uu_Max_CL = 0; vv_Max_CL = 0; ww_Max_CL = 0
itersum=0
errsum=0


step     = 0
time_sec = 0
den      = den_w

 CONTINUE !} end implicit none

!01/24/04
KdifY=KdifX

  CALL     READINPUT
   dt=dt_o

  CALL     INIT_VOL(voll,volsum) !OUT
  CALL     INIT3D(ff_initial,ff_strati,ff,ffX,ffy,ffZ,NX,NY,NZ,NXS_start,NXS_end,NYS,NZS,NZS2,dx,dy,dz)
  CALL     DENSITYSAL_R(den,ff)
    
  CALL     OUTPUT_INIT(voll,volsum,uu,vv,ww,pp,ff,den)
  IF (ff_initial>0) CALL  TOTAL(ff,voll,totalsalin_ori) ! for conservation, one-domain only
  CALL     CPU_TIME(time_begin)

                   
dt  = dt_o 

LOOP_TimeStepping_OneDomain: DO step=1,NT
  
 CALL ACC_2
 
!---------------------------------------------------------------------------------------
uu_0        = uu         ;vv_0        = vv       ;ww_0        = ww
ff_0        = ff         ;ffX_0       = ffX      ;ffY_0       = ffY ;ffZ_0       = ffZ
den_0       = den
volsum_0    = volsum     ;voll_0      = voll
ff_ori = ff
!---------------------------------------------------------------------------------------


!***************************************************************************************************


LOOP_SplitTimeStep_OneDomain: DO step_1=1,1

time_i = time_i + 1
time_sec = time_i * dt_o


 SELECTCASE(TIME_ORDER)
 CASE(1) ! Euler
!---------------------------------------------------------(
 CALL VELOCITY(uu,   vv,   ww,   pp,   volsum,   den,   ff,   ffX,   ffY,   ffZ,     & !OUT
               uu_0, vv_0, ww_0,       volsum_0,        ff_0, ffX_0, ffY_0, ffZ_0,   & !IN old
               uu_0, vv_0, ww_0, pp_0, volsum_0, den_0,                              & !IN m
               1)                                      
!---------------------------------------------------------)

 CASE(2) ! Predictor-Corrector 
!---------------------------------------------------------(
 dt=0.5*dt_o
 CALL VELOCITY(uu_m1,vv_m1,ww_m1,pp_m1,volsum_m1,den_m1,ff_m1,ffX_m1,ffY_m1,ffZ_m1,  & !OUT
               uu_0, vv_0, ww_0,       volsum_0,        ff_0, ffX_0, ffY_0, ffZ_0,   & !IN old
               uu_0, vv_0, ww_0, pp_0, volsum_0, den_0,                              & !IN m
               pp_solve_1)                                                          !IN
             
 dt=dt_o
 CALL VELOCITY(uu,   vv,   ww,   pp,   volsum,   den,   ff,   ffX,   ffY,   ffZ,     & !OUT
               uu_0, vv_0, ww_0,       volsum_0,        ff_0, ffX_0, ffY_0, ffZ_0,   & !IN old
               uu_m1,vv_m1,ww_m1,pp_m1,volsum_m1,den_m1,                              & !IN m
               1)                                                          !IN
!---------------------------------------------------------)

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
 ! CALL ASSIGNVOL_R(voll,topcell,volsum_0)
  CALL VELOCITYONLY(uu_m1 ,vv_m1 ,ww_m1 ,pp_m1,                     & !OUT
                     uu_0  ,vv_0  ,ww_0  ,                           & !IN old 
                     uu_0  ,vv_0  ,ww_0  ,pp_0  ,volsum_0  ,den_0  , & !IN m
                     pp_solve_1)
 !============================================================================                     
  CALL ASSIGNVOL_R(voll,volsum_0)
  CALL HEIGHTONLY_R(volsum_m1,voll,uu_0  ,vv_0  ,volsum_0  )
  IF (height_smoothing==1) CALL HEIGHT_SMOOTH(volsum_m1)
 !============================================================================ 
  CALL BDC3D_R(ff_0  ,ffX_0  ,ffY_0  ,ffZ_0  )
  CALL TRANSPORT_CIP_R3D(ff_m1,  ffX_m1,  ffY_m1,  ffZ_m1,        & ! OUT
                         ff_0  , ffX_0  , ffY_0  , ffZ_0  ,       & ! IN
                        uu_0  ,  vv_0  ,  ww_0  )          ! IN m
  CALL DENSITYSAL_R(den_m1,ff_m1)
 !============================================================================ 
 !U 2
dt=0.5*dt_o 
 !============================================================================                    
 ! CALL ASSIGNVOL_R(voll,topcell,volsum_m1)
  CALL VELOCITYONLY(uu_m2, vv_m2, ww_m2, pp_m2,                     & !OUT
                     uu_m1 ,vv_m1 ,ww_m1 ,                           & !IN old
                     uu_m1 ,vv_m1 ,ww_m1 ,pp_m1, volsum_m1, den_m1,  & !IN m
                     pp_solve_2)
 !============================================================================                     
  CALL ASSIGNVOL_R(voll,volsum_m1)
  CALL HEIGHTONLY_R(volsum_m2,voll,uu_m1,vv_m1,volsum_m1  )
  IF (height_smoothing==1) CALL HEIGHT_SMOOTH(volsum_m2)
 !============================================================================ 
 ! CALL BDC3D_R(ff_0  ,ffX_0  ,ffY_0  ,ffZ_0  )
  CALL TRANSPORT_CIP_R3D(ff_m2,  ffX_m2,  ffY_m2,  ffZ_m2,        & ! OUT 
                         ff_m1,  ffX_m1,  ffY_m1,  ffZ_m1,        & ! IN
                        uu_m1,  vv_m1,  ww_m1)          ! IN m
  CALL DENSITYSAL_R(den_m2,ff_m2)
 !============================================================================ 
 !U 3
dt=1./6*dt_o
 uu_I     = 2./3*uu_0     + 1./3*uu_m2
 vv_I     = 2./3*vv_0     + 1./3*vv_m2
 ww_I     = 2./3*ww_0     + 1./3*ww_m2
 pp_I     = pp_m2   
 den_I    = 2./3*den_0    + 1./3*den_m2   
 volsum_I = 2./3*volsum_0 + 1./3*volsum_m2
 ff_I     = 2./3*ff_0     + 1./3*ff_m2      
 ffX_I    = 2./3*ffX_0    + 1./3*ffX_m2
 ffY_I    = 2./3*ffY_0    + 1./3*ffY_m2
 ffZ_I    = 2./3*ffZ_0    + 1./3*ffZ_m2
 !============================================================================                     
 ! CALL ASSIGNVOL_R(voll,topcell,volsum_I)
  CALL VELOCITYONLY(uu_m3, vv_m3, ww_m3, pp_m3,                     & !OUT
                     uu_I  ,vv_I  ,ww_I,                            & !IN old 
                     uu_m2, vv_m2, ww_m2, pp_m2, volsum_m2, den_m2,  & !IN m
                     pp_solve_3)                     
 !============================================================================                     
  CALL ASSIGNVOL_R(voll,volsum_I  )
  CALL HEIGHTONLY_R(volsum_m3,voll,uu_m2,vv_m2,volsum_I  )
  IF (height_smoothing==1) CALL HEIGHT_SMOOTH(volsum_m3)
 !============================================================================ 
  CALL BDC3D_R(ff_I  ,ffX_I  ,ffY_I  ,ffZ_I  )
  CALL TRANSPORT_CIP_R3D(ff_m3,  ffX_m3,  ffY_m3,  ffZ_m3,        & ! OUT 
                         ff_I  , ffX_I  , ffY_I  , ffZ_I  ,       & ! IN
                        uu_m2,  vv_m2,  ww_m2)          ! IN m
  CALL DENSITYSAL_R(den_m3,ff_m3)
 !============================================================================ 
 !U 4 - final
 dt=0.5*dt_o
 !============================================================================                     
 ! CALL ASSIGNVOL_R(voll,topcell,volsum_m3  )
  CALL VELOCITYONLY(uu  ,   vv,     ww     , pp,                     & !OUT
                    uu_m3,  vv_m3,  ww_m3  ,                           & !IN old 
                    uu_m3,  vv_m3,  ww_m3, pp_m3, volsum_m3, den_m3,  & !IN m
                    1)                     
 !============================================================================                     
  CALL ASSIGNVOL_R(voll,volsum_m3  )
  CALL HEIGHTONLY_R(volsum,voll,uu_m3,vv_m3,volsum_m3  )
  IF (height_smoothing==1) CALL HEIGHT_SMOOTH(volsum)
 !============================================================================ 
  CALL BDC3D_R(ff_m3  ,ffX_m3  ,ffY_m3  ,ffZ_m3  )
  CALL TRANSPORT_CIP_R3D(ff,    ffX,    ffY,    ffZ,           & ! OUT !Recycle ff_m1
                         ff_m3, ffX_m3, ffY_m3, ffZ_m3,        & ! IN
                        uu_m3,  vv_m3,  ww_m3)           ! IN m
  CALL DENSITYSAL_R(den,ff)
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
  CALL VELOCITYONLY(uu_m2, vv_m2, ww_m2, pp_m2,                     & !OUT
                     uu_0  ,vv_0  ,ww_0  ,                           & !IN old
                     uu_m1, vv_m1, ww_m1, pp_m1, volsum_m1, den_m1,  & !IN m
                     pp_solve_2)
 !============================================================================                     
  CALL ASSIGNVOL_R(voll,volsum_0  )
  CALL HEIGHTONLY_R(volsum_m2,voll,uu_m1,vv_m1,volsum_0  )
  IF (height_smoothing==1) CALL HEIGHT_SMOOTH(volsum_m2)
 !============================================================================ 
  CALL BDC3D_R(ff_0  ,ffX_0  ,ffY_0  ,ffZ_0  )
  CALL TRANSPORT_CIP_R3D(ff_m1,  ffX_m1,  ffY_m1,  ffZ_m1,        & ! OUT !Recycle ff_m1
                         ff_0  , ffX_0  , ffY_0  , ffZ_0  ,       & ! IN
                        uu_m1,  vv_m1,  ww_m1)          ! IN m
  CALL DENSITYSAL_R(den_m2,ff_m1)
 !============================================================================ 
                     
 !K 3
 dt=0.5*dt_o
 !============================================================================                     
  CALL VELOCITYONLY(uu_m3, vv_m3, ww_m3, pp_m3,                     & !OUT
                     uu_0  ,vv_0  ,ww_0  ,                           & !IN old 
                     uu_m2, vv_m2, ww_m2, pp_m2, volsum_m2, den_m2,  & !IN m
                     pp_solve_3)                     
 !============================================================================                     
 ! CALL ASSIGNVOL_R(voll,topcell,volsum_0  )
  CALL HEIGHTONLY_R(volsum_m3,voll,uu_m2,vv_m2,volsum_0  )
  IF (height_smoothing==1) CALL HEIGHT_SMOOTH(volsum_m3)
 !============================================================================ 
 ! CALL BDC3D_R(ff_0  ,ffX_0  ,ffY_0  ,ffZ_0  )
  CALL TRANSPORT_CIP_R3D(ff_m1,  ffX_m1,  ffY_m1,  ffZ_m1,        & ! OUT !Recycle ff_m1
                         ff_0  , ffX_0  , ffY_0  , ffZ_0  ,       & ! IN
                        uu_m2,  vv_m2,  ww_m2)          ! IN m
  CALL DENSITYSAL_R(den_m3,ff_m1)
 !============================================================================ 
                     
 !K 4
 dt=dt_o
 !============================================================================                     
  CALL VELOCITYONLY(uu_m4, vv_m4, ww_m4, pp_m4,                     & !OUT
                     uu_0  ,vv_0  ,ww_0  ,                           & !IN old 
                     uu_m3, vv_m3, ww_m3, pp_m3, volsum_m3, den_m3,  & !IN m
                     pp_solve_4)                     
 !============================================================================                     
 ! CALL ASSIGNVOL_R(voll,topcell,volsum_0  )
  CALL HEIGHTONLY_R(volsum_m4,voll,uu_m3,vv_m3,volsum_0  )
  IF (height_smoothing==1) CALL HEIGHT_SMOOTH(volsum_m4)
 !============================================================================ 
 ! CALL BDC3D_R(ff_0  ,ffX_0  ,ffY_0  ,ffZ_0  )
  CALL TRANSPORT_CIP_R3D(ff_m1,  ffX_m1,  ffY_m1,  ffZ_m1,        & ! OUT !Recycle ff_m1
                         ff_0  , ffX_0  , ffY_0  , ffZ_0  ,       & ! IN
                        uu_m3,  vv_m3,  ww_m3)          ! IN m
  CALL DENSITYSAL_R(den_m4,ff_m1)
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
  CALL VELOCITYONLY(uu    ,vv    ,ww,    pp,                        & !OUT
                     uu_0  ,vv_0  ,ww_0  ,                           & !IN old 
                     uu_m5, vv_m5, ww_m5, pp_m5, volsum_m5, den_m5,  & !IN m
                     1)                     
 !============================================================================                     
 ! CALL ASSIGNVOL_R(voll,topcell,volsum_0  )
  CALL HEIGHTONLY_R(volsum,voll,uu_m5,vv_m5,volsum_0  )
  IF (height_smoothing==1) CALL HEIGHT_SMOOTH(volsum)
 !============================================================================ 
 ! CALL BDC3D_R(ff_0  ,ffX_0  ,ffY_0  ,ffZ_0  )
  CALL TRANSPORT_CIP_R3D(ff    ,  ffX   ,  ffY   ,  ffZ   ,        & ! OUT
                         ff_0  , ffX_0  , ffY_0  , ffZ_0  ,        & ! IN
                        uu_m5,  vv_m5,  ww_m5)              ! IN m
  CALL DENSITYSAL_R(den,ff)
 !============================================================================
 
  CASE DEFAULT
   WRITE(*,*) ' Time_Order Error! '; STOP   
!---------------------------------------------------------)
 ENDSELECT
dt=dt_o


!---------------------------------------------------------(
uu_0 = uu ; vv_0 = vv ; ww_0 = ww
ff_0 = ff ; ffX_0 = ffX ; ffY_0 = ffY ; ffZ_0 = ffZ
den_0 = den ; volsum_0 = volsum 
pp_0 = pp
!---------------------------------------------------------)

!----------------------------------------------------------------------------------------)
ENDDO LOOP_SplitTimeStep_OneDomain

!***************************************************************************************************


 CALL ASSIGNVOL(voll,volsum)

 CALL STABILITY(uu,vv,ww) 

IF (ff_initial>0)  THEN
 CALL TOTAL(ff,voll,totalsalin)
 IF(conservation==1) THEN
  CALL RESET_FF(totalsalin_ori, totalsalin, ff,ff_ori)
 ENDIF
ENDIF 

  CALL CPU_TIME(time_end); PRINT *, 'CPU_TIME =', time_end-time_begin; WRITE(18,*) 'CPU_TIME =', time_end-time_begin

IF(uu_Max_CL > uu_Max_CL_T) THEN
  uu_Max_CL_T = uu_Max_CL
ENDIF 

IF(vv_Max_CL > vv_Max_CL_T) THEN
  vv_Max_CL_T = vv_Max_CL
ENDIF 

IF(ww_Max_CL > ww_Max_CL_T) THEN
  ww_Max_CL_T = ww_Max_CL
ENDIF

itersum = itersum + iSOR * ii_SOR_Pos
errsum  = errsum  + err
!------------------------------------------------------------------( 
  IF (MOD(step,1*ipltn) == 0) THEN

    CALL OUTPUT(voll,volsum,uu,vv,ww,pp,ff,den)
    IF (ff_initial>0) write(182,'(F8.3,A)') totalsalin/totalsalin_ori*100,'%' 

    WRITE(181,1011) step*dt, time_end-time_begin, REAL(itersum)/step, REAL(errsum)/step, uu_Max_CL_T, vv_Max_CL_T, ww_Max_CL_T

  ENDIF
!-------------------------------------------------------------------)

ENDDO LOOP_TimeStepping_OneDomain

1011   format(3(F11.4,2X), E11.4,2X    ,3(F11.4,2X) )
1012   format(4(F11.4,2X), 2(E11.4,2X) ,3(F11.4,2X) )
ENDSUBROUTINE EX_R