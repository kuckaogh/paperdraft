SUBROUTINE EX_colocated(NX,NY,NZ,NXm,NYm,NZm)

USE COMDAT
INTEGER(4)                      ,INTENT(IN)    :: NX,NY,NZ,NXm,NYm,NZm
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1)    :: ff,ffX,ffY,ffZ,ff_ori
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1)    :: ff_0,ffX_0,ffY_0,ffZ_0
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1)    :: den,den_0  
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1)    :: uu_0,uu
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1)    :: vv_0,vv
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1)    :: ww_0,ww
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1)    :: pp_0,pp
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1)           :: volsum,volsum_0 
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1)    :: voll

REAL(8)    ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1)  :: ff_r  ,ffX_r  ,ffY_r  ,ffZ_r
REAL(8)    ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1)  :: ff_r_0,ffX_r_0,ffY_r_0,ffZ_r_0
REAL(8)    ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1)  :: den_r,den_r_0  
REAL(8)    ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1)  :: uu_r_0,uu_r
REAL(8)    ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1)  :: vv_r_0,vv_r
REAL(8)    ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1)  :: ww_r_0,ww_r
REAL(8)    ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1)  :: pp_r_0,pp_r
REAL(8)    ,DIMENSION(0:NXm+1,0:NYm+1)          :: volsum_r,volsum_r_0 
REAL(8)    ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1)  :: voll_r

REAL(8)                                         :: totalsalin_ori=0, totalsalin=0
INTEGER(4)                                      :: time_i=0

!REAL(8)    ,DIMENSION(:,:,:), Allocatable    :: test

ff=0;ffX=0;ffY=0;ffZ=0
ff_0=0;ffX_0=0;ffY_0=0;ffZ_0=0
den=0;den_0=0 
uu_0=0;vv_0=0;ww_0=0
uu=0;vv=0;ww=0 
pp=0;pp_0=0
volsum=0;volsum_0=0;voll=0
uu_Max_CL_T = 0; vv_Max_CL_T = 0; ww_Max_CL_T = 0; uu_Max_CL = 0; vv_Max_CL = 0; ww_Max_CL = 0
itersum=0;errsum=0

step     = 0
time_sec = 0
den      = den_w

 CONTINUE !} end implicit none

!01/24/04
KdifY=KdifX

  CALL     READINPUT
  CALL     INIT_VOL(voll,volsum,NX,NY,NZ)  !OUT
  CALL     INIT3D(ff_initial,ff_strati,ff,ffX,ffy,ffZ,NX,NY,NZ,NXS_start,NXS_end,NYS,NZS,NZS2,dx,dy,dz)
  CALL     DENSITYSAL_R(den,ff,NX,NY,NZ)
  CALL     OUTPUT_INIT_colocated(voll,volsum,uu,vv,ww,pp,ff,den,NX,NY,NZ)
  IF (ff_initial>0) THEN
    CALL  TOTAL(ff,voll,totalsalin_ori,NX,NY,NZ)  ! for conservation, one-domain only
  ENDIF
  CALL     CPU_TIME(time_begin)


!---------------------------------------------------------------------------------------
uu_r_0  (0:NXm+1,:,:) = uu  (0:NXm+1,:,:)
vv_r_0  (0:NXm+1,:,:) = vv  (0:NXm+1,:,:)     
ww_r_0  (0:NXm+1,:,:) = ww  (0:NXm+1,:,:)
ff_r_0  (0:NXm+1,:,:) = ff  (0:NXm+1,:,:) 
ffX_r_0 (0:NXm+1,:,:) = ffX (0:NXm+1,:,:)   
ffY_r_0 (0:NXm+1,:,:) = ffY (0:NXm+1,:,:) 
ffZ_r_0 (0:NXm+1,:,:) = ffZ (0:NXm+1,:,:)
den_r_0 (0:NXm+1,:,:) = den (0:NXm+1,:,:)
 pp_r_0 (0:NXm+1,:,:) = pp  (0:NXm+1,:,:)  
voll_r  (0:NXm+1,:,:) = voll(0:NXm+1,:,:)
volsum_r_0(0:NXm+1,:) = volsum(0:NXm+1,:) 
!; ff_ori = ff !OneDomain Only
!---------------------------------------------------------------------------------------

TimeStep: DO step=1,NT
  
 CALL ACC_2
 
!***************************************************************************************************

time_i = time_i + 1
time_sec = time_i * dt_o


!CALL TimeIntegration(uu_r  ,vv_r  ,ww_r,  pp_r  ,volsum_r  ,den_r  ,ff_r  ,ffX_r  ,ffY_r  ,ffZ_r,   & !OUT
!                     uu_r_0,vv_r_0,ww_r_0,pp_r_0,volsum_r_0,den_r_0,ff_r_0,ffX_r_0,ffY_r_0,ffZ_r_0, & !IN
!                     NXm ,NYm ,NZm, 3)                                            !IN

CALL TimeIntegration_colocated(uu_r  ,vv_r  ,ww_r,  pp_r  ,volsum_r  ,den_r  ,ff_r  ,ffX_r  ,ffY_r  ,ffZ_r,   & !OUT
                     uu_r_0,vv_r_0,ww_r_0,pp_r_0,volsum_r_0,den_r_0,ff_r_0,ffX_r_0,ffY_r_0,ffZ_r_0, & !IN
                     NXm ,NYm ,NZm, 3) 

!---------------------------------------------------------------------------------------
uu_r_0  = uu_r  ; vv_r_0     = vv_r     ; ww_r_0   = ww_r
ff_r_0  = ff_r  ; ffX_r_0    = ffX_r    ; ffY_r_0  = ffY_r ; ffZ_r_0  = ffZ_r
den_r_0 = den_r ; volsum_r_0 = volsum_r      !; ff_ori = ff !Onedomain only
pp_r_0  = pp_r
!---------------------------------------------------------------------------------------

!***************************************************************************************************
 CALL ASSIGNVOL(voll_r,volsum_r,NXm ,NYm ,NZm)
 CALL STABILITY_colocated(uu_r,vv_r,ww_r,NXm,NYm,NZm) 

IF (ff_initial>0)  THEN ! One Domain Only
 !CALL TOTAL(ff,voll,totalsalin,NX,NY,NZ)
 CALL TOTAL(ff_r,voll,totalsalin,NX,NY,NZ) 
 IF(conservation==1) THEN
  CALL RESET_FF(totalsalin_ori, totalsalin, ff,ff_ori,NX,NY,NZ) 
 ENDIF
ENDIF 

  CALL CPU_TIME(time_end); PRINT *, 'CPU_TIME =', time_end-time_begin; WRITE(18,*) 'CPU_TIME =', time_end-time_begin

IF(uu_Max_CL > uu_Max_CL_T) uu_Max_CL_T = uu_Max_CL
IF(vv_Max_CL > vv_Max_CL_T) vv_Max_CL_T = vv_Max_CL
IF(ww_Max_CL > ww_Max_CL_T) ww_Max_CL_T = ww_Max_CL

itersum = itersum + iSOR * ii_SOR_Pos
errsum  = errsum  + err
!------------------------------------------------------------------( 
  IF (MOD(step,1*ipltn) == 0) THEN
    !CALL OUTPUT(voll,volsum,uu,vv,ww,pp,ff,den,NXm,NYm,NZm)
    CALL OUTPUT_INIT_colocated(voll_r,volsum_r,uu_r,vv_r,ww_r,pp_r,ff_r,den_r,NXm,NYm,NZm)
    IF (ff_initial>0) write(182,'(F8.3,A)') totalsalin/totalsalin_ori*100,'%' 
    WRITE(181,1011) step*dt, time_end-time_begin, REAL(itersum)/step, REAL(errsum)/step, uu_Max_CL_T, vv_Max_CL_T, ww_Max_CL_T
  ENDIF
!-------------------------------------------------------------------)

ENDDO TimeStep

1011   format(3(F11.4,2X), E11.4,2X    ,3(F11.4,2X) )
1012   format(4(F11.4,2X), 2(E11.4,2X) ,3(F11.4,2X) )
ENDSUBROUTINE EX_colocated