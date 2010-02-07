SUBROUTINE RESTART_READ(ff,ffX,ffY,ffZ  &
& ,ff_0,ffX_0,ffY_0,ffZ_0                &
& ,ff_OUT,ffX_OUT,ffY_OUT,ffZ_OUT        &
& ,ff_IN,ffX_IN,ffY_IN,ffZ_IN            &
& ,den,den_0,den_IN                      &
& ,den_OUT                               &
& ,uu_0                                  &
& ,uu,uu_old,uu_OUT,uu_IN                &
& ,vv_0                                  &
& ,vv,vv_old,vv_OUT,vv_IN                &
& ,ww_0                                  &
& ,ww,ww_old,ww_OUT,ww_IN                &
& ,pp,pp_part1,pp_part2                  &
& ,voll,voll_0,voll_IN                   &
& ,pp_OUT                                &
& ,volsum,volsum_0,volsum_IN             &
& ,volsum_OUT                            &
& ,w_kinetic_OUT                         &
& ,w_kinetic,w_kinetic_0,w_kinetic_IN    &
& ,topcell,topcell_0,topcell_IN          &
& ,domain                                &
& ,totalsalin_ori, totalsalin)

USE COMDAT

IMPLICIT NONE !{
!INTEGER :: NXBdum=0,NXEdum=0
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1)    :: ff,ffX,ffY,ffZ
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1)    :: ff_0,ffX_0,ffY_0,ffZ_0
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1)    :: ff_OUT,ffX_OUT,ffY_OUT,ffZ_OUT
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1)    :: ff_IN,ffX_IN,ffY_IN,ffZ_IN
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1)    :: den,den_0,den_IN 
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1)    :: den_OUT
REAL(8)    ,DIMENSION(0:NX  ,0:NY+1,0:NZ+1)    :: uu_0
REAL(8)    ,DIMENSION(0:NX  ,0:NY+1,0:NZ+1)    :: uu,uu_old,uu_OUT,uu_IN !uu_R
REAL(8)    ,DIMENSION(0:NX+1,0:NY  ,0:NZ+1)    :: vv_0
REAL(8)    ,DIMENSION(0:NX+1,0:NY  ,0:NZ+1)    :: vv,vv_old,vv_OUT,vv_IN !vv_R
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ  )    :: ww_0
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ  )    :: ww,ww_old,ww_OUT,ww_IN !ww_R
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1)    :: pp,pp_part1,pp_part2 
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1)    :: voll,voll_0,voll_IN 
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1)    :: pp_OUT 
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1)           :: volsum,volsum_0,volsum_IN
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1)           :: volsum_OUT 
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1)           :: w_kinetic_OUT
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1)           :: w_kinetic,w_kinetic_0,w_kinetic_IN 
INTEGER(4) ,DIMENSION(0:NX+1,0:NY+1)           :: topcell,topcell_0,topcell_IN 
INTEGER(2)                                     :: domain
 
REAL(8)                                        :: totalsalin_ori, totalsalin

READ(16,*)    
READ(16,*)     step
READ(16,901)  ff,ffX,ffY,ffZ
READ(16,901)   ff_0,ffX_0,ffY_0,ffZ_0
READ(16,901)   ff_OUT,ffX_OUT,ffY_OUT,ffZ_OUT
READ(16,901)   ff_IN,ffX_IN,ffY_IN,ffZ_IN
READ(16,901)   den,den_0,den_IN 
READ(16,901)   den_OUT
READ(16,901)   uu_0
READ(16,901)  uu,uu_old,uu_OUT,uu_IN 
READ(16,901)  vv_0
READ(16,901)  vv,vv_old,vv_OUT,vv_IN 
READ(16,901)  ww_0
READ(16,901)  ww,ww_old,ww_OUT,ww_IN 
READ(16,901)  pp,pp_part1,pp_part2 
READ(16,901)  voll,voll_0,voll_IN 
READ(16,901)  pp_OUT 
READ(16,901)  volsum,volsum_0,volsum_IN
READ(16,901)  volsum_OUT 
READ(16,901)  w_kinetic_OUT
READ(16,901)  w_kinetic,w_kinetic_0,w_kinetic_IN 
READ(16,*)  topcell,topcell_0,topcell_IN 
READ(16,*)  domain
READ(16,901)  totalsalin_ori, totalsalin, time_sec, dt_o


901   format(E17.9)

ENDSUBROUTINE