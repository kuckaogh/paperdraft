SUBROUTINE STABILITY(uu,vv,ww,NXm,NYm,NZm) 

USE COMDAT
INTEGER(4)   ,INTENT(IN)           :: NXm,NYm,NZm
REAL(8)    ,DIMENSION(0:NXm  ,0:NYm+1,0:NZm+1) ,INTENT(IN)    :: uu
REAL(8)    ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm  ) ,INTENT(IN)    :: ww
REAL(8)    ,DIMENSION(0:NXm+1,0:NYm  ,0:NZm+1) ,INTENT(IN)    :: vv

uu_Max_CL = MAXVAL(ABS(uu))*dt_o/dx
vv_Max_CL = MAXVAL(ABS(vv))*dt_o/dy
ww_Max_CL = MAXVAL(ABS(ww))*dt_o/dz

Write(18,*) '----------------------------------------------'
Write(18,'(" Max UU Corant# =         ",E11.4)') uu_Max_CL
Write(18,'(" Max VV Corant# =         ",E11.4)') vv_Max_CL
Write(18,'(" Max WW Corant# =         ",E11.4)') ww_Max_CL
Write(18,'(" Max (6D*dt/Min(dx,dy,dz)^2) =       ",E11.4)') 6*Visc*dt_o/MIN(dx,dy,dz)**2
Write(18,*) '----------------------------------------------'  
Write(*,*) '----------------------------------------------'
Write(*,'(" Max UU Corant# =         ",E11.4)') uu_Max_CL
Write(*,'(" Max VV Corant# =         ",E11.4)') vv_Max_CL
Write(*,'(" Max WW Corant# =         ",E11.4)') ww_Max_CL
Write(*,'(" Max (6D*dt/Min(dx,dy,dz)^2) =       ",E11.4)') 6*Visc*dt_o/MIN(dx,dy,dz)**2
Write(*,*) '----------------------------------------------' 

ENDSUBROUTINE