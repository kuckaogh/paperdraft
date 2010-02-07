SUBROUTINE READINPUT

USE COMDAT


READ(14,*) dummy
READ(14,*) NXS_start,NXS_end, NYS, NZS , NZS2
READ(14,*) dummy

READ(14,*) dx, dy, visc, gravity, den_w
READ(14,*) dummy
READ(14,*) dt_o, NT, ipltn, restart, irestartn
READ(14,*) dummy 
READ(14,*) ff_initial,ff_strati
READ(14,*) dummy
READ(14,*) Am0, h0
READ(14,*) dummy
READ(14,*) ii_SOR, ii_SOR_Min, ii_SOR_Pos, beta,  errmax
READ(14,*) dummy
READ(14,*) uuwall, Height_field, height_smoothing
READ(14,*) dummy 
READ(14,*) KdifX,KdifZ, transp_x, transp_y, transp_z
READ(14,*) dummy 
READ(14,*) WBaC, VBC(1), VBC(2), VBC(3), VBC(4), VBC(5)
READ(14,*) dummy
READ(14,*) Bossi,i_dfdxN,i_dfdx_vN
READ(14,*) dummy
READ(14,*) ANH, CIP2D
READ(14,*) dummy
READ(14,*) ISOR_ErrMax
READ(14,*) dummy
READ(14,*) conservation, conservation_power
READ(14,*) dummy
READ(14,*) TIME_ORDER, pp_solve_1, pp_solve_2, pp_solve_3, pp_solve_4
READ(14,*) dummy
READ(14,*) OPENZ, SurfaceCorrection
READ(14,*) dummy
READ(14,*) MultiLevel,i_Multi1,i_Multi2,i_Multi3,i_Multi4,i_Multi5
READ(14,*) dummy
READ(14,*) errLevel1,errLevel2,errLevel3,errLevel4
!*********************************************************************************


!Initial check

IF((ANH==0).and.(WBAC/=1)) THEN
  WRITE(*,*) 'check ANH and WBAC'
  STOP
ELSEIF(ii_SOR < ii_SOR_Min) THEN
  WRITE(*,*) 'check ii_SOR_Min'
  STOP
ENDIF   



ENDSUBROUTINE READINPUT