PROGRAM RecedingBoundary

USE COMDAT
INTEGER ii
INTEGER(4)                                  :: standwav_ana

OPEN( unit = 14, FILE = 'input.dat')
OPEN( unit = 15, file = 'out.plt'  , status = 'unknown', form = 'formatted' )
OPEN( unit = 151, file = 'out_D1.plt'  , status = 'unknown', form = 'formatted' )
OPEN( unit = 152, file = 'out_D2.plt'  , status = 'unknown', form = 'formatted' )
OPEN( unit = 16, file = 'restart_read.txt' )
OPEN( unit = 17, file = 'restart_write.txt' )
OPEN( unit = 18, file = 'log.txt' )
OPEN( unit = 181, file = 'log-SOR-CFL.txt' )
OPEN( unit = 182, file = 'log-total-ff.txt' )

READ(14,*) dummy
READ(14,*) NX,NY, NZ, dz,  standwav_ana
WRITE(181,'(A)') '    Time   CPU   Ave.Iters(2)  err(2)   CFL-U   CFL-V   CFL-W   '



        NXm=NX; NYm=NY; NZm=NZ
        NXm2=FLOOR(REAL(NXm)/2); NYm2=MAX(FLOOR(REAL(NYm)/2),1); NZm2=FLOOR(REAL(NZm)/2)
        NXm3=FLOOR(REAL(NXm2)/2);NYm3=MAX(FLOOR(REAL(NYm2)/2),1);NZm3=FLOOR(REAL(NZm2)/2)
        NXm4=FLOOR(REAL(NXm3)/2);NYm4=MAX(FLOOR(REAL(NYm3)/2),1);NZm4=FLOOR(REAL(NZm3)/2)
        NXm5=FLOOR(REAL(NXm4)/2);NYm5=MAX(FLOOR(REAL(NYm4)/2),1);NZm5=FLOOR(REAL(NZm4)/2)
 CALL EX_R

IF (standwav_ana>0) THEN
OPEN( unit = 153, file = 'out-analyt.plt'  , status = 'unknown', form = 'formatted' )
OPEN( unit = 156, file = 'out-analyt-surf.plt'  , status = 'unknown', form = 'formatted' )
 DO ii=0,NT,ipltn
 CALL OUTPUT_ANALYTICAL_INIT(NX,NY,NZ,dx,dz,dt,ii)
 ENDDO
ENDIF


ENDPROGRAM