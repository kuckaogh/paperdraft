PROGRAM RecedingBoundary

USE COMDAT
INTEGER ii
INTEGER(4)                                  :: standwav_ana
INTEGER(4)                                  :: NX,NY,NZ
INTEGER(4)                                  :: NXm,NYm,NZm


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
        !NXm=NX/2
        
! CALL EX_R(NX,NY,NZ,NXm,NYm,NZm)
 CALL EX_colocated(NX,NY,NZ,NXm,NYm,NZm)
 
IF (standwav_ana>0) THEN
OPEN( unit = 153, file = 'out-analyt.plt'  , status = 'unknown', form = 'formatted' )
OPEN( unit = 156, file = 'out-analyt-surf.plt'  , status = 'unknown', form = 'formatted' )
 DO ii=0,NT,ipltn
 CALL OUTPUT_ANALYTICAL_INIT(NX,NY,NZ,dx,dz,dt,ii)
 ENDDO
ENDIF


ENDPROGRAM