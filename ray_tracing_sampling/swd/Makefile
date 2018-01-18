include ../Makefile.compiler
## Mac OSX:
#FC = mpif90 -fc=pgfortran -C
#FC = mpif90 -fc=pgfortran -fast -O4 -fastsse
#FC = mpif90 -fc=pgfortran -fast -O4 -fastsse -Mfprelaxed -Mipa=fast,inline
## Kispiox:
#FC = mpif90 -f90=/srv/pgi/linux86-64/9.0/bin/pgfortran -tp=shanghai-64 -C
#FC = mpif90 -f90=/srv/pgi/linux86-64/9.0/bin/pgfortran -tp=shanghai-64 -fast -O4 -fastsse -Mfprelaxed -Mipa=fast,inline

FO=libswd.a

MODS = 
OBJS = dispersion.o lovmrx.o raymrx.o lovdsp.o raydsp.o four1.o
                  
$(FO): ${MODS} ${OBJS}
	ar urv libswd.a ${MODS} ${OBJS}
	cp *.a ../obj/

dispersion.o: dispersion.f90; $(FCSWD) -c dispersion.f90
lovmrx.o: lovmrx.f; $(FCSWD) -c lovmrx.f
raymrx.o: raymrx.f; $(FCSWD) -c raymrx.f
lovdsp.o: lovdsp.f; $(FCSWD) -c lovdsp.f
raydsp.o: raydsp.f; $(FCSWD) -c raydsp.f
four1.o: four1.f; $(FCSWD) -c four1.f

clean:	
	/bin/rm -f *.o *.mod *.l core *.a
