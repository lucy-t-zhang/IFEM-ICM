#  Makefile 
F90 = /usr/bin/gfortran

OBJ1= sboc.o sbpress.o scalfp.o scalfu.o scalkpp.o scalkup.o scalkuu.o 
OBJ2= scauchy.o smaterj.o spiola.o spress.o sstif.o sstrain.o sstress.o
OBJ3= stang.o stoxc.o eleinit.o element.o surcon.o print9.o acca.o
OBJ4= acc.o assbf.o asmpr.o bdpd.o bload.o jacob.o print.o 
OBJ5= fluid.o init.o inter.o rebound.o readinit.o gauss.o
OBJ6= load.o main.o mater.o nodalf.o press.o inputdat.o belement.o 
OBJ7= antiskew.o skew.o strut.o timefun.o surfree.o eleleng.o  
OBJ8= nonass.o energy.o concont.o convel.o force.o fluida.o
OBJ9= surfluid.o surele.o suracc.o smooth.o conraw.o conset.o
OBJ10= sur1d.o surdif.o delta.o map.o ghost.o stanga.o sstifa.o

run: $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) $(OBJ5) $(OBJ6) $(OBJ7) $(OBJ8) $(OBJ9) $(OBJ10)
		$(F90) -o run  $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) $(OBJ5) $(OBJ6) $(OBJ7) $(OBJ8) $(OBJ9) $(OBJ10)


main.o: main.f90

		$(F90) -c main.f90

surfree.o: surfree.f90

		$(F90) -c surfree.f90

strut.o: strut.f90

		$(F90) -c strut.f90

timefun.o: timefun.f90

		$(F90) -c timefun.f90

force.o: force.f90

		$(F90) -c force.f90

nonass.o: nonass.f90

		$(F90) -c nonass.f90

inter.o: inter.f90

		$(F90) -c inter.f90

assbf.o: assbf.f90

		$(F90) -c assbf.f90

energy.o: energy.f90

		$(F90) -c energy.f90

sstifa.o: sstifa.f90

		$(F90) -c sstifa.f90

scalkup.o: scalkup.f90

		$(F90) -c scalkup.f90

fluida.o: fluida.f90

		$(F90) -c fluida.f90

surele.o: surele.f90

		$(F90) -c surele.f90

skew.o: skew.f90

		$(F90) -c skew.f90

rebound.o: rebound.f90

		$(F90) -c rebound.f90

delta.o: delta.f90

		$(F90) -c delta.f90

surfluid.o: surfluid.f90

		$(F90) -c surfluid.f90

stoxc.o: stoxc.f90

		$(F90) -c stoxc.f90

nodalf.o: nodalf.f90

		$(F90) -c nodalf.f90

inputdat.o: inputdat.f90

		$(F90) -c inputdat.f90

asmpr.o: asmpr.f90

		$(F90) -c asmpr.f90

sstif.o: sstif.f90

		$(F90) -c sstif.f90

sstifa.o: sstifa.f90

		$(F90) -c sstifa.f90

scalkpp.o: scalkpp.f90

		$(F90) -c scalkpp.f90

convel.o: convel.f90

		$(F90) -c convel.f90

map.o: map.f90

		$(F90) -c map.f90

print9.o: print9.f90

		$(F90) -c print9.f90

fluid.o: fluid.f90

		$(F90) -c fluid.f90

fluida.o: fluida.f90

		$(F90) -c fluida.f90

surdif.o: surdif.f90

		$(F90) -c surdif.f90

readinit.o: readinit.f90

		$(F90) -c readinit.f90

concont.o: concont.f90

		$(F90) -c concont.f90

stanga.o: stanga.f90

		$(F90) -c stanga.f90

mater.o: mater.f90

		$(F90) -c mater.f90

init.o: init.f90

		$(F90) -c init.f90

antiskew.o: antiskew.f90

		$(F90) -c antiskew.f90

spress.o: spress.f90

		$(F90) -c spress.f90

scalfu.o: scalfu.f90

		$(F90) -c scalfu.f90

energy.o: energy.f90

		$(F90) -c energy.f90

surcon.o: surcon.f90

		$(F90) -c surcon.f90

print9.o: print9.f90

	$(F90) -c print9.f90

conset.o: conset.f90

	$(F90) -c conset.f90

bload.o: bload.f90

	$(F90) -c bload.f90

stang.o: stang.f90

	$(F90) -c stang.f90

map.o: map.f90

	$(F90) -c map.f90

acc.o: acc.f90

	$(F90) -c acc.f90


acca.o: acca.f90

	$(F90) -c acca.f90

spiola.o: spiola.f90

	$(F90) -c spiola.f90

scalfp.o: scalfp.f90

	$(F90) -c scalfp.f90

element.o: element.f90

	$(F90) -c element.f90

tridag.o: tridag.f90

	$(F90) -c tridag.f90

suracc.o: suracc.f90

	$(F90) -c suracc.f90

print.o: print.f90

	$(F90) -c print.f90

load.o: load.f90

	$(F90) -c load.f90

conraw.o: conraw.f90

	$(F90) -c conraw.f90

belement.o: belement.f90

	$(F90) -c belement.f90

sstress.o: sstress.f90

	$(F90) -c sstress.f90

scauchy.o: scauchy.f90

	$(F90) -c scauchy.f90

gauss.o: gauss.f90

	$(F90) -c gauss.f90

smooth.o: smooth.f90

	$(F90) -c smooth.f90

sbpress.o: sbpress.f90

	$(F90) -c sbpress.f90

eleleng.o: eleleng.f90

	$(F90) -c eleleng.f90

timefun.o: timefun.f90

	$(F90) -c timefun.f90

sur1d.o: sur1d.f90

	$(F90) -c sur1d.f90

press.o: press.f90

	$(F90) -c press.f90

ghost.o: ghost.f90

	$(F90) -c ghost.f90

jacob.o: jacob.f90

	$(F90) -c jacob.f90

bdpd.o: bdpd.f90

	$(F90) -c bdpd.f90

sstrain.o: sstrain.f90

	$(F90) -c sstrain.f90

scalkuu.o: scalkuu.f90

	$(F90) -c scalkuu.f90

force.o: force.f90

	$(F90) -c force.f90

smaterj.o: smaterj.f90

	$(F90) -c smaterj.f90

sboc.o: sboc.f90

	$(F90) -c sboc.f90

eleinit.o: eleinit.f90

	$(F90) -c eleinit.f90

clean:
	rm  *.o run   
	
