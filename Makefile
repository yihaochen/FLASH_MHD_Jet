# additional files required by the jet problem

Simulation += Simulation_data.o\
			  Simulation_jiggle.o\
			  Simulation_jetNozzleUpdate.o

Hydro += hy_uhd_electricNozzle.o\
		 hy_uhd_jetNozzleGeometry.o\
		 hy_uhd_jetNozzleGeometryOld.o\
		 hy_uhd_getA.o\
		 hy_uhd_getAOld.o

Heat += Heat_data.o\
		Heat_fillnozzle.o\
		ht_getValueAtPoint.o

Grid += gr_markJet.o

Particles += pt_getRandomPos.o
