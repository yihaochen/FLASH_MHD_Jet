# additional files required by the jet problem

Simulation += Simulation_data.o\
			  Simulation_jetNozzleUpdate.o

Hydro += hy_uhd_electricNozzle.o\
		 hy_uhd_jetNozzleGeometry.o\
		 hy_uhd_jetNozzleGeometryOld.o\
		 hy_uhd_getA.o\
		 hy_uhd_getAOld.o

Heat += Heat_fillnozzle.o

Grid += gr_markJet.o
