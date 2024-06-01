CC = g++
CFLAGS = -Wall -O3 --std=c++11

NVE_NVT.exe : MD_MC.o random.o
	$(CC) random.o MD_MC.o -o NVE_NVT.exe
MD_MC.o : MD_MC.cpp
	$(CC) -c MD_MC.cpp -o MD_MC.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o NVE_NVT.exe
cleansolid :
	rm output_ekin_solid.dat output_epot_solid.dat output_etot_solid.dat output_pres_solid.dat output_temp_solid.dat 
cleanliquid :
	rm output_ekin_liquid.dat output_epot_liquid.dat output_etot_liquid.dat output_pres_liquid.dat output_temp_liquid.dat 
cleangas :
	rm output_ekin_gas.dat output_epot_gas.dat output_etot_gas.dat output_pres_gas.dat output_temp_gas.dat 
