// Includes for 2D version of the code
#include "palabos2D.h"
#include "palabos2D.hh"

typedef double T;
using namespace plb;

#include <ctime>

// Data processors
#include "../src/io/InputOutput.h"
#include "../src/physics/ApplyNoFlux2DFunctional.h"
#include "../src/physics/ApplyNoFlux2DFunctional.hh"
#include "./ApplyPeriodicBC2DFunctional.h"
#include "./ApplyPeriodicBC2DFunctional.hh"
#include "../src/physics/solveSBMProcessor2D.h"
#include "../src/physics/solveSBMProcessor2D.hh"
#include "../src/physics/dataWrapper2D.h"


void applyNoFlux2Dfunction(MultiScalarField2D<T>& mu, MultiScalarField2D<T>& C) {

     const plint nx = C.getNx();
     const plint ny = C.getNy();

     Box2D left(0, 1, 0, ny-1);
     Box2D right(nx-2, nx-1, 0, ny-1);
     Box2D bottom(0, nx-1, 0, 1);
     Box2D top(0, nx-1, ny-2, ny-1);

     // Domain for periodic BC
     // The idea here is to create two shifted boxes
     // The first (leftPeriodic) will serve to assign the left ghost layer
     // to the last known layer
     // The second (rightPeriodic) does the opposite and serve to 
     // assigning the right ghost layer 
     
     // NOTE!!! *******************************************************************
     // With the trick of creating these domains we avoid the sitations
     // of performing non local operations (e.g domain.x0 = domain.x1 -1)
     // The -1 operation is non local, because we attempt to access a cell 
     // which can be placed in an another block (when the program runs in parallel)
     // This is dangerous in parallel, mainly when we perform operations on the 
     // envelope of a block (always the case when imposing BCs)
     // ***************************************************************************
     Box2D leftPeriodic(0, nx-2, 0, ny -1);
     Box2D rightPeriodic(1, nx-1, 0, ny -1);
     
     // No flux for mu
     // Preserve conservation of mass
     applyProcessingFunctional(new applyNoFlux2D<T>(0,-1), left, mu);
     applyProcessingFunctional(new applyNoFlux2D<T>(0, 1), right, mu);
     applyProcessingFunctional(new applyNoFlux2D<T>(1,-1), bottom, mu);
     applyProcessingFunctional(new applyNoFlux2D<T>(1, 1), top, mu);

     // No-flux for C on y direction (dir = 1)
     applyProcessingFunctional(new applyNoFlux2D<T>(1,-1), bottom, C);
     applyProcessingFunctional(new applyNoFlux2D<T>(1, 1), top, C);

     // Here we apply the periodic bounary conditions on x (dir = 0)
     // The trick is to pass the right domains to the dataprocessor
     applyProcessingFunctional(new applyPeriodicBC2D<T>(0,-1), leftPeriodic, C);
     applyProcessingFunctional(new applyPeriodicBC2D<T>(0, 1), rightPeriodic, C);

    
}

int main(int argc, char *argv[]) {

    // Init simulation here
    plbInit(&argc, &argv);

    if (argc != 2) {

	pcout << "Wrong number of input parameters" << std::endl;
	return -1;

    }

    SimulationParameter2D<T> param(argv[1]);
    param.print();


    // Allocate needed scalar field
    MultiScalarField2D<T> mu(param.nx, param.ny);
    MultiScalarField2D<T> C(param.nx,  param.ny);
    MultiScalarField2D<T> Psi(param.nx, param.ny);
    MultiScalarField2D<T> geom(param.nx, param.ny);

    // Read geometries
    readGeom2D(param.fileCName, C);
    readGeom2D(param.filePsiName, Psi);

    // Test
    /*pcout << "Warning: Test for Mass" << std::endl;
    pcout << "Psi == 1 everywhere"    << std::endl;
    pcout << "Remove the section for real studies" << std::endl;
    setToConstant(Psi, Psi.getBoundingBox(), 1.);*/	

    // Write Inital vtks
    std::vector<MultiBlock2D*> blocksForOutput;
    blocksForOutput.push_back(&geom);
    blocksForOutput.push_back(&C);
    blocksForOutput.push_back(&Psi);
    
    // Blocks for simulations
    std::vector<MultiBlock2D*> blocksForSim;
    blocksForSim.push_back(&mu);
    blocksForSim.push_back(&C);
    blocksForSim.push_back(&Psi);

    // Only for test
    writeSingleVTK2D<T>(Psi, "PsiIni", 0, 0);


    // Domain of the simulation
    Box2D dmn(1, param.nx -2, 1, param.ny -2);
    
    clock_t begin = clock();
    for (plint iTT=param.initIter; iTT < param.maxIter; ++iTT) {

        pcout << "Iteration:  " << iTT << std::endl;
	    
	if (iTT%param.freqOut == 0) {
            	
		// Output
            	pcout << "Printing output:  " << iTT << std::endl;
        	writeDat2D<T>(C, param.outDatName, iTT, 15);
        	writeFullVTK2D<T>(blocksForOutput, C.getBoundingBox() , param.outVTKName, iTT, 15);
	   }

	// Cahn-Hillard step
	compute_SBM_CH_step2D<T>(blocksForSim, dmn, param);
        
	// Boundary conditions
        applyNoFlux2Dfunction(mu, C);
    }
    
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    pcout << "Elapsed time:   " << elapsed_secs << std::endl;

    pcout << "Hello palabos!!!" << std::endl;
    return 0;
}
