// Includes for 3D version of the code
#include "palabos3D.h"
#include "palabos3D.hh"
#include "math.h"

typedef double T;
using namespace plb;

#include <ctime>

// Data processors
#include "../src/io/InputOutputDavis.h"
#include "../src/physics/ApplyNoFlux3DFunctional.h"
#include "../src/physics/ApplyNoFlux3DFunctional.hh"
#include "../src/physics/solveSBMProcessor3D_Davis.h"
#include "../src/physics/solveSBMProcessor3D_Davis.hh"
#include "../src/physics/dataWrapper3D_Davis.h"

void applyNoFlux3Dfunction(MultiScalarField3D<T>& mu, MultiScalarField3D<T>& C) {

     const plint nx = C.getNx();
     const plint ny = C.getNy();
     const plint nz = C.getNz();

     Box3D left(0, 1, 0, ny-1, 0, nz-1);
     Box3D right(nx-2, nx-1, 0, ny-1, 0, nz-1);
     Box3D bottom(0, nx-1, 0, 1, 0, nz-1);
     Box3D top(0, nx-1, ny-2, ny-1, 0, nz-1);
     Box3D frontP(0, nx-1, 0, ny-1, 0, 1);
     Box3D backP(0, nx-1, 0, ny-1, nz-2, nz-1);

     // No flux for mu
     applyProcessingFunctional(new applyNoFlux3D<T>(0,-1), left, mu);
     applyProcessingFunctional(new applyNoFlux3D<T>(0, 1), right, mu);
     applyProcessingFunctional(new applyNoFlux3D<T>(1,-1), bottom, mu);
     applyProcessingFunctional(new applyNoFlux3D<T>(1, 1), top, mu);
     applyProcessingFunctional(new applyNoFlux3D<T>(2,-1), frontP, mu);
     applyProcessingFunctional(new applyNoFlux3D<T>(2, 1), backP, mu);

     // No-flux for C
     applyProcessingFunctional(new applyNoFlux3D<T>(0,-1), left, C);
     applyProcessingFunctional(new applyNoFlux3D<T>(0, 1), right, C);
     applyProcessingFunctional(new applyNoFlux3D<T>(1,-1), bottom, C);
     applyProcessingFunctional(new applyNoFlux3D<T>(1, 1), top, C);
     applyProcessingFunctional(new applyNoFlux3D<T>(2,-1), frontP, C);
     applyProcessingFunctional(new applyNoFlux3D<T>(2, 1), backP, C);
}

int main(int argc, char *argv[]) {

    // Init simulation here
    plbInit(&argc, &argv);

    if (argc != 2) {

	pcout << "Wrong number of input parameters" << std::endl;
	return -1;

    }

    SimulationParameter3D<T> param(argv[1]);
    param.print();


    // Allocate needed scalar field
    MultiScalarField3D<T> mu(param.nx, param.ny, param.nz);
    MultiScalarField3D<T> C(param.nx,  param.ny , param.nz);
    MultiScalarField3D<T> Psi(param.nx, param.ny, param.nz);
    MultiScalarField3D<T> geom(param.nx, param.ny, param.nz);
    MultiScalarField3D<T> Mc(param.nx, param.ny, param.nz);

    // Read geometries
    readGeom3D(param.fileCName, C);
    readGeom3D(param.filePsiName, Psi);

    // Test
    /*pcout << "Warning: Test for Mass" << std::endl;
    pcout << "Psi == 1 everywhere"    << std::endl;
    pcout << "Remove the section for real studies" << std::endl;
    setToConstant(Psi, Psi.getBoundingBox(), 1.);*/	

    // Write Inital vtks
    std::vector<MultiBlock3D*> blocksForOutput;
    blocksForOutput.push_back(&geom);
    blocksForOutput.push_back(&C);
    blocksForOutput.push_back(&Psi);
    
    // Blocks for simulations
    std::vector<MultiBlock3D*> blocksForSim;
    blocksForSim.push_back(&mu);
    blocksForSim.push_back(&C);
    blocksForSim.push_back(&Psi);
    blocksForSim.push_back(&Mc);

    // Only for test
    // writeSingleVTK<T>(Psi, "PsiIni", 0, 0);


    // Domain of the simulation
    Box3D dmn(1, param.nx -2, 1, param.ny -2, 1, param.nz -2);
    
    clock_t begin = clock();
    for (plint iTT=param.initIter; iTT <= param.maxIter; ++iTT) {
    
      if (iTT%1000 == 0) {
	
	    T mass = computeSum(C,C.getBoundingBox());	
            pcout << "Iteration:  " << iTT << "Mass: " << mass << std::endl;
	
       }

	if (iTT%param.freqOut == 0) {
            	
		// Output
            	pcout << "Printing output:  " << iTT << std::endl;
        	writeDat3D<T>(C, param.outDatName, iTT, 15);
        	writeFullVTK<T>(blocksForOutput, C.getBoundingBox() , param.outVTKName, iTT, 15);
	   }

	// Cahn-Hillard step
	compute_SBM_CH_step_Davis_dM<T>(blocksForSim, dmn, param);
        
	// Boundary conditions
        applyNoFlux3Dfunction(mu, C);
    }
    
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    pcout << "Elapsed time:   " << elapsed_secs << std::endl;

    pcout << "Hello palabos!!!" << std::endl;
    return 0;
}
