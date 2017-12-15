#ifndef INPUTOUTPUT_H_INCLUDED
#define INPUTOUTPUT_H_INCLUDED

#include "../physics/ComputeSegmentation.h"
#include "../physics/ComputeSegmentation.hh"

/* Header file with classes and functions for input-output */ 

namespace plb {


// =================================================================================================
// =================== I N P U T ===================================================================

// Class of simulation parameters
template <typename T>
class SimulationParameter3D {
public:

      SimulationParameter3D(std::string xmlFileName) {

            /// Reading xml file
            XMLreader xmlFile(xmlFileName.c_str());
            xmlFile["Geometry"]["inputCFile"].read(fileCName);
            xmlFile["Geometry"]["inputPsiFile"].read(filePsiName);
            xmlFile["Geometry"]["size"]["nx"].read(nx);
            xmlFile["Geometry"]["size"]["ny"].read(ny);
            xmlFile["Geometry"]["size"]["nz"].read(nz);
            xmlFile["Simulation"]["Q"].read(Q);
            xmlFile["Simulation"]["Epsilon"].read(eps);
            xmlFile["Simulation"]["Mobility"].read(M);
            xmlFile["Simulation"]["Theta"].read(theta);
            xmlFile["Numerics"]["h"].read(h);
            xmlFile["Numerics"]["RedFactor"].read(red_factor);
            xmlFile["Numerics"]["SurfMobility"].read(SMobility);
            xmlFile["Output"]["OutputDatFile"].read(outDatName);
            xmlFile["Output"]["OutputVTKFile"].read(outVTKName);
            xmlFile["Output"]["MaxIter"].read(maxIter);
            xmlFile["Output"]["FreqOutput"].read(freqOut);




            /// Assign remaining parameters
            const T dim_domain = 3;
	    const T MaxM = 1.;
            dt = red_factor*std::pow(h,4)/(MaxM*std::pow(2, 2*dim_domain +1));

            const T pi = 3.14159265;
            cosTheta = std::cos((theta*pi)/180.);

            }

        void print() {

            /// Print simulation parameters to screen
            pcout << "\n===========================================================" << std::endl;
            pcout << "=====  S I M U L A T I O N    P A R A M E T E R S  ========" << std::endl;
            pcout << "" << std::endl;
            pcout << "*      Input C File:     " << this->fileCName << std::endl;
            pcout << "*      Input Psi File:   " << this->filePsiName << std::endl;
            pcout << "*      Nx:               " << this->nx << std::endl;
            pcout << "*      Ny:               " << this->ny << std::endl;
            pcout << "*      Nz:               " << this->nz << std::endl;
            pcout << "*      Q:                " << this->Q << std::endl;
            pcout << "*      Eps:              " << this->eps << std::endl;
            pcout << "*      M:                " << this->M << std::endl;
            pcout << "*      Theta [deg]:      " << this->theta << std::endl;
            pcout << "*      h:                " << this->h << std::endl;
            pcout << "*      dt:               " << this->dt << std::endl;
            pcout << "*      VTK output name:  " << this->outVTKName << std::endl;
            pcout << "*      DAT output name:  " << this->outDatName << std::endl;
            pcout << "*      MaxIter:          " << this->maxIter << std::endl;
            pcout << "*      FreqOut:          " << this->freqOut << std::endl;
            pcout << "*      SurfMobility:     " << this->SMobility << std::endl;
            pcout << "" << std::endl;
            pcout << "===========================================================" << std::endl;
            pcout << "===========================================================" << std::endl;
        }

public:

	// XML input file name
	std::string fileCName, filePsiName;
    	// Dimension of geometry 2D
    	plint nx, ny, nz;
    	// Integration and derivation constants^M
    	T h, red_factor, dt;
    	// Parameters related to phase field simulations
    	T Q, eps, M, theta, cosTheta;
	bool SMobility;
    	// Output files names
    	std::string outDatName, outVTKName;
    	// For variables
    	plint maxIter, freqOut;
};

// Read geometries from file
template <typename T>
void readGeom3D(std::string filename, MultiScalarField3D<T> &scalarField) {

    plb_ifstream geometryFile(filename.c_str());
    if (!geometryFile.is_open()) {
            pcout << "Error: could not open geometry file " << filename << std::endl;
            exit(EXIT_FAILURE); }

    geometryFile >> scalarField;

}

// ===============================================================================================
// ========================== O U T P U T ========================================================


/// Functions which writes a VTI file ***
// Single parameter, assign name to parameter via
template <typename T>
void writeSingleVTK(MultiScalarField3D<T>& field, std::string filename, plint iter, plint num_zeros)

{
        VtkImageOutput3D<T> vtkOut(createFileName(filename, iter, num_zeros), 1.);
        vtkOut.template writeData<float>(field, "orderParameter", 1.);
}


// Write a full VTK with segmented structure in it
template <typename T>
void writeFullVTK(std::vector<MultiBlock3D* > &fields, Box3D domain, std::string fileName, plint iter, plint num_zeros) {

    pcout << "Writing VTK file..." << std::endl;

    // Populate the geom scalarfield
    applyProcessingFunctional(new ComputeSegmentation<T>(), domain, fields);

    // Cast back types
    MultiScalarField3D<T>& geom = *dynamic_cast<MultiScalarField3D<T>* >(fields[0]);
    MultiScalarField3D<T>& C = *dynamic_cast<MultiScalarField3D<T>* >(fields[1]);

    // Write VTK
    VtkImageOutput3D<T> vtkOut(createFileName(fileName, iter, num_zeros), 1.);
    vtkOut.template writeData<float>(C, "orderParameter", 1.);
    vtkOut.template writeData<float>(geom, "segmentedGeom", 1.);

}

/// Functions which write a DAT file ***
template <typename T>
void writeDat3D(MultiScalarField3D<T> &field, std::string filename, plint iter, plint num_images) {


    pcout << "Writing DAT file..." << std::endl;
    std::string outputName = createFileName(filename, iter, num_images);

    // Stream outfile in a text file
    plb_ofstream outFile(outputName.c_str());
    outFile  << field;

};



} // end namspace plb



#endif // INPUTOUTPUT_H_INCLUDED
