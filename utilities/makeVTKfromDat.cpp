#include "palabos3D.h"
#include "palabos3D.hh"

#include "../src/io/InputOutput.h"

using namespace plb;
typedef double T;

int main(int argc, char* argv[]) {

    // Init palabos simulatiion
    plbInit(&argc, &argv);

    // Read from std input
    if (argc != 5) {
       
       std::cout << "Wrong number of input paramters" << std::endl;
       std::cout << "File name and dimensions  required!!!" << std::endl;
       return -1;

    }
    
    // File name
    std::string fileName = argv[1];
    const plint nx = atoi(argv[2]);
    const plint ny = atoi(argv[3]);
    const plint nz = atoi(argv[4]);


    // Read the .dat file
    MultiScalarField3D<T> field(nx, ny, nz);
    readGeom3D(fileName, field);

    // Write VTK
    fileName.erase(fileName.find_last_of("."), std::string::npos);
    writeGeneralVTK3D<T>(field, fileName, "content"); 



    return 0;

}
