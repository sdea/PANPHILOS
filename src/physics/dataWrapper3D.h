#include "./ComputeSegmentation.h"
#include "./ComputeSegmentation.h"
#include "./solveSBMProcessor3D.h"
#include "./solveSBMProcessor3D.hh"
#include "./ApplyNoFlux3DFunctional.h"
#include "./ApplyNoFlux3DFunctional.hh"
#include "../io/InputOutput.h"

namespace plb {
// === List of wrapper functions around data processors 
template<typename T>
void compute_SBM_CH_step(std::vector<MultiBlock3D* >& blocksForSim, Box3D dmn, SimulationParameter3D<T>& param) {
	
     // Apply processing functional on the blocks
     T eps2 = param.eps*param.eps;
     T invH2 = 1./(param.h*param.h);
     applyProcessingFunctional(new solveSBM3D<T>(eps2, param.Q, invH2, param.dt, param.M , param.h, param.cosTheta, param.SMobility), dmn , blocksForSim);

}



}
