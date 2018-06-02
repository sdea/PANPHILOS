#include "./ComputeSegmentation.h"
#include "./ComputeSegmentation.h"
#include "./solveSBMProcessor3D_Davis.h"
#include "./solveSBMProcessor3D_Davis.hh"
#include "./ApplyNoFlux3DFunctional.h"
#include "./ApplyNoFlux3DFunctional.hh"
#include "../io/InputOutputDavis.h"

namespace plb {
// === List of wrapper functions around data processors 
template<typename T>
void compute_SBM_CH_step_Davis_dM(std::vector<MultiBlock3D* >& blocksForSim, Box3D dmn, SimulationParameter3D<T>& param) {
	
     // Apply processing functional on the blocks
     T eps2 = param.eps*param.eps;
     T invH2 = 1./(param.h*param.h);
     applyProcessingFunctional(new solveSBM3D_Davis_dM<T>(eps2, param.Q, invH2, param.dt, param.M , param.h, param.cosTheta, param.a1, param.a2, param.a3, param.SMobility, param.useG), dmn , blocksForSim);

}



}
