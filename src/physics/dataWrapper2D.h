#include "./ComputeSegmentation2D.h"
#include "./ComputeSegmentation2D.h"
#include "./solveSBMProcessor2D.h"
#include "./solveSBMProcessor2D.hh"
#include "./ApplyNoFlux2DFunctional.h"
#include "./ApplyNoFlux2DFunctional.hh"
#include "../io/InputOutput.h"

namespace plb {
// === List of wrapper functions around data processors 
template<typename T>
void compute_SBM_CH_step2D(std::vector<MultiBlock2D* >& blocksForSim, Box2D dmn, SimulationParameter2D<T>& param) {
	
     // Apply processing functional on the blocks
     T eps2 = param.eps*param.eps;
     T invH2 = 1./(param.h*param.h);
     applyProcessingFunctional(new solveSBM2D<T>(eps2, param.Q, invH2, param.dt, param.M , param.h, param.cosTheta, param.SMobility, param.useG), dmn , blocksForSim);

}



}
