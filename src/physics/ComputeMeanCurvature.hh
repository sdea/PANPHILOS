#ifndef COMPUTEMEANCURVATURE_HH_INCLUDED
#define COMPUTEMEANCURVATURE_HH_INCLUDED

#include "ComputeMeanCurvature.h"

namespace plb {

template<typename T>
ComputeMeanCurvature<T>::ComputeMeanCurvature(T h_): h(h_)
                    { }

// Override generic process method
template<typename T>
void ComputeMeanCurvature<T>::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D* > atomicBlocks) {

     // Casting back the blocks to their specific type
     ScalarField3D<T>& MeanCurv = *dynamic_cast<ScalarField3D<T>* >(atomicBlocks[0]);
     ScalarField3D<T>& C = *dynamic_cast<ScalarField3D<T>* >(atomicBlocks[1]);

    // Compute displacments
     Dot3D offset_C = computeRelativeDisplacement(MeanCurv, C);

     // Loop over the indexes of the blocks
     for (plint iX=domain.x0; iX<=domain.x1; ++iX) {

         // X Displacements for block
         plint iX_C = iX + offset_C.x;

         for (plint iY=domain.y0; iY<=domain.y1; ++iY) {

            // Y Displacements for block C
            plint iY_C = iY + offset_C.y;

            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {

                // Z Displacements for block C
                plint iZ_C = iZ + offset_C.z;


		if (C.get(iX_C, iY_C, iZ_C) > 0.1 && C.get(iX_C, iY_C, iZ_C) < 0.9){
		// Compute C gradient
	        T dxC = (C.get(iX_C + 1, iY_C, iZ_C) - C.get(iX_C - 1, iY_C, iZ_C))/(2.*h);
		T dyC = (C.get(iX_C, iY_C + 1, iZ_C) - C.get(iX_C, iY_C -1, iZ_C))/(2.*h);
		T dzC = (C.get(iX_C, iY_C, iZ_C + 1) - C.get(iX_C, iY_C, iZ_C - 1))/(2.*h);

		// Compute ModGradC
		T ModGradC = std::sqrt(dxC*dxC + dyC*dyC + dzC*dzC);
		
		// Compute gradC/ModGradC
		T dCx = dxC/ModGradC;
		T dCy = dyC/ModGradC;
		T dCz = dzC/ModGradC;

		// Compute Mean Curv
		MeanCurv.get(iX, iY, iZ) = dCx; // (dCx*dCx + dCy*dCy + dCz*dCz);		
		
		pcout << "Curv:  " << dCx*dCx + dCy*dCy + dCz*dCz << std::endl;

		}
	   }
	}
     } // End of for

} // End of process function

template<typename T>
ComputeMeanCurvature<T>*
        ComputeMeanCurvature<T>::clone() const {

    return new ComputeMeanCurvature<T>(*this);
}

template<typename T>
void ComputeMeanCurvature<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {

    modified[0] = modif::allVariables;
    modified[1] = modif::nothing;

}

} // End namespace

#endif // COMPUTEMEANCURVATURE_HH_INCLUDED
