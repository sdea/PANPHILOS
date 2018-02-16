#ifndef COMPUTESEGMENTATION2D_HH_INCLUDED
#define COMPUTESEGMENTATION2D_HH_INCLUDED

#include "ComputeSegmentation2D.h"

namespace plb {

template<typename T>
ComputeSegmentation2D<T>::ComputeSegmentation2D()
                    { }

// Override generic process method
template<typename T>
void ComputeSegmentation2D<T>::processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D* > atomicBlocks) {

     // Casting back the blocks to their specific type
     ScalarField2D<T>& geom = *dynamic_cast<ScalarField2D<T>* >(atomicBlocks[0]);
     ScalarField2D<T>& C = *dynamic_cast<ScalarField2D<T>* >(atomicBlocks[1]);
     ScalarField2D<T>& Psi = *dynamic_cast<ScalarField2D<T>* >(atomicBlocks[2]);

    // Compute displacments
     Dot2D offset_C = computeRelativeDisplacement(geom, C);
     Dot2D offset_Psi = computeRelativeDisplacement(geom, Psi);

     // Loop over the indexes of the blocks
     for (plint iX=domain.x0; iX<=domain.x1; ++iX) {

         // X Displacements for block C and dPsiOvPsi
         plint iX_C = iX + offset_C.x;
         plint iX_Psi = iX + offset_Psi.x;

         for (plint iY=domain.y0; iY<=domain.y1; ++iY) {

            // Y Displacements for block C and dPsiOvPsi
            plint iY_C = iY + offset_C.y;
            plint iY_Psi = iY + offset_Psi.y;

		// Perform segmetation
		if(Psi.get(iX_Psi, iY_Psi) <= 0.5) {

			geom.get(iX, iY) = 1;
		}
		else {
			if (C.get(iX_C, iY_C) >= 0.5) {

				geom.get(iX, iY) = 2;

			}
			else {

				geom.get(iX, iY) = 0;
			}

		}
	}
     } // End of for

} // End of process function

template<typename T>
ComputeSegmentation2D<T>*
        ComputeSegmentation2D<T>::clone() const {

    return new ComputeSegmentation2D<T>(*this);
}

template<typename T>
void ComputeSegmentation2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {

    modified[0] = modif::allVariables;
    modified[1] = modif::nothing;
    modified[2] = modif::nothing;

}

} // End namespace

#endif // COMPUTESEGMENTATION2D_HH_INCLUDED
