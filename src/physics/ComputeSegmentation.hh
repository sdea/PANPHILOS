#ifndef COMPUTESEGMENTATION_HH_INCLUDED
#define COMPUTESEGMENTATION_HH_INCLUDED

#include "ComputeSegmentation.h"

namespace plb {

template<typename T>
ComputeSegmentation<T>::ComputeSegmentation()
                    { }

// Override generic process method
template<typename T>
void ComputeSegmentation<T>::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D* > atomicBlocks) {

     // Casting back the blocks to their specific type
     ScalarField3D<T>& geom = *dynamic_cast<ScalarField3D<T>* >(atomicBlocks[0]);
     ScalarField3D<T>& C = *dynamic_cast<ScalarField3D<T>* >(atomicBlocks[1]);
     ScalarField3D<T>& Psi = *dynamic_cast<ScalarField3D<T>* >(atomicBlocks[2]);

    // Compute displacments
     Dot3D offset_C = computeRelativeDisplacement(geom, C);
     Dot3D offset_Psi = computeRelativeDisplacement(geom, Psi);

     // Loop over the indexes of the blocks
     for (plint iX=domain.x0; iX<=domain.x1; ++iX) {

         // X Displacements for block C and dPsiOvPsi
         plint iX_C = iX + offset_C.x;
         plint iX_Psi = iX + offset_Psi.x;

         for (plint iY=domain.y0; iY<=domain.y1; ++iY) {

            // Y Displacements for block C and dPsiOvPsi
            plint iY_C = iY + offset_C.y;
            plint iY_Psi = iY + offset_Psi.y;

            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {

                // Z Displacements for block C and dPsiOvPsi
                plint iZ_C = iZ + offset_C.z;
                plint iZ_Psi = iZ + offset_Psi.z;

		// Perform segmetation
		if(Psi.get(iX_Psi, iY_Psi, iZ_Psi) >= 0.5) {

			geom.get(iX, iY, iZ) = 1;
		}
		else {
			if (C.get(iX_C, iY_C, iZ_C) >= 0.5) {

				geom.get(iX, iY, iZ) = 2;

			}
			else {

				geom.get(iX, iY, iZ) = 0;
			}

		}


	   }
	}
     } // End of for

} // End of process function

template<typename T>
ComputeSegmentation<T>*
        ComputeSegmentation<T>::clone() const {

    return new ComputeSegmentation<T>(*this);
}

template<typename T>
void ComputeSegmentation<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {

    modified[0] = modif::allVariables;
    modified[1] = modif::nothing;
    modified[2] = modif::nothing;

}

} // End namespace

#endif // COMPUTESEGMENTATION_HH_INCLUDED
