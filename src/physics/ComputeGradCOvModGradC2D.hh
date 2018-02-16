#ifndef COMPUTEGRADCOVMODGRADC2D_HH_INCLUDED
#define COMPUTEGRADCOVMODGRADC2D_HH_INCLUDED

#include "ComputeGradCOvModGradC2D.h"

namespace plb {

template<typename T1, typename T2, int nDim>
ComputeGradCOvModGradC2D<T1, T2, nDim>::ComputeGradCOvModGradC2D(T h_): h(h_)
                    { }

// Override generic process method
template<typename T1, typename T2, int nDim>
void ComputeGradCOvModGradC2D<T1, T2, nDim>::process(Box2D domain, ScalarField2D<T1>& C, TensorField2D<T2,nDim>& GradCOvModGradC) {

    // Compute displacments
     Dot2D offset_C = computeRelativeDisplacement(GradCOvModGradC, C);

     // Loop over the indexes of the blocks
     for (plint iX=domain.x0; iX<=domain.x1; ++iX) {

         // X Displacements for block
         plint iX_C = iX + offset_C.x;

         for (plint iY=domain.y0; iY<=domain.y1; ++iY) {

            // Y Displacements for block C
            plint iY_C = iY + offset_C.y;


		// Compute GradC/ModGradC only at the interface
		if (C.get(iX_C, iY_C) > 0.1 && C.get(iX_C, iY_C) < 0.9){


		// Compute C gradient
	        T dxC = (C.get(iX_C + 1, iY_C) - C.get(iX_C - 1, iY_C))/(2.*h);
		T dyC = (C.get(iX_C, iY_C + 1) - C.get(iX_C, iY_C - 1))/(2.*h);

		// Compute ModGradC
		T ModGradC = std::sqrt(dxC*dxC + dyC*dyC);
		
		// Compute gradC/ModGradC
		GradCOvModGradC.get(iX, iY)[0] = dxC/ModGradC;
		GradCOvModGradC.get(iX, iY)[1] = dyC/ModGradC;
				
     	   }
	}
     } // End of for

} // End of process function

template<typename T1, typename T2, int nDim>
ComputeGradCOvModGradC2D<T1, T2, nDim>*
        ComputeGradCOvModGradC2D<T1, T2, nDim>::clone() const {

    return new ComputeGradCOvModGradC2D<T1, T2, nDim>(*this);
}

template<typename T1, typename T2, int nDim>
void ComputeGradCOvModGradC2D<T1, T2, nDim>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {

    modified[0] = modif::allVariables;
    modified[1] = modif::nothing;

}

} // End namespace

#endif // COMPUTEGRADCOVMODGRADC2D_HH_INCLUDED
