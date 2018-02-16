#ifndef SOLVESBMPROCESSOR2D_HH_INCLUDED
#define SOLVESBMPROCESSOR2D_HH_INCLUDED

#include "solveSBMProcessor2D.h"

template<typename T>
solveSBM2D<T>::solveSBM2D(T eps2_, T Q_, T invh2_, T dt_, T M_, T h_, T cosTheta_, bool SM_, bool useG_):
                    eps2(eps2_), Q(Q_), invh2(invh2_), dt(dt_), M(M_), h(h_), cosTheta(cosTheta_), SM(SM_),
		    useG(useG_)	{ }

// Override generic process method
template<typename T>
void solveSBM2D<T>::processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D* > atomicBlocks) {

     // Casting back the blocks to their specific type
     ScalarField2D<T>& mu = *dynamic_cast<ScalarField2D<T>* >(atomicBlocks[0]);
     ScalarField2D<T>& C = *dynamic_cast<ScalarField2D<T>* >(atomicBlocks[1]);
     ScalarField2D<T>& Psi = *dynamic_cast<ScalarField2D<T>* >(atomicBlocks[2]);

    // Compute displacments
     Dot2D offset_C = computeRelativeDisplacement(mu, C);
     Dot2D offset_Psi = computeRelativeDisplacement(mu, Psi);

     // Loop over the indexes of the blocks
     for (plint iX=domain.x0; iX<=domain.x1; ++iX) {

         // X Displacements for block C and dPsiOvPsi
         plint iX_C = iX + offset_C.x;
         plint iX_Psi = iX + offset_Psi.x;


         for (plint iY=domain.y0; iY<=domain.y1; ++iY) {

            // Y Displacements for block C and dPsiOvPsi
            plint iY_C = iY + offset_C.y;
            plint iY_Psi = iY + offset_Psi.y;


                // Compute dPsiOvPsi
                T dxPsi = (Psi.get(iX_Psi +1, iY_Psi) - Psi.get(iX_Psi -1, iY_Psi))/(2.*h);
                T dyPsi = (Psi.get(iX_Psi, iY_Psi +1) - Psi.get(iX_Psi, iY_Psi -1))/(2.*h);

		T dPx = dxPsi/Psi.get(iX_Psi,iY_Psi);
		T dPy = dyPsi/Psi.get(iX_Psi,iY_Psi);

                // Compute modGPsiOvPsi
                T  modGPsiOvPsi = std::sqrt(dPx*dPx + dPy*dPy);

                // Update Mu ******************************************************************
                T laplacianC = invh2*(C.get(iX_C +1, iY_C) + C.get(iX_C -1,iY_C)   +
                                      C.get(iX_C ,iY_C +1) + C.get(iX_C, iY_C -1)  -
                                      4.*C.get(iX_C, iY_C));


                // Get derivative of C
                T dxC = (C.get(iX_C +1, iY_C) - C.get(iX_C -1,iY_C))/(2.*h);
                T dyC = (C.get(iX_C ,iY_C +1) - C.get(iX_C, iY_C -1))/(2.*h);


                T dotdC_dPsiOvPsi = (dPx*dxC + dPy*dyC);

                // Derivative of energy funcitonal
                T F = (Q/4.)*(C.get(iX_C,iY_C))*(C.get(iX_C,iY_C))*
                             (1 - C.get(iX_C,iY_C))*(1 - C.get(iX_C,iY_C));
                T dF = (Q/2.)*C.get(iX_C,iY_C)*(1. - C.get(iX_C,iY_C))*
                                                                (1. - 2.*C.get(iX_C,iY_C));

                mu.get(iX,iY) = dF- eps2*(dotdC_dPsiOvPsi + laplacianC) -
                                        std::sqrt(eps2)*
                                        modGPsiOvPsi*
                                        std::sqrt(2*F)*cosTheta;



                // *****************************************************************************

                // Get derivatives of mu
                T dxMu = (mu.get(iX +1, iY) - mu.get(iX -1,iY))/(2.*h);
                T dyMu = (mu.get(iX ,iY +1) - mu.get(iX, iY -1))/(2.*h);


                T laplacianMu = invh2*(mu.get(iX +1, iY) + mu.get(iX -1,iY)   +
                                      mu.get(iX ,iY +1) + mu.get(iX, iY -1) -
                                      4.*mu.get(iX, iY));

                T dotdMu_dPsiOvPsi = (dPx*dxMu + dPy*dyMu);

                // Update C
                if(SM == false){
                	
			C.get(iX_C,iY_C) += M*dt*(dotdMu_dPsiOvPsi + laplacianMu);
		}
		else {

			// Mobility and its derivative
			T Mc = M*(C.get(iX_C,iY_C))*(C.get(iX_C,iY_C))*
				(1 - C.get(iX_C,iY_C))*(1 - C.get(iX_C,iY_C));
			T dMc = 2*M*C.get(iX_C,iY_C)*(1. - C.get(iX_C,iY_C))*
				(1. - 2.*C.get(iX_C,iY_C));

			// Default (useG == false)
			T G = 1.;
			T dG = 0.;

			if (useG == true) {
				
				// Compute the polinomial function G and its deivatives
				T psi = Psi.get(iX_Psi,iY_Psi);
				T psi6 = std::pow(psi, 6);
				T psi5 = std::pow(psi, 5);
				T psi2 = psi*psi;

				// G function
				G = psi6*(10.*psi2 - 15.*psi + 6.);
				dG = psi5*(80.*psi2 - 105.*psi + 36.);

			} 

			// Derivatives of surface mobility
			T dxM = G*dMc*dxC + Mc*dG*dxPsi;
			T dyM = G*dMc*dyC + Mc*dG*dyPsi;

			// Dot produt dM * dMu
			T dotdMdMu = dxM*dxMu + dyM*dyMu;

			// Update C
			C.get(iX_C,iY_C) += dt*(Mc*G*dotdMu_dPsiOvPsi + dotdMdMu +
								Mc*G*laplacianMu);
                }
         }
     } // end for
}

template<typename T>
solveSBM2D<T>*
        solveSBM2D<T>::clone() const {

    return new solveSBM2D<T>(*this);
}

 
template<typename T>
void solveSBM2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {

    modified[0] = modif::staticVariables;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
    modified[3] = modif::nothing;

}


#endif // SOLVESBMPROCESSOR2D_HH_INCLUDED
