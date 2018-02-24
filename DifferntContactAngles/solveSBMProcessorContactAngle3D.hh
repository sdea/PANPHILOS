#ifndef SOLVESBMPROCESSORCONTACTANGLE3D_HH_INCLUDED
#define SOLVESBMPROCESSORCONTACTANGLE3D_HH_INCLUDED

#include "solveSBMProcessorContactAngle3D.h"

template<typename T>
solveSBMContactAngle3D<T>::solveSBMContactAngle3D(T eps2_, T Q_, T invh2_, T dt_, T M_, T h_, T cosTheta_, bool SM_, bool useG_):
                    eps2(eps2_), Q(Q_), invh2(invh2_), dt(dt_), M(M_), h(h_), cosTheta(cosTheta_), SM(SM_),
		    useG(useG_)	{ }

 
// Override generic process method
template<typename T>
void solveSBMContactAngle3D<T>::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D* > atomicBlocks) {

     // Casting back the blocks to their specific type
     ScalarField3D<T>& mu = *dynamic_cast<ScalarField3D<T>* >(atomicBlocks[0]);
     ScalarField3D<T>& C = *dynamic_cast<ScalarField3D<T>* >(atomicBlocks[1]);
     ScalarField3D<T>& Psi = *dynamic_cast<ScalarField3D<T>* >(atomicBlocks[2]);
     ScalarField3D<T>& CosTheta = *dynamic_cast<ScalarField3D<T>* >(atomicBlocks[3]);

    // Compute displacments
     Dot3D offset_C = computeRelativeDisplacement(mu, C);
     Dot3D offset_Psi = computeRelativeDisplacement(mu, Psi);
     Dot3D offset_CosTheta = computeRelativeDisplacement(mu, CosTheta);

     // Loop over the indexes of the blocks
     for (plint iX=domain.x0; iX<=domain.x1; ++iX) {

         // X Displacements for block C and dPsiOvPsi
         plint iX_C = iX + offset_C.x;
         plint iX_Psi = iX + offset_Psi.x;
	 plint iX_CosTheta = iX + offset_CosTheta.x;

         for (plint iY=domain.y0; iY<=domain.y1; ++iY) {

            // Y Displacements for block C and dPsiOvPsi
            plint iY_C = iY + offset_C.y;
            plint iY_Psi = iY + offset_Psi.y;
	    plint iY_CosTheta = iY + offset_CosTheta.y;


            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {

                // Z Displacements for block C and dPsiOvPsi
                plint iZ_C = iZ + offset_C.z;
                plint iZ_Psi = iZ + offset_Psi.z;
		plint iZ_CosTheta = iZ + offset_CosTheta.z;

                // Compute dPsiOvPsi
                T dxPsi = ((Psi.get(iX_Psi +1, iY_Psi, iZ_Psi) - Psi.get(iX_Psi -1, iY_Psi, iZ_Psi))/(2.*h));
                T dyPsi = ((Psi.get(iX_Psi, iY_Psi +1, iZ_Psi) - Psi.get(iX_Psi, iY_Psi -1, iZ_Psi))/(2.*h));
                T dzPsi = ((Psi.get(iX_Psi, iY_Psi, iZ_Psi +1) - Psi.get(iX_Psi, iY_Psi, iZ_Psi -1))/(2.*h));

		T dPx = dxPsi/Psi.get(iX_Psi,iY_Psi,iZ_Psi);
		T dPy = dyPsi/Psi.get(iX_Psi,iY_Psi,iZ_Psi);
		T dPz = dzPsi/Psi.get(iX_Psi,iY_Psi,iZ_Psi);
                // Compute modGPsiOvPsi
                T  modGPsiOvPsi = std::sqrt(dPx*dPx + dPy*dPy + dPz*dPz);

                // Update Mu ******************************************************************
                T laplacianC = invh2*(C.get(iX_C +1, iY_C, iZ_C) + C.get(iX_C -1,iY_C,iZ_C)   +
                                      C.get(iX_C ,iY_C +1, iZ_C) + C.get(iX_C, iY_C -1, iZ_C) +
                                      C.get(iX_C, iY_C, iZ_C +1) + C.get(iX_C, iY_C, iZ_C -1) -
                                      6.*C.get(iX_C, iY_C, iZ_C));


                // Get derivative of C
                T dxC = (C.get(iX_C +1, iY_C, iZ_C) - C.get(iX_C -1,iY_C,iZ_C))/(2.*h);
                T dyC = (C.get(iX_C ,iY_C +1, iZ_C) - C.get(iX_C, iY_C -1, iZ_C))/(2.*h);
                T dzC = (C.get(iX_C, iY_C, iZ_C +1) - C.get(iX_C, iY_C, iZ_C -1))/(2.*h);


                T dotdC_dPsiOvPsi = (dPx*dxC + dPy*dyC + dPz*dzC);

                // Derivative of energy funcitonal
                T F = (Q/4.)*(C.get(iX_C,iY_C,iZ_C))*(C.get(iX_C,iY_C,iZ_C))*
                             (1 - C.get(iX_C,iY_C,iZ_C))*(1 - C.get(iX_C,iY_C,iZ_C));
                T dF = (Q/2.)*C.get(iX_C,iY_C,iZ_C)*(1. - C.get(iX_C,iY_C,iZ_C))*
                                                                (1. - 2.*C.get(iX_C,iY_C,iZ_C));

                mu.get(iX,iY,iZ) = dF- eps2*(dotdC_dPsiOvPsi + laplacianC) -
                                        std::sqrt(eps2)*
                                        (modGPsiOvPsi*
                                        std::sqrt(2*F)*CosTheta.get(iX_CosTheta,iY_CosTheta,iZ_CosTheta));


                // *****************************************************************************

                // Get derivatives of mu
                T dxMu = (mu.get(iX +1, iY, iZ) - mu.get(iX -1,iY,iZ))/(2.*h);
                T dyMu = (mu.get(iX ,iY +1, iZ) - mu.get(iX, iY -1, iZ))/(2.*h);
                T dzMu = (mu.get(iX, iY, iZ +1) - mu.get(iX, iY, iZ -1))/(2.*h);


                T laplacianMu = invh2*(mu.get(iX +1, iY, iZ) + mu.get(iX -1,iY,iZ)   +
                                      mu.get(iX ,iY +1, iZ) + mu.get(iX, iY -1, iZ) +
                                      mu.get(iX, iY, iZ +1) + mu.get(iX, iY, iZ -1) -
                                      6.*mu.get(iX, iY, iZ));

                T dotdMu_dPsiOvPsi = (dPx*dxMu + dPy*dyMu + dPz*dzMu);

                // Update C
                if(SM == false){
                	
			C.get(iX_C,iY_C,iZ_C) += M*dt*(dotdMu_dPsiOvPsi + laplacianMu);
		}
		else {

			// Mobility and its derivative
			T Mc = M*(C.get(iX_C,iY_C,iZ_C))*(C.get(iX_C,iY_C,iZ_C))*
				(1 - C.get(iX_C,iY_C,iZ_C))*(1 - C.get(iX_C,iY_C,iZ_C));
			T dMc = 2*M*C.get(iX_C,iY_C,iZ_C)*(1. - C.get(iX_C,iY_C,iZ_C))*
				(1. - 2.*C.get(iX_C,iY_C,iZ_C));

			// Default (useG == false)
			T G = 1.;
			T dG = 0.;

			if (useG == true) {
				
				// Compute the polinomial function G and its deivatives
				T psi = Psi.get(iX_Psi,iY_Psi,iZ_Psi);
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
			T dzM = G*dMc*dzC + Mc*dG*dzPsi;

			// Dot produt dM * dMu
			T dotdMdMu = dxM*dxMu + dyM*dyMu + dzM*dzMu;

			// Update C
			C.get(iX_C,iY_C,iZ_C) += dt*(Mc*G*dotdMu_dPsiOvPsi + dotdMdMu +
								Mc*G*laplacianMu);
                }
	    }
         }
     } // end for
}

template<typename T>
solveSBMContactAngle3D<T>*
        solveSBMContactAngle3D<T>::clone() const {

    return new solveSBMContactAngle3D<T>(*this);
}

 
template<typename T>
void solveSBMContactAngle3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {

    modified[0] = modif::staticVariables;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
    modified[3] = modif::nothing;

}


#endif // SOLVESBMPROCESSORCONTACTANGLE3D_HH_INCLUDED
