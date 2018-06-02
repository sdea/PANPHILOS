#ifndef SOLVESBMPROCESSOR2D_IMPLICIT_HH_INCLUDED
#define SOLVESBMPROCESSOR2D_IMPLICIT_HH_INCLUDED

#include "solveSBMProcessor2DImplicit.h"

template<typename T>
solveSBM2DImplicit<T>::solveSBM2DImplicit(T eps2_, T Q_, T invh2_, T dt_, T M_, T h_, T cosTheta_, bool SM_, bool useG_):
                    eps2(eps2_), Q(Q_), invh2(invh2_), dt(dt_), M(M_), h(h_), cosTheta(cosTheta_), SM(SM_),
		    useG(useG_)	{ }

// Override generic process method
template<typename T>
void solveSBM2DImplicit<T>::processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D* > atomicBlocks) {

     // Casting back the blocks to their specific type
     ScalarField2D<T>& mu = *dynamic_cast<ScalarField2D<T>* >(atomicBlocks[0]);
     ScalarField2D<T>& C = *dynamic_cast<ScalarField2D<T>* >(atomicBlocks[1]);
     ScalarField2D<T>& Psi = *dynamic_cast<ScalarField2D<T>* >(atomicBlocks[2]);
     ScalarField2D<T>& M0 = *dynamic_cast<ScalarField2D<T>* >(atomicBlocks[3]);
     ScalarField2D<T>& Mu0 = *dynamic_cast<ScalarField2D<T>* >(atomicBlocks[4]);
     ScalarField2D<T>& C0 = *dynamic_cast<ScalarField2D<T>* >(atomicBlocks[5]);

    // Compute displacments
     Dot2D offset_C = computeRelativeDisplacement(mu, C);
     Dot2D offset_Psi = computeRelativeDisplacement(mu, Psi);
     Dot2D offset_M0 = computeRelativeDisplacement(mu, M0);
     Dot2D offset_Mu0 = computeRelativeDisplacement(mu, Mu0);
     Dot2D offset_C0 = computeRelativeDisplacement(mu, C0);

     // Loop over the indexes of the blocks
     for (plint iX=domain.x0; iX<=domain.x1; ++iX) {

         // X Displacements for block C and dPsiOvPsi
         plint iX_C = iX + offset_C.x;
         plint iX_Psi = iX + offset_Psi.x;
	 plint iX_M0 = iX + offset_M0.x;
	 plint iX_Mu0 = iX + offset_Mu0.x;
	 plint iX_C0 = iX + offset_C0.x;


         for (plint iY=domain.y0; iY<=domain.y1; ++iY) {

            // Y Displacements for block C and dPsiOvPsi
            plint iY_C = iY + offset_C.y;
            plint iY_Psi = iY + offset_Psi.y;
	    plint iY_M0 = iY + offset_M0.y;
	    plint iY_Mu0 = iY + offset_Mu0.y;
	    plint iY_C0 = iY + offset_C0.y;

		
		// compute invariant terms		
		// derivatives of Psi
		T dPsi_dx = (Psi.get(iX_Psi + 1, iY_Psi) - Psi.get(iX_Psi - 1, iY_Psi))/(2.*h);
		T dPsi_dy = (Psi.get(iX_Psi, iY_Psi - 1) - Psi.get(iX_Psi, iY_Psi - 1))/(2.*h);

		// compute abs of psi gradient
		T abs_psi = std::sqrt(dPsi_dx*dPsi_dx + dPsi_dy*dPsi_dy);

		// compute grad_psi_term
		T grad_psi_term = ((2.*std::sqrt(eps2)*abs_psi)*cosTheta)/Psi.get(iX_Psi, iY_Psi);

		// compute psi at half point
		T Psi_x_plus_half = (Psi.get(iX_Psi + 0.5, iY_Psi) + Psi.get(iX_Psi, iY_Psi))/2.; 
		T Psi_x_minus_half = (Psi.get(iX_Psi, iY_Psi) + Psi.get(iX_Psi - 0.5, iY_Psi))/2.;
		T Psi_y_plus_half = (Psi.get(iX_Psi, iY_Psi + 0.5) + Psi.get(iX_Psi, iY_Psi))/2.;
		T Psi_y_minus_half = (Psi.get(iX_Psi, iY_Psi) + Psi.get(iX_Psi, iY_Psi - 0.5))/2.;

		// compute psi sum term
		T sum_psi = ((2.*eps2)/(Psi.get(iX_Psi, iY_Psi)*h*h))*(Psi_x_plus_half +
			      Psi_x_minus_half + Psi_y_plus_half + Psi_y_minus_half);

		// compute M0 elements for computing M0 at half points
		T M0_x_y = M0.get(iX_M0, iY_M0);
                T M0_x_plus_y = M0.get(iX_M0 + 1, iY_M0);
                T M0_x_y_plus = M0.get(iX_M0, iY_M0 + 1);
                T M0_x_minus_y = M0.get(iX_M0 - 1, iY_M0);
                T M0_x_y_minus = M0.get(iX_M0, iY_M0 - 1);

		// compute M0_psi at haf points
		T M0_psi_x_plus_half = (M0_x_plus_y*Psi.get(iX_Psi + 1, iY_Psi) +
                                       M0_x_y*Psi.get(iX_Psi , iY_Psi))/2.;
                T M0_psi_x_minus_half = (M0_x_y*Psi.get(iX_Psi, iY_Psi) +
                                       M0_x_minus_y*Psi.get(iX_Psi - 1, iY_Psi))/2.;
                T M0_psi_y_plus_half = (M0_x_y_plus*Psi.get(iX_Psi, iY_Psi + 1) +
                                       M0_x_y*Psi.get(iX_Psi , iY_Psi))/2.;
                T M0_psi_y_minus_half = (M0_x_y*Psi.get(iX_Psi, iY_Psi) +
                                       M0_x_y_minus*Psi.get(iX_Psi, iY_Psi - 1))/2.;

		// compute sum_M0_psi_grad_Mu0
		T sum_M0_psi_gradMu = M0_psi_x_plus_half*(Mu0.get(iX_Mu0 + 1, iY_Mu0) - Mu0.get(iX_Mu0, iY_Mu0)) - 
				      M0_psi_x_minus_half*(Mu0.get(iX_Mu0, iY_Mu0) - Mu0.get(iX_Mu0 - 1, iY_Mu0)) +
				      M0_psi_y_plus_half*(Mu0.get(iX_Mu0, iY_Mu0 + 1) - Mu0.get(iX_Mu0, iY_Mu0)) -
				      M0_psi_y_minus_half*(Mu0.get(iX_Mu0, iY_Mu0) - Mu0.get(iX_Mu0, iY_Mu0 - 1));

		// compute time_n_terms
		T time_n_terms = C0.get(iX_C0, iY_C0) + (dt/(2.*Psi.get(iX_Psi, iY_Psi)))*sum_M0_psi_gradMu;
		
		// compute double well potential F
		T F = (Q/4.)*(C.get(iX_C, iY_C))*(C.get(iX_C, iY_C))*
                              (1 - C.get(iX_C, iY_C))*(1. - C.get(iX_C, iY_C));

		// compute derivative o F
		T dF_dC = (Q/2.)*C.get(iX_C, iY_C)*(1. - C.get(iX_C, iY_C))*(1. -2.*C.get(iX_C, iY_C));
		
		// compute second derivative of F (from Taylor  expansion)
		T d2F_dC2 = (Q/2.)*(1 - C.get(iX_C, iY_C))*(1 - C.get(iX_C, iY_C))*
			      (1. - 2.*C.get(iX_C, iY_C));

		// compute M elements for computing M at the half point
		T M_x_y = M*(C.get(iX_C, iY_C))*(C.get(iX_C, iY_C))*
                                (1. - C.get(iX_C, iY_C))*(1. - C.get(iX_C, iY_C));
		T M_x_plus_y = M*(C.get(iX_C + 1, iY_C))*(C.get(iX_C + 1, iY_C))*
                                (1. - C.get(iX_C + 1, iY_C))*(1. - C.get(iX_C + 1, iY_C));
		T M_x_y_plus = M*(C.get(iX_C, iY_C + 1))*(C.get(iX_C, iY_C + 1))*
                                (1. - C.get(iX_C, iY_C + 1))*(1. - C.get(iX_C, iY_C + 1));
		T M_x_minus_y = M*(C.get(iX_C - 1, iY_C))*(C.get(iX_C - 1, iY_C))*
                                (1. - C.get(iX_C - 1, iY_C))*(1. - C.get(iX_C - 1, iY_C));
		T M_x_y_minus = M*(C.get(iX_C, iY_C - 1))*(C.get(iX_C, iY_C - 1))*
                                (1 - C.get(iX_C, iY_C - 1))*(1 - C.get(iX_C, iY_C - 1));

		// compute M_psi at half points
		T M_psi_x_plus_half = (M_x_plus_y*Psi.get(iX_Psi + 1, iY_Psi) + 
				       M_x_y*Psi.get(iX_Psi , iY_Psi))/2.; 
		T M_psi_x_minus_half = (M_x_y*Psi.get(iX_Psi, iY_Psi) +
                                       M_x_minus_y*Psi.get(iX_Psi - 1, iY_Psi))/2.;
		T M_psi_y_plus_half = (M_x_y_plus*Psi.get(iX_Psi, iY_Psi + 1) +
                                       M_x_y*Psi.get(iX_Psi , iY_Psi))/2.;
		T M_psi_y_minus_half = (M_x_y*Psi.get(iX_Psi, iY_Psi) +
                                       M_x_y_minus*Psi.get(iX_Psi, iY_Psi - 1))/2.;
		
		// compute sum_M_psi
		T sum_M_psi = (dt/(2.*Psi.get(iX_Psi, iY_Psi)*h*h))*
			      (M_psi_x_plus_half + M_psi_x_minus_half +
			       M_psi_y_plus_half + M_psi_y_minus_half);

		// compute sum_M_psi_mu
		T sum_M_psi_mu = (dt/(2.*Psi.get(iX_Psi, iY_Psi)*h*h))*
				 (M_psi_x_plus_half*mu.get(iX + 1, iY) + M_psi_x_minus_half*mu.get(iX - 1, iY) +
				  M_psi_y_plus_half*mu.get(iX, iY + 1) + M_psi_y_minus_half*mu.get(iX, iY - 1));

		// compute psi at half at half points
		T psi_x_plus_half_y = (Psi.get(iX_Psi + 1, iY_Psi) + Psi.get(iX_Psi, iY_Psi))/2.;
		T psi_x_minus_half_y = (Psi.get(iX_Psi, iY_Psi) + Psi.get(iX_Psi - 1, iY_Psi))/2.;
		T psi_x_y_plus_half = (Psi.get(iX_Psi, iY_Psi + 1) + Psi.get(iX_Psi, iY_Psi))/2.;
		T psi_x_y_minus_half = (Psi.get(iX_Psi, iY_Psi) + Psi.get(iX_Psi, iY_Psi - 1))/2.;

		// compute sum_psi_C
		T sum_psi_C = ((2.*eps2)/(Psi.get(iX_Psi, iY_Psi)))*(psi_x_plus_half_y*C.get(iX_C + 1, iY_C) +
								     psi_x_minus_half_y*C.get(iX_C - 1, iY_C) + 
								     psi_x_y_plus_half*C.get(iX_C, iY_C + 1) + 
								     psi_x_y_minus_half*C.get(iX_C, iY_C - 1));

		// compute determinant of matrix A
		T det = 1 + sum_M_psi*d2F_dC2 + sum_M_psi*sum_psi;
		
		// compute A22*R1 = R1
		T A22_R1 = time_n_terms*sum_M0_psi_gradMu;		

		// compute A12*R2
		T A12_R2 = (- d2F_dC2 - sum_psi)*(dF_dC - d2F_dC2*C.get(iX_C, iY_C) -
						    grad_psi_term*std::sqrt(F) - sum_psi_C);

		// compute A21*R1
		T A21_R1 = (- d2F_dC2 - sum_psi)*(time_n_terms + sum_M_psi_mu);

		// compute A11*R2 = R2
		T A11_R2 = dF_dC - d2F_dC2*C.get(iX_C, iY_C) - 
			   grad_psi_term*std::sqrt(F) - sum_psi_C; 

		// update C and mu
		C.get(iX_C, iY_C) = (A22_R1 - A12_R2)/det;
		mu.get(iX, iY) = (A11_R2 - A21_R1)/det;
		
	}
    }
}

template<typename T>
solveSBM2DImplicit<T>*
        solveSBM2DImplicit<T>::clone() const {

    return new solveSBM2DImplicit<T>(*this);
}

 
template<typename T>
void solveSBM2DImplicit<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {

    modified[0] = modif::staticVariables;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
    modified[3] = modif::nothing;

}


#endif // SOLVESBMPROCESSOR2D_IMPLICIT_HH_INCLUDED
