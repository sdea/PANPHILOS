#ifndef APPLYNOFLUXANDPERIODIC2D_HH_INCLUDED
#define APPLYNOFLUXANDPERIODIC2D_HH_INCLUDED

#include "ApplyNoFluxAndPeriodic2D.h"
namespace plb {

template<typename T>
applyNoFluxAndPeriodic2D<T>::applyNoFluxAndPeriodic2D(plint dir_, plint normal_):
                                    dir(dir_), normal(normal_) { }

template<typename T>
void applyNoFluxAndPeriodic2D<T>::process(Box2D domain, ScalarField2D<T>& field) {

    for(plint iX=domain.x0; iX <= domain.x1; ++iX) {
            for(plint iY=domain.y0; iY <= domain.y1; ++iY) {

                switch (dir) {
                case 0:
	            
		      // Here we apply periodic boundary conditions for x
		      field.get(domain.x0, iY) = field.get(domain.x1 -1, iY);
		      field.get(domain.x1, iY) = field.get(domain.x0 +1, iY);
                    }
                break;

                case 1:

                    // Here we apply periodic boudary conditions for y
		    field.get(iX, domain.y0) = field.get(iX, domain.y1 -1);
		    field.get(iX, domain.y1) = field.get(iX, domain.y0 +1);
 
		    
                break;

                
                default:
                break;


                } // end switch

          }
        } // end for
}

template<typename T>
applyNoFluxAndPeriodic2D<T>* applyNoFluxAndPeriodic2D<T>::clone() const {

    return new applyNoFluxAndPeriodic2D<T>(*this);
}

template<typename T>
void applyNoFluxAndPeriodic2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {

     modified[0] = modif::allVariables;

}

}



#endif // APPLYNOFLUXANDPERIODIC2D_HH_INCLUDED
