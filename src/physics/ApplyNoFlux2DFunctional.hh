#ifndef APPLYNOFLUX2DFUNCTIONAL_HH_INCLUDED
#define APPLYNOFLUX2DFUNCTIONAL_HH_INCLUDED

#include "ApplyNoFlux2DFunctional.h"
namespace plb {

template<typename T>
applyNoFlux2D<T>::applyNoFlux2D(plint dir_, plint normal_):
                                    dir(dir_), normal(normal_) { }

template<typename T>
void applyNoFlux2D<T>::process(Box2D domain, ScalarField2D<T>& field) {

    for(plint iX=domain.x0; iX <= domain.x1; ++iX) {
            for(plint iY=domain.y0; iY <= domain.y1; ++iY) {

                switch (dir) {
                case 0:

                    // Boundary conditions on X
                    if(normal == -1) {

                        field.get(domain.x0, iY) = field.get(domain.x1, iY);
                    }
                    else {
                        field.get(domain.x1, iY) = field.get(domain.x0, iY);
                    }
                break;

                case 1:

                    // Boundary conditions on Y
                    if(normal == -1) {

                        field.get(iX, domain.y0) = field.get(iX, domain.y1);
                    }
                    else {
                        field.get(iX, domain.y1) = field.get(iX, domain.y0);
                    }
                break;

                
                default:
                break;


                } // end switch

          }
        } // end for
}

template<typename T>
applyNoFlux2D<T>* applyNoFlux2D<T>::clone() const {

    return new applyNoFlux2D<T>(*this);
}

template<typename T>
void applyNoFlux2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {

     modified[0] = modif::allVariables;

}

}



#endif // APPLYNOFLUX2DFUNCTIONAL_HH_INCLUDED
