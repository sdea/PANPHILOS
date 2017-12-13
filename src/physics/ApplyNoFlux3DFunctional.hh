#ifndef APPLYNOFLUX3DFUNCTIONAL_HH_INCLUDED
#define APPLYNOFLUX3DFUNCTIONAL_HH_INCLUDED

#include "ApplyNoFlux3DFunctional.h"
namespace plb {

template<typename T>
applyNoFlux3D<T>::applyNoFlux3D(plint dir_, plint normal_):
                                    dir(dir_), normal(normal_) { }

template<typename T>
void applyNoFlux3D<T>::process(Box3D domain, ScalarField3D<T>& field) {

    for(plint iX=domain.x0; iX <= domain.x1; ++iX) {
            for(plint iY=domain.y0; iY <= domain.y1; ++iY) {
                for(plint iZ=domain.z0; iZ <= domain.z1; ++iZ) {

                switch (dir) {
                case 0:

                    // Boundary conditions on X
                    if(normal == -1) {

                        field.get(domain.x0, iY, iZ) = field.get(domain.x1, iY, iZ);
                    }
                    else {
                        field.get(domain.x1, iY,iZ) = field.get(domain.x0, iY, iZ);
                    }
                break;

                case 1:

                    // Boundary conditions on Y
                    if(normal == -1) {

                        field.get(iX, domain.y0, iZ) = field.get(iX, domain.y1, iZ);
                    }
                    else {
                        field.get(iX, domain.y1, iZ) = field.get(iX, domain.y0, iZ);
                    }
                break;

                case 2:
                    // Boundary conditions on Z
                    if(normal == -1) {

                        field.get(iX, iY, domain.z0) = field.get(iX, iY, domain.z1);
                    }
                    else {
                        field.get(iX, iY, domain.z1) = field.get(iX, iY, domain.z0);
                    }
                break;

                default:
                break;


                } // end switch
            }
          }
        } // end for
}

template<typename T>
applyNoFlux3D<T>* applyNoFlux3D<T>::clone() const {

    return new applyNoFlux3D<T>(*this);
}

template<typename T>
void applyNoFlux3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {

     modified[0] = modif::allVariables;

}

}



#endif // APPLYNOFLUX3DFUNCTIONAL_HH_INCLUDED
