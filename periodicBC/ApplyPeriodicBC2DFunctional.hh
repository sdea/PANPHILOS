#ifndef APPLYPERIODIC2DBCFUNCTIONAL_HH_INCLUDED
#define APPLYPERIODIC2DBCFUNCTIONAL_HH_INCLUDED

#include "./ApplyPeriodicBC2DFunctional.h"
namespace plb {

template<typename T>
applyPeriodicBC2D<T>::applyPeriodicBC2D(plint dir_, plint normal_):
                                    dir(dir_), normal(normal_) { }

template<typename T>
void applyPeriodicBC2D<T>::process(Box2D domain, ScalarField2D<T>& field) {

    for(plint iX=domain.x0; iX <= domain.x1; ++iX) {
            for(plint iY=domain.y0; iY <= domain.y1; ++iY) {

                switch (dir) {
                case 0:
		      
			if (normal == -1) {
	            
		      	   // Domain MUST BE leftPeriodic, when the normal
			   // is -1, left boundary 
		           field.get(domain.x0, iY) = field.get(domain.x1, iY);
			}
			else {

		           // Domain MUST BE rightPeriodic, when the normal
			   // is +1, right boundary
	                   field.get(domain.x1, iY) = field.get(domain.x0, iY);

			}

                    
                break;
		case 1:

		    if (normal == -1) {
	            
		      	   // Domain MUST BE downPeriodic, same as above but in y 
		           field.get(iX, domain.y0) = field.get(iX, domain.y1);
			}
			else {

		           // Domain MUST BE upPeriodic, same as above but in y
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
applyPeriodicBC2D<T>* applyPeriodicBC2D<T>::clone() const {

    return new applyPeriodicBC2D<T>(*this);
}

template<typename T>
void applyPeriodicBC2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {

     modified[0] = modif::allVariables;

}

}



#endif // APPLYNOFLUXANDPERIODIC2D_HH_INCLUDED
