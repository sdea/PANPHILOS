#ifndef APPLYNOFLUX3DFUNCTIONAL_H_INCLUDED
#define APPLYNOFLUX3DFUNCTIONAL_H_INCLUDED


namespace plb {

template<typename T>
class applyNoFlux3D: public BoxProcessingFunctional3D_S<T> {

public:

        applyNoFlux3D(plint dir_, plint normal_);
        virtual applyNoFlux3D<T>* clone() const;
        virtual void process(Box3D domain, ScalarField3D<T>& field);
        virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:

    plint dir;
    plint normal;

};

}

#endif // APPLYNOFLUX3DFUNCTIONAL_H_INCLUDED
