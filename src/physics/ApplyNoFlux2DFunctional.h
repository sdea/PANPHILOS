#ifndef APPLYNOFLUX2DFUNCTIONAL_H_INCLUDED
#define APPLYNOFLUX2DFUNCTIONAL_H_INCLUDED


namespace plb {

template<typename T>
class applyNoFlux2D: public BoxProcessingFunctional2D_S<T> {

public:

        applyNoFlux2D(plint dir_, plint normal_);
        virtual applyNoFlux2D<T>* clone() const;
        virtual void process(Box2D domain, ScalarField2D<T>& field);
        virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:

    plint dir;
    plint normal;

};

}

#endif // APPLYNOFLUX2DFUNCTIONAL_H_INCLUDED
