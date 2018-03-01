#ifndef APPLYPERIODICBC2DFUNCTIONAL_H_INCLUDED
#define APPLYPERIODICBC2DFUNCTIONAL_H_INCLUDED


namespace plb {

template<typename T>
class applyPeriodicBC2D: public BoxProcessingFunctional2D_S<T> {

public:

        applyPeriodicBC2D(plint dir_, plint normal_);
        virtual applyPeriodicBC2D<T>* clone() const;
        virtual void process(Box2D domain, ScalarField2D<T>& field);
        virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:

    plint dir;
    plint normal;

};

}

#endif // APPLYNOFLUXANDPERIODIC2D_H_INCLUDED
