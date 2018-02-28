#ifndef APPLYNOFLUXANDPERIODIC2D_H_INCLUDED
#define APPLYNOFLUXANDPERIODIC2D_H_INCLUDED


namespace plb {

template<typename T>
class applyNoFluxAndPeriodic2D: public BoxProcessingFunctional2D_S<T> {

public:

        applyNoFluxAndPeriodic2D(plint dir_, plint normal_);
        virtual applyNoFluxAndPeriodic2D<T>* clone() const;
        virtual void process(Box2D domain, ScalarField2D<T>& field);
        virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:

    plint dir;
    plint normal;

};

}

#endif // APPLYNOFLUXANDPERIODIC2D_H_INCLUDED
