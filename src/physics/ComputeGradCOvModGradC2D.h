#ifndef COMPUTEGRADCOVMODGRADC2D_H_INCLUDED
#define COMPUTEGRADCOVMODGRADC2D_H_INCLUDED

namespace plb {

template<typename T1, typename T2, int nDim>
class ComputeGradCOvModGradC2D : public BoxProcessingFunctional2D_ST<T1, T2, nDim> {

public:

    ComputeGradCOvModGradC2D(T h_);
    virtual void process(Box2D domain, ScalarField2D<T1>& C, TensorField2D<T2,nDim>& GradCOvModGradC);
    virtual ComputeGradCOvModGradC2D<T1, T2, nDim> * clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;

private:

   T h;
};

}
#endif // COMPUTEGRADCOVMODGRADC2D_H_INCLUDED
