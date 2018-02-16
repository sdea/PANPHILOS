#ifndef COMPUTEGRADCOVMODGRADC_H_INCLUDED
#define COMPUTEGRADCOVMODGRADC_H_INCLUDED

namespace plb {

template<typename T1, typename T2, int nDim>
class ComputeGradCOvModGradC : public BoxProcessingFunctional3D_ST<T1, T2, nDim> {

public:

    ComputeGradCOvModGradC(T h_);
    virtual void process(Box3D domain, ScalarField3D<T1>& C, TensorField3D<T2,nDim>& GradCOvModGradC);
    virtual ComputeGradCOvModGradC<T1, T2, nDim> * clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;

private:

   T h;
};

}
#endif // COMPUTEGRADCOVMODGRADC_H_INCLUDED
