#ifndef DUNEURO_ANALYTIC_SOLUTION_HH
#define DUNEURO_ANALYTIC_SOLUTION_HH

#include <dune/common/fvector.hh>
#include <duneuro/common/dipole.hh>

namespace duneuro {

  // implements the analytic MEG forwad solution for multilayer sphere models in 3 dimensions
  // We assume layer wise isotropic conductivity
  template<class FieldType>
  class AnalyticSolutionMEG
  {
  public:
    static constexpr size_t dim = 3;
    using Coordinate = Dune::FieldVector<FieldType, dim>;
  
    // constructor
    AnalyticSolutionMEG(const Coordinate& sphereCenter, FieldType scalingFactor = 1.0)
      : sphereCenter_(sphereCenter)
      , scalingFactor_(scalingFactor)
    {
    }
  
    void bind(const Dipole<FieldType, dim>& dipole)
    {
      R_0 = dipole.position() - sphereCenter_;
      moment_ = dipole.moment();
    }
    
    //////////////////////////////////
    // we define methods to compute the total B-field, the primary B-field and the secondary B-field
    //////////////////////////////////
  
    // compute analytical solution of the total magnetic Field as described in 
    // Basic mathematical and electromagnetic concepts of the biomagnetic inverse problem, Jukka Sarvas, 1987, §4
    Coordinate totalField(const Coordinate& coilPos)
    {
      Coordinate R = coilPos - sphereCenter_;
      Coordinate A = R - R_0;
      FieldType r = R.two_norm();
      FieldType a = A.two_norm();
      
      FieldType F = a * (r * a + r * r - R_0 * R);
      
      Coordinate grad_F = (a*a/r + A * R/a + 2*(a + r)) * R - (a + 2*r + A * R/a)* R_0; 
      
      return scalingFactor_ * (F * crossProduct(moment_, R_0) - (crossProduct(moment_, R_0) * R) * grad_F) / (F * F);
    }
    
    FieldType totalField(const Coordinate& coilPos, const Coordinate& direction)
    {
      return totalField(coilPos) * direction;
    }
    
    // compute primary field
    Coordinate primaryField(const Coordinate& coilPos)
    {
      Coordinate R = coilPos - sphereCenter_;
      Coordinate diff = R - R_0;
      FieldType diffNorm = diff.two_norm();
      diff /= (diffNorm * diffNorm * diffNorm);
      return scalingFactor_ * crossProduct(moment_, diff);
    }
    
    FieldType primaryField(const Coordinate& coilPos, const Coordinate& direction)
    {
      return primaryField(coilPos) * direction;
    }
    
    // compute secondary field
    Coordinate secondaryField(const Coordinate& coilPos)
    {
      return primaryField(coilPos) - totalField(coilPos);
    }
    
    FieldType secondaryField(const Coordinate& coilPos, const Coordinate& direction)
    {
      return primaryField(coilPos, direction) - totalField(coilPos, direction);
    }
    
  private:
    //set in constructor
    Coordinate sphereCenter_;
    FieldType scalingFactor_;
    
    // set later on
    Coordinate moment_;
    Coordinate R_0;
    
    // compute the cross prodcut of two vectors. Copied from dune/pdelab/common/crossproduct.hh
    Coordinate crossProduct(const Coordinate& vec_1, const Coordinate& vec_2)
    {
      Coordinate crossProd;
      for(size_t i = 0; i < 3; ++i) {
        size_t j = (i + 1) % 3;
        size_t k = (i + 2) % 3;
        
        crossProd[i] = vec_1[j] * vec_2[k] - vec_1[k] * vec_2[j];
      }
      
      return crossProd;
    }
  }; // end class AnalyticSolutionMEG 

} // end namespace duneuro
#endif // DUNEURO_ANALYTIC_SOLUTION_HH
