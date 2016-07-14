// $Id$
//==============================================================================
//!
//! \file CubicMinimum.h
//!
//! \date Jul 13 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class finding the minimum of a cubic hermite interpolant.
//!
//==============================================================================

#include <cstdlib>
#include <vector>


class CubicMinimum {
public:
  static bool Find(double& alpha,
                   const std::vector<double>& params,
                   const std::vector<double>& vals,
                   const std::vector<double>& tgts);
};
