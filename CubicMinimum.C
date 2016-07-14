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

#include "CubicMinimum.h"
#include <GoTools/geometry/HermiteInterpolator.h>
#include <GoTools/geometry/SplineCurve.h>
#include <map>


bool CubicMinimum::Find(double& alpha,
                        const std::vector<double>& params,
                        const std::vector<double>& vals,
                        const std::vector<double>& tgts)
{
  if (vals.empty() || tgts.size() != vals.size())
    return false;

  // stick data in gotools structure
  std::vector<Go::Point> samples;
  Go::Point val(1);
  for (size_t i = 0; i < vals.size(); ++i) {
    val[0] = vals[i];
    samples.push_back(val);
    val[0] = tgts[i];
    samples.push_back(val);
  }

  // interpolate
  Go::HermiteInterpolator interp;
  std::vector<double> coefs;
  interp.interpolate(samples, params, coefs);

  // create curve
  Go::SplineCurve crv(interp.basis(), coefs.begin(), 1);

  // grab tangent curve
  std::unique_ptr<Go::SplineCurve> dcrv(crv.derivCurve(1));

  // for each knotspan find point closest to 0 and mark as an extremum if close enough.
  Go::Point nullpt(1);
  nullpt[0] = 0.0;
  std::vector<double> extrema;
  for (size_t i = 0; i < vals.size()-1; ++i) {
    Go::Point a, b;
    Go::Point loc_pt;
    double loc_alpha, loc_dist;
    dcrv->closestPoint(nullpt, params[i], params[i+1], loc_alpha, loc_pt, loc_dist);
    if ((loc_dist < 1e-5))
      extrema.push_back(loc_alpha);
  }

  // no extrema found
  if (extrema.empty())
    return false;

  // choose acceptable solution among values
  std::unique_ptr<Go::SplineCurve> ddcrv(crv.derivCurve(2));
  alpha = 0;
  double minVal = 1e100;
  for (auto& it : extrema) {
    Go::Point pt_dd, pt;
    crv.point(pt, it);
    ddcrv->point(pt_dd, it);
    if (pt_dd[0] >= 0.0 && pt[0] < minVal) {
      alpha = it;
      minVal = pt[0];
    }
  }

  return true;
}
