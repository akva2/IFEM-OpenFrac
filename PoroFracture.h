// $Id$
//==============================================================================
//!
//! \file PoroFracture.h
//!
//! \date Apr 15 2016
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand implementations for poroelasticity problems with fracture.
//!
//==============================================================================

#ifndef _PORO_FRACTURE_H
#define _PORO_FRACTURE_H

#include "PoroElasticity.h"

class FractureElasticity;


/*!
  \brief Class representing the integrand of poroelasticity with fracture.
  \details This class inherits PoroElasticity and uses elements from
  FractureElasticity in addition through a private member.
*/

class PoroFracture : public PoroElasticity
{
public:
  //! \brief The constructor allocates the internal FractureElasticy object.
  //! \param[in] n Number of spatial dimensions
  PoroFracture(unsigned short int n);
  //! \brief The destructor deletes the internal FractureElasticy object.
  virtual ~PoroFracture();

  //! \brief Parses a data section from an XML-element.
  virtual bool parse(const TiXmlElement* elem);

  using PoroElasticity::parseMatProp;
  //! \brief Parses material properties from an XML-element.
  virtual Material* parseMatProp(const TiXmlElement* elem, bool);

  //! Defines the material properties.
  virtual void setMaterial(Material* mat);

  //! \brief Defines the solution mode before the element assembly is started.
  //! \param[in] mode The solution mode to use
  virtual void setMode(SIM::SolutionMode mode);

  //! \brief Initializes the integrand with the number of integration points.
  //! \param[in] nGp Total number of interior integration points
  //! \param[in] nBp Total number of boundary integration points
  virtual void initIntegration(size_t nGp, size_t nBp);

  //! \brief Returns a local integral contribution object for the given element.
  //! \param[in] nen Number of nodes on element
  //! \param[in] neumann Whether or not we are assembling Neumann BCs
  virtual LocalIntegral* getLocalIntegral(size_t nen,
                                          size_t, bool neumann) const;
  //! \brief Returns a local integral contribution object for the given element.
  //! \param[in] nen Number of nodes on element for each basis
  //! \param[in] neumann Whether or not we are assembling Neumann BCs
  virtual LocalIntegral* getLocalIntegral(const std::vector<size_t>& nen,
                                          size_t, bool neumann) const;

  using PoroElasticity::initElement;
  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param elmInt The local integral object for current element
  virtual bool initElement(const std::vector<int>& MNPC, LocalIntegral& elmInt);
  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param[in] elem_sizes Size of each basis on the element
  //! \param[in] basis_sizes Size of each basis on the patch level
  //! \param elmInt The local integral object for current element
  virtual bool initElement(const std::vector<int>& MNPC,
                           const std::vector<size_t>& elem_sizes,
                           const std::vector<size_t>& basis_sizes,
                           LocalIntegral& elmInt);

  //! \brief Returns a pointer to the Gauss-point tensile energy array.
  virtual const RealArray* getTensileEnergy() const;

protected:
  //! \brief Computes the elasticity matrices at a quadrature point.
  //! \param elMat The element matrix object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalElasticityMatrices(ElmMats& elMat, const Matrix&,
                                      const FiniteElement& fe,
                                      const Vec3& X) const;

private:
  FractureElasticity* fracEl; //!< Integrand for tangent stiffness evaluation
};

#endif
