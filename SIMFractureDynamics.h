// $Id$
//==============================================================================
//!
//! \file SIMFractureDynamics.h
//!
//! \date Jul 13 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Driver for fracture-dynamic problems.
//!
//==============================================================================

#ifndef _SIM_FRACTURE_DYNAMICS_H_
#define _SIM_FRACTURE_DYNAMICS_H_

#include "ProcessAdm.h"
#ifdef HAS_LRSPLINE
#include "ASMu2D.h"
#include "LRSpline/LRSplineSurface.h"
#endif
#include <fstream>
#include <numeric>


/*!
  \brief Driver class for fracture dynamics simulators.
  \details A fracture dynamics simulator is a coupling between
  a dynamic elasticity solver and a phase field solver.
*/

template<class SolidSolver, class PhaseSolver,
         template<class S1, class S2> class Coupling>
class SIMFracture : public Coupling<SolidSolver,PhaseSolver>
{
  //! Convenience type
  typedef Coupling<SolidSolver,PhaseSolver> CoupledSIM;

public:
  //! \brief The constructor initializes the references to the two solvers.
  SIMFracture(SolidSolver& s1, PhaseSolver& s2, const std::string& inputfile)
    : CoupledSIM(s1,s2), infile(inputfile), aMin(0.0)
  {
    doStop = false;
    irfStop = 0;
    E0 = Ec = Ep = stopVal = 0.0;
  }

  //! \brief Empty destructor.
  virtual ~SIMFracture() {}

  //! \brief Initializes and sets up field dependencies.
  virtual void setupDependencies()
  {
    this->S1.registerDependency(&this->S2,"phasefield",1);
    // The tensile energy is defined on integration points and not nodal points.
    // It is a global buffer array across all patches in the model.
    // Use an explicit call instead of normal couplings for this.
    this->S2.setTensileEnergy(this->S1.getTensileEnergy());
  }

  //! \brief Advances the time step one step forward.
  virtual bool advanceStep(TimeStep& tp)
  {
    return this->CoupledSIM::advanceStep(tp) && !doStop;
  }

  //! \brief Computes the solution for the current time step.
  virtual bool solveStep(TimeStep& tp, bool firstS1 = true)
  {
    if (tp.step == 1 && this->S1.haveCrackPressure())
      // Start the initial step by solving the phase-field first
      if (!this->S2.solveStep(tp,false))
        return false;

    return this->CoupledSIM::solveStep(tp,firstS1);
  }

  //! \brief Saves the converged results to VTF-file of a given time step.
  //! \details It also writes global energy quantities to file for plotting.
  virtual bool saveStep(const TimeStep& tp, int& nBlock)
  {
    Vector RF;
    this->S1.getBoundaryReactions(RF);

    if (!energFile.empty() && tp.step > 0 &&
        this->S1.getProcessAdm().getProcId() == 0)
    {
      std::ofstream os(energFile, tp.step == 1 ? std::ios::out : std::ios::app);

      Vector BF;
      this->S1.getBoundaryForce(BF,this->S1.getSolutions(),tp);

      if (tp.step == 1)
      {
        size_t i;
        os <<"#t eps_e external_energy eps+ eps- eps_b |c|"
           <<" eps_d-eps_d(0) eps_d";
        for (i = 0; i < BF.size(); i++)
          os <<" load_"<< char('X'+i);
        for (i = 0; i < RF.size(); i++)
          os <<" react_"<< char('X'+i);
        os << std::endl;
      }

      os << std::setprecision(11) << std::setw(6) << std::scientific
         << tp.time.t;
      for (double n1 : this->S1.getGlobalNorms()) os <<" "<< n1;
      const Vector& n2 = this->S2.getGlobalNorms();
      os <<" "<< (n2.size() > 2 ? n2[1] : 0.0);
      os <<" "<< (n2.size() > 1 ? n2[n2.size()-2] : 0.0);
      os <<" "<< (n2.size() > 0 ? n2.back() : 0.0);
      for (double f : BF) os <<" "<< utl::trunc(f);
      for (double f : RF) os <<" "<< utl::trunc(f);
      os << std::endl;
    }

    // Check stop criterion
    if (tp.step > 1 && irfStop > 0 && irfStop <= RF.size())
      if ((doStop = fabs(RF(irfStop)) < stopVal))
        IFEM::cout <<"\n >>> Terminating simulation due to stop criterion |RF("
                   << irfStop <<")| = "<< fabs(RF(irfStop)) <<" < "<< stopVal
                   << std::endl;

    return (this->S2.saveStep(tp,nBlock) && this->S1.saveStep(tp,nBlock) &&
            this->S2.saveResidual(tp,residual,nBlock));
  }

  //! \brief Parses staggering parameters from an XML element.
  virtual void parseStaggering(const TiXmlElement* elem)
  {
    const TiXmlElement* child = elem->FirstChildElement("stop");
    if (child)
    {
      utl::getAttribute(child,"rcomp",irfStop);
      utl::getAttribute(child,"force",stopVal);
    }
  }

  //! \brief Assigns the file name for global energy output.
  void setEnergyFile(const char* fName)
  {
    if (fName)
    {
      energFile = fName;
      IFEM::cout <<"\tFile for global energy output: "<< energFile << std::endl;
    }
  }

  //! \brief Stores current solution state in an internal buffer.
  void saveState()
  {
    sols = this->S1.getSolutions();
    sols.push_back(this->S2.getSolution());
    hsol = this->S2.getHistoryField();
  }

  //! \brief Refines the mesh on the initial configuration.
  bool initialRefine(double beta, double min_frac, int nrefinements)
  {
    if (this->S2.getInitRefine() >= nrefinements)
      return true; // Grid is sufficiently refined during input parsing
    else if (this->S2.hasIC("phasefield"))
      return true; // No initial refinement when specified initial phase field

    TimeStep step0;
    int newElements = 1;
    for (step0.iter = 0; newElements > 0; step0.iter++)
      if (!this->S2.solveStep(step0))
        return false;
      else
        newElements = this->adaptMesh(beta,min_frac,nrefinements);

    return newElements == 0;
  }

  //! \brief Refines the mesh with transfer of solution onto the new mesh.
  int adaptMesh(double beta, double min_frac, int nrefinements)
  {
#ifdef HAS_LRSPLINE
    // TODO: Add multi-patch support
    ASMu2D* pch = dynamic_cast<ASMu2D*>(this->S1.getPatch(1));
    if (!pch) return -999; // Logic error, should not happen

    if (aMin <= 0.0) // maximum refinements per element
    {
      double redMax = pow(2.0,nrefinements);
      aMin = pch->getBasis()->getElement(0)->area()/(redMax*redMax);
    }

    // Fetch element norms to use as refinement criteria
    Vector eNorm;
    double gNorm = this->S2.getNorm(eNorm,3);
    if (eNorm.empty())
    {
      std::cerr <<" *** SIMFractureDynamics:adaptMesh: Missing refinement"
                <<" indicators, expected as the 3rd element norm."<< std::endl;
      return -1;
    }

    // Sort element indices based on comparing values in eNorm
    IntVec idx(eNorm.size());
    std::iota(idx.begin(),idx.end(),0);
    std::sort(idx.begin(),idx.end(),
              [&eNorm](size_t i1, size_t i2) { return eNorm[i1] < eNorm[i2]; });

    double eMin = min_frac < 0.0 ? -min_frac*gNorm/sqrt(idx.size()) : min_frac;
    size_t eMax = beta < 0.0 ? idx.size() : idx.size()*beta/100.0;
    IFEM::cout <<"\n  Lowest element: "<< std::setw(8) << idx.front()
               <<"    |c| = "<< eNorm[idx.front()]
               <<"\n  Highest element:"<< std::setw(8) << idx.back()
               <<"    |c| = "<< eNorm[idx.back()]
               <<"\n  Minimum |c|-value for refinement: "<< eMin
               <<"\n  Minimum element area: "<< aMin << std::endl;

    IntVec elements; // Find the elements to refine
    elements.reserve(eMax);
    for (int eid : idx)
      if (eNorm[eid] > eMin || elements.size() >= eMax)
        break;
      else if (pch->getBasis()->getElement(eid)->area() > aMin+1.0e-12)
        elements.push_back(eid);

    if (elements.empty())
      return 0;

    IFEM::cout <<"  Elements to refine: "<< elements.size()
               <<" (|c| = ["<< eNorm[elements.front()]
               <<","<< eNorm[elements.back()] <<"])\n"<< std::endl;

    LR::LRSplineSurface* oldBasis = nullptr;
    if (!hsol.empty()) oldBasis = pch->getBasis()->copy();

    // Do the mesh refinement
    LR::RefineData prm;
    prm.options = { 10, 1, 2, 0, 1 };
    prm.elements = pch->getFunctionsForElements(elements);
    if (!this->S1.refine(prm,sols) || !this->S2.refine(prm))
      return -2;

    // Re-initialize the simulators for the new mesh
    this->S1.clearProperties();
    this->S2.clearProperties();
    if (!this->S1.read(infile.c_str()) || !this->S2.read(infile.c_str()))
      return -3;

    if (!this->preprocess())
      return -4;

    if (!this->init(TimeStep()))
      return -5;

    if (!this->S1.initSystem(this->S1.opt.solver,1,1,0,true) ||
        !this->S2.initSystem(this->S2.opt.solver))
      return -6;

    // Transfer solution variables onto the new mesh
    if (!sols.empty())
    {
      IFEM::cout <<"\nTransferring "<< sols.size()-1 <<"x"<< sols.front().size()
                 <<" solution variables to new mesh for "<< this->S1.getName();
      this->S1.setSolutions(sols);
      IFEM::cout <<"\nTransferring "<< sols.back().size()
                 <<" solution variables to new mesh for "<< this->S2.getName();
      this->S2.setSolution(sols.back());
    }
    if (!hsol.empty())
    {
      IFEM::cout <<"\nTransferring "<< hsol.size()
                 <<" history variables to new mesh for "<< this->S2.getName()
                 << std::endl;
      this->S2.transferHistory2D(hsol,oldBasis);
      delete oldBasis;
    }

    return elements.size();
#else
    std::cerr <<" *** SIMFractureDynamics:adaptMesh: No LR-spline support.\n";
    return -1;
#endif
  }

  //! \brief Dumps the current grid to specified file.
  bool dumpMesh(const char* fileName)
  {
    std::ofstream os(fileName);
    return this->S2.dumpGeometry(os);
  }

protected:
  //! \brief Calculates and prints the solution and residual norms.
  double calcResidual(const TimeStep& tp, bool cycles = false)
  {
    // Compute residual of the elasticity equation
    this->S1.setMode(SIM::RHS_ONLY);
    if (!this->S1.assembleSystem(tp.time,this->S1.getSolutions(),false))
      return -1.0;

    if (!this->S1.extractLoadVec(residual))
      return -1.0;

    double rNorm1 = residual.norm2();
    double eNorm1 = this->S1.extractScalar();

    // Compute residual of the phase-field equation
    if (!this->S2.setMode(SIM::INT_FORCES))
      return -2.0;

    Vectors sol2(1,this->S2.getSolution());
    if (!this->S2.assembleSystem(tp.time,sol2,false))
      return -2.0;

    if (!this->S2.extractLoadVec(residual))
      return -2.0;

    double rNorm2 = residual.norm2();
    double eNorm2 = this->S2.extractScalar();

    double rConv = rNorm1 + rNorm2;
    double eConv = eNorm1 + eNorm2;
    if (cycles)
    {
      IFEM::cout <<"  cycle "<< tp.iter
                 <<": Res = "<< rNorm1 <<" + "<< rNorm2 <<" = "<< rConv;
      if (eConv > 0.0)
        IFEM::cout <<"  E = "<< eNorm1 <<" + "<< eNorm2 <<" = "<< eConv;
      if (tp.iter == 0)
        E0 = eConv;
      else
      {
        Ep = tp.iter > 1 ? Ec : E0;
        Ec = eConv;
        if (eConv > 0.0)
          IFEM::cout <<"  beta="<< atan2(tp.iter*(Ep-Ec),E0-Ec) * 180.0/M_PI;
      }
    }
    else
    {
      IFEM::cout <<"  Res = "<< rNorm1 <<" + "<< rNorm2 <<" = "<< rConv;
      if (eConv > 0)
        IFEM::cout <<"\n    E = "<< eNorm1 <<" + "<< eNorm2 <<" = "<< eConv;
    }
    IFEM::cout << std::endl;
    return rConv;
  }

private:
  std::string energFile; //!< File name for global energy output
  std::string infile;    //!< Input file parsed

  double    aMin; //!< Minimum element area
  Vectors   sols; //!< Solution state to transfer onto refined mesh
  RealArray hsol; //!< History field to transfer onto refined mesh

  size_t irfStop; //!< Reaction force component to use as stop criterion
  double stopVal; //!< Stop simulation when less that this value
  bool   doStop;  //!< If \e true, terminate simulation due to stop criterion

  double E0; //!< Energy norm of initial staggering cycle
  double Ec; //!< Energy norm of current staggering cycle
  double Ep; //!< Energy norm of previous staggering cycle

  Vector residual; //!< Residual force vector (of the phase field equation)
};

#endif
