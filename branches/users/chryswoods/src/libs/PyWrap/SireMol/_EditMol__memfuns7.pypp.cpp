// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "_EditMol__memfuns7.pypp.hpp"
#include "boost/python.hpp"
#include "siremol_headers.h"
#include "SireMol/moleculedata.h"
#include "SireMaths/angle.h"
#include "SireMaths/quaternion.h"
#include "SireMaths/matrix.h"
#include "SireMaths/triangle.h"
#include "SireMaths/line.h"
#include "SireMaths/torsion.h"

namespace bp = boost::python;

void register_EditMol_memfuns7( EditMol_exposer_t& EditMol_exposer ){

    { //::SireMol::EditMol::measure
    
        typedef double ( ::SireMol::EditMol::*measure_function_type )( ::SireMol::Bond const & ) const;
        
        EditMol_exposer.def( 
            "measure"
            , measure_function_type( &::SireMol::EditMol::measure )
            , ( bp::arg("bnd") ) );
    
    }

    { //::SireMol::EditMol::measure
    
        typedef ::SireMaths::Angle ( ::SireMol::EditMol::*measure_function_type )( ::SireMol::Angle const & ) const;
        
        EditMol_exposer.def( 
            "measure"
            , measure_function_type( &::SireMol::EditMol::measure )
            , ( bp::arg("ang") ) );
    
    }

    { //::SireMol::EditMol::measure
    
        typedef ::SireMaths::Angle ( ::SireMol::EditMol::*measure_function_type )( ::SireMol::Dihedral const & ) const;
        
        EditMol_exposer.def( 
            "measure"
            , measure_function_type( &::SireMol::EditMol::measure )
            , ( bp::arg("dih") ) );
    
    }

    { //::SireMol::EditMol::measure
    
        typedef ::SireMaths::Angle ( ::SireMol::EditMol::*measure_function_type )( ::SireMol::Improper const & ) const;
        
        EditMol_exposer.def( 
            "measure"
            , measure_function_type( &::SireMol::EditMol::measure )
            , ( bp::arg("improper") ) );
    
    }

    { //::SireMol::EditMol::nAtoms
    
        typedef int ( ::SireMol::EditMol::*nAtoms_function_type )(  ) const;
        
        EditMol_exposer.def( 
            "nAtoms"
            , nAtoms_function_type( &::SireMol::EditMol::nAtoms ) );
    
    }

    { //::SireMol::EditMol::nAtoms
    
        typedef int ( ::SireMol::EditMol::*nAtoms_function_type )( ::SireMol::ResNum ) const;
        
        EditMol_exposer.def( 
            "nAtoms"
            , nAtoms_function_type( &::SireMol::EditMol::nAtoms )
            , ( bp::arg("resnum") ) );
    
    }

    { //::SireMol::EditMol::nAtoms
    
        typedef int ( ::SireMol::EditMol::*nAtoms_function_type )( ::SireMol::ResID ) const;
        
        EditMol_exposer.def( 
            "nAtoms"
            , nAtoms_function_type( &::SireMol::EditMol::nAtoms )
            , ( bp::arg("resid") ) );
    
    }

    { //::SireMol::EditMol::nAtoms
    
        typedef int ( ::SireMol::EditMol::*nAtoms_function_type )( ::SireMol::CutGroupID ) const;
        
        EditMol_exposer.def( 
            "nAtoms"
            , nAtoms_function_type( &::SireMol::EditMol::nAtoms )
            , ( bp::arg("id") ) );
    
    }

    { //::SireMol::EditMol::nAtoms
    
        typedef int ( ::SireMol::EditMol::*nAtoms_function_type )( ::SireMol::CutGroupNum ) const;
        
        EditMol_exposer.def( 
            "nAtoms"
            , nAtoms_function_type( &::SireMol::EditMol::nAtoms )
            , ( bp::arg("cgnum") ) );
    
    }

    { //::SireMol::EditMol::nBonds
    
        typedef int ( ::SireMol::EditMol::*nBonds_function_type )(  ) const;
        
        EditMol_exposer.def( 
            "nBonds"
            , nBonds_function_type( &::SireMol::EditMol::nBonds ) );
    
    }

    { //::SireMol::EditMol::nBonds
    
        typedef int ( ::SireMol::EditMol::*nBonds_function_type )( ::SireMol::ResNum ) const;
        
        EditMol_exposer.def( 
            "nBonds"
            , nBonds_function_type( &::SireMol::EditMol::nBonds )
            , ( bp::arg("resnum") ) );
    
    }

    { //::SireMol::EditMol::nBonds
    
        typedef int ( ::SireMol::EditMol::*nBonds_function_type )( ::SireMol::ResID ) const;
        
        EditMol_exposer.def( 
            "nBonds"
            , nBonds_function_type( &::SireMol::EditMol::nBonds )
            , ( bp::arg("resid") ) );
    
    }

    { //::SireMol::EditMol::nCutGroups
    
        typedef int ( ::SireMol::EditMol::*nCutGroups_function_type )(  ) const;
        
        EditMol_exposer.def( 
            "nCutGroups"
            , nCutGroups_function_type( &::SireMol::EditMol::nCutGroups ) );
    
    }

    { //::SireMol::EditMol::nInterBonds
    
        typedef int ( ::SireMol::EditMol::*nInterBonds_function_type )(  ) const;
        
        EditMol_exposer.def( 
            "nInterBonds"
            , nInterBonds_function_type( &::SireMol::EditMol::nInterBonds ) );
    
    }

    { //::SireMol::EditMol::nInterBonds
    
        typedef int ( ::SireMol::EditMol::*nInterBonds_function_type )( ::SireMol::ResNum ) const;
        
        EditMol_exposer.def( 
            "nInterBonds"
            , nInterBonds_function_type( &::SireMol::EditMol::nInterBonds )
            , ( bp::arg("resnum") ) );
    
    }

    { //::SireMol::EditMol::nInterBonds
    
        typedef int ( ::SireMol::EditMol::*nInterBonds_function_type )( ::SireMol::ResID ) const;
        
        EditMol_exposer.def( 
            "nInterBonds"
            , nInterBonds_function_type( &::SireMol::EditMol::nInterBonds )
            , ( bp::arg("resid") ) );
    
    }

    { //::SireMol::EditMol::nIntraBonds
    
        typedef int ( ::SireMol::EditMol::*nIntraBonds_function_type )(  ) const;
        
        EditMol_exposer.def( 
            "nIntraBonds"
            , nIntraBonds_function_type( &::SireMol::EditMol::nIntraBonds ) );
    
    }

    { //::SireMol::EditMol::nIntraBonds
    
        typedef int ( ::SireMol::EditMol::*nIntraBonds_function_type )( ::SireMol::ResNum ) const;
        
        EditMol_exposer.def( 
            "nIntraBonds"
            , nIntraBonds_function_type( &::SireMol::EditMol::nIntraBonds )
            , ( bp::arg("resnum") ) );
    
    }

    { //::SireMol::EditMol::nIntraBonds
    
        typedef int ( ::SireMol::EditMol::*nIntraBonds_function_type )( ::SireMol::ResID ) const;
        
        EditMol_exposer.def( 
            "nIntraBonds"
            , nIntraBonds_function_type( &::SireMol::EditMol::nIntraBonds )
            , ( bp::arg("resid") ) );
    
    }

    { //::SireMol::EditMol::nResidues
    
        typedef int ( ::SireMol::EditMol::*nResidues_function_type )(  ) const;
        
        EditMol_exposer.def( 
            "nResidues"
            , nResidues_function_type( &::SireMol::EditMol::nResidues ) );
    
    }

}
