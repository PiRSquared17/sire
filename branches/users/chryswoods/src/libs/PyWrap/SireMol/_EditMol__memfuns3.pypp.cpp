// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "_EditMol__memfuns3.pypp.hpp"
#include "boost/python.hpp"
#include "siremol_headers.h"
#include "SireMaths/angle.h"
#include "SireMaths/quaternion.h"
#include "SireMaths/matrix.h"
#include "SireMaths/triangle.h"
#include "SireMaths/line.h"
#include "SireMaths/torsion.h"

namespace bp = boost::python;

void register_EditMol_memfuns3( EditMol_exposer_t& EditMol_exposer ){

    { //::SireMol::EditMol::commit
    
        typedef ::SireMol::detail::MolData ( ::SireMol::EditMol::*commit_function_type )(  ) const;
        
        EditMol_exposer.def( 
            "commit"
            , commit_function_type( &::SireMol::EditMol::commit ) );
    
    }

    { //::SireMol::EditMol::connectivity
    
        typedef ::SireMol::MoleculeBonds ( ::SireMol::EditMol::*connectivity_function_type )(  ) const;
        
        EditMol_exposer.def( 
            "connectivity"
            , connectivity_function_type( &::SireMol::EditMol::connectivity ) );
    
    }

    { //::SireMol::EditMol::connectivity
    
        typedef ::SireMol::ResidueBonds ( ::SireMol::EditMol::*connectivity_function_type )( ::SireMol::ResNum ) const;
        
        EditMol_exposer.def( 
            "connectivity"
            , connectivity_function_type( &::SireMol::EditMol::connectivity )
            , ( bp::arg("resnum") ) );
    
    }

    { //::SireMol::EditMol::connectivity
    
        typedef ::SireMol::ResidueBonds ( ::SireMol::EditMol::*connectivity_function_type )( ::SireMol::ResID ) const;
        
        EditMol_exposer.def( 
            "connectivity"
            , connectivity_function_type( &::SireMol::EditMol::connectivity )
            , ( bp::arg("resid") ) );
    
    }

    { //::SireMol::EditMol::contains
    
        typedef bool ( ::SireMol::EditMol::*contains_function_type )( ::SireMol::CutGroupID ) const;
        
        EditMol_exposer.def( 
            "contains"
            , contains_function_type( &::SireMol::EditMol::contains )
            , ( bp::arg("cgid") ) );
    
    }

    { //::SireMol::EditMol::contains
    
        typedef bool ( ::SireMol::EditMol::*contains_function_type )( ::SireMol::CutGroupNum ) const;
        
        EditMol_exposer.def( 
            "contains"
            , contains_function_type( &::SireMol::EditMol::contains )
            , ( bp::arg("cgnum") ) );
    
    }

    { //::SireMol::EditMol::contains
    
        typedef bool ( ::SireMol::EditMol::*contains_function_type )( ::SireMol::ResNum ) const;
        
        EditMol_exposer.def( 
            "contains"
            , contains_function_type( &::SireMol::EditMol::contains )
            , ( bp::arg("resnum") ) );
    
    }

    { //::SireMol::EditMol::contains
    
        typedef bool ( ::SireMol::EditMol::*contains_function_type )( ::SireMol::ResID ) const;
        
        EditMol_exposer.def( 
            "contains"
            , contains_function_type( &::SireMol::EditMol::contains )
            , ( bp::arg("resid") ) );
    
    }

    { //::SireMol::EditMol::contains
    
        typedef bool ( ::SireMol::EditMol::*contains_function_type )( ::SireMol::ResNum,::QString const & ) const;
        
        EditMol_exposer.def( 
            "contains"
            , contains_function_type( &::SireMol::EditMol::contains )
            , ( bp::arg("resnum"), bp::arg("atomname") ) );
    
    }

    { //::SireMol::EditMol::contains
    
        typedef bool ( ::SireMol::EditMol::*contains_function_type )( ::SireMol::AtomIndex const & ) const;
        
        EditMol_exposer.def( 
            "contains"
            , contains_function_type( &::SireMol::EditMol::contains )
            , ( bp::arg("atm") ) );
    
    }

    { //::SireMol::EditMol::contains
    
        typedef bool ( ::SireMol::EditMol::*contains_function_type )( ::SireMol::CutGroupID,::SireMol::AtomID ) const;
        
        EditMol_exposer.def( 
            "contains"
            , contains_function_type( &::SireMol::EditMol::contains )
            , ( bp::arg("cgid"), bp::arg("atomid") ) );
    
    }

    { //::SireMol::EditMol::contains
    
        typedef bool ( ::SireMol::EditMol::*contains_function_type )( ::SireMol::CutGroupNum,::SireMol::AtomID ) const;
        
        EditMol_exposer.def( 
            "contains"
            , contains_function_type( &::SireMol::EditMol::contains )
            , ( bp::arg("cgnum"), bp::arg("atomid") ) );
    
    }

    { //::SireMol::EditMol::contains
    
        typedef bool ( ::SireMol::EditMol::*contains_function_type )( ::SireMol::CGAtomID const & ) const;
        
        EditMol_exposer.def( 
            "contains"
            , contains_function_type( &::SireMol::EditMol::contains )
            , ( bp::arg("cgatomid") ) );
    
    }

    { //::SireMol::EditMol::contains
    
        typedef bool ( ::SireMol::EditMol::*contains_function_type )( ::SireMol::CGNumAtomID const & ) const;
        
        EditMol_exposer.def( 
            "contains"
            , contains_function_type( &::SireMol::EditMol::contains )
            , ( bp::arg("cgatomid") ) );
    
    }

    { //::SireMol::EditMol::contains
    
        typedef bool ( ::SireMol::EditMol::*contains_function_type )( ::SireMol::ResNum,::SireMol::AtomID ) const;
        
        EditMol_exposer.def( 
            "contains"
            , contains_function_type( &::SireMol::EditMol::contains )
            , ( bp::arg("resnum"), bp::arg("atomid") ) );
    
    }

    { //::SireMol::EditMol::contains
    
        typedef bool ( ::SireMol::EditMol::*contains_function_type )( ::SireMol::ResNumAtomID const & ) const;
        
        EditMol_exposer.def( 
            "contains"
            , contains_function_type( &::SireMol::EditMol::contains )
            , ( bp::arg("resatomid") ) );
    
    }

    { //::SireMol::EditMol::contains
    
        typedef bool ( ::SireMol::EditMol::*contains_function_type )( ::SireMol::ResID,::SireMol::AtomID ) const;
        
        EditMol_exposer.def( 
            "contains"
            , contains_function_type( &::SireMol::EditMol::contains )
            , ( bp::arg("resid"), bp::arg("atomid") ) );
    
    }

    { //::SireMol::EditMol::contains
    
        typedef bool ( ::SireMol::EditMol::*contains_function_type )( ::SireMol::ResIDAtomID const & ) const;
        
        EditMol_exposer.def( 
            "contains"
            , contains_function_type( &::SireMol::EditMol::contains )
            , ( bp::arg("resatomid") ) );
    
    }

    { //::SireMol::EditMol::contains
    
        typedef bool ( ::SireMol::EditMol::*contains_function_type )( ::SireMol::Bond const & ) const;
        
        EditMol_exposer.def( 
            "contains"
            , contains_function_type( &::SireMol::EditMol::contains )
            , ( bp::arg("bond") ) );
    
    }

    { //::SireMol::EditMol::coordGroup
    
        typedef ::SireVol::CoordGroup ( ::SireMol::EditMol::*coordGroup_function_type )( ::SireMol::CutGroupID ) const;
        
        EditMol_exposer.def( 
            "coordGroup"
            , coordGroup_function_type( &::SireMol::EditMol::coordGroup )
            , ( bp::arg("id") ) );
    
    }

}
