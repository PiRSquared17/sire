// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "_Residue__memfuns5.pypp.hpp"
#include "boost/python.hpp"
#include "siremol_headers.h"
#include "SireMaths/angle.h"
#include "SireMaths/quaternion.h"
#include "SireMaths/matrix.h"
#include "SireMaths/triangle.h"
#include "SireMaths/line.h"
#include "SireMaths/torsion.h"

namespace bp = boost::python;

void register_Residue_memfuns5( Residue_exposer_t& Residue_exposer ){

    { //::SireMol::Residue::setCoordinates
    
        typedef void ( ::SireMol::Residue::*setCoordinates_function_type )( ::QHash<SireMol::AtomID,SireMaths::Vector> const & ) ;
        
        Residue_exposer.def( 
            "setCoordinates"
            , setCoordinates_function_type( &::SireMol::Residue::setCoordinates )
            , ( bp::arg("newcoords") ) );
    
    }

    { //::SireMol::Residue::setCoordinates
    
        typedef void ( ::SireMol::Residue::*setCoordinates_function_type )( ::QString const &,::SireMaths::Vector const & ) ;
        
        Residue_exposer.def( 
            "setCoordinates"
            , setCoordinates_function_type( &::SireMol::Residue::setCoordinates )
            , ( bp::arg("atomname"), bp::arg("newcoords") ) );
    
    }

    { //::SireMol::Residue::setCoordinates
    
        typedef void ( ::SireMol::Residue::*setCoordinates_function_type )( ::QHash<QString,SireMaths::Vector> const & ) ;
        
        Residue_exposer.def( 
            "setCoordinates"
            , setCoordinates_function_type( &::SireMol::Residue::setCoordinates )
            , ( bp::arg("newcoords") ) );
    
    }

    { //::SireMol::Residue::setCoordinates
    
        typedef void ( ::SireMol::Residue::*setCoordinates_function_type )( ::SireMol::AtomIndex const &,::SireMaths::Vector const & ) ;
        
        Residue_exposer.def( 
            "setCoordinates"
            , setCoordinates_function_type( &::SireMol::Residue::setCoordinates )
            , ( bp::arg("atom"), bp::arg("newcoords") ) );
    
    }

    { //::SireMol::Residue::setCoordinates
    
        typedef void ( ::SireMol::Residue::*setCoordinates_function_type )( ::QHash<SireMol::AtomIndex,SireMaths::Vector> const & ) ;
        
        Residue_exposer.def( 
            "setCoordinates"
            , setCoordinates_function_type( &::SireMol::Residue::setCoordinates )
            , ( bp::arg("newcoords") ) );
    
    }

    { //::SireMol::Residue::translate
    
        typedef void ( ::SireMol::Residue::*translate_function_type )( ::SireMaths::Vector const & ) ;
        
        Residue_exposer.def( 
            "translate"
            , translate_function_type( &::SireMol::Residue::translate )
            , ( bp::arg("delta") ) );
    
    }

    { //::SireMol::Residue::translate
    
        typedef void ( ::SireMol::Residue::*translate_function_type )( ::SireMol::AtomID,::SireMaths::Vector const & ) ;
        
        Residue_exposer.def( 
            "translate"
            , translate_function_type( &::SireMol::Residue::translate )
            , ( bp::arg("atomid"), bp::arg("delta") ) );
    
    }

    { //::SireMol::Residue::translate
    
        typedef void ( ::SireMol::Residue::*translate_function_type )( ::QSet<SireMol::AtomID> const &,::SireMaths::Vector const & ) ;
        
        Residue_exposer.def( 
            "translate"
            , translate_function_type( &::SireMol::Residue::translate )
            , ( bp::arg("atomids"), bp::arg("delta") ) );
    
    }

    { //::SireMol::Residue::translate
    
        typedef void ( ::SireMol::Residue::*translate_function_type )( ::QString const &,::SireMaths::Vector const & ) ;
        
        Residue_exposer.def( 
            "translate"
            , translate_function_type( &::SireMol::Residue::translate )
            , ( bp::arg("atom"), bp::arg("delta") ) );
    
    }

    { //::SireMol::Residue::translate
    
        typedef void ( ::SireMol::Residue::*translate_function_type )( ::QStringList const &,::SireMaths::Vector const & ) ;
        
        Residue_exposer.def( 
            "translate"
            , translate_function_type( &::SireMol::Residue::translate )
            , ( bp::arg("atoms"), bp::arg("delta") ) );
    
    }

    { //::SireMol::Residue::translate
    
        typedef void ( ::SireMol::Residue::*translate_function_type )( ::SireMol::AtomIndex const &,::SireMaths::Vector const & ) ;
        
        Residue_exposer.def( 
            "translate"
            , translate_function_type( &::SireMol::Residue::translate )
            , ( bp::arg("atom"), bp::arg("delta") ) );
    
    }

    { //::SireMol::Residue::translate
    
        typedef void ( ::SireMol::Residue::*translate_function_type )( ::QSet<SireMol::AtomIndex> const &,::SireMaths::Vector const & ) ;
        
        Residue_exposer.def( 
            "translate"
            , translate_function_type( &::SireMol::Residue::translate )
            , ( bp::arg("atoms"), bp::arg("delta") ) );
    
    }

    { //::SireMol::Residue::version
    
        typedef ::SireMol::MoleculeVersion const & ( ::SireMol::Residue::*version_function_type )(  ) const;
        
        Residue_exposer.def( 
            "version"
            , version_function_type( &::SireMol::Residue::version )
            , bp::return_value_policy< bp::copy_const_reference >() );
    
    }

}
