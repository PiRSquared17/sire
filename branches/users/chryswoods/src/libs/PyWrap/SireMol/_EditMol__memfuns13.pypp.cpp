// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "_EditMol__memfuns13.pypp.hpp"
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

void register_EditMol_memfuns13( EditMol_exposer_t& EditMol_exposer ){

    { //::SireMol::EditMol::setCoordinates
    
        typedef void ( ::SireMol::EditMol::*setCoordinates_function_type )( ::SireMol::CGAtomID const &,::SireMaths::Vector const & ) ;
        
        EditMol_exposer.def( 
            "setCoordinates"
            , setCoordinates_function_type( &::SireMol::EditMol::setCoordinates )
            , ( bp::arg("cgatomid"), bp::arg("newcoords") ) );
    
    }

    { //::SireMol::EditMol::setCoordinates
    
        typedef void ( ::SireMol::EditMol::*setCoordinates_function_type )( ::QHash<SireMol::CGAtomID,SireMaths::Vector> const & ) ;
        
        EditMol_exposer.def( 
            "setCoordinates"
            , setCoordinates_function_type( &::SireMol::EditMol::setCoordinates )
            , ( bp::arg("newcoords") ) );
    
    }

    { //::SireMol::EditMol::setCoordinates
    
        typedef void ( ::SireMol::EditMol::*setCoordinates_function_type )( ::SireMol::ResNumAtomID const &,::SireMaths::Vector const & ) ;
        
        EditMol_exposer.def( 
            "setCoordinates"
            , setCoordinates_function_type( &::SireMol::EditMol::setCoordinates )
            , ( bp::arg("resatomid"), bp::arg("newcoords") ) );
    
    }

    { //::SireMol::EditMol::setCoordinates
    
        typedef void ( ::SireMol::EditMol::*setCoordinates_function_type )( ::QHash<SireMol::ResNumAtomID,SireMaths::Vector> const & ) ;
        
        EditMol_exposer.def( 
            "setCoordinates"
            , setCoordinates_function_type( &::SireMol::EditMol::setCoordinates )
            , ( bp::arg("newcoords") ) );
    
    }

    { //::SireMol::EditMol::setCoordinates
    
        typedef void ( ::SireMol::EditMol::*setCoordinates_function_type )( ::SireMol::ResIDAtomID const &,::SireMaths::Vector const & ) ;
        
        EditMol_exposer.def( 
            "setCoordinates"
            , setCoordinates_function_type( &::SireMol::EditMol::setCoordinates )
            , ( bp::arg("resatomid"), bp::arg("newcoords") ) );
    
    }

    { //::SireMol::EditMol::setCoordinates
    
        typedef void ( ::SireMol::EditMol::*setCoordinates_function_type )( ::QHash<SireMol::ResIDAtomID,SireMaths::Vector> const & ) ;
        
        EditMol_exposer.def( 
            "setCoordinates"
            , setCoordinates_function_type( &::SireMol::EditMol::setCoordinates )
            , ( bp::arg("newcoords") ) );
    
    }

    { //::SireMol::EditMol::setName
    
        typedef void ( ::SireMol::EditMol::*setName_function_type )( ::QString const & ) ;
        
        EditMol_exposer.def( 
            "setName"
            , setName_function_type( &::SireMol::EditMol::setName )
            , ( bp::arg("name") ) );
    
    }

    { //::SireMol::EditMol::setResidueName
    
        typedef void ( ::SireMol::EditMol::*setResidueName_function_type )( ::SireMol::ResNum,::QString const & ) ;
        
        EditMol_exposer.def( 
            "setResidueName"
            , setResidueName_function_type( &::SireMol::EditMol::setResidueName )
            , ( bp::arg("resnum"), bp::arg("name") ) );
    
    }

    { //::SireMol::EditMol::translate
    
        typedef void ( ::SireMol::EditMol::*translate_function_type )( ::SireMaths::Vector const & ) ;
        
        EditMol_exposer.def( 
            "translate"
            , translate_function_type( &::SireMol::EditMol::translate )
            , ( bp::arg("delta") ) );
    
    }

    { //::SireMol::EditMol::translate
    
        typedef void ( ::SireMol::EditMol::*translate_function_type )( ::SireMol::AtomIDGroup const &,::SireMaths::Vector const & ) ;
        
        EditMol_exposer.def( 
            "translate"
            , translate_function_type( &::SireMol::EditMol::translate )
            , ( bp::arg("group"), bp::arg("delta") ) );
    
    }

    { //::SireMol::EditMol::translate
    
        typedef void ( ::SireMol::EditMol::*translate_function_type )( ::SireMol::AtomIndex const &,::SireMaths::Vector const & ) ;
        
        EditMol_exposer.def( 
            "translate"
            , translate_function_type( &::SireMol::EditMol::translate )
            , ( bp::arg("atom"), bp::arg("delta") ) );
    
    }

    { //::SireMol::EditMol::translate
    
        typedef void ( ::SireMol::EditMol::*translate_function_type )( ::QSet<SireMol::AtomIndex> const &,::SireMaths::Vector const & ) ;
        
        EditMol_exposer.def( 
            "translate"
            , translate_function_type( &::SireMol::EditMol::translate )
            , ( bp::arg("atoms"), bp::arg("delta") ) );
    
    }

    { //::SireMol::EditMol::translate
    
        typedef void ( ::SireMol::EditMol::*translate_function_type )( ::SireMol::ResNum,::QStringList const &,::SireMaths::Vector const & ) ;
        
        EditMol_exposer.def( 
            "translate"
            , translate_function_type( &::SireMol::EditMol::translate )
            , ( bp::arg("resnum"), bp::arg("atoms"), bp::arg("delta") ) );
    
    }

    { //::SireMol::EditMol::translate
    
        typedef void ( ::SireMol::EditMol::*translate_function_type )( ::SireMol::ResNum,::SireMaths::Vector const & ) ;
        
        EditMol_exposer.def( 
            "translate"
            , translate_function_type( &::SireMol::EditMol::translate )
            , ( bp::arg("resnum"), bp::arg("delta") ) );
    
    }

    { //::SireMol::EditMol::translate
    
        typedef void ( ::SireMol::EditMol::*translate_function_type )( ::QSet<SireMol::ResNum> const &,::SireMaths::Vector const & ) ;
        
        EditMol_exposer.def( 
            "translate"
            , translate_function_type( &::SireMol::EditMol::translate )
            , ( bp::arg("resnums"), bp::arg("delta") ) );
    
    }

    { //::SireMol::EditMol::translate
    
        typedef void ( ::SireMol::EditMol::*translate_function_type )( ::SireMol::ResID,::QStringList const &,::SireMaths::Vector const & ) ;
        
        EditMol_exposer.def( 
            "translate"
            , translate_function_type( &::SireMol::EditMol::translate )
            , ( bp::arg("resid"), bp::arg("atoms"), bp::arg("delta") ) );
    
    }

    { //::SireMol::EditMol::translate
    
        typedef void ( ::SireMol::EditMol::*translate_function_type )( ::SireMol::ResID,::SireMaths::Vector const & ) ;
        
        EditMol_exposer.def( 
            "translate"
            , translate_function_type( &::SireMol::EditMol::translate )
            , ( bp::arg("resid"), bp::arg("delta") ) );
    
    }

    { //::SireMol::EditMol::translate
    
        typedef void ( ::SireMol::EditMol::*translate_function_type )( ::QSet<SireMol::ResID> const &,::SireMaths::Vector const & ) ;
        
        EditMol_exposer.def( 
            "translate"
            , translate_function_type( &::SireMol::EditMol::translate )
            , ( bp::arg("resids"), bp::arg("delta") ) );
    
    }

    { //::SireMol::EditMol::translate
    
        typedef void ( ::SireMol::EditMol::*translate_function_type )( ::SireMol::CutGroupNum,::SireMaths::Vector const & ) ;
        
        EditMol_exposer.def( 
            "translate"
            , translate_function_type( &::SireMol::EditMol::translate )
            , ( bp::arg("cgnum"), bp::arg("delta") ) );
    
    }

    { //::SireMol::EditMol::translate
    
        typedef void ( ::SireMol::EditMol::*translate_function_type )( ::QSet<SireMol::CutGroupNum> const &,::SireMaths::Vector const & ) ;
        
        EditMol_exposer.def( 
            "translate"
            , translate_function_type( &::SireMol::EditMol::translate )
            , ( bp::arg("cgnums"), bp::arg("delta") ) );
    
    }

}
