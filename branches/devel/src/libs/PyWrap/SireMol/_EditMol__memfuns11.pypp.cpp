// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "_EditMol__memfuns11.pypp.hpp"
#include "boost/python.hpp"
#include "siremol_headers.h"
#include "SireMaths/angle.h"
#include "SireMaths/quaternion.h"
#include "SireMaths/matrix.h"
#include "SireMaths/triangle.h"
#include "SireMaths/line.h"
#include "SireMaths/torsion.h"

namespace bp = boost::python;

void register_EditMol_memfuns11( EditMol_exposer_t& EditMol_exposer ){

    { //::SireMol::EditMol::rotate
    
        typedef void ( ::SireMol::EditMol::*rotate_function_type )( ::SireMol::AtomIDGroup const &,::SireMaths::Matrix const &,::SireMaths::Vector const & ) ;
        
        EditMol_exposer.def( 
            "rotate"
            , rotate_function_type( &::SireMol::EditMol::rotate )
            , ( bp::arg("group"), bp::arg("matrix"), bp::arg("point") ) );
    
    }

    { //::SireMol::EditMol::rotate
    
        typedef void ( ::SireMol::EditMol::*rotate_function_type )( ::SireMol::AtomIndex const &,::SireMaths::Matrix const &,::SireMaths::Vector const & ) ;
        
        EditMol_exposer.def( 
            "rotate"
            , rotate_function_type( &::SireMol::EditMol::rotate )
            , ( bp::arg("atom"), bp::arg("matrix"), bp::arg("point") ) );
    
    }

    { //::SireMol::EditMol::rotate
    
        typedef void ( ::SireMol::EditMol::*rotate_function_type )( ::QSet<SireMol::AtomIndex> const &,::SireMaths::Matrix const &,::SireMaths::Vector const & ) ;
        
        EditMol_exposer.def( 
            "rotate"
            , rotate_function_type( &::SireMol::EditMol::rotate )
            , ( bp::arg("atoms"), bp::arg("matrix"), bp::arg("point") ) );
    
    }

    { //::SireMol::EditMol::rotate
    
        typedef void ( ::SireMol::EditMol::*rotate_function_type )( ::SireMol::ResNum,::QStringList const &,::SireMaths::Matrix const &,::SireMaths::Vector const & ) ;
        
        EditMol_exposer.def( 
            "rotate"
            , rotate_function_type( &::SireMol::EditMol::rotate )
            , ( bp::arg("resnum"), bp::arg("atoms"), bp::arg("matrix"), bp::arg("point") ) );
    
    }

    { //::SireMol::EditMol::rotate
    
        typedef void ( ::SireMol::EditMol::*rotate_function_type )( ::SireMol::ResNum,::SireMaths::Matrix const &,::SireMaths::Vector const & ) ;
        
        EditMol_exposer.def( 
            "rotate"
            , rotate_function_type( &::SireMol::EditMol::rotate )
            , ( bp::arg("resnum"), bp::arg("matrix"), bp::arg("point") ) );
    
    }

    { //::SireMol::EditMol::rotate
    
        typedef void ( ::SireMol::EditMol::*rotate_function_type )( ::QSet<SireMol::ResNum> const &,::SireMaths::Matrix const &,::SireMaths::Vector const & ) ;
        
        EditMol_exposer.def( 
            "rotate"
            , rotate_function_type( &::SireMol::EditMol::rotate )
            , ( bp::arg("resnums"), bp::arg("matrix"), bp::arg("point") ) );
    
    }

    { //::SireMol::EditMol::rotate
    
        typedef void ( ::SireMol::EditMol::*rotate_function_type )( ::SireMol::ResID,::QStringList const &,::SireMaths::Matrix const &,::SireMaths::Vector const & ) ;
        
        EditMol_exposer.def( 
            "rotate"
            , rotate_function_type( &::SireMol::EditMol::rotate )
            , ( bp::arg("resid"), bp::arg("atoms"), bp::arg("matrix"), bp::arg("point") ) );
    
    }

    { //::SireMol::EditMol::rotate
    
        typedef void ( ::SireMol::EditMol::*rotate_function_type )( ::SireMol::ResID,::SireMaths::Matrix const &,::SireMaths::Vector const & ) ;
        
        EditMol_exposer.def( 
            "rotate"
            , rotate_function_type( &::SireMol::EditMol::rotate )
            , ( bp::arg("resid"), bp::arg("matrix"), bp::arg("point") ) );
    
    }

    { //::SireMol::EditMol::rotate
    
        typedef void ( ::SireMol::EditMol::*rotate_function_type )( ::QSet<SireMol::ResID> const &,::SireMaths::Matrix const &,::SireMaths::Vector const & ) ;
        
        EditMol_exposer.def( 
            "rotate"
            , rotate_function_type( &::SireMol::EditMol::rotate )
            , ( bp::arg("resids"), bp::arg("matrix"), bp::arg("point") ) );
    
    }

    { //::SireMol::EditMol::rotate
    
        typedef void ( ::SireMol::EditMol::*rotate_function_type )( ::SireMol::CutGroupNum,::SireMaths::Matrix const &,::SireMaths::Vector const & ) ;
        
        EditMol_exposer.def( 
            "rotate"
            , rotate_function_type( &::SireMol::EditMol::rotate )
            , ( bp::arg("cgnum"), bp::arg("matrix"), bp::arg("point") ) );
    
    }

    { //::SireMol::EditMol::rotate
    
        typedef void ( ::SireMol::EditMol::*rotate_function_type )( ::QSet<SireMol::CutGroupNum> const &,::SireMaths::Matrix const &,::SireMaths::Vector const & ) ;
        
        EditMol_exposer.def( 
            "rotate"
            , rotate_function_type( &::SireMol::EditMol::rotate )
            , ( bp::arg("cgnums"), bp::arg("matrix"), bp::arg("point") ) );
    
    }

    { //::SireMol::EditMol::rotate
    
        typedef void ( ::SireMol::EditMol::*rotate_function_type )( ::SireMol::CutGroupID,::SireMaths::Matrix const &,::SireMaths::Vector const & ) ;
        
        EditMol_exposer.def( 
            "rotate"
            , rotate_function_type( &::SireMol::EditMol::rotate )
            , ( bp::arg("cgid"), bp::arg("matrix"), bp::arg("point") ) );
    
    }

    { //::SireMol::EditMol::rotate
    
        typedef void ( ::SireMol::EditMol::*rotate_function_type )( ::QSet<SireMol::CutGroupID> const &,::SireMaths::Matrix const &,::SireMaths::Vector const & ) ;
        
        EditMol_exposer.def( 
            "rotate"
            , rotate_function_type( &::SireMol::EditMol::rotate )
            , ( bp::arg("cgids"), bp::arg("matrix"), bp::arg("point") ) );
    
    }

    { //::SireMol::EditMol::set
    
        typedef void ( ::SireMol::EditMol::*set_function_type )( ::SireMol::Bond const &,double,::QSet<SireMol::AtomIndex> const & ) ;
        
        EditMol_exposer.def( 
            "set"
            , set_function_type( &::SireMol::EditMol::set )
            , ( bp::arg("bond"), bp::arg("lgth"), bp::arg("anchors")=::QSet<SireMol::AtomIndex>( ) ) );
    
    }

    { //::SireMol::EditMol::set
    
        typedef void ( ::SireMol::EditMol::*set_function_type )( ::SireMol::Bond const &,double,::SireMol::WeightFunction const &,::QSet<SireMol::AtomIndex> const & ) ;
        
        EditMol_exposer.def( 
            "set"
            , set_function_type( &::SireMol::EditMol::set )
            , ( bp::arg("bond"), bp::arg("lgth"), bp::arg("func"), bp::arg("anchors")=::QSet<SireMol::AtomIndex>( ) ) );
    
    }

    { //::SireMol::EditMol::set
    
        typedef void ( ::SireMol::EditMol::*set_function_type )( ::SireMol::Angle const &,::SireMaths::Angle const &,::QSet<SireMol::AtomIndex> const & ) ;
        
        EditMol_exposer.def( 
            "set"
            , set_function_type( &::SireMol::EditMol::set )
            , ( bp::arg("angle"), bp::arg("ang"), bp::arg("anchors")=::QSet<SireMol::AtomIndex>( ) ) );
    
    }

    { //::SireMol::EditMol::set
    
        typedef void ( ::SireMol::EditMol::*set_function_type )( ::SireMol::Angle const &,::SireMaths::Angle const &,::SireMol::WeightFunction const &,::QSet<SireMol::AtomIndex> const & ) ;
        
        EditMol_exposer.def( 
            "set"
            , set_function_type( &::SireMol::EditMol::set )
            , ( bp::arg("angle"), bp::arg("ang"), bp::arg("func"), bp::arg("anchors")=::QSet<SireMol::AtomIndex>( ) ) );
    
    }

    { //::SireMol::EditMol::set
    
        typedef void ( ::SireMol::EditMol::*set_function_type )( ::SireMol::Dihedral const &,::SireMaths::Angle const &,::QSet<SireMol::AtomIndex> const & ) ;
        
        EditMol_exposer.def( 
            "set"
            , set_function_type( &::SireMol::EditMol::set )
            , ( bp::arg("dihedral"), bp::arg("ang"), bp::arg("anchors")=::QSet<SireMol::AtomIndex>( ) ) );
    
    }

    { //::SireMol::EditMol::set
    
        typedef void ( ::SireMol::EditMol::*set_function_type )( ::SireMol::Dihedral const &,::SireMaths::Angle const &,::SireMol::WeightFunction const &,::QSet<SireMol::AtomIndex> const & ) ;
        
        EditMol_exposer.def( 
            "set"
            , set_function_type( &::SireMol::EditMol::set )
            , ( bp::arg("dihedral"), bp::arg("ang"), bp::arg("func"), bp::arg("anchors")=::QSet<SireMol::AtomIndex>( ) ) );
    
    }

    { //::SireMol::EditMol::set
    
        typedef void ( ::SireMol::EditMol::*set_function_type )( ::SireMol::Improper const &,::SireMaths::Angle const &,::QSet<SireMol::AtomIndex> const & ) ;
        
        EditMol_exposer.def( 
            "set"
            , set_function_type( &::SireMol::EditMol::set )
            , ( bp::arg("improper"), bp::arg("ang"), bp::arg("anchors")=::QSet<SireMol::AtomIndex>( ) ) );
    
    }

}
