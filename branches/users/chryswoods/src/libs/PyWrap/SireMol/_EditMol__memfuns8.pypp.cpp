// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "_EditMol__memfuns8.pypp.hpp"
#include "boost/python.hpp"
#include "siremol_headers.h"
#include "SireMaths/angle.h"
#include "SireMaths/quaternion.h"
#include "SireMaths/matrix.h"
#include "SireMaths/triangle.h"
#include "SireMaths/line.h"
#include "SireMaths/torsion.h"

namespace bp = boost::python;

void register_EditMol_memfuns8( EditMol_exposer_t& EditMol_exposer ){

    { //::SireMol::EditMol::name
    
        typedef ::QString ( ::SireMol::EditMol::*name_function_type )(  ) const;
        
        EditMol_exposer.def( 
            "name"
            , name_function_type( &::SireMol::EditMol::name ) );
    
    }

    { //::SireMol::EditMol::operator[]
    
        typedef ::SireMol::CutGroup ( ::SireMol::EditMol::*__getitem___function_type )( ::SireMol::CutGroupID ) const;
        
        EditMol_exposer.def( 
            "__getitem__"
            , __getitem___function_type( &::SireMol::EditMol::operator[] )
            , ( bp::arg("cgid") ) );
    
    }

    { //::SireMol::EditMol::operator[]
    
        typedef ::SireMol::CutGroup ( ::SireMol::EditMol::*__getitem___function_type )( ::SireMol::CutGroupNum ) const;
        
        EditMol_exposer.def( 
            "__getitem__"
            , __getitem___function_type( &::SireMol::EditMol::operator[] )
            , ( bp::arg("cgnum") ) );
    
    }

    { //::SireMol::EditMol::operator[]
    
        typedef ::SireMol::EditRes ( ::SireMol::EditMol::*__getitem___function_type )( ::SireMol::ResID ) const;
        
        EditMol_exposer.def( 
            "__getitem__"
            , __getitem___function_type( &::SireMol::EditMol::operator[] )
            , ( bp::arg("resid") ) );
    
    }

    { //::SireMol::EditMol::operator[]
    
        typedef ::SireMol::EditRes ( ::SireMol::EditMol::*__getitem___function_type )( ::SireMol::ResNum ) const;
        
        EditMol_exposer.def( 
            "__getitem__"
            , __getitem___function_type( &::SireMol::EditMol::operator[] )
            , ( bp::arg("resnum") ) );
    
    }

    { //::SireMol::EditMol::operator[]
    
        typedef ::SireMol::Atom ( ::SireMol::EditMol::*__getitem___function_type )( ::SireMol::AtomID ) const;
        
        EditMol_exposer.def( 
            "__getitem__"
            , __getitem___function_type( &::SireMol::EditMol::operator[] )
            , ( bp::arg("atomid") ) );
    
    }

    { //::SireMol::EditMol::operator[]
    
        typedef ::SireMol::Atom ( ::SireMol::EditMol::*__getitem___function_type )( ::SireMol::CGAtomID const & ) const;
        
        EditMol_exposer.def( 
            "__getitem__"
            , __getitem___function_type( &::SireMol::EditMol::operator[] )
            , ( bp::arg("cgatomid") ) );
    
    }

    { //::SireMol::EditMol::operator[]
    
        typedef ::SireMol::Atom ( ::SireMol::EditMol::*__getitem___function_type )( ::SireMol::CGNumAtomID const & ) const;
        
        EditMol_exposer.def( 
            "__getitem__"
            , __getitem___function_type( &::SireMol::EditMol::operator[] )
            , ( bp::arg("cgatomid") ) );
    
    }

    { //::SireMol::EditMol::operator[]
    
        typedef ::SireMol::Atom ( ::SireMol::EditMol::*__getitem___function_type )( ::SireMol::ResNumAtomID const & ) const;
        
        EditMol_exposer.def( 
            "__getitem__"
            , __getitem___function_type( &::SireMol::EditMol::operator[] )
            , ( bp::arg("resatomid") ) );
    
    }

    { //::SireMol::EditMol::operator[]
    
        typedef ::SireMol::Atom ( ::SireMol::EditMol::*__getitem___function_type )( ::SireMol::ResIDAtomID const & ) const;
        
        EditMol_exposer.def( 
            "__getitem__"
            , __getitem___function_type( &::SireMol::EditMol::operator[] )
            , ( bp::arg("resatomid") ) );
    
    }

    { //::SireMol::EditMol::operator[]
    
        typedef ::SireMol::Atom ( ::SireMol::EditMol::*__getitem___function_type )( ::SireMol::AtomIndex const & ) const;
        
        EditMol_exposer.def( 
            "__getitem__"
            , __getitem___function_type( &::SireMol::EditMol::operator[] )
            , ( bp::arg("atm") ) );
    
    }

    { //::SireMol::EditMol::remove
    
        typedef void ( ::SireMol::EditMol::*remove_function_type )( ::SireMol::AtomIndex const & ) ;
        
        EditMol_exposer.def( 
            "remove"
            , remove_function_type( &::SireMol::EditMol::remove )
            , ( bp::arg("atom") ) );
    
    }

    { //::SireMol::EditMol::remove
    
        typedef void ( ::SireMol::EditMol::*remove_function_type )( ::SireMol::Bond const & ) ;
        
        EditMol_exposer.def( 
            "remove"
            , remove_function_type( &::SireMol::EditMol::remove )
            , ( bp::arg("bond") ) );
    
    }

    { //::SireMol::EditMol::remove
    
        typedef void ( ::SireMol::EditMol::*remove_function_type )( ::SireMol::ResNum ) ;
        
        EditMol_exposer.def( 
            "remove"
            , remove_function_type( &::SireMol::EditMol::remove )
            , ( bp::arg("resnum") ) );
    
    }

    { //::SireMol::EditMol::remove
    
        typedef void ( ::SireMol::EditMol::*remove_function_type )( ::SireMol::ResID ) ;
        
        EditMol_exposer.def( 
            "remove"
            , remove_function_type( &::SireMol::EditMol::remove )
            , ( bp::arg("resid") ) );
    
    }

    { //::SireMol::EditMol::removeAllBonds
    
        typedef void ( ::SireMol::EditMol::*removeAllBonds_function_type )( ::SireMol::ResNum ) ;
        
        EditMol_exposer.def( 
            "removeAllBonds"
            , removeAllBonds_function_type( &::SireMol::EditMol::removeAllBonds )
            , ( bp::arg("resnum") ) );
    
    }

    { //::SireMol::EditMol::removeAllBonds
    
        typedef void ( ::SireMol::EditMol::*removeAllBonds_function_type )( ::SireMol::ResID ) ;
        
        EditMol_exposer.def( 
            "removeAllBonds"
            , removeAllBonds_function_type( &::SireMol::EditMol::removeAllBonds )
            , ( bp::arg("resnum") ) );
    
    }

    { //::SireMol::EditMol::removeAllBonds
    
        typedef void ( ::SireMol::EditMol::*removeAllBonds_function_type )( ::SireMol::AtomIndex const & ) ;
        
        EditMol_exposer.def( 
            "removeAllBonds"
            , removeAllBonds_function_type( &::SireMol::EditMol::removeAllBonds )
            , ( bp::arg("atom") ) );
    
    }

    { //::SireMol::EditMol::removeAllBonds
    
        typedef void ( ::SireMol::EditMol::*removeAllBonds_function_type )(  ) ;
        
        EditMol_exposer.def( 
            "removeAllBonds"
            , removeAllBonds_function_type( &::SireMol::EditMol::removeAllBonds ) );
    
    }

    { //::SireMol::EditMol::removeAllInterBonds
    
        typedef void ( ::SireMol::EditMol::*removeAllInterBonds_function_type )( ::SireMol::ResNum ) ;
        
        EditMol_exposer.def( 
            "removeAllInterBonds"
            , removeAllInterBonds_function_type( &::SireMol::EditMol::removeAllInterBonds )
            , ( bp::arg("resnum") ) );
    
    }

}
