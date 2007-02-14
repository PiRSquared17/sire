// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "_Residue__memfuns3.pypp.hpp"
#include "boost/python.hpp"
#include "siremol_headers.h"
#include "SireMaths/angle.h"
#include "SireMaths/quaternion.h"
#include "SireMaths/matrix.h"
#include "SireMaths/triangle.h"
#include "SireMaths/line.h"
#include "SireMaths/torsion.h"

namespace bp = boost::python;

void register_Residue_memfuns3( Residue_exposer_t& Residue_exposer ){

    { //::SireMol::Residue::nInterBonds
    
        typedef int ( ::SireMol::Residue::*nInterBonds_function_type )(  ) const;
        
        Residue_exposer.def( 
            "nInterBonds"
            , nInterBonds_function_type( &::SireMol::Residue::nInterBonds ) );
    
    }

    { //::SireMol::Residue::nIntraBonds
    
        typedef int ( ::SireMol::Residue::*nIntraBonds_function_type )(  ) const;
        
        Residue_exposer.def( 
            "nIntraBonds"
            , nIntraBonds_function_type( &::SireMol::Residue::nIntraBonds ) );
    
    }

    { //::SireMol::Residue::name
    
        typedef ::QString ( ::SireMol::Residue::*name_function_type )(  ) const;
        
        Residue_exposer.def( 
            "name"
            , name_function_type( &::SireMol::Residue::name ) );
    
    }

    { //::SireMol::Residue::number
    
        typedef ::SireMol::ResNum ( ::SireMol::Residue::*number_function_type )(  ) const;
        
        Residue_exposer.def( 
            "number"
            , number_function_type( &::SireMol::Residue::number ) );
    
    }

    { //::SireMol::Residue::operator[]
    
        typedef ::SireMol::Atom ( ::SireMol::Residue::*__getitem___function_type )( ::SireMol::AtomID ) const;
        
        Residue_exposer.def( 
            "__getitem__"
            , __getitem___function_type( &::SireMol::Residue::operator[] )
            , ( bp::arg("i") ) );
    
    }

    { //::SireMol::Residue::operator[]
    
        typedef ::SireMol::Atom ( ::SireMol::Residue::*__getitem___function_type )( ::QString const & ) const;
        
        Residue_exposer.def( 
            "__getitem__"
            , __getitem___function_type( &::SireMol::Residue::operator[] )
            , ( bp::arg("atomname") ) );
    
    }

    { //::SireMol::Residue::operator[]
    
        typedef ::SireMol::Atom ( ::SireMol::Residue::*__getitem___function_type )( ::SireMol::AtomIndex const & ) const;
        
        Residue_exposer.def( 
            "__getitem__"
            , __getitem___function_type( &::SireMol::Residue::operator[] )
            , ( bp::arg("atom") ) );
    
    }

    { //::SireMol::Residue::operator[]
    
        typedef ::SireMol::CutGroup ( ::SireMol::Residue::*__getitem___function_type )( ::SireMol::CutGroupID ) const;
        
        Residue_exposer.def( 
            "__getitem__"
            , __getitem___function_type( &::SireMol::Residue::operator[] )
            , ( bp::arg("cgid") ) );
    
    }

    { //::SireMol::Residue::operator[]
    
        typedef ::SireMol::CutGroup ( ::SireMol::Residue::*__getitem___function_type )( ::SireMol::CutGroupNum ) const;
        
        Residue_exposer.def( 
            "__getitem__"
            , __getitem___function_type( &::SireMol::Residue::operator[] )
            , ( bp::arg("cgnum") ) );
    
    }

    { //::SireMol::Residue::resName
    
        typedef ::QString ( ::SireMol::Residue::*resName_function_type )(  ) const;
        
        Residue_exposer.def( 
            "resName"
            , resName_function_type( &::SireMol::Residue::resName ) );
    
    }

    { //::SireMol::Residue::resNum
    
        typedef ::SireMol::ResNum ( ::SireMol::Residue::*resNum_function_type )(  ) const;
        
        Residue_exposer.def( 
            "resNum"
            , resNum_function_type( &::SireMol::Residue::resNum ) );
    
    }

    { //::SireMol::Residue::residuesBondedTo
    
        typedef ::QHash<SireMol::ResNum,SireMol::Residue> ( ::SireMol::Residue::*residuesBondedTo_function_type )( ::SireMol::AtomID ) const;
        
        Residue_exposer.def( 
            "residuesBondedTo"
            , residuesBondedTo_function_type( &::SireMol::Residue::residuesBondedTo )
            , ( bp::arg("atom") ) );
    
    }

    { //::SireMol::Residue::residuesBondedTo
    
        typedef ::QHash<SireMol::ResNum,SireMol::Residue> ( ::SireMol::Residue::*residuesBondedTo_function_type )( ::QString const & ) const;
        
        Residue_exposer.def( 
            "residuesBondedTo"
            , residuesBondedTo_function_type( &::SireMol::Residue::residuesBondedTo )
            , ( bp::arg("atomname") ) );
    
    }

    { //::SireMol::Residue::residuesBondedTo
    
        typedef ::QHash<SireMol::ResNum,SireMol::Residue> ( ::SireMol::Residue::*residuesBondedTo_function_type )( ::SireMol::AtomIndex const & ) const;
        
        Residue_exposer.def( 
            "residuesBondedTo"
            , residuesBondedTo_function_type( &::SireMol::Residue::residuesBondedTo )
            , ( bp::arg("atom") ) );
    
    }

    { //::SireMol::Residue::rotate
    
        typedef void ( ::SireMol::Residue::*rotate_function_type )( ::SireMaths::Quaternion const &,::SireMaths::Vector const & ) ;
        
        Residue_exposer.def( 
            "rotate"
            , rotate_function_type( &::SireMol::Residue::rotate )
            , ( bp::arg("quat"), bp::arg("point") ) );
    
    }

    { //::SireMol::Residue::rotate
    
        typedef void ( ::SireMol::Residue::*rotate_function_type )( ::SireMol::AtomID,::SireMaths::Quaternion const &,::SireMaths::Vector const & ) ;
        
        Residue_exposer.def( 
            "rotate"
            , rotate_function_type( &::SireMol::Residue::rotate )
            , ( bp::arg("atomid"), bp::arg("quat"), bp::arg("point") ) );
    
    }

    { //::SireMol::Residue::rotate
    
        typedef void ( ::SireMol::Residue::*rotate_function_type )( ::QSet<SireMol::AtomID> const &,::SireMaths::Quaternion const &,::SireMaths::Vector const & ) ;
        
        Residue_exposer.def( 
            "rotate"
            , rotate_function_type( &::SireMol::Residue::rotate )
            , ( bp::arg("atomids"), bp::arg("quat"), bp::arg("point") ) );
    
    }

    { //::SireMol::Residue::rotate
    
        typedef void ( ::SireMol::Residue::*rotate_function_type )( ::QString const &,::SireMaths::Quaternion const &,::SireMaths::Vector const & ) ;
        
        Residue_exposer.def( 
            "rotate"
            , rotate_function_type( &::SireMol::Residue::rotate )
            , ( bp::arg("atom"), bp::arg("quat"), bp::arg("point") ) );
    
    }

    { //::SireMol::Residue::rotate
    
        typedef void ( ::SireMol::Residue::*rotate_function_type )( ::QStringList const &,::SireMaths::Quaternion const &,::SireMaths::Vector const & ) ;
        
        Residue_exposer.def( 
            "rotate"
            , rotate_function_type( &::SireMol::Residue::rotate )
            , ( bp::arg("atoms"), bp::arg("quat"), bp::arg("point") ) );
    
    }

    { //::SireMol::Residue::rotate
    
        typedef void ( ::SireMol::Residue::*rotate_function_type )( ::SireMol::AtomIndex const &,::SireMaths::Quaternion const &,::SireMaths::Vector const & ) ;
        
        Residue_exposer.def( 
            "rotate"
            , rotate_function_type( &::SireMol::Residue::rotate )
            , ( bp::arg("atom"), bp::arg("quat"), bp::arg("point") ) );
    
    }

}
