// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "_EditRes__memfuns0.pypp.hpp"
#include "boost/python.hpp"
#include "siremol_headers.h"
#include "SireMol/moleculedata.h"
#include "SireMaths/angle.h"
#include "SireMaths/quaternion.h"
#include "SireMaths/matrix.h"
#include "SireMaths/triangle.h"
#include "SireMaths/line.h"
#include "SireMaths/torsion.h"
#include "SireVol/space.h"

namespace bp = boost::python;

void register_EditRes_memfuns0( EditRes_exposer_t& EditRes_exposer ){

    { //::SireMol::EditRes::add
    
        typedef void ( ::SireMol::EditRes::*add_function_type )( ::QString const & ) ;
        
        EditRes_exposer.def( 
            "add"
            , add_function_type( &::SireMol::EditRes::add )
            , ( bp::arg("atm") ) );
    
    }

    { //::SireMol::EditRes::add
    
        typedef void ( ::SireMol::EditRes::*add_function_type )( ::SireMol::Atom const & ) ;
        
        EditRes_exposer.def( 
            "add"
            , add_function_type( &::SireMol::EditRes::add )
            , ( bp::arg("atm") ) );
    
    }

    { //::SireMol::EditRes::add
    
        typedef void ( ::SireMol::EditRes::*add_function_type )( ::QString const &,::SireMol::CutGroupNum ) ;
        
        EditRes_exposer.def( 
            "add"
            , add_function_type( &::SireMol::EditRes::add )
            , ( bp::arg("atm"), bp::arg("cgnum") ) );
    
    }

    { //::SireMol::EditRes::add
    
        typedef void ( ::SireMol::EditRes::*add_function_type )( ::SireMol::Atom const &,::SireMol::CutGroupNum ) ;
        
        EditRes_exposer.def( 
            "add"
            , add_function_type( &::SireMol::EditRes::add )
            , ( bp::arg("atm"), bp::arg("cgnum") ) );
    
    }

    { //::SireMol::EditRes::add
    
        typedef void ( ::SireMol::EditRes::*add_function_type )( ::SireMol::Bond const & ) ;
        
        EditRes_exposer.def( 
            "add"
            , add_function_type( &::SireMol::EditRes::add )
            , ( bp::arg("bond") ) );
    
    }

    { //::SireMol::EditRes::addAutoBonds
    
        typedef void ( ::SireMol::EditRes::*addAutoBonds_function_type )(  ) ;
        
        EditRes_exposer.def( 
            "addAutoBonds"
            , addAutoBonds_function_type( &::SireMol::EditRes::addAutoBonds ) );
    
    }

    { //::SireMol::EditRes::addAutoBonds
    
        typedef void ( ::SireMol::EditRes::*addAutoBonds_function_type )( ::SireMol::BondAddingFunction const & ) ;
        
        EditRes_exposer.def( 
            "addAutoBonds"
            , addAutoBonds_function_type( &::SireMol::EditRes::addAutoBonds )
            , ( bp::arg("bondfunc") ) );
    
    }

    { //::SireMol::EditRes::addBond
    
        typedef void ( ::SireMol::EditRes::*addBond_function_type )( ::QString const &,::QString const & ) ;
        
        EditRes_exposer.def( 
            "addBond"
            , addBond_function_type( &::SireMol::EditRes::addBond )
            , ( bp::arg("atmnam0"), bp::arg("atmnam1") ) );
    
    }

    { //::SireMol::EditRes::angle
    
        typedef ::SireMaths::Triangle ( ::SireMol::EditRes::*angle_function_type )( ::SireMol::Angle const & ) const;
        
        EditRes_exposer.def( 
            "angle"
            , angle_function_type( &::SireMol::EditRes::angle )
            , ( bp::arg("ang") ) );
    
    }

    { //::SireMol::EditRes::applyTemplate
    
        typedef void ( ::SireMol::EditRes::*applyTemplate_function_type )( ::SireMol::EditRes const &,::SireMol::TemplateFunction const & ) ;
        
        EditRes_exposer.def( 
            "applyTemplate"
            , applyTemplate_function_type( &::SireMol::EditRes::applyTemplate )
            , ( bp::arg("tmpl"), bp::arg("tmplfunc") ) );
    
    }

    { //::SireMol::EditRes::assertSameResidue
    
        typedef void ( ::SireMol::EditRes::*assertSameResidue_function_type )( ::SireMol::AtomIndex const & ) const;
        
        EditRes_exposer.def( 
            "assertSameResidue"
            , assertSameResidue_function_type( &::SireMol::EditRes::assertSameResidue )
            , ( bp::arg("atom") ) );
    
    }

    { //::SireMol::EditRes::assertSameResidue
    
        typedef void ( ::SireMol::EditRes::*assertSameResidue_function_type )( ::QSet<SireMol::AtomIndex> const & ) const;
        
        EditRes_exposer.def( 
            "assertSameResidue"
            , assertSameResidue_function_type( &::SireMol::EditRes::assertSameResidue )
            , ( bp::arg("atoms") ) );
    
    }

    { //::SireMol::EditRes::assertSameResidue
    
        typedef void ( ::SireMol::EditRes::*assertSameResidue_function_type )( ::SireMol::CutGroupNum ) const;
        
        EditRes_exposer.def( 
            "assertSameResidue"
            , assertSameResidue_function_type( &::SireMol::EditRes::assertSameResidue )
            , ( bp::arg("cgnum") ) );
    
    }

    { //::SireMol::EditRes::assertSameResidue
    
        typedef void ( ::SireMol::EditRes::*assertSameResidue_function_type )( ::SireMol::CutGroupID ) const;
        
        EditRes_exposer.def( 
            "assertSameResidue"
            , assertSameResidue_function_type( &::SireMol::EditRes::assertSameResidue )
            , ( bp::arg("cgid") ) );
    
    }

    { //::SireMol::EditRes::assertSameResidue
    
        typedef void ( ::SireMol::EditRes::*assertSameResidue_function_type )( ::SireMol::Bond const & ) const;
        
        EditRes_exposer.def( 
            "assertSameResidue"
            , assertSameResidue_function_type( &::SireMol::EditRes::assertSameResidue )
            , ( bp::arg("bond") ) );
    
    }

    { //::SireMol::EditRes::at
    
        typedef ::SireMol::Atom ( ::SireMol::EditRes::*at_function_type )( ::SireMol::AtomID ) const;
        
        EditRes_exposer.def( 
            "at"
            , at_function_type( &::SireMol::EditRes::at )
            , ( bp::arg("i") ) );
    
    }

    { //::SireMol::EditRes::at
    
        typedef ::SireMol::Atom ( ::SireMol::EditRes::*at_function_type )( ::QString const & ) const;
        
        EditRes_exposer.def( 
            "at"
            , at_function_type( &::SireMol::EditRes::at )
            , ( bp::arg("atomname") ) );
    
    }

    { //::SireMol::EditRes::at
    
        typedef ::SireMol::Atom ( ::SireMol::EditRes::*at_function_type )( ::SireMol::AtomIndex const & ) const;
        
        EditRes_exposer.def( 
            "at"
            , at_function_type( &::SireMol::EditRes::at )
            , ( bp::arg("atom") ) );
    
    }

    { //::SireMol::EditRes::at
    
        typedef ::SireMol::CutGroup ( ::SireMol::EditRes::*at_function_type )( ::SireMol::CutGroupID ) const;
        
        EditRes_exposer.def( 
            "at"
            , at_function_type( &::SireMol::EditRes::at )
            , ( bp::arg("cgid") ) );
    
    }

    { //::SireMol::EditRes::at
    
        typedef ::SireMol::CutGroup ( ::SireMol::EditRes::*at_function_type )( ::SireMol::CutGroupNum ) const;
        
        EditRes_exposer.def( 
            "at"
            , at_function_type( &::SireMol::EditRes::at )
            , ( bp::arg("cgnum") ) );
    
    }

}
