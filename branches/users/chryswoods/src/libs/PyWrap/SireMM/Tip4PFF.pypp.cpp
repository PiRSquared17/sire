// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "Tip4PFF.pypp.hpp"
#include "boost/python.hpp"
#include "siremm_headers.h"
#include "SireMol/molecule.h"
#include "SireMol/residue.h"
#include "SireMol/newatom.h"
#include "SireMol/atom.h"

namespace bp = boost::python;

void register_Tip4PFF_class(){

    { //::SireMM::Tip4PFF
        typedef bp::class_< SireMM::Tip4PFF, bp::bases< SireMM::CLJFF > > Tip4PFF_exposer_t;
        Tip4PFF_exposer_t Tip4PFF_exposer = Tip4PFF_exposer_t( "Tip4PFF" );
        bp::scope Tip4PFF_scope( Tip4PFF_exposer );
        bp::class_< SireMM::Tip4PFF::Components, bp::bases< SireMM::CLJFF::Components > >( "Components" )    
            .def( bp::init< >() )    
            .def( bp::init< SireFF::FFBase const &, SireCAS::Symbols const & >(( bp::arg("ffbase"), bp::arg("symbols") )) );
        bp::class_< SireMM::Tip4PFF::Groups, bp::bases< SireMM::CLJFF::Groups > >( "Groups" )    
            .def( bp::init< >() );
        bp::class_< SireMM::Tip4PFF::Parameters, bp::bases< SireMM::CLJFF::Parameters > >( "Parameters" )    
            .def( bp::init< >() );
        Tip4PFF_exposer.def( bp::init< >() );
        Tip4PFF_exposer.def( bp::init< SireVol::Space const &, SireMM::SwitchingFunction const & >(( bp::arg("space"), bp::arg("switchfunc") )) );
        { //::SireMM::Tip4PFF::add
        
            typedef bool ( ::SireMM::Tip4PFF::*add_function_type )( ::SireMol::Molecule const &,::SireFF::ParameterMap const & ) ;
            
            Tip4PFF_exposer.def( 
                "add"
                , add_function_type( &::SireMM::Tip4PFF::add )
                , ( bp::arg("molecule"), bp::arg("map")=::SireFF::ParameterMap( ) ) );
        
        }
        { //::SireMM::Tip4PFF::change
        
            typedef bool ( ::SireMM::Tip4PFF::*change_function_type )( ::SireMol::Molecule const & ) ;
            
            Tip4PFF_exposer.def( 
                "change"
                , change_function_type( &::SireMM::Tip4PFF::change )
                , ( bp::arg("molecule") ) );
        
        }
        { //::SireMM::Tip4PFF::change
        
            typedef bool ( ::SireMM::Tip4PFF::*change_function_type )( ::SireMol::Residue const & ) ;
            
            Tip4PFF_exposer.def( 
                "change"
                , change_function_type( &::SireMM::Tip4PFF::change )
                , ( bp::arg("residue") ) );
        
        }
        { //::SireMM::Tip4PFF::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            
            Tip4PFF_exposer.def( 
                "typeName"
                , typeName_function_type( &::SireMM::Tip4PFF::typeName ) );
        
        }
        { //::SireMM::Tip4PFF::what
        
            typedef char const * ( ::SireMM::Tip4PFF::*what_function_type )(  ) const;
            
            Tip4PFF_exposer.def( 
                "what"
                , what_function_type( &::SireMM::Tip4PFF::what ) );
        
        }
        Tip4PFF_exposer.staticmethod( "typeName" );
    }

}
