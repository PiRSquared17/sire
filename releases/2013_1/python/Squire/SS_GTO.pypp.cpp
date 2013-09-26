// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "SS_GTO.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireMaths/boys.h"

#include "SireMaths/maths.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "pointcharge.h"

#include "pointdipole.h"

#include "sgto.h"

#include "sgto.h"

Squire::SS_GTO __copy__(const Squire::SS_GTO &other){ return Squire::SS_GTO(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_SS_GTO_class(){

    { //::Squire::SS_GTO
        typedef bp::class_< Squire::SS_GTO, bp::bases< Squire::GTOPair, Squire::ShellPair, SireBase::Property > > SS_GTO_exposer_t;
        SS_GTO_exposer_t SS_GTO_exposer = SS_GTO_exposer_t( "SS_GTO", bp::init< >() );
        bp::scope SS_GTO_scope( SS_GTO_exposer );
        SS_GTO_exposer.def( bp::init< SireMaths::Vector const &, Squire::S_GTO const &, SireMaths::Vector const &, Squire::S_GTO const & >(( bp::arg("A"), bp::arg("a"), bp::arg("B"), bp::arg("b") )) );
        SS_GTO_exposer.def( bp::init< Squire::SS_GTO const & >(( bp::arg("other") )) );
        { //::Squire::SS_GTO::angularMomentum0
        
            typedef int ( ::Squire::SS_GTO::*angularMomentum0_function_type )(  ) const;
            angularMomentum0_function_type angularMomentum0_function_value( &::Squire::SS_GTO::angularMomentum0 );
            
            SS_GTO_exposer.def( 
                "angularMomentum0"
                , angularMomentum0_function_value );
        
        }
        { //::Squire::SS_GTO::angularMomentum1
        
            typedef int ( ::Squire::SS_GTO::*angularMomentum1_function_type )(  ) const;
            angularMomentum1_function_type angularMomentum1_function_value( &::Squire::SS_GTO::angularMomentum1 );
            
            SS_GTO_exposer.def( 
                "angularMomentum1"
                , angularMomentum1_function_value );
        
        }
        { //::Squire::SS_GTO::nOrbitals0
        
            typedef int ( ::Squire::SS_GTO::*nOrbitals0_function_type )(  ) const;
            nOrbitals0_function_type nOrbitals0_function_value( &::Squire::SS_GTO::nOrbitals0 );
            
            SS_GTO_exposer.def( 
                "nOrbitals0"
                , nOrbitals0_function_value );
        
        }
        { //::Squire::SS_GTO::nOrbitals1
        
            typedef int ( ::Squire::SS_GTO::*nOrbitals1_function_type )(  ) const;
            nOrbitals1_function_type nOrbitals1_function_value( &::Squire::SS_GTO::nOrbitals1 );
            
            SS_GTO_exposer.def( 
                "nOrbitals1"
                , nOrbitals1_function_value );
        
        }
        SS_GTO_exposer.def( bp::self != bp::self );
        { //::Squire::SS_GTO::operator=
        
            typedef ::Squire::SS_GTO & ( ::Squire::SS_GTO::*assign_function_type )( ::Squire::SS_GTO const & ) ;
            assign_function_type assign_function_value( &::Squire::SS_GTO::operator= );
            
            SS_GTO_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        SS_GTO_exposer.def( bp::self == bp::self );
        { //::Squire::SS_GTO::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::Squire::SS_GTO::typeName );
            
            SS_GTO_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        SS_GTO_exposer.staticmethod( "typeName" );
        SS_GTO_exposer.def( "__copy__", &__copy__);
        SS_GTO_exposer.def( "__deepcopy__", &__copy__);
        SS_GTO_exposer.def( "clone", &__copy__);
        SS_GTO_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::Squire::SS_GTO >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        SS_GTO_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::Squire::SS_GTO >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        SS_GTO_exposer.def( "__str__", &__str__< ::Squire::SS_GTO > );
        SS_GTO_exposer.def( "__repr__", &__str__< ::Squire::SS_GTO > );
    }

}
