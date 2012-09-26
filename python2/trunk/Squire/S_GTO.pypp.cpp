// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "S_GTO.pypp.hpp"

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

Squire::S_GTO __copy__(const Squire::S_GTO &other){ return Squire::S_GTO(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_S_GTO_class(){

    { //::Squire::S_GTO
        typedef bp::class_< Squire::S_GTO, bp::bases< Squire::GTO, Squire::OrbitalShell, Squire::Orbital, SireBase::Property > > S_GTO_exposer_t;
        S_GTO_exposer_t S_GTO_exposer = S_GTO_exposer_t( "S_GTO", bp::init< >() );
        bp::scope S_GTO_scope( S_GTO_exposer );
        S_GTO_exposer.def( bp::init< double, bp::optional< double > >(( bp::arg("alpha"), bp::arg("scale")=1 )) );
        S_GTO_exposer.def( bp::init< Squire::S_GTO const & >(( bp::arg("other") )) );
        { //::Squire::S_GTO::angularMomentum
        
            typedef int ( ::Squire::S_GTO::*angularMomentum_function_type )(  ) const;
            angularMomentum_function_type angularMomentum_function_value( &::Squire::S_GTO::angularMomentum );
            
            S_GTO_exposer.def( 
                "angularMomentum"
                , angularMomentum_function_value );
        
        }
        { //::Squire::S_GTO::nOrbitals
        
            typedef int ( ::Squire::S_GTO::*nOrbitals_function_type )(  ) const;
            nOrbitals_function_type nOrbitals_function_value( &::Squire::S_GTO::nOrbitals );
            
            S_GTO_exposer.def( 
                "nOrbitals"
                , nOrbitals_function_value );
        
        }
        S_GTO_exposer.def( bp::self != bp::self );
        { //::Squire::S_GTO::operator=
        
            typedef ::Squire::S_GTO & ( ::Squire::S_GTO::*assign_function_type )( ::Squire::S_GTO const & ) ;
            assign_function_type assign_function_value( &::Squire::S_GTO::operator= );
            
            S_GTO_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        S_GTO_exposer.def( bp::self == bp::self );
        { //::Squire::S_GTO::toString
        
            typedef ::QString ( ::Squire::S_GTO::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::Squire::S_GTO::toString );
            
            S_GTO_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::Squire::S_GTO::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::Squire::S_GTO::typeName );
            
            S_GTO_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        S_GTO_exposer.staticmethod( "typeName" );
        S_GTO_exposer.def( "__copy__", &__copy__);
        S_GTO_exposer.def( "__deepcopy__", &__copy__);
        S_GTO_exposer.def( "clone", &__copy__);
        S_GTO_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::Squire::S_GTO >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        S_GTO_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::Squire::S_GTO >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        S_GTO_exposer.def( "__str__", &__str__< ::Squire::S_GTO > );
        S_GTO_exposer.def( "__repr__", &__str__< ::Squire::S_GTO > );
    }

}
