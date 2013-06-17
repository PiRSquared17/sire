// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "NullPerturbation.pypp.hpp"

namespace bp = boost::python;

#include "SireCAS/identities.h"

#include "SireCAS/values.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "geometryperturbation.h"

#include "molecule.h"

#include "moleditor.h"

#include "mover.hpp"

#include "perturbation.h"

#include "perturbation.h"

SireMol::NullPerturbation __copy__(const SireMol::NullPerturbation &other){ return SireMol::NullPerturbation(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_NullPerturbation_class(){

    { //::SireMol::NullPerturbation
        typedef bp::class_< SireMol::NullPerturbation, bp::bases< SireMol::Perturbation, SireBase::Property > > NullPerturbation_exposer_t;
        NullPerturbation_exposer_t NullPerturbation_exposer = NullPerturbation_exposer_t( "NullPerturbation", bp::init< >() );
        bp::scope NullPerturbation_scope( NullPerturbation_exposer );
        NullPerturbation_exposer.def( bp::init< SireMol::NullPerturbation const & >(( bp::arg("other") )) );
        NullPerturbation_exposer.def( bp::self != bp::self );
        { //::SireMol::NullPerturbation::operator=
        
            typedef ::SireMol::NullPerturbation & ( ::SireMol::NullPerturbation::*assign_function_type )( ::SireMol::NullPerturbation const & ) ;
            assign_function_type assign_function_value( &::SireMol::NullPerturbation::operator= );
            
            NullPerturbation_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        NullPerturbation_exposer.def( bp::self == bp::self );
        { //::SireMol::NullPerturbation::requiredProperties
        
            typedef ::QSet< QString > ( ::SireMol::NullPerturbation::*requiredProperties_function_type )(  ) const;
            requiredProperties_function_type requiredProperties_function_value( &::SireMol::NullPerturbation::requiredProperties );
            
            NullPerturbation_exposer.def( 
                "requiredProperties"
                , requiredProperties_function_value );
        
        }
        { //::SireMol::NullPerturbation::requiredSymbols
        
            typedef ::QSet< SireCAS::Symbol > ( ::SireMol::NullPerturbation::*requiredSymbols_function_type )(  ) const;
            requiredSymbols_function_type requiredSymbols_function_value( &::SireMol::NullPerturbation::requiredSymbols );
            
            NullPerturbation_exposer.def( 
                "requiredSymbols"
                , requiredSymbols_function_value );
        
        }
        { //::SireMol::NullPerturbation::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::NullPerturbation::typeName );
            
            NullPerturbation_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireMol::NullPerturbation::wouldChange
        
            typedef bool ( ::SireMol::NullPerturbation::*wouldChange_function_type )( ::SireMol::Molecule const &,::SireCAS::Values const & ) const;
            wouldChange_function_type wouldChange_function_value( &::SireMol::NullPerturbation::wouldChange );
            
            NullPerturbation_exposer.def( 
                "wouldChange"
                , wouldChange_function_value
                , ( bp::arg("molecule"), bp::arg("values") ) );
        
        }
        NullPerturbation_exposer.staticmethod( "typeName" );
        NullPerturbation_exposer.def( "__copy__", &__copy__);
        NullPerturbation_exposer.def( "__deepcopy__", &__copy__);
        NullPerturbation_exposer.def( "clone", &__copy__);
        NullPerturbation_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::NullPerturbation >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        NullPerturbation_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::NullPerturbation >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        NullPerturbation_exposer.def( "__str__", &__str__< ::SireMol::NullPerturbation > );
        NullPerturbation_exposer.def( "__repr__", &__str__< ::SireMol::NullPerturbation > );
    }

}
