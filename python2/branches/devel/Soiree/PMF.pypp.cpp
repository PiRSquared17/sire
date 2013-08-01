// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "PMF.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireMaths/maths.h"

#include "SireStream/shareddatastream.h"

#include "fep.h"

#include "tostring.h"

#include "fep.h"

Soiree::PMF __copy__(const Soiree::PMF &other){ return Soiree::PMF(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_PMF_class(){

    { //::Soiree::PMF
        typedef bp::class_< Soiree::PMF, bp::bases< SireBase::Property > > PMF_exposer_t;
        PMF_exposer_t PMF_exposer = PMF_exposer_t( "PMF", bp::init< >() );
        bp::scope PMF_scope( PMF_exposer );
        PMF_exposer.def( bp::init< QVector< Soiree::DataPoint > const & >(( bp::arg("values") )) );
        PMF_exposer.def( bp::init< Soiree::PMF const & >(( bp::arg("other") )) );
        { //::Soiree::PMF::deltaG
        
            typedef double ( ::Soiree::PMF::*deltaG_function_type )(  ) const;
            deltaG_function_type deltaG_function_value( &::Soiree::PMF::deltaG );
            
            PMF_exposer.def( 
                "deltaG"
                , deltaG_function_value );
        
        }
        { //::Soiree::PMF::error
        
            typedef double ( ::Soiree::PMF::*error_function_type )(  ) const;
            error_function_type error_function_value( &::Soiree::PMF::error );
            
            PMF_exposer.def( 
                "error"
                , error_function_value );
        
        }
        { //::Soiree::PMF::isEmpty
        
            typedef bool ( ::Soiree::PMF::*isEmpty_function_type )(  ) const;
            isEmpty_function_type isEmpty_function_value( &::Soiree::PMF::isEmpty );
            
            PMF_exposer.def( 
                "isEmpty"
                , isEmpty_function_value );
        
        }
        PMF_exposer.def( bp::self != bp::self );
        { //::Soiree::PMF::operator=
        
            typedef ::Soiree::PMF & ( ::Soiree::PMF::*assign_function_type )( ::Soiree::PMF const & ) ;
            assign_function_type assign_function_value( &::Soiree::PMF::operator= );
            
            PMF_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        PMF_exposer.def( bp::self == bp::self );
        { //::Soiree::PMF::rangeMax
        
            typedef double ( ::Soiree::PMF::*rangeMax_function_type )(  ) const;
            rangeMax_function_type rangeMax_function_value( &::Soiree::PMF::rangeMax );
            
            PMF_exposer.def( 
                "rangeMax"
                , rangeMax_function_value );
        
        }
        { //::Soiree::PMF::rangeMin
        
            typedef double ( ::Soiree::PMF::*rangeMin_function_type )(  ) const;
            rangeMin_function_type rangeMin_function_value( &::Soiree::PMF::rangeMin );
            
            PMF_exposer.def( 
                "rangeMin"
                , rangeMin_function_value );
        
        }
        { //::Soiree::PMF::toString
        
            typedef ::QString ( ::Soiree::PMF::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::Soiree::PMF::toString );
            
            PMF_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::Soiree::PMF::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::Soiree::PMF::typeName );
            
            PMF_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::Soiree::PMF::values
        
            typedef ::QVector< Soiree::DataPoint > ( ::Soiree::PMF::*values_function_type )(  ) const;
            values_function_type values_function_value( &::Soiree::PMF::values );
            
            PMF_exposer.def( 
                "values"
                , values_function_value );
        
        }
        { //::Soiree::PMF::what
        
            typedef char const * ( ::Soiree::PMF::*what_function_type )(  ) const;
            what_function_type what_function_value( &::Soiree::PMF::what );
            
            PMF_exposer.def( 
                "what"
                , what_function_value );
        
        }
        PMF_exposer.staticmethod( "typeName" );
        PMF_exposer.def( "__copy__", &__copy__);
        PMF_exposer.def( "__deepcopy__", &__copy__);
        PMF_exposer.def( "clone", &__copy__);
        PMF_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::Soiree::PMF >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        PMF_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::Soiree::PMF >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        PMF_exposer.def( "__str__", &__str__< ::Soiree::PMF > );
        PMF_exposer.def( "__repr__", &__str__< ::Soiree::PMF > );
    }

}
