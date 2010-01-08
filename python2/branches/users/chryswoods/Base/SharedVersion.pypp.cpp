// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "SharedVersion.pypp.hpp"

namespace bp = boost::python;

#include "sharedversion.h"

#include "versionregistry.hpp"

#include "sharedversion.h"

#include "Siren/logger.h"

SireBase::SharedVersion __copy__(const SireBase::SharedVersion &other){ return SireBase::SharedVersion(other); }

#include "Siren/str.hpp"

void register_SharedVersion_class(){

    { //::SireBase::SharedVersion
        typedef bp::class_< SireBase::SharedVersion, bp::bases< Siren::Handle > > SharedVersion_exposer_t;
        SharedVersion_exposer_t SharedVersion_exposer = SharedVersion_exposer_t( "SharedVersion" );
        bp::scope SharedVersion_scope( SharedVersion_exposer );
        SharedVersion_exposer.def( bp::init< >() );
        SharedVersion_exposer.def( bp::init< quint64, quint64 >(( bp::arg("vmaj"), bp::arg("vmin") )) );
        SharedVersion_exposer.def( bp::init< SireBase::SharedVersion const & >(( bp::arg("other") )) );
        { //::SireBase::SharedVersion::hashCode
        
            typedef ::uint ( ::SireBase::SharedVersion::*hashCode_function_type )(  ) const;
            hashCode_function_type hashCode_function_value( &::SireBase::SharedVersion::hashCode );
            
            SharedVersion_exposer.def( 
                "hashCode"
                , hashCode_function_value );
        
        }
        { //::SireBase::SharedVersion::incrementMajor
        
            typedef void ( ::SireBase::SharedVersion::*incrementMajor_function_type )(  ) ;
            incrementMajor_function_type incrementMajor_function_value( &::SireBase::SharedVersion::incrementMajor );
            
            SharedVersion_exposer.def( 
                "incrementMajor"
                , incrementMajor_function_value );
        
        }
        { //::SireBase::SharedVersion::incrementMinor
        
            typedef void ( ::SireBase::SharedVersion::*incrementMinor_function_type )(  ) ;
            incrementMinor_function_type incrementMinor_function_value( &::SireBase::SharedVersion::incrementMinor );
            
            SharedVersion_exposer.def( 
                "incrementMinor"
                , incrementMinor_function_value );
        
        }
        { //::SireBase::SharedVersion::majorVersion
        
            typedef ::quint64 ( ::SireBase::SharedVersion::*majorVersion_function_type )(  ) const;
            majorVersion_function_type majorVersion_function_value( &::SireBase::SharedVersion::majorVersion );
            
            SharedVersion_exposer.def( 
                "majorVersion"
                , majorVersion_function_value );
        
        }
        { //::SireBase::SharedVersion::minorVersion
        
            typedef ::quint64 ( ::SireBase::SharedVersion::*minorVersion_function_type )(  ) const;
            minorVersion_function_type minorVersion_function_value( &::SireBase::SharedVersion::minorVersion );
            
            SharedVersion_exposer.def( 
                "minorVersion"
                , minorVersion_function_value );
        
        }
        SharedVersion_exposer.def( bp::self != bp::self );
        { //::SireBase::SharedVersion::operator=
        
            typedef ::SireBase::SharedVersion & ( ::SireBase::SharedVersion::*assign_function_type )( ::SireBase::SharedVersion const & ) ;
            assign_function_type assign_function_value( &::SireBase::SharedVersion::operator= );
            
            SharedVersion_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        SharedVersion_exposer.def( bp::self == bp::self );
        { //::SireBase::SharedVersion::toString
        
            typedef ::QString ( ::SireBase::SharedVersion::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireBase::SharedVersion::toString );
            
            SharedVersion_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireBase::SharedVersion::version
        
            typedef ::SireBase::Version const & ( ::SireBase::SharedVersion::*version_function_type )(  ) const;
            version_function_type version_function_value( &::SireBase::SharedVersion::version );
            
            SharedVersion_exposer.def( 
                "version"
                , version_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        SharedVersion_exposer.def( "__copy__", &__copy__);
        SharedVersion_exposer.def( "__deepcopy__", &__copy__);
        SharedVersion_exposer.def( "clone", &__copy__);
        SharedVersion_exposer.def( "__str__", &__str__< ::SireBase::SharedVersion > );
        SharedVersion_exposer.def( "__repr__", &__str__< ::SireBase::SharedVersion > );
    }

}
