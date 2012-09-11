// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "SysIdx.pypp.hpp"

namespace bp = boost::python;

#include "sysidx.h"

#include "sysidx.h"

SireSystem::SysIdx __copy__(const SireSystem::SysIdx &other){ return SireSystem::SysIdx(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_SysIdx_class(){

    { //::SireSystem::SysIdx
        typedef bp::class_< SireSystem::SysIdx, bp::bases< SireSystem::SysID, SireID::ID, SireID::IndexBase > > SysIdx_exposer_t;
        SysIdx_exposer_t SysIdx_exposer = SysIdx_exposer_t( "SysIdx", bp::init< >() );
        bp::scope SysIdx_scope( SysIdx_exposer );
        SysIdx_exposer.def( bp::init< qint32 >(( bp::arg("idx") )) );
        SysIdx_exposer.def( bp::init< SireSystem::SysIdx const & >(( bp::arg("other") )) );
        { //::SireSystem::SysIdx::hash
        
            typedef ::uint ( ::SireSystem::SysIdx::*hash_function_type )(  ) const;
            hash_function_type hash_function_value( &::SireSystem::SysIdx::hash );
            
            SysIdx_exposer.def( 
                "hash"
                , hash_function_value );
        
        }
        { //::SireSystem::SysIdx::isNull
        
            typedef bool ( ::SireSystem::SysIdx::*isNull_function_type )(  ) const;
            isNull_function_type isNull_function_value( &::SireSystem::SysIdx::isNull );
            
            SysIdx_exposer.def( 
                "isNull"
                , isNull_function_value );
        
        }
        { //::SireSystem::SysIdx::map
        
            typedef ::QList< SireSystem::SysIdx > ( ::SireSystem::SysIdx::*map_function_type )( ::SireSystem::Systems const & ) const;
            map_function_type map_function_value( &::SireSystem::SysIdx::map );
            
            SysIdx_exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("systems") ) );
        
        }
        { //::SireSystem::SysIdx::null
        
            typedef ::SireSystem::SysIdx ( *null_function_type )(  );
            null_function_type null_function_value( &::SireSystem::SysIdx::null );
            
            SysIdx_exposer.def( 
                "null"
                , null_function_value );
        
        }
        { //::SireSystem::SysIdx::operator=
        
            typedef ::SireSystem::SysIdx & ( ::SireSystem::SysIdx::*assign_function_type )( ::SireSystem::SysIdx const & ) ;
            assign_function_type assign_function_value( &::SireSystem::SysIdx::operator= );
            
            SysIdx_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        { //::SireSystem::SysIdx::toString
        
            typedef ::QString ( ::SireSystem::SysIdx::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireSystem::SysIdx::toString );
            
            SysIdx_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireSystem::SysIdx::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireSystem::SysIdx::typeName );
            
            SysIdx_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireSystem::SysIdx::what
        
            typedef char const * ( ::SireSystem::SysIdx::*what_function_type )(  ) const;
            what_function_type what_function_value( &::SireSystem::SysIdx::what );
            
            SysIdx_exposer.def( 
                "what"
                , what_function_value );
        
        }
        SysIdx_exposer.staticmethod( "null" );
        SysIdx_exposer.staticmethod( "typeName" );
        SysIdx_exposer.def( "__copy__", &__copy__);
        SysIdx_exposer.def( "__deepcopy__", &__copy__);
        SysIdx_exposer.def( "clone", &__copy__);
        SysIdx_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireSystem::SysIdx >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        SysIdx_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireSystem::SysIdx >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        SysIdx_exposer.def( "__str__", &__str__< ::SireSystem::SysIdx > );
        SysIdx_exposer.def( "__repr__", &__str__< ::SireSystem::SysIdx > );
    }

}
