// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "SupraMove.pypp.hpp"

namespace bp = boost::python;

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "supramove.h"

#include "suprasystem.h"

#include "supramove.h"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_SupraMove_class(){

    { //::SireMove::SupraMove
        typedef bp::class_< SireMove::SupraMove, bp::bases< SireBase::Property >, boost::noncopyable > SupraMove_exposer_t;
        SupraMove_exposer_t SupraMove_exposer = SupraMove_exposer_t( "SupraMove", bp::no_init );
        bp::scope SupraMove_scope( SupraMove_exposer );
        { //::SireMove::SupraMove::clearStatistics
        
            typedef void ( ::SireMove::SupraMove::*clearStatistics_function_type )(  ) ;
            clearStatistics_function_type clearStatistics_function_value( &::SireMove::SupraMove::clearStatistics );
            
            SupraMove_exposer.def( 
                "clearStatistics"
                , clearStatistics_function_value );
        
        }
        { //::SireMove::SupraMove::move
        
            typedef void ( ::SireMove::SupraMove::*move_function_type )( ::SireMove::SupraSystem &,int,bool ) ;
            move_function_type move_function_value( &::SireMove::SupraMove::move );
            
            SupraMove_exposer.def( 
                "move"
                , move_function_value
                , ( bp::arg("system"), bp::arg("nmoves"), bp::arg("record_stats")=(bool)(true) ) );
        
        }
        { //::SireMove::SupraMove::nMoves
        
            typedef int ( ::SireMove::SupraMove::*nMoves_function_type )(  ) const;
            nMoves_function_type nMoves_function_value( &::SireMove::SupraMove::nMoves );
            
            SupraMove_exposer.def( 
                "nMoves"
                , nMoves_function_value );
        
        }
        { //::SireMove::SupraMove::null
        
            typedef ::SireMove::NullSupraMove const & ( *null_function_type )(  );
            null_function_type null_function_value( &::SireMove::SupraMove::null );
            
            SupraMove_exposer.def( 
                "null"
                , null_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMove::SupraMove::toString
        
            typedef ::QString ( ::SireMove::SupraMove::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireMove::SupraMove::toString );
            
            SupraMove_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireMove::SupraMove::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMove::SupraMove::typeName );
            
            SupraMove_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        SupraMove_exposer.staticmethod( "null" );
        SupraMove_exposer.staticmethod( "typeName" );
        SupraMove_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMove::SupraMove >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        SupraMove_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMove::SupraMove >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        SupraMove_exposer.def( "__str__", &__str__< ::SireMove::SupraMove > );
        SupraMove_exposer.def( "__repr__", &__str__< ::SireMove::SupraMove > );
    }

}
