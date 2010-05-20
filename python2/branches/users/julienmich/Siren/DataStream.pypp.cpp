// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "DataStream.pypp.hpp"

namespace bp = boost::python;

#include "Siren/errors.h"

#include "datastream.h"

#include "hanref.h"

#include "logger.h"

#include "tester.h"

#include <QHash>

#include <QMutex>

#include "datastream.h"

Siren::DataStream __copy__(const Siren::DataStream &other){ return Siren::DataStream(other); }

#include "Helpers/str.hpp"

void register_DataStream_class(){

    { //::Siren::DataStream
        typedef bp::class_< Siren::DataStream, bp::bases< Siren::Stream, Siren::Handle > > DataStream_exposer_t;
        DataStream_exposer_t DataStream_exposer = DataStream_exposer_t( "DataStream" );
        bp::scope DataStream_scope( DataStream_exposer );
        DataStream_exposer.def( bp::init< >() );
        DataStream_exposer.def( bp::init< QIODevice * >(( bp::arg("d") )) );
        DataStream_exposer.def( bp::init< QByteArray *, QFlags< QIODevice::OpenModeFlag > >(( bp::arg("a"), bp::arg("mode") )) );
        DataStream_exposer.def( bp::init< QByteArray const & >(( bp::arg("a") )) );
        DataStream_exposer.def( bp::init< QDataStream & >(( bp::arg("ds") )) );
        DataStream_exposer.def( bp::init< Siren::DataStream const & >(( bp::arg("other") )) );
        { //::Siren::DataStream::hashCode
        
            typedef ::uint ( ::Siren::DataStream::*hashCode_function_type )(  ) const;
            hashCode_function_type hashCode_function_value( &::Siren::DataStream::hashCode );
            
            DataStream_exposer.def( 
                "hashCode"
                , hashCode_function_value );
        
        }
        DataStream_exposer.def( bp::self != bp::self );
        { //::Siren::DataStream::operator=
        
            typedef ::Siren::DataStream & ( ::Siren::DataStream::*assign_function_type )( ::Siren::DataStream const & ) ;
            assign_function_type assign_function_value( &::Siren::DataStream::operator= );
            
            DataStream_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        DataStream_exposer.def( bp::self == bp::self );
        { //::Siren::DataStream::toString
        
            typedef ::QString ( ::Siren::DataStream::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::Siren::DataStream::toString );
            
            DataStream_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        DataStream_exposer.def( "__copy__", &__copy__);
        DataStream_exposer.def( "__deepcopy__", &__copy__);
        DataStream_exposer.def( "clone", &__copy__);
        DataStream_exposer.def( "__str__", &__str__< ::Siren::DataStream > );
        DataStream_exposer.def( "__repr__", &__str__< ::Siren::DataStream > );
    }

}
