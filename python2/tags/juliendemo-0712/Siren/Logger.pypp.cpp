// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Logger.pypp.hpp"

namespace bp = boost::python;

#include "Siren/errors.h"

#include "hanref.h"

#include "logger.h"

#include <QBuffer>

#include <QFile>

#include <QProcess>

#include "logger.h"

Siren::Logger __copy__(const Siren::Logger &other){ return Siren::Logger(other); }

#include "Helpers/str.hpp"

void register_Logger_class(){

    { //::Siren::Logger
        typedef bp::class_< Siren::Logger, bp::bases< Siren::Handle > > Logger_exposer_t;
        Logger_exposer_t Logger_exposer = Logger_exposer_t( "Logger" );
        bp::scope Logger_scope( Logger_exposer );
        Logger_exposer.def( bp::init< >() );
        Logger_exposer.def( bp::init< QTextStream * >(( bp::arg("stream") )) );
        Logger_exposer.def( bp::init< Siren::Logger const & >(( bp::arg("other") )) );
        { //::Siren::Logger::hashCode
        
            typedef ::uint ( ::Siren::Logger::*hashCode_function_type )(  ) const;
            hashCode_function_type hashCode_function_value( &::Siren::Logger::hashCode );
            
            Logger_exposer.def( 
                "hashCode"
                , hashCode_function_value );
        
        }
        Logger_exposer.def( bp::self != bp::self );
        { //::Siren::Logger::operator=
        
            typedef ::Siren::Logger & ( ::Siren::Logger::*assign_function_type )( ::Siren::Logger const & ) ;
            assign_function_type assign_function_value( &::Siren::Logger::operator= );
            
            Logger_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        Logger_exposer.def( bp::self == bp::self );
        { //::Siren::Logger::toString
        
            typedef ::QString ( ::Siren::Logger::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::Siren::Logger::toString );
            
            Logger_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::Siren::Logger::write
        
            typedef void ( ::Siren::Logger::*write_function_type )( ::QString const & ) ;
            write_function_type write_function_value( &::Siren::Logger::write );
            
            Logger_exposer.def( 
                "write"
                , write_function_value
                , ( bp::arg("text") ) );
        
        }
        { //::Siren::Logger::writeHeader
        
            typedef void ( ::Siren::Logger::*writeHeader_function_type )( ::QString const &,int ) ;
            writeHeader_function_type writeHeader_function_value( &::Siren::Logger::writeHeader );
            
            Logger_exposer.def( 
                "writeHeader"
                , writeHeader_function_value
                , ( bp::arg("text"), bp::arg("level") ) );
        
        }
        Logger_exposer.def( "__copy__", &__copy__);
        Logger_exposer.def( "__deepcopy__", &__copy__);
        Logger_exposer.def( "clone", &__copy__);
        Logger_exposer.def( "__str__", &__str__< ::Siren::Logger > );
        Logger_exposer.def( "__repr__", &__str__< ::Siren::Logger > );
    }

}
