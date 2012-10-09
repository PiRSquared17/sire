// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Options.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireSim/option.h"

#include "option.h"

#include "tostring.h"

#include <QDebug>

#include <QDir>

#include <QDomDocument>

#include <QDomElement>

#include <QDomNode>

#include <QFile>

#include <QFileInfo>

#include <QHash>

#include <QList>

#include "option.h"

SireSim::Options __copy__(const SireSim::Options &other){ return SireSim::Options(other); }

#include "Helpers/str.hpp"

void register_Options_class(){

    { //::SireSim::Options
        typedef bp::class_< SireSim::Options, bp::bases< SireSim::Value > > Options_exposer_t;
        Options_exposer_t Options_exposer = Options_exposer_t( "Options", bp::init< >() );
        bp::scope Options_scope( Options_exposer );
        Options_exposer.def( bp::init< QString, bp::optional< QStringList > >(( bp::arg("xml"), bp::arg("path")=::QStringList( ) )) );
        Options_exposer.def( bp::init< QDomElement, bp::optional< QStringList > >(( bp::arg("elem"), bp::arg("path")=::QStringList( ) )) );
        Options_exposer.def( bp::init< SireSim::Option const & >(( bp::arg("option") )) );
        Options_exposer.def( bp::init< QList< SireSim::Option > const &, bp::optional< bool > >(( bp::arg("options"), bp::arg("mutually_exclusive")=(bool)(false) )) );
        Options_exposer.def( bp::init< SireSim::Options const & >(( bp::arg("other") )) );
        { //::SireSim::Options::add
        
            typedef ::SireSim::Options ( ::SireSim::Options::*add_function_type )( ::SireSim::Options const & ) const;
            add_function_type add_function_value( &::SireSim::Options::add );
            
            Options_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("other") ) );
        
        }
        { //::SireSim::Options::clear
        
            typedef ::SireSim::ValuePtr ( ::SireSim::Options::*clear_function_type )(  ) const;
            clear_function_type clear_function_value( &::SireSim::Options::clear );
            
            Options_exposer.def( 
                "clear"
                , clear_function_value );
        
        }
        { //::SireSim::Options::fromConfig
        
            typedef ::SireSim::Options ( ::SireSim::Options::*fromConfig_function_type )( ::QString ) const;
            fromConfig_function_type fromConfig_function_value( &::SireSim::Options::fromConfig );
            
            Options_exposer.def( 
                "fromConfig"
                , fromConfig_function_value
                , ( bp::arg("text") ) );
        
        }
        { //::SireSim::Options::fromConfigFile
        
            typedef ::SireSim::Options ( ::SireSim::Options::*fromConfigFile_function_type )( ::QString ) const;
            fromConfigFile_function_type fromConfigFile_function_value( &::SireSim::Options::fromConfigFile );
            
            Options_exposer.def( 
                "fromConfigFile"
                , fromConfigFile_function_value
                , ( bp::arg("configfile") ) );
        
        }
        { //::SireSim::Options::fromXML
        
            typedef ::SireSim::Options ( *fromXML_function_type )( ::QString,::QStringList );
            fromXML_function_type fromXML_function_value( &::SireSim::Options::fromXML );
            
            Options_exposer.def( 
                "fromXML"
                , fromXML_function_value
                , ( bp::arg("xml"), bp::arg("path")=::QStringList( ) ) );
        
        }
        { //::SireSim::Options::fromXMLFile
        
            typedef ::SireSim::Options ( *fromXMLFile_function_type )( ::QString,::QStringList );
            fromXMLFile_function_type fromXMLFile_function_value( &::SireSim::Options::fromXMLFile );
            
            Options_exposer.def( 
                "fromXMLFile"
                , fromXMLFile_function_value
                , ( bp::arg("xmlfile"), bp::arg("path")=::QStringList( ) ) );
        
        }
        { //::SireSim::Options::getNestedOption
        
            typedef ::SireSim::Option ( ::SireSim::Options::*getNestedOption_function_type )( ::QString ) const;
            getNestedOption_function_type getNestedOption_function_value( &::SireSim::Options::getNestedOption );
            
            Options_exposer.def( 
                "getNestedOption"
                , getNestedOption_function_value
                , ( bp::arg("key") ) );
        
        }
        { //::SireSim::Options::getNestedValue
        
            typedef ::SireSim::ValuePtr ( ::SireSim::Options::*getNestedValue_function_type )( ::QString ) const;
            getNestedValue_function_type getNestedValue_function_value( &::SireSim::Options::getNestedValue );
            
            Options_exposer.def( 
                "getNestedValue"
                , getNestedValue_function_value
                , ( bp::arg("key") ) );
        
        }
        { //::SireSim::Options::hasSubOptions
        
            typedef bool ( ::SireSim::Options::*hasSubOptions_function_type )(  ) const;
            hasSubOptions_function_type hasSubOptions_function_value( &::SireSim::Options::hasSubOptions );
            
            Options_exposer.def( 
                "hasSubOptions"
                , hasSubOptions_function_value );
        
        }
        { //::SireSim::Options::hasValue
        
            typedef bool ( ::SireSim::Options::*hasValue_function_type )( ::QString ) const;
            hasValue_function_type hasValue_function_value( &::SireSim::Options::hasValue );
            
            Options_exposer.def( 
                "hasValue"
                , hasValue_function_value
                , ( bp::arg("key") ) );
        
        }
        { //::SireSim::Options::hasValue
        
            typedef bool ( ::SireSim::Options::*hasValue_function_type )(  ) const;
            hasValue_function_type hasValue_function_value( &::SireSim::Options::hasValue );
            
            Options_exposer.def( 
                "hasValue"
                , hasValue_function_value );
        
        }
        { //::SireSim::Options::keys
        
            typedef ::QStringList ( ::SireSim::Options::*keys_function_type )(  ) const;
            keys_function_type keys_function_value( &::SireSim::Options::keys );
            
            Options_exposer.def( 
                "keys"
                , keys_function_value );
        
        }
        { //::SireSim::Options::keysWithValue
        
            typedef ::QStringList ( ::SireSim::Options::*keysWithValue_function_type )(  ) const;
            keysWithValue_function_type keysWithValue_function_value( &::SireSim::Options::keysWithValue );
            
            Options_exposer.def( 
                "keysWithValue"
                , keysWithValue_function_value );
        
        }
        Options_exposer.def( bp::self != bp::self );
        Options_exposer.def( bp::self + bp::self );
        { //::SireSim::Options::operator=
        
            typedef ::SireSim::Options & ( ::SireSim::Options::*assign_function_type )( ::SireSim::Options const & ) ;
            assign_function_type assign_function_value( &::SireSim::Options::operator= );
            
            Options_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        Options_exposer.def( bp::self == bp::self );
        { //::SireSim::Options::operator[]
        
            typedef ::SireSim::Option ( ::SireSim::Options::*__getitem___function_type )( ::QString ) const;
            __getitem___function_type __getitem___function_value( &::SireSim::Options::operator[] );
            
            Options_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("key") ) );
        
        }
        { //::SireSim::Options::options
        
            typedef ::QList< SireSim::Option > ( ::SireSim::Options::*options_function_type )(  ) const;
            options_function_type options_function_value( &::SireSim::Options::options );
            
            Options_exposer.def( 
                "options"
                , options_function_value );
        
        }
        { //::SireSim::Options::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireSim::Options::typeName );
            
            Options_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireSim::Options::what
        
            typedef char const * ( ::SireSim::Options::*what_function_type )(  ) const;
            what_function_type what_function_value( &::SireSim::Options::what );
            
            Options_exposer.def( 
                "what"
                , what_function_value );
        
        }
        Options_exposer.staticmethod( "fromXML" );
        Options_exposer.staticmethod( "fromXMLFile" );
        Options_exposer.staticmethod( "typeName" );
        Options_exposer.def( "__copy__", &__copy__);
        Options_exposer.def( "__deepcopy__", &__copy__);
        Options_exposer.def( "clone", &__copy__);
        Options_exposer.def( "__str__", &__str__< ::SireSim::Options > );
        Options_exposer.def( "__repr__", &__str__< ::SireSim::Options > );
    }

}