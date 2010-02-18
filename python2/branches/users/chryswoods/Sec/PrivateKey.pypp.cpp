// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "PrivateKey.pypp.hpp"

namespace bp = boost::python;

#include "SireSec/errors.h"

#include "Siren/errors.h"

#include "Siren/mutex.h"

#include "Siren/stream.h"

#include "privatekey.h"

#include "publickey.h"

#include "pubprilock.h"

#include "signaturelock.h"

#include <QDebug>

#include <QHash>

#include <QThreadStorage>

#include <QUuid>

#include <boost/bind.hpp>

#include "privatekey.h"

#include "Siren/logger.h"

SireSec::PrivateKey __copy__(const SireSec::PrivateKey &other){ return SireSec::PrivateKey(other); }

#include "Siren/str.hpp"

void register_PrivateKey_class(){

    { //::SireSec::PrivateKey
        typedef bp::class_< SireSec::PrivateKey, bp::bases< SireSec::Key, Siren::Object > > PrivateKey_exposer_t;
        PrivateKey_exposer_t PrivateKey_exposer = PrivateKey_exposer_t( "PrivateKey" );
        bp::scope PrivateKey_scope( PrivateKey_exposer );
        PrivateKey_exposer.def( bp::init< >() );
        PrivateKey_exposer.def( bp::init< SireSec::PrivateKey const & >(( bp::arg("other") )) );
        { //::SireSec::PrivateKey::availableToThisThread
        
            typedef bool ( ::SireSec::PrivateKey::*availableToThisThread_function_type )(  ) const;
            availableToThisThread_function_type availableToThisThread_function_value( &::SireSec::PrivateKey::availableToThisThread );
            
            PrivateKey_exposer.def( 
                "availableToThisThread"
                , availableToThisThread_function_value );
        
        }
        { //::SireSec::PrivateKey::generate
        
            typedef ::boost::tuples::tuple< SireSec::PublicKey, SireSec::PrivateKey, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type > ( *generate_function_type )( ::QString,::SireSec::KeyTypes::KeyType,int );
            generate_function_type generate_function_value( &::SireSec::PrivateKey::generate );
            
            PrivateKey_exposer.def( 
                "generate"
                , generate_function_value
                , ( bp::arg("label"), bp::arg("keytype")=::SireSec::KeyTypes::DEFAULT, bp::arg("keylength")=(int)(0) ) );
        
        }
        { //::SireSec::PrivateKey::generate
        
            typedef ::boost::tuples::tuple< SireSec::PublicKey, SireSec::PrivateKey, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type > ( *generate_function_type )( ::SireSec::KeyTypes::KeyType,int );
            generate_function_type generate_function_value( &::SireSec::PrivateKey::generate );
            
            PrivateKey_exposer.def( 
                "generate"
                , generate_function_value
                , ( bp::arg("keytype")=::SireSec::KeyTypes::DEFAULT, bp::arg("keylength")=(int)(0) ) );
        
        }
        { //::SireSec::PrivateKey::generate
        
            typedef ::boost::tuples::tuple< SireSec::PublicKey, SireSec::PrivateKey, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type > ( *generate_function_type )( ::QDateTime const &,::QString,::SireSec::KeyTypes::KeyType,int );
            generate_function_type generate_function_value( &::SireSec::PrivateKey::generate );
            
            PrivateKey_exposer.def( 
                "generate"
                , generate_function_value
                , ( bp::arg("expiry"), bp::arg("label"), bp::arg("keytype")=::SireSec::KeyTypes::DEFAULT, bp::arg("keylength")=(int)(0) ) );
        
        }
        { //::SireSec::PrivateKey::generate
        
            typedef ::boost::tuples::tuple< SireSec::PublicKey, SireSec::PrivateKey, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type > ( *generate_function_type )( ::QFlags< SireSec::Key::Option >,::QString,::SireSec::KeyTypes::KeyType,int );
            generate_function_type generate_function_value( &::SireSec::PrivateKey::generate );
            
            PrivateKey_exposer.def( 
                "generate"
                , generate_function_value
                , ( bp::arg("keyoptions"), bp::arg("label"), bp::arg("keytype")=::SireSec::KeyTypes::DEFAULT, bp::arg("keylength")=(int)(0) ) );
        
        }
        { //::SireSec::PrivateKey::generate
        
            typedef ::boost::tuples::tuple< SireSec::PublicKey, SireSec::PrivateKey, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type > ( *generate_function_type )( ::QFlags< SireSec::Key::Option >,::QDateTime const &,::QString,::SireSec::KeyTypes::KeyType,int );
            generate_function_type generate_function_value( &::SireSec::PrivateKey::generate );
            
            PrivateKey_exposer.def( 
                "generate"
                , generate_function_value
                , ( bp::arg("keyoptions"), bp::arg("expiry"), bp::arg("label"), bp::arg("keytype")=::SireSec::KeyTypes::DEFAULT, bp::arg("keylength")=(int)(0) ) );
        
        }
        { //::SireSec::PrivateKey::generate
        
            typedef ::boost::tuples::tuple< SireSec::PublicKey, SireSec::PrivateKey, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type > ( *generate_function_type )( ::QDateTime const &,::SireSec::KeyTypes::KeyType,int );
            generate_function_type generate_function_value( &::SireSec::PrivateKey::generate );
            
            PrivateKey_exposer.def( 
                "generate"
                , generate_function_value
                , ( bp::arg("expiry"), bp::arg("keytype")=::SireSec::KeyTypes::DEFAULT, bp::arg("keylength")=(int)(0) ) );
        
        }
        { //::SireSec::PrivateKey::generate
        
            typedef ::boost::tuples::tuple< SireSec::PublicKey, SireSec::PrivateKey, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type > ( *generate_function_type )( ::QFlags< SireSec::Key::Option >,::SireSec::KeyTypes::KeyType,int );
            generate_function_type generate_function_value( &::SireSec::PrivateKey::generate );
            
            PrivateKey_exposer.def( 
                "generate"
                , generate_function_value
                , ( bp::arg("keyoptions"), bp::arg("keytype")=::SireSec::KeyTypes::DEFAULT, bp::arg("keylength")=(int)(0) ) );
        
        }
        { //::SireSec::PrivateKey::generate
        
            typedef ::boost::tuples::tuple< SireSec::PublicKey, SireSec::PrivateKey, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type > ( *generate_function_type )( ::QFlags< SireSec::Key::Option >,::QDateTime const &,::SireSec::KeyTypes::KeyType,int );
            generate_function_type generate_function_value( &::SireSec::PrivateKey::generate );
            
            PrivateKey_exposer.def( 
                "generate"
                , generate_function_value
                , ( bp::arg("keyoptions"), bp::arg("expiry"), bp::arg("keytype")=::SireSec::KeyTypes::DEFAULT, bp::arg("keylength")=(int)(0) ) );
        
        }
        { //::SireSec::PrivateKey::hashCode
        
            typedef ::uint ( ::SireSec::PrivateKey::*hashCode_function_type )(  ) const;
            hashCode_function_type hashCode_function_value( &::SireSec::PrivateKey::hashCode );
            
            PrivateKey_exposer.def( 
                "hashCode"
                , hashCode_function_value );
        
        }
        { //::SireSec::PrivateKey::isValid
        
            typedef bool ( ::SireSec::PrivateKey::*isValid_function_type )(  ) const;
            isValid_function_type isValid_function_value( &::SireSec::PrivateKey::isValid );
            
            PrivateKey_exposer.def( 
                "isValid"
                , isValid_function_value );
        
        }
        PrivateKey_exposer.def( bp::self != bp::self );
        { //::SireSec::PrivateKey::operator=
        
            typedef ::SireSec::PrivateKey & ( ::SireSec::PrivateKey::*assign_function_type )( ::SireSec::PrivateKey const & ) ;
            assign_function_type assign_function_value( &::SireSec::PrivateKey::operator= );
            
            PrivateKey_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        PrivateKey_exposer.def( bp::self == bp::self );
        { //::SireSec::PrivateKey::stream
        
            typedef void ( ::SireSec::PrivateKey::*stream_function_type )( ::Siren::Stream & ) ;
            stream_function_type stream_function_value( &::SireSec::PrivateKey::stream );
            
            PrivateKey_exposer.def( 
                "stream"
                , stream_function_value
                , ( bp::arg("s") ) );
        
        }
        { //::SireSec::PrivateKey::toString
        
            typedef ::QString ( ::SireSec::PrivateKey::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireSec::PrivateKey::toString );
            
            PrivateKey_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        PrivateKey_exposer.staticmethod( "generate" );
        PrivateKey_exposer.def( "__copy__", &__copy__);
        PrivateKey_exposer.def( "__deepcopy__", &__copy__);
        PrivateKey_exposer.def( "clone", &__copy__);
        PrivateKey_exposer.def( "__str__", &__str__< ::SireSec::PrivateKey > );
        PrivateKey_exposer.def( "__repr__", &__str__< ::SireSec::PrivateKey > );
    }

}
