// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Tinker.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireIO/errors.h"

#include "SireMaths/vector.h"

#include "SireMol/element.h"

#include "SireMol/moleculegroup.h"

#include "SireMol/molecules.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "tinker.h"

#include <QFile>

#include <QRegExp>

#include <QString>

#include <QStringList>

#include <QTextStream>

#include "tinker.h"

SireIO::Tinker __copy__(const SireIO::Tinker &other){ return SireIO::Tinker(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_Tinker_class(){

    { //::SireIO::Tinker
        typedef bp::class_< SireIO::Tinker, bp::bases< SireIO::IOBase, SireBase::Property > > Tinker_exposer_t;
        Tinker_exposer_t Tinker_exposer = Tinker_exposer_t( "Tinker" );
        bp::scope Tinker_scope( Tinker_exposer );
        Tinker_exposer.def( bp::init< >() );
        Tinker_exposer.def( bp::init< SireIO::Tinker const & >(( bp::arg("other") )) );
        { //::SireIO::Tinker::loadParameters
        
            typedef void ( ::SireIO::Tinker::*loadParameters_function_type )( ::QString const & ) ;
            loadParameters_function_type loadParameters_function_value( &::SireIO::Tinker::loadParameters );
            
            Tinker_exposer.def( 
                "loadParameters"
                , loadParameters_function_value
                , ( bp::arg("prmfile") ) );
        
        }
        Tinker_exposer.def( bp::self != bp::self );
        { //::SireIO::Tinker::operator=
        
            typedef ::SireIO::Tinker & ( ::SireIO::Tinker::*assign_function_type )( ::SireIO::Tinker const & ) ;
            assign_function_type assign_function_value( &::SireIO::Tinker::operator= );
            
            Tinker_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        Tinker_exposer.def( bp::self == bp::self );
        { //::SireIO::Tinker::parameters
        
            typedef ::SireIO::TinkerParameters const & ( *parameters_function_type )(  );
            parameters_function_type parameters_function_value( &::SireIO::Tinker::parameters );
            
            Tinker_exposer.def( 
                "parameters"
                , parameters_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireIO::Tinker::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireIO::Tinker::typeName );
            
            Tinker_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        Tinker_exposer.staticmethod( "parameters" );
        Tinker_exposer.staticmethod( "typeName" );
        Tinker_exposer.def( "__copy__", &__copy__);
        Tinker_exposer.def( "__deepcopy__", &__copy__);
        Tinker_exposer.def( "clone", &__copy__);
        Tinker_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireIO::Tinker >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Tinker_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireIO::Tinker >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Tinker_exposer.def( "__str__", &__str__< ::SireIO::Tinker > );
        Tinker_exposer.def( "__repr__", &__str__< ::SireIO::Tinker > );
    }

}
