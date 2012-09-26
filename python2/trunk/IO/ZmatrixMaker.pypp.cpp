// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "ZmatrixMaker.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireIO/errors.h"

#include "SireMol/bondid.h"

#include "SireMol/connectivity.h"

#include "SireMol/errors.h"

#include "SireMol/molecule.h"

#include "SireMol/moleculedata.h"

#include "SireMol/molecules.h"

#include "SireMol/moleditor.h"

#include "SireMol/reseditor.h"

#include "SireMol/selector.hpp"

#include "SireMove/zmatrix.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireUnits/convert.h"

#include "SireUnits/units.h"

#include "zmatrixmaker.h"

#include <QFile>

#include <QTextStream>

#include "zmatrixmaker.h"

SireIO::ZmatrixMaker __copy__(const SireIO::ZmatrixMaker &other){ return SireIO::ZmatrixMaker(other); }

#include "Qt/qdatastream.hpp"

const char* pvt_get_name(const SireIO::ZmatrixMaker&){ return "SireIO::ZmatrixMaker";}

void register_ZmatrixMaker_class(){

    { //::SireIO::ZmatrixMaker
        typedef bp::class_< SireIO::ZmatrixMaker > ZmatrixMaker_exposer_t;
        ZmatrixMaker_exposer_t ZmatrixMaker_exposer = ZmatrixMaker_exposer_t( "ZmatrixMaker", bp::init< >() );
        bp::scope ZmatrixMaker_scope( ZmatrixMaker_exposer );
        { //::SireIO::ZmatrixMaker::applyTemplates
        
            typedef ::SireMol::Molecule ( ::SireIO::ZmatrixMaker::*applyTemplates_function_type )( ::SireMol::Molecule & ) ;
            applyTemplates_function_type applyTemplates_function_value( &::SireIO::ZmatrixMaker::applyTemplates );
            
            ZmatrixMaker_exposer.def( 
                "applyTemplates"
                , applyTemplates_function_value
                , ( bp::arg("molecule") ) );
        
        }
        { //::SireIO::ZmatrixMaker::loadTemplates
        
            typedef void ( ::SireIO::ZmatrixMaker::*loadTemplates_function_type )( ::QString const & ) ;
            loadTemplates_function_type loadTemplates_function_value( &::SireIO::ZmatrixMaker::loadTemplates );
            
            ZmatrixMaker_exposer.def( 
                "loadTemplates"
                , loadTemplates_function_value
                , ( bp::arg("templatesfile") ) );
        
        }
        { //::SireIO::ZmatrixMaker::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireIO::ZmatrixMaker::typeName );
            
            ZmatrixMaker_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireIO::ZmatrixMaker::what
        
            typedef char const * ( ::SireIO::ZmatrixMaker::*what_function_type )(  ) const;
            what_function_type what_function_value( &::SireIO::ZmatrixMaker::what );
            
            ZmatrixMaker_exposer.def( 
                "what"
                , what_function_value );
        
        }
        ZmatrixMaker_exposer.staticmethod( "typeName" );
        ZmatrixMaker_exposer.def( "__copy__", &__copy__);
        ZmatrixMaker_exposer.def( "__deepcopy__", &__copy__);
        ZmatrixMaker_exposer.def( "clone", &__copy__);
        ZmatrixMaker_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireIO::ZmatrixMaker >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ZmatrixMaker_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireIO::ZmatrixMaker >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ZmatrixMaker_exposer.def( "__str__", &pvt_get_name);
        ZmatrixMaker_exposer.def( "__repr__", &pvt_get_name);
    }

}
