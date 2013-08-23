// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "BeadID.pypp.hpp"

namespace bp = boost::python;

#include "SireStream/datastream.h"

#include "beadid.h"

#include "beadidx.h"

#include "beadnum.h"

#include "beadid.h"

#include "Helpers/str.hpp"

void register_BeadID_class(){

    { //::SireMol::BeadID
        typedef bp::class_< SireMol::BeadID, bp::bases< SireID::ID >, boost::noncopyable > BeadID_exposer_t;
        BeadID_exposer_t BeadID_exposer = BeadID_exposer_t( "BeadID", bp::no_init );
        bp::scope BeadID_scope( BeadID_exposer );
        { //::SireMol::BeadID::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::BeadID::typeName );
            
            BeadID_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        BeadID_exposer.staticmethod( "typeName" );
        BeadID_exposer.def( "__str__", &__str__< ::SireMol::BeadID > );
        BeadID_exposer.def( "__repr__", &__str__< ::SireMol::BeadID > );
    }

}