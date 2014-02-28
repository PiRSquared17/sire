// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "NullCLJFunction.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/errors.h"

#include "SireBase/lengthproperty.h"

#include "SireBase/numberproperty.h"

#include "SireBase/properties.h"

#include "SireBase/stringproperty.h"

#include "SireError/errors.h"

#include "SireMaths/multidouble.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireVol/cartesian.h"

#include "SireVol/periodicbox.h"

#include "cljboxes.h"

#include "cljdelta.h"

#include "cljfunction.h"

#include "gridinfo.h"

#include "switchingfunction.h"

#include "tbb/blocked_range.h"

#include "tbb/parallel_for.h"

#include "tostring.h"

#include <QElapsedTimer>

#include "cljfunction.h"

SireMM::NullCLJFunction __copy__(const SireMM::NullCLJFunction &other){ return SireMM::NullCLJFunction(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_NullCLJFunction_class(){

    { //::SireMM::NullCLJFunction
        typedef bp::class_< SireMM::NullCLJFunction, bp::bases< SireMM::CLJFunction, SireBase::Property > > NullCLJFunction_exposer_t;
        NullCLJFunction_exposer_t NullCLJFunction_exposer = NullCLJFunction_exposer_t( "NullCLJFunction", bp::init< >() );
        bp::scope NullCLJFunction_scope( NullCLJFunction_exposer );
        NullCLJFunction_exposer.def( bp::init< SireMM::NullCLJFunction const & >(( bp::arg("other") )) );
        NullCLJFunction_exposer.def( bp::self != bp::self );
        { //::SireMM::NullCLJFunction::operator=
        
            typedef ::SireMM::NullCLJFunction & ( ::SireMM::NullCLJFunction::*assign_function_type )( ::SireMM::NullCLJFunction const & ) ;
            assign_function_type assign_function_value( &::SireMM::NullCLJFunction::operator= );
            
            NullCLJFunction_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        NullCLJFunction_exposer.def( bp::self == bp::self );
        { //::SireMM::NullCLJFunction::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMM::NullCLJFunction::typeName );
            
            NullCLJFunction_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireMM::NullCLJFunction::what
        
            typedef char const * ( ::SireMM::NullCLJFunction::*what_function_type )(  ) const;
            what_function_type what_function_value( &::SireMM::NullCLJFunction::what );
            
            NullCLJFunction_exposer.def( 
                "what"
                , what_function_value );
        
        }
        NullCLJFunction_exposer.staticmethod( "typeName" );
        NullCLJFunction_exposer.def( "__copy__", &__copy__);
        NullCLJFunction_exposer.def( "__deepcopy__", &__copy__);
        NullCLJFunction_exposer.def( "clone", &__copy__);
        NullCLJFunction_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::NullCLJFunction >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        NullCLJFunction_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::NullCLJFunction >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        NullCLJFunction_exposer.def( "__str__", &__str__< ::SireMM::NullCLJFunction > );
        NullCLJFunction_exposer.def( "__repr__", &__str__< ::SireMM::NullCLJFunction > );
    }

}
