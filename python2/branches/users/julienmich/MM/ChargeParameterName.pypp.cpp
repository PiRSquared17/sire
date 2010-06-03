// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "ChargeParameterName.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/countflops.h"

#include "SireBase/errors.h"

#include "SireBase/sparsematrix.hpp"

#include "SireError/errors.h"

#include "SireFF/errors.h"

#include "SireMaths/maths.h"

#include "SireMol/atomcoords.h"

#include "SireMol/mover.hpp"

#include "SireStream/datastream.h"

#include "SireUnits/units.h"

#include "SireVol/cartesian.h"

#include "coulombpotential.h"

#include "switchingfunction.h"

#include <QDebug>

#include "coulombpotential.h"

SireMM::ChargeParameterName __copy__(const SireMM::ChargeParameterName &other){ return SireMM::ChargeParameterName(other); }

const char* pvt_get_name(const SireMM::ChargeParameterName&){ return "SireMM::ChargeParameterName";}

void register_ChargeParameterName_class(){

    { //::SireMM::ChargeParameterName
        typedef bp::class_< SireMM::ChargeParameterName > ChargeParameterName_exposer_t;
        ChargeParameterName_exposer_t ChargeParameterName_exposer = ChargeParameterName_exposer_t( "ChargeParameterName" );
        bp::scope ChargeParameterName_scope( ChargeParameterName_exposer );
        ChargeParameterName_exposer.def( bp::init< >() );
        { //::SireMM::ChargeParameterName::charge
        
            typedef ::QString const & ( ::SireMM::ChargeParameterName::*charge_function_type )(  ) const;
            charge_function_type charge_function_value( &::SireMM::ChargeParameterName::charge );
            
            ChargeParameterName_exposer.def( 
                "charge"
                , charge_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        ChargeParameterName_exposer.def( "__copy__", &__copy__);
        ChargeParameterName_exposer.def( "__deepcopy__", &__copy__);
        ChargeParameterName_exposer.def( "clone", &__copy__);
        ChargeParameterName_exposer.def( "__str__", &pvt_get_name);
        ChargeParameterName_exposer.def( "__repr__", &pvt_get_name);
    }

}
