// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "CovalentBondHunterParameters.pypp.hpp"

namespace bp = boost::python;

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireVol/coordgroup.h"

#include "atom.h"

#include "atomcoords.h"

#include "atomelements.h"

#include "atomselection.h"

#include "bondhunter.h"

#include "connectivity.h"

#include "molecule.h"

#include "moleculedata.h"

#include "moleculeinfodata.h"

#include "moleculeview.h"

#include "mover.hpp"

#include "selector.hpp"

#include <QDebug>

#include <QMutex>

#include "bondhunter.h"

SireMol::CovalentBondHunterParameters __copy__(const SireMol::CovalentBondHunterParameters &other){ return SireMol::CovalentBondHunterParameters(other); }

const char* pvt_get_name(const SireMol::CovalentBondHunterParameters&){ return "SireMol::CovalentBondHunterParameters";}

void register_CovalentBondHunterParameters_class(){

    { //::SireMol::CovalentBondHunterParameters
        typedef bp::class_< SireMol::CovalentBondHunterParameters > CovalentBondHunterParameters_exposer_t;
        CovalentBondHunterParameters_exposer_t CovalentBondHunterParameters_exposer = CovalentBondHunterParameters_exposer_t( "CovalentBondHunterParameters" );
        bp::scope CovalentBondHunterParameters_scope( CovalentBondHunterParameters_exposer );
        CovalentBondHunterParameters_exposer.def( bp::init< >() );
        { //::SireMol::CovalentBondHunterParameters::coordinates
        
            typedef ::SireBase::PropertyName const & ( *coordinates_function_type )(  );
            coordinates_function_type coordinates_function_value( &::SireMol::CovalentBondHunterParameters::coordinates );
            
            CovalentBondHunterParameters_exposer.def( 
                "coordinates"
                , coordinates_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMol::CovalentBondHunterParameters::element
        
            typedef ::SireBase::PropertyName const & ( *element_function_type )(  );
            element_function_type element_function_value( &::SireMol::CovalentBondHunterParameters::element );
            
            CovalentBondHunterParameters_exposer.def( 
                "element"
                , element_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        CovalentBondHunterParameters_exposer.staticmethod( "coordinates" );
        CovalentBondHunterParameters_exposer.staticmethod( "element" );
        CovalentBondHunterParameters_exposer.def( "__copy__", &__copy__);
        CovalentBondHunterParameters_exposer.def( "__deepcopy__", &__copy__);
        CovalentBondHunterParameters_exposer.def( "clone", &__copy__);
        CovalentBondHunterParameters_exposer.def( "__str__", &pvt_get_name);
        CovalentBondHunterParameters_exposer.def( "__repr__", &pvt_get_name);
    }

}
