// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License


#include "boost/python.hpp"

#include "Celsius.pypp.hpp"

#include "Fahrenheit.pypp.hpp"

#include "GeneralUnit.pypp.hpp"

#include "GeneralUnitBase.pypp.hpp"

#include "GeneralUnitObject.pypp.hpp"

#include "Temperature.pypp.hpp"

#include "Unit.pypp.hpp"

#include "_Units_free_functions.pypp.hpp"

namespace bp = boost::python;

#include "SireUnits_containers.h"

#include "SireUnits/temperature.h"

#include "SireUnits/generalunit.h"

#include "sireunits_dimensions.h"

#include "_Units_global_variables.pyman.hpp"

BOOST_PYTHON_MODULE(_Units){
    register_SireUnits_containers();

    register_Temperature_class();

    register_Celsius_class();

    register_Unit_class();

    register_GeneralUnitBase_class();

    register_GeneralUnit_class();

    register_Fahrenheit_class();

    register_GeneralUnitObject_class();

    bp::implicitly_convertible< SireUnits::Dimension::detail::Temperature, SireUnits::Dimension::Temperature >();

    bp::implicitly_convertible< SireUnits::Dimension::detail::Temperature, double >();

    register_SireUnits_dimensions();

    register_man_global_variables();

    register_free_functions();
}

