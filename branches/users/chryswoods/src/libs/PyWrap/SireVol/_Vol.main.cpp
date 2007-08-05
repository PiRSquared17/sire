// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License


#include "boost/python.hpp"

#include "sirevol_headers.h"

#include "SireMaths/vector.h"

#include "SireMaths/quaternion.h"

#include "SireMaths/matrix.h"

#include "AABox.pypp.hpp"

#include "Cartesian.pypp.hpp"

#include "CoordGroup.pypp.hpp"

#include "CoordGroupBase.pypp.hpp"

#include "CoordGroupEditor.pypp.hpp"

#include "MapFromCartesianFunction.pypp.hpp"

#include "MapFromSelfFunction.pypp.hpp"

#include "MappingFunctionBase.pypp.hpp"

#include "PeriodicBox.pypp.hpp"

#include "SpaceBase.pypp.hpp"

#include "sirevol_containers.h"

namespace bp = boost::python;

BOOST_PYTHON_MODULE(_Vol){
    register_SireVol_containers();

    register_AABox_class();

    register_SpaceBase_class();

    register_Cartesian_class();

    register_CoordGroupBase_class();

    register_CoordGroup_class();

    register_CoordGroupEditor_class();

    register_MappingFunctionBase_class();

    register_MapFromCartesianFunction_class();

    register_MapFromSelfFunction_class();

    bp::implicitly_convertible< QVector<SireMaths::Vector>, SireVol::CoordGroup >();

    bp::implicitly_convertible< const SireVol::MappingFunctionBase&, SireVol::MappingFunction >();

    register_PeriodicBox_class();
}

