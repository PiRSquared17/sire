// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License


#include "boost/python.hpp"

#include "Cube.pypp.hpp"

#include "IOBase.pypp.hpp"

#include "IOParametersBase.pypp.hpp"

#include "NullIO.pypp.hpp"

#include "PDB.pypp.hpp"

#include "PDBParameters.pypp.hpp"

#include "ProtoMS.pypp.hpp"

#include "ProtoMSParameters.pypp.hpp"

#include "Tinker.pypp.hpp"

#include "TinkerParameters.pypp.hpp"

#include "TrajectoryMonitor.pypp.hpp"

namespace bp = boost::python;

#include "SireIO_registrars.h"

#include "SireIO_properties.h"

BOOST_PYTHON_MODULE(_IO){
    register_SireIO_objects();

    register_Cube_class();

    register_IOBase_class();

    register_IOParametersBase_class();

    register_NullIO_class();

    register_PDB_class();

    register_PDBParameters_class();

    register_ProtoMS_class();

    register_ProtoMSParameters_class();

    register_Tinker_class();

    register_TinkerParameters_class();

    register_SireIO_properties();

    register_TrajectoryMonitor_class();
}

