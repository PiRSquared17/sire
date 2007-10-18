// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License


#include "boost/python.hpp"

#include "siresystem_headers.h"

#include "SireMol/molecule.h"

#include "SireMol/molecules.h"

#include "SireMol/residue.h"

#include "SireMol/newatom.h"

#include "SireFF/ffgroupid.h"

#include "CheckPoint.pypp.hpp"

#include "LocalSimSystem.pypp.hpp"

#include "MonitorEnergy.pypp.hpp"

#include "MoveBase.pypp.hpp"

#include "MovesBase.pypp.hpp"

#include "QuerySystem.pypp.hpp"

#include "RDF.pypp.hpp"

#include "RDFMonitor.pypp.hpp"

#include "SameMoves.pypp.hpp"

#include "SimSystem.pypp.hpp"

#include "System.pypp.hpp"

#include "SystemData.pypp.hpp"

#include "SystemID.pypp.hpp"

#include "SystemMonitorBase.pypp.hpp"

#include "SystemMonitors.pypp.hpp"

#include "siresystem_containers.h"

namespace bp = boost::python;

BOOST_PYTHON_MODULE(_System){
    register_SireSystem_containers();

    register_CheckPoint_class();

    register_QuerySystem_class();

    register_SimSystem_class();

    register_LocalSimSystem_class();

    register_SystemMonitorBase_class();

    register_MonitorEnergy_class();

    register_MoveBase_class();

    register_MovesBase_class();

    register_RDF_class();

    register_RDFMonitor_class();

    register_SameMoves_class();

    register_SystemMonitors_class();

    register_System_class();

    register_SystemData_class();

    bp::implicitly_convertible< const SireSystem::MoveBase&, SireSystem::Move >();

    bp::implicitly_convertible< SireSystem::Move, SireSystem::Moves >();

    bp::implicitly_convertible< const SireSystem::MoveBase&, SireSystem::Moves >();

    bp::implicitly_convertible< const SireSystem::MovesBase&, SireSystem::Moves >();

    register_SystemID_class();
}

