// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License


#include "boost/python.hpp"

#include "Helpers/clone_const_reference.hpp"

#include "AngleComponent.pypp.hpp"

#include "ChargeConstraint.pypp.hpp"

#include "CheckPoint.pypp.hpp"

#include "CloseMols.pypp.hpp"

#include "ComponentConstraint.pypp.hpp"

#include "Constraint.pypp.hpp"

#include "Constraints.pypp.hpp"

#include "DihedralComponent.pypp.hpp"

#include "DistanceComponent.pypp.hpp"

#include "DoubleDistanceComponent.pypp.hpp"

#include "EnergyMonitor.pypp.hpp"

#include "GeometryComponent.pypp.hpp"

#include "IDAndSet_MonitorID_.pypp.hpp"

#include "IDAndSet_SysID_.pypp.hpp"

#include "IDAssigner.pypp.hpp"

#include "IDOrSet_MonitorID_.pypp.hpp"

#include "IDOrSet_SysID_.pypp.hpp"

#include "IdentityConstraint.pypp.hpp"

#include "MoleculeConstraint.pypp.hpp"

#include "MonitorComponent.pypp.hpp"

#include "MonitorComponents.pypp.hpp"

#include "MonitorID.pypp.hpp"

#include "MonitorIdx.pypp.hpp"

#include "MonitorMonitor.pypp.hpp"

#include "MonitorName.pypp.hpp"

#include "MonitorProperty.pypp.hpp"

#include "NullConstraint.pypp.hpp"

#include "NullMonitor.pypp.hpp"

#include "PerturbationConstraint.pypp.hpp"

#include "PolariseCharges.pypp.hpp"

#include "PolariseChargesFF.pypp.hpp"

#include "PropertyConstraint.pypp.hpp"

#include "SpaceWrapper.pypp.hpp"

#include "Specify_MonitorID_.pypp.hpp"

#include "Specify_SysID_.pypp.hpp"

#include "SysID.pypp.hpp"

#include "SysIdx.pypp.hpp"

#include "SysName.pypp.hpp"

#include "System.pypp.hpp"

#include "SystemMonitor.pypp.hpp"

#include "SystemMonitors.pypp.hpp"

#include "TripleDistanceComponent.pypp.hpp"

#include "WindowedComponent.pypp.hpp"

namespace bp = boost::python;

#include "SireSystem_containers.h"

#include "SireSystem_registrars.h"

#include "SireSystem_properties.h"

BOOST_PYTHON_MODULE(_System){
    register_SireSystem_objects();

    register_SireSystem_containers();

    register_MonitorID_class();

    register_IDAndSet_MonitorID__class();

    register_SysID_class();

    register_IDAndSet_SysID__class();

    register_IDOrSet_MonitorID__class();

    register_IDOrSet_SysID__class();

    register_Specify_MonitorID__class();

    register_Specify_SysID__class();

    register_Constraint_class();

    register_GeometryComponent_class();

    register_AngleComponent_class();

    register_MoleculeConstraint_class();

    register_ChargeConstraint_class();

    register_CheckPoint_class();

    register_CloseMols_class();

    register_ComponentConstraint_class();

    register_Constraints_class();

    register_DihedralComponent_class();

    register_DistanceComponent_class();

    register_DoubleDistanceComponent_class();

    register_SystemMonitor_class();

    register_EnergyMonitor_class();

    register_IDAssigner_class();

    register_IdentityConstraint_class();

    register_MonitorComponent_class();

    register_MonitorComponents_class();

    register_MonitorIdx_class();

    register_MonitorMonitor_class();

    register_MonitorName_class();

    register_MonitorProperty_class();

    register_NullConstraint_class();

    register_NullMonitor_class();

    register_PerturbationConstraint_class();

    register_PolariseCharges_class();

    register_PolariseChargesFF_class();

    register_PropertyConstraint_class();

    register_SpaceWrapper_class();

    register_SysIdx_class();

    register_SysName_class();

    register_System_class();

    register_SystemMonitors_class();

    register_TripleDistanceComponent_class();

    register_SireSystem_properties();

    register_WindowedComponent_class();
}

