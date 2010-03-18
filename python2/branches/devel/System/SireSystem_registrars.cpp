//WARNING - AUTOGENERATED FILE - CONTENTS WILL BE OVERWRITTEN!
#include <Python.h>

#include "SireSystem_registrars.h"

#include "sysidentifier.h"
#include "perturbationconstraint.h"
#include "rdfmonitor.h"
#include "distancecomponent.h"
#include "monitormonitor.h"
#include "spacewrapper.h"
#include "monitoridx.h"
#include "constraints.h"
#include "monitorname.h"
#include "sysidx.h"
#include "monitorcomponent.h"
#include "systemmonitor.h"
#include "sysid.h"
#include "constraint.h"
#include "checkpoint.h"
#include "closemols.h"
#include "delta.h"
#include "monitoridentifier.h"
#include "systemmonitors.h"
#include "system.h"
#include "monitorcomponents.h"
#include "identityconstraint.h"
#include "sysname.h"
#include "monitorid.h"

#include "Helpers/objectregistry.hpp"

void register_SireSystem_objects()
{

    ObjectRegistry::registerConverterFor< SireSystem::SysIdentifier >();
    ObjectRegistry::registerConverterFor< SireSystem::PerturbationConstraint >();
    ObjectRegistry::registerConverterFor< SireSystem::RDF >();
    ObjectRegistry::registerConverterFor< SireSystem::RDFMonitor >();
    ObjectRegistry::registerConverterFor< SireSystem::DistanceComponent >();
    ObjectRegistry::registerConverterFor< SireSystem::DoubleDistanceComponent >();
    ObjectRegistry::registerConverterFor< SireSystem::TripleDistanceComponent >();
    ObjectRegistry::registerConverterFor< SireSystem::MonitorMonitor >();
    ObjectRegistry::registerConverterFor< SireSystem::SpaceWrapper >();
    ObjectRegistry::registerConverterFor< SireSystem::MonitorIdx >();
    ObjectRegistry::registerConverterFor< SireSystem::Constraints >();
    ObjectRegistry::registerConverterFor< SireSystem::MonitorName >();
    ObjectRegistry::registerConverterFor< SireSystem::SysIdx >();
    ObjectRegistry::registerConverterFor< SireSystem::MonitorComponent >();
    ObjectRegistry::registerConverterFor< SireSystem::NullMonitor >();
    ObjectRegistry::registerConverterFor< SireID::Specify<SireSystem::SysID> >();
    ObjectRegistry::registerConverterFor< SireID::IDAndSet<SireSystem::SysID> >();
    ObjectRegistry::registerConverterFor< SireID::IDOrSet<SireSystem::SysID> >();
    ObjectRegistry::registerConverterFor< SireSystem::NullConstraint >();
    ObjectRegistry::registerConverterFor< SireSystem::PropertyConstraint >();
    ObjectRegistry::registerConverterFor< SireSystem::ComponentConstraint >();
    ObjectRegistry::registerConverterFor< SireSystem::WindowedComponent >();
    ObjectRegistry::registerConverterFor< SireSystem::CheckPoint >();
    ObjectRegistry::registerConverterFor< SireSystem::CloseMols >();
    ObjectRegistry::registerConverterFor< SireSystem::Delta >();
    ObjectRegistry::registerConverterFor< SireSystem::MonitorIdentifier >();
    ObjectRegistry::registerConverterFor< SireSystem::SystemMonitors >();
    ObjectRegistry::registerConverterFor< SireSystem::System >();
    ObjectRegistry::registerConverterFor< SireSystem::MonitorComponents >();
    ObjectRegistry::registerConverterFor< SireSystem::IdentityConstraint >();
    ObjectRegistry::registerConverterFor< SireSystem::SysName >();
    ObjectRegistry::registerConverterFor< SireID::Specify<SireSystem::MonitorID> >();
    ObjectRegistry::registerConverterFor< SireID::IDAndSet<SireSystem::MonitorID> >();
    ObjectRegistry::registerConverterFor< SireID::IDOrSet<SireSystem::MonitorID> >();

}

