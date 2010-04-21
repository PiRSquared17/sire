#include <Python.h>
#include <boost/python.hpp>

#include "Base/convertproperty.hpp"
#include "SireVol_properties.h"

#include "SireError/errors.h"
#include "SireMaths/rotate.h"
#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "SireUnits/units.h"
#include "grid.h"
#include "grid.h"
#include "SireError/errors.h"
#include "SireMaths/rangenerator.h"
#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "cartesian.h"
#include "space.h"
#include <QMutex>
#include "space.h"
void register_SireVol_properties()
{
    register_property_container< SireVol::GridPtr, SireVol::Grid >();
    register_property_container< SireVol::SpacePtr, SireVol::Space >();
}
