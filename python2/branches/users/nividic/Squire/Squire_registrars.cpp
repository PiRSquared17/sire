//WARNING - AUTOGENERATED FILE - CONTENTS WILL BE OVERWRITTEN!
#include <Python.h>

#include "Squire_registrars.h"

#include "pgto.h"
#include "qmprogram.h"
#include "qmff.h"
#include "am1bcc.h"
#include "qmchargeconstraint.h"
#include "sgto.h"
#include "qmmmff.h"
#include "mopac.h"
#include "qmchargecalculator.h"
#include "qmpotential.h"
#include "pointdipole.h"
#include "molpro.h"
#include "pointcharge.h"

#include "Helpers/objectregistry.hpp"

void register_Squire_objects()
{

    ObjectRegistry::registerConverterFor< Squire::P_GTO >();
    ObjectRegistry::registerConverterFor< Squire::PS_GTO >();
    ObjectRegistry::registerConverterFor< Squire::PP_GTO >();
    ObjectRegistry::registerConverterFor< Squire::NullQM >();
    ObjectRegistry::registerConverterFor< Squire::QMFF >();
    ObjectRegistry::registerConverterFor< Squire::AM1BCC >();
    ObjectRegistry::registerConverterFor< Squire::QMChargeConstraint >();
    ObjectRegistry::registerConverterFor< Squire::S_GTO >();
    ObjectRegistry::registerConverterFor< Squire::SS_GTO >();
    ObjectRegistry::registerConverterFor< Squire::QMMMFF >();
    ObjectRegistry::registerConverterFor< Squire::Mopac >();
    ObjectRegistry::registerConverterFor< Squire::NullQMChargeCalculator >();
    ObjectRegistry::registerConverterFor< Squire::QMComponent >();
    ObjectRegistry::registerConverterFor< Squire::PointDipole >();
    ObjectRegistry::registerConverterFor< Squire::Molpro >();
    ObjectRegistry::registerConverterFor< Squire::PointCharge >();

}

