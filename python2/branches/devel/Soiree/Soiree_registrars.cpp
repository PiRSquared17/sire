//WARNING - AUTOGENERATED FILE - CONTENTS WILL BE OVERWRITTEN!
#include <Python.h>

#include "Soiree_registrars.h"

#include "fep.h"
#include "ti.h"

#include "Helpers/objectregistry.hpp"

void register_Soiree_objects()
{

    ObjectRegistry::registerConverterFor< Soiree::FEP >();
    ObjectRegistry::registerConverterFor< Soiree::FEPDeltas >();
    ObjectRegistry::registerConverterFor< Soiree::DataPoint >();
    ObjectRegistry::registerConverterFor< Soiree::PMF >();
    ObjectRegistry::registerConverterFor< Soiree::Gradients >();
    ObjectRegistry::registerConverterFor< Soiree::TI >();
    ObjectRegistry::registerConverterFor< Soiree::TIPMF >();

}

