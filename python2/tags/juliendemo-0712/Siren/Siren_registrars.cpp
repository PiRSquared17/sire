//WARNING - AUTOGENERATED FILE - CONTENTS WILL BE OVERWRITTEN!
#include <Python.h>

#include "Siren_registrars.h"

#include "datastream.h"
#include "tester.h"
#include "logger.h"
#include "primitive.h"
#include "objptr.hpp"
#include "md5sum.h"
#include "object.h"
#include "class.h"
#include "xmlstream.h"

#include "Helpers/objectregistry.hpp"

void register_Siren_objects()
{

    ObjectRegistry::registerConverterFor< Siren::DataStream >();
    ObjectRegistry::registerConverterFor< Siren::Tester >();
    ObjectRegistry::registerConverterFor< Siren::Logger >();
    ObjectRegistry::registerConverterFor< Siren::String >();
    ObjectRegistry::registerConverterFor< Siren::Number >();
    ObjectRegistry::registerConverterFor< Siren::ObjectPtr >();
    ObjectRegistry::registerConverterFor< Siren::GlobalObjectPtr >();
    ObjectRegistry::registerConverterFor< Siren::MD5Sum >();
    ObjectRegistry::registerConverterFor< Siren::None >();
    ObjectRegistry::registerConverterFor< Siren::Class >();
    ObjectRegistry::registerConverterFor< Siren::XMLStream >();

}
