//WARNING - AUTOGENERATED FILE - CONTENTS WILL BE OVERWRITTEN!
#include <Python.h>

#include "SireBase_registrars.h"

#include "property.h"
#include "properties.h"
#include "propertymap.h"
#include "majorminorversion.h"
#include "stringmangler.h"
#include "linktoproperty.h"

#include "Helpers/objectregistry.hpp"

void register_SireBase_objects()
{

    ObjectRegistry::registerConverterFor< SireBase::NullProperty >();
    ObjectRegistry::registerConverterFor< SireBase::VariantProperty >();
    ObjectRegistry::registerConverterFor< SireBase::Properties >();
    ObjectRegistry::registerConverterFor< SireBase::PropertyName >();
    ObjectRegistry::registerConverterFor< SireBase::PropertyMap >();
    ObjectRegistry::registerConverterFor< SireBase::MajorMinorVersion >();
    ObjectRegistry::registerConverterFor< SireBase::Version >();
    ObjectRegistry::registerConverterFor< SireBase::NoMangling >();
    ObjectRegistry::registerConverterFor< SireBase::TrimString >();
    ObjectRegistry::registerConverterFor< SireBase::UpperCaseString >();
    ObjectRegistry::registerConverterFor< SireBase::LowerCaseString >();
    ObjectRegistry::registerConverterFor< SireBase::LinkToProperty >();

}

