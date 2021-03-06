//WARNING - AUTOGENERATED FILE - CONTENTS WILL BE OVERWRITTEN!
#include <Python.h>

#include "SireMaths_registrars.h"

#include "vectorproperty.h"
#include "freeenergyaverage.h"
#include "quaternion.h"
#include "histogram.h"
#include "n4matrix.h"
#include "nvector.h"
#include "sphere.h"
#include "accumulator.h"
#include "nmatrix.h"
#include "line.h"
#include "vector.h"
#include "matrix.h"
#include "distvector.h"
#include "triangle.h"
#include "align.h"
#include "trigmatrix.h"
#include "torsion.h"
#include "rangenerator.h"
#include "axisset.h"
#include "plane.h"
#include "complex.h"

#include "Helpers/objectregistry.hpp"

void register_SireMaths_objects()
{

    ObjectRegistry::registerConverterFor< SireMaths::VectorProperty >();
    ObjectRegistry::registerConverterFor< SireMaths::VectorArrayProperty >();
    ObjectRegistry::registerConverterFor< SireMaths::FreeEnergyAverage >();
    ObjectRegistry::registerConverterFor< SireMaths::BennettsFreeEnergyAverage >();
    ObjectRegistry::registerConverterFor< SireMaths::Quaternion >();
    ObjectRegistry::registerConverterFor< SireMaths::Histogram >();
    ObjectRegistry::registerConverterFor< SireMaths::N4Matrix >();
    ObjectRegistry::registerConverterFor< SireMaths::NVector >();
    ObjectRegistry::registerConverterFor< SireMaths::Sphere >();
    ObjectRegistry::registerConverterFor< SireMaths::NullAccumulator >();
    ObjectRegistry::registerConverterFor< SireMaths::Average >();
    ObjectRegistry::registerConverterFor< SireMaths::AverageAndStddev >();
    ObjectRegistry::registerConverterFor< SireMaths::ExpAverage >();
    ObjectRegistry::registerConverterFor< SireMaths::Median >();
    ObjectRegistry::registerConverterFor< SireMaths::RecordValues >();
    ObjectRegistry::registerConverterFor< SireMaths::NMatrix >();
    ObjectRegistry::registerConverterFor< SireMaths::Line >();
    ObjectRegistry::registerConverterFor< SireMaths::Vector >();
    ObjectRegistry::registerConverterFor< SireMaths::Matrix >();
    ObjectRegistry::registerConverterFor< SireMaths::DistVector >();
    ObjectRegistry::registerConverterFor< SireMaths::Triangle >();
    ObjectRegistry::registerConverterFor< SireMaths::Transform >();
    ObjectRegistry::registerConverterFor< SireMaths::TrigMatrix >();
    ObjectRegistry::registerConverterFor< SireMaths::Torsion >();
    ObjectRegistry::registerConverterFor< SireMaths::RanGenerator >();
    ObjectRegistry::registerConverterFor< SireMaths::AxisSet >();
    ObjectRegistry::registerConverterFor< SireMaths::Plane >();
    ObjectRegistry::registerConverterFor< SireMaths::Complex >();

}

