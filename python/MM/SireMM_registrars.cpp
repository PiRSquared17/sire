//WARNING - AUTOGENERATED FILE - CONTENTS WILL BE OVERWRITTEN!
#include <Python.h>

#include "SireMM_registrars.h"

#include "cljrffunction.h"
#include "dihedralrestraint.h"
#include "ljparameter.h"
#include "softcljcomponent.h"
#include "cljshiftfunction.h"
#include "ljpair.h"
#include "restraintcomponent.h"
#include "intragroupff.h"
#include "gridff.h"
#include "anglerestraint.h"
#include "intersoftcljff.h"
#include "cljboxes.h"
#include "interljff.h"
#include "cljnbpairs.h"
#include "cljprobe.h"
#include "intrasoftcljff.h"
#include "switchingfunction.h"
#include "cljfunction.h"
#include "intergroupff.h"
#include "multicljcomponent.h"
#include "cljcomponent.h"
#include "clj14group.h"
#include "atomljs.h"
#include "distancerestraint.h"
#include "fouratomfunctions.h"
#include "cljcalculator.h"
#include "cljgrid.h"
#include "intracljff.h"
#include "cljdelta.h"
#include "cljextractor.h"
#include "internalparameters.h"
#include "internalcomponent.h"
#include "threeatomfunctions.h"
#include "restraintff.h"
#include "interff.h"
#include "cljworkspace.h"
#include "intercoulombff.h"
#include "intraff.h"
#include "cljatoms.h"
#include "intraljff.h"
#include "internalperturbation.h"
#include "intercljff.h"
#include "ljperturbation.h"
#include "twoatomfunctions.h"
#include "intracoulombff.h"
#include "internalff.h"
#include "cljparam.h"
#include "restraint.h"
#include "cljgroup.h"
#include "gridff2.h"

#include "Helpers/objectregistry.hpp"

void register_SireMM_objects()
{

    ObjectRegistry::registerConverterFor< SireMM::CLJRFFunction >();
    ObjectRegistry::registerConverterFor< SireMM::CLJIntraRFFunction >();
    ObjectRegistry::registerConverterFor< SireMM::CLJSoftRFFunction >();
    ObjectRegistry::registerConverterFor< SireMM::CLJSoftIntraRFFunction >();
    ObjectRegistry::registerConverterFor< SireMM::DihedralRestraint >();
    ObjectRegistry::registerConverterFor< SireMM::LJParameter >();
    ObjectRegistry::registerConverterFor< SireMM::SoftCLJComponent >();
    ObjectRegistry::registerConverterFor< SireMM::CLJShiftFunction >();
    ObjectRegistry::registerConverterFor< SireMM::CLJIntraShiftFunction >();
    ObjectRegistry::registerConverterFor< SireMM::CLJSoftShiftFunction >();
    ObjectRegistry::registerConverterFor< SireMM::CLJSoftIntraShiftFunction >();
    ObjectRegistry::registerConverterFor< SireMM::LJPair >();
    ObjectRegistry::registerConverterFor< SireMM::RestraintComponent >();
    ObjectRegistry::registerConverterFor< SireMM::IntraGroupFF >();
    ObjectRegistry::registerConverterFor< SireMM::GridFF >();
    ObjectRegistry::registerConverterFor< SireMM::AngleRestraint >();
    ObjectRegistry::registerConverterFor< SireMM::InterSoftCLJFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::InterSoftCLJFF >();
    ObjectRegistry::registerConverterFor< SireMM::InterGroupSoftCLJFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::InterGroupSoftCLJFF >();
    ObjectRegistry::registerConverterFor< SireMM::CLJBox >();
    ObjectRegistry::registerConverterFor< SireMM::CLJBoxIndex >();
    ObjectRegistry::registerConverterFor< SireMM::CLJBoxes >();
    ObjectRegistry::registerConverterFor< SireMM::CLJBoxDistance >();
    ObjectRegistry::registerConverterFor< SireMM::InterLJFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::InterLJFF >();
    ObjectRegistry::registerConverterFor< SireMM::InterGroupLJFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::InterGroupLJFF >();
    ObjectRegistry::registerConverterFor< SireMM::CoulombScaleFactor >();
    ObjectRegistry::registerConverterFor< SireMM::LJScaleFactor >();
    ObjectRegistry::registerConverterFor< SireMM::CLJScaleFactor >();
    ObjectRegistry::registerConverterFor< SireMM::CLJNBPairs >();
    ObjectRegistry::registerConverterFor< SireMM::CoulombNBPairs >();
    ObjectRegistry::registerConverterFor< SireMM::LJNBPairs >();
    ObjectRegistry::registerConverterFor< SireMM::CoulombProbe >();
    ObjectRegistry::registerConverterFor< SireMM::LJProbe >();
    ObjectRegistry::registerConverterFor< SireMM::CLJProbe >();
    ObjectRegistry::registerConverterFor< SireMM::IntraSoftCLJFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::IntraSoftCLJFF >();
    ObjectRegistry::registerConverterFor< SireMM::IntraGroupSoftCLJFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::IntraGroupSoftCLJFF >();
    ObjectRegistry::registerConverterFor< SireMM::NoCutoff >();
    ObjectRegistry::registerConverterFor< SireMM::HarmonicSwitchingFunction >();
    ObjectRegistry::registerConverterFor< SireMM::CHARMMSwitchingFunction >();
    ObjectRegistry::registerConverterFor< SireMM::NullCLJFunction >();
    ObjectRegistry::registerConverterFor< SireMM::InterGroupFF >();
    ObjectRegistry::registerConverterFor< SireMM::MultiCLJComponent >();
    ObjectRegistry::registerConverterFor< SireMM::CoulombComponent >();
    ObjectRegistry::registerConverterFor< SireMM::LJComponent >();
    ObjectRegistry::registerConverterFor< SireMM::CLJComponent >();
    ObjectRegistry::registerConverterFor< SireMM::CLJ14Group >();
    ObjectRegistry::registerConverterFor< SireMM::AtomLJs >();
    ObjectRegistry::registerConverterFor< SireMM::DistanceRestraint >();
    ObjectRegistry::registerConverterFor< SireMM::DoubleDistanceRestraint >();
    ObjectRegistry::registerConverterFor< SireMM::TripleDistanceRestraint >();
    ObjectRegistry::registerConverterFor< SireMM::HarmonicDistanceForceConstant >();
    ObjectRegistry::registerConverterFor< SireMM::FourAtomFunctions >();
    ObjectRegistry::registerConverterFor< SireMM::CLJCalculator >();
    ObjectRegistry::registerConverterFor< SireMM::CLJGrid >();
    ObjectRegistry::registerConverterFor< SireMM::IntraCLJFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::IntraCLJFF >();
    ObjectRegistry::registerConverterFor< SireMM::IntraGroupCLJFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::IntraGroupCLJFF >();
    ObjectRegistry::registerConverterFor< SireMM::CLJDelta >();
    ObjectRegistry::registerConverterFor< SireMM::CLJExtractor >();
    ObjectRegistry::registerConverterFor< SireMM::InternalParameters >();
    ObjectRegistry::registerConverterFor< SireMM::InternalParameters3D >();
    ObjectRegistry::registerConverterFor< SireMM::BondComponent >();
    ObjectRegistry::registerConverterFor< SireMM::AngleComponent >();
    ObjectRegistry::registerConverterFor< SireMM::DihedralComponent >();
    ObjectRegistry::registerConverterFor< SireMM::ImproperComponent >();
    ObjectRegistry::registerConverterFor< SireMM::UreyBradleyComponent >();
    ObjectRegistry::registerConverterFor< SireMM::StretchStretchComponent >();
    ObjectRegistry::registerConverterFor< SireMM::StretchBendComponent >();
    ObjectRegistry::registerConverterFor< SireMM::BendBendComponent >();
    ObjectRegistry::registerConverterFor< SireMM::StretchBendTorsionComponent >();
    ObjectRegistry::registerConverterFor< SireMM::Intra14CoulombComponent >();
    ObjectRegistry::registerConverterFor< SireMM::Intra14LJComponent >();
    ObjectRegistry::registerConverterFor< SireMM::Intra14Component >();
    ObjectRegistry::registerConverterFor< SireMM::InternalComponent >();
    ObjectRegistry::registerConverterFor< SireMM::ThreeAtomFunctions >();
    ObjectRegistry::registerConverterFor< SireMM::RestraintFF >();
    ObjectRegistry::registerConverterFor< SireMM::InterFF >();
    ObjectRegistry::registerConverterFor< SireMM::CLJWorkspace >();
    ObjectRegistry::registerConverterFor< SireMM::InterCoulombFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::InterCoulombFF >();
    ObjectRegistry::registerConverterFor< SireMM::InterGroupCoulombFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::InterGroupCoulombFF >();
    ObjectRegistry::registerConverterFor< SireMM::IntraFF >();
    ObjectRegistry::registerConverterFor< SireMM::CLJAtom >();
    ObjectRegistry::registerConverterFor< SireMM::CLJAtoms >();
    ObjectRegistry::registerConverterFor< SireMM::IntraLJFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::IntraLJFF >();
    ObjectRegistry::registerConverterFor< SireMM::IntraGroupLJFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::IntraGroupLJFF >();
    ObjectRegistry::registerConverterFor< SireMM::TwoAtomPerturbation >();
    ObjectRegistry::registerConverterFor< SireMM::ThreeAtomPerturbation >();
    ObjectRegistry::registerConverterFor< SireMM::FourAtomPerturbation >();
    ObjectRegistry::registerConverterFor< SireMM::InterCLJFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::InterCLJFF >();
    ObjectRegistry::registerConverterFor< SireMM::InterGroupCLJFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::InterGroupCLJFF >();
    ObjectRegistry::registerConverterFor< SireMM::LJPerturbation >();
    ObjectRegistry::registerConverterFor< SireMM::TwoAtomFunctions >();
    ObjectRegistry::registerConverterFor< SireMM::IntraCoulombFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::IntraCoulombFF >();
    ObjectRegistry::registerConverterFor< SireMM::IntraGroupCoulombFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::IntraGroupCoulombFF >();
    ObjectRegistry::registerConverterFor< SireMM::InternalFF >();
    ObjectRegistry::registerConverterFor< SireMM::CLJParams >();
    ObjectRegistry::registerConverterFor< SireMM::CLJParamsArray >();
    ObjectRegistry::registerConverterFor< SireMM::NullRestraint >();
    ObjectRegistry::registerConverterFor< SireMM::CLJGroup >();
    ObjectRegistry::registerConverterFor< SireMM::GridFF2 >();

}

