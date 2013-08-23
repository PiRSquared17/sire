// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License


#include "boost/python.hpp"

#include "Helpers/clone_const_reference.hpp"

#include "AM1BCC.pypp.hpp"

#include "ChargeElementParameterNames3D.pypp.hpp"

#include "ElementParameterName.pypp.hpp"

#include "ElementParameterName3D.pypp.hpp"

#include "GTO.pypp.hpp"

#include "GTOPair.pypp.hpp"

#include "HF.pypp.hpp"

#include "Molpro.pypp.hpp"

#include "Mopac.pypp.hpp"

#include "NullQM.pypp.hpp"

#include "NullQMChargeCalculator.pypp.hpp"

#include "Orbital.pypp.hpp"

#include "OrbitalShell.pypp.hpp"

#include "PP_GTO.pypp.hpp"

#include "PP_GTOs.pypp.hpp"

#include "PS_GTO.pypp.hpp"

#include "PS_GTOs.pypp.hpp"

#include "P_GTO.pypp.hpp"

#include "PointCharge.pypp.hpp"

#include "PointDipole.pypp.hpp"

#include "QMChargeCalculator.pypp.hpp"

#include "QMChargeConstraint.pypp.hpp"

#include "QMComponent.pypp.hpp"

#include "QMFF.pypp.hpp"

#include "QMMMFF.pypp.hpp"

#include "QMProgram.pypp.hpp"

#include "SS_GTO.pypp.hpp"

#include "SS_GTOs.pypp.hpp"

#include "S_GTO.pypp.hpp"

#include "ShellPair.pypp.hpp"

#include "_Squire_free_functions.pypp.hpp"

namespace bp = boost::python;

#include "Squire_containers.h"

#include "Squire_registrars.h"

#include "Squire_properties.h"

BOOST_PYTHON_MODULE(_Squire){
    register_Squire_objects();

    register_Squire_containers();

    register_QMChargeCalculator_class();

    register_AM1BCC_class();

    register_ElementParameterName_class();

    register_ChargeElementParameterNames3D_class();

    register_ElementParameterName3D_class();

    register_Orbital_class();

    register_OrbitalShell_class();

    register_GTO_class();

    register_ShellPair_class();

    register_GTOPair_class();

    register_HF_class();

    register_QMProgram_class();

    register_Molpro_class();

    register_Mopac_class();

    register_NullQM_class();

    register_NullQMChargeCalculator_class();

    register_PP_GTO_class();

    register_PP_GTOs_class();

    register_PS_GTO_class();

    register_PS_GTOs_class();

    register_P_GTO_class();

    register_PointCharge_class();

    register_PointDipole_class();

    register_QMChargeConstraint_class();

    register_QMComponent_class();

    register_QMFF_class();

    register_QMMMFF_class();

    register_SS_GTO_class();

    register_SS_GTOs_class();

    register_S_GTO_class();

    register_Squire_properties();

    register_free_functions();
}

