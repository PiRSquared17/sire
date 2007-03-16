// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "boost/python.hpp"
#include "siremol_headers.h"
#include "SireMaths/angle.h"
#include "SireMaths/quaternion.h"
#include "SireMaths/matrix.h"
#include "SireMaths/triangle.h"
#include "SireMaths/line.h"
#include "SireMaths/torsion.h"
#include "_EditMol__enums.pypp.hpp"
#include "_EditMol__unnamed_enums.pypp.hpp"
#include "_EditMol__classes.pypp.hpp"
#include "_EditMol__memfuns0.pypp.hpp"
#include "_EditMol__memfuns1.pypp.hpp"
#include "_EditMol__memfuns2.pypp.hpp"
#include "_EditMol__memfuns3.pypp.hpp"
#include "_EditMol__memfuns4.pypp.hpp"
#include "_EditMol__memfuns5.pypp.hpp"
#include "_EditMol__memfuns6.pypp.hpp"
#include "_EditMol__memfuns7.pypp.hpp"
#include "_EditMol__memfuns8.pypp.hpp"
#include "_EditMol__memfuns9.pypp.hpp"
#include "_EditMol__memfuns10.pypp.hpp"
#include "_EditMol__memfuns11.pypp.hpp"
#include "_EditMol__memfuns12.pypp.hpp"
#include "_EditMol__memfuns13.pypp.hpp"
#include "_EditMol__memfuns14.pypp.hpp"
#include "_EditMol__memfuns_virtual.pypp.hpp"
#include "_EditMol__memfuns_pvirtual.pypp.hpp"
#include "_EditMol__protected_memfuns.pypp.hpp"

namespace bp = boost::python;

SireMol::EditMol __copy__(const SireMol::EditMol &other){ return SireMol::EditMol(other); }

#include "SireQt/qdatastream.hpp"

const char* pvt_get_name(const SireMol::EditMol&){ return "SireMol::EditMol";}

void register_EditMol_class(){

{ //::SireMol::EditMol
    typedef bp::class_< SireMol::EditMol > EditMol_exposer_t;
    EditMol_exposer_t EditMol_exposer = EditMol_exposer_t( "EditMol" );
    bp::scope EditMol_scope( EditMol_exposer );
    EditMol_exposer.def( bp::init< >() );
    EditMol_exposer.def( bp::init< QString const & >(( bp::arg("name") )) );
    EditMol_exposer.def( bp::init< QString const &, SireMol::CuttingFunction const & >(( bp::arg("name"), bp::arg("cuttingfunc") )) );
    EditMol_exposer.def( bp::init< SireMol::CuttingFunction const & >(( bp::arg("cuttingfunc") )) );
    EditMol_exposer.def( bp::init< SireMol::EditRes const & >(( bp::arg("other") )) );
    register_EditMol_memfuns0(EditMol_exposer);
    register_EditMol_memfuns1(EditMol_exposer);
    register_EditMol_memfuns2(EditMol_exposer);
    register_EditMol_memfuns3(EditMol_exposer);
    register_EditMol_memfuns4(EditMol_exposer);
    register_EditMol_memfuns5(EditMol_exposer);
    register_EditMol_memfuns6(EditMol_exposer);
    register_EditMol_memfuns7(EditMol_exposer);
    register_EditMol_memfuns8(EditMol_exposer);
    EditMol_exposer.def( bp::self != bp::self );
    EditMol_exposer.def( bp::self == bp::self );
    register_EditMol_memfuns9(EditMol_exposer);
    register_EditMol_memfuns10(EditMol_exposer);
    register_EditMol_memfuns11(EditMol_exposer);
    register_EditMol_memfuns12(EditMol_exposer);
    register_EditMol_memfuns13(EditMol_exposer);
    register_EditMol_memfuns14(EditMol_exposer);
    EditMol_exposer.def( "__copy__", &__copy__);
    EditMol_exposer.def( "__rlshift__", &SireQt::__rlshift__QDataStream< ::SireMol::EditMol >,
                        bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
    EditMol_exposer.def( "__rrshift__", &SireQt::__rrshift__QDataStream< ::SireMol::EditMol >,
                        bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
    EditMol_exposer.def( "__str__", &pvt_get_name);
}

}
