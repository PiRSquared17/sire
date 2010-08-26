// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "DofID.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/angleid.h"

#include "SireMol/atomidx.h"

#include "SireMol/bondid.h"

#include "SireMol/connectivity.h"

#include "SireMol/dihedralid.h"

#include "SireMol/molecule.h"

#include "SireMol/moleditor.h"

#include "SireMol/mover.hpp"

#include "SireMol/partialmolecule.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireSystem/system.h"

#include "SireUnits/dimensions.h"

#include "SireUnits/temperature.h"

#include "SireUnits/units.h"

#include "ensemble.h"

#include "movermove.h"

#include <QDebug>

#include <QTime>

#include "movermove.h"

SireMove::DofID __copy__(const SireMove::DofID &other){ return SireMove::DofID(other); }

#include "Qt/qdatastream.hpp"

const char* pvt_get_name(const SireMove::DofID&){ return "SireMove::DofID";}

void register_DofID_class(){

    { //::SireMove::DofID
        typedef bp::class_< SireMove::DofID > DofID_exposer_t;
        DofID_exposer_t DofID_exposer = DofID_exposer_t( "DofID", bp::init< bp::optional< qint32, qint32, qint32, qint32 > >(( bp::arg("index0")=(int)(-0x000000001), bp::arg("index1")=(int)(-0x000000001), bp::arg("index2")=(int)(-0x000000001), bp::arg("index3")=(int)(-0x000000001) )) );
        bp::scope DofID_scope( DofID_exposer );
        DofID_exposer.def( bp::init< SireMol::AtomID const &, SireMol::AtomID const & >(( bp::arg("atom0"), bp::arg("atom1") )) );
        DofID_exposer.def( bp::init< SireMol::AtomID const &, SireMol::AtomID const &, SireMol::AtomID const & >(( bp::arg("atom0"), bp::arg("atom1"), bp::arg("atom2") )) );
        DofID_exposer.def( bp::init< SireMol::AtomID const &, SireMol::AtomID const &, SireMol::AtomID const &, SireMol::AtomID const & >(( bp::arg("atom0"), bp::arg("atom1"), bp::arg("atom2"), bp::arg("atom3") )) );
        DofID_exposer.def( bp::init< SireMol::AtomIdx const &, SireMol::AtomIdx const & >(( bp::arg("atom0"), bp::arg("atom1") )) );
        DofID_exposer.def( bp::init< SireMol::AtomIdx const &, SireMol::AtomIdx const &, SireMol::AtomIdx const & >(( bp::arg("atom0"), bp::arg("atom1"), bp::arg("atom2") )) );
        DofID_exposer.def( bp::init< SireMol::AtomIdx const &, SireMol::AtomIdx const &, SireMol::AtomIdx const &, SireMol::AtomIdx const & >(( bp::arg("atom0"), bp::arg("atom1"), bp::arg("atom2"), bp::arg("atom3") )) );
        DofID_exposer.def( bp::init< SireMol::BondID const & >(( bp::arg("bond") )) );
        DofID_exposer.def( bp::init< SireMol::AngleID const & >(( bp::arg("angle") )) );
        DofID_exposer.def( bp::init< SireMol::DihedralID const & >(( bp::arg("dihedral") )) );
        DofID_exposer.def( bp::init< SireMove::DofID const & >(( bp::arg("other") )) );
        { //::SireMove::DofID::atom0
        
            typedef ::SireMol::AtomIdx ( ::SireMove::DofID::*atom0_function_type )(  ) const;
            atom0_function_type atom0_function_value( &::SireMove::DofID::atom0 );
            
            DofID_exposer.def( 
                "atom0"
                , atom0_function_value );
        
        }
        { //::SireMove::DofID::atom1
        
            typedef ::SireMol::AtomIdx ( ::SireMove::DofID::*atom1_function_type )(  ) const;
            atom1_function_type atom1_function_value( &::SireMove::DofID::atom1 );
            
            DofID_exposer.def( 
                "atom1"
                , atom1_function_value );
        
        }
        { //::SireMove::DofID::atom2
        
            typedef ::SireMol::AtomIdx ( ::SireMove::DofID::*atom2_function_type )(  ) const;
            atom2_function_type atom2_function_value( &::SireMove::DofID::atom2 );
            
            DofID_exposer.def( 
                "atom2"
                , atom2_function_value );
        
        }
        { //::SireMove::DofID::atom3
        
            typedef ::SireMol::AtomIdx ( ::SireMove::DofID::*atom3_function_type )(  ) const;
            atom3_function_type atom3_function_value( &::SireMove::DofID::atom3 );
            
            DofID_exposer.def( 
                "atom3"
                , atom3_function_value );
        
        }
        DofID_exposer.def( bp::self != bp::self );
        { //::SireMove::DofID::operator=
        
            typedef ::SireMove::DofID & ( ::SireMove::DofID::*assign_function_type )( ::SireMove::DofID const & ) ;
            assign_function_type assign_function_value( &::SireMove::DofID::operator= );
            
            DofID_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        DofID_exposer.def( bp::self == bp::self );
        { //::SireMove::DofID::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMove::DofID::typeName );
            
            DofID_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        DofID_exposer.staticmethod( "typeName" );
        DofID_exposer.def( "__copy__", &__copy__);
        DofID_exposer.def( "__deepcopy__", &__copy__);
        DofID_exposer.def( "clone", &__copy__);
        DofID_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMove::DofID >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        DofID_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMove::DofID >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        DofID_exposer.def( "__str__", &pvt_get_name);
        DofID_exposer.def( "__repr__", &pvt_get_name);
    }

}
