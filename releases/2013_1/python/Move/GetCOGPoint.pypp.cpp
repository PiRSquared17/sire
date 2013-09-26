// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "GetCOGPoint.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireMaths/vector.h"

#include "SireMol/atom.h"

#include "SireMol/atomcoords.h"

#include "SireMol/evaluator.h"

#include "SireMol/moleculeview.h"

#include "SireMol/mover.hpp"

#include "SireMol/selector.hpp"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "getpoint.h"

#include <QDebug>

#include "getpoint.h"

SireMove::GetCOGPoint __copy__(const SireMove::GetCOGPoint &other){ return SireMove::GetCOGPoint(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_GetCOGPoint_class(){

    { //::SireMove::GetCOGPoint
        typedef bp::class_< SireMove::GetCOGPoint, bp::bases< SireMove::GetPoint, SireBase::Property > > GetCOGPoint_exposer_t;
        GetCOGPoint_exposer_t GetCOGPoint_exposer = GetCOGPoint_exposer_t( "GetCOGPoint", bp::init< >() );
        bp::scope GetCOGPoint_scope( GetCOGPoint_exposer );
        GetCOGPoint_exposer.def( bp::init< SireMol::AtomID const & >(( bp::arg("atomid") )) );
        GetCOGPoint_exposer.def( bp::init< SireMol::AtomID const &, SireMol::AtomID const & >(( bp::arg("atomid0"), bp::arg("atomid1") )) );
        GetCOGPoint_exposer.def( bp::init< QList< SireMol::AtomIdentifier > const & >(( bp::arg("atomids") )) );
        GetCOGPoint_exposer.def( bp::init< SireMove::GetCOGPoint const & >(( bp::arg("other") )) );
        { //::SireMove::GetCOGPoint::atomID
        
            typedef ::SireMol::AtomID const & ( ::SireMove::GetCOGPoint::*atomID_function_type )(  ) const;
            atomID_function_type atomID_function_value( &::SireMove::GetCOGPoint::atomID );
            
            GetCOGPoint_exposer.def( 
                "atomID"
                , atomID_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireMove::GetCOGPoint::getPoint
        
            typedef ::SireMaths::Vector ( ::SireMove::GetCOGPoint::*getPoint_function_type )( ::SireMol::MoleculeView const &,::SireBase::PropertyMap const & ) const;
            getPoint_function_type getPoint_function_value( &::SireMove::GetCOGPoint::getPoint );
            
            GetCOGPoint_exposer.def( 
                "getPoint"
                , getPoint_function_value
                , ( bp::arg("molecule"), bp::arg("map")=SireBase::PropertyMap() ) );
        
        }
        GetCOGPoint_exposer.def( bp::self != bp::self );
        { //::SireMove::GetCOGPoint::operator=
        
            typedef ::SireMove::GetCOGPoint & ( ::SireMove::GetCOGPoint::*assign_function_type )( ::SireMove::GetCOGPoint const & ) ;
            assign_function_type assign_function_value( &::SireMove::GetCOGPoint::operator= );
            
            GetCOGPoint_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        GetCOGPoint_exposer.def( bp::self == bp::self );
        { //::SireMove::GetCOGPoint::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMove::GetCOGPoint::typeName );
            
            GetCOGPoint_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        GetCOGPoint_exposer.staticmethod( "typeName" );
        GetCOGPoint_exposer.def( "__copy__", &__copy__);
        GetCOGPoint_exposer.def( "__deepcopy__", &__copy__);
        GetCOGPoint_exposer.def( "clone", &__copy__);
        GetCOGPoint_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMove::GetCOGPoint >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        GetCOGPoint_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMove::GetCOGPoint >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        GetCOGPoint_exposer.def( "__str__", &__str__< ::SireMove::GetCOGPoint > );
        GetCOGPoint_exposer.def( "__repr__", &__str__< ::SireMove::GetCOGPoint > );
    }

}
