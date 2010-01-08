// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "CoordGroup.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/quickcopy.hpp"

#include "SireMaths/axisset.h"

#include "SireMaths/matrix.h"

#include "SireMaths/quaternion.h"

#include "SireMaths/rotate.h"

#include "Siren/errors.h"

#include "Siren/stream.h"

#include "coordgroup.h"

#include <QDebug>

#include "coordgroup.h"

#include "Siren/logger.h"

SireVol::CoordGroup __copy__(const SireVol::CoordGroup &other){ return SireVol::CoordGroup(other); }

#include "Siren/str.hpp"

void register_CoordGroup_class(){

    { //::SireVol::CoordGroup
        typedef bp::class_< SireVol::CoordGroup, bp::bases< SireVol::CoordGroupBase, Siren::Object > > CoordGroup_exposer_t;
        CoordGroup_exposer_t CoordGroup_exposer = CoordGroup_exposer_t( "CoordGroup" );
        bp::scope CoordGroup_scope( CoordGroup_exposer );
        CoordGroup_exposer.def( bp::init< >() );
        CoordGroup_exposer.def( bp::init< quint32 >(( bp::arg("size") )) );
        CoordGroup_exposer.def( bp::init< quint32, SireMaths::Vector const & >(( bp::arg("size"), bp::arg("value") )) );
        CoordGroup_exposer.def( bp::init< quint32, SireMaths::Vector const * >(( bp::arg("size"), bp::arg("values") )) );
        CoordGroup_exposer.def( bp::init< SireVol::CoordGroupArray const & >(( bp::arg("cgarray") )) );
        CoordGroup_exposer.def( bp::init< SireVol::CoordGroupArrayArray const & >(( bp::arg("cgarrays") )) );
        CoordGroup_exposer.def( bp::init< QVector< SireMaths::Vector > const & >(( bp::arg("points") )) );
        CoordGroup_exposer.def( bp::init< SireVol::CoordGroupEditor const & >(( bp::arg("editor") )) );
        CoordGroup_exposer.def( bp::init< SireVol::CoordGroup const & >(( bp::arg("other") )) );
        { //::SireVol::CoordGroup::edit
        
            typedef ::SireVol::CoordGroupEditor ( ::SireVol::CoordGroup::*edit_function_type )(  ) const;
            edit_function_type edit_function_value( &::SireVol::CoordGroup::edit );
            
            CoordGroup_exposer.def( 
                "edit"
                , edit_function_value );
        
        }
        { //::SireVol::CoordGroup::hashCode
        
            typedef ::uint ( ::SireVol::CoordGroup::*hashCode_function_type )(  ) const;
            hashCode_function_type hashCode_function_value( &::SireVol::CoordGroup::hashCode );
            
            CoordGroup_exposer.def( 
                "hashCode"
                , hashCode_function_value );
        
        }
        { //::SireVol::CoordGroup::operator=
        
            typedef ::SireVol::CoordGroup & ( ::SireVol::CoordGroup::*assign_function_type )( ::SireVol::CoordGroup const & ) ;
            assign_function_type assign_function_value( &::SireVol::CoordGroup::operator= );
            
            CoordGroup_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        { //::SireVol::CoordGroup::operator=
        
            typedef ::SireVol::CoordGroup & ( ::SireVol::CoordGroup::*assign_function_type )( ::SireVol::CoordGroupEditor const & ) ;
            assign_function_type assign_function_value( &::SireVol::CoordGroup::operator= );
            
            CoordGroup_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        { //::SireVol::CoordGroup::stream
        
            typedef void ( ::SireVol::CoordGroup::*stream_function_type )( ::Siren::Stream & ) ;
            stream_function_type stream_function_value( &::SireVol::CoordGroup::stream );
            
            CoordGroup_exposer.def( 
                "stream"
                , stream_function_value
                , ( bp::arg("s") ) );
        
        }
        { //::SireVol::CoordGroup::toString
        
            typedef ::QString ( ::SireVol::CoordGroup::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireVol::CoordGroup::toString );
            
            CoordGroup_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        CoordGroup_exposer.def( "__copy__", &__copy__);
        CoordGroup_exposer.def( "__deepcopy__", &__copy__);
        CoordGroup_exposer.def( "clone", &__copy__);
        CoordGroup_exposer.def( "__str__", &__str__< ::SireVol::CoordGroup > );
        CoordGroup_exposer.def( "__repr__", &__str__< ::SireVol::CoordGroup > );
    }

}
