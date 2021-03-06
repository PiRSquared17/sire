// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "MGIDsAndMaps.pypp.hpp"

namespace bp = boost::python;

#include "SireID/idorset.hpp"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "mgidsandmaps.h"

#include "mgidsandmaps.h"

SireMol::MGIDsAndMaps __copy__(const SireMol::MGIDsAndMaps &other){ return SireMol::MGIDsAndMaps(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/len.hpp"

void register_MGIDsAndMaps_class(){

    { //::SireMol::MGIDsAndMaps
        typedef bp::class_< SireMol::MGIDsAndMaps, bp::bases< SireBase::Property > > MGIDsAndMaps_exposer_t;
        MGIDsAndMaps_exposer_t MGIDsAndMaps_exposer = MGIDsAndMaps_exposer_t( "MGIDsAndMaps", bp::init< >() );
        bp::scope MGIDsAndMaps_scope( MGIDsAndMaps_exposer );
        MGIDsAndMaps_exposer.def( bp::init< SireMol::MoleculeGroup const & >(( bp::arg("mgroup") )) );
        MGIDsAndMaps_exposer.def( bp::init< SireMol::MoleculeGroup const &, SireBase::PropertyMap const & >(( bp::arg("mgroup"), bp::arg("map") )) );
        MGIDsAndMaps_exposer.def( bp::init< boost::tuples::tuple< SireBase::PropPtr< SireMol::MoleculeGroup >, SireBase::PropertyMap, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type > >(( bp::arg("group_and_map") )) );
        MGIDsAndMaps_exposer.def( bp::init< QList< SireBase::PropPtr< SireMol::MoleculeGroup > > const & >(( bp::arg("mgroups") )) );
        MGIDsAndMaps_exposer.def( bp::init< QList< boost::tuples::tuple< SireBase::PropPtr< SireMol::MoleculeGroup >, SireBase::PropertyMap, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type > > const & >(( bp::arg("groups_and_maps") )) );
        MGIDsAndMaps_exposer.def( bp::init< QList< SireBase::PropPtr< SireMol::MoleculeGroup > > const &, SireBase::PropertyMap const & >(( bp::arg("mgroups"), bp::arg("map") )) );
        MGIDsAndMaps_exposer.def( bp::init< boost::tuples::tuple< QList< SireBase::PropPtr< SireMol::MoleculeGroup > >, SireBase::PropertyMap, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type > const & >(( bp::arg("groups_and_maps") )) );
        MGIDsAndMaps_exposer.def( bp::init< SireMol::MGID const & >(( bp::arg("mgid") )) );
        MGIDsAndMaps_exposer.def( bp::init< SireMol::MGID const &, SireBase::PropertyMap const & >(( bp::arg("mgid"), bp::arg("map") )) );
        MGIDsAndMaps_exposer.def( bp::init< boost::tuples::tuple< SireMol::MGIdentifier, SireBase::PropertyMap, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type > >(( bp::arg("mgid_and_map") )) );
        MGIDsAndMaps_exposer.def( bp::init< QList< SireMol::MGIdentifier > const & >(( bp::arg("mgids") )) );
        MGIDsAndMaps_exposer.def( bp::init< QList< boost::tuples::tuple< SireMol::MGIdentifier, SireBase::PropertyMap, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type > > const & >(( bp::arg("mgids_and_maps") )) );
        MGIDsAndMaps_exposer.def( bp::init< QList< SireMol::MGIdentifier > const &, SireBase::PropertyMap const & >(( bp::arg("mgids"), bp::arg("map") )) );
        MGIDsAndMaps_exposer.def( bp::init< boost::tuples::tuple< QList< SireMol::MGIdentifier >, SireBase::PropertyMap, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type > const & >(( bp::arg("mgids_and_maps") )) );
        MGIDsAndMaps_exposer.def( bp::init< QList< SireMol::MGIDsAndMaps > const & >(( bp::arg("mgids_and_maps") )) );
        MGIDsAndMaps_exposer.def( bp::init< SireMol::MGIDsAndMaps const & >(( bp::arg("other") )) );
        { //::SireMol::MGIDsAndMaps::count
        
            typedef int ( ::SireMol::MGIDsAndMaps::*count_function_type )(  ) const;
            count_function_type count_function_value( &::SireMol::MGIDsAndMaps::count );
            
            MGIDsAndMaps_exposer.def( 
                "count"
                , count_function_value );
        
        }
        { //::SireMol::MGIDsAndMaps::isEmpty
        
            typedef bool ( ::SireMol::MGIDsAndMaps::*isEmpty_function_type )(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMol::MGIDsAndMaps::isEmpty );
            
            MGIDsAndMaps_exposer.def( 
                "isEmpty"
                , isEmpty_function_value );
        
        }
        { //::SireMol::MGIDsAndMaps::mgIDs
        
            typedef ::QVector< SireMol::MGIdentifier > const & ( ::SireMol::MGIDsAndMaps::*mgIDs_function_type )(  ) const;
            mgIDs_function_type mgIDs_function_value( &::SireMol::MGIDsAndMaps::mgIDs );
            
            MGIDsAndMaps_exposer.def( 
                "mgIDs"
                , mgIDs_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        MGIDsAndMaps_exposer.def( bp::self != bp::self );
        { //::SireMol::MGIDsAndMaps::operator=
        
            typedef ::SireMol::MGIDsAndMaps & ( ::SireMol::MGIDsAndMaps::*assign_function_type )( ::SireMol::MGIDsAndMaps const & ) ;
            assign_function_type assign_function_value( &::SireMol::MGIDsAndMaps::operator= );
            
            MGIDsAndMaps_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        MGIDsAndMaps_exposer.def( bp::self == bp::self );
        { //::SireMol::MGIDsAndMaps::propertyMaps
        
            typedef ::QVector< SireBase::PropertyMap > const & ( ::SireMol::MGIDsAndMaps::*propertyMaps_function_type )(  ) const;
            propertyMaps_function_type propertyMaps_function_value( &::SireMol::MGIDsAndMaps::propertyMaps );
            
            MGIDsAndMaps_exposer.def( 
                "propertyMaps"
                , propertyMaps_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMol::MGIDsAndMaps::toString
        
            typedef ::QString ( ::SireMol::MGIDsAndMaps::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireMol::MGIDsAndMaps::toString );
            
            MGIDsAndMaps_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireMol::MGIDsAndMaps::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::MGIDsAndMaps::typeName );
            
            MGIDsAndMaps_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        MGIDsAndMaps_exposer.staticmethod( "typeName" );
        MGIDsAndMaps_exposer.def( "__copy__", &__copy__);
        MGIDsAndMaps_exposer.def( "__deepcopy__", &__copy__);
        MGIDsAndMaps_exposer.def( "clone", &__copy__);
        MGIDsAndMaps_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::MGIDsAndMaps >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        MGIDsAndMaps_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::MGIDsAndMaps >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        MGIDsAndMaps_exposer.def( "__str__", &__str__< ::SireMol::MGIDsAndMaps > );
        MGIDsAndMaps_exposer.def( "__repr__", &__str__< ::SireMol::MGIDsAndMaps > );
        MGIDsAndMaps_exposer.def( "__len__", &__len_count< ::SireMol::MGIDsAndMaps > );
    }

}
