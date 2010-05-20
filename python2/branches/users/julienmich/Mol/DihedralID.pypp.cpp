// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "DihedralID.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/property.h"

#include "SireMaths/torsion.h"

#include "SireMaths/vector.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "atomcoords.h"

#include "dihedralid.h"

#include "moleculedata.h"

#include "moleculeinfodata.h"

#include "dihedralid.h"

SireMol::DihedralID __copy__(const SireMol::DihedralID &other){ return SireMol::DihedralID(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_DihedralID_class(){

    { //::SireMol::DihedralID
        typedef bp::class_< SireMol::DihedralID, bp::bases< SireID::ID > > DihedralID_exposer_t;
        DihedralID_exposer_t DihedralID_exposer = DihedralID_exposer_t( "DihedralID" );
        bp::scope DihedralID_scope( DihedralID_exposer );
        DihedralID_exposer.def( bp::init< >() );
        DihedralID_exposer.def( bp::init< SireMol::AtomID const &, SireMol::AtomID const &, SireMol::AtomID const &, SireMol::AtomID const & >(( bp::arg("atom0"), bp::arg("atom1"), bp::arg("atom2"), bp::arg("atom3") )) );
        DihedralID_exposer.def( bp::init< SireMol::DihedralID const & >(( bp::arg("other") )) );
        { //::SireMol::DihedralID::atom0
        
            typedef ::SireMol::AtomID const & ( ::SireMol::DihedralID::*atom0_function_type )(  ) const;
            atom0_function_type atom0_function_value( &::SireMol::DihedralID::atom0 );
            
            DihedralID_exposer.def( 
                "atom0"
                , atom0_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireMol::DihedralID::atom1
        
            typedef ::SireMol::AtomID const & ( ::SireMol::DihedralID::*atom1_function_type )(  ) const;
            atom1_function_type atom1_function_value( &::SireMol::DihedralID::atom1 );
            
            DihedralID_exposer.def( 
                "atom1"
                , atom1_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireMol::DihedralID::atom2
        
            typedef ::SireMol::AtomID const & ( ::SireMol::DihedralID::*atom2_function_type )(  ) const;
            atom2_function_type atom2_function_value( &::SireMol::DihedralID::atom2 );
            
            DihedralID_exposer.def( 
                "atom2"
                , atom2_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireMol::DihedralID::atom3
        
            typedef ::SireMol::AtomID const & ( ::SireMol::DihedralID::*atom3_function_type )(  ) const;
            atom3_function_type atom3_function_value( &::SireMol::DihedralID::atom3 );
            
            DihedralID_exposer.def( 
                "atom3"
                , atom3_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireMol::DihedralID::hash
        
            typedef ::uint ( ::SireMol::DihedralID::*hash_function_type )(  ) const;
            hash_function_type hash_function_value( &::SireMol::DihedralID::hash );
            
            DihedralID_exposer.def( 
                "hash"
                , hash_function_value );
        
        }
        { //::SireMol::DihedralID::isNull
        
            typedef bool ( ::SireMol::DihedralID::*isNull_function_type )(  ) const;
            isNull_function_type isNull_function_value( &::SireMol::DihedralID::isNull );
            
            DihedralID_exposer.def( 
                "isNull"
                , isNull_function_value );
        
        }
        { //::SireMol::DihedralID::map
        
            typedef ::boost::tuples::tuple< SireMol::AtomIdx, SireMol::AtomIdx, SireMol::AtomIdx, SireMol::AtomIdx, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type > ( ::SireMol::DihedralID::*map_function_type )( ::SireMol::MoleculeInfoData const & ) const;
            map_function_type map_function_value( &::SireMol::DihedralID::map );
            
            DihedralID_exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("molinfo") ) );
        
        }
        { //::SireMol::DihedralID::map
        
            typedef ::boost::tuples::tuple< SireMol::AtomIdx, SireMol::AtomIdx, SireMol::AtomIdx, SireMol::AtomIdx, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type > ( ::SireMol::DihedralID::*map_function_type )( ::SireMol::MoleculeInfoData const &,::SireMol::MoleculeInfoData const &,::SireMol::MoleculeInfoData const &,::SireMol::MoleculeInfoData const & ) const;
            map_function_type map_function_value( &::SireMol::DihedralID::map );
            
            DihedralID_exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("mol0info"), bp::arg("mol1info"), bp::arg("mol2info"), bp::arg("mol3info") ) );
        
        }
        DihedralID_exposer.def( bp::self != bp::self );
        { //::SireMol::DihedralID::operator=
        
            typedef ::SireMol::DihedralID & ( ::SireMol::DihedralID::*assign_function_type )( ::SireMol::DihedralID const & ) ;
            assign_function_type assign_function_value( &::SireMol::DihedralID::operator= );
            
            DihedralID_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        DihedralID_exposer.def( bp::self == bp::other< SireID::ID >() );
        DihedralID_exposer.def( bp::self == bp::self );
        { //::SireMol::DihedralID::size
        
            typedef ::SireUnits::Dimension::Angle ( ::SireMol::DihedralID::*size_function_type )( ::SireMol::MoleculeData const &,::SireBase::PropertyMap const & ) const;
            size_function_type size_function_value( &::SireMol::DihedralID::size );
            
            DihedralID_exposer.def( 
                "size"
                , size_function_value
                , ( bp::arg("moldata"), bp::arg("map")=SireBase::PropertyMap() ) );
        
        }
        { //::SireMol::DihedralID::size
        
            typedef ::SireUnits::Dimension::Angle ( ::SireMol::DihedralID::*size_function_type )( ::SireMol::MoleculeData const &,::SireMol::MoleculeData const &,::SireMol::MoleculeData const &,::SireMol::MoleculeData const &,::SireBase::PropertyMap const & ) const;
            size_function_type size_function_value( &::SireMol::DihedralID::size );
            
            DihedralID_exposer.def( 
                "size"
                , size_function_value
                , ( bp::arg("mol0data"), bp::arg("mol1data"), bp::arg("mol2data"), bp::arg("mol3data"), bp::arg("map")=SireBase::PropertyMap() ) );
        
        }
        { //::SireMol::DihedralID::size
        
            typedef ::SireUnits::Dimension::Angle ( ::SireMol::DihedralID::*size_function_type )( ::SireMol::MoleculeData const &,::SireBase::PropertyMap const &,::SireMol::MoleculeData const &,::SireBase::PropertyMap const &,::SireMol::MoleculeData const &,::SireBase::PropertyMap const &,::SireMol::MoleculeData const &,::SireBase::PropertyMap const & ) const;
            size_function_type size_function_value( &::SireMol::DihedralID::size );
            
            DihedralID_exposer.def( 
                "size"
                , size_function_value
                , ( bp::arg("mol0data"), bp::arg("map0"), bp::arg("mol1data"), bp::arg("map1"), bp::arg("mol2data"), bp::arg("map2"), bp::arg("mol3data"), bp::arg("map3") ) );
        
        }
        { //::SireMol::DihedralID::toString
        
            typedef ::QString ( ::SireMol::DihedralID::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireMol::DihedralID::toString );
            
            DihedralID_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireMol::DihedralID::torsion
        
            typedef ::SireMaths::Torsion ( ::SireMol::DihedralID::*torsion_function_type )( ::SireMol::MoleculeData const &,::SireBase::PropertyMap const & ) const;
            torsion_function_type torsion_function_value( &::SireMol::DihedralID::torsion );
            
            DihedralID_exposer.def( 
                "torsion"
                , torsion_function_value
                , ( bp::arg("moldata"), bp::arg("map")=SireBase::PropertyMap() ) );
        
        }
        { //::SireMol::DihedralID::torsion
        
            typedef ::SireMaths::Torsion ( ::SireMol::DihedralID::*torsion_function_type )( ::SireMol::MoleculeData const &,::SireMol::MoleculeData const &,::SireMol::MoleculeData const &,::SireMol::MoleculeData const &,::SireBase::PropertyMap const & ) const;
            torsion_function_type torsion_function_value( &::SireMol::DihedralID::torsion );
            
            DihedralID_exposer.def( 
                "torsion"
                , torsion_function_value
                , ( bp::arg("mol0data"), bp::arg("mol1data"), bp::arg("mol2data"), bp::arg("mol3data"), bp::arg("map")=SireBase::PropertyMap() ) );
        
        }
        { //::SireMol::DihedralID::torsion
        
            typedef ::SireMaths::Torsion ( ::SireMol::DihedralID::*torsion_function_type )( ::SireMol::MoleculeData const &,::SireBase::PropertyMap const &,::SireMol::MoleculeData const &,::SireBase::PropertyMap const &,::SireMol::MoleculeData const &,::SireBase::PropertyMap const &,::SireMol::MoleculeData const &,::SireBase::PropertyMap const & ) const;
            torsion_function_type torsion_function_value( &::SireMol::DihedralID::torsion );
            
            DihedralID_exposer.def( 
                "torsion"
                , torsion_function_value
                , ( bp::arg("mol0data"), bp::arg("map0"), bp::arg("mol1data"), bp::arg("map1"), bp::arg("mol2data"), bp::arg("map2"), bp::arg("mol3data"), bp::arg("map3") ) );
        
        }
        { //::SireMol::DihedralID::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::DihedralID::typeName );
            
            DihedralID_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireMol::DihedralID::what
        
            typedef char const * ( ::SireMol::DihedralID::*what_function_type )(  ) const;
            what_function_type what_function_value( &::SireMol::DihedralID::what );
            
            DihedralID_exposer.def( 
                "what"
                , what_function_value );
        
        }
        DihedralID_exposer.staticmethod( "typeName" );
        DihedralID_exposer.def( "__copy__", &__copy__);
        DihedralID_exposer.def( "__deepcopy__", &__copy__);
        DihedralID_exposer.def( "clone", &__copy__);
        DihedralID_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::DihedralID >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        DihedralID_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::DihedralID >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        DihedralID_exposer.def( "__str__", &__str__< ::SireMol::DihedralID > );
        DihedralID_exposer.def( "__repr__", &__str__< ::SireMol::DihedralID > );
    }

}
