// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Flexibility.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireMol/angleid.h"

#include "SireMol/atomidx.h"

#include "SireMol/bondid.h"

#include "SireMol/dihedralid.h"

#include "SireMol/errors.h"

#include "SireMol/molecule.h"

#include "SireMol/mover.hpp"

#include "SireMol/partialmolecule.h"

#include "SireMove/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireUnits/convert.h"

#include "SireUnits/units.h"

#include "flexibility.h"

#include <QList>

#include "flexibility.h"

SireMove::Flexibility __copy__(const SireMove::Flexibility &other){ return SireMove::Flexibility(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_Flexibility_class(){

    { //::SireMove::Flexibility
        typedef bp::class_< SireMove::Flexibility, bp::bases< SireMol::MoleculeProperty, SireMol::MolViewProperty, SireBase::Property > > Flexibility_exposer_t;
        Flexibility_exposer_t Flexibility_exposer = Flexibility_exposer_t( "Flexibility", bp::init< >() );
        bp::scope Flexibility_scope( Flexibility_exposer );
        Flexibility_exposer.def( bp::init< SireMol::MoleculeData const & >(( bp::arg("molecule") )) );
        Flexibility_exposer.def( bp::init< SireMove::Flexibility const & >(( bp::arg("other") )) );
        { //::SireMove::Flexibility::add
        
            typedef void ( ::SireMove::Flexibility::*add_function_type )( ::SireMol::BondID const &,::SireUnits::Dimension::Length const & ) ;
            add_function_type add_function_value( &::SireMove::Flexibility::add );
            
            Flexibility_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("bond"), bp::arg("delta") ) );
        
        }
        { //::SireMove::Flexibility::add
        
            typedef void ( ::SireMove::Flexibility::*add_function_type )( ::SireMol::AngleID const &,::SireUnits::Dimension::Angle const & ) ;
            add_function_type add_function_value( &::SireMove::Flexibility::add );
            
            Flexibility_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("angle"), bp::arg("delta") ) );
        
        }
        { //::SireMove::Flexibility::add
        
            typedef void ( ::SireMove::Flexibility::*add_function_type )( ::SireMol::DihedralID const &,::SireUnits::Dimension::Angle const & ) ;
            add_function_type add_function_value( &::SireMove::Flexibility::add );
            
            Flexibility_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("dihedral"), bp::arg("delta") ) );
        
        }
        { //::SireMove::Flexibility::contains
        
            typedef bool ( ::SireMove::Flexibility::*contains_function_type )( ::SireMol::BondID const & ) const;
            contains_function_type contains_function_value( &::SireMove::Flexibility::contains );
            
            Flexibility_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("bond") ) );
        
        }
        { //::SireMove::Flexibility::contains
        
            typedef bool ( ::SireMove::Flexibility::*contains_function_type )( ::SireMol::AngleID const & ) const;
            contains_function_type contains_function_value( &::SireMove::Flexibility::contains );
            
            Flexibility_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("angle") ) );
        
        }
        { //::SireMove::Flexibility::contains
        
            typedef bool ( ::SireMove::Flexibility::*contains_function_type )( ::SireMol::DihedralID const & ) const;
            contains_function_type contains_function_value( &::SireMove::Flexibility::contains );
            
            Flexibility_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("dihedral") ) );
        
        }
        { //::SireMove::Flexibility::delta
        
            typedef ::SireUnits::Dimension::Length ( ::SireMove::Flexibility::*delta_function_type )( ::SireMol::BondID const & ) const;
            delta_function_type delta_function_value( &::SireMove::Flexibility::delta );
            
            Flexibility_exposer.def( 
                "delta"
                , delta_function_value
                , ( bp::arg("bond") ) );
        
        }
        { //::SireMove::Flexibility::delta
        
            typedef ::SireUnits::Dimension::Angle ( ::SireMove::Flexibility::*delta_function_type )( ::SireMol::AngleID const & ) const;
            delta_function_type delta_function_value( &::SireMove::Flexibility::delta );
            
            Flexibility_exposer.def( 
                "delta"
                , delta_function_value
                , ( bp::arg("angle") ) );
        
        }
        { //::SireMove::Flexibility::delta
        
            typedef ::SireUnits::Dimension::Angle ( ::SireMove::Flexibility::*delta_function_type )( ::SireMol::DihedralID const & ) const;
            delta_function_type delta_function_value( &::SireMove::Flexibility::delta );
            
            Flexibility_exposer.def( 
                "delta"
                , delta_function_value
                , ( bp::arg("dihedral") ) );
        
        }
        { //::SireMove::Flexibility::flexibleAngles
        
            typedef ::QList< SireMol::AngleID > ( ::SireMove::Flexibility::*flexibleAngles_function_type )(  ) const;
            flexibleAngles_function_type flexibleAngles_function_value( &::SireMove::Flexibility::flexibleAngles );
            
            Flexibility_exposer.def( 
                "flexibleAngles"
                , flexibleAngles_function_value );
        
        }
        { //::SireMove::Flexibility::flexibleBonds
        
            typedef ::QList< SireMol::BondID > ( ::SireMove::Flexibility::*flexibleBonds_function_type )(  ) const;
            flexibleBonds_function_type flexibleBonds_function_value( &::SireMove::Flexibility::flexibleBonds );
            
            Flexibility_exposer.def( 
                "flexibleBonds"
                , flexibleBonds_function_value );
        
        }
        { //::SireMove::Flexibility::flexibleDihedrals
        
            typedef ::QList< SireMol::DihedralID > ( ::SireMove::Flexibility::*flexibleDihedrals_function_type )(  ) const;
            flexibleDihedrals_function_type flexibleDihedrals_function_value( &::SireMove::Flexibility::flexibleDihedrals );
            
            Flexibility_exposer.def( 
                "flexibleDihedrals"
                , flexibleDihedrals_function_value );
        
        }
        { //::SireMove::Flexibility::info
        
            typedef ::SireMol::MoleculeInfoData const & ( ::SireMove::Flexibility::*info_function_type )(  ) const;
            info_function_type info_function_value( &::SireMove::Flexibility::info );
            
            Flexibility_exposer.def( 
                "info"
                , info_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMove::Flexibility::isCompatibleWith
        
            typedef bool ( ::SireMove::Flexibility::*isCompatibleWith_function_type )( ::SireMol::MoleculeInfoData const & ) const;
            isCompatibleWith_function_type isCompatibleWith_function_value( &::SireMove::Flexibility::isCompatibleWith );
            
            Flexibility_exposer.def( 
                "isCompatibleWith"
                , isCompatibleWith_function_value
                , ( bp::arg("molinfo") ) );
        
        }
        { //::SireMove::Flexibility::maximumAngleVar
        
            typedef int ( ::SireMove::Flexibility::*maximumAngleVar_function_type )(  ) const;
            maximumAngleVar_function_type maximumAngleVar_function_value( &::SireMove::Flexibility::maximumAngleVar );
            
            Flexibility_exposer.def( 
                "maximumAngleVar"
                , maximumAngleVar_function_value );
        
        }
        { //::SireMove::Flexibility::maximumBondVar
        
            typedef int ( ::SireMove::Flexibility::*maximumBondVar_function_type )(  ) const;
            maximumBondVar_function_type maximumBondVar_function_value( &::SireMove::Flexibility::maximumBondVar );
            
            Flexibility_exposer.def( 
                "maximumBondVar"
                , maximumBondVar_function_value );
        
        }
        { //::SireMove::Flexibility::maximumDihedralVar
        
            typedef int ( ::SireMove::Flexibility::*maximumDihedralVar_function_type )(  ) const;
            maximumDihedralVar_function_type maximumDihedralVar_function_value( &::SireMove::Flexibility::maximumDihedralVar );
            
            Flexibility_exposer.def( 
                "maximumDihedralVar"
                , maximumDihedralVar_function_value );
        
        }
        Flexibility_exposer.def( bp::self != bp::self );
        { //::SireMove::Flexibility::operator=
        
            typedef ::SireMove::Flexibility & ( ::SireMove::Flexibility::*assign_function_type )( ::SireMove::Flexibility const & ) ;
            assign_function_type assign_function_value( &::SireMove::Flexibility::operator= );
            
            Flexibility_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        Flexibility_exposer.def( bp::self == bp::self );
        { //::SireMove::Flexibility::remove
        
            typedef void ( ::SireMove::Flexibility::*remove_function_type )( ::SireMol::BondID const & ) ;
            remove_function_type remove_function_value( &::SireMove::Flexibility::remove );
            
            Flexibility_exposer.def( 
                "remove"
                , remove_function_value
                , ( bp::arg("bond") ) );
        
        }
        { //::SireMove::Flexibility::remove
        
            typedef void ( ::SireMove::Flexibility::*remove_function_type )( ::SireMol::AngleID const & ) ;
            remove_function_type remove_function_value( &::SireMove::Flexibility::remove );
            
            Flexibility_exposer.def( 
                "remove"
                , remove_function_value
                , ( bp::arg("angle") ) );
        
        }
        { //::SireMove::Flexibility::remove
        
            typedef void ( ::SireMove::Flexibility::*remove_function_type )( ::SireMol::DihedralID const & ) ;
            remove_function_type remove_function_value( &::SireMove::Flexibility::remove );
            
            Flexibility_exposer.def( 
                "remove"
                , remove_function_value
                , ( bp::arg("dihedral") ) );
        
        }
        { //::SireMove::Flexibility::rotation
        
            typedef ::SireUnits::Dimension::Angle ( ::SireMove::Flexibility::*rotation_function_type )(  ) const;
            rotation_function_type rotation_function_value( &::SireMove::Flexibility::rotation );
            
            Flexibility_exposer.def( 
                "rotation"
                , rotation_function_value );
        
        }
        { //::SireMove::Flexibility::setDelta
        
            typedef void ( ::SireMove::Flexibility::*setDelta_function_type )( ::SireMol::BondID const &,::SireUnits::Dimension::Length const & ) ;
            setDelta_function_type setDelta_function_value( &::SireMove::Flexibility::setDelta );
            
            Flexibility_exposer.def( 
                "setDelta"
                , setDelta_function_value
                , ( bp::arg("bond"), bp::arg("delta") ) );
        
        }
        { //::SireMove::Flexibility::setDelta
        
            typedef void ( ::SireMove::Flexibility::*setDelta_function_type )( ::SireMol::AngleID const &,::SireUnits::Dimension::Angle const & ) ;
            setDelta_function_type setDelta_function_value( &::SireMove::Flexibility::setDelta );
            
            Flexibility_exposer.def( 
                "setDelta"
                , setDelta_function_value
                , ( bp::arg("angle"), bp::arg("delta") ) );
        
        }
        { //::SireMove::Flexibility::setDelta
        
            typedef void ( ::SireMove::Flexibility::*setDelta_function_type )( ::SireMol::DihedralID const &,::SireUnits::Dimension::Angle const & ) ;
            setDelta_function_type setDelta_function_value( &::SireMove::Flexibility::setDelta );
            
            Flexibility_exposer.def( 
                "setDelta"
                , setDelta_function_value
                , ( bp::arg("dihedral"), bp::arg("delta") ) );
        
        }
        { //::SireMove::Flexibility::setMaximumAngleVar
        
            typedef void ( ::SireMove::Flexibility::*setMaximumAngleVar_function_type )( int ) ;
            setMaximumAngleVar_function_type setMaximumAngleVar_function_value( &::SireMove::Flexibility::setMaximumAngleVar );
            
            Flexibility_exposer.def( 
                "setMaximumAngleVar"
                , setMaximumAngleVar_function_value
                , ( bp::arg("maxvar") ) );
        
        }
        { //::SireMove::Flexibility::setMaximumBondVar
        
            typedef void ( ::SireMove::Flexibility::*setMaximumBondVar_function_type )( int ) ;
            setMaximumBondVar_function_type setMaximumBondVar_function_value( &::SireMove::Flexibility::setMaximumBondVar );
            
            Flexibility_exposer.def( 
                "setMaximumBondVar"
                , setMaximumBondVar_function_value
                , ( bp::arg("maxvar") ) );
        
        }
        { //::SireMove::Flexibility::setMaximumDihedralVar
        
            typedef void ( ::SireMove::Flexibility::*setMaximumDihedralVar_function_type )( int ) ;
            setMaximumDihedralVar_function_type setMaximumDihedralVar_function_value( &::SireMove::Flexibility::setMaximumDihedralVar );
            
            Flexibility_exposer.def( 
                "setMaximumDihedralVar"
                , setMaximumDihedralVar_function_value
                , ( bp::arg("maxvar") ) );
        
        }
        { //::SireMove::Flexibility::setRotation
        
            typedef void ( ::SireMove::Flexibility::*setRotation_function_type )( ::SireUnits::Dimension::Angle const & ) ;
            setRotation_function_type setRotation_function_value( &::SireMove::Flexibility::setRotation );
            
            Flexibility_exposer.def( 
                "setRotation"
                , setRotation_function_value
                , ( bp::arg("rotation") ) );
        
        }
        { //::SireMove::Flexibility::setTranslation
        
            typedef void ( ::SireMove::Flexibility::*setTranslation_function_type )( ::SireUnits::Dimension::Length const & ) ;
            setTranslation_function_type setTranslation_function_value( &::SireMove::Flexibility::setTranslation );
            
            Flexibility_exposer.def( 
                "setTranslation"
                , setTranslation_function_value
                , ( bp::arg("translation") ) );
        
        }
        { //::SireMove::Flexibility::toString
        
            typedef ::QString ( ::SireMove::Flexibility::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireMove::Flexibility::toString );
            
            Flexibility_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireMove::Flexibility::translation
        
            typedef ::SireUnits::Dimension::Length ( ::SireMove::Flexibility::*translation_function_type )(  ) const;
            translation_function_type translation_function_value( &::SireMove::Flexibility::translation );
            
            Flexibility_exposer.def( 
                "translation"
                , translation_function_value );
        
        }
        { //::SireMove::Flexibility::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMove::Flexibility::typeName );
            
            Flexibility_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        Flexibility_exposer.staticmethod( "typeName" );
        Flexibility_exposer.def( "__copy__", &__copy__);
        Flexibility_exposer.def( "__deepcopy__", &__copy__);
        Flexibility_exposer.def( "clone", &__copy__);
        Flexibility_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMove::Flexibility >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Flexibility_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMove::Flexibility >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Flexibility_exposer.def( "__str__", &__str__< ::SireMove::Flexibility > );
        Flexibility_exposer.def( "__repr__", &__str__< ::SireMove::Flexibility > );
    }

}
