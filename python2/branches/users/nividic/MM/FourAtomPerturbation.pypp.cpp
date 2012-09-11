// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "FourAtomPerturbation.pypp.hpp"

namespace bp = boost::python;

#include "SireCAS/values.h"

#include "SireMol/molecule.h"

#include "SireMol/moleditor.h"

#include "SireMol/mover.hpp"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "fouratomfunctions.h"

#include "internalperturbation.h"

#include "threeatomfunctions.h"

#include "twoatomfunctions.h"

#include "internalperturbation.h"

SireMM::FourAtomPerturbation __copy__(const SireMM::FourAtomPerturbation &other){ return SireMM::FourAtomPerturbation(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_FourAtomPerturbation_class(){

    { //::SireMM::FourAtomPerturbation
        typedef bp::class_< SireMM::FourAtomPerturbation, bp::bases< SireMM::InternalPerturbation, SireMol::Perturbation, SireBase::Property > > FourAtomPerturbation_exposer_t;
        FourAtomPerturbation_exposer_t FourAtomPerturbation_exposer = FourAtomPerturbation_exposer_t( "FourAtomPerturbation", bp::init< >() );
        bp::scope FourAtomPerturbation_scope( FourAtomPerturbation_exposer );
        FourAtomPerturbation_exposer.def( bp::init< SireMol::AtomID const &, SireMol::AtomID const &, SireMol::AtomID const &, SireMol::AtomID const &, SireCAS::Expression const &, SireCAS::Expression const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("atom0"), bp::arg("atom1"), bp::arg("atom2"), bp::arg("atom3"), bp::arg("initial_form"), bp::arg("final_form"), bp::arg("map")=SireBase::PropertyMap() )) );
        FourAtomPerturbation_exposer.def( bp::init< SireMol::AtomID const &, SireMol::AtomID const &, SireMol::AtomID const &, SireMol::AtomID const &, SireCAS::Expression const &, SireCAS::Expression const &, SireCAS::Expression const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("atom0"), bp::arg("atom1"), bp::arg("atom2"), bp::arg("atom3"), bp::arg("initial_form"), bp::arg("final_form"), bp::arg("mapping_function"), bp::arg("map")=SireBase::PropertyMap() )) );
        FourAtomPerturbation_exposer.def( bp::init< SireMol::AtomID const &, SireMol::AtomID const &, SireMol::AtomID const &, SireMol::AtomID const &, SireCAS::Expression const &, SireCAS::Identities const &, SireCAS::Identities const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("atom0"), bp::arg("atom1"), bp::arg("atom2"), bp::arg("atom3"), bp::arg("base_expression"), bp::arg("initial_forms"), bp::arg("final_forms"), bp::arg("map")=SireBase::PropertyMap() )) );
        FourAtomPerturbation_exposer.def( bp::init< SireMol::AtomID const &, SireMol::AtomID const &, SireMol::AtomID const &, SireMol::AtomID const &, SireCAS::Expression const &, SireCAS::Identities const &, SireCAS::Identities const &, SireCAS::Expression const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("atom0"), bp::arg("atom1"), bp::arg("atom2"), bp::arg("atom3"), bp::arg("base_expression"), bp::arg("initial_forms"), bp::arg("final_forms"), bp::arg("mapping_function"), bp::arg("map")=SireBase::PropertyMap() )) );
        FourAtomPerturbation_exposer.def( bp::init< SireMM::FourAtomPerturbation const & >(( bp::arg("other") )) );
        { //::SireMM::FourAtomPerturbation::atom0
        
            typedef ::SireMol::AtomID const & ( ::SireMM::FourAtomPerturbation::*atom0_function_type )(  ) const;
            atom0_function_type atom0_function_value( &::SireMM::FourAtomPerturbation::atom0 );
            
            FourAtomPerturbation_exposer.def( 
                "atom0"
                , atom0_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireMM::FourAtomPerturbation::atom1
        
            typedef ::SireMol::AtomID const & ( ::SireMM::FourAtomPerturbation::*atom1_function_type )(  ) const;
            atom1_function_type atom1_function_value( &::SireMM::FourAtomPerturbation::atom1 );
            
            FourAtomPerturbation_exposer.def( 
                "atom1"
                , atom1_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireMM::FourAtomPerturbation::atom2
        
            typedef ::SireMol::AtomID const & ( ::SireMM::FourAtomPerturbation::*atom2_function_type )(  ) const;
            atom2_function_type atom2_function_value( &::SireMM::FourAtomPerturbation::atom2 );
            
            FourAtomPerturbation_exposer.def( 
                "atom2"
                , atom2_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireMM::FourAtomPerturbation::atom3
        
            typedef ::SireMol::AtomID const & ( ::SireMM::FourAtomPerturbation::*atom3_function_type )(  ) const;
            atom3_function_type atom3_function_value( &::SireMM::FourAtomPerturbation::atom3 );
            
            FourAtomPerturbation_exposer.def( 
                "atom3"
                , atom3_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        FourAtomPerturbation_exposer.def( bp::self != bp::self );
        { //::SireMM::FourAtomPerturbation::operator=
        
            typedef ::SireMM::FourAtomPerturbation & ( ::SireMM::FourAtomPerturbation::*assign_function_type )( ::SireMM::FourAtomPerturbation const & ) ;
            assign_function_type assign_function_value( &::SireMM::FourAtomPerturbation::operator= );
            
            FourAtomPerturbation_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        FourAtomPerturbation_exposer.def( bp::self == bp::self );
        { //::SireMM::FourAtomPerturbation::requiredProperties
        
            typedef ::QSet< QString > ( ::SireMM::FourAtomPerturbation::*requiredProperties_function_type )(  ) const;
            requiredProperties_function_type requiredProperties_function_value( &::SireMM::FourAtomPerturbation::requiredProperties );
            
            FourAtomPerturbation_exposer.def( 
                "requiredProperties"
                , requiredProperties_function_value );
        
        }
        { //::SireMM::FourAtomPerturbation::toString
        
            typedef ::QString ( ::SireMM::FourAtomPerturbation::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireMM::FourAtomPerturbation::toString );
            
            FourAtomPerturbation_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireMM::FourAtomPerturbation::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMM::FourAtomPerturbation::typeName );
            
            FourAtomPerturbation_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireMM::FourAtomPerturbation::wouldChange
        
            typedef bool ( ::SireMM::FourAtomPerturbation::*wouldChange_function_type )( ::SireMol::Molecule const &,::SireCAS::Values const & ) const;
            wouldChange_function_type wouldChange_function_value( &::SireMM::FourAtomPerturbation::wouldChange );
            
            FourAtomPerturbation_exposer.def( 
                "wouldChange"
                , wouldChange_function_value
                , ( bp::arg("molecule"), bp::arg("values") ) );
        
        }
        FourAtomPerturbation_exposer.staticmethod( "typeName" );
        FourAtomPerturbation_exposer.def( "__copy__", &__copy__);
        FourAtomPerturbation_exposer.def( "__deepcopy__", &__copy__);
        FourAtomPerturbation_exposer.def( "clone", &__copy__);
        FourAtomPerturbation_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::FourAtomPerturbation >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        FourAtomPerturbation_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::FourAtomPerturbation >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        FourAtomPerturbation_exposer.def( "__str__", &__str__< ::SireMM::FourAtomPerturbation > );
        FourAtomPerturbation_exposer.def( "__repr__", &__str__< ::SireMM::FourAtomPerturbation > );
    }

}
