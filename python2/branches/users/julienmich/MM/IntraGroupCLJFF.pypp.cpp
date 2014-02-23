// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "IntraGroupCLJFF.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/mover.hpp"

#include "SireMol/partialmolecule.h"

#include "intracljff.h"

#include "intracljff.h"

SireFF::Intra2B2G3DFF<SireMM::CLJPotentialInterface<SireMM::IntraCLJPotential> > __copy__(const SireFF::Intra2B2G3DFF<SireMM::CLJPotentialInterface<SireMM::IntraCLJPotential> > &other){ return SireFF::Intra2B2G3DFF<SireMM::CLJPotentialInterface<SireMM::IntraCLJPotential> >(other); }

#include "Helpers/str.hpp"

void register_IntraGroupCLJFF_class(){

    { //::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >
        typedef bp::class_< SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >, bp::bases< SireFF::FF3D, SireFF::Intra2B2GFF<SireMM::CLJPotentialInterface<SireMM::IntraCLJPotential> >, SireMM::CLJPotentialInterface<SireMM::IntraCLJPotential>, SireFF::G2FF, SireFF::FF, SireMol::MolGroupsBase, SireBase::Property > > IntraGroupCLJFF_exposer_t;
        IntraGroupCLJFF_exposer_t IntraGroupCLJFF_exposer = IntraGroupCLJFF_exposer_t( "IntraGroupCLJFF" );
        bp::scope IntraGroupCLJFF_scope( IntraGroupCLJFF_exposer );
        IntraGroupCLJFF_exposer.def( bp::init< >() );
        IntraGroupCLJFF_exposer.def( bp::init< QString const & >(( bp::arg("name") )) );
        IntraGroupCLJFF_exposer.def( bp::init< SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > > const & >(( bp::arg("other") )) );
        { //::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::energy
        
            typedef SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::*energy_function_type )( ::SireFF::EnergyTable &,double ) ;
            energy_function_type energy_function_value( &::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::energy );
            
            IntraGroupCLJFF_exposer.def( 
                "energy"
                , energy_function_value
                , ( bp::arg("energytable"), bp::arg("scale_energy")=1 ) );
        
        }
        { //::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::energy
        
            typedef SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::*energy_function_type )( ::SireFF::EnergyTable &,::SireCAS::Symbol const &,double ) ;
            energy_function_type energy_function_value( &::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::energy );
            
            IntraGroupCLJFF_exposer.def( 
                "energy"
                , energy_function_value
                , ( bp::arg("energytable"), bp::arg("symbol"), bp::arg("scale_energy")=1 ) );
        
        }
        { //::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::field
        
            typedef SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::*field_function_type )( ::SireFF::FieldTable &,double ) ;
            field_function_type field_function_value( &::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::field );
            
            IntraGroupCLJFF_exposer.def( 
                "field"
                , field_function_value
                , ( bp::arg("fieldtable"), bp::arg("scale_field")=1 ) );
        
        }
        { //::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::field
        
            typedef SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::*field_function_type )( ::SireFF::FieldTable &,::SireCAS::Symbol const &,double ) ;
            field_function_type field_function_value( &::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::field );
            
            IntraGroupCLJFF_exposer.def( 
                "field"
                , field_function_value
                , ( bp::arg("fieldtable"), bp::arg("component"), bp::arg("scale_field")=1 ) );
        
        }
        { //::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::field
        
            typedef SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::*field_function_type )( ::SireFF::FieldTable &,::SireFF::Probe const &,double ) ;
            field_function_type field_function_value( &::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::field );
            
            IntraGroupCLJFF_exposer.def( 
                "field"
                , field_function_value
                , ( bp::arg("fieldtable"), bp::arg("probe"), bp::arg("scale_field")=1 ) );
        
        }
        { //::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::field
        
            typedef SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::*field_function_type )( ::SireFF::FieldTable &,::SireCAS::Symbol const &,::SireFF::Probe const &,double ) ;
            field_function_type field_function_value( &::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::field );
            
            IntraGroupCLJFF_exposer.def( 
                "field"
                , field_function_value
                , ( bp::arg("fieldtable"), bp::arg("component"), bp::arg("probe"), bp::arg("scale_field")=1 ) );
        
        }
        { //::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::force
        
            typedef SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::*force_function_type )( ::SireFF::ForceTable &,double ) ;
            force_function_type force_function_value( &::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::force );
            
            IntraGroupCLJFF_exposer.def( 
                "force"
                , force_function_value
                , ( bp::arg("forcetable"), bp::arg("scale_force")=1 ) );
        
        }
        { //::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::force
        
            typedef SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::*force_function_type )( ::SireFF::ForceTable &,::SireCAS::Symbol const &,double ) ;
            force_function_type force_function_value( &::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::force );
            
            IntraGroupCLJFF_exposer.def( 
                "force"
                , force_function_value
                , ( bp::arg("forcetable"), bp::arg("symbol"), bp::arg("scale_force")=1 ) );
        
        }
        IntraGroupCLJFF_exposer.def( bp::self != bp::self );
        { //::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::operator=
        
            typedef SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > > exported_class_t;
            typedef ::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > > & ( ::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::*assign_function_type )( ::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > > const & ) ;
            assign_function_type assign_function_value( &::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::operator= );
            
            IntraGroupCLJFF_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        IntraGroupCLJFF_exposer.def( bp::self == bp::self );
        { //::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::potential
        
            typedef SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::*potential_function_type )( ::SireFF::PotentialTable &,double ) ;
            potential_function_type potential_function_value( &::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::potential );
            
            IntraGroupCLJFF_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("potentialtable"), bp::arg("scale_potential")=1 ) );
        
        }
        { //::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::potential
        
            typedef SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::*potential_function_type )( ::SireFF::PotentialTable &,::SireCAS::Symbol const &,double ) ;
            potential_function_type potential_function_value( &::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::potential );
            
            IntraGroupCLJFF_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("potentialtable"), bp::arg("component"), bp::arg("scale_potential")=1 ) );
        
        }
        { //::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::potential
        
            typedef SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::*potential_function_type )( ::SireFF::PotentialTable &,::SireFF::Probe const &,double ) ;
            potential_function_type potential_function_value( &::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::potential );
            
            IntraGroupCLJFF_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("potentialtable"), bp::arg("probe"), bp::arg("scale_potential")=1 ) );
        
        }
        { //::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::potential
        
            typedef SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::*potential_function_type )( ::SireFF::PotentialTable &,::SireCAS::Symbol const &,::SireFF::Probe const &,double ) ;
            potential_function_type potential_function_value( &::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::potential );
            
            IntraGroupCLJFF_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("potentialtable"), bp::arg("component"), bp::arg("probe"), bp::arg("scale_potential")=1 ) );
        
        }
        { //::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::typeName
        
            typedef SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::typeName );
            
            IntraGroupCLJFF_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::what
        
            typedef SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > > exported_class_t;
            typedef char const * ( ::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::*what_function_type )(  ) const;
            what_function_type what_function_value( &::SireFF::Intra2B2G3DFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::what );
            
            IntraGroupCLJFF_exposer.def( 
                "what"
                , what_function_value );
        
        }
        IntraGroupCLJFF_exposer.staticmethod( "typeName" );
        IntraGroupCLJFF_exposer.def( "__copy__", &__copy__);
        IntraGroupCLJFF_exposer.def( "__deepcopy__", &__copy__);
        IntraGroupCLJFF_exposer.def( "clone", &__copy__);
        IntraGroupCLJFF_exposer.def( "__str__", &__str__< ::SireFF::Intra2B2G3DFF<SireMM::CLJPotentialInterface<SireMM::IntraCLJPotential> > > );
        IntraGroupCLJFF_exposer.def( "__repr__", &__str__< ::SireFF::Intra2B2G3DFF<SireMM::CLJPotentialInterface<SireMM::IntraCLJPotential> > > );
    }

}
