// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "InterSoftCLJFF.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/mover.hpp"

#include "SireMol/partialmolecule.h"

#include "intersoftcljff.h"

#include "intersoftcljff.h"

SireFF::Inter2B3DFF<SireMM::SoftCLJPotentialInterface<SireMM::InterSoftCLJPotential> > __copy__(const SireFF::Inter2B3DFF<SireMM::SoftCLJPotentialInterface<SireMM::InterSoftCLJPotential> > &other){ return SireFF::Inter2B3DFF<SireMM::SoftCLJPotentialInterface<SireMM::InterSoftCLJPotential> >(other); }

#include "Helpers/str.hpp"

void register_InterSoftCLJFF_class(){

    { //::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >
        typedef bp::class_< SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >, bp::bases< SireFF::FF3D, SireFF::Inter2BFF<SireMM::SoftCLJPotentialInterface<SireMM::InterSoftCLJPotential> >, SireMM::SoftCLJPotentialInterface<SireMM::InterSoftCLJPotential>, SireMM::CLJPotentialInterface<SireMM::InterSoftCLJPotential>, SireFF::G1FF, SireFF::FF, SireMol::MolGroupsBase, SireBase::Property > > InterSoftCLJFF_exposer_t;
        InterSoftCLJFF_exposer_t InterSoftCLJFF_exposer = InterSoftCLJFF_exposer_t( "InterSoftCLJFF" );
        bp::scope InterSoftCLJFF_scope( InterSoftCLJFF_exposer );
        InterSoftCLJFF_exposer.def( bp::init< >() );
        InterSoftCLJFF_exposer.def( bp::init< QString const & >(( bp::arg("name") )) );
        InterSoftCLJFF_exposer.def( bp::init< SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > const & >(( bp::arg("other") )) );
        { //::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::energy
        
            typedef SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::*energy_function_type )( ::SireFF::EnergyTable &,double ) ;
            energy_function_type energy_function_value( &::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::energy );
            
            InterSoftCLJFF_exposer.def( 
                "energy"
                , energy_function_value
                , ( bp::arg("energytable"), bp::arg("scale_energy")=1 ) );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::energy
        
            typedef SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::*energy_function_type )( ::SireFF::EnergyTable &,::SireCAS::Symbol const &,double ) ;
            energy_function_type energy_function_value( &::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::energy );
            
            InterSoftCLJFF_exposer.def( 
                "energy"
                , energy_function_value
                , ( bp::arg("energytable"), bp::arg("symbol"), bp::arg("scale_energy")=1 ) );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::field
        
            typedef SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::*field_function_type )( ::SireFF::FieldTable &,double ) ;
            field_function_type field_function_value( &::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::field );
            
            InterSoftCLJFF_exposer.def( 
                "field"
                , field_function_value
                , ( bp::arg("fieldtable"), bp::arg("scale_field")=1 ) );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::field
        
            typedef SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::*field_function_type )( ::SireFF::FieldTable &,::SireCAS::Symbol const &,double ) ;
            field_function_type field_function_value( &::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::field );
            
            InterSoftCLJFF_exposer.def( 
                "field"
                , field_function_value
                , ( bp::arg("fieldtable"), bp::arg("component"), bp::arg("scale_field")=1 ) );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::field
        
            typedef SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::*field_function_type )( ::SireFF::FieldTable &,::SireFF::Probe const &,double ) ;
            field_function_type field_function_value( &::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::field );
            
            InterSoftCLJFF_exposer.def( 
                "field"
                , field_function_value
                , ( bp::arg("fieldtable"), bp::arg("probe"), bp::arg("scale_field")=1 ) );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::field
        
            typedef SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::*field_function_type )( ::SireFF::FieldTable &,::SireCAS::Symbol const &,::SireFF::Probe const &,double ) ;
            field_function_type field_function_value( &::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::field );
            
            InterSoftCLJFF_exposer.def( 
                "field"
                , field_function_value
                , ( bp::arg("fieldtable"), bp::arg("component"), bp::arg("probe"), bp::arg("scale_field")=1 ) );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::force
        
            typedef SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::*force_function_type )( ::SireFF::ForceTable &,double ) ;
            force_function_type force_function_value( &::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::force );
            
            InterSoftCLJFF_exposer.def( 
                "force"
                , force_function_value
                , ( bp::arg("forcetable"), bp::arg("scale_force")=1 ) );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::force
        
            typedef SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::*force_function_type )( ::SireFF::ForceTable &,::SireCAS::Symbol const &,double ) ;
            force_function_type force_function_value( &::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::force );
            
            InterSoftCLJFF_exposer.def( 
                "force"
                , force_function_value
                , ( bp::arg("forcetable"), bp::arg("symbol"), bp::arg("scale_force")=1 ) );
        
        }
        InterSoftCLJFF_exposer.def( bp::self != bp::self );
        { //::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::operator=
        
            typedef SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > exported_class_t;
            typedef ::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > & ( ::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::*assign_function_type )( ::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > const & ) ;
            assign_function_type assign_function_value( &::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::operator= );
            
            InterSoftCLJFF_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        InterSoftCLJFF_exposer.def( bp::self == bp::self );
        { //::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::packCoordinates
        
            typedef SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::*packCoordinates_function_type )(  ) ;
            packCoordinates_function_type packCoordinates_function_value( &::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::packCoordinates );
            
            InterSoftCLJFF_exposer.def( 
                "packCoordinates"
                , packCoordinates_function_value );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::potential
        
            typedef SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::*potential_function_type )( ::SireFF::PotentialTable &,double ) ;
            potential_function_type potential_function_value( &::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::potential );
            
            InterSoftCLJFF_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("potentialtable"), bp::arg("scale_potential")=1 ) );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::potential
        
            typedef SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::*potential_function_type )( ::SireFF::PotentialTable &,::SireCAS::Symbol const &,double ) ;
            potential_function_type potential_function_value( &::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::potential );
            
            InterSoftCLJFF_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("potentialtable"), bp::arg("component"), bp::arg("scale_potential")=1 ) );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::potential
        
            typedef SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::*potential_function_type )( ::SireFF::PotentialTable &,::SireFF::Probe const &,double ) ;
            potential_function_type potential_function_value( &::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::potential );
            
            InterSoftCLJFF_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("potentialtable"), bp::arg("probe"), bp::arg("scale_potential")=1 ) );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::potential
        
            typedef SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::*potential_function_type )( ::SireFF::PotentialTable &,::SireCAS::Symbol const &,::SireFF::Probe const &,double ) ;
            potential_function_type potential_function_value( &::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::potential );
            
            InterSoftCLJFF_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("potentialtable"), bp::arg("component"), bp::arg("probe"), bp::arg("scale_potential")=1 ) );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::typeName
        
            typedef SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::typeName );
            
            InterSoftCLJFF_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::what
        
            typedef SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > exported_class_t;
            typedef char const * ( ::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::*what_function_type )(  ) const;
            what_function_type what_function_value( &::SireFF::Inter2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::what );
            
            InterSoftCLJFF_exposer.def( 
                "what"
                , what_function_value );
        
        }
        InterSoftCLJFF_exposer.staticmethod( "typeName" );
        InterSoftCLJFF_exposer.def( "__copy__", &__copy__);
        InterSoftCLJFF_exposer.def( "__deepcopy__", &__copy__);
        InterSoftCLJFF_exposer.def( "clone", &__copy__);
        InterSoftCLJFF_exposer.def( "__str__", &__str__< ::SireFF::Inter2B3DFF<SireMM::SoftCLJPotentialInterface<SireMM::InterSoftCLJPotential> > > );
        InterSoftCLJFF_exposer.def( "__repr__", &__str__< ::SireFF::Inter2B3DFF<SireMM::SoftCLJPotentialInterface<SireMM::InterSoftCLJPotential> > > );
    }

}
