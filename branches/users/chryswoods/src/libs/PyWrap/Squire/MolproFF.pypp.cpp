// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "MolproFF.pypp.hpp"
#include "boost/python.hpp"
#include "squire_headers.h"
#include "SireMol/molecule.h"
#include "SireMol/residue.h"
#include "SireMol/newatom.h"

namespace bp = boost::python;

#include "SireQt/qdatastream.hpp"

#include "SirePy/str.hpp"

void register_MolproFF_class(){

    { //::Squire::MolproFF
        typedef bp::class_< Squire::MolproFF, bp::bases< SireFF::FFBase > > MolproFF_exposer_t;
        MolproFF_exposer_t MolproFF_exposer = MolproFF_exposer_t( "MolproFF" );
        bp::scope MolproFF_scope( MolproFF_exposer );
        bp::class_< Squire::MolproFF::Components, bp::bases< SireFF::FFBase::Components > >( "Components" )    
            .def( bp::init< >() )    
            .def( bp::init< SireFF::FFBase const &, bp::optional< SireCAS::Symbols const & > >(( bp::arg("ffbase"), bp::arg("symbols")=::SireCAS::Symbols( ) )) )    
            .def( 
                "describe_qm"
                , (::QString (*)(  ))( &::Squire::MolproFF::Components::describe_qm ) )    
            .def( 
                "qm"
                , (::SireFF::FFComponent const & ( ::Squire::MolproFF::Components::* )(  ) const)( &::Squire::MolproFF::Components::qm )
                , bp::return_value_policy< bp::copy_const_reference >() )    
            .staticmethod( "describe_qm" );
        bp::class_< Squire::MolproFF::Parameters, bp::bases< SireFF::FFBase::Parameters > >( "Parameters" )    
            .def( bp::init< >() )    
            .def( 
                "coulomb"
                , (::SireFF::ParameterName const & ( ::Squire::MolproFF::Parameters::* )(  ) const)( &::Squire::MolproFF::Parameters::coulomb )
                , bp::return_value_policy< bp::copy_const_reference >() )    
            .def_readonly( "default_sources", Squire::MolproFF::Parameters::default_sources );
        MolproFF_exposer.def( bp::init< >() );
        MolproFF_exposer.def( bp::init< SireVol::Space const &, SireMM::SwitchingFunction const & >(( bp::arg("space"), bp::arg("switchfunc") )) );
        { //::Squire::MolproFF::addToMM
        
            typedef bool ( ::Squire::MolproFF::*addToMM_function_type )( ::SireMol::Molecule const &,::SireFF::ParameterMap const & ) ;
            
            MolproFF_exposer.def( 
                "addToMM"
                , addToMM_function_type( &::Squire::MolproFF::addToMM )
                , ( bp::arg("molecule"), bp::arg("map")=::SireFF::ParameterMap( ) ) );
        
        }
        { //::Squire::MolproFF::addToMM
        
            typedef bool ( ::Squire::MolproFF::*addToMM_function_type )( ::QList<SireMol::Molecule> const &,::SireFF::ParameterMap const & ) ;
            
            MolproFF_exposer.def( 
                "addToMM"
                , addToMM_function_type( &::Squire::MolproFF::addToMM )
                , ( bp::arg("molecules"), bp::arg("map")=::SireFF::ParameterMap( ) ) );
        
        }
        { //::Squire::MolproFF::addToQM
        
            typedef bool ( ::Squire::MolproFF::*addToQM_function_type )( ::SireMol::Molecule const & ) ;
            
            MolproFF_exposer.def( 
                "addToQM"
                , addToQM_function_type( &::Squire::MolproFF::addToQM )
                , ( bp::arg("molecule") ) );
        
        }
        { //::Squire::MolproFF::addToQM
        
            typedef bool ( ::Squire::MolproFF::*addToQM_function_type )( ::QList<SireMol::Molecule> const & ) ;
            
            MolproFF_exposer.def( 
                "addToQM"
                , addToQM_function_type( &::Squire::MolproFF::addToQM )
                , ( bp::arg("molecules") ) );
        
        }
        { //::Squire::MolproFF::change
        
            typedef bool ( ::Squire::MolproFF::*change_function_type )( ::SireMol::Molecule const & ) ;
            
            MolproFF_exposer.def( 
                "change"
                , change_function_type( &::Squire::MolproFF::change )
                , ( bp::arg("molecule") ) );
        
        }
        { //::Squire::MolproFF::change
        
            typedef bool ( ::Squire::MolproFF::*change_function_type )( ::SireMol::Residue const & ) ;
            
            MolproFF_exposer.def( 
                "change"
                , change_function_type( &::Squire::MolproFF::change )
                , ( bp::arg("residue") ) );
        
        }
        { //::Squire::MolproFF::change
        
            typedef bool ( ::Squire::MolproFF::*change_function_type )( ::SireMol::NewAtom const & ) ;
            
            MolproFF_exposer.def( 
                "change"
                , change_function_type( &::Squire::MolproFF::change )
                , ( bp::arg("atom") ) );
        
        }
        { //::Squire::MolproFF::components
        
            typedef ::Squire::MolproFF::Components const & ( ::Squire::MolproFF::*components_function_type )(  ) const;
            
            MolproFF_exposer.def( 
                "components"
                , components_function_type( &::Squire::MolproFF::components )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::Squire::MolproFF::mmCoordsAndChargesArray
        
            typedef ::QVector<double> const & ( ::Squire::MolproFF::*mmCoordsAndChargesArray_function_type )(  ) const;
            
            MolproFF_exposer.def( 
                "mmCoordsAndChargesArray"
                , mmCoordsAndChargesArray_function_type( &::Squire::MolproFF::mmCoordsAndChargesArray )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::Squire::MolproFF::molproCommandInput
        
            typedef ::QString ( ::Squire::MolproFF::*molproCommandInput_function_type )(  ) ;
            
            MolproFF_exposer.def( 
                "molproCommandInput"
                , molproCommandInput_function_type( &::Squire::MolproFF::molproCommandInput ) );
        
        }
        { //::Squire::MolproFF::molproExe
        
            typedef ::QFileInfo const & ( ::Squire::MolproFF::*molproExe_function_type )(  ) const;
            
            MolproFF_exposer.def( 
                "molproExe"
                , molproExe_function_type( &::Squire::MolproFF::molproExe )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::Squire::MolproFF::molproTempDir
        
            typedef ::QDir const & ( ::Squire::MolproFF::*molproTempDir_function_type )(  ) const;
            
            MolproFF_exposer.def( 
                "molproTempDir"
                , molproTempDir_function_type( &::Squire::MolproFF::molproTempDir )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::Squire::MolproFF::parameters
        
            typedef ::Squire::MolproFF::Parameters const & ( ::Squire::MolproFF::*parameters_function_type )(  ) const;
            
            MolproFF_exposer.def( 
                "parameters"
                , parameters_function_type( &::Squire::MolproFF::parameters )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::Squire::MolproFF::qmCoordsArray
        
            typedef ::QVector<double> const & ( ::Squire::MolproFF::*qmCoordsArray_function_type )(  ) const;
            
            MolproFF_exposer.def( 
                "qmCoordsArray"
                , qmCoordsArray_function_type( &::Squire::MolproFF::qmCoordsArray )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::Squire::MolproFF::qmVersion
        
            typedef ::quint32 ( ::Squire::MolproFF::*qmVersion_function_type )(  ) const;
            
            MolproFF_exposer.def( 
                "qmVersion"
                , qmVersion_function_type( &::Squire::MolproFF::qmVersion ) );
        
        }
        { //::Squire::MolproFF::remove
        
            typedef bool ( ::Squire::MolproFF::*remove_function_type )( ::SireMol::Molecule const & ) ;
            
            MolproFF_exposer.def( 
                "remove"
                , remove_function_type( &::Squire::MolproFF::remove )
                , ( bp::arg("molecule") ) );
        
        }
        { //::Squire::MolproFF::setEnergyOrigin
        
            typedef void ( ::Squire::MolproFF::*setEnergyOrigin_function_type )( double ) ;
            
            MolproFF_exposer.def( 
                "setEnergyOrigin"
                , setEnergyOrigin_function_type( &::Squire::MolproFF::setEnergyOrigin )
                , ( bp::arg("nrg_zero") ) );
        
        }
        { //::Squire::MolproFF::setMolproExe
        
            typedef void ( ::Squire::MolproFF::*setMolproExe_function_type )( ::QFileInfo const & ) ;
            
            MolproFF_exposer.def( 
                "setMolproExe"
                , setMolproExe_function_type( &::Squire::MolproFF::setMolproExe )
                , ( bp::arg("molpro") ) );
        
        }
        { //::Squire::MolproFF::setMolproTempDir
        
            typedef void ( ::Squire::MolproFF::*setMolproTempDir_function_type )( ::QDir const & ) ;
            
            MolproFF_exposer.def( 
                "setMolproTempDir"
                , setMolproTempDir_function_type( &::Squire::MolproFF::setMolproTempDir )
                , ( bp::arg("tempdir") ) );
        
        }
        { //::Squire::MolproFF::space
        
            typedef ::SireVol::Space const & ( ::Squire::MolproFF::*space_function_type )(  ) const;
            
            MolproFF_exposer.def( 
                "space"
                , space_function_type( &::Squire::MolproFF::space )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::Squire::MolproFF::switchingFunction
        
            typedef ::SireMM::SwitchingFunction const & ( ::Squire::MolproFF::*switchingFunction_function_type )(  ) const;
            
            MolproFF_exposer.def( 
                "switchingFunction"
                , switchingFunction_function_type( &::Squire::MolproFF::switchingFunction )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::Squire::MolproFF::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            
            MolproFF_exposer.def( 
                "typeName"
                , typeName_function_type( &::Squire::MolproFF::typeName ) );
        
        }
        { //::Squire::MolproFF::what
        
            typedef char const * ( ::Squire::MolproFF::*what_function_type )(  ) const;
            
            MolproFF_exposer.def( 
                "what"
                , what_function_type( &::Squire::MolproFF::what ) );
        
        }
        MolproFF_exposer.staticmethod( "typeName" );
        MolproFF_exposer.def( "__rlshift__", &SireQt::__rlshift__QDataStream< ::Squire::MolproFF >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        MolproFF_exposer.def( "__rrshift__", &SireQt::__rrshift__QDataStream< ::Squire::MolproFF >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        MolproFF_exposer.def( "__str__", &SirePy::__str__< ::Squire::MolproFF > );
    }

}
