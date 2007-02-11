// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "_Molecule__memfuns0.pypp.hpp"
#include "boost/python.hpp"
#include "siremol_headers.h"
#include "SireMaths/angle.h"
#include "SireMaths/quaternion.h"
#include "SireMaths/matrix.h"
#include "SireMaths/triangle.h"
#include "SireMaths/line.h"
#include "SireMaths/torsion.h"

namespace bp = boost::python;

void register_Molecule_memfuns0( Molecule_exposer_t& Molecule_exposer ){

    { //::SireMol::Molecule::ID
    
        typedef ::SireMol::MoleculeID ( ::SireMol::Molecule::*ID_function_type )(  ) const;
        
        Molecule_exposer.def( 
            "ID"
            , ID_function_type( &::SireMol::Molecule::ID ) );
    
    }

    { //::SireMol::Molecule::addProperty
    
        typedef void ( ::SireMol::Molecule::*addProperty_function_type )( ::QString const &,::SireMol::PropertyBase const & ) ;
        
        Molecule_exposer.def( 
            "addProperty"
            , addProperty_function_type( &::SireMol::Molecule::addProperty )
            , ( bp::arg("name"), bp::arg("value") ) );
    
    }

    { //::SireMol::Molecule::addProperty
    
        typedef void ( ::SireMol::Molecule::*addProperty_function_type )( ::QString const &,::SireMol::Property const & ) ;
        
        Molecule_exposer.def( 
            "addProperty"
            , addProperty_function_type( &::SireMol::Molecule::addProperty )
            , ( bp::arg("name"), bp::arg("value") ) );
    
    }

    { //::SireMol::Molecule::addProperty
    
        typedef void ( ::SireMol::Molecule::*addProperty_function_type )( ::QString const &,::QVariant const & ) ;
        
        Molecule_exposer.def( 
            "addProperty"
            , addProperty_function_type( &::SireMol::Molecule::addProperty )
            , ( bp::arg("name"), bp::arg("value") ) );
    
    }

    { //::SireMol::Molecule::angle
    
        typedef ::SireMaths::Triangle ( ::SireMol::Molecule::*angle_function_type )( ::SireMol::Angle const & ) const;
        
        Molecule_exposer.def( 
            "angle"
            , angle_function_type( &::SireMol::Molecule::angle )
            , ( bp::arg("ang") ) );
    
    }

    { //::SireMol::Molecule::assertSameMajorVersion
    
        typedef void ( ::SireMol::Molecule::*assertSameMajorVersion_function_type )( ::SireMol::Molecule const & ) const;
        
        Molecule_exposer.def( 
            "assertSameMajorVersion"
            , assertSameMajorVersion_function_type( &::SireMol::Molecule::assertSameMajorVersion )
            , ( bp::arg("other") ) );
    
    }

    { //::SireMol::Molecule::assertSameMolecule
    
        typedef void ( ::SireMol::Molecule::*assertSameMolecule_function_type )( ::SireMol::Molecule const & ) const;
        
        Molecule_exposer.def( 
            "assertSameMolecule"
            , assertSameMolecule_function_type( &::SireMol::Molecule::assertSameMolecule )
            , ( bp::arg("other") ) );
    
    }

    { //::SireMol::Molecule::assertSameVersion
    
        typedef void ( ::SireMol::Molecule::*assertSameVersion_function_type )( ::SireMol::Molecule const & ) const;
        
        Molecule_exposer.def( 
            "assertSameVersion"
            , assertSameVersion_function_type( &::SireMol::Molecule::assertSameVersion )
            , ( bp::arg("other") ) );
    
    }

    { //::SireMol::Molecule::at
    
        typedef ::SireMol::Residue ( ::SireMol::Molecule::*at_function_type )( ::SireMol::ResNum ) const;
        
        Molecule_exposer.def( 
            "at"
            , at_function_type( &::SireMol::Molecule::at )
            , ( bp::arg("resnum") ) );
    
    }

    { //::SireMol::Molecule::at
    
        typedef ::SireMol::Residue ( ::SireMol::Molecule::*at_function_type )( ::SireMol::ResID ) const;
        
        Molecule_exposer.def( 
            "at"
            , at_function_type( &::SireMol::Molecule::at )
            , ( bp::arg("i") ) );
    
    }

    { //::SireMol::Molecule::at
    
        typedef ::SireMol::CutGroup ( ::SireMol::Molecule::*at_function_type )( ::SireMol::CutGroupID ) const;
        
        Molecule_exposer.def( 
            "at"
            , at_function_type( &::SireMol::Molecule::at )
            , ( bp::arg("cgid") ) );
    
    }

    { //::SireMol::Molecule::at
    
        typedef ::SireMol::Atom ( ::SireMol::Molecule::*at_function_type )( ::SireMol::CGAtomID const & ) const;
        
        Molecule_exposer.def( 
            "at"
            , at_function_type( &::SireMol::Molecule::at )
            , ( bp::arg("cgatomid") ) );
    
    }

    { //::SireMol::Molecule::at
    
        typedef ::SireMol::Atom ( ::SireMol::Molecule::*at_function_type )( ::SireMol::IDMolAtom const & ) const;
        
        Molecule_exposer.def( 
            "at"
            , at_function_type( &::SireMol::Molecule::at )
            , ( bp::arg("atomid") ) );
    
    }

    { //::SireMol::Molecule::atom
    
        typedef ::SireMol::Atom ( ::SireMol::Molecule::*atom_function_type )( ::SireMol::CGAtomID const & ) const;
        
        Molecule_exposer.def( 
            "atom"
            , atom_function_type( &::SireMol::Molecule::atom )
            , ( bp::arg("cgatmid") ) );
    
    }

    { //::SireMol::Molecule::atom
    
        typedef ::SireMol::Atom ( ::SireMol::Molecule::*atom_function_type )( ::SireMol::IDMolAtom const & ) const;
        
        Molecule_exposer.def( 
            "atom"
            , atom_function_type( &::SireMol::Molecule::atom )
            , ( bp::arg("atomid") ) );
    
    }

    { //::SireMol::Molecule::atomNames
    
        typedef ::QStringList ( ::SireMol::Molecule::*atomNames_function_type )( ::SireMol::ResNum ) const;
        
        Molecule_exposer.def( 
            "atomNames"
            , atomNames_function_type( &::SireMol::Molecule::atomNames )
            , ( bp::arg("resnum") ) );
    
    }

    { //::SireMol::Molecule::atomNames
    
        typedef ::QStringList ( ::SireMol::Molecule::*atomNames_function_type )( ::SireMol::ResID ) const;
        
        Molecule_exposer.def( 
            "atomNames"
            , atomNames_function_type( &::SireMol::Molecule::atomNames )
            , ( bp::arg("resid") ) );
    
    }

    { //::SireMol::Molecule::atoms
    
        typedef ::QVector<SireMol::Atom> ( ::SireMol::Molecule::*atoms_function_type )(  ) const;
        
        Molecule_exposer.def( 
            "atoms"
            , atoms_function_type( &::SireMol::Molecule::atoms ) );
    
    }

    { //::SireMol::Molecule::atoms
    
        typedef ::QHash<SireMol::AtomID,SireMol::Atom> ( ::SireMol::Molecule::*atoms_function_type )( ::QSet<SireMol::AtomID> const & ) const;
        
        Molecule_exposer.def( 
            "atoms"
            , atoms_function_type( &::SireMol::Molecule::atoms )
            , ( bp::arg("atomids") ) );
    
    }

    { //::SireMol::Molecule::atoms
    
        typedef ::QHash<SireMol::CGAtomID,SireMol::Atom> ( ::SireMol::Molecule::*atoms_function_type )( ::QSet<SireMol::CGAtomID> const & ) const;
        
        Molecule_exposer.def( 
            "atoms"
            , atoms_function_type( &::SireMol::Molecule::atoms )
            , ( bp::arg("cgatomids") ) );
    
    }

}
