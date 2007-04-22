// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "_Molecule__memfuns1.pypp.hpp"
#include "boost/python.hpp"
#include "siremol_headers.h"
#include "SireMol/moleculedata.h"
#include "SireMaths/angle.h"
#include "SireMaths/quaternion.h"
#include "SireMaths/matrix.h"
#include "SireMaths/triangle.h"
#include "SireMaths/line.h"
#include "SireMaths/torsion.h"
#include "SireVol/space.h"

namespace bp = boost::python;

void register_Molecule_memfuns1( Molecule_exposer_t& Molecule_exposer ){

    { //::SireMol::Molecule::atoms
    
        typedef ::QHash<SireMol::ResNumAtomID,SireMol::Atom> ( ::SireMol::Molecule::*atoms_function_type )( ::QSet<SireMol::ResNumAtomID> const & ) const;
        
        Molecule_exposer.def( 
            "atoms"
            , atoms_function_type( &::SireMol::Molecule::atoms )
            , ( bp::arg("resatomids") ) );
    
    }

    { //::SireMol::Molecule::atoms
    
        typedef ::QHash<SireMol::ResIDAtomID,SireMol::Atom> ( ::SireMol::Molecule::*atoms_function_type )( ::QSet<SireMol::ResIDAtomID> const & ) const;
        
        Molecule_exposer.def( 
            "atoms"
            , atoms_function_type( &::SireMol::Molecule::atoms )
            , ( bp::arg("resatomids") ) );
    
    }

    { //::SireMol::Molecule::atoms
    
        typedef ::QHash<SireMol::AtomIndex,SireMol::Atom> ( ::SireMol::Molecule::*atoms_function_type )( ::QSet<SireMol::AtomIndex> const & ) const;
        
        Molecule_exposer.def( 
            "atoms"
            , atoms_function_type( &::SireMol::Molecule::atoms )
            , ( bp::arg("atoms") ) );
    
    }

    { //::SireMol::Molecule::atoms
    
        typedef ::QVector<SireMol::Atom> ( ::SireMol::Molecule::*atoms_function_type )( ::SireMol::CutGroupID ) const;
        
        Molecule_exposer.def( 
            "atoms"
            , atoms_function_type( &::SireMol::Molecule::atoms )
            , ( bp::arg("cgid") ) );
    
    }

    { //::SireMol::Molecule::atoms
    
        typedef ::QHash<SireMol::CutGroupID,QVector<SireMol::Atom> > ( ::SireMol::Molecule::*atoms_function_type )( ::QSet<SireMol::CutGroupID> const & ) const;
        
        Molecule_exposer.def( 
            "atoms"
            , atoms_function_type( &::SireMol::Molecule::atoms )
            , ( bp::arg("cgids") ) );
    
    }

    { //::SireMol::Molecule::atoms
    
        typedef ::QVector<SireMol::Atom> ( ::SireMol::Molecule::*atoms_function_type )( ::SireMol::ResNum ) const;
        
        Molecule_exposer.def( 
            "atoms"
            , atoms_function_type( &::SireMol::Molecule::atoms )
            , ( bp::arg("resnum") ) );
    
    }

    { //::SireMol::Molecule::atoms
    
        typedef ::QHash<SireMol::ResNum,QVector<SireMol::Atom> > ( ::SireMol::Molecule::*atoms_function_type )( ::QSet<SireMol::ResNum> const & ) const;
        
        Molecule_exposer.def( 
            "atoms"
            , atoms_function_type( &::SireMol::Molecule::atoms )
            , ( bp::arg("resnums") ) );
    
    }

    { //::SireMol::Molecule::atoms
    
        typedef ::QVector<SireMol::Atom> ( ::SireMol::Molecule::*atoms_function_type )( ::SireMol::ResID ) const;
        
        Molecule_exposer.def( 
            "atoms"
            , atoms_function_type( &::SireMol::Molecule::atoms )
            , ( bp::arg("resid") ) );
    
    }

    { //::SireMol::Molecule::atoms
    
        typedef ::QHash<SireMol::ResID,QVector<SireMol::Atom> > ( ::SireMol::Molecule::*atoms_function_type )( ::QSet<SireMol::ResID> const & ) const;
        
        Molecule_exposer.def( 
            "atoms"
            , atoms_function_type( &::SireMol::Molecule::atoms )
            , ( bp::arg("resids") ) );
    
    }

    { //::SireMol::Molecule::bond
    
        typedef ::SireMaths::Line ( ::SireMol::Molecule::*bond_function_type )( ::SireMol::Bond const & ) const;
        
        Molecule_exposer.def( 
            "bond"
            , bond_function_type( &::SireMol::Molecule::bond )
            , ( bp::arg("bnd") ) );
    
    }

    { //::SireMol::Molecule::change
    
        typedef void ( ::SireMol::Molecule::*change_function_type )( ::SireMol::Bond const &,double,::QSet<SireMol::AtomIndex> const & ) ;
        
        Molecule_exposer.def( 
            "change"
            , change_function_type( &::SireMol::Molecule::change )
            , ( bp::arg("bond"), bp::arg("delta"), bp::arg("anchors")=::QSet<SireMol::AtomIndex>( ) ) );
    
    }

    { //::SireMol::Molecule::change
    
        typedef void ( ::SireMol::Molecule::*change_function_type )( ::SireMol::Bond const &,double,::SireMol::WeightFunction const &,::QSet<SireMol::AtomIndex> const & ) ;
        
        Molecule_exposer.def( 
            "change"
            , change_function_type( &::SireMol::Molecule::change )
            , ( bp::arg("bond"), bp::arg("delta"), bp::arg("func"), bp::arg("anchors")=::QSet<SireMol::AtomIndex>( ) ) );
    
    }

    { //::SireMol::Molecule::change
    
        typedef void ( ::SireMol::Molecule::*change_function_type )( ::SireMol::Angle const &,::SireMaths::Angle const &,::QSet<SireMol::AtomIndex> const & ) ;
        
        Molecule_exposer.def( 
            "change"
            , change_function_type( &::SireMol::Molecule::change )
            , ( bp::arg("angle"), bp::arg("delta"), bp::arg("anchors")=::QSet<SireMol::AtomIndex>( ) ) );
    
    }

    { //::SireMol::Molecule::change
    
        typedef void ( ::SireMol::Molecule::*change_function_type )( ::SireMol::Angle const &,::SireMaths::Angle const &,::SireMol::WeightFunction const &,::QSet<SireMol::AtomIndex> const & ) ;
        
        Molecule_exposer.def( 
            "change"
            , change_function_type( &::SireMol::Molecule::change )
            , ( bp::arg("angle"), bp::arg("delta"), bp::arg("func"), bp::arg("anchors")=::QSet<SireMol::AtomIndex>( ) ) );
    
    }

    { //::SireMol::Molecule::change
    
        typedef void ( ::SireMol::Molecule::*change_function_type )( ::SireMol::Dihedral const &,::SireMaths::Angle const &,::QSet<SireMol::AtomIndex> const & ) ;
        
        Molecule_exposer.def( 
            "change"
            , change_function_type( &::SireMol::Molecule::change )
            , ( bp::arg("dihedral"), bp::arg("delta"), bp::arg("anchors")=::QSet<SireMol::AtomIndex>( ) ) );
    
    }

    { //::SireMol::Molecule::change
    
        typedef void ( ::SireMol::Molecule::*change_function_type )( ::SireMol::Dihedral const &,::SireMaths::Angle const &,::SireMol::WeightFunction const &,::QSet<SireMol::AtomIndex> const & ) ;
        
        Molecule_exposer.def( 
            "change"
            , change_function_type( &::SireMol::Molecule::change )
            , ( bp::arg("dihedral"), bp::arg("delta"), bp::arg("func"), bp::arg("anchors")=::QSet<SireMol::AtomIndex>( ) ) );
    
    }

    { //::SireMol::Molecule::change
    
        typedef void ( ::SireMol::Molecule::*change_function_type )( ::SireMol::Bond const &,::SireMaths::Angle const &,::QSet<SireMol::AtomIndex> const & ) ;
        
        Molecule_exposer.def( 
            "change"
            , change_function_type( &::SireMol::Molecule::change )
            , ( bp::arg("bond"), bp::arg("delta"), bp::arg("anchors")=::QSet<SireMol::AtomIndex>( ) ) );
    
    }

    { //::SireMol::Molecule::change
    
        typedef void ( ::SireMol::Molecule::*change_function_type )( ::SireMol::Bond const &,::SireMaths::Angle const &,::SireMol::WeightFunction const &,::QSet<SireMol::AtomIndex> const & ) ;
        
        Molecule_exposer.def( 
            "change"
            , change_function_type( &::SireMol::Molecule::change )
            , ( bp::arg("bond"), bp::arg("delta"), bp::arg("func"), bp::arg("anchors")=::QSet<SireMol::AtomIndex>( ) ) );
    
    }

    { //::SireMol::Molecule::change
    
        typedef void ( ::SireMol::Molecule::*change_function_type )( ::SireMol::Improper const &,::SireMaths::Angle const &,::QSet<SireMol::AtomIndex> const & ) ;
        
        Molecule_exposer.def( 
            "change"
            , change_function_type( &::SireMol::Molecule::change )
            , ( bp::arg("improper"), bp::arg("delta"), bp::arg("anchors")=::QSet<SireMol::AtomIndex>( ) ) );
    
    }

    { //::SireMol::Molecule::change
    
        typedef void ( ::SireMol::Molecule::*change_function_type )( ::SireMol::Improper const &,::SireMaths::Angle const &,::SireMol::WeightFunction const &,::QSet<SireMol::AtomIndex> const & ) ;
        
        Molecule_exposer.def( 
            "change"
            , change_function_type( &::SireMol::Molecule::change )
            , ( bp::arg("improper"), bp::arg("delta"), bp::arg("func"), bp::arg("anchors")=::QSet<SireMol::AtomIndex>( ) ) );
    
    }

}
