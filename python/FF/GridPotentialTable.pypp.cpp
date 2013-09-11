// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "GridPotentialTable.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireID/index.h"

#include "SireMol/atomselection.h"

#include "SireMol/errors.h"

#include "SireMol/moleculegroup.h"

#include "SireMol/moleculeview.h"

#include "SireMol/mover.hpp"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "potentialtable.h"

#include "potentialtable.h"

SireFF::GridPotentialTable __copy__(const SireFF::GridPotentialTable &other){ return SireFF::GridPotentialTable(other); }

#include "Qt/qdatastream.hpp"

const char* pvt_get_name(const SireFF::GridPotentialTable&){ return "SireFF::GridPotentialTable";}

void register_GridPotentialTable_class(){

    { //::SireFF::GridPotentialTable
        typedef bp::class_< SireFF::GridPotentialTable > GridPotentialTable_exposer_t;
        GridPotentialTable_exposer_t GridPotentialTable_exposer = GridPotentialTable_exposer_t( "GridPotentialTable", bp::init< >() );
        bp::scope GridPotentialTable_scope( GridPotentialTable_exposer );
        GridPotentialTable_exposer.def( bp::init< SireVol::Grid const & >(( bp::arg("grid") )) );
        GridPotentialTable_exposer.def( bp::init< SireFF::GridPotentialTable const & >(( bp::arg("other") )) );
        { //::SireFF::GridPotentialTable::add
        
            typedef void ( ::SireFF::GridPotentialTable::*add_function_type )( int,::SireUnits::Dimension::MolarEnergy const & ) ;
            add_function_type add_function_value( &::SireFF::GridPotentialTable::add );
            
            GridPotentialTable_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("ipoint"), bp::arg("potential") ) );
        
        }
        { //::SireFF::GridPotentialTable::add
        
            typedef void ( ::SireFF::GridPotentialTable::*add_function_type )( ::SireFF::GridPotentialTable const & ) ;
            add_function_type add_function_value( &::SireFF::GridPotentialTable::add );
            
            GridPotentialTable_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("other") ) );
        
        }
        { //::SireFF::GridPotentialTable::add
        
            typedef void ( ::SireFF::GridPotentialTable::*add_function_type )( ::SireUnits::Dimension::MolarEnergy const & ) ;
            add_function_type add_function_value( &::SireFF::GridPotentialTable::add );
            
            GridPotentialTable_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("potential") ) );
        
        }
        { //::SireFF::GridPotentialTable::at
        
            typedef ::SireUnits::Dimension::MolarEnergy const & ( ::SireFF::GridPotentialTable::*at_function_type )( int ) const;
            at_function_type at_function_value( &::SireFF::GridPotentialTable::at );
            
            GridPotentialTable_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireFF::GridPotentialTable::count
        
            typedef int ( ::SireFF::GridPotentialTable::*count_function_type )(  ) const;
            count_function_type count_function_value( &::SireFF::GridPotentialTable::count );
            
            GridPotentialTable_exposer.def( 
                "count"
                , count_function_value );
        
        }
        { //::SireFF::GridPotentialTable::divide
        
            typedef void ( ::SireFF::GridPotentialTable::*divide_function_type )( double ) ;
            divide_function_type divide_function_value( &::SireFF::GridPotentialTable::divide );
            
            GridPotentialTable_exposer.def( 
                "divide"
                , divide_function_value
                , ( bp::arg("value") ) );
        
        }
        { //::SireFF::GridPotentialTable::grid
        
            typedef ::SireVol::Grid const & ( ::SireFF::GridPotentialTable::*grid_function_type )(  ) const;
            grid_function_type grid_function_value( &::SireFF::GridPotentialTable::grid );
            
            GridPotentialTable_exposer.def( 
                "grid"
                , grid_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireFF::GridPotentialTable::initialise
        
            typedef void ( ::SireFF::GridPotentialTable::*initialise_function_type )(  ) ;
            initialise_function_type initialise_function_value( &::SireFF::GridPotentialTable::initialise );
            
            GridPotentialTable_exposer.def( 
                "initialise"
                , initialise_function_value );
        
        }
        { //::SireFF::GridPotentialTable::multiply
        
            typedef void ( ::SireFF::GridPotentialTable::*multiply_function_type )( double ) ;
            multiply_function_type multiply_function_value( &::SireFF::GridPotentialTable::multiply );
            
            GridPotentialTable_exposer.def( 
                "multiply"
                , multiply_function_value
                , ( bp::arg("value") ) );
        
        }
        { //::SireFF::GridPotentialTable::nPoints
        
            typedef int ( ::SireFF::GridPotentialTable::*nPoints_function_type )(  ) const;
            nPoints_function_type nPoints_function_value( &::SireFF::GridPotentialTable::nPoints );
            
            GridPotentialTable_exposer.def( 
                "nPoints"
                , nPoints_function_value );
        
        }
        GridPotentialTable_exposer.def( bp::self != bp::self );
        GridPotentialTable_exposer.def( bp::self * bp::other< double >() );
        GridPotentialTable_exposer.def( bp::self + bp::self );
        GridPotentialTable_exposer.def( bp::self + bp::other< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >() );
        GridPotentialTable_exposer.def( bp::self - bp::self );
        GridPotentialTable_exposer.def( bp::self - bp::other< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >() );
        GridPotentialTable_exposer.def( -bp::self );
        GridPotentialTable_exposer.def( bp::self / bp::other< double >() );
        { //::SireFF::GridPotentialTable::operator=
        
            typedef ::SireFF::GridPotentialTable & ( ::SireFF::GridPotentialTable::*assign_function_type )( ::SireFF::GridPotentialTable const & ) ;
            assign_function_type assign_function_value( &::SireFF::GridPotentialTable::operator= );
            
            GridPotentialTable_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        { //::SireFF::GridPotentialTable::operator=
        
            typedef ::SireFF::GridPotentialTable & ( ::SireFF::GridPotentialTable::*assign_function_type )( ::SireUnits::Dimension::MolarEnergy const & ) ;
            assign_function_type assign_function_value( &::SireFF::GridPotentialTable::operator= );
            
            GridPotentialTable_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("potential") )
                , bp::return_self< >() );
        
        }
        GridPotentialTable_exposer.def( bp::self == bp::self );
        { //::SireFF::GridPotentialTable::operator[]
        
            typedef ::SireUnits::Dimension::MolarEnergy & ( ::SireFF::GridPotentialTable::*__getitem___function_type )( int ) ;
            __getitem___function_type __getitem___function_value( &::SireFF::GridPotentialTable::operator[] );
            
            GridPotentialTable_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , bp::return_internal_reference< >() );
        
        }
        { //::SireFF::GridPotentialTable::operator[]
        
            typedef ::SireUnits::Dimension::MolarEnergy const & ( ::SireFF::GridPotentialTable::*__getitem___function_type )( int ) const;
            __getitem___function_type __getitem___function_value( &::SireFF::GridPotentialTable::operator[] );
            
            GridPotentialTable_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireFF::GridPotentialTable::setAll
        
            typedef void ( ::SireFF::GridPotentialTable::*setAll_function_type )( ::SireUnits::Dimension::MolarEnergy const & ) ;
            setAll_function_type setAll_function_value( &::SireFF::GridPotentialTable::setAll );
            
            GridPotentialTable_exposer.def( 
                "setAll"
                , setAll_function_value
                , ( bp::arg("potential") ) );
        
        }
        { //::SireFF::GridPotentialTable::subtract
        
            typedef void ( ::SireFF::GridPotentialTable::*subtract_function_type )( int,::SireUnits::Dimension::MolarEnergy const & ) ;
            subtract_function_type subtract_function_value( &::SireFF::GridPotentialTable::subtract );
            
            GridPotentialTable_exposer.def( 
                "subtract"
                , subtract_function_value
                , ( bp::arg("ipoint"), bp::arg("potential") ) );
        
        }
        { //::SireFF::GridPotentialTable::subtract
        
            typedef void ( ::SireFF::GridPotentialTable::*subtract_function_type )( ::SireFF::GridPotentialTable const & ) ;
            subtract_function_type subtract_function_value( &::SireFF::GridPotentialTable::subtract );
            
            GridPotentialTable_exposer.def( 
                "subtract"
                , subtract_function_value
                , ( bp::arg("other") ) );
        
        }
        { //::SireFF::GridPotentialTable::subtract
        
            typedef void ( ::SireFF::GridPotentialTable::*subtract_function_type )( ::SireUnits::Dimension::MolarEnergy const & ) ;
            subtract_function_type subtract_function_value( &::SireFF::GridPotentialTable::subtract );
            
            GridPotentialTable_exposer.def( 
                "subtract"
                , subtract_function_value
                , ( bp::arg("potential") ) );
        
        }
        { //::SireFF::GridPotentialTable::toVector
        
            typedef ::QVector< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > > ( ::SireFF::GridPotentialTable::*toVector_function_type )(  ) const;
            toVector_function_type toVector_function_value( &::SireFF::GridPotentialTable::toVector );
            
            GridPotentialTable_exposer.def( 
                "toVector"
                , toVector_function_value );
        
        }
        { //::SireFF::GridPotentialTable::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireFF::GridPotentialTable::typeName );
            
            GridPotentialTable_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        GridPotentialTable_exposer.staticmethod( "typeName" );
        GridPotentialTable_exposer.def( "__copy__", &__copy__);
        GridPotentialTable_exposer.def( "__deepcopy__", &__copy__);
        GridPotentialTable_exposer.def( "clone", &__copy__);
        GridPotentialTable_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireFF::GridPotentialTable >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        GridPotentialTable_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireFF::GridPotentialTable >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        GridPotentialTable_exposer.def( "__str__", &pvt_get_name);
        GridPotentialTable_exposer.def( "__repr__", &pvt_get_name);
    }

}
