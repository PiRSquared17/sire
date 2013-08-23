// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "PDB.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/stringmangler.h"

#include "SireError/errors.h"

#include "SireIO/errors.h"

#include "SireMol/atomcharges.h"

#include "SireMol/atomcoords.h"

#include "SireMol/atomeditor.h"

#include "SireMol/atomelements.h"

#include "SireMol/cgeditor.h"

#include "SireMol/chain.h"

#include "SireMol/chaineditor.h"

#include "SireMol/connectivity.h"

#include "SireMol/cuttingfunction.h"

#include "SireMol/element.h"

#include "SireMol/molecule.h"

#include "SireMol/moleculegroup.h"

#include "SireMol/moleditor.h"

#include "SireMol/molidentifier.h"

#include "SireMol/molidx.h"

#include "SireMol/mover.hpp"

#include "SireMol/reseditor.h"

#include "SireMol/residue.h"

#include "SireMol/segeditor.h"

#include "SireMol/segment.h"

#include "SireMol/selector.hpp"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireStream/sharestrings.h"

#include "SireUnits/units.h"

#include "pdb.h"

#include "tostring.h"

#include <QByteArray>

#include <QDebug>

#include <QMap>

#include <QString>

#include <QStringList>

#include <QTextStream>

#include "pdb.h"

SireIO::PDB __copy__(const SireIO::PDB &other){ return SireIO::PDB(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_PDB_class(){

    { //::SireIO::PDB
        typedef bp::class_< SireIO::PDB, bp::bases< SireIO::IOBase, SireBase::Property > > PDB_exposer_t;
        PDB_exposer_t PDB_exposer = PDB_exposer_t( "PDB", bp::init< >() );
        bp::scope PDB_scope( PDB_exposer );
        PDB_exposer.def( bp::init< SireIO::PDB const & >(( bp::arg("other") )) );
        PDB_exposer.def( bp::self != bp::self );
        { //::SireIO::PDB::operator=
        
            typedef ::SireIO::PDB & ( ::SireIO::PDB::*assign_function_type )( ::SireIO::PDB const & ) ;
            assign_function_type assign_function_value( &::SireIO::PDB::operator= );
            
            PDB_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        PDB_exposer.def( bp::self == bp::self );
        { //::SireIO::PDB::parameters
        
            typedef ::SireIO::PDBParameters const & ( *parameters_function_type )(  );
            parameters_function_type parameters_function_value( &::SireIO::PDB::parameters );
            
            PDB_exposer.def( 
                "parameters"
                , parameters_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireIO::PDB::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireIO::PDB::typeName );
            
            PDB_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireIO::PDB::what
        
            typedef char const * ( ::SireIO::PDB::*what_function_type )(  ) const;
            what_function_type what_function_value( &::SireIO::PDB::what );
            
            PDB_exposer.def( 
                "what"
                , what_function_value );
        
        }
        PDB_exposer.staticmethod( "parameters" );
        PDB_exposer.staticmethod( "typeName" );
        PDB_exposer.def( "__copy__", &__copy__);
        PDB_exposer.def( "__deepcopy__", &__copy__);
        PDB_exposer.def( "clone", &__copy__);
        PDB_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireIO::PDB >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        PDB_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireIO::PDB >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        PDB_exposer.def( "__str__", &__str__< ::SireIO::PDB > );
        PDB_exposer.def( "__repr__", &__str__< ::SireIO::PDB > );
    }

}
