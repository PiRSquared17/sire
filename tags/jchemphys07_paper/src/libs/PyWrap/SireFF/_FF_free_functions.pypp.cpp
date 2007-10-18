// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "_FF_free_functions.pypp.hpp"
#include "boost/python.hpp"
#include "sireff_headers.h"
#include "SireMol/molecule.h"
#include "SireMol/molecules.h"
#include "SireMol/partialmolecule.h"
#include "SireMol/residue.h"
#include "SireMol/newatom.h"
#include "SireMol/atom.h"
#include "SireMol/cgatomid.h"
#include "SireMol/cutgroupnum.h"
#include "SireMol/resnumatomid.h"
#include "SireMol/resid.h"
#include "SireMol/moleculeid.h"
#include "SireBase/property.h"

namespace bp = boost::python;

void register_free_functions(){

    { //::SireFF::e_total
    
        typedef ::SireCAS::Symbol const & ( *e_total_function_type )(  );
        
        bp::def( 
            "e_total"
            , e_total_function_type( &::SireFF::e_total )
            , bp::return_value_policy< bp::copy_const_reference >() );
    
    }

}
