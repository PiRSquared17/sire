// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "assign_internals_MolBondInfo_.pypp.hpp"
#include "boost/python.hpp"
#include "siremm_headers.h"
#include "SireMol/molecule.h"
#include "SireMol/residue.h"
#include "SireMol/newatom.h"
#include "SireMol/atom.h"

namespace bp = boost::python;

void register_assign_internals_MolBondInfo__class(){

    bp::class_< SireMM::assign_internals<SireMM::MolBondInfo>, bp::bases< SireDB::AssignBase >, boost::noncopyable >( "assign_internals_MolBondInfo_", bp::no_init );

}
