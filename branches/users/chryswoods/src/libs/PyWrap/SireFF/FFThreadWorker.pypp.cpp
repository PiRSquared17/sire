// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "FFThreadWorker.pypp.hpp"
#include "boost/python.hpp"
#include "sireff_headers.h"
#include "SireMol/molecule.h"
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

const char* pvt_get_name(const SireFF::FFThreadWorker&){ return "SireFF::FFThreadWorker";}

void register_FFThreadWorker_class(){

    bp::class_< SireFF::FFThreadWorker, bp::bases< SireFF::FFLocalWorker, SireCluster::ThreadWorker >, boost::noncopyable >( "FFThreadWorker", bp::no_init )    
        .def( "__str__", &pvt_get_name);

}
