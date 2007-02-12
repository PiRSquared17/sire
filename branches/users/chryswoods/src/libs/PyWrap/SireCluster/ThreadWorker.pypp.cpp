// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License

#include "ThreadWorker.pypp.hpp"
#include "boost/python.hpp"
#include "sirecluster_headers.h"

namespace bp = boost::python;

void register_ThreadWorker_class(){

    bp::class_< SireCluster::ThreadWorker, bp::bases< SireCluster::WorkerBase >, boost::noncopyable >( "ThreadWorker", bp::no_init );

}
