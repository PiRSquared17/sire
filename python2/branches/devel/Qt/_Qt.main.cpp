// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License


#include "boost/python.hpp"

#include "__array_1.pypp.hpp"

#include <QByteArray>

#include <QDateTime>

#include <QUuid>

#include <QBitArray>

#include "QByteArray.pypp.hpp"

#include "QDate.pypp.hpp"

#include "QDateTime.pypp.hpp"

#include "QTime.pypp.hpp"

#include "QUuid.pypp.hpp"

namespace bp = boost::python;

#include "sireqt_containers.h"

void autoconvert_QChar();

void autoconvert_QString();

BOOST_PYTHON_MODULE(_Qt){

    register_SireQt_containers();

    register_QByteArray_class();
   
    register_QDate_class();
 
   /* register_QDateTime_class(); */

    register_QTime_class();

    autoconvert_QChar();
    
    autoconvert_QString();

    register_QUuid_class(); 
}

