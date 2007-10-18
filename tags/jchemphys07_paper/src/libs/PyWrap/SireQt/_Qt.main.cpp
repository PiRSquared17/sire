// This file has been generated by Py++.

// (C) Christopher Woods, GPL >=2 License


#include "boost/python.hpp"

#include "sireqt_headers.h"

#include "Filters.pypp.hpp"

#include "NumberFlags.pypp.hpp"

#include "OpenMode.pypp.hpp"

#include "Permissions.pypp.hpp"

#include "QByteArray.pypp.hpp"

#include "QChar.pypp.hpp"

#include "QDataStream.pypp.hpp"

#include "QDate.pypp.hpp"

#include "QDateTime.pypp.hpp"

#include "QDir.pypp.hpp"

#include "QFile.pypp.hpp"

#include "QFileInfo.pypp.hpp"

#include "QIODevice.pypp.hpp"

#include "QTextStream.pypp.hpp"

#include "QTime.pypp.hpp"

#include "QVariant.pypp.hpp"

#include "SortFlags.pypp.hpp"

#include "QString.py.h"

#include "QtGlobal.py.h"

#include "sireqt_containers.h"

namespace bp = boost::python;

BOOST_PYTHON_MODULE(_Qt){
    export_QString();

    export_QtGlobal();

    register_SireQt_containers();

    register_QByteArray_class();

    register_QChar_class();

    register_QDataStream_class();

    register_QDate_class();

    register_QDateTime_class();

    register_Filters_class();

    bp::implicitly_convertible< QFlags<QDir::Filter>, int >();

    register_SortFlags_class();

    bp::implicitly_convertible< QFlags<QDir::SortFlag>, int >();

    register_QDir_class();

    register_QIODevice_class();

    register_QFile_class();

    register_QFileInfo_class();

    register_Permissions_class();

    bp::implicitly_convertible< QFlags<QFile::Permission>, int >();

    register_OpenMode_class();

    bp::implicitly_convertible< QFlags<QIODevice::OpenModeFlag>, int >();

    register_QTextStream_class();

    register_NumberFlags_class();

    bp::implicitly_convertible< QFlags<QTextStream::NumberFlag>, int >();

    register_QTime_class();

    bp::implicitly_convertible< QDir::SortFlag, QFlags<QDir::SortFlag> >();

    bp::implicitly_convertible< QDir::Filter, QFlags<QDir::Filter> >();

    bp::implicitly_convertible< QFile::Permission, QFlags<QFile::Permission> >();

    bp::implicitly_convertible< QIODevice::OpenModeFlag, QFlags<QIODevice::OpenModeFlag> >();

    bp::implicitly_convertible< QTextStream::NumberFlag, QFlags<QTextStream::NumberFlag> >();

    bp::implicitly_convertible< QString, QFileInfo >();

    bp::implicitly_convertible< QString, QFile >();

    bp::implicitly_convertible< QString, QDir >();

    register_QVariant_class();
}

