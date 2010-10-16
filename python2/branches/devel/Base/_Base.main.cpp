// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License


#include "boost/python.hpp"

#include "boost/python/suite/indexing/vector_indexing_suite.hpp"

#include "Helpers/clone_const_reference.hpp"

#include "Array2DBase.pypp.hpp"

#include "Array2D_double_.pypp.hpp"

#include "ChunkedVector_double_.pypp.hpp"

#include "CombineProperties.pypp.hpp"

#include "FlopsMark.pypp.hpp"

#include "Incremint.pypp.hpp"

#include "LinkToProperty.pypp.hpp"

#include "LowerCaseString.pypp.hpp"

#include "MemInfo.pypp.hpp"

#include "NoMangling.pypp.hpp"

#include "NullProperty.pypp.hpp"

#include "PackedArray2D_QString_.pypp.hpp"

#include "PackedArray2D_QString_Array.pypp.hpp"

#include "PackedArray2D_double_.pypp.hpp"

#include "PackedArray2D_double_Array.pypp.hpp"

#include "PackedArray2D_int_.pypp.hpp"

#include "PackedArray2D_int_Array.pypp.hpp"

#include "Process.pypp.hpp"

#include "Properties.pypp.hpp"

#include "Property.pypp.hpp"

#include "PropertyMap.pypp.hpp"

#include "PropertyName.pypp.hpp"

#include "StringMangler.pypp.hpp"

#include "TempDir.pypp.hpp"

#include "TrigArray2DBase.pypp.hpp"

#include "TrigArray2D_double_.pypp.hpp"

#include "TrimString.pypp.hpp"

#include "UpperCaseString.pypp.hpp"

#include "VariantProperty.pypp.hpp"

#include "Version.pypp.hpp"

#include "_Base_free_functions.pypp.hpp"

#include "vector_less__double__grate_.pypp.hpp"

namespace bp = boost::python;

#include "SireBase_containers.h"

#include "SireBase_registrars.h"

#include "SireBase_properties.h"

#include "SireBase/propertymap.h"

#include <QString>

BOOST_PYTHON_MODULE(_Base){
    register_SireBase_objects();

    register_SireBase_containers();

    register_vector_less__double__grate__class();

    register_Array2DBase_class();

    register_Array2D_double__class();

    register_ChunkedVector_double__class();

    register_Property_class();

    register_CombineProperties_class();

    register_FlopsMark_class();

    register_Incremint_class();

    register_LinkToProperty_class();

    register_StringMangler_class();

    register_LowerCaseString_class();

    register_MemInfo_class();

    register_NoMangling_class();

    register_NullProperty_class();

    register_PackedArray2D_QString__class();

    register_PackedArray2D_double__class();

    register_PackedArray2D_int__class();

    register_Process_class();

    register_Properties_class();

    register_PropertyMap_class();

    register_PropertyName_class();

    register_TempDir_class();

    register_TrigArray2DBase_class();

    register_TrigArray2D_double__class();

    register_TrimString_class();

    register_UpperCaseString_class();

    register_VariantProperty_class();

    register_Version_class();

    register_PackedArray2D_QString_Array_class();

    register_PackedArray2D_double_Array_class();

    register_PackedArray2D_int_Array_class();

    register_SireBase_properties();

    bp::implicitly_convertible< QString, SireBase::PropertyName >();

    bp::implicitly_convertible< QString, SireBase::VariantProperty >();

    bp::implicitly_convertible< double, SireBase::VariantProperty >();

    bp::implicitly_convertible< QString, SireBase::PropertyName >();

    bp::implicitly_convertible< SireBase::Property, SireBase::PropertyName >();

    bp::implicitly_convertible< QHash<QString,SireBase::PropertyName>, SireBase::PropertyMap >();

    register_free_functions();
}

