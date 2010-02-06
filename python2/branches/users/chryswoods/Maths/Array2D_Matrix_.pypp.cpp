// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Array2D_Matrix_.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/array2d.hpp"

#include "SireBase/trigarray2d.hpp"

#include "Siren/errors.h"

#include "Siren/stream.h"

#include "Siren/streamqt.h"

#include "matrix.h"

#include "nmatrix.h"

#include "nvector.h"

#include "sire_blas.h"

#include "sire_lapack.h"

#include "sire_linpack.h"

#include "trigmatrix.h"

#include "vector.h"

#include "nmatrix.h"

#include "Siren/logger.h"

SireBase::Array2D<SireMaths::Matrix> __copy__(const SireBase::Array2D<SireMaths::Matrix> &other){ return SireBase::Array2D<SireMaths::Matrix>(other); }

#include "Siren/str.hpp"

void register_Array2D_Matrix__class(){

    { //::SireBase::Array2D< SireMaths::Matrix >
        typedef bp::class_< SireBase::Array2D< SireMaths::Matrix >, bp::bases< Siren::Mutable, SireBase::Array2DBase, Siren::Object > > Array2D_Matrix__exposer_t;
        Array2D_Matrix__exposer_t Array2D_Matrix__exposer = Array2D_Matrix__exposer_t( "Array2D_Matrix_" );
        bp::scope Array2D_Matrix__scope( Array2D_Matrix__exposer );
        Array2D_Matrix__exposer.def( bp::init< >() );
        Array2D_Matrix__exposer.def( bp::init< int, int >(( bp::arg("nrows"), bp::arg("ncolumns") )) );
        Array2D_Matrix__exposer.def( bp::init< int, int, SireMaths::Matrix const & >(( bp::arg("nrows"), bp::arg("ncolumns"), bp::arg("default_value") )) );
        Array2D_Matrix__exposer.def( bp::init< SireBase::Array2D< SireMaths::Matrix > const & >(( bp::arg("other") )) );
        { //::SireBase::Array2D< SireMaths::Matrix >::at
        
            typedef SireBase::Array2D< SireMaths::Matrix > exported_class_t;
            typedef ::SireMaths::Matrix const & ( ::SireBase::Array2D< SireMaths::Matrix >::*at_function_type )( int,int ) const;
            at_function_type at_function_value( &::SireBase::Array2D< SireMaths::Matrix >::at );
            
            Array2D_Matrix__exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("i"), bp::arg("j") )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireBase::Array2D< SireMaths::Matrix >::get
        
            typedef SireBase::Array2D< SireMaths::Matrix > exported_class_t;
            typedef ::SireMaths::Matrix const & ( ::SireBase::Array2D< SireMaths::Matrix >::*get_function_type )( int,int ) const;
            get_function_type get_function_value( &::SireBase::Array2D< SireMaths::Matrix >::get );
            
            Array2D_Matrix__exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("i"), bp::arg("j") )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        Array2D_Matrix__exposer.def( bp::self != bp::self );
        { //::SireBase::Array2D< SireMaths::Matrix >::operator()
        
            typedef SireBase::Array2D< SireMaths::Matrix > exported_class_t;
            typedef ::SireMaths::Matrix const & ( ::SireBase::Array2D< SireMaths::Matrix >::*__call___function_type )( int,int ) const;
            __call___function_type __call___function_value( &::SireBase::Array2D< SireMaths::Matrix >::operator() );
            
            Array2D_Matrix__exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("i"), bp::arg("j") )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireBase::Array2D< SireMaths::Matrix >::operator=
        
            typedef SireBase::Array2D< SireMaths::Matrix > exported_class_t;
            typedef ::SireBase::Array2D< SireMaths::Matrix > & ( ::SireBase::Array2D< SireMaths::Matrix >::*assign_function_type )( ::SireBase::Array2D< SireMaths::Matrix > const & ) ;
            assign_function_type assign_function_value( &::SireBase::Array2D< SireMaths::Matrix >::operator= );
            
            Array2D_Matrix__exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        Array2D_Matrix__exposer.def( bp::self == bp::self );
        { //::SireBase::Array2D< SireMaths::Matrix >::redimension
        
            typedef SireBase::Array2D< SireMaths::Matrix > exported_class_t;
            typedef void ( ::SireBase::Array2D< SireMaths::Matrix >::*redimension_function_type )( int,int ) ;
            redimension_function_type redimension_function_value( &::SireBase::Array2D< SireMaths::Matrix >::redimension );
            
            Array2D_Matrix__exposer.def( 
                "redimension"
                , redimension_function_value
                , ( bp::arg("nrows"), bp::arg("ncolumns") ) );
        
        }
        { //::SireBase::Array2D< SireMaths::Matrix >::set
        
            typedef SireBase::Array2D< SireMaths::Matrix > exported_class_t;
            typedef void ( ::SireBase::Array2D< SireMaths::Matrix >::*set_function_type )( int,int,::SireMaths::Matrix const & ) ;
            set_function_type set_function_value( &::SireBase::Array2D< SireMaths::Matrix >::set );
            
            Array2D_Matrix__exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("i"), bp::arg("j"), bp::arg("value") ) );
        
        }
        { //::SireBase::Array2D< SireMaths::Matrix >::setAll
        
            typedef SireBase::Array2D< SireMaths::Matrix > exported_class_t;
            typedef void ( ::SireBase::Array2D< SireMaths::Matrix >::*setAll_function_type )( ::SireMaths::Matrix const & ) ;
            setAll_function_type setAll_function_value( &::SireBase::Array2D< SireMaths::Matrix >::setAll );
            
            Array2D_Matrix__exposer.def( 
                "setAll"
                , setAll_function_value
                , ( bp::arg("value") ) );
        
        }
        { //::SireBase::Array2D< SireMaths::Matrix >::stream
        
            typedef SireBase::Array2D< SireMaths::Matrix > exported_class_t;
            typedef void ( ::SireBase::Array2D< SireMaths::Matrix >::*stream_function_type )( ::Siren::Stream & ) ;
            stream_function_type stream_function_value( &::SireBase::Array2D< SireMaths::Matrix >::stream );
            
            Array2D_Matrix__exposer.def( 
                "stream"
                , stream_function_value
                , ( bp::arg("s") ) );
        
        }
        { //::SireBase::Array2D< SireMaths::Matrix >::toString
        
            typedef SireBase::Array2D< SireMaths::Matrix > exported_class_t;
            typedef ::QString ( ::SireBase::Array2D< SireMaths::Matrix >::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireBase::Array2D< SireMaths::Matrix >::toString );
            
            Array2D_Matrix__exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireBase::Array2D< SireMaths::Matrix >::transpose
        
            typedef SireBase::Array2D< SireMaths::Matrix > exported_class_t;
            typedef ::SireBase::Array2D< SireMaths::Matrix > ( ::SireBase::Array2D< SireMaths::Matrix >::*transpose_function_type )(  ) const;
            transpose_function_type transpose_function_value( &::SireBase::Array2D< SireMaths::Matrix >::transpose );
            
            Array2D_Matrix__exposer.def( 
                "transpose"
                , transpose_function_value );
        
        }
        { //::SireBase::Array2D< SireMaths::Matrix >::typeName
        
            typedef SireBase::Array2D< SireMaths::Matrix > exported_class_t;
            typedef ::QString ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireBase::Array2D< SireMaths::Matrix >::typeName );
            
            Array2D_Matrix__exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        Array2D_Matrix__exposer.staticmethod( "typeName" );
        Array2D_Matrix__exposer.def( "__copy__", &__copy__);
        Array2D_Matrix__exposer.def( "__deepcopy__", &__copy__);
        Array2D_Matrix__exposer.def( "clone", &__copy__);
        Array2D_Matrix__exposer.def( "__str__", &__str__< ::SireBase::Array2D<SireMaths::Matrix> > );
        Array2D_Matrix__exposer.def( "__repr__", &__str__< ::SireBase::Array2D<SireMaths::Matrix> > );
    }

}