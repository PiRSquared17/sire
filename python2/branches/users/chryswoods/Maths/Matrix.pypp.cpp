// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Matrix.pypp.hpp"

namespace bp = boost::python;

#include "SireMaths/errors.h"

#include "Siren/errors.h"

#include "Siren/stream.h"

#include "maths.h"

#include "matrix.h"

#include "vector.h"

#include <QString>

#include <boost/scoped_array.hpp>

#include <gsl/gsl_eigen.h>

#include <gsl/gsl_matrix.h>

#include <gsl/gsl_vector.h>

#include "matrix.h"

#include "Siren/logger.h"

SireMaths::Matrix __copy__(const SireMaths::Matrix &other){ return SireMaths::Matrix(other); }

#include "Siren/str.hpp"

void register_Matrix_class(){

    { //::SireMaths::Matrix
        typedef bp::class_< SireMaths::Matrix, bp::bases< Siren::Primitive<SireMaths::Matrix> > > Matrix_exposer_t;
        Matrix_exposer_t Matrix_exposer = Matrix_exposer_t( "Matrix" );
        bp::scope Matrix_scope( Matrix_exposer );
        Matrix_exposer.def( bp::init< >() );
        Matrix_exposer.def( bp::init< double >(( bp::arg("diagonal_value") )) );
        Matrix_exposer.def( bp::init< double, double, double, double, double, double, double, double, double >(( bp::arg("xx"), bp::arg("xy"), bp::arg("xz"), bp::arg("yx"), bp::arg("yy"), bp::arg("yz"), bp::arg("zx"), bp::arg("zy"), bp::arg("zz") )) );
        Matrix_exposer.def( bp::init< SireMaths::Vector const &, SireMaths::Vector const &, SireMaths::Vector const & >(( bp::arg("r1"), bp::arg("r2"), bp::arg("r3") )) );
        Matrix_exposer.def( bp::init< boost::tuples::tuple< SireMaths::Vector, SireMaths::Vector, SireMaths::Vector, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type > const & >(( bp::arg("rows") )) );
        Matrix_exposer.def( bp::init< SireMaths::Matrix const & >(( bp::arg("m") )) );
        { //::SireMaths::Matrix::checkedOffset
        
            typedef int ( ::SireMaths::Matrix::*checkedOffset_function_type )( int,int ) const;
            checkedOffset_function_type checkedOffset_function_value( &::SireMaths::Matrix::checkedOffset );
            
            Matrix_exposer.def( 
                "checkedOffset"
                , checkedOffset_function_value
                , ( bp::arg("i"), bp::arg("j") ) );
        
        }
        { //::SireMaths::Matrix::column0
        
            typedef ::SireMaths::Vector ( ::SireMaths::Matrix::*column0_function_type )(  ) const;
            column0_function_type column0_function_value( &::SireMaths::Matrix::column0 );
            
            Matrix_exposer.def( 
                "column0"
                , column0_function_value );
        
        }
        { //::SireMaths::Matrix::column1
        
            typedef ::SireMaths::Vector ( ::SireMaths::Matrix::*column1_function_type )(  ) const;
            column1_function_type column1_function_value( &::SireMaths::Matrix::column1 );
            
            Matrix_exposer.def( 
                "column1"
                , column1_function_value );
        
        }
        { //::SireMaths::Matrix::column2
        
            typedef ::SireMaths::Vector ( ::SireMaths::Matrix::*column2_function_type )(  ) const;
            column2_function_type column2_function_value( &::SireMaths::Matrix::column2 );
            
            Matrix_exposer.def( 
                "column2"
                , column2_function_value );
        
        }
        { //::SireMaths::Matrix::determinant
        
            typedef double ( ::SireMaths::Matrix::*determinant_function_type )(  ) const;
            determinant_function_type determinant_function_value( &::SireMaths::Matrix::determinant );
            
            Matrix_exposer.def( 
                "determinant"
                , determinant_function_value );
        
        }
        { //::SireMaths::Matrix::enforceSymmetric
        
            typedef void ( ::SireMaths::Matrix::*enforceSymmetric_function_type )(  ) ;
            enforceSymmetric_function_type enforceSymmetric_function_value( &::SireMaths::Matrix::enforceSymmetric );
            
            Matrix_exposer.def( 
                "enforceSymmetric"
                , enforceSymmetric_function_value );
        
        }
        { //::SireMaths::Matrix::getPrincipalAxes
        
            typedef ::SireMaths::Matrix ( ::SireMaths::Matrix::*getPrincipalAxes_function_type )(  ) const;
            getPrincipalAxes_function_type getPrincipalAxes_function_value( &::SireMaths::Matrix::getPrincipalAxes );
            
            Matrix_exposer.def( 
                "getPrincipalAxes"
                , getPrincipalAxes_function_value );
        
        }
        { //::SireMaths::Matrix::hashCode
        
            typedef ::uint ( ::SireMaths::Matrix::*hashCode_function_type )(  ) const;
            hashCode_function_type hashCode_function_value( &::SireMaths::Matrix::hashCode );
            
            Matrix_exposer.def( 
                "hashCode"
                , hashCode_function_value );
        
        }
        { //::SireMaths::Matrix::identity
        
            typedef ::SireMaths::Matrix ( *identity_function_type )(  );
            identity_function_type identity_function_value( &::SireMaths::Matrix::identity );
            
            Matrix_exposer.def( 
                "identity"
                , identity_function_value );
        
        }
        { //::SireMaths::Matrix::inverse
        
            typedef ::SireMaths::Matrix ( ::SireMaths::Matrix::*inverse_function_type )(  ) const;
            inverse_function_type inverse_function_value( &::SireMaths::Matrix::inverse );
            
            Matrix_exposer.def( 
                "inverse"
                , inverse_function_value );
        
        }
        { //::SireMaths::Matrix::isIdentity
        
            typedef bool ( ::SireMaths::Matrix::*isIdentity_function_type )(  ) const;
            isIdentity_function_type isIdentity_function_value( &::SireMaths::Matrix::isIdentity );
            
            Matrix_exposer.def( 
                "isIdentity"
                , isIdentity_function_value );
        
        }
        { //::SireMaths::Matrix::isSymmetric
        
            typedef bool ( ::SireMaths::Matrix::*isSymmetric_function_type )(  ) const;
            isSymmetric_function_type isSymmetric_function_value( &::SireMaths::Matrix::isSymmetric );
            
            Matrix_exposer.def( 
                "isSymmetric"
                , isSymmetric_function_value );
        
        }
        { //::SireMaths::Matrix::offset
        
            typedef int ( ::SireMaths::Matrix::*offset_function_type )( int,int ) const;
            offset_function_type offset_function_value( &::SireMaths::Matrix::offset );
            
            Matrix_exposer.def( 
                "offset"
                , offset_function_value
                , ( bp::arg("i"), bp::arg("j") ) );
        
        }
        { //::SireMaths::Matrix::operator()
        
            typedef double const & ( ::SireMaths::Matrix::*__call___function_type )( int,int ) const;
            __call___function_type __call___function_value( &::SireMaths::Matrix::operator() );
            
            Matrix_exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("i"), bp::arg("j") )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        Matrix_exposer.def( bp::self == bp::self );
        { //::SireMaths::Matrix::row0
        
            typedef ::SireMaths::Vector ( ::SireMaths::Matrix::*row0_function_type )(  ) const;
            row0_function_type row0_function_value( &::SireMaths::Matrix::row0 );
            
            Matrix_exposer.def( 
                "row0"
                , row0_function_value );
        
        }
        { //::SireMaths::Matrix::row1
        
            typedef ::SireMaths::Vector ( ::SireMaths::Matrix::*row1_function_type )(  ) const;
            row1_function_type row1_function_value( &::SireMaths::Matrix::row1 );
            
            Matrix_exposer.def( 
                "row1"
                , row1_function_value );
        
        }
        { //::SireMaths::Matrix::row2
        
            typedef ::SireMaths::Vector ( ::SireMaths::Matrix::*row2_function_type )(  ) const;
            row2_function_type row2_function_value( &::SireMaths::Matrix::row2 );
            
            Matrix_exposer.def( 
                "row2"
                , row2_function_value );
        
        }
        { //::SireMaths::Matrix::setToIdentity
        
            typedef void ( ::SireMaths::Matrix::*setToIdentity_function_type )(  ) ;
            setToIdentity_function_type setToIdentity_function_value( &::SireMaths::Matrix::setToIdentity );
            
            Matrix_exposer.def( 
                "setToIdentity"
                , setToIdentity_function_value );
        
        }
        { //::SireMaths::Matrix::stream
        
            typedef void ( ::SireMaths::Matrix::*stream_function_type )( ::Siren::Stream & ) ;
            stream_function_type stream_function_value( &::SireMaths::Matrix::stream );
            
            Matrix_exposer.def( 
                "stream"
                , stream_function_value
                , ( bp::arg("s") ) );
        
        }
        { //::SireMaths::Matrix::toString
        
            typedef ::QString ( ::SireMaths::Matrix::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireMaths::Matrix::toString );
            
            Matrix_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireMaths::Matrix::trace
        
            typedef ::SireMaths::Vector ( ::SireMaths::Matrix::*trace_function_type )(  ) const;
            trace_function_type trace_function_value( &::SireMaths::Matrix::trace );
            
            Matrix_exposer.def( 
                "trace"
                , trace_function_value );
        
        }
        { //::SireMaths::Matrix::transpose
        
            typedef ::SireMaths::Matrix ( ::SireMaths::Matrix::*transpose_function_type )(  ) const;
            transpose_function_type transpose_function_value( &::SireMaths::Matrix::transpose );
            
            Matrix_exposer.def( 
                "transpose"
                , transpose_function_value );
        
        }
        { //::SireMaths::Matrix::xx
        
            typedef double ( ::SireMaths::Matrix::*xx_function_type )(  ) const;
            xx_function_type xx_function_value( &::SireMaths::Matrix::xx );
            
            Matrix_exposer.def( 
                "xx"
                , xx_function_value );
        
        }
        { //::SireMaths::Matrix::xy
        
            typedef double ( ::SireMaths::Matrix::*xy_function_type )(  ) const;
            xy_function_type xy_function_value( &::SireMaths::Matrix::xy );
            
            Matrix_exposer.def( 
                "xy"
                , xy_function_value );
        
        }
        { //::SireMaths::Matrix::xz
        
            typedef double ( ::SireMaths::Matrix::*xz_function_type )(  ) const;
            xz_function_type xz_function_value( &::SireMaths::Matrix::xz );
            
            Matrix_exposer.def( 
                "xz"
                , xz_function_value );
        
        }
        { //::SireMaths::Matrix::yx
        
            typedef double ( ::SireMaths::Matrix::*yx_function_type )(  ) const;
            yx_function_type yx_function_value( &::SireMaths::Matrix::yx );
            
            Matrix_exposer.def( 
                "yx"
                , yx_function_value );
        
        }
        { //::SireMaths::Matrix::yy
        
            typedef double ( ::SireMaths::Matrix::*yy_function_type )(  ) const;
            yy_function_type yy_function_value( &::SireMaths::Matrix::yy );
            
            Matrix_exposer.def( 
                "yy"
                , yy_function_value );
        
        }
        { //::SireMaths::Matrix::yz
        
            typedef double ( ::SireMaths::Matrix::*yz_function_type )(  ) const;
            yz_function_type yz_function_value( &::SireMaths::Matrix::yz );
            
            Matrix_exposer.def( 
                "yz"
                , yz_function_value );
        
        }
        { //::SireMaths::Matrix::zero
        
            typedef ::SireMaths::Matrix ( *zero_function_type )(  );
            zero_function_type zero_function_value( &::SireMaths::Matrix::zero );
            
            Matrix_exposer.def( 
                "zero"
                , zero_function_value );
        
        }
        { //::SireMaths::Matrix::zx
        
            typedef double ( ::SireMaths::Matrix::*zx_function_type )(  ) const;
            zx_function_type zx_function_value( &::SireMaths::Matrix::zx );
            
            Matrix_exposer.def( 
                "zx"
                , zx_function_value );
        
        }
        { //::SireMaths::Matrix::zy
        
            typedef double ( ::SireMaths::Matrix::*zy_function_type )(  ) const;
            zy_function_type zy_function_value( &::SireMaths::Matrix::zy );
            
            Matrix_exposer.def( 
                "zy"
                , zy_function_value );
        
        }
        { //::SireMaths::Matrix::zz
        
            typedef double ( ::SireMaths::Matrix::*zz_function_type )(  ) const;
            zz_function_type zz_function_value( &::SireMaths::Matrix::zz );
            
            Matrix_exposer.def( 
                "zz"
                , zz_function_value );
        
        }
        Matrix_exposer.staticmethod( "identity" );
        Matrix_exposer.staticmethod( "zero" );
        Matrix_exposer.def( bp::other< double >() * bp::self );
        Matrix_exposer.def( bp::self * bp::other< double >() );
        Matrix_exposer.def( bp::self * bp::self );
        Matrix_exposer.def( bp::self * bp::other< SireMaths::Vector >() );
        Matrix_exposer.def( bp::self + bp::self );
        Matrix_exposer.def( bp::self - bp::self );
        Matrix_exposer.def( "__copy__", &__copy__);
        Matrix_exposer.def( "__deepcopy__", &__copy__);
        Matrix_exposer.def( "clone", &__copy__);
        Matrix_exposer.def( "__str__", &__str__< ::SireMaths::Matrix > );
        Matrix_exposer.def( "__repr__", &__str__< ::SireMaths::Matrix > );
    }

}
