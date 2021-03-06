// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "TrigMatrix.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/array2d.hpp"

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "matrix.h"

#include "nmatrix.h"

#include "nvector.h"

#include "sire_blas.h"

#include "sire_lapack.h"

#include "trigmatrix.h"

#include "vector.h"

#include "trigmatrix.h"

SireMaths::TrigMatrix __copy__(const SireMaths::TrigMatrix &other){ return SireMaths::TrigMatrix(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/len.hpp"

void register_TrigMatrix_class(){

    { //::SireMaths::TrigMatrix
        typedef bp::class_< SireMaths::TrigMatrix > TrigMatrix_exposer_t;
        TrigMatrix_exposer_t TrigMatrix_exposer = TrigMatrix_exposer_t( "TrigMatrix", bp::init< >() );
        bp::scope TrigMatrix_scope( TrigMatrix_exposer );
        TrigMatrix_exposer.def( bp::init< int >(( bp::arg("dimension") )) );
        TrigMatrix_exposer.def( bp::init< int, double >(( bp::arg("dimension"), bp::arg("initial_value") )) );
        TrigMatrix_exposer.def( bp::init< SireMaths::NMatrix const &, bp::optional< bool > >(( bp::arg("matrix"), bp::arg("take_upper")=(bool)(true) )) );
        TrigMatrix_exposer.def( bp::init< SireMaths::TrigMatrix const & >(( bp::arg("other") )) );
        { //::SireMaths::TrigMatrix::assertNColumns
        
            typedef void ( ::SireMaths::TrigMatrix::*assertNColumns_function_type )( int ) const;
            assertNColumns_function_type assertNColumns_function_value( &::SireMaths::TrigMatrix::assertNColumns );
            
            TrigMatrix_exposer.def( 
                "assertNColumns"
                , assertNColumns_function_value
                , ( bp::arg("ncolumns") ) );
        
        }
        { //::SireMaths::TrigMatrix::assertNRows
        
            typedef void ( ::SireMaths::TrigMatrix::*assertNRows_function_type )( int ) const;
            assertNRows_function_type assertNRows_function_value( &::SireMaths::TrigMatrix::assertNRows );
            
            TrigMatrix_exposer.def( 
                "assertNRows"
                , assertNRows_function_value
                , ( bp::arg("nrows") ) );
        
        }
        { //::SireMaths::TrigMatrix::assertSquare
        
            typedef void ( ::SireMaths::TrigMatrix::*assertSquare_function_type )(  ) const;
            assertSquare_function_type assertSquare_function_value( &::SireMaths::TrigMatrix::assertSquare );
            
            TrigMatrix_exposer.def( 
                "assertSquare"
                , assertSquare_function_value );
        
        }
        { //::SireMaths::TrigMatrix::assertValidColumn
        
            typedef void ( ::SireMaths::TrigMatrix::*assertValidColumn_function_type )( int ) const;
            assertValidColumn_function_type assertValidColumn_function_value( &::SireMaths::TrigMatrix::assertValidColumn );
            
            TrigMatrix_exposer.def( 
                "assertValidColumn"
                , assertValidColumn_function_value
                , ( bp::arg("j") ) );
        
        }
        { //::SireMaths::TrigMatrix::assertValidIndex
        
            typedef void ( ::SireMaths::TrigMatrix::*assertValidIndex_function_type )( int,int ) const;
            assertValidIndex_function_type assertValidIndex_function_value( &::SireMaths::TrigMatrix::assertValidIndex );
            
            TrigMatrix_exposer.def( 
                "assertValidIndex"
                , assertValidIndex_function_value
                , ( bp::arg("i"), bp::arg("j") ) );
        
        }
        { //::SireMaths::TrigMatrix::assertValidRow
        
            typedef void ( ::SireMaths::TrigMatrix::*assertValidRow_function_type )( int ) const;
            assertValidRow_function_type assertValidRow_function_value( &::SireMaths::TrigMatrix::assertValidRow );
            
            TrigMatrix_exposer.def( 
                "assertValidRow"
                , assertValidRow_function_value
                , ( bp::arg("i") ) );
        
        }
        { //::SireMaths::TrigMatrix::checkedOffset
        
            typedef int ( ::SireMaths::TrigMatrix::*checkedOffset_function_type )( int,int ) const;
            checkedOffset_function_type checkedOffset_function_value( &::SireMaths::TrigMatrix::checkedOffset );
            
            TrigMatrix_exposer.def( 
                "checkedOffset"
                , checkedOffset_function_value
                , ( bp::arg("i"), bp::arg("j") ) );
        
        }
        { //::SireMaths::TrigMatrix::column
        
            typedef ::SireMaths::NVector ( ::SireMaths::TrigMatrix::*column_function_type )( int ) const;
            column_function_type column_function_value( &::SireMaths::TrigMatrix::column );
            
            TrigMatrix_exposer.def( 
                "column"
                , column_function_value
                , ( bp::arg("j") ) );
        
        }
        { //::SireMaths::TrigMatrix::count
        
            typedef int ( ::SireMaths::TrigMatrix::*count_function_type )(  ) const;
            count_function_type count_function_value( &::SireMaths::TrigMatrix::count );
            
            TrigMatrix_exposer.def( 
                "count"
                , count_function_value );
        
        }
        { //::SireMaths::TrigMatrix::determinant
        
            typedef double ( ::SireMaths::TrigMatrix::*determinant_function_type )(  ) const;
            determinant_function_type determinant_function_value( &::SireMaths::TrigMatrix::determinant );
            
            TrigMatrix_exposer.def( 
                "determinant"
                , determinant_function_value );
        
        }
        { //::SireMaths::TrigMatrix::diagonal
        
            typedef ::SireMaths::NVector ( ::SireMaths::TrigMatrix::*diagonal_function_type )(  ) const;
            diagonal_function_type diagonal_function_value( &::SireMaths::TrigMatrix::diagonal );
            
            TrigMatrix_exposer.def( 
                "diagonal"
                , diagonal_function_value );
        
        }
        { //::SireMaths::TrigMatrix::diagonalise
        
            typedef ::std::pair< SireMaths::NVector, SireMaths::NMatrix > ( ::SireMaths::TrigMatrix::*diagonalise_function_type )(  ) const;
            diagonalise_function_type diagonalise_function_value( &::SireMaths::TrigMatrix::diagonalise );
            
            TrigMatrix_exposer.def( 
                "diagonalise"
                , diagonalise_function_value );
        
        }
        { //::SireMaths::TrigMatrix::fullTranspose
        
            typedef ::SireMaths::TrigMatrix ( ::SireMaths::TrigMatrix::*fullTranspose_function_type )(  ) const;
            fullTranspose_function_type fullTranspose_function_value( &::SireMaths::TrigMatrix::fullTranspose );
            
            TrigMatrix_exposer.def( 
                "fullTranspose"
                , fullTranspose_function_value );
        
        }
        { //::SireMaths::TrigMatrix::inverse
        
            typedef ::SireMaths::TrigMatrix ( ::SireMaths::TrigMatrix::*inverse_function_type )(  ) const;
            inverse_function_type inverse_function_value( &::SireMaths::TrigMatrix::inverse );
            
            TrigMatrix_exposer.def( 
                "inverse"
                , inverse_function_value );
        
        }
        { //::SireMaths::TrigMatrix::isTransposed
        
            typedef bool ( ::SireMaths::TrigMatrix::*isTransposed_function_type )(  ) const;
            isTransposed_function_type isTransposed_function_value( &::SireMaths::TrigMatrix::isTransposed );
            
            TrigMatrix_exposer.def( 
                "isTransposed"
                , isTransposed_function_value );
        
        }
        { //::SireMaths::TrigMatrix::memory
        
            typedef ::QVector< double > ( ::SireMaths::TrigMatrix::*memory_function_type )(  ) const;
            memory_function_type memory_function_value( &::SireMaths::TrigMatrix::memory );
            
            TrigMatrix_exposer.def( 
                "memory"
                , memory_function_value );
        
        }
        { //::SireMaths::TrigMatrix::nColumns
        
            typedef int ( ::SireMaths::TrigMatrix::*nColumns_function_type )(  ) const;
            nColumns_function_type nColumns_function_value( &::SireMaths::TrigMatrix::nColumns );
            
            TrigMatrix_exposer.def( 
                "nColumns"
                , nColumns_function_value );
        
        }
        { //::SireMaths::TrigMatrix::nRows
        
            typedef int ( ::SireMaths::TrigMatrix::*nRows_function_type )(  ) const;
            nRows_function_type nRows_function_value( &::SireMaths::TrigMatrix::nRows );
            
            TrigMatrix_exposer.def( 
                "nRows"
                , nRows_function_value );
        
        }
        { //::SireMaths::TrigMatrix::offset
        
            typedef int ( ::SireMaths::TrigMatrix::*offset_function_type )( int,int ) const;
            offset_function_type offset_function_value( &::SireMaths::TrigMatrix::offset );
            
            TrigMatrix_exposer.def( 
                "offset"
                , offset_function_value
                , ( bp::arg("i"), bp::arg("j") ) );
        
        }
        TrigMatrix_exposer.def( bp::self != bp::self );
        { //::SireMaths::TrigMatrix::operator()
        
            typedef double const & ( ::SireMaths::TrigMatrix::*__call___function_type )( int,int ) const;
            __call___function_type __call___function_value( &::SireMaths::TrigMatrix::operator() );
            
            TrigMatrix_exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("i"), bp::arg("j") )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        TrigMatrix_exposer.def( bp::self * bp::self );
        TrigMatrix_exposer.def( bp::self * bp::other< double >() );
        TrigMatrix_exposer.def( bp::self * bp::other< SireMaths::NVector >() );
        TrigMatrix_exposer.def( bp::self * bp::other< SireMaths::Vector >() );
        TrigMatrix_exposer.def( bp::self + bp::self );
        TrigMatrix_exposer.def( -bp::self );
        TrigMatrix_exposer.def( bp::self - bp::self );
        TrigMatrix_exposer.def( bp::self / bp::self );
        TrigMatrix_exposer.def( bp::self / bp::other< double >() );
        { //::SireMaths::TrigMatrix::operator=
        
            typedef ::SireMaths::TrigMatrix & ( ::SireMaths::TrigMatrix::*assign_function_type )( ::SireMaths::TrigMatrix const & ) ;
            assign_function_type assign_function_value( &::SireMaths::TrigMatrix::operator= );
            
            TrigMatrix_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        TrigMatrix_exposer.def( bp::self == bp::self );
        { //::SireMaths::TrigMatrix::redimension
        
            typedef void ( ::SireMaths::TrigMatrix::*redimension_function_type )( int ) ;
            redimension_function_type redimension_function_value( &::SireMaths::TrigMatrix::redimension );
            
            TrigMatrix_exposer.def( 
                "redimension"
                , redimension_function_value
                , ( bp::arg("dimension") ) );
        
        }
        { //::SireMaths::TrigMatrix::row
        
            typedef ::SireMaths::NVector ( ::SireMaths::TrigMatrix::*row_function_type )( int ) const;
            row_function_type row_function_value( &::SireMaths::TrigMatrix::row );
            
            TrigMatrix_exposer.def( 
                "row"
                , row_function_value
                , ( bp::arg("i") ) );
        
        }
        { //::SireMaths::TrigMatrix::set
        
            typedef void ( ::SireMaths::TrigMatrix::*set_function_type )( int,int,double ) ;
            set_function_type set_function_value( &::SireMaths::TrigMatrix::set );
            
            TrigMatrix_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("i"), bp::arg("j"), bp::arg("value") ) );
        
        }
        { //::SireMaths::TrigMatrix::setAll
        
            typedef void ( ::SireMaths::TrigMatrix::*setAll_function_type )( double ) ;
            setAll_function_type setAll_function_value( &::SireMaths::TrigMatrix::setAll );
            
            TrigMatrix_exposer.def( 
                "setAll"
                , setAll_function_value
                , ( bp::arg("value") ) );
        
        }
        { //::SireMaths::TrigMatrix::setColumn
        
            typedef void ( ::SireMaths::TrigMatrix::*setColumn_function_type )( int,double ) ;
            setColumn_function_type setColumn_function_value( &::SireMaths::TrigMatrix::setColumn );
            
            TrigMatrix_exposer.def( 
                "setColumn"
                , setColumn_function_value
                , ( bp::arg("j"), bp::arg("value") ) );
        
        }
        { //::SireMaths::TrigMatrix::setColumn
        
            typedef void ( ::SireMaths::TrigMatrix::*setColumn_function_type )( int,::SireMaths::NVector const & ) ;
            setColumn_function_type setColumn_function_value( &::SireMaths::TrigMatrix::setColumn );
            
            TrigMatrix_exposer.def( 
                "setColumn"
                , setColumn_function_value
                , ( bp::arg("j"), bp::arg("column") ) );
        
        }
        { //::SireMaths::TrigMatrix::setRow
        
            typedef void ( ::SireMaths::TrigMatrix::*setRow_function_type )( int,double ) ;
            setRow_function_type setRow_function_value( &::SireMaths::TrigMatrix::setRow );
            
            TrigMatrix_exposer.def( 
                "setRow"
                , setRow_function_value
                , ( bp::arg("i"), bp::arg("value") ) );
        
        }
        { //::SireMaths::TrigMatrix::setRow
        
            typedef void ( ::SireMaths::TrigMatrix::*setRow_function_type )( int,::SireMaths::NVector const & ) ;
            setRow_function_type setRow_function_value( &::SireMaths::TrigMatrix::setRow );
            
            TrigMatrix_exposer.def( 
                "setRow"
                , setRow_function_value
                , ( bp::arg("i"), bp::arg("row") ) );
        
        }
        { //::SireMaths::TrigMatrix::size
        
            typedef int ( ::SireMaths::TrigMatrix::*size_function_type )(  ) const;
            size_function_type size_function_value( &::SireMaths::TrigMatrix::size );
            
            TrigMatrix_exposer.def( 
                "size"
                , size_function_value );
        
        }
        { //::SireMaths::TrigMatrix::toString
        
            typedef ::QString ( ::SireMaths::TrigMatrix::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireMaths::TrigMatrix::toString );
            
            TrigMatrix_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireMaths::TrigMatrix::trace
        
            typedef double ( ::SireMaths::TrigMatrix::*trace_function_type )(  ) const;
            trace_function_type trace_function_value( &::SireMaths::TrigMatrix::trace );
            
            TrigMatrix_exposer.def( 
                "trace"
                , trace_function_value );
        
        }
        { //::SireMaths::TrigMatrix::transpose
        
            typedef ::SireMaths::TrigMatrix ( ::SireMaths::TrigMatrix::*transpose_function_type )(  ) const;
            transpose_function_type transpose_function_value( &::SireMaths::TrigMatrix::transpose );
            
            TrigMatrix_exposer.def( 
                "transpose"
                , transpose_function_value );
        
        }
        { //::SireMaths::TrigMatrix::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMaths::TrigMatrix::typeName );
            
            TrigMatrix_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireMaths::TrigMatrix::what
        
            typedef char const * ( ::SireMaths::TrigMatrix::*what_function_type )(  ) const;
            what_function_type what_function_value( &::SireMaths::TrigMatrix::what );
            
            TrigMatrix_exposer.def( 
                "what"
                , what_function_value );
        
        }
        TrigMatrix_exposer.staticmethod( "typeName" );
        TrigMatrix_exposer.def( "__copy__", &__copy__);
        TrigMatrix_exposer.def( "__deepcopy__", &__copy__);
        TrigMatrix_exposer.def( "clone", &__copy__);
        TrigMatrix_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMaths::TrigMatrix >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        TrigMatrix_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMaths::TrigMatrix >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        TrigMatrix_exposer.def( "__str__", &__str__< ::SireMaths::TrigMatrix > );
        TrigMatrix_exposer.def( "__repr__", &__str__< ::SireMaths::TrigMatrix > );
        TrigMatrix_exposer.def( "__len__", &__len_size< ::SireMaths::TrigMatrix > );
    }

}
