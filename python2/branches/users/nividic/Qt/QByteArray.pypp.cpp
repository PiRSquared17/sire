// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include <QString>
#include <QByteArray>
#include <QFile>
#include <QFileInfo>
#include <QDir>
#include <QTextStream>
#include <QDateTime>
#include <QLocale>
#include <QUuid>
#include <qnamespace.h>
#include <QVariant>
#include <QUrl>
#include <QBitArray>
#include "QByteArray.pypp.hpp"

namespace bp = boost::python;

QByteArray __copy__(const QByteArray &other){ return QByteArray(other); }

#include "Qt/qdatastream.hpp"

void register_QByteArray_class(){

    { //::QByteArray
        typedef bp::class_< QByteArray > QByteArray_exposer_t;
        QByteArray_exposer_t QByteArray_exposer = QByteArray_exposer_t( "QByteArray" );
        bp::scope QByteArray_scope( QByteArray_exposer );
        QByteArray_exposer.def( bp::init< >() );
        QByteArray_exposer.def( bp::init< char const * >(( bp::arg("arg0") )) );
        QByteArray_exposer.def( bp::init< char const *, int >(( bp::arg("arg0"), bp::arg("size") )) );
        QByteArray_exposer.def( bp::init< int, char >(( bp::arg("size"), bp::arg("c") )) );
        QByteArray_exposer.def( bp::init< QByteArray const & >(( bp::arg("arg0") )) );
        { //::QByteArray::append
        
            typedef ::QByteArray & ( ::QByteArray::*append_function_type )( char ) ;
            append_function_type append_function_value( &::QByteArray::append );
            
            QByteArray_exposer.def( 
                "append"
                , append_function_value
                , ( bp::arg("c") )
                , bp::return_self< >() );
        
        }
        { //::QByteArray::append
        
            typedef ::QByteArray & ( ::QByteArray::*append_function_type )( char const * ) ;
            append_function_type append_function_value( &::QByteArray::append );
            
            QByteArray_exposer.def( 
                "append"
                , append_function_value
                , ( bp::arg("s") )
                , bp::return_self< >() );
        
        }
        { //::QByteArray::append
        
            typedef ::QByteArray & ( ::QByteArray::*append_function_type )( ::QByteArray const & ) ;
            append_function_type append_function_value( &::QByteArray::append );
            
            QByteArray_exposer.def( 
                "append"
                , append_function_value
                , ( bp::arg("a") )
                , bp::return_self< >() );
        
        }
        { //::QByteArray::capacity
        
            typedef int ( ::QByteArray::*capacity_function_type )(  ) const;
            capacity_function_type capacity_function_value( &::QByteArray::capacity );
            
            QByteArray_exposer.def( 
                "capacity"
                , capacity_function_value );
        
        }
        { //::QByteArray::chop
        
            typedef void ( ::QByteArray::*chop_function_type )( int ) ;
            chop_function_type chop_function_value( &::QByteArray::chop );
            
            QByteArray_exposer.def( 
                "chop"
                , chop_function_value
                , ( bp::arg("n") ) );
        
        }
        { //::QByteArray::clear
        
            typedef void ( ::QByteArray::*clear_function_type )(  ) ;
            clear_function_type clear_function_value( &::QByteArray::clear );
            
            QByteArray_exposer.def( 
                "clear"
                , clear_function_value );
        
        }
        { //::QByteArray::constData
        
            typedef char const * ( ::QByteArray::*constData_function_type )(  ) const;
            constData_function_type constData_function_value( &::QByteArray::constData );
            
            QByteArray_exposer.def( 
                "constData"
                , constData_function_value );
        
        }
        { //::QByteArray::contains
        
            typedef ::QBool ( ::QByteArray::*contains_function_type )( char ) const;
            contains_function_type contains_function_value( &::QByteArray::contains );
            
            QByteArray_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("c") ) );
        
        }
        { //::QByteArray::contains
        
            typedef ::QBool ( ::QByteArray::*contains_function_type )( char const * ) const;
            contains_function_type contains_function_value( &::QByteArray::contains );
            
            QByteArray_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("a") ) );
        
        }
        { //::QByteArray::contains
        
            typedef ::QBool ( ::QByteArray::*contains_function_type )( ::QByteArray const & ) const;
            contains_function_type contains_function_value( &::QByteArray::contains );
            
            QByteArray_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("a") ) );
        
        }
        { //::QByteArray::count
        
            typedef int ( ::QByteArray::*count_function_type )( char ) const;
            count_function_type count_function_value( &::QByteArray::count );
            
            QByteArray_exposer.def( 
                "count"
                , count_function_value
                , ( bp::arg("c") ) );
        
        }
        { //::QByteArray::count
        
            typedef int ( ::QByteArray::*count_function_type )( char const * ) const;
            count_function_type count_function_value( &::QByteArray::count );
            
            QByteArray_exposer.def( 
                "count"
                , count_function_value
                , ( bp::arg("a") ) );
        
        }
        { //::QByteArray::count
        
            typedef int ( ::QByteArray::*count_function_type )( ::QByteArray const & ) const;
            count_function_type count_function_value( &::QByteArray::count );
            
            QByteArray_exposer.def( 
                "count"
                , count_function_value
                , ( bp::arg("a") ) );
        
        }
        { //::QByteArray::count
        
            typedef int ( ::QByteArray::*count_function_type )(  ) const;
            count_function_type count_function_value( &::QByteArray::count );
            
            QByteArray_exposer.def( 
                "count"
                , count_function_value );
        
        }
        { //::QByteArray::data
        
            typedef char const * ( ::QByteArray::*data_function_type )(  ) const;
            data_function_type data_function_value( &::QByteArray::data );
            
            QByteArray_exposer.def( 
                "data"
                , data_function_value );
        
        }
        { //::QByteArray::detach
        
            typedef void ( ::QByteArray::*detach_function_type )(  ) ;
            detach_function_type detach_function_value( &::QByteArray::detach );
            
            QByteArray_exposer.def( 
                "detach"
                , detach_function_value );
        
        }
        { //::QByteArray::endsWith
        
            typedef bool ( ::QByteArray::*endsWith_function_type )( ::QByteArray const & ) const;
            endsWith_function_type endsWith_function_value( &::QByteArray::endsWith );
            
            QByteArray_exposer.def( 
                "endsWith"
                , endsWith_function_value
                , ( bp::arg("a") ) );
        
        }
        { //::QByteArray::endsWith
        
            typedef bool ( ::QByteArray::*endsWith_function_type )( char ) const;
            endsWith_function_type endsWith_function_value( &::QByteArray::endsWith );
            
            QByteArray_exposer.def( 
                "endsWith"
                , endsWith_function_value
                , ( bp::arg("c") ) );
        
        }
        { //::QByteArray::endsWith
        
            typedef bool ( ::QByteArray::*endsWith_function_type )( char const * ) const;
            endsWith_function_type endsWith_function_value( &::QByteArray::endsWith );
            
            QByteArray_exposer.def( 
                "endsWith"
                , endsWith_function_value
                , ( bp::arg("c") ) );
        
        }
        { //::QByteArray::fill
        
            typedef ::QByteArray & ( ::QByteArray::*fill_function_type )( char,int ) ;
            fill_function_type fill_function_value( &::QByteArray::fill );
            
            QByteArray_exposer.def( 
                "fill"
                , fill_function_value
                , ( bp::arg("c"), bp::arg("size")=(int)(-0x00000000000000001) )
                , bp::return_self< >() );
        
        }
        { //::QByteArray::toBase64
        
            typedef ::QByteArray ( ::QByteArray::*toBase64_function_type )( ) const;
            toBase64_function_type toBase64_function_value( &::QByteArray::toBase64 );
            
            QByteArray_exposer.def( 
                "toBase64"
                , toBase64_function_value
                , ( bp::arg("toBase64") ) );
        
        }
        { //::QByteArray::fromBase64

            typedef ::QByteArray ( *fromBase64_function_type )( ::QByteArray const & );
            fromBase64_function_type fromBase64_function_value( &::QByteArray::fromBase64 );

            QByteArray_exposer.def(
                "fromBase64"
                , fromBase64_function_value
                , ( bp::arg("base64") ) );

        }
        { //::QByteArray::fromRawData
        
            typedef ::QByteArray ( *fromRawData_function_type )( char const *,int );
            fromRawData_function_type fromRawData_function_value( &::QByteArray::fromRawData );
            
            QByteArray_exposer.def( 
                "fromRawData"
                , fromRawData_function_value
                , ( bp::arg("arg0"), bp::arg("size") ) );
        
        }
        { //::QByteArray::indexOf
        
            typedef int ( ::QByteArray::*indexOf_function_type )( char,int ) const;
            indexOf_function_type indexOf_function_value( &::QByteArray::indexOf );
            
            QByteArray_exposer.def( 
                "indexOf"
                , indexOf_function_value
                , ( bp::arg("c"), bp::arg("from")=(int)(0) ) );
        
        }
        { //::QByteArray::indexOf
        
            typedef int ( ::QByteArray::*indexOf_function_type )( char const *,int ) const;
            indexOf_function_type indexOf_function_value( &::QByteArray::indexOf );
            
            QByteArray_exposer.def( 
                "indexOf"
                , indexOf_function_value
                , ( bp::arg("c"), bp::arg("from")=(int)(0) ) );
        
        }
        { //::QByteArray::indexOf
        
            typedef int ( ::QByteArray::*indexOf_function_type )( ::QByteArray const &,int ) const;
            indexOf_function_type indexOf_function_value( &::QByteArray::indexOf );
            
            QByteArray_exposer.def( 
                "indexOf"
                , indexOf_function_value
                , ( bp::arg("a"), bp::arg("from")=(int)(0) ) );
        
        }
        { //::QByteArray::insert
        
            typedef ::QByteArray & ( ::QByteArray::*insert_function_type )( int,char ) ;
            insert_function_type insert_function_value( &::QByteArray::insert );
            
            QByteArray_exposer.def( 
                "insert"
                , insert_function_value
                , ( bp::arg("i"), bp::arg("c") )
                , bp::return_self< >() );
        
        }
        { //::QByteArray::insert
        
            typedef ::QByteArray & ( ::QByteArray::*insert_function_type )( int,char const * ) ;
            insert_function_type insert_function_value( &::QByteArray::insert );
            
            QByteArray_exposer.def( 
                "insert"
                , insert_function_value
                , ( bp::arg("i"), bp::arg("s") )
                , bp::return_self< >() );
        
        }
        { //::QByteArray::insert
        
            typedef ::QByteArray & ( ::QByteArray::*insert_function_type )( int,::QByteArray const & ) ;
            insert_function_type insert_function_value( &::QByteArray::insert );
            
            QByteArray_exposer.def( 
                "insert"
                , insert_function_value
                , ( bp::arg("i"), bp::arg("a") )
                , bp::return_self< >() );
        
        }
        { //::QByteArray::isDetached
        
            typedef bool ( ::QByteArray::*isDetached_function_type )(  ) const;
            isDetached_function_type isDetached_function_value( &::QByteArray::isDetached );
            
            QByteArray_exposer.def( 
                "isDetached"
                , isDetached_function_value );
        
        }
        { //::QByteArray::isEmpty
        
            typedef bool ( ::QByteArray::*isEmpty_function_type )(  ) const;
            isEmpty_function_type isEmpty_function_value( &::QByteArray::isEmpty );
            
            QByteArray_exposer.def( 
                "isEmpty"
                , isEmpty_function_value );
        
        }
        { //::QByteArray::isNull
        
            typedef bool ( ::QByteArray::*isNull_function_type )(  ) const;
            isNull_function_type isNull_function_value( &::QByteArray::isNull );
            
            QByteArray_exposer.def( 
                "isNull"
                , isNull_function_value );
        
        }
        { //::QByteArray::lastIndexOf
        
            typedef int ( ::QByteArray::*lastIndexOf_function_type )( char,int ) const;
            lastIndexOf_function_type lastIndexOf_function_value( &::QByteArray::lastIndexOf );
            
            QByteArray_exposer.def( 
                "lastIndexOf"
                , lastIndexOf_function_value
                , ( bp::arg("c"), bp::arg("from")=(int)(-0x00000000000000001) ) );
        
        }
        { //::QByteArray::lastIndexOf
        
            typedef int ( ::QByteArray::*lastIndexOf_function_type )( char const *,int ) const;
            lastIndexOf_function_type lastIndexOf_function_value( &::QByteArray::lastIndexOf );
            
            QByteArray_exposer.def( 
                "lastIndexOf"
                , lastIndexOf_function_value
                , ( bp::arg("c"), bp::arg("from")=(int)(-0x00000000000000001) ) );
        
        }
        { //::QByteArray::lastIndexOf
        
            typedef int ( ::QByteArray::*lastIndexOf_function_type )( ::QByteArray const &,int ) const;
            lastIndexOf_function_type lastIndexOf_function_value( &::QByteArray::lastIndexOf );
            
            QByteArray_exposer.def( 
                "lastIndexOf"
                , lastIndexOf_function_value
                , ( bp::arg("a"), bp::arg("from")=(int)(-0x00000000000000001) ) );
        
        }
        { //::QByteArray::left
        
            typedef ::QByteArray ( ::QByteArray::*left_function_type )( int ) const;
            left_function_type left_function_value( &::QByteArray::left );
            
            QByteArray_exposer.def( 
                "left"
                , left_function_value
                , ( bp::arg("len") ) );
        
        }
        { //::QByteArray::leftJustified
        
            typedef ::QByteArray ( ::QByteArray::*leftJustified_function_type )( int,char,bool ) const;
            leftJustified_function_type leftJustified_function_value( &::QByteArray::leftJustified );
            
            QByteArray_exposer.def( 
                "leftJustified"
                , leftJustified_function_value
                , ( bp::arg("width"), bp::arg("fill")=(char)(' '), bp::arg("truncate")=(bool)(false) ) );
        
        }
        { //::QByteArray::length
        
            typedef int ( ::QByteArray::*length_function_type )(  ) const;
            length_function_type length_function_value( &::QByteArray::length );
            
            QByteArray_exposer.def( 
                "length"
                , length_function_value );
        
        }
        { //::QByteArray::mid
        
            typedef ::QByteArray ( ::QByteArray::*mid_function_type )( int,int ) const;
            mid_function_type mid_function_value( &::QByteArray::mid );
            
            QByteArray_exposer.def( 
                "mid"
                , mid_function_value
                , ( bp::arg("index"), bp::arg("len")=(int)(-0x00000000000000001) ) );
        
        }
        { //::QByteArray::number
        
            typedef ::QByteArray ( *number_function_type )( int,int );
            number_function_type number_function_value( &::QByteArray::number );
            
            QByteArray_exposer.def( 
                "number"
                , number_function_value
                , ( bp::arg("arg0"), bp::arg("base")=(int)(10) ) );
        
        }
        { //::QByteArray::number
        
            typedef ::QByteArray ( *number_function_type )( ::uint,int );
            number_function_type number_function_value( &::QByteArray::number );
            
            QByteArray_exposer.def( 
                "number"
                , number_function_value
                , ( bp::arg("arg0"), bp::arg("base")=(int)(10) ) );
        
        }
        { //::QByteArray::number
        
            typedef ::QByteArray ( *number_function_type )( double,char,int );
            number_function_type number_function_value( &::QByteArray::number );
            
            QByteArray_exposer.def( 
                "number"
                , number_function_value
                , ( bp::arg("arg0"), bp::arg("f")=(char)('g'), bp::arg("prec")=(int)(6) ) );
        
        }
        { //::QByteArray::operator=
        
            typedef ::QByteArray & ( ::QByteArray::*assign_function_type )( ::QByteArray const & ) ;
            assign_function_type assign_function_value( &::QByteArray::operator= );
            
            QByteArray_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        { //::QByteArray::prepend
        
            typedef ::QByteArray & ( ::QByteArray::*prepend_function_type )( char ) ;
            prepend_function_type prepend_function_value( &::QByteArray::prepend );
            
            QByteArray_exposer.def( 
                "prepend"
                , prepend_function_value
                , ( bp::arg("c") )
                , bp::return_self< >() );
        
        }
        { //::QByteArray::prepend
        
            typedef ::QByteArray & ( ::QByteArray::*prepend_function_type )( char const * ) ;
            prepend_function_type prepend_function_value( &::QByteArray::prepend );
            
            QByteArray_exposer.def( 
                "prepend"
                , prepend_function_value
                , ( bp::arg("s") )
                , bp::return_self< >() );
        
        }
        { //::QByteArray::prepend
        
            typedef ::QByteArray & ( ::QByteArray::*prepend_function_type )( ::QByteArray const & ) ;
            prepend_function_type prepend_function_value( &::QByteArray::prepend );
            
            QByteArray_exposer.def( 
                "prepend"
                , prepend_function_value
                , ( bp::arg("a") )
                , bp::return_self< >() );
        
        }
        { //::QByteArray::push_back
        
            typedef void ( ::QByteArray::*push_back_function_type )( char ) ;
            push_back_function_type push_back_function_value( &::QByteArray::push_back );
            
            QByteArray_exposer.def( 
                "push_back"
                , push_back_function_value
                , ( bp::arg("c") ) );
        
        }
        { //::QByteArray::push_back
        
            typedef void ( ::QByteArray::*push_back_function_type )( char const * ) ;
            push_back_function_type push_back_function_value( &::QByteArray::push_back );
            
            QByteArray_exposer.def( 
                "push_back"
                , push_back_function_value
                , ( bp::arg("c") ) );
        
        }
        { //::QByteArray::push_back
        
            typedef void ( ::QByteArray::*push_back_function_type )( ::QByteArray const & ) ;
            push_back_function_type push_back_function_value( &::QByteArray::push_back );
            
            QByteArray_exposer.def( 
                "push_back"
                , push_back_function_value
                , ( bp::arg("a") ) );
        
        }
        { //::QByteArray::push_front
        
            typedef void ( ::QByteArray::*push_front_function_type )( char ) ;
            push_front_function_type push_front_function_value( &::QByteArray::push_front );
            
            QByteArray_exposer.def( 
                "push_front"
                , push_front_function_value
                , ( bp::arg("c") ) );
        
        }
        { //::QByteArray::push_front
        
            typedef void ( ::QByteArray::*push_front_function_type )( char const * ) ;
            push_front_function_type push_front_function_value( &::QByteArray::push_front );
            
            QByteArray_exposer.def( 
                "push_front"
                , push_front_function_value
                , ( bp::arg("c") ) );
        
        }
        { //::QByteArray::push_front
        
            typedef void ( ::QByteArray::*push_front_function_type )( ::QByteArray const & ) ;
            push_front_function_type push_front_function_value( &::QByteArray::push_front );
            
            QByteArray_exposer.def( 
                "push_front"
                , push_front_function_value
                , ( bp::arg("a") ) );
        
        }
        { //::QByteArray::remove
        
            typedef ::QByteArray & ( ::QByteArray::*remove_function_type )( int,int ) ;
            remove_function_type remove_function_value( &::QByteArray::remove );
            
            QByteArray_exposer.def( 
                "remove"
                , remove_function_value
                , ( bp::arg("index"), bp::arg("len") )
                , bp::return_self< >() );
        
        }
        { //::QByteArray::replace
        
            typedef ::QByteArray & ( ::QByteArray::*replace_function_type )( int,int,char const * ) ;
            replace_function_type replace_function_value( &::QByteArray::replace );
            
            QByteArray_exposer.def( 
                "replace"
                , replace_function_value
                , ( bp::arg("index"), bp::arg("len"), bp::arg("s") )
                , bp::return_self< >() );
        
        }
        { //::QByteArray::replace
        
            typedef ::QByteArray & ( ::QByteArray::*replace_function_type )( int,int,::QByteArray const & ) ;
            replace_function_type replace_function_value( &::QByteArray::replace );
            
            QByteArray_exposer.def( 
                "replace"
                , replace_function_value
                , ( bp::arg("index"), bp::arg("len"), bp::arg("s") )
                , bp::return_self< >() );
        
        }
        { //::QByteArray::replace
        
            typedef ::QByteArray & ( ::QByteArray::*replace_function_type )( char,char const * ) ;
            replace_function_type replace_function_value( &::QByteArray::replace );
            
            QByteArray_exposer.def( 
                "replace"
                , replace_function_value
                , ( bp::arg("before"), bp::arg("after") )
                , bp::return_self< >() );
        
        }
        { //::QByteArray::replace
        
            typedef ::QByteArray & ( ::QByteArray::*replace_function_type )( char,::QByteArray const & ) ;
            replace_function_type replace_function_value( &::QByteArray::replace );
            
            QByteArray_exposer.def( 
                "replace"
                , replace_function_value
                , ( bp::arg("before"), bp::arg("after") )
                , bp::return_self< >() );
        
        }
        { //::QByteArray::replace
        
            typedef ::QByteArray & ( ::QByteArray::*replace_function_type )( char const *,char const * ) ;
            replace_function_type replace_function_value( &::QByteArray::replace );
            
            QByteArray_exposer.def( 
                "replace"
                , replace_function_value
                , ( bp::arg("before"), bp::arg("after") )
                , bp::return_self< >() );
        
        }
        { //::QByteArray::replace
        
            typedef ::QByteArray & ( ::QByteArray::*replace_function_type )( ::QByteArray const &,::QByteArray const & ) ;
            replace_function_type replace_function_value( &::QByteArray::replace );
            
            QByteArray_exposer.def( 
                "replace"
                , replace_function_value
                , ( bp::arg("before"), bp::arg("after") )
                , bp::return_self< >() );
        
        }
        { //::QByteArray::replace
        
            typedef ::QByteArray & ( ::QByteArray::*replace_function_type )( ::QByteArray const &,char const * ) ;
            replace_function_type replace_function_value( &::QByteArray::replace );
            
            QByteArray_exposer.def( 
                "replace"
                , replace_function_value
                , ( bp::arg("before"), bp::arg("after") )
                , bp::return_self< >() );
        
        }
        { //::QByteArray::replace
        
            typedef ::QByteArray & ( ::QByteArray::*replace_function_type )( char const *,::QByteArray const & ) ;
            replace_function_type replace_function_value( &::QByteArray::replace );
            
            QByteArray_exposer.def( 
                "replace"
                , replace_function_value
                , ( bp::arg("before"), bp::arg("after") )
                , bp::return_self< >() );
        
        }
        { //::QByteArray::replace
        
            typedef ::QByteArray & ( ::QByteArray::*replace_function_type )( char,char ) ;
            replace_function_type replace_function_value( &::QByteArray::replace );
            
            QByteArray_exposer.def( 
                "replace"
                , replace_function_value
                , ( bp::arg("before"), bp::arg("after") )
                , bp::return_self< >() );
        
        }
        { //::QByteArray::reserve
        
            typedef void ( ::QByteArray::*reserve_function_type )( int ) ;
            reserve_function_type reserve_function_value( &::QByteArray::reserve );
            
            QByteArray_exposer.def( 
                "reserve"
                , reserve_function_value
                , ( bp::arg("size") ) );
        
        }
        { //::QByteArray::resize
        
            typedef void ( ::QByteArray::*resize_function_type )( int ) ;
            resize_function_type resize_function_value( &::QByteArray::resize );
            
            QByteArray_exposer.def( 
                "resize"
                , resize_function_value
                , ( bp::arg("size") ) );
        
        }
        { //::QByteArray::right
        
            typedef ::QByteArray ( ::QByteArray::*right_function_type )( int ) const;
            right_function_type right_function_value( &::QByteArray::right );
            
            QByteArray_exposer.def( 
                "right"
                , right_function_value
                , ( bp::arg("len") ) );
        
        }
        { //::QByteArray::rightJustified
        
            typedef ::QByteArray ( ::QByteArray::*rightJustified_function_type )( int,char,bool ) const;
            rightJustified_function_type rightJustified_function_value( &::QByteArray::rightJustified );
            
            QByteArray_exposer.def( 
                "rightJustified"
                , rightJustified_function_value
                , ( bp::arg("width"), bp::arg("fill")=(char)(' '), bp::arg("truncate")=(bool)(false) ) );
        
        }
        { //::QByteArray::setNum
        
            typedef ::QByteArray & ( ::QByteArray::*setNum_function_type )( short int,int ) ;
            setNum_function_type setNum_function_value( &::QByteArray::setNum );
            
            QByteArray_exposer.def( 
                "setNum"
                , setNum_function_value
                , ( bp::arg("arg0"), bp::arg("base")=(int)(10) )
                , bp::return_self< >() );
        
        }
        { //::QByteArray::setNum
        
            typedef ::QByteArray & ( ::QByteArray::*setNum_function_type )( ::ushort,int ) ;
            setNum_function_type setNum_function_value( &::QByteArray::setNum );
            
            QByteArray_exposer.def( 
                "setNum"
                , setNum_function_value
                , ( bp::arg("arg0"), bp::arg("base")=(int)(10) )
                , bp::return_self< >() );
        
        }
        { //::QByteArray::setNum
        
            typedef ::QByteArray & ( ::QByteArray::*setNum_function_type )( int,int ) ;
            setNum_function_type setNum_function_value( &::QByteArray::setNum );
            
            QByteArray_exposer.def( 
                "setNum"
                , setNum_function_value
                , ( bp::arg("arg0"), bp::arg("base")=(int)(10) )
                , bp::return_self< >() );
        
        }
        { //::QByteArray::setNum
        
            typedef ::QByteArray & ( ::QByteArray::*setNum_function_type )( ::uint,int ) ;
            setNum_function_type setNum_function_value( &::QByteArray::setNum );
            
            QByteArray_exposer.def( 
                "setNum"
                , setNum_function_value
                , ( bp::arg("arg0"), bp::arg("base")=(int)(10) )
                , bp::return_self< >() );
        
        }
        { //::QByteArray::setNum
        
            typedef ::QByteArray & ( ::QByteArray::*setNum_function_type )( float,char,int ) ;
            setNum_function_type setNum_function_value( &::QByteArray::setNum );
            
            QByteArray_exposer.def( 
                "setNum"
                , setNum_function_value
                , ( bp::arg("arg0"), bp::arg("f")=(char)('g'), bp::arg("prec")=(int)(6) )
                , bp::return_self< >() );
        
        }
        { //::QByteArray::setNum
        
            typedef ::QByteArray & ( ::QByteArray::*setNum_function_type )( double,char,int ) ;
            setNum_function_type setNum_function_value( &::QByteArray::setNum );
            
            QByteArray_exposer.def( 
                "setNum"
                , setNum_function_value
                , ( bp::arg("arg0"), bp::arg("f")=(char)('g'), bp::arg("prec")=(int)(6) )
                , bp::return_self< >() );
        
        }
        { //::QByteArray::simplified
        
            typedef ::QByteArray ( ::QByteArray::*simplified_function_type )(  ) const;
            simplified_function_type simplified_function_value( &::QByteArray::simplified );
            
            QByteArray_exposer.def( 
                "simplified"
                , simplified_function_value );
        
        }
        { //::QByteArray::size
        
            typedef int ( ::QByteArray::*size_function_type )(  ) const;
            size_function_type size_function_value( &::QByteArray::size );
            
            QByteArray_exposer.def( 
                "size"
                , size_function_value );
        
        }
        { //::QByteArray::split
        
            typedef ::QList< QByteArray > ( ::QByteArray::*split_function_type )( char ) const;
            split_function_type split_function_value( &::QByteArray::split );
            
            QByteArray_exposer.def( 
                "split"
                , split_function_value
                , ( bp::arg("sep") ) );
        
        }
        { //::QByteArray::squeeze
        
            typedef void ( ::QByteArray::*squeeze_function_type )(  ) ;
            squeeze_function_type squeeze_function_value( &::QByteArray::squeeze );
            
            QByteArray_exposer.def( 
                "squeeze"
                , squeeze_function_value );
        
        }
        { //::QByteArray::startsWith
        
            typedef bool ( ::QByteArray::*startsWith_function_type )( ::QByteArray const & ) const;
            startsWith_function_type startsWith_function_value( &::QByteArray::startsWith );
            
            QByteArray_exposer.def( 
                "startsWith"
                , startsWith_function_value
                , ( bp::arg("a") ) );
        
        }
        { //::QByteArray::startsWith
        
            typedef bool ( ::QByteArray::*startsWith_function_type )( char ) const;
            startsWith_function_type startsWith_function_value( &::QByteArray::startsWith );
            
            QByteArray_exposer.def( 
                "startsWith"
                , startsWith_function_value
                , ( bp::arg("c") ) );
        
        }
        { //::QByteArray::startsWith
        
            typedef bool ( ::QByteArray::*startsWith_function_type )( char const * ) const;
            startsWith_function_type startsWith_function_value( &::QByteArray::startsWith );
            
            QByteArray_exposer.def( 
                "startsWith"
                , startsWith_function_value
                , ( bp::arg("c") ) );
        
        }
        { //::QByteArray::toLower
        
            typedef ::QByteArray ( ::QByteArray::*toLower_function_type )(  ) const;
            toLower_function_type toLower_function_value( &::QByteArray::toLower );
            
            QByteArray_exposer.def( 
                "toLower"
                , toLower_function_value );
        
        }
        { //::QByteArray::toUpper
        
            typedef ::QByteArray ( ::QByteArray::*toUpper_function_type )(  ) const;
            toUpper_function_type toUpper_function_value( &::QByteArray::toUpper );
            
            QByteArray_exposer.def( 
                "toUpper"
                , toUpper_function_value );
        
        }
        { //::QByteArray::trimmed
        
            typedef ::QByteArray ( ::QByteArray::*trimmed_function_type )(  ) const;
            trimmed_function_type trimmed_function_value( &::QByteArray::trimmed );
            
            QByteArray_exposer.def( 
                "trimmed"
                , trimmed_function_value );
        
        }
        { //::QByteArray::truncate
        
            typedef void ( ::QByteArray::*truncate_function_type )( int ) ;
            truncate_function_type truncate_function_value( &::QByteArray::truncate );
            
            QByteArray_exposer.def( 
                "truncate"
                , truncate_function_value
                , ( bp::arg("pos") ) );
        
        }
        QByteArray_exposer.staticmethod( "fromBase64" );
        QByteArray_exposer.staticmethod( "fromRawData" );
        QByteArray_exposer.staticmethod( "number" );
        QByteArray_exposer.def( "__copy__", &__copy__);
        QByteArray_exposer.def( "__deepcopy__", &__copy__);
        QByteArray_exposer.def( "clone", &__copy__);
        QByteArray_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::QByteArray >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        QByteArray_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::QByteArray >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
    }

}