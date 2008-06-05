#ifndef SIREGLOBAL_H
#define SIREGLOBAL_H

/** This file contains some global macros that are used
    to improve the portability of the code. */

//include qglobal.h - this includes all of the definitions
//that are needed for these macros
#include <qconfig.h>
#include <qglobal.h>

#ifdef __cplusplus

#include <qmetatype.h>

#include <boost/current_function.hpp>
//macro used to return the current file and line as a QString
#ifdef CODELOC
#undef CODELOC
#endif

#define CODELOC QString("FILE: %1, LINE: %2, FUNCTION: %3").arg(__FILE__).arg(__LINE__).arg(BOOST_CURRENT_FUNCTION)

//copy the QT_BEGIN_HEADER and QT_END_HEADER
//to SIRE_BEGIN_HEADER and SIRE_END_HEADER. This
//will allow me to change their definition at some
//future point in time without having to mess with
//Qt
#define SIRE_BEGIN_HEADER QT_BEGIN_HEADER
#define SIRE_END_HEADER QT_END_HEADER

//create keywords that control whether classes or functions are
//exposed to a scripting language
#define SIRE_EXPOSE_FUNCTION(f)  /* Exposing function #1 */
#define SIRE_EXPOSE_CLASS(c)     /* Exposing class #1 */
#define SIRE_EXPOSE_ALIAS(c,a)   /* Exposing class #1 as alias #2 */

//create the keyword used to export a symbol - this
//is a copy of Q_DECL_EXPORT, which will definitely
//be set to the correct value
#ifndef SIRE_NO_VISIBILITY_AVAILABLE
#define SIRE_EXPORT Q_DECL_EXPORT

//create the keyword to import a symbol - this is copied
//from Q_DECL_IMPORT, which will definitely be set to
//the correct value
#define SIRE_IMPORT Q_DECL_IMPORT
#else
#define SIRE_EXPORT
#define SIRE_IMPORT
#endif

//create the keyword to fix symbol visibility problems for out-of-line
//template functions
#define SIRE_OUTOFLINE_TEMPLATE Q_OUTOFLINE_TEMPLATE

//do the same for inline template functions
#define SIRE_INLINE_TEMPLATE Q_INLINE_TEMPLATE

//Sire is much more picky about what the compiler can support
//than Qt - use the results of the Qt portability tests
//to check that the C++ compiler is up to scratch

#ifdef Q_BROKEN_TEMPLATE_SPECIALIZATION
#error Sire requires a C++ compiler that can handle templates in \
       a standards compliant manner. Please use a more recent version \
       of your compiler, or try gcc >= 3.0
#endif

#ifdef QT_NO_PARTIAL_TEMPLATE_SPECIALIZATION
#error Sire requires that the C++ compiler supports \
       partial template specialization. Please use a more \
     recent version of your compiler, or try gcc >= 3.0
#endif

#ifdef Q_NO_USING_KEYWORD
#error Sire requires that the C++ compiler supports the 'using' \
       keyword. Please use a more recent version of your compiler \
     or try gcc >= 3.0
#endif

#ifdef QT_NO_MEMBER_TEMPLATES
#error Sire requires that the C++ compiler supports the use of \
       member templates. Please use a more recent version of your \
     compiler or try gcc >= 3.0
#endif

#ifdef QT_NO_THREAD
#error Sire requires that Qt was built with threading enabled. \
       Please recompile Qt with threading enabled.
#endif

#ifdef QT_NO_TEXTSTREAM
#error Sire requires that Qt was built with QTextStream enabled \
       Please recompile Qt with QTextStream enabled.
#endif

//sire requires that at least the QtCore and QtSql modules
//are available
#if  !(QT_EDITION & QT_MODULE_CORE)
#error Sire requires that the QtCore module is enabled.
#endif

#if !(QT_EDITION & QT_MODULE_SQL)
#error Sire requires that the QtSql module is enabled.
#endif

/** These functions convert a pointer into an integer.

    toInt(void*) : Returns an integer from the pointer. Returns
                   an integer of the same size as the pointer.

    toUInt(void*) : Returns an unsigned integer of the same size as
                    the pointer.

    toInt32(void*)   : Returns a 32bit integer from a pointer. This
                       will truncate if the pointer is larger than 32bits.

    toUInt32(void*)  : Same as toInt32, but returns the pointer as an
                       unsigned 32bit integer.

    toInt64(void*) : Returns a 64bit integer

    toUInt64(void*) : Returns a 64bit unsigned integer.

    @author Christopher Woods
*/
#ifdef QT_POINTER_SIZE
  #if QT_POINTER_SIZE == 4

    // Functions used if we have 32bit pointers
    inline qint32 toInt(const void *ptr)
    {
        return qint32(ptr);
    }

    inline quint32 toUInt(const void *ptr)
    {
        return quint32(ptr);
    }

    inline qint32 toInt32(const void *ptr)
    {
        return toInt(ptr);
    }

    inline quint32 toUInt32(const void *ptr)
    {
        return toUInt(ptr);
    }

    inline qint64 toInt64(const void *ptr)
    {
        return qint64( toInt(ptr) );
    }

    inline quint64 toUInt64(const void *ptr)
    {
        return quint64( toUInt(ptr) );
    }

  #elif QT_POINTER_SIZE == 8

    // Functions used if we have 64bit pointers
    inline qint64 toInt(const void *ptr)
    {
        return qint64(ptr);
    }

    inline quint64 toUInt(const void *ptr)
    {
        return quint64(ptr);
    }

    inline qint32 toInt32(const void *ptr)
    {
        return qint32( toInt(ptr) );
    }

    inline quint32 toUInt32(const void *ptr)
    {
        return qint32( toUInt(ptr) );
    }

    inline qint64 toInt64(const void *ptr)
    {
        return toInt(ptr);
    }

    inline quint64 toUInt64(const void *ptr)
    {
        return toUInt(ptr);
    }

  #else
    #fatal Invalid QT_POINTER_SIZE value stored!
  #endif
#else
  #fatal No QT_POINTER_SIZE macro defined!
#endif

//Add functions that are used to register all public
//types (this allows them to be streamed to a binary file and
//created dynamically)
#include <QMetaType>

namespace Sire
{
typedef quint32 MagicID;

MagicID getMagic(const char *type_name);

/** Enum used to register with only the magic ID */
enum MagicOnlyEnum{ MAGIC_ONLY = 1 };
/** Enum used to register without using the streaming operators */
enum NoStreamEnum{ NO_STREAM = 2 };

class RegisterMetaTypeBase
{
public:
    RegisterMetaTypeBase(MagicID magic, const char *type_name, int ID)
          : magicid(magic), typnam(type_name), id(ID)
    {}

    ~RegisterMetaTypeBase()
    {}

    MagicID magicID() const
    {
        return magicid;
    }

    const char* typeName() const
    {
        return typnam;
    }

    int ID() const
    {
        return id;
    }

private:
    MagicID magicid;
    const char *typnam;
    int id;
};

template<class T>
class MetaType
{
public:
    MetaType()
    {
        magicid = getMagic(this->typeName());
    }

    ~MetaType()
    {}

    static MagicID magicID()
    {
        return magicid;
    }

    static int ID()
    {
        return qMetaTypeId<T>();
    }

    static const char* typeName()
    {
        return QMetaType::typeName( MetaType<T>::ID() );
    }

private:
    static MagicID magicid;
};

template<class T>
MagicID MetaType<T>::magicid(0);

/** This is used to register the type 'T' - this
    needs to be called once for each public type.

    @author Christopher Woods
*/
template<class T>
class RegisterMetaType : public RegisterMetaTypeBase
{
public:

    /** Use this constructor to register a concrete class */
    RegisterMetaType()
        : RegisterMetaTypeBase( MetaType<T>::magicID(),
                                MetaType<T>::typeName(),
                                MetaType<T>::ID() )
    {
        qRegisterMetaType<T>(this->typeName());
        qRegisterMetaTypeStreamOperators<T>(this->typeName());
    }

    /** Use this construct to register a pure-virtual class */
    RegisterMetaType(MagicOnlyEnum, const char *type_name)
        : RegisterMetaTypeBase( getMagic(type_name), type_name, 0 )
    {}

    /** Use this constructor to register a class that cannot be streamed */
    RegisterMetaType(NoStreamEnum )
        : RegisterMetaTypeBase( MetaType<T>::magicID(),
                                MetaType<T>::typeName(),
                                MetaType<T>::ID() )
    {
        qRegisterMetaType<T>(this->typeName());
    }

    /** Destructor */
    ~RegisterMetaType()
    {}
};

} // namespace Sire

namespace SireStream
{
class XMLStream; // This class is used to stream objects to and from XML files

} // namespace SireStream

using Sire::MetaType;
using Sire::RegisterMetaTypeBase;
using Sire::RegisterMetaType;
using Sire::MAGIC_ONLY;
using Sire::NO_STREAM;
using Sire::MagicID;
using Sire::getMagic;
using SireStream::XMLStream;

class QDataStream;
class QTextStream;

#else  // else #ifdef __cplusplus

//copied directly from qglobal.h /////////////////
#ifndef Q_DECL_EXPORT
#  ifdef Q_OS_WIN
#    define Q_DECL_EXPORT __declspec(dllexport)
#  elif defined(QT_VISIBILITY_AVAILABLE)
#    define Q_DECL_EXPORT __attribute__((visibility("default")))
#  endif
#  ifndef Q_DECL_EXPORT
#    define Q_DECL_EXPORT
#  endif
#endif
#ifndef Q_DECL_IMPORT
#  ifdef Q_OS_WIN
#    define Q_DECL_IMPORT __declspec(dllimport)
#  else
#    define Q_DECL_IMPORT
#  endif
#endif
///////////////////////////////////////////////////

#ifndef SIRE_NO_VISIBILITY_AVAILABLE
#define SIRE_EXPORT Q_DECL_EXPORT
#define SIRE_IMPORT Q_DECL_IMPORT
#else
#define SIRE_EXPORT
#define SIRE_IMPORT
#endif

#endif // #ifdef __cplusplus

//I now define seperate SIRE_EXPORT macros for each of the different Sire libraries.
//This allows me to use 'grep' to find all of the exported symbols for each library
//(by doing 'grep -r SIREBASE_EXPORT *'
#define SIREBASE_EXPORT SIRE_EXPORT
#define SIRECAS_EXPORT SIRE_EXPORT
#define SIRECLUSTER_EXPORT SIRE_EXPORT
#define SIREDB_EXPORT SIRE_EXPORT
#define SIREERROR_EXPORT SIRE_EXPORT
#define SIREFF_EXPORT SIRE_EXPORT
#define SIREID_EXPORT SIRE_EXPORT
#define SIREIO_EXPORT SIRE_EXPORT
#define SIREMM_EXPORT SIRE_EXPORT
#define SIREMATHS_EXPORT SIRE_EXPORT
#define SIREMOL_EXPORT SIRE_EXPORT
#define SIREMOVE_EXPORT SIRE_EXPORT
#define SIREPY_EXPORT SIRE_EXPORT
#define SIREQT_EXPORT SIRE_EXPORT
#define SIRESIM_EXPORT SIRE_EXPORT
#define SIRESTREAM_EXPORT SIRE_EXPORT
#define SIRESYSTEM_EXPORT SIRE_EXPORT
#define SIRETEST_EXPORT SIRE_EXPORT
#define SIREUNITTEST_EXPORT SIRE_EXPORT
#define SIREUNITS_EXPORT SIRE_EXPORT
#define SIREVOL_EXPORT SIRE_EXPORT
#define SPIER_EXPORT SIRE_EXPORT
#define SQUIRE_EXPORT SIRE_EXPORT

//whether or not to build shared python wrappers
#define BUILD_SHARED_WRAPPERS 1

#endif // SIREGLOBAL_H
